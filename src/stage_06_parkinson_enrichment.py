"""Stage 6: matched gene-set permutation tests for Parkinson spatial signal.

The primary result is the Stage 5 L2G-weighted score tested against 10,000
matched random gene sets. Broad and stringent scores are tested in parallel as
prespecified robustness analyses. This is a gene-set null, not a spatial null.
"""

from __future__ import annotations

import json
import os
import platform
import resource
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
from nilearn import datasets, plotting
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
GENES = ROOT / "data" / "genes"
DATA_RESULTS = ROOT / "data" / "results"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
BRAIN_MAPS = ROOT / "results" / "brain_maps"
REPORTS = ROOT / "reports"
ATLAS = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.nii.gz"

N_PERMUTATIONS = 10_000
RANDOM_SEED = 20260913
NEIGHBOR_POOL_SIZE = 200
BATCH_SIZE = 250
FDR_THRESHOLD = 0.05
PRIMARY_METHOD = "weighted"
MATCH_VARIABLES = ("mean_expression", "log10_expression_variance", "log1p_probe_count")


def load_and_audit_stage_05() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame
]:
    """Load inputs and enforce the Stage 5 scale, formulas, and coverage."""
    expression = pd.read_csv(
        PROCESSED / "brain_region_gene_expression.csv", index_col="region_id"
    )
    expression.index = expression.index.astype(int)
    regions = pd.read_csv(PROCESSED / "region_metadata.csv").set_index("region_id")
    gene_metadata = pd.read_csv(PROCESSED / "gene_metadata.csv")
    gene_sets = {
        name: pd.read_csv(GENES / f"parkinson_genes_{name}.csv")
        for name in ("broad", "stringent", "weighted")
    }
    stage5 = pd.read_csv(DATA_RESULTS / "parkinson_regional_raw_scores.csv").set_index("region_id")
    donor = pd.read_csv(DATA_RESULTS / "parkinson_donor_regional_scores.csv")

    errors: list[str] = []
    if expression.shape != (138, 15_632):
        errors.append(f"Unexpected AHBA matrix dimensions: {expression.shape}")
    if set(stage5.get("analysis_scale", [])) != {"raw_scale_analysis"}:
        errors.append("Stage 5 primary scale is not labelled raw_scale_analysis")
    if set(expression.index) != set(stage5.index):
        errors.append("Stage 5 and AHBA region identifiers differ")
    if set(expression.columns) != set(gene_metadata["gene_symbol"]):
        errors.append("AHBA matrix and Stage 2 gene metadata differ")
    if donor["donor"].nunique() != 6:
        errors.append("Stage 5 donor scores do not contain six donors")

    available = set(expression.columns)
    matched = {
        name: [gene for gene in frame["gene"] if gene in available]
        for name, frame in gene_sets.items()
    }
    if matched["broad"] != matched["weighted"]:
        errors.append("Broad and weighted represented genes are not identically ordered")
    recalculated = pd.DataFrame(index=expression.index)
    recalculated["broad"] = expression[matched["broad"]].mean(axis=1)
    recalculated["stringent"] = expression[matched["stringent"]].mean(axis=1)
    weights = gene_sets["weighted"].set_index("gene").loc[
        matched["weighted"], "gene_weight"
    ].astype(float)
    recalculated["weighted"] = (
        expression[matched["weighted"]].mul(weights, axis=1).sum(axis=1)
        / weights.sum()
    )
    stage5_columns = {
        "broad": "broad_mean_expression",
        "stringent": "stringent_mean_expression",
        "weighted": "weighted_mean_expression",
    }
    for method, column in stage5_columns.items():
        if not np.allclose(recalculated[method], stage5[column]):
            errors.append(f"Stage 5 {method} formula could not be reproduced")
    if errors:
        raise RuntimeError("Stage 5 audit failed:\n- " + "\n- ".join(errors))
    return expression, regions, gene_metadata, gene_sets, stage5


def gene_covariates(
    expression: pd.DataFrame, gene_metadata: pd.DataFrame, broad_genes: list[str]
) -> pd.DataFrame:
    """Compute candidate technical matching variables for every AHBA gene."""
    metadata = gene_metadata.set_index("gene_symbol")
    covariates = pd.DataFrame(index=expression.columns)
    covariates.index.name = "gene"
    covariates["mean_expression"] = expression.mean(axis=0)
    variance = expression.var(axis=0, ddof=1)
    covariates["expression_variance"] = variance
    covariates["log10_expression_variance"] = np.log10(variance + 1e-12)
    covariates["probe_count"] = metadata.loc[covariates.index, "n_reannotated_candidate_probes"].astype(float)
    covariates["log1p_probe_count"] = np.log1p(covariates["probe_count"])
    covariates["measurable_region_count"] = expression.notna().sum(axis=0)
    covariates["is_parkinson_broad"] = covariates.index.isin(broad_genes)
    if covariates[list(MATCH_VARIABLES)].isna().any().any():
        raise RuntimeError("Matching covariates contain missing values")
    return covariates


def nearest_neighbors(
    disease_genes: list[str], background_genes: list[str], covariates: pd.DataFrame
) -> tuple[np.ndarray, pd.DataFrame]:
    """Rank background genes by Euclidean distance in standardized covariates."""
    all_values = covariates.loc[:, MATCH_VARIABLES]
    center = all_values.mean(axis=0)
    scale = all_values.std(axis=0, ddof=0)
    if (scale <= 0).any():
        raise RuntimeError(f"Constant requested matching variable: {list(scale.index[scale <= 0])}")
    disease = ((all_values.loc[disease_genes] - center) / scale).to_numpy()
    background = ((all_values.loc[background_genes] - center) / scale).to_numpy()
    order = np.argsort(cdist(disease, background, metric="euclidean"), axis=1)
    standardization = pd.DataFrame({"mean": center, "population_sd": scale})
    return order, standardization


def sample_unique_sets(
    neighbor_order: np.ndarray,
    n_background: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample one neighbor per disease gene without replacement within each set."""
    n_genes = neighbor_order.shape[0]
    sampled = np.empty((N_PERMUTATIONS, n_genes), dtype=np.int32)
    used = np.zeros(n_background, dtype=bool)
    for permutation in range(N_PERMUTATIONS):
        used.fill(False)
        for slot in rng.permutation(n_genes):
            pool = neighbor_order[slot, : min(NEIGHBOR_POOL_SIZE, n_background)]
            available = pool[~used[pool]]
            if available.size == 0:
                full = neighbor_order[slot]
                available = full[~used[full]]
            chosen = int(available[rng.integers(available.size)])
            sampled[permutation, slot] = chosen
            used[chosen] = True
    if not np.all(np.diff(np.sort(sampled, axis=1), axis=1) > 0):
        raise RuntimeError("Duplicate gene found within a matched random set")
    return sampled


def null_scores(
    background_expression: np.ndarray,
    sampled: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Calculate region × permutation scores in batches without a 3D result."""
    n_regions = background_expression.shape[0]
    output = np.empty((N_PERMUTATIONS, n_regions), dtype=np.float32)
    denominator = float(weights.sum()) if weights is not None else None
    for start in range(0, N_PERMUTATIONS, BATCH_SIZE):
        stop = min(start + BATCH_SIZE, N_PERMUTATIONS)
        values = background_expression[:, sampled[start:stop]]
        if weights is None:
            batch = values.mean(axis=2)
        else:
            batch = np.einsum("rbg,g->rb", values, weights, optimize=True) / denominator
        output[start:stop] = batch.T
    return output


def enrichment_table(
    observed: pd.Series,
    random_scores: np.ndarray,
    regions: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    """Calculate one-sided empirical gene-set permutation statistics."""
    obs = observed.to_numpy(dtype=float)
    random_mean = random_scores.mean(axis=0, dtype=np.float64)
    random_std = random_scores.std(axis=0, ddof=1, dtype=np.float64)
    z = np.divide(obs - random_mean, random_std, out=np.full_like(obs, np.nan), where=random_std > 0)
    n_ge = (random_scores >= obs[None, :]).sum(axis=0)
    empirical = (n_ge + 1) / (N_PERMUTATIONS + 1)
    fdr = multipletests(empirical, alpha=FDR_THRESHOLD, method="fdr_bh")[1]
    output = pd.DataFrame(
        {
            "region_id": observed.index,
            "region_name": regions.loc[observed.index, "region_name"].to_numpy(),
            "gene_set_version": method,
            "observed_score": obs,
            "random_mean": random_mean,
            "random_std": random_std,
            "z_score": z,
            "number_random_equal_or_greater": n_ge,
            "empirical_p": empirical,
            "fdr_p": fdr,
            "effect_size": obs - random_mean,
            "n_permutations": N_PERMUTATIONS,
            "inference_scope": "gene-set permutation significance; not a spatially independent test",
        }
    )
    return output.sort_values(["fdr_p", "z_score"], ascending=[True, False]).reset_index(drop=True)


def method_comparison(results: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    pairs = (("broad", "stringent"), ("broad", "weighted"), ("stringent", "weighted"))
    for first, second in pairs:
        a = results[first].set_index("region_id")["z_score"].sort_index()
        b = results[second].set_index("region_id")["z_score"].sort_index()
        rows.append(
            {
                "method_1": first,
                "method_2": second,
                "pearson_r": pearsonr(a, b).statistic,
                "spearman_rho": spearmanr(a, b).statistic,
                "n_regions": len(a),
            }
        )
    return pd.DataFrame(rows)


def matching_balance(
    covariates: pd.DataFrame,
    disease_genes: list[str],
    background_genes: list[str],
    sampled: np.ndarray,
    gene_set_version: str,
) -> pd.DataFrame:
    flat = sampled.ravel()
    records = []
    for variable in (*MATCH_VARIABLES, "measurable_region_count"):
        background = covariates.loc[background_genes, variable].to_numpy()
        disease = covariates.loc[disease_genes, variable].to_numpy()
        matched = background[flat]
        sd = background.std(ddof=0)
        records.append(
            {
                "gene_set_version": gene_set_version,
                "variable": variable,
                "parkinson_mean": disease.mean(),
                "all_eligible_background_mean": background.mean(),
                "matched_draw_mean": matched.mean(),
                "standardized_mean_difference_before": (disease.mean() - background.mean()) / sd if sd else np.nan,
                "standardized_mean_difference_after": (disease.mean() - matched.mean()) / sd if sd else np.nan,
                "used_for_matching": variable in MATCH_VARIABLES,
            }
        )
    return pd.DataFrame(records)


def make_matching_figures(
    covariates: pd.DataFrame,
    disease_genes: list[str],
    background_genes: list[str],
    sampled: np.ndarray,
) -> None:
    rng = np.random.default_rng(RANDOM_SEED + 1)
    flat = sampled.ravel()
    subset = flat[rng.choice(flat.size, size=min(100_000, flat.size), replace=False)]
    background = covariates.loc[background_genes]
    disease = covariates.loc[disease_genes]
    matched = background.iloc[subset]
    sns.set_theme(style="whitegrid", context="talk")

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, variable, label in (
        (axes[0], "mean_expression", "Mean expression across 138 regions"),
        (axes[1], "expression_variance", "Expression variance across 138 regions"),
    ):
        for frame, name, color in (
            (background, "Eligible AHBA background", "#777777"),
            (disease, "Parkinson broad", "#d95f02"),
            (matched, "Matched random draws", "#1b9e77"),
        ):
            sns.ecdfplot(data=frame, x=variable, label=name, color=color, ax=ax)
        ax.set(xlabel=label, ylabel="Empirical cumulative proportion")
    axes[0].legend(fontsize=10)
    if axes[1].get_legend() is not None:
        axes[1].get_legend().remove()
    fig.suptitle("Stage 6 Gene Matching: Expression Properties")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_gene_matching_expression.png", dpi=250)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for frame, name, color in (
        (background, "Eligible AHBA background", "#777777"),
        (disease, "Parkinson broad", "#d95f02"),
        (matched, "Matched random draws", "#1b9e77"),
    ):
        sns.ecdfplot(data=frame, x="probe_count", label=name, color=color, ax=axes[0])
    axes[0].set(xlabel="Reannotated candidate probe count", ylabel="Empirical cumulative proportion")
    coverage = pd.DataFrame(
        {
            "group": ["AHBA background", "Parkinson broad", "Matched draws"],
            "measurable_regions": [
                int(background["measurable_region_count"].median()),
                int(disease["measurable_region_count"].median()),
                int(matched["measurable_region_count"].median()),
            ],
        }
    )
    sns.barplot(data=coverage, x="group", y="measurable_regions", color="#4c78a8", ax=axes[1])
    axes[1].set(
        xlabel="", ylabel="Regions with a finite value",
        title="Coverage is constant after Stage 2 complete-case QC",
        ylim=(0, expression_region_count(covariates) * 1.08),
    )
    axes[1].tick_params(axis="x", rotation=20)
    fig.suptitle("Stage 6 Gene Matching: Technical Coverage")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_gene_matching_coverage.png", dpi=250)
    plt.close(fig)


def expression_region_count(covariates: pd.DataFrame) -> int:
    return int(covariates["measurable_region_count"].max())


def atlas_image(values: pd.Series) -> nib.Nifti1Image:
    atlas = nib.load(ATLAS)
    labels = np.asanyarray(atlas.dataobj).astype(np.int16)
    output = np.zeros(labels.shape, dtype=np.float32)
    for region_id, value in values.items():
        output[labels == int(region_id)] = float(value)
    return nib.Nifti1Image(output, atlas.affine, atlas.header)


def make_result_figures(
    primary: pd.DataFrame,
    primary_null: np.ndarray,
    all_results: dict[str, pd.DataFrame],
    comparison: pd.DataFrame,
    donor: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    by_region = primary.set_index("region_id")
    ordered_z = primary.sort_values("z_score", ascending=False)
    representative_ids = [
        int(ordered_z.iloc[0]["region_id"]),
        int(ordered_z.iloc[len(ordered_z) // 2]["region_id"]),
        int(ordered_z.iloc[-1]["region_id"]),
    ]
    labels = ["highest Z", "median Z rank", "lowest Z"]
    region_order = list(pd.read_csv(DATA_RESULTS / "parkinson_regional_raw_scores.csv")["region_id"])
    column_for_region = {int(region): index for index, region in enumerate(region_order)}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    for ax, region_id, rule in zip(axes, representative_ids, labels):
        row = by_region.loc[region_id]
        values = primary_null[:, column_for_region[region_id]]
        sns.histplot(values, bins=40, color="#4c78a8", ax=ax)
        ax.axvline(row["observed_score"], color="#d62728", linewidth=3, label="Observed")
        ax.set(
            title=f"{row['region_name']}\n{rule}; Z={row['z_score']:.2f}",
            xlabel="Matched-random regional score",
            ylabel="Permutations",
        )
        ax.legend()
    fig.suptitle("Parkinson Regional Null Distributions")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_parkinson_null_distributions.png", dpi=250)
    plt.close(fig)

    template = datasets.load_mni152_template(resolution=2)
    z_values = by_region["z_score"]
    max_abs = max(abs(float(z_values.min())), abs(float(z_values.max())))
    display = plotting.plot_stat_map(
        atlas_image(z_values), bg_img=template, display_mode="ortho", cut_coords=(0, -20, 5),
        cmap="RdBu_r", threshold=1e-8, vmin=-max_abs, vmax=max_abs,
        symmetric_cbar=True, colorbar=True, draw_cross=False,
        title="Parkinson Genetic Enrichment Across the Healthy Brain\nGene-set permutation Z scores",
        figure=plt.figure(figsize=(12, 5)),
    )
    display.savefig(BRAIN_MAPS / "stage_06_parkinson_zscore_map.png", dpi=250)
    display.close()

    passing = by_region["fdr_p"] < FDR_THRESHOLD
    if passing.any():
        masked = z_values.where(passing, 0.0)
        display = plotting.plot_stat_map(
            atlas_image(masked), bg_img=template, display_mode="ortho", cut_coords=(0, -20, 5),
            cmap="Reds", threshold=1e-8, vmin=0, vmax=max(float(masked.max()), 1e-8),
            symmetric_cbar=False, colorbar=True, draw_cross=False,
            title="Parkinson FDR-Significant Regional Enrichment\nGene-set permutation q < 0.05",
            figure=plt.figure(figsize=(12, 5)),
        )
    else:
        figure = plt.figure(figsize=(12, 5))
        display = plotting.plot_anat(
            template, display_mode="ortho", cut_coords=(0, -20, 5), draw_cross=False,
            title=None, figure=figure,
        )
        figure.text(
            0.01, 0.98,
            "Parkinson FDR-Significant Regional Enrichment\nNo regions passed gene-set permutation q < 0.05",
            ha="left", va="top", fontsize=15, color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 3},
        )
    display.savefig(BRAIN_MAPS / "stage_06_parkinson_fdr_map.png", dpi=250)
    display.close()

    fig, ax = plt.subplots(figsize=(9, 8))
    plot = primary.copy()
    plot["passes_fdr"] = plot["fdr_p"] < FDR_THRESHOLD
    sns.scatterplot(
        data=plot, x="random_mean", y="observed_score", hue="passes_fdr",
        palette={False: "#777777", True: "#d62728"}, s=70, ax=ax,
    )
    limits = [min(plot["random_mean"].min(), plot["observed_score"].min()),
              max(plot["random_mean"].max(), plot["observed_score"].max())]
    ax.plot(limits, limits, linestyle="--", color="black", label="Identity")
    ax.set(
        title="Observed vs Random Regional Expression",
        xlabel="Matched-random mean", ylabel="Observed Parkinson weighted score",
    )
    ax.legend(title="Gene-set q < 0.05")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_observed_vs_random.png", dpi=250)
    plt.close(fig)

    top = primary.nlargest(20, "z_score")
    colors = np.where(top["fdr_p"] < FDR_THRESHOLD, "#d62728", "#7f8c8d")
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.barh(top["region_name"], top["z_score"], color=colors)
    ax.invert_yaxis()
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set(title="Regional Enrichment Ranking", xlabel="Gene-set permutation Z score", ylabel="AAL3 region")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#d62728", label="q < 0.05"), Patch(color="#7f8c8d", label="q ≥ 0.05")])
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_parkinson_enrichment_ranking.png", dpi=250)
    plt.close(fig)

    methods = ["broad", "stringent", "weighted"]
    pearson = pd.DataFrame(np.eye(3), index=methods, columns=methods)
    spearman = pearson.copy()
    for _, row in comparison.iterrows():
        a, b = row["method_1"], row["method_2"]
        pearson.loc[a, b] = pearson.loc[b, a] = row["pearson_r"]
        spearman.loc[a, b] = spearman.loc[b, a] = row["spearman_rho"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for matrix, title, ax in ((pearson, "Pearson", axes[0]), (spearman, "Spearman", axes[1])):
        sns.heatmap(matrix, annot=True, fmt=".3f", cmap="vlag", center=0, vmin=-1, vmax=1,
                    square=True, cbar_kws={"label": "Correlation"}, ax=ax)
        ax.set_title(title)
    fig.suptitle("Stage 6 Null-Model Method Robustness")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_method_robustness_heatmap.png", dpi=250)
    plt.close(fig)

    top_ids = list(primary.nlargest(10, "z_score")["region_id"])
    names = by_region.loc[top_ids, "region_name"]
    selected = donor.loc[donor["region_id"].isin(top_ids)]
    heat = selected.pivot(index="donor", columns="region_id", values="weighted_score").reindex(columns=top_ids)
    heat.columns = [names.loc[x] for x in heat.columns]
    fig, ax = plt.subplots(figsize=(16, 7))
    sns.heatmap(heat, cmap="viridis", annot=True, fmt=".2f", linewidths=0.3,
                cbar_kws={"label": "Observed donor weighted score"}, ax=ax)
    ax.set(
        title="Donor Robustness for Top Gene-Set Z-Score Regions",
        xlabel="Regions ordered by pooled Z score", ylabel="Donor",
    )
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_06_donor_robustness.png", dpi=250)
    plt.close(fig)


def markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for _, row in frame.iterrows():
        cells = []
        for column in columns:
            value = row[column]
            cells.append(f"{value:.4g}" if isinstance(value, (float, np.floating)) else str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_reports(
    balance: pd.DataFrame,
    results: dict[str, pd.DataFrame],
    elapsed: float,
    peak_mb: float,
    counts: dict[str, int],
) -> None:
    bias = f"""# Stage 6 gene-set bias assessment

## Candidate confounders

The Stage 5 matrix contains scaled-robust-sigmoid normalized expression for 15,632
genes across 138 retained AAL3 regions. Before random sampling, Parkinson broad-set
genes were compared with the eligible AHBA background on global mean expression,
regional expression variance, reannotated candidate probe count, and the number of
regions with a finite measurement.

The matching model uses mean expression, log10 variance, and log1p probe count.
Each variable is standardized over the AHBA universe and given equal unit-variance
weight in Euclidean distance. These variables address baseline abundance, spatial
variability, and microarray representation. They summarize measurement properties
without matching the region-by-region pattern being tested.

Coverage was inspected but not used because Stage 2 complete-case filtering makes
it exactly 138 regions for every eligible gene. Gene length is not available in the
existing Stage 2 metadata and was not inferred or silently substituted. Candidate
probe count is an imperfect technical proxy: `abagen` exposes the retained gene,
while Stage 2 records reannotated candidate probes rather than claiming the exact
selected probe identity.

## Background and matching

All 123 AHBA-represented broad Parkinson genes were excluded from the common
background, leaving 15,509 genes. For each disease gene, the algorithm computes its
200 nearest background neighbors, randomizes disease-gene processing order within
each permutation, and draws an unused candidate. It expands beyond 200 only if a
pool is exhausted. This yields unique genes within every set while keeping the
procedure transparent. The same broad matched sets are reused for the weighted null.

## Balance results

{markdown_table(balance)}

Standardized mean differences after matching are diagnostics rather than outcome
tests. The figures show complete empirical distributions and a fixed 100,000-draw
subsample of the matched distribution for plotting efficiency.
"""
    (REPORTS / "stage_06_gene_set_bias_assessment.md").write_text(bias, encoding="utf-8")

    primary = results[PRIMARY_METHOD]
    n_pass = int((primary["fdr_p"] < FDR_THRESHOLD).sum())
    interpretation = f"""# Stage 6 interpretation

## What this analysis means

For each region, the observed Stage 5 L2G-weighted Parkinson score was compared with
{N_PERMUTATIONS:,} unique-gene random sets matched on global expression mean,
expression variance, and candidate probe count. A positive effect size or Z score
means the Parkinson-prioritised genes have higher expression than those matched
sets in that region. The one-sided empirical p-value includes the prespecified +1
correction, and Benjamini-Hochberg correction is applied across 138 regions at
q < {FDR_THRESHOLD:.2f}. {n_pass} regions pass that threshold in the primary weighted analysis.

This result is labelled **gene-set permutation significance**. It does not show that
a region causes Parkinson disease, initiates pathology, is specific to Parkinson,
or will be affected in an individual. L2G scores remain prioritisation evidence,
and healthy adult post-mortem expression is not a direct measurement of disease.

## Assumptions and limitations

- The selected matching summaries adequately control major measurable gene-level
  abundance, variance, and probe-representation differences.
- Matching does not control unrecorded properties such as gene length, GC content,
  cell-type specificity, network degree, or all probe-design effects.
- Excluding Parkinson genes from the background avoids direct contamination but
  slightly changes the eligible-gene universe.
- Fixed Parkinson L2G weights are carried to their gene-specific matched replacements,
  preserving the exact weight vector without inventing random weights.
- The six AHBA donors are few and unevenly sampled; four are predominantly left-sided.
- Benjamini-Hochberg correction addresses multiple regional tests under its standard
  assumptions, but the gene-set permutations do not remove spatial autocorrelation
  among neighboring atlas parcels. A spatially informed null remains a separate
  sensitivity analysis for later work.
- Results depend on the AAL3 parcellation, Stage 2 structural normalization, Stage 4
  gene thresholds, and Open Targets release.
"""
    (REPORTS / "stage_06_interpretation.md").write_text(interpretation, encoding="utf-8")

    performance = f"""# Stage 6 performance

- Regions: 138
- AHBA genes: 15,632
- Eligible non-Parkinson background genes: 15,509
- Represented genes: broad={counts['broad']}, stringent={counts['stringent']}, weighted={counts['weighted']}
- Permutations per analysis: {N_PERMUTATIONS:,}
- Total runtime: {elapsed:.2f} seconds
- Approximate peak resident memory: {peak_mb:.1f} MB
- Batch size: {BATCH_SIZE} permutations
- Memory strategy: random indices are stored as compact int32 matrices; null scores
  are stored as three two-dimensional permutation × region float32 arrays. Temporary
  region × batch × gene arrays are released each batch; no full three-dimensional
  permutation object is retained.
- Python: {platform.python_version()}
- Platform: {platform.platform()}
- Processor: {platform.processor() or 'not reported'}
- Logical CPU count: {os.cpu_count()}
- Random seed: {RANDOM_SEED}
"""
    (REPORTS / "stage_06_performance.md").write_text(performance, encoding="utf-8")


def peak_memory_mb() -> float:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (1024 * 1024) if sys.platform == "darwin" else value / 1024


def print_summary(
    results: dict[str, pd.DataFrame], comparison: pd.DataFrame, counts: dict[str, int],
    donor: pd.DataFrame,
) -> None:
    print("\nPARKINSON ENRICHMENT SUMMARY")
    print("Gene-set versions: broad, stringent, weighted (primary)")
    print(f"Number of genes: broad={counts['broad']}, stringent={counts['stringent']}, weighted={counts['weighted']}")
    print(f"Number of permutations: {N_PERMUTATIONS}")
    print("Matching strategy: random draw from 200 nearest neighbors on standardized mean expression, log variance, and log probe count; unique within set")
    for method, table in results.items():
        top_z = ", ".join(table.nlargest(5, "z_score")["region_name"])
        top_fdr = ", ".join(table.sort_values(["fdr_p", "z_score"], ascending=[True, False]).head(5)["region_name"])
        print(f"{method} top regions by Z score: {top_z}")
        print(f"{method} top regions by FDR: {top_fdr}")
        print(f"{method} regions p < 0.05: {(table['empirical_p'] < 0.05).sum()}")
        print(f"{method} regions FDR < 0.05: {(table['fdr_p'] < 0.05).sum()}")
    print("Broad/stringent/weighted spatial correlations:")
    for _, row in comparison.iterrows():
        print(f"  {row['method_1']} vs {row['method_2']}: Pearson={row['pearson_r']:.3f}, Spearman={row['spearman_rho']:.3f}")
    top_ids = set(results[PRIMARY_METHOD].nlargest(10, "z_score")["region_id"])
    subset = donor.loc[donor["region_id"].isin(top_ids)]
    observed = subset.groupby("region_id")["weighted_score"].count()
    print(f"Donor robustness summary: top regions have {observed.min()}–{observed.max()} observed donors; no donor-specific inference performed")
    print("STOP: no pathway, cross-disease, prediction, pathology-training, or spatial-null analysis performed.")


def main() -> None:
    started = time.perf_counter()
    for directory in (DATA_RESULTS, TABLES, FIGURES, BRAIN_MAPS, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    expression, regions, gene_metadata, gene_sets, stage5 = load_and_audit_stage_05()
    available = set(expression.columns)
    genes = {
        name: [gene for gene in frame["gene"] if gene in available]
        for name, frame in gene_sets.items()
    }
    covariates = gene_covariates(expression, gene_metadata, genes["broad"])
    background_genes = [gene for gene in expression.columns if gene not in set(genes["broad"])]
    background_expression = expression[background_genes].to_numpy(dtype=np.float32)

    broad_neighbors, standardization = nearest_neighbors(genes["broad"], background_genes, covariates)
    stringent_neighbors, _ = nearest_neighbors(genes["stringent"], background_genes, covariates)
    rng = np.random.default_rng(RANDOM_SEED)
    broad_sets = sample_unique_sets(broad_neighbors, len(background_genes), rng)
    stringent_sets = sample_unique_sets(stringent_neighbors, len(background_genes), rng)
    weights = gene_sets["weighted"].set_index("gene").loc[genes["broad"], "gene_weight"].to_numpy(dtype=np.float32)

    nulls = {
        "broad": null_scores(background_expression, broad_sets),
        "stringent": null_scores(background_expression, stringent_sets),
        "weighted": null_scores(background_expression, broad_sets, weights=weights),
    }
    stage5_columns = {
        "broad": "broad_mean_expression",
        "stringent": "stringent_mean_expression",
        "weighted": "weighted_mean_expression",
    }
    results = {
        method: enrichment_table(stage5[column], nulls[method], regions, method)
        for method, column in stage5_columns.items()
    }
    all_results = pd.concat(results.values(), ignore_index=True)
    primary = results[PRIMARY_METHOD]
    comparison = method_comparison(results)
    broad_balance = matching_balance(
        covariates, genes["broad"], background_genes, broad_sets, "broad"
    )
    stringent_balance = matching_balance(
        covariates, genes["stringent"], background_genes, stringent_sets, "stringent"
    )
    weighted_balance = broad_balance.copy()
    weighted_balance["gene_set_version"] = "weighted"
    balance = pd.concat(
        [broad_balance, stringent_balance, weighted_balance], ignore_index=True
    )
    donor = pd.read_csv(DATA_RESULTS / "parkinson_donor_regional_scores.csv")

    primary.to_csv(DATA_RESULTS / "parkinson_regional_enrichment.csv", index=False)
    all_results.to_csv(DATA_RESULTS / "parkinson_regional_enrichment_all_methods.csv", index=False)
    comparison.to_csv(DATA_RESULTS / "parkinson_null_model_method_comparison.csv", index=False)
    balance.to_csv(DATA_RESULTS / "stage_06_gene_matching_balance.csv", index=False)
    primary.head(25).to_csv(TABLES / "parkinson_top_enriched_regions.csv", index=False)

    make_matching_figures(covariates, genes["broad"], background_genes, broad_sets)
    make_result_figures(primary, nulls[PRIMARY_METHOD], results, comparison, donor)

    parameters = {
        "random_seed": RANDOM_SEED,
        "n_permutations_per_gene_set": N_PERMUTATIONS,
        "gene_set_versions": ["broad", "stringent", "weighted"],
        "primary_gene_set_version": PRIMARY_METHOD,
        "genes_per_set": {name: len(value) for name, value in genes.items()},
        "matching_variables": list(MATCH_VARIABLES),
        "variables_inspected_but_not_matched": {
            "measurable_region_count": "constant at 138 after Stage 2 complete-case QC",
            "gene_length": "not available in existing project metadata",
        },
        "matching_method": "Euclidean nearest neighbors after unit-variance standardization",
        "neighbor_pool_size": NEIGHBOR_POOL_SIZE,
        "sampling_algorithm": "random disease-gene order; random unused candidate from nearest-neighbor pool; without replacement within permutation",
        "background": "all AHBA genes excluding all represented broad Parkinson genes",
        "background_gene_count": len(background_genes),
        "weighted_null": "reuse broad matched replacements and carry the corresponding fixed Parkinson L2G weight to each replacement",
        "batch_size": BATCH_SIZE,
        "empirical_test": "one-sided greater-or-equal with +1 numerator and denominator correction",
        "multiple_testing": "Benjamini-Hochberg across 138 regions separately for each gene-set version",
        "fdr_threshold": FDR_THRESHOLD,
        "standardization_parameters": {
            variable: {
                "mean": float(standardization.loc[variable, "mean"]),
                "population_sd": float(standardization.loc[variable, "population_sd"]),
            }
            for variable in MATCH_VARIABLES
        },
        "uniqueness_verified": {
            "broad_and_weighted": bool(np.all(np.diff(np.sort(broad_sets, axis=1), axis=1) > 0)),
            "stringent": bool(np.all(np.diff(np.sort(stringent_sets, axis=1), axis=1) > 0)),
        },
    }
    (DATA_RESULTS / "stage_06_permutation_parameters.json").write_text(
        json.dumps(parameters, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    elapsed = time.perf_counter() - started
    write_reports(balance, results, elapsed, peak_memory_mb(), {k: len(v) for k, v in genes.items()})
    print_summary(results, comparison, {k: len(v) for k, v in genes.items()}, donor)


if __name__ == "__main__":
    main()
