"""Stage 5: observed Parkinson-associated gene expression across healthy AHBA.

This stage calculates regional expression profiles from the broad, stringent,
and L2G-weighted Parkinson gene sets. It deliberately performs no random-set
analysis, null modelling, empirical inference, FDR, or cross-disease analysis.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
from nilearn import datasets, plotting
from scipy.stats import pearsonr


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
GENES = ROOT / "data" / "genes"
INTERMEDIATE = ROOT / "data" / "intermediate"
DATA_RESULTS = ROOT / "data" / "results"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
BRAIN_MAPS = ROOT / "results" / "brain_maps"
REPORTS = ROOT / "reports"
ATLAS = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.nii.gz"

EXPRESSION = PROCESSED / "brain_region_gene_expression.csv"
REGION_METADATA = PROCESSED / "region_metadata.csv"
GENE_METADATA = PROCESSED / "gene_metadata.csv"
GENE_FILES = {
    "broad": GENES / "parkinson_genes_broad.csv",
    "stringent": GENES / "parkinson_genes_stringent.csv",
    "weighted": GENES / "parkinson_genes_weighted.csv",
}
TOP_REGIONS = 25


def load_and_validate() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame
]:
    """Load inputs, verify identifier compatibility, and record every match."""
    expression = pd.read_csv(EXPRESSION, index_col="region_id")
    expression.index = expression.index.astype(int)
    regions = pd.read_csv(REGION_METADATA).set_index("region_id")
    gene_metadata = pd.read_csv(GENE_METADATA)
    gene_sets = {name: pd.read_csv(path) for name, path in GENE_FILES.items()}

    errors: list[str] = []
    if not expression.index.is_unique or not expression.columns.is_unique:
        errors.append("AHBA matrix identifiers are not unique")
    if not np.isfinite(expression.to_numpy()).all():
        errors.append("Pooled AHBA matrix contains non-finite values")
    if set(expression.columns) != set(gene_metadata["gene_symbol"]):
        errors.append("Stage 2 gene metadata does not match the AHBA matrix columns")
    retained = regions.loc[regions["retained_in_main_matrix"].astype(bool)]
    if set(expression.index) != set(retained.index):
        errors.append("Retained region metadata does not match the AHBA matrix rows")
    for name, frame in gene_sets.items():
        if frame["gene"].isna().any() or not frame["gene"].is_unique:
            errors.append(f"{name} Parkinson gene identifiers are missing or duplicated")
        if name == "weighted" and (
            frame["gene_weight"].isna().any() or not frame["gene_weight"].between(0, 1).all()
        ):
            errors.append("Weighted gene set has invalid L2G weights")
    if errors:
        raise RuntimeError("Stage 5 input audit failed:\n- " + "\n- ".join(errors))

    available = set(expression.columns)
    audit_rows: list[dict[str, object]] = []
    for set_name, frame in gene_sets.items():
        for _, row in frame.iterrows():
            audit_rows.append(
                {
                    "gene_set": set_name,
                    "gene": row["gene"],
                    "ensembl_id": row.get("ensembl_id"),
                    "present_in_ahba": row["gene"] in available,
                    "match_method": "exact approved gene symbol",
                    "gene_weight": row.get("gene_weight", 1.0),
                    "missing_handling": (
                        "included" if row["gene"] in available else "excluded with identifier retained in audit"
                    ),
                }
            )
    audit = pd.DataFrame(audit_rows)
    return expression, retained, gene_metadata, gene_sets, audit


def present_genes(frame: pd.DataFrame, expression: pd.DataFrame) -> list[str]:
    """Return gene-set order restricted to exact symbols present in AHBA."""
    available = set(expression.columns)
    return [gene for gene in frame["gene"] if gene in available]


def calculate_scores(
    expression: pd.DataFrame, gene_sets: dict[str, pd.DataFrame]
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Calculate the three prespecified regional profiles on a supplied scale."""
    matched = {name: present_genes(frame, expression) for name, frame in gene_sets.items()}
    broad = expression[matched["broad"]].mean(axis=1)
    stringent = expression[matched["stringent"]].mean(axis=1)

    weight_table = gene_sets["weighted"].set_index("gene")
    weights = weight_table.loc[matched["weighted"], "gene_weight"].astype(float)
    weighted_values = expression[matched["weighted"]].mul(weights, axis=1).sum(axis=1) / weights.sum()
    scores = pd.DataFrame(
        {
            "broad_mean_expression": broad,
            "stringent_mean_expression": stringent,
            "weighted_mean_expression": weighted_values,
            "n_broad_genes_present": len(matched["broad"]),
            "n_stringent_genes_present": len(matched["stringent"]),
            "n_weighted_genes_present": len(matched["weighted"]),
        },
        index=expression.index,
    )
    scores.index.name = "region_id"
    return scores, matched


def add_region_names(scores: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    output = scores.join(regions[["region_name", "atlas_label", "broad_system"]], how="left")
    output["analysis_scale"] = "raw_scale_analysis"
    columns = ["region_name"] + [c for c in output.columns if c != "region_name"]
    return output[columns].reset_index()


def standardized_expression(expression: pd.DataFrame) -> pd.DataFrame:
    """Standardize each gene across the 138 pooled regions (population SD)."""
    means = expression.mean(axis=0)
    scales = expression.std(axis=0, ddof=0)
    if (scales <= 0).any():
        bad = list(scales.index[scales <= 0])
        raise RuntimeError(f"Cannot standardize zero-variance AHBA genes: {bad[:10]}")
    return (expression - means) / scales


def ranking_table(
    regional: pd.DataFrame, score_column: str, gene_count: int
) -> pd.DataFrame:
    ranked = regional.sort_values(score_column, ascending=False).head(TOP_REGIONS).copy()
    output = pd.DataFrame(
        {
            "rank": np.arange(1, len(ranked) + 1),
            "region": ranked["region_name"].to_numpy(),
            "region_id": ranked["region_id"].to_numpy(),
            "score": ranked[score_column].to_numpy(),
            "number_of_genes": gene_count,
            "interpretation": "observed Parkinson-associated gene expression",
        }
    )
    return output


def method_correlations(scores: pd.DataFrame) -> pd.DataFrame:
    columns = {
        "broad": "broad_mean_expression",
        "stringent": "stringent_mean_expression",
        "weighted": "weighted_mean_expression",
    }
    return scores[list(columns.values())].rename(columns={v: k for k, v in columns.items()}).corr()


def sensitivity_table(
    raw: pd.DataFrame, standardized: pd.DataFrame, regions: pd.DataFrame
) -> pd.DataFrame:
    columns = {
        "broad": "broad_mean_expression",
        "stringent": "stringent_mean_expression",
        "weighted": "weighted_mean_expression",
    }
    records: list[dict[str, object]] = []
    for method, column in columns.items():
        correlation = pearsonr(raw[column], standardized[column]).statistic
        for region_id in raw.index:
            records.append(
                {
                    "region_id": region_id,
                    "region_name": regions.loc[region_id, "region_name"],
                    "method": method,
                    "raw_scale_score": raw.loc[region_id, column],
                    "gene_standardized_score": standardized.loc[region_id, column],
                    "raw_scale_label": "raw_scale_analysis",
                    "standardized_scale_label": "gene_standardized_analysis",
                    "raw_vs_gene_standardized_pearson_r": correlation,
                    "primary_or_sensitivity": "primary raw-scale plus standardized sensitivity",
                }
            )
    return pd.DataFrame(records)


def gene_contributions(
    expression: pd.DataFrame,
    regional: pd.DataFrame,
    regions: pd.DataFrame,
    weighted_genes: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate additive weighted-score contributions in the top ten regions."""
    genes = present_genes(weighted_genes, expression)
    weights = weighted_genes.set_index("gene").loc[genes, "gene_weight"].astype(float)
    top_ids = regional.nlargest(10, "weighted_mean_expression").index
    records: list[dict[str, object]] = []
    for region_id in top_ids:
        values = expression.loc[region_id, genes].astype(float)
        additive = values * weights / weights.sum()
        share = additive / additive.sum()
        order = additive.sort_values(ascending=False)
        for rank, gene in enumerate(order.index, start=1):
            records.append(
                {
                    "region_id": region_id,
                    "region": regions.loc[region_id, "region_name"],
                    "regional_weighted_score": regional.loc[region_id, "weighted_mean_expression"],
                    "contribution_rank": rank,
                    "gene": gene,
                    "expression": values[gene],
                    "l2g_score": weights[gene],
                    "weighted_contribution": additive[gene],
                    "normalized_contribution": share[gene],
                }
            )
    return pd.DataFrame(records)


def donor_scores(
    gene_sets: dict[str, pd.DataFrame], regions: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate donor profiles, retaining unobserved donor-region rows as NaN."""
    paths = sorted(INTERMEDIATE.glob("donor_*_expression.csv"))
    if not paths:
        raise FileNotFoundError(
            "Stage 2 donor matrices are required; rerun src/stage_02_preprocess.py"
        )
    frames: list[pd.DataFrame] = []
    for path in paths:
        donor = path.stem.removeprefix("donor_").removesuffix("_expression")
        matrix = pd.read_csv(path, index_col="region_id")
        matrix.index = matrix.index.astype(int)
        matched = {name: present_genes(frame, matrix) for name, frame in gene_sets.items()}
        broad = matrix[matched["broad"]].mean(axis=1, skipna=True)
        stringent = matrix[matched["stringent"]].mean(axis=1, skipna=True)
        weights = gene_sets["weighted"].set_index("gene").loc[matched["weighted"], "gene_weight"].astype(float)
        values = matrix[matched["weighted"]]
        valid_weight = values.notna().mul(weights, axis=1).sum(axis=1)
        weighted = values.mul(weights, axis=1).sum(axis=1, min_count=1) / valid_weight.replace(0, np.nan)
        frame = pd.DataFrame(
            {
                "donor": donor,
                "region_id": matrix.index,
                "region": regions.loc[matrix.index, "region_name"].to_numpy(),
                "broad_score": broad.to_numpy(),
                "stringent_score": stringent.to_numpy(),
                "weighted_score": weighted.to_numpy(),
                "n_broad_genes_observed": matrix[matched["broad"]].notna().sum(axis=1).to_numpy(),
                "n_stringent_genes_observed": matrix[matched["stringent"]].notna().sum(axis=1).to_numpy(),
                "n_weighted_genes_observed": matrix[matched["weighted"]].notna().sum(axis=1).to_numpy(),
            }
        )
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)

    records: list[dict[str, object]] = []
    donors = sorted(combined["donor"].unique())
    for i, donor_a in enumerate(donors):
        left = combined.loc[combined["donor"] == donor_a].set_index("region_id")
        for donor_b in donors[i + 1 :]:
            right = combined.loc[combined["donor"] == donor_b].set_index("region_id")
            for method in ("broad", "stringent", "weighted"):
                column = f"{method}_score"
                paired = pd.concat([left[column], right[column]], axis=1, keys=["a", "b"]).dropna()
                value = pearsonr(paired["a"], paired["b"]).statistic if len(paired) >= 3 else np.nan
                records.append(
                    {
                        "donor_1": donor_a,
                        "donor_2": donor_b,
                        "score_method": method,
                        "pearson_r": value,
                        "n_common_regions": len(paired),
                    }
                )
    return combined, pd.DataFrame(records)


def atlas_score_image(scores: pd.Series) -> nib.Nifti1Image:
    atlas = nib.load(ATLAS)
    labels = np.asanyarray(atlas.dataobj).astype(np.int16)
    output = np.zeros(labels.shape, dtype=np.float32)
    for region_id, score in scores.items():
        output[labels == int(region_id)] = float(score)
    return nib.Nifti1Image(output, atlas.affine, atlas.header)


def make_brain_maps(regional: pd.DataFrame) -> None:
    mappings = [
        ("broad_mean_expression", "Broad Parkinson gene set", "stage_05_parkinson_broad_expression.png"),
        ("stringent_mean_expression", "Stringent Parkinson gene set", "stage_05_parkinson_stringent_expression.png"),
        ("weighted_mean_expression", "L2G-weighted Parkinson gene set", "stage_05_parkinson_weighted_expression.png"),
    ]
    values = regional[[x[0] for x in mappings]].to_numpy()
    vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
    threshold = max(1e-8, vmin - 1e-8)
    anatomical_template = datasets.load_mni152_template(resolution=2)
    for column, subtitle, filename in mappings:
        image = atlas_score_image(regional.set_index("region_id")[column])
        display = plotting.plot_stat_map(
            image,
            bg_img=anatomical_template,
            display_mode="ortho",
            cut_coords=(0, -20, 5),
            cmap="viridis",
            threshold=threshold,
            vmin=vmin,
            vmax=vmax,
            symmetric_cbar=False,
            colorbar=True,
            annotate=True,
            draw_cross=False,
            black_bg=False,
            title=f"Parkinson-Associated Gene Expression Across the Healthy Brain\n{subtitle}",
            figure=plt.figure(figsize=(12, 5)),
        )
        display.savefig(BRAIN_MAPS / filename, dpi=250)
        display.close()


def make_figures(
    regional: pd.DataFrame,
    method_corr: pd.DataFrame,
    donor: pd.DataFrame,
    donor_corr: pd.DataFrame,
    contributions: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    top = regional.nlargest(25, "weighted_mean_expression")
    fig, ax = plt.subplots(figsize=(11, 10))
    sns.barplot(data=top, x="weighted_mean_expression", y="region_name", color="#4c78a8", ax=ax)
    ax.set(
        title="Parkinson Regional Expression Ranking",
        xlabel="Observed L2G-weighted mean expression",
        ylabel="AAL3 region",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_05_parkinson_top_regions.png", dpi=250)
    plt.close(fig)

    weighted_wide = donor.pivot(index="donor", columns="region", values="weighted_score")
    region_order = regional.sort_values("weighted_mean_expression", ascending=False)["region_name"]
    weighted_wide = weighted_wide.reindex(columns=region_order)
    fig, axes = plt.subplots(2, 1, figsize=(18, 10), gridspec_kw={"height_ratios": [2.2, 1]})
    sns.heatmap(
        weighted_wide,
        cmap="viridis",
        xticklabels=False,
        cbar_kws={"label": "Observed weighted expression"},
        ax=axes[0],
    )
    axes[0].set(title="Parkinson Spatial Signal: Donor Concordance", xlabel="Regions ordered by pooled weighted score", ylabel="Donor")
    weighted_pairs = donor_corr.loc[donor_corr["score_method"] == "weighted"]
    donors = sorted(donor["donor"].unique())
    matrix = pd.DataFrame(np.eye(len(donors)), index=donors, columns=donors)
    for _, row in weighted_pairs.iterrows():
        matrix.loc[row["donor_1"], row["donor_2"]] = row["pearson_r"]
        matrix.loc[row["donor_2"], row["donor_1"]] = row["pearson_r"]
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="vlag", center=0, vmin=-1, vmax=1,
                cbar_kws={"label": "Pearson r"}, ax=axes[1])
    axes[1].set(title="Pairwise donor correlation: weighted regional profile", xlabel="Donor", ylabel="Donor")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_05_parkinson_donor_concordance.png", dpi=250)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(method_corr, annot=True, fmt=".3f", cmap="viridis", vmin=0, vmax=1,
                square=True, cbar_kws={"label": "Pearson r"}, ax=ax)
    ax.set(title="Parkinson Spatial Signal: Scoring-Method Comparison")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_05_parkinson_method_comparison.png", dpi=250)
    plt.close(fig)

    top_five = list(regional.nlargest(5, "weighted_mean_expression")["region_name"])
    subset = contributions.loc[contributions["region"].isin(top_five)]
    top_genes = list(
        subset.groupby("gene")["normalized_contribution"].max().nlargest(20).index
    )
    heat = subset.loc[subset["gene"].isin(top_genes)].pivot(
        index="region", columns="gene", values="normalized_contribution"
    ).reindex(top_five)
    heat = heat.reindex(columns=heat.max().sort_values(ascending=False).index)
    fig, ax = plt.subplots(figsize=(16, 6))
    sns.heatmap(heat, cmap="mako", cbar_kws={"label": "Fraction of weighted regional score"}, ax=ax)
    ax.set(title="Top Region Gene Contribution", xlabel="Top contributing Parkinson genes", ylabel="Top weighted-score regions")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_05_parkinson_gene_contributions.png", dpi=250)
    plt.close(fig)


def write_methods(
    expression: pd.DataFrame,
    regions: pd.DataFrame,
    gene_sets: dict[str, pd.DataFrame],
    audit: pd.DataFrame,
    donor_corr: pd.DataFrame,
    sensitivity: pd.DataFrame,
) -> None:
    counts = audit.groupby("gene_set")["present_in_ahba"].agg(["size", "sum"])
    absent = {
        name: sorted(audit.loc[(audit["gene_set"] == name) & ~audit["present_in_ahba"], "gene"])
        for name in gene_sets
    }
    sensitivity_r = sensitivity.groupby("method")["raw_vs_gene_standardized_pearson_r"].first()
    mean_donor = donor_corr.groupby("score_method")["pearson_r"].mean()
    content = rf"""# Stage 5 methods: Parkinson × AHBA observed spatial signal

## Scope and interpretation

This stage describes the observed spatial expression of genes prioritised for
Parkinson disease in the healthy adult AHBA reference. A high score means only
that the selected Parkinson-associated genes have higher expression on the Stage 2
scale in that region. No random gene sets, null distributions, p-values, FDR,
spatial permutations, pathways, or other diseases were analysed.

## Inputs and identifier matching

- AHBA matrix: {expression.shape[0]} retained AAL3 regions × {expression.shape[1]:,} genes.
- Broad set: {counts.loc['broad', 'size']} total, {counts.loc['broad', 'sum']} exact-symbol matches.
- Stringent set: {counts.loc['stringent', 'size']} total, {counts.loc['stringent', 'sum']} exact-symbol matches.
- Weighted set: {counts.loc['weighted', 'size']} total, {counts.loc['weighted', 'sum']} exact-symbol matches.
- Matching used the approved gene-symbol column from Stage 4 against the columns of
  `brain_region_gene_expression.csv`. `gene_metadata.csv` was checked to have the
  same {expression.shape[1]:,}-symbol universe. No identifier conversion was needed.
- Every included and missing identifier is recorded in
  `data/results/parkinson_gene_match_audit.csv`; nothing was silently discarded.

Missing broad/weighted genes ({len(absent['broad'])}): {', '.join(absent['broad'])}.

Missing stringent genes ({len(absent['stringent'])}): {', '.join(absent['stringent'])}.

## Expression scale

Stage 2 already applied scaled robust sigmoid normalization across genes within each
sample and then across matched samples for each gene within donor, separately for
cortex, subcortex/brainstem, and cerebellum. Samples were averaged within region and
donor, then available donors were equally averaged. Values therefore lie in [0, 1]
and are neither raw microarray intensities nor z-scores. The primary analysis,
labelled `raw_scale_analysis`, preserves this Stage 2 scale without another
transformation.

The separate `gene_standardized_analysis` sensitivity analysis transforms each gene
across the {expression.shape[0]} pooled regions as

$$Z_{{rg}} = (E_{{rg}} - \bar E_g) / \sigma_g,$$

using the population standard deviation (`ddof=0`). It is not mixed with or used to
replace the primary score.

## Regional scores

For region $r$, broad genes $G_B$, and stringent genes $G_S$:

$$S_B(r) = \frac{{1}}{{|G_B|}} \sum_{{g \in G_B}} E_{{rg}},$$

$$S_S(r) = \frac{{1}}{{|G_S|}} \sum_{{g \in G_S}} E_{{rg}}.$$

For weighted genes $G_W$ and their unmodified Stage 4 maximum L2G score $w_g$:

$$S_W(r) = \frac{{\sum_{{g \in G_W}} E_{{rg}} w_g}}{{\sum_{{g \in G_W}} w_g}}.$$

Weights were not rescaled, normalized, thresholded again, or exponentiated. Missing
genes are excluded before both numerator and denominator are formed. Regions are
ranked independently for each score; ranking is descriptive.

## Gene contributions

For the ten highest weighted-score regions, the additive contribution of gene $g$ is
$E_{{rg}}w_g / \sum_g w_g$. These contributions sum to the region's weighted score.
The heatmap additionally divides each additive contribution by the regional score,
so its displayed fractions sum to one within each region.

## Donor robustness

The six Stage 2 donor matrices were scored with the same matched genes and weights.
Unobserved donor-region combinations remain missing; they were not imputed.
Pairwise Pearson correlations use only regions observed in both donors. Mean
pairwise correlations were broad={mean_donor['broad']:.3f},
stringent={mean_donor['stringent']:.3f}, and weighted={mean_donor['weighted']:.3f}.

## Sensitivity analysis

Raw-scale versus gene-standardized regional Pearson correlations were
broad={sensitivity_r['broad']:.3f}, stringent={sensitivity_r['stringent']:.3f}, and
weighted={sensitivity_r['weighted']:.3f}. The long-form sensitivity file preserves
both profiles for every region and method.

## Limitations

The six post-mortem donors are few and unevenly sampled, four predominantly in the
left hemisphere. The pooled matrix retains only regions meeting the Stage 2 coverage
rule and is an adult healthy-brain reference. Exact-symbol matching leaves genes
without reliable AHBA representation out of the score. Broad anatomical-structure
normalization limits interpretation of absolute offsets between cortex,
subcortex/brainstem, and cerebellum. L2G is prioritisation evidence rather than a
definitive gene assignment. These observed rankings require an appropriate matched
null model before inferential interpretation.
"""
    (REPORTS / "stage_05_methods.md").write_text(content, encoding="utf-8")


def print_summary(
    expression: pd.DataFrame,
    gene_sets: dict[str, pd.DataFrame],
    matched: dict[str, list[str]],
    rankings: dict[str, pd.DataFrame],
    correlations: pd.DataFrame,
    donor_correlations: pd.DataFrame,
) -> None:
    print("\nPARKINSON × AHBA RAW SPATIAL SIGNAL")
    print(f"AHBA regions: {expression.shape[0]}")
    print(f"AHBA genes: {expression.shape[1]}")
    print(f"Broad Parkinson genes: {len(gene_sets['broad'])}")
    print(f"Broad genes represented: {len(matched['broad'])}")
    print(f"Stringent Parkinson genes: {len(gene_sets['stringent'])}")
    print(f"Stringent genes represented: {len(matched['stringent'])}")
    print(f"Weighted genes represented: {len(matched['weighted'])}")
    for method in ("broad", "stringent", "weighted"):
        names = ", ".join(rankings[method].head(10)["region"])
        print(f"Top 10 regions by {method} score: {names}")
    print("Correlation:")
    print(f"  broad vs stringent: {correlations.loc['broad', 'stringent']:.4f}")
    print(f"  broad vs weighted: {correlations.loc['broad', 'weighted']:.4f}")
    print(f"  stringent vs weighted: {correlations.loc['stringent', 'weighted']:.4f}")
    weighted = donor_correlations.loc[donor_correlations["score_method"] == "weighted", "pearson_r"]
    print(f"Donor concordance mean correlation (weighted): {weighted.mean():.4f}")
    print("STOP: observed spatial expression only; no null model or inferential analysis was performed.")


def main() -> None:
    for directory in (DATA_RESULTS, TABLES, FIGURES, BRAIN_MAPS, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    expression, regions, _, gene_sets, audit = load_and_validate()
    raw_scores, matched = calculate_scores(expression, gene_sets)
    standardized_scores, _ = calculate_scores(standardized_expression(expression), gene_sets)
    regional = add_region_names(raw_scores, regions)
    method_corr = method_correlations(raw_scores)
    sensitivity = sensitivity_table(raw_scores, standardized_scores, regions)
    contributions = gene_contributions(expression, raw_scores, regions, gene_sets["weighted"])
    donor, donor_corr = donor_scores(gene_sets, regions)

    rankings = {
        "broad": ranking_table(regional, "broad_mean_expression", len(matched["broad"])),
        "stringent": ranking_table(regional, "stringent_mean_expression", len(matched["stringent"])),
        "weighted": ranking_table(regional, "weighted_mean_expression", len(matched["weighted"])),
    }

    regional.to_csv(DATA_RESULTS / "parkinson_regional_raw_scores.csv", index=False)
    audit.to_csv(DATA_RESULTS / "parkinson_gene_match_audit.csv", index=False)
    contributions.to_csv(DATA_RESULTS / "parkinson_top_region_gene_contributions.csv", index=False)
    donor.to_csv(DATA_RESULTS / "parkinson_donor_regional_scores.csv", index=False)
    donor_corr.to_csv(DATA_RESULTS / "parkinson_donor_pattern_correlations.csv", index=False)
    method_corr.to_csv(DATA_RESULTS / "parkinson_score_method_correlations.csv", index_label="method")
    sensitivity.to_csv(DATA_RESULTS / "parkinson_expression_scale_sensitivity.csv", index=False)
    for method, table in rankings.items():
        table.to_csv(TABLES / f"parkinson_top_regions_{method}.csv", index=False)

    make_brain_maps(regional)
    make_figures(regional, method_corr, donor, donor_corr, contributions)
    write_methods(expression, regions, gene_sets, audit, donor_corr, sensitivity)
    print_summary(expression, gene_sets, matched, rankings, method_corr, donor_corr)


if __name__ == "__main__":
    main()
