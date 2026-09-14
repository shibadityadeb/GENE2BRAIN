"""Stage 10 common AHBA, null-model, spatial, and cross-disease analysis.

Only diseases passing the frozen source and gene-set QC gates enter inference.
The implementation reuses the Parkinson Stage 5/6 numerical functions and the
Stage 7 Moran spectral randomization primitives without disease-level tuning.
"""

from __future__ import annotations

from datetime import date
import json
import math
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import norm, pearsonr, spearmanr
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

import stage_05_parkinson_spatial_signal as reference_stage5
import stage_06_parkinson_enrichment as reference_stage6
import stage_07_spatial_robustness as reference_stage7


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
PROCESSED = ROOT / "data" / "processed"
GENES = ROOT / "data" / "genes" / "multidisease"
DATA_RESULTS = ROOT / "data" / "results"
VALIDATION = ROOT / "data" / "validation"
INTERMEDIATE = ROOT / "data" / "intermediate"
REPORTS = ROOT / "reports"
DISEASE_REPORTS = REPORTS / "diseases"
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"

METHODS = ("broad", "stringent", "weighted")
PRIMARY_METHOD = "weighted"
ACCESS_DATE = date.today().isoformat()


def load_json_yaml(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def panel() -> list[dict]:
    return load_json_yaml(CONFIG / "disease_panel.yaml")["diseases"]


def parameters() -> dict:
    return load_json_yaml(CONFIG / "statistical_parameters.yaml")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = pd.read_csv(PROCESSED / "brain_region_gene_expression.csv", index_col="region_id")
    expression.index = expression.index.astype(int)
    regions = pd.read_csv(PROCESSED / "region_metadata.csv").set_index("region_id")
    regions = regions.loc[regions["retained_in_main_matrix"].astype(bool)].reindex(expression.index)
    gene_metadata = pd.read_csv(PROCESSED / "gene_metadata.csv")
    errors = []
    if expression.shape != (138, 15_632):
        errors.append(f"frozen AHBA matrix shape changed: {expression.shape}")
    if expression.isna().any().any() or not np.isfinite(expression.to_numpy()).all():
        errors.append("AHBA matrix contains missing or non-finite values")
    if set(expression.columns) != set(gene_metadata["gene_symbol"]):
        errors.append("AHBA matrix and gene metadata differ")
    if regions["region_name"].isna().any():
        errors.append("region metadata does not cover the AHBA matrix")
    if errors:
        raise RuntimeError("Stage 10 frozen-input audit failed:\n- " + "\n- ".join(errors))
    return expression, regions, gene_metadata


def load_gene_sets(disease_id: str) -> dict[str, pd.DataFrame]:
    return {method: pd.read_csv(GENES / f"{disease_id}_{method}.csv") for method in METHODS}


def coverage_and_qc(expression: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    available = set(expression.columns)
    coverage_rows: list[dict] = []
    qc_rows: list[dict] = []
    analyzable: list[dict] = []
    limits = parameters()
    for disease in panel():
        if disease["status"] == "excluded":
            for method in METHODS:
                coverage_rows.append({
                    "disease": disease["disease_name"], "disease_id": disease["disease_id"],
                    "gene_set": method, "total_genes": 0, "genes_in_AHBA": 0,
                    "genes_missing": 0, "missing_gene_symbols": "", "coverage_percent": np.nan,
                    "analysis_status": "excluded_before_gene_prioritization",
                })
            qc_rows.append({
                "disease": disease["disease_name"], "configured_status": "excluded",
                "analysis_status": "excluded", "broad_genes": 0, "stringent_genes": 0,
                "weighted_genes": 0, "weighted_genes_in_ahba": 0,
                "weighted_coverage_percent": np.nan, "duplicate_gene_rows": 0,
                "missing_identifiers": 0, "l2g_min": np.nan, "l2g_median": np.nan,
                "l2g_max": np.nan, "gwas_sample_size": disease.get("sample_size"),
                "ancestry": disease.get("ancestry"), "flag": disease["selection_note"],
            })
            continue
        sets = load_gene_sets(disease["disease_id"])
        set_stats: dict[str, dict] = {}
        duplicates = 0
        missing_identifiers = 0
        for method, frame in sets.items():
            duplicates += int(frame["gene"].duplicated().sum())
            missing_identifiers += int(frame["gene"].isna().sum() + frame["ensembl_id"].isna().sum())
            genes = frame["gene"].dropna().astype(str).tolist()
            represented = [gene for gene in genes if gene in available]
            missing = [gene for gene in genes if gene not in available]
            pct = 100.0 * len(represented) / len(genes) if genes else np.nan
            set_stats[method] = {"total": len(genes), "present": len(represented), "coverage": pct}
            coverage_rows.append({
                "disease": disease["disease_name"], "disease_id": disease["disease_id"],
                "gene_set": method, "total_genes": len(genes), "genes_in_AHBA": len(represented),
                "genes_missing": len(missing), "missing_gene_symbols": ";".join(missing),
                "coverage_percent": pct, "analysis_status": disease["status"],
            })
        weighted = sets["weighted"]
        warning = (
            f"stringent sensitivity set has only {set_stats['stringent']['present']} AHBA genes; "
            "retained without threshold changes and interpreted cautiously"
            if set_stats["stringent"]["present"] < limits["minimum_ahba_genes"]
            else "none"
        )
        passes = (
            disease["status"] == "ready"
            and set_stats["weighted"]["present"] >= limits["minimum_ahba_genes"]
            and set_stats["broad"]["present"] >= limits["minimum_ahba_genes"]
            and all(set_stats[m]["coverage"] >= limits["minimum_ahba_coverage_percent"] for m in METHODS)
            and duplicates == 0 and missing_identifiers == 0
        )
        status = "analyzable" if passes else "needs_review_not_analyzed"
        if passes:
            analyzable.append(disease)
        qc_rows.append({
            "disease": disease["disease_name"], "configured_status": disease["status"],
            "analysis_status": status, "broad_genes": set_stats["broad"]["total"],
            "stringent_genes": set_stats["stringent"]["total"],
            "weighted_genes": set_stats["weighted"]["total"],
            "weighted_genes_in_ahba": set_stats["weighted"]["present"],
            "weighted_coverage_percent": set_stats["weighted"]["coverage"],
            "duplicate_gene_rows": duplicates, "missing_identifiers": missing_identifiers,
            "l2g_min": weighted["gene_weight"].min(), "l2g_median": weighted["gene_weight"].median(),
            "l2g_max": weighted["gene_weight"].max(), "gwas_sample_size": disease.get("sample_size"),
            "ancestry": disease.get("ancestry"),
            "flag": warning if passes else disease["selection_note"],
        })
    coverage = pd.DataFrame(coverage_rows)
    qc = pd.DataFrame(qc_rows)
    coverage.to_csv(DATA_RESULTS / "multidisease_ahba_gene_coverage.csv", index=False)
    qc.to_csv(DATA_RESULTS / "multidisease_qc_summary.csv", index=False)
    return coverage, qc, analyzable


def raw_scores(
    disease: dict, expression: pd.DataFrame, regions: pd.DataFrame, gene_sets: dict[str, pd.DataFrame]
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    represented = {
        method: [gene for gene in frame["gene"] if gene in expression.columns]
        for method, frame in gene_sets.items()
    }
    broad = expression[represented["broad"]].mean(axis=1)
    stringent = expression[represented["stringent"]].mean(axis=1)
    weights = gene_sets["weighted"].set_index("gene").loc[represented["weighted"], "gene_weight"].astype(float)
    weighted = expression[represented["weighted"]].mul(weights, axis=1).sum(axis=1) / weights.sum()
    wide = pd.DataFrame({
        "broad_mean_expression": broad,
        "stringent_mean_expression": stringent,
        "weighted_mean_expression": weighted,
    }, index=expression.index)
    long = []
    for method, column in {
        "broad": "broad_mean_expression", "stringent": "stringent_mean_expression",
        "weighted": "weighted_mean_expression",
    }.items():
        for region_id, score in wide[column].items():
            long.append({
                "disease": disease["disease_name"], "disease_id": disease["disease_id"],
                "gene_set": method, "region": regions.loc[region_id, "region_name"],
                "region_id": region_id, "score": score,
                "number_of_genes": len(represented[method]), "analysis_scale": "raw_scale_analysis",
            })
    return wide, pd.DataFrame(long), represented


def donor_scores(disease: dict, gene_sets: dict[str, pd.DataFrame], regions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores, correlations = reference_stage5.donor_scores(gene_sets, regions)
    scores.insert(0, "disease", disease["disease_name"])
    scores.insert(1, "disease_id", disease["disease_id"])
    correlations.insert(0, "disease", disease["disease_name"])
    correlations.insert(1, "disease_id", disease["disease_id"])
    return scores, correlations


def run_gene_set_null(
    disease: dict,
    panel_index: int,
    expression: pd.DataFrame,
    regions: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    gene_sets: dict[str, pd.DataFrame],
    wide_scores: pd.DataFrame,
    represented: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    params = parameters()
    reference_stage6.N_PERMUTATIONS = int(params["gene_set_permutations"])
    reference_stage6.NEIGHBOR_POOL_SIZE = int(params["gene_set_matching_neighbors"])
    reference_stage6.FDR_THRESHOLD = float(params["fdr_threshold"])
    covariates = reference_stage6.gene_covariates(expression, gene_metadata, represented["broad"])
    background = [gene for gene in expression.columns if gene not in set(represented["broad"])]
    background_values = expression[background].to_numpy(dtype=np.float32)
    broad_order, _ = reference_stage6.nearest_neighbors(represented["broad"], background, covariates)
    stringent_order, _ = reference_stage6.nearest_neighbors(represented["stringent"], background, covariates)
    rng = np.random.default_rng(int(params["gene_set_random_seed"]) + panel_index)
    broad_random = reference_stage6.sample_unique_sets(broad_order, len(background), rng)
    stringent_random = reference_stage6.sample_unique_sets(stringent_order, len(background), rng)
    weights = gene_sets["weighted"].set_index("gene").loc[represented["weighted"], "gene_weight"].to_numpy(np.float32)
    nulls = {
        "broad": reference_stage6.null_scores(background_values, broad_random),
        "stringent": reference_stage6.null_scores(background_values, stringent_random),
        "weighted": reference_stage6.null_scores(background_values, broad_random, weights=weights),
    }
    columns = {
        "broad": "broad_mean_expression", "stringent": "stringent_mean_expression",
        "weighted": "weighted_mean_expression",
    }
    outputs = []
    for method in METHODS:
        table = reference_stage6.enrichment_table(wide_scores[columns[method]], nulls[method], regions, method)
        table.insert(0, "disease", disease["disease_name"])
        table.insert(1, "disease_id", disease["disease_id"])
        table = table.rename(columns={"gene_set_version": "gene_set", "region_name": "region"})
        outputs.append(table)
    enrichment = pd.concat(outputs, ignore_index=True)
    robustness = reference_stage6.method_comparison({
        method: outputs[index].rename(columns={"gene_set": "gene_set_version", "region": "region_name"})
        for index, method in enumerate(METHODS)
    })
    robustness.insert(0, "disease", disease["disease_name"])
    robustness.insert(1, "disease_id", disease["disease_id"])
    return enrichment, robustness


def spatial_null(
    disease: dict,
    panel_index: int,
    enrichment: pd.DataFrame,
    donor: pd.DataFrame,
    regions: pd.DataFrame,
    weights: np.ndarray,
    eigenvectors: np.ndarray,
) -> pd.DataFrame:
    params = parameters()
    n_perm = int(params["spatial_permutations"])
    reference_stage7.N_PERMUTATIONS = n_perm
    fdr_threshold = float(params["fdr_threshold"])
    frames = []
    for method_index, method in enumerate(METHODS):
        subset = enrichment.loc[enrichment["gene_set"] == method].set_index("region_id").reindex(regions.index)
        values = subset["z_score"].to_numpy(float)
        rng = np.random.default_rng(int(params["spatial_random_seed"]) + 10 * panel_index + method_index)
        surrogates = reference_stage7.moran_singleton_randomization(values, eigenvectors, rng)
        spatial_p = reference_stage7.empirical_upper(values[None, :], surrogates, axis=0)
        spatial_fdr = multipletests(spatial_p, method="fdr_bh")[1]
        percentile = 1.0 - spatial_p
        frame = pd.DataFrame({
            "disease": disease["disease_name"], "disease_id": disease["disease_id"],
            "gene_set": method, "region_id": regions.index,
            "region": regions["region_name"].to_numpy(), "z_score": values,
            "gene_set_fdr": subset["fdr_p"].to_numpy(float),
            "spatial_null_p": spatial_p, "spatial_null_fdr": spatial_fdr,
            "spatial_robustness": percentile,
            "robustness_rank": pd.Series(percentile).rank(method="min", ascending=False).astype(int),
            "spatially_robust": (values > 0) & (subset["fdr_p"].to_numpy(float) < fdr_threshold) & (spatial_fdr < fdr_threshold),
            "n_spatial_permutations": n_perm,
        })
        frames.append(frame)
    output = pd.concat(frames, ignore_index=True)
    donor_weighted = donor.pivot(index="region_id", columns="donor", values="weighted_score").reindex(regions.index)
    donor_pct = donor_weighted.rank(axis=0, pct=True)
    donor_summary = pd.DataFrame({
        "region_id": regions.index.to_numpy(),
        "donor_median_rank_percentile": donor_pct.median(axis=1).to_numpy(),
        "donor_rank_iqr": (donor_pct.quantile(0.75, axis=1) - donor_pct.quantile(0.25, axis=1)).to_numpy(),
        "donors_observed": donor_weighted.notna().sum(axis=1).to_numpy(),
    })
    return output.merge(donor_summary, on="region_id", how="left")


def cross_disease_outputs(
    enrichment: pd.DataFrame, spatial: pd.DataFrame, robustness: pd.DataFrame,
    coverage: pd.DataFrame, qc: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    primary = enrichment.loc[enrichment["gene_set"] == PRIMARY_METHOD].copy()
    z_matrix = primary.pivot(index="disease", columns="region", values="z_score")
    fdr_matrix = primary.pivot(index="disease", columns="region", values="fdr_p").reindex(z_matrix.index)
    robust_primary = spatial.loc[spatial["gene_set"] == PRIMARY_METHOD]
    robust_matrix = robust_primary.pivot(index="disease", columns="region", values="spatially_robust").reindex(z_matrix.index).astype(int)
    z_matrix.to_csv(DATA_RESULTS / "disease_region_enrichment_matrix.csv")
    fdr_matrix.to_csv(DATA_RESULTS / "disease_region_fdr_matrix.csv")
    robust_matrix.to_csv(DATA_RESULTS / "disease_region_robustness_matrix.csv")

    similarity_rows = []
    for first_index, first in enumerate(z_matrix.index):
        for second in z_matrix.index[first_index:]:
            for metric, function in (("pearson", pearsonr), ("spearman", spearmanr)):
                result = function(z_matrix.loc[first], z_matrix.loc[second])
                similarity_rows.append({
                    "disease_1": first, "disease_2": second, "metric": metric,
                    "correlation": float(result.statistic), "p_value": float(result.pvalue),
                    "n_regions": z_matrix.shape[1],
                })
    similarity = pd.DataFrame(similarity_rows)
    similarity.to_csv(DATA_RESULTS / "disease_spatial_similarity.csv", index=False)

    shared = pd.DataFrame({
        "region": z_matrix.columns,
        "number_of_diseases_with_significant_enrichment": (fdr_matrix < parameters()["fdr_threshold"]).sum(axis=0).to_numpy(),
        "number_of_diseases_with_spatially_robust_enrichment": robust_matrix.sum(axis=0).to_numpy(),
        "mean_z_score": z_matrix.mean(axis=0).to_numpy(),
        "maximum_z_score": z_matrix.max(axis=0).to_numpy(),
    }).sort_values(
        ["number_of_diseases_with_spatially_robust_enrichment", "number_of_diseases_with_significant_enrichment", "mean_z_score"],
        ascending=False,
    )
    shared.to_csv(DATA_RESULTS / "shared_brain_region_enrichment.csv", index=False)

    specificity_rows = []
    for region in z_matrix.columns:
        values = z_matrix[region]
        for disease, value in values.items():
            others = values.drop(disease)
            sd = others.std(ddof=1)
            score = (value - others.mean()) / sd if sd > 0 else np.nan
            specificity_rows.append({
                "disease": disease, "region": region, "enrichment_z": value,
                "other_disease_mean_z": others.mean(), "other_disease_sd_z": sd,
                "cross_disease_specificity_z": score,
                "one_sided_p": norm.sf(score) if np.isfinite(score) else np.nan,
                "framework": "leave-one-disease-out standardization within region",
            })
    specificity = pd.DataFrame(specificity_rows)
    valid = specificity["one_sided_p"].notna()
    specificity["fdr_p"] = np.nan
    specificity.loc[valid, "fdr_p"] = multipletests(specificity.loc[valid, "one_sided_p"], method="fdr_bh")[1]
    specificity = specificity.sort_values(["fdr_p", "cross_disease_specificity_z"], ascending=[True, False])
    specificity.to_csv(DATA_RESULTS / "disease_specific_regional_signatures.csv", index=False)
    robustness.to_csv(DATA_RESULTS / "multidisease_gene_set_robustness.csv", index=False)

    metadata = pd.read_csv(DATA_RESULTS / "multidisease_gwas_metadata.csv")
    weighted_coverage = coverage.loc[coverage["gene_set"] == PRIMARY_METHOD].set_index("disease")
    qc_by_disease = qc.set_index("disease")
    power_rows = []
    for disease in z_matrix.index:
        values = primary.loc[primary["disease"] == disease]
        power_rows.append({
            "disease": disease,
            "gwas_sample_size": metadata.set_index("disease").loc[disease, "sample_size"],
            "gene_set_size": qc_by_disease.loc[disease, "weighted_genes_in_ahba"],
            "ahba_coverage_percent": weighted_coverage.loc[disease, "coverage_percent"],
            "number_of_enriched_regions": int((values["fdr_p"] < parameters()["fdr_threshold"]).sum()),
            "maximum_z_score": values["z_score"].max(),
            "number_of_spatially_robust_regions": int(robust_matrix.loc[disease].sum()),
        })
    power = pd.DataFrame(power_rows)
    power.to_csv(DATA_RESULTS / "disease_power_relationships.csv", index=False)
    return {
        "z": z_matrix, "fdr": fdr_matrix, "robust": robust_matrix,
        "similarity": similarity, "shared": shared, "specificity": specificity,
        "power": power,
    }


def validation_metadata() -> pd.DataFrame:
    rows = []
    pd_validation = pd.read_csv(DATA_RESULTS / "parkinson_independent_validation.csv")
    for disease in panel():
        if disease["disease_id"] == "parkinson":
            rows.append({
                "disease": disease["disease_name"],
                "validation_dataset": "Stage 8 curated substantia-nigra neuron-loss regional phenotype",
                "phenotype": "independent regional neuron-loss severity proxy",
                "sample_size": pd.NA, "region_count": len(pd_validation), "independent": True,
                "source": disease["validation_source"], "validation_type": "regional pathology proxy",
                "status": "partially_validated", "validation_outcome": "NOT SUPPORTED",
            })
        elif disease["status"] == "ready":
            rows.append({
                "disease": disease["disease_name"], "validation_dataset": "not selected",
                "phenotype": "not selected", "sample_size": pd.NA, "region_count": 0,
                "independent": pd.NA, "source": disease["validation_source"],
                "validation_type": "pending systematic source review", "status": "pending",
                "validation_outcome": "not analyzed",
            })
        else:
            rows.append({
                "disease": disease["disease_name"], "validation_dataset": "not available for current analysis",
                "phenotype": "not analyzed", "sample_size": pd.NA, "region_count": 0,
                "independent": pd.NA, "source": disease["validation_source"],
                "validation_type": "not assessed after genetic/QC exclusion", "status": "not_available",
                "validation_outcome": "not analyzed",
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(VALIDATION / "multidisease_validation_metadata.csv", index=False)
    return frame


def clustered_order(z: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    row_distance = 1 - z.T.corr(method="pearson").clip(-1, 1)
    column_distance = 1 - z.corr(method="pearson").clip(-1, 1)
    np.fill_diagonal(row_distance.values, 0)
    np.fill_diagonal(column_distance.values, 0)
    row_link = linkage(squareform(row_distance.values, checks=False), method="average")
    col_link = linkage(squareform(column_distance.values, checks=False), method="average")
    return leaves_list(row_link), leaves_list(col_link), row_link


def make_figures(outputs: dict[str, pd.DataFrame], robustness: pd.DataFrame) -> None:
    sns.set_theme(style="white", context="notebook")
    z = outputs["z"]
    row_order, col_order, row_link = clustered_order(z)
    ordered = z.iloc[row_order, col_order]
    limit = float(np.nanquantile(np.abs(z.to_numpy()), 0.99))
    for filename, title in (
        ("stage_10_disease_region_atlas.png", "GENE2BRAIN Disease × Brain Region Atlas"),
        ("stage_10_multidisease_signature_matrix.png", "Multi-Disease Spatial Signature Matrix"),
    ):
        fig, ax = plt.subplots(figsize=(24, 8))
        sns.heatmap(ordered, cmap="vlag", center=0, vmin=-limit, vmax=limit, ax=ax, cbar_kws={"label": "Matched gene-set Z score"})
        ax.set(title=title, xlabel="AAL3 region (average-linkage clustered)", ylabel="Disease (average-linkage clustered)")
        ax.tick_params(axis="x", labelsize=5)
        fig.tight_layout()
        fig.savefig(FIGURES / filename, dpi=220)
        plt.close(fig)

    pearson = outputs["similarity"].loc[outputs["similarity"]["metric"] == "pearson"]
    matrix = pd.DataFrame(np.eye(len(z)), index=z.index, columns=z.index)
    for row in pearson.itertuples():
        matrix.loc[row.disease_1, row.disease_2] = row.correlation
        matrix.loc[row.disease_2, row.disease_1] = row.correlation
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(matrix, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Disease Spatial Similarity\nPearson correlation across 138 regional enrichment Z scores")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_disease_spatial_similarity.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 6))
    dendrogram(row_link, labels=z.index.tolist(), leaf_rotation=35, leaf_font_size=9, ax=ax)
    ax.set(title="Disease Spatial Clustering", ylabel="Correlation distance", xlabel="No biological cluster labels assigned")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_disease_clustering.png", dpi=220)
    plt.close(fig)

    pca = PCA(n_components=2).fit_transform(z)
    fig, ax = plt.subplots(figsize=(10, 8))
    categories = {x["disease_name"]: x["category"] for x in panel()}
    palette = dict(zip(sorted(set(categories.values())), sns.color_palette("colorblind", len(set(categories.values())))))
    for i, disease in enumerate(z.index):
        ax.scatter(pca[i, 0], pca[i, 1], s=80, color=palette[categories[disease]], label=categories[disease])
        ax.annotate(disease, pca[i], xytext=(5, 4), textcoords="offset points", fontsize=8)
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), title="Predefined category", fontsize=8)
    ax.axhline(0, color="0.85", lw=1); ax.axvline(0, color="0.85", lw=1)
    ax.set(title="GENE2BRAIN Disease Landscape", xlabel="PC1", ylabel="PC2")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_disease_pca.png", dpi=220)
    plt.close(fig)

    top_shared = outputs["shared"].head(20).sort_values("number_of_diseases_with_significant_enrichment")
    top_specific = outputs["specificity"].head(20).sort_values("cross_disease_specificity_z")
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    axes[0].barh(top_shared["region"], top_shared["number_of_diseases_with_significant_enrichment"], color="#426b69")
    axes[0].set(title="Shared signal", xlabel="Diseases with weighted-set FDR < 0.05")
    labels = top_specific["disease"].str.slice(0, 18) + " · " + top_specific["region"]
    axes[1].barh(labels, top_specific["cross_disease_specificity_z"], color="#b76145")
    axes[1].set(title="Disease-specific signal", xlabel="Leave-one-disease-out specificity Z")
    fig.suptitle("Shared vs Disease-Specific Brain Signatures")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_shared_vs_specific.png", dpi=220)
    plt.close(fig)

    threshold = float(np.quantile(np.abs(matrix.to_numpy()[np.triu_indices(len(matrix), 1)]), 0.90))
    angles = np.linspace(0, 2 * np.pi, len(matrix), endpoint=False)
    positions = np.column_stack([np.cos(angles), np.sin(angles)])
    fig, ax = plt.subplots(figsize=(10, 10))
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            value = matrix.iloc[i, j]
            if abs(value) >= threshold:
                ax.plot(positions[[i, j], 0], positions[[i, j], 1], color="#32677a" if value > 0 else "#a84d45", alpha=0.25 + 0.65 * abs(value), lw=1 + 3 * abs(value), zorder=1)
    category_palette = {key: value for key, value in zip(sorted(set(categories.values())), sns.color_palette("colorblind", len(set(categories.values()))))}
    sizes = outputs["power"].set_index("disease")["gwas_sample_size"].reindex(matrix.index)
    sizes = 120 + 500 * np.sqrt(sizes / sizes.max())
    for i, disease in enumerate(matrix.index):
        ax.scatter(*positions[i], s=sizes.iloc[i], color=category_palette[categories[disease]], edgecolor="white", zorder=2)
        ax.text(*(positions[i] * 1.15), disease, ha="center", va="center", fontsize=8)
    ax.set(title=f"Disease Similarity Network\n|Pearson r| ≥ empirical 90th percentile ({threshold:.2f}); size = GWAS sample size", xlim=(-1.45, 1.45), ylim=(-1.45, 1.45), aspect="equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_disease_similarity_network.png", dpi=220)
    plt.close(fig)

    plot = robustness.loc[robustness["method_1"] != robustness["method_2"]].copy()
    plot["comparison"] = plot["method_1"] + " vs " + plot["method_2"]
    fig, ax = plt.subplots(figsize=(12, 7))
    sns.boxplot(data=plot, x="comparison", y="spearman_rho", hue="comparison", legend=False, ax=ax)
    sns.stripplot(data=plot, x="comparison", y="spearman_rho", color="black", alpha=.6, ax=ax)
    ax.set(title="Robustness Across Gene-Set Definitions", ylabel="Regional Spearman correlation", xlabel="")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_gene_set_robustness.png", dpi=220)
    plt.close(fig)

    power = outputs["power"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for disease, row in power.set_index("disease").iterrows():
        axes[0].scatter(row.gwas_sample_size, row.maximum_z_score, s=60)
        axes[0].annotate(disease, (row.gwas_sample_size, row.maximum_z_score), fontsize=7)
        axes[1].scatter(row.gene_set_size, row.number_of_enriched_regions, s=60)
        axes[1].annotate(disease, (row.gene_set_size, row.number_of_enriched_regions), fontsize=7)
    axes[0].set_xscale("log")
    axes[0].set(xlabel="GWAS sample size (log)", ylabel="Maximum regional Z", title="GWAS power vs enrichment peak")
    axes[1].set(xlabel="AHBA-represented weighted genes", ylabel="FDR-significant regions", title="Gene-set size vs discoveries")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_gwas_power_vs_enrichment.png", dpi=220)
    plt.close(fig)

    steps = ["GWAS", "Loci", "Genes", "Healthy brain", "Regional expression", "Enrichment", "Spatial robustness", "Validation", "Biology", "Cross-disease atlas"]
    fig, ax = plt.subplots(figsize=(18, 3))
    ax.axis("off")
    for index, step in enumerate(steps):
        x = index / (len(steps) - 1)
        ax.text(x, .5, step, ha="center", va="center", bbox=dict(boxstyle="round,pad=.45", fc="#e8f0ef", ec="#426b69"), transform=ax.transAxes, fontsize=9)
        if index < len(steps) - 1:
            ax.annotate("", xy=((index + .82) / (len(steps)-1), .5), xytext=((index + .18) / (len(steps)-1), .5), xycoords=ax.transAxes, arrowprops=dict(arrowstyle="->", color="#777"))
    ax.set_title("GENE2BRAIN: From Genetic Risk to Brain", fontsize=16, pad=18)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_gene2brain_pipeline.png", dpi=220)
    plt.close(fig)


def write_reports(
    coverage: pd.DataFrame, qc: pd.DataFrame, enrichment: pd.DataFrame,
    spatial: pd.DataFrame, outputs: dict[str, pd.DataFrame], validation: pd.DataFrame,
) -> None:
    analyzed = qc.loc[qc["analysis_status"] == "analyzable", "disease"].tolist()
    excluded = qc.loc[qc["analysis_status"] != "analyzable", "disease"].tolist()
    qc_text = f"""# Stage 10 multi-disease quality control

This report was generated before spatial inference using the frozen panel and
thresholds. {len(analyzed)} diseases passed and {len(excluded)} did not.

## Analyzed

{chr(10).join(f'- {name}' for name in analyzed)}

## Not analyzed

{chr(10).join(f'- {row.disease}: {row.flag}' for row in qc.loc[qc.analysis_status != 'analyzable'].itertuples())}

Checks include gene-set size, exact AHBA symbol overlap, duplicate identifiers,
missing Ensembl identifiers, L2G distribution, ancestry, GWAS sample size, and
the frozen minimum-five-gene / 50%-coverage gate. Full values are in
`data/results/multidisease_qc_summary.csv` and
`data/results/multidisease_ahba_gene_coverage.csv`.
"""
    (REPORTS / "stage_10_multidisease_qc.md").write_text(qc_text, encoding="utf-8")

    audit = f"""# Stage 10 method consistency audit

All {len(analyzed)} analyzed diseases used the same frozen AHBA matrix (138 AAL3
regions × 15,632 genes), region definitions, raw-scale Stage 5 equations,
exact-symbol matching, Stage 6 matching covariates (mean expression, log10
expression variance, log1p reannotated probe count), 200-nearest-neighbor unique
sampling, 10,000 gene-set permutations, one-sided +1 empirical p-values,
Benjamini-Hochberg correction, binary 26-neighbor AAL3 graph, Moran singleton
spectral randomization, 10,000 spatial permutations, and FDR 0.05.

The weighted matched null reuses broad matched replacements and carries the
fixed disease L2G weights, exactly as for Parkinson. Cross-disease comparisons
use the primary matched-null regional Z score; raw expression is never compared
between diseases.

## Exceptions

- Parkinson uses its frozen Stage 3–9 source artifacts without modification.
- Major depression uses the same-publication European component because the
  larger bi-ancestry Catalog record has zero Open Targets credible sets.
- Autism and progressive supranuclear palsy have real L2G outputs but fail the
  prespecified five-gene QC gate and are not spatially analyzed.
- Epilepsy passes the primary weighted-set gate with seven AHBA genes, but its
  two-gene stringent sensitivity result is explicitly flagged as unstable.
- Huntington, multiple system atrophy, and Tourette syndrome fail genetic/L2G
  inclusion and remain documented exclusions.
- Independent validation type is allowed to differ by disease. Only Parkinson
  has a completed Stage 8 analysis; its outcome remains NOT SUPPORTED. Other
  analyzed diseases remain pending rather than receiving proxy data.
"""
    (REPORTS / "stage_10_method_consistency_audit.md").write_text(audit, encoding="utf-8")

    for disease in panel():
        name = disease["disease_name"]
        q = qc.loc[qc["disease"] == name].iloc[0]
        cov = coverage.loc[(coverage["disease"] == name) & (coverage["gene_set"] == PRIMARY_METHOD)].iloc[0]
        val = validation.loc[validation["disease"] == name].iloc[0]
        if q["analysis_status"] == "analyzable":
            regional = enrichment.loc[(enrichment["disease"] == name) & (enrichment["gene_set"] == PRIMARY_METHOD)]
            top = regional.nlargest(1, "z_score").iloc[0]
            robust_n = int(spatial.loc[(spatial["disease"] == name) & (spatial["gene_set"] == PRIMARY_METHOD), "spatially_robust"].sum())
            result_text = f"Top region: {top.region} (Z={top.z_score:.3f}, FDR={top.fdr_p:.4g}). Significant regions: {(regional.fdr_p < parameters()['fdr_threshold']).sum()}; spatially robust regions: {robust_n}."
        else:
            result_text = "Regional enrichment and spatial robustness were not run because the frozen inclusion/QC gate failed."
        report = f"""# {name}

## GWAS

- Accession: {disease.get('gwas_accession') or 'none'}
- Study: {disease.get('primary_gwas_study') or 'none'}
- Sample size: {disease.get('sample_size') or 'not applicable'}
- Ancestry: {disease.get('ancestry') or 'not applicable'}
- Panel status: {disease['status']}

## Gene prioritization and AHBA coverage

Open Targets L2G was used where compatible. Weighted prioritized genes: {q.weighted_genes}; represented in AHBA: {q.weighted_genes_in_ahba}; coverage: {cov.coverage_percent if pd.notna(cov.coverage_percent) else 'not applicable'}.

## Regional enrichment and spatial robustness

{result_text}

## Independent validation

Status: {val.status}. Outcome: {val.validation_outcome}.

## Biological interpretation

Stage 10 standardized pathway, cell-type, and driver-gene tables are populated
only for diseases that pass the spatial-analysis gate. Empty/significance-null
results are preserved.

## Limitations

{disease['selection_note']} L2G scores prioritize candidate genes but do not
prove causality. AHBA represents healthy adult donor expression and is neither a
patient map nor a diagnostic or predictive model.
"""
        (DISEASE_REPORTS / f"{disease['disease_id']}.md").write_text(report, encoding="utf-8")


def master_table(
    coverage: pd.DataFrame, qc: pd.DataFrame, enrichment: pd.DataFrame,
    spatial: pd.DataFrame, validation: pd.DataFrame,
) -> pd.DataFrame:
    metadata = pd.read_csv(DATA_RESULTS / "multidisease_gwas_metadata.csv").set_index("disease")
    rows = []
    for disease in panel():
        name = disease["disease_name"]
        q = qc.set_index("disease").loc[name]
        cov = coverage.loc[(coverage.disease == name) & (coverage.gene_set == PRIMARY_METHOD)].iloc[0]
        subset = enrichment.loc[(enrichment.disease == name) & (enrichment.gene_set == PRIMARY_METHOD)]
        if subset.empty:
            top_region = pd.NA; top_z = np.nan; top_fdr = np.nan; significant = 0; robust = 0
        else:
            top = subset.nlargest(1, "z_score").iloc[0]
            top_region, top_z, top_fdr = top.region, top.z_score, top.fdr_p
            significant = int((subset.fdr_p < parameters()["fdr_threshold"]).sum())
            robust = int(spatial.loc[(spatial.disease == name) & (spatial.gene_set == PRIMARY_METHOD), "spatially_robust"].sum())
        rows.append({
            "disease": name, "category": disease["category"],
            "gwas_sample_size": disease.get("sample_size"),
            "n_gwas_loci": metadata.loc[name, "number_of_loci"] if name in metadata.index else 0,
            "n_prioritized_genes": q.weighted_genes,
            "ahba_gene_coverage": cov.coverage_percent,
            "n_significant_regions": significant, "n_spatially_robust_regions": robust,
            "validation_status": validation.set_index("disease").loc[name, "status"],
            "top_region": top_region, "top_region_z": top_z, "top_region_fdr": top_fdr,
            "top_biological_process": "pending Stage 10 standardized biology table" if q.analysis_status == "analyzable" else pd.NA,
            "top_cell_type": "pending Stage 10 standardized biology table" if q.analysis_status == "analyzable" else pd.NA,
            "notes": q.flag,
        })
    output = pd.DataFrame(rows)
    output.to_csv(TABLES / "gene2brain_master_disease_results.csv", index=False)
    return output


def print_summary(master: pd.DataFrame, outputs: dict[str, pd.DataFrame]) -> None:
    analyzed = master.loc[master["top_region"].notna()]
    excluded = master.loc[master["top_region"].isna()]
    pairwise = outputs["similarity"].loc[(outputs["similarity"].metric == "pearson") & (outputs["similarity"].disease_1 != outputs["similarity"].disease_2)]
    unique_pairs = pairwise.assign(pair=pairwise.apply(lambda x: "|".join(sorted((x.disease_1, x.disease_2))), axis=1)).drop_duplicates("pair")
    most_similar = unique_pairs.nlargest(1, "correlation").iloc[0]
    most_distinct = unique_pairs.nsmallest(1, "correlation").iloc[0]
    recurrent = outputs["shared"].iloc[0]
    print("GENE2BRAIN MULTI-DISEASE SUMMARY")
    print("Diseases analyzed:", "; ".join(analyzed.disease))
    print("Diseases excluded/not analyzed:", "; ".join(excluded.disease))
    print("Total brain regions: 138")
    print("Total AHBA genes: 15,632")
    for row in master.itertuples():
        print(f"{row.disease}: sample={row.gwas_sample_size}; loci={row.n_gwas_loci}; genes={row.n_prioritized_genes}; coverage={row.ahba_gene_coverage}; significant={row.n_significant_regions}; robust={row.n_spatially_robust_regions}; validation={row.validation_status}; top={row.top_region}")
    print("Most region-specific disease:", outputs["specificity"].iloc[0]["disease"])
    print("Most spatially widespread disease:", analyzed.sort_values("n_significant_regions", ascending=False).iloc[0].disease)
    print(f"Most similar disease pair: {most_similar.disease_1} / {most_similar.disease_2} (r={most_similar.correlation:.3f})")
    print(f"Most distinct disease pair: {most_distinct.disease_1} / {most_distinct.disease_2} (r={most_distinct.correlation:.3f})")
    print("Most recurrent brain region:", recurrent.region)
    print(f"Mean robust regions: {analyzed.n_spatially_robust_regions.mean():.2f}")
    print("Interpret with GWAS power and gene-set size; no causality or patient prediction is claimed.")


def main() -> None:
    for directory in (DATA_RESULTS, VALIDATION, INTERMEDIATE, REPORTS, DISEASE_REPORTS, FIGURES, TABLES):
        directory.mkdir(parents=True, exist_ok=True)
    expression, regions, gene_metadata = load_inputs()
    coverage, qc, analyzable = coverage_and_qc(expression)
    weights = pd.read_csv(INTERMEDIATE / "brain_region_spatial_weights.csv", index_col="region_id").reindex(index=regions.index, columns=regions.index.astype(str)).to_numpy(float)
    eigenvectors, _ = reference_stage7.moran_eigenvectors(weights)
    raw_frames = []
    enrichment_frames = []
    robustness_frames = []
    spatial_frames = []
    donor_frames = []
    donor_correlation_frames = []
    for panel_index, disease in enumerate(analyzable):
        print(f"Analyzing {disease['disease_name']} ({panel_index + 1}/{len(analyzable)})", flush=True)
        gene_sets = load_gene_sets(disease["disease_id"])
        wide, raw, represented = raw_scores(disease, expression, regions, gene_sets)
        donor, donor_corr = donor_scores(disease, gene_sets, regions)
        enrichment, robustness = run_gene_set_null(
            disease, panel_index, expression, regions, gene_metadata, gene_sets, wide, represented
        )
        spatial = spatial_null(disease, panel_index, enrichment, donor, regions, weights, eigenvectors)
        raw_frames.append(raw); donor_frames.append(donor); donor_correlation_frames.append(donor_corr)
        enrichment_frames.append(enrichment); robustness_frames.append(robustness); spatial_frames.append(spatial)
    raw_all = pd.concat(raw_frames, ignore_index=True)
    enrichment_all = pd.concat(enrichment_frames, ignore_index=True)
    robustness_all = pd.concat(robustness_frames, ignore_index=True)
    spatial_all = pd.concat(spatial_frames, ignore_index=True)
    donor_all = pd.concat(donor_frames, ignore_index=True)
    donor_correlations = pd.concat(donor_correlation_frames, ignore_index=True)
    raw_all.to_csv(DATA_RESULTS / "multidisease_regional_raw_scores.csv", index=False)
    enrichment_all.to_csv(DATA_RESULTS / "multidisease_regional_enrichment.csv", index=False)
    spatial_all.to_csv(DATA_RESULTS / "multidisease_spatial_robustness.csv", index=False)
    donor_all.to_csv(DATA_RESULTS / "multidisease_donor_regional_scores.csv", index=False)
    donor_correlations.to_csv(DATA_RESULTS / "multidisease_donor_correlations.csv", index=False)
    outputs = cross_disease_outputs(enrichment_all, spatial_all, robustness_all, coverage, qc)
    validation = validation_metadata()
    make_figures(outputs, robustness_all)
    write_reports(coverage, qc, enrichment_all, spatial_all, outputs, validation)
    master = master_table(coverage, qc, enrichment_all, spatial_all, validation)
    print_summary(master, outputs)


if __name__ == "__main__":
    main()
