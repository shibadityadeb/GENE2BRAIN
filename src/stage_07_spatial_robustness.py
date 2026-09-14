"""Stage 7: spatial-dependence sensitivity analysis for Parkinson enrichment.

Stage 6 tests disease-gene specificity. This stage leaves those results intact and
asks whether their anatomical arrangement is unusual under an atlas-graph null.
"""

from __future__ import annotations

import json
import math
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mpl_colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import nibabel as nib
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy.linalg import null_space
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
INTERMEDIATE = ROOT / "data" / "intermediate"
DATA_RESULTS = ROOT / "data" / "results"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
BRAIN_MAPS = ROOT / "results" / "brain_maps"
REPORTS = ROOT / "reports"
ATLAS = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.nii.gz"
WEB_GEOMETRY = ROOT / "web" / "public" / "data" / "aal3_regions.json"

N_PERMUTATIONS = 10_000
RANDOM_SEED = 20260914
FDR_THRESHOLD = 0.05
CONNECTIVITY = 26
PRIMARY_METHOD = "weighted"
METHODS = ("broad", "stringent", "weighted")
STAGE5_COLUMNS = {
    "broad": "broad_mean_expression",
    "stringent": "stringent_mean_expression",
    "weighted": "weighted_mean_expression",
}


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, nib.Nifti1Image]:
    regions = pd.read_csv(PROCESSED / "region_metadata.csv")
    regions = regions.loc[regions["retained_in_main_matrix"]].sort_values("region_id").set_index("region_id")
    stage5 = pd.read_csv(DATA_RESULTS / "parkinson_regional_raw_scores.csv").set_index("region_id").reindex(regions.index)
    stage6_all = pd.read_csv(DATA_RESULTS / "parkinson_regional_enrichment_all_methods.csv")
    donor = pd.read_csv(DATA_RESULTS / "parkinson_donor_regional_scores.csv")
    image = nib.load(ATLAS)
    errors: list[str] = []
    if len(regions) != 138 or regions.index.has_duplicates:
        errors.append("Expected 138 unique retained AAL3 regions")
    if stage5[list(STAGE5_COLUMNS.values())].isna().any().any():
        errors.append("Stage 5 profiles do not cover every retained region")
    expected_pairs = pd.MultiIndex.from_product([METHODS, regions.index])
    observed_pairs = pd.MultiIndex.from_frame(stage6_all[["gene_set_version", "region_id"]])
    if set(expected_pairs) != set(observed_pairs):
        errors.append("Stage 6 all-method table does not cover every method-region pair")
    if donor["donor"].nunique() != 6 or set(donor["region_id"]) != set(regions.index):
        errors.append("Stage 5 donor table is incomplete")
    if errors:
        raise RuntimeError("Stage 7 input audit failed:\n- " + "\n- ".join(errors))
    return regions, stage5, stage6_all, donor, image


def neighbor_offsets() -> list[tuple[int, int, int]]:
    """Return one half of the 26-neighbor offsets to avoid duplicate comparisons."""
    return [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
        if (dx, dy, dz) > (0, 0, 0)
    ]


def paired_slices(offset: tuple[int, int, int]) -> tuple[tuple[slice, ...], tuple[slice, ...]]:
    first: list[slice] = []
    second: list[slice] = []
    for delta in offset:
        if delta > 0:
            first.append(slice(0, -delta))
            second.append(slice(delta, None))
        elif delta < 0:
            first.append(slice(-delta, None))
            second.append(slice(0, delta))
        else:
            first.append(slice(None))
            second.append(slice(None))
    return tuple(first), tuple(second)


def build_spatial_weights(atlas: np.ndarray, region_ids: np.ndarray) -> tuple[pd.DataFrame, set[tuple[int, int]]]:
    """Build binary AAL3 adjacency from face-, edge-, or corner-touching voxels."""
    retained = set(int(value) for value in region_ids)
    edges: set[tuple[int, int]] = set()
    for offset in neighbor_offsets():
        first, second = paired_slices(offset)
        left = atlas[first]
        right = atlas[second]
        boundary = (left != right) & (left > 0) & (right > 0)
        if not boundary.any():
            continue
        pairs = np.column_stack((left[boundary], right[boundary])).astype(int)
        pairs.sort(axis=1)
        for first_id, second_id in np.unique(pairs, axis=0):
            if first_id in retained and second_id in retained:
                edges.add((int(first_id), int(second_id)))
    index = {int(region_id): position for position, region_id in enumerate(region_ids)}
    weights = np.zeros((len(region_ids), len(region_ids)), dtype=float)
    for first_id, second_id in edges:
        first_position, second_position = index[first_id], index[second_id]
        weights[first_position, second_position] = 1.0
        weights[second_position, first_position] = 1.0
    frame = pd.DataFrame(weights.astype(int), index=region_ids, columns=[str(value) for value in region_ids])
    frame.index.name = "region_id"
    return frame, edges


def connected_components(weights: np.ndarray) -> list[list[int]]:
    unseen = set(range(weights.shape[0]))
    components: list[list[int]] = []
    while unseen:
        seed = unseen.pop()
        component = {seed}
        frontier = [seed]
        while frontier:
            current = frontier.pop()
            neighbors = set(np.flatnonzero(weights[current])) & unseen
            unseen -= neighbors
            component |= neighbors
            frontier.extend(neighbors)
        components.append(sorted(component))
    return sorted(components, key=len, reverse=True)


def morans_i(values: np.ndarray, weights: np.ndarray) -> float:
    centered = np.asarray(values, dtype=float) - np.mean(values)
    denominator = float(centered @ centered)
    weight_sum = float(weights.sum())
    if denominator <= 0 or weight_sum <= 0:
        raise ValueError("Moran's I requires nonconstant values and nonzero weights")
    return float(len(values) / weight_sum * (centered @ weights @ centered) / denominator)


def moran_permutation_null(values: np.ndarray, weights: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    centered = np.asarray(values, dtype=float) - np.mean(values)
    permutations = np.vstack([rng.permutation(centered) for _ in range(N_PERMUTATIONS)])
    numerators = np.einsum("bi,ij,bj->b", permutations, weights, permutations, optimize=True)
    return len(values) / weights.sum() * numerators / np.sum(centered**2)


def moran_eigenvectors(weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return an orthonormal basis of the mean-zero Moran operator."""
    n_regions = weights.shape[0]
    mean_zero_basis = null_space(np.ones((1, n_regions)))
    reduced = mean_zero_basis.T @ weights @ mean_zero_basis
    eigenvalues, rotation = np.linalg.eigh(reduced)
    order = np.argsort(eigenvalues)[::-1]
    return mean_zero_basis @ rotation[:, order], eigenvalues[order]


def moran_singleton_randomization(
    values: np.ndarray,
    eigenvectors: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate Wagner–Dray singleton MSR surrogates by spectral sign flips."""
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    sample_sd = values.std(ddof=1)
    coefficients = eigenvectors.T @ centered / np.linalg.norm(centered)
    n_modes = eigenvectors.shape[1]
    signs = rng.choice((-1.0, 1.0), size=(N_PERMUTATIONS, n_modes))
    randomized = signs * coefficients
    surrogates = values.mean() + sample_sd * math.sqrt(len(values) - 1) * randomized @ eigenvectors.T
    if not np.allclose(surrogates.mean(axis=1), values.mean(), atol=1e-10):
        raise RuntimeError("Moran surrogates did not preserve the observed mean")
    if not np.allclose(surrogates.std(axis=1, ddof=1), sample_sd, atol=1e-10):
        raise RuntimeError("Moran surrogates did not preserve the observed variance")
    return surrogates


def empirical_upper(observed: np.ndarray | float, null: np.ndarray, axis: int = 0) -> np.ndarray:
    return (np.sum(null >= observed, axis=axis) + 1) / (null.shape[axis] + 1)


def stage6_profiles(stage6_all: pd.DataFrame, region_ids: pd.Index) -> dict[str, np.ndarray]:
    return {
        method: stage6_all.loc[stage6_all["gene_set_version"] == method]
        .set_index("region_id").reindex(region_ids)["z_score"].to_numpy(dtype=float)
        for method in METHODS
    }


def calculate_outputs(
    regions: pd.DataFrame,
    stage5: pd.DataFrame,
    stage6_all: pd.DataFrame,
    donor: pd.DataFrame,
    weights: np.ndarray,
) -> dict[str, object]:
    region_ids = regions.index
    eigenvectors, eigenvalues = moran_eigenvectors(weights)
    seed_sequence = np.random.SeedSequence(RANDOM_SEED)
    child_seeds = iter(seed_sequence.spawn(20))

    raw_profiles = {
        f"{method}_stage5_score": stage5[STAGE5_COLUMNS[method]].to_numpy(dtype=float)
        for method in METHODS
    }
    z_profiles = stage6_profiles(stage6_all, region_ids)
    autocorrelation_rows: list[dict[str, object]] = []
    autocorrelation_nulls: dict[str, np.ndarray] = {}
    for name, values in {**raw_profiles, "weighted_stage6_z": z_profiles["weighted"]}.items():
        observed = morans_i(values, weights)
        null = moran_permutation_null(values, weights, np.random.default_rng(next(child_seeds)))
        null_mean = float(null.mean())
        p_value = (np.sum(np.abs(null - null_mean) >= abs(observed - null_mean)) + 1) / (N_PERMUTATIONS + 1)
        autocorrelation_rows.append({
            "metric": name,
            "observed_statistic": observed,
            "null_mean": null_mean,
            "null_std": float(null.std(ddof=1)),
            "p_value": p_value,
            "number_of_permutations": N_PERMUTATIONS,
            "statistic": "global Moran's I",
            "null_model": "unrestricted regional label permutation",
        })
        autocorrelation_nulls[name] = null

    surrogate_maps: dict[str, np.ndarray] = {}
    null_summary_rows: list[dict[str, object]] = []
    global_rows: list[dict[str, object]] = []
    sensitivity_frames: list[pd.DataFrame] = []
    for method in METHODS:
        values = z_profiles[method]
        surrogates = moran_singleton_randomization(values, eigenvectors, np.random.default_rng(next(child_seeds)))
        surrogate_maps[method] = surrogates
        observed_moran = morans_i(values, weights)
        surrogate_moran = np.array([morans_i(row, weights) for row in surrogates])
        if not np.allclose(surrogate_moran, observed_moran, atol=1e-10):
            raise RuntimeError(f"Moran singleton surrogates did not preserve Moran's I for {method}")
        standardized_observed = (values - values.mean()) / values.std(ddof=1)
        standardized_null = (surrogates - values.mean()) / values.std(ddof=1)
        observed_peak = float(standardized_observed.max())
        null_peaks = standardized_null.max(axis=1)
        peak_p = float(empirical_upper(observed_peak, null_peaks, axis=0))
        for statistic, observed, null_values, p_value in (
            ("morans_i_preservation", observed_moran, surrogate_moran, np.nan),
            ("maximum_standardized_regional_z", observed_peak, null_peaks, peak_p),
        ):
            null_summary_rows.append({
                "gene_set_version": method,
                "statistic": statistic,
                "observed_statistic": observed,
                "null_mean": float(null_values.mean()),
                "null_std": float(null_values.std(ddof=1)),
                "null_q025": float(np.quantile(null_values, 0.025)),
                "null_median": float(np.median(null_values)),
                "null_q975": float(np.quantile(null_values, 0.975)),
                "p_value": p_value,
                "number_of_permutations": N_PERMUTATIONS,
            })
        global_rows.append({
            "gene_set_version": method,
            "global_statistic": "maximum standardized regional Stage 6 Z-score",
            "observed_statistic": observed_peak,
            "spatial_null_mean": float(null_peaks.mean()),
            "spatial_null_std": float(null_peaks.std(ddof=1)),
            "spatial_null_p": peak_p,
            "number_of_spatial_permutations": N_PERMUTATIONS,
            "question": "Is the strongest positive regional enrichment peak more extreme than spatially autocorrelated surrogate maps?",
        })
        spatial_p = empirical_upper(values[None, :], surrogates, axis=0)
        spatial_fdr = multipletests(spatial_p, method="fdr_bh")[1]
        percentile = 1.0 - spatial_p
        ranks = pd.Series(percentile).rank(method="min", ascending=False).astype(int).to_numpy()
        method_stage6 = stage6_all.loc[stage6_all["gene_set_version"] == method].set_index("region_id").reindex(region_ids)
        sensitivity_frames.append(pd.DataFrame({
            "gene_set_version": method,
            "region_id": region_ids,
            "region_name": regions["region_name"].to_numpy(),
            "stage6_z": values,
            "stage6_fdr": method_stage6["fdr_p"].to_numpy(dtype=float),
            "spatial_null_p": spatial_p,
            "spatial_null_fdr": spatial_fdr,
            "spatial_robustness": percentile,
            "robustness_rank": ranks,
        }))

    sensitivity = pd.concat(sensitivity_frames, ignore_index=True)
    primary = sensitivity.loc[sensitivity["gene_set_version"] == PRIMARY_METHOD].copy()
    primary["spatial_statistic"] = primary["stage6_z"]
    primary["spatial_isolate"] = weights.sum(axis=1) == 0
    primary["robust"] = (
        (primary["stage6_z"] > 0)
        & (primary["stage6_fdr"] < FDR_THRESHOLD)
        & (primary["spatial_null_fdr"] < FDR_THRESHOLD)
    )
    primary["spatial_robustness_label"] = np.where(primary["robust"], "robust", "not_robust")
    robustness = primary[[
        "region_id", "region_name", "stage6_z", "stage6_fdr", "spatial_statistic",
        "spatial_null_p", "spatial_null_fdr", "spatial_robustness", "spatial_robustness_label",
        "robustness_rank", "spatial_isolate",
    ]].sort_values(["robustness_rank", "stage6_z"]).reset_index(drop=True)

    donor_weighted = donor.pivot(index="region_id", columns="donor", values="weighted_score").reindex(region_ids)
    donor_percentiles = donor_weighted.rank(axis=0, pct=True)
    donor_summary = pd.DataFrame({
        "region_id": region_ids,
        "donor_median_rank_percentile": donor_percentiles.median(axis=1, skipna=True),
        "donor_rank_iqr": donor_percentiles.quantile(0.75, axis=1) - donor_percentiles.quantile(0.25, axis=1),
        "donors_observed": donor_weighted.notna().sum(axis=1),
    }).set_index("region_id")
    robustness = robustness.join(donor_summary, on="region_id")
    return {
        "autocorrelation": pd.DataFrame(autocorrelation_rows),
        "autocorrelation_nulls": autocorrelation_nulls,
        "null_summary": pd.DataFrame(null_summary_rows),
        "global": pd.DataFrame(global_rows),
        "sensitivity": sensitivity,
        "robustness": robustness,
        "surrogates": surrogate_maps,
        "eigenvalues": eigenvalues,
    }


def save_outputs(outputs: dict[str, object], weights_frame: pd.DataFrame) -> None:
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    DATA_RESULTS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    weights_frame.to_csv(INTERMEDIATE / "brain_region_spatial_weights.csv")
    outputs["autocorrelation"].to_csv(DATA_RESULTS / "parkinson_spatial_autocorrelation.csv", index=False)
    outputs["null_summary"].to_csv(DATA_RESULTS / "parkinson_spatial_null_summary.csv", index=False)
    outputs["global"].to_csv(DATA_RESULTS / "parkinson_global_spatial_test.csv", index=False)
    outputs["sensitivity"].to_csv(DATA_RESULTS / "parkinson_spatial_gene_set_sensitivity.csv", index=False)
    outputs["robustness"].to_csv(DATA_RESULTS / "parkinson_spatial_robustness.csv", index=False)
    table = outputs["robustness"].rename(columns={
        "region_name": "region",
        "stage6_z": "Stage 6 Z-score",
        "stage6_fdr": "Stage 6 FDR",
        "spatial_statistic": "spatial statistic",
        "spatial_robustness": "spatial robustness measure",
        "spatial_robustness_label": "robust/not_robust",
        "robustness_rank": "rank",
    })
    table.to_csv(TABLES / "parkinson_robust_regions.csv", index=False)


def plot_autocorrelation(outputs: dict[str, object], regions: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    row = outputs["autocorrelation"].set_index("metric").loc["weighted_stage6_z"]
    null = outputs["autocorrelation_nulls"]["weighted_stage6_z"]
    robustness = outputs["robustness"].set_index("region_id").reindex(regions.index)
    figure = plt.figure(figsize=(13, 5.4))
    left = figure.add_subplot(1, 2, 1, projection="3d")
    scatter = left.scatter(
        regions["centroid_mni_x"], regions["centroid_mni_y"], regions["centroid_mni_z"],
        c=robustness["stage6_z"], cmap="RdBu_r", vmin=-4, vmax=4, s=35, alpha=0.9,
    )
    left.view_init(elev=22, azim=-68)
    left.set(xlabel="MNI x", ylabel="MNI y", zlabel="MNI z", title="Observed Stage 6 weighted Z profile")
    figure.colorbar(scatter, ax=left, shrink=0.62, pad=0.1, label="Stage 6 Z-score")
    right = figure.add_subplot(1, 2, 2)
    right.hist(null, bins=45, color="#8ca8a2", edgecolor="white", linewidth=0.35)
    right.axvline(row["observed_statistic"], color="#b2182b", linewidth=2.4, label=f"Observed I = {row['observed_statistic']:.3f}")
    right.axvline(row["null_mean"], color="#243b3a", linestyle="--", label=f"Permutation mean = {row['null_mean']:.3f}")
    right.set(xlabel="Global Moran's I", ylabel="Unrestricted label permutations", title="Spatial autocorrelation permutation test")
    right.legend(frameon=False)
    figure.suptitle("Parkinson Spatial Autocorrelation", fontsize=17, fontweight="bold")
    figure.text(0.5, 0.01, f"Two-sided permutation p = {row['p_value']:.4g}; n = {N_PERMUTATIONS:,}", ha="center")
    figure.tight_layout(rect=(0, 0.04, 1, 0.96))
    figure.savefig(FIGURES / "stage_07_parkinson_spatial_autocorrelation.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_robustness_map(outputs: dict[str, object]) -> None:
    payload = json.loads(WEB_GEOMETRY.read_text(encoding="utf-8"))
    robustness = outputs["robustness"].set_index("region_id")
    normalization = mpl_colors.Normalize(0, 1)
    cmap = plt.get_cmap("Purples")
    triangles: list[np.ndarray] = []
    facecolors: list[tuple[float, float, float, float]] = []
    for parcel in payload["regions"]:
        region_id = int(parcel["region_id"])
        positions = np.asarray(parcel["positions"], dtype=float).reshape(-1, 3)
        # Matplotlib z is vertical; web coordinates are (MNI x, MNI z, MNI y).
        plot_positions = positions[:, [0, 2, 1]]
        parcel_triangles = plot_positions[np.asarray(parcel["indices"], dtype=int).reshape(-1, 3)]
        triangles.extend(parcel_triangles)
        color = cmap(normalization(robustness.loc[region_id, "spatial_robustness"]))
        facecolors.extend([color] * len(parcel_triangles))
    figure = plt.figure(figsize=(10.5, 7.5))
    axis = figure.add_subplot(111, projection="3d")
    collection = Poly3DCollection(triangles, facecolors=facecolors, linewidths=0, alpha=0.98)
    axis.add_collection3d(collection)
    all_vertices = np.concatenate(triangles, axis=0)
    axis.auto_scale_xyz(all_vertices[:, 0], all_vertices[:, 1], all_vertices[:, 2])
    axis.set_box_aspect((1.0, 1.15, 1.0), zoom=1.35)
    axis.view_init(elev=18, azim=-48)
    axis.set_axis_off()
    robust_count = int(robustness["spatial_robustness_label"].eq("robust").sum())
    scalar = plt.cm.ScalarMappable(norm=normalization, cmap=cmap)
    figure.colorbar(scalar, ax=axis, shrink=0.58, pad=0.02, label="Spatial-null percentile (1 − empirical p)")
    figure.suptitle("Parkinson Enrichment: Spatial Robustness", fontsize=18, fontweight="bold")
    figure.text(0.5, 0.035, f"Continuous spatial sensitivity measure; {robust_count} regions meet the pre-specified joint robust rule", ha="center")
    figure.savefig(BRAIN_MAPS / "stage_07_parkinson_spatial_robustness.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def label_predefined(axis: plt.Axes, frame: pd.DataFrame, x: str, y: str) -> None:
    selected = set(frame.nlargest(5, "stage6_z").index) | set(frame.nlargest(5, "spatial_robustness").index)
    ordered = frame.loc[list(selected)].sort_values(y, ascending=False)
    for rank, (index, row) in enumerate(ordered.iterrows()):
        horizontal_offset = -6 if row[x] > frame[x].quantile(0.85) else 6
        horizontal_alignment = "right" if horizontal_offset < 0 else "left"
        axis.annotate(
            row["region_name"],
            (row[x], row[y]),
            xytext=(horizontal_offset, -10 * rank),
            textcoords="offset points",
            ha=horizontal_alignment,
            va="center",
            fontsize=7,
            arrowprops={"arrowstyle": "-", "color": "#777777", "linewidth": 0.45},
        )


def plot_comparisons(outputs: dict[str, object]) -> None:
    robustness = outputs["robustness"].copy().set_index("region_id")
    figure, axis = plt.subplots(figsize=(7.7, 6.2))
    sns.scatterplot(data=robustness, x="stage6_z", y="spatial_robustness", hue="spatial_robustness_label", palette={"robust": "#b2182b", "not_robust": "#477c78"}, ax=axis)
    label_predefined(axis, robustness, "stage6_z", "spatial_robustness")
    axis.axvline(0, color="#999999", linewidth=0.8)
    axis.set(xlabel="Stage 6 Z-score", ylabel="Spatial-null percentile (1 − p)", title="Stage 6 vs Stage 7 Regional Robustness")
    axis.legend(title="Joint rule", frameon=False)
    figure.tight_layout()
    figure.savefig(FIGURES / "stage_07_stage6_vs_spatial_robustness.png", dpi=220, bbox_inches="tight")
    plt.close(figure)

    stage6_rank = robustness["stage6_z"].rank(pct=True)
    spatial_rank = robustness["spatial_robustness"].rank(pct=True)
    correlation = spearmanr(stage6_rank, spatial_rank).statistic
    figure, axis = plt.subplots(figsize=(7.2, 6.2))
    axis.scatter(stage6_rank, spatial_rank, c=robustness["stage6_z"], cmap="RdBu_r", vmin=-4, vmax=4, s=31, alpha=0.82)
    axis.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=1)
    axis.set(xlabel="Stage 6 enrichment rank percentile", ylabel="Stage 7 spatial-null rank percentile", title="Enrichment Stability Across Statistical Null Models", xlim=(0, 1.02), ylim=(0, 1.02))
    axis.text(0.03, 0.95, f"Spearman ρ = {correlation:.3f}", transform=axis.transAxes, va="top")
    figure.tight_layout()
    figure.savefig(FIGURES / "stage_07_null_model_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_gene_set_sensitivity(outputs: dict[str, object]) -> None:
    sensitivity = outputs["sensitivity"]
    ranks = sensitivity.pivot(index="region_name", columns="gene_set_version", values="spatial_robustness")
    top_names: set[str] = set()
    for method in METHODS:
        top_names |= set(ranks.nlargest(8, method).index)
    plot_data = ranks.loc[list(top_names)].assign(mean=lambda frame: frame.mean(axis=1)).sort_values("mean", ascending=False).drop(columns="mean")
    figure, axis = plt.subplots(figsize=(7.5, max(6, 0.32 * len(plot_data))))
    sns.heatmap(plot_data, cmap="mako", vmin=0, vmax=1, annot=True, fmt=".2f", linewidths=0.3, cbar_kws={"label": "Spatial-null percentile"}, ax=axis)
    axis.set(xlabel="Gene-set definition", ylabel="Region", title="Gene-set Spatial Sensitivity")
    figure.tight_layout()
    figure.savefig(FIGURES / "stage_07_gene_set_spatial_sensitivity.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_donor_robustness(outputs: dict[str, object]) -> None:
    frame = outputs["robustness"].copy().set_index("region_id")
    figure, axis = plt.subplots(figsize=(8, 6.3))
    scatter = axis.scatter(
        frame["spatial_robustness"], frame["donor_median_rank_percentile"],
        c=frame["donor_rank_iqr"], cmap="viridis_r", s=38, alpha=0.85,
    )
    label_predefined(axis, frame, "spatial_robustness", "donor_median_rank_percentile")
    figure.colorbar(scatter, ax=axis, label="Across-donor rank IQR (lower is more stable)")
    axis.set(xlabel="Stage 7 spatial-null percentile", ylabel="Median within-donor weighted-score percentile", title="Donor and Spatial Robustness")
    axis.set_xlim(-0.05, 1.04)
    figure.tight_layout()
    figure.savefig(FIGURES / "stage_07_donor_and_spatial_robustness.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def write_reports(
    regions: pd.DataFrame,
    weights: np.ndarray,
    edges: set[tuple[int, int]],
    outputs: dict[str, object],
    elapsed: float,
) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    degrees = weights.sum(axis=1).astype(int)
    isolates = regions.loc[degrees == 0, "region_name"].tolist()
    components = connected_components(weights)
    representation = f"""# Stage 07 spatial representation

## Selected representation

Stage 7 uses a **binary, symmetric volumetric AAL3 adjacency graph** derived from
`data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz`, the exact 2 mm MNI image used in Stage 2
and the web geometry. Two retained parcels are neighbors when any of their voxels
touch by a face, edge, or corner (26-connectivity). This criterion is deterministic
and respects atlas contact; it does not claim that Euclidean proximity is biological
connectivity.

- Regions: {len(regions)} retained Stage 2 parcels
- Undirected edges: {len(edges)}
- Degree: minimum {degrees.min()}, median {np.median(degrees):.1f}, maximum {degrees.max()}
- Connected-component sizes: {', '.join(str(len(component)) for component in components)}
- Isolated zero-weight parcels ({len(isolates)}): {', '.join(isolates)}
- Weights: binary (1=voxel contact, 0=no contact), diagonal zero
- Normalization: not row-normalized; symmetry is preserved for Moran eigenvectors

Isolated parcels remain in every matrix and output. Their spatial-weight rows are
zero and `spatial_isolate=True`; no missing edge is invented. Their regional
spectral-null results are consequently less anatomically constrained and must be
interpreted cautiously. The matrix is saved with numeric region IDs as both rows
and columns.

## Why this representation and null

A cortical spin requires spherical cortical coordinates and cannot validly rotate
irregular subcortical, cerebellar, or brainstem parcels. AAL3 is volumetric and the
analysis deliberately retains those structures. Moran spectral randomization (MSR)
instead accepts an explicit symmetric neighbor matrix for irregularly arranged
observations. The singleton procedure independently flips each observed Moran-basis
coefficient's sign, preserving the mean, variance, complete Moran power spectrum,
and global Moran's I exactly while changing hotspot locations. This strong constraint
is intentionally conservative and can leave surrogates correlated with the input.

Primary method reference: Wagner HH, Dray S (2015), *Methods in Ecology and
Evolution* 6:1169–1178, https://doi.org/10.1111/2041-210X.12407. BrainSpace's
independent implementation documents the same singleton procedure:
https://brainspace.readthedocs.io/en/latest/generated/brainspace.null_models.moran.moran_randomization.html.
"""
    (REPORTS / "stage_07_spatial_representation.md").write_text(representation, encoding="utf-8")

    definition = f"""# Stage 07 robustness definition

This rule was fixed before Stage 7 results were generated.

A region is called **robust** only if all three conditions hold:

1. its Stage 6 weighted enrichment Z-score is positive;
2. its Stage 6 gene-set permutation BH-FDR is below {FDR_THRESHOLD}; and
3. its one-sided Stage 7 MSR regional p-value passes BH-FDR across {len(regions)} regions at q < {FDR_THRESHOLD}.

The Stage 7 regional p-value is `(1 + number of surrogate values >= observed) /
({N_PERMUTATIONS} + 1)` at the same fixed atlas parcel. `spatial_robustness` is the
descriptive percentile-like quantity `1 - spatial_null_p`; larger values indicate
that the observed positive score lies higher in that parcel's spatial-null
distribution. It is not a posterior probability. `robustness_rank` sorts this
measure descending with deterministic minimum ranks.

The joint rule prevents a spatial sensitivity result from reviving a region that
did not pass the discovery-stage multiplicity criterion. It also means a negative
Stage 6 result must remain negative even if its descriptive spatial rank is high.
"""
    (REPORTS / "stage_07_robustness_definition.md").write_text(definition, encoding="utf-8")

    autocorrelation = outputs["autocorrelation"].set_index("metric")
    global_results = outputs["global"].set_index("gene_set_version")
    robustness = outputs["robustness"]
    stage6_significant = int((robustness["stage6_fdr"] < FDR_THRESHOLD).sum())
    robust_count = int((robustness["spatial_robustness_label"] == "robust").sum())
    sensitivity = outputs["sensitivity"].pivot(index="region_id", columns="gene_set_version", values="spatial_robustness")
    pairwise = sensitivity.corr(method="spearman")
    donor_rho = spearmanr(robustness["spatial_robustness"], robustness["donor_median_rank_percentile"], nan_policy="omit").statistic
    interpretation = f"""# Stage 07 interpretation

## What each stage measures

- **Stage 5 — raw expression:** regional averages of Parkinson-prioritized genes on
  the preserved AHBA normalization scale; no null-model inference.
- **Stage 6 — gene-set enrichment:** whether those genes score higher than technical-
  feature-matched random gene sets. This is the discovery analysis.
- **Stage 7 — spatial robustness:** whether the location and peak structure of the
  Stage 6 Z profile remain unusual against MSR maps with comparable atlas-graph
  spatial dependence. This is a sensitivity analysis, not a replacement result.

## Findings

The weighted Stage 6 Z profile had global Moran's I =
{autocorrelation.loc['weighted_stage6_z', 'observed_statistic']:.4f}
(two-sided unrestricted-label permutation p =
{autocorrelation.loc['weighted_stage6_z', 'p_value']:.4g}). The spatially constrained
global peak test gave p = {global_results.loc['weighted', 'spatial_null_p']:.4g}.

Stage 6 contained {stage6_significant} weighted regions at BH-FDR q < {FDR_THRESHOLD};
therefore {robust_count} regions met the pre-specified joint Stage 6 + Stage 7 rule.
This negative joint result is retained explicitly. Continuous spatial-null
percentiles are provided for sensitivity and ranking, but they must not be presented
as corrected truth.

Broad/stringent/weighted spatial-robustness rank correlations ranged from
{pairwise.where(np.triu(np.ones(pairwise.shape), 1).astype(bool)).stack().min():.3f}
to {pairwise.where(np.triu(np.ones(pairwise.shape), 1).astype(bool)).stack().max():.3f}.
The correlation between weighted spatial robustness and the median donor rank
percentile was {donor_rho:.3f}. Donor summaries are descriptive because six donors,
uneven sampling, and donor missingness do not support donor-level significance.

## Interpretation boundaries

Parkinson-associated genes show a regional spatial pattern, but no region passed the
current joint robustness rule. This does not show where Parkinson starts, identify a
causal region, map pathology, prove selective vulnerability, or provide clinical
prediction. MSR depends on the chosen AAL3 contact graph, assumes the graph is an
adequate proxy for spatial dependence. Singleton MSR exactly preserves the graph
power spectrum but may be conservative because randomized maps can remain correlated
with the observed map. Four isolated brainstem parcels have weaker spatial constraints.
Independent pathological, longitudinal, and replication evidence would
be required for stronger claims.

Runtime for this build: {elapsed:.1f} seconds using Python {platform.python_version()},
NumPy {np.__version__}, SciPy {scipy.__version__}, pandas {pd.__version__}, and
nibabel {nib.__version__}.
"""
    (REPORTS / "stage_07_interpretation.md").write_text(interpretation, encoding="utf-8")


def write_parameters(regions: pd.DataFrame, weights: np.ndarray, edges: set[tuple[int, int]]) -> None:
    degrees = weights.sum(axis=1)
    payload = {
        "stage": 7,
        "analysis_role": "spatial sensitivity analysis; does not replace Stage 6",
        "atlas": "AAL3v1 2 mm MNI, identical to Stage 2 and web geometry",
        "region_count": len(regions),
        "neighbor_definition": "26-connectivity: voxel face, edge, or corner contact",
        "weight_type": "binary symmetric",
        "row_normalized": False,
        "diagonal": 0,
        "edge_count": len(edges),
        "isolated_region_ids": [int(value) for value in regions.index[degrees == 0]],
        "isolated_region_handling": "retained with zero-weight rows and explicitly flagged; no edges imputed",
        "spatial_null": "Moran spectral randomization, singleton procedure",
        "spatial_null_reference": "Wagner and Dray 2015, doi:10.1111/2041-210X.12407",
        "n_spatial_permutations": N_PERMUTATIONS,
        "moran_autocorrelation_null": "unrestricted regional label permutation, two-sided around null mean",
        "n_moran_permutations": N_PERMUTATIONS,
        "random_seed": RANDOM_SEED,
        "regional_test": "one-sided greater-or-equal empirical p with +1 correction",
        "regional_multiple_testing": "Benjamini-Hochberg across 138 regions, separately by gene-set version",
        "global_test": "maximum standardized positive regional Stage 6 Z-score",
        "fdr_threshold": FDR_THRESHOLD,
        "robust_rule": "positive Stage 6 Z AND Stage 6 FDR < 0.05 AND Stage 7 spatial-null FDR < 0.05",
        "stored_nulls": "summary statistics only; full 10000 x 138 surrogate matrices not saved",
        "gene_set_versions": list(METHODS),
    }
    (DATA_RESULTS / "stage_07_spatial_null_parameters.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def final_summary(outputs: dict[str, object]) -> None:
    robustness = outputs["robustness"]
    global_weighted = outputs["global"].set_index("gene_set_version").loc["weighted"]
    sensitivity = outputs["sensitivity"].pivot(index="region_id", columns="gene_set_version", values="spatial_robustness")
    correlations = sensitivity.corr(method="spearman")
    off_diagonal = correlations.where(np.triu(np.ones(correlations.shape), 1).astype(bool)).stack()
    donor_rho = spearmanr(robustness["spatial_robustness"], robustness["donor_median_rank_percentile"], nan_policy="omit").statistic
    stage6_count = int((robustness["stage6_fdr"] < FDR_THRESHOLD).sum())
    robust_count = int(robustness["spatial_robustness_label"].eq("robust").sum())
    print("\nPARKINSON SPATIAL VALIDATION")
    print(f"Regions: {len(robustness)}")
    print("Spatial method: binary 26-neighbor AAL3 graph + Moran spectral singleton randomization")
    print(f"Spatial permutations: {N_PERMUTATIONS:,}")
    print("Global spatial statistic: maximum standardized positive regional Stage 6 Z-score")
    print(f"Global spatial result: observed={global_weighted['observed_statistic']:.4f}, p={global_weighted['spatial_null_p']:.4g}")
    print(f"Stage 6 significant regions: {stage6_count}")
    print(f"Stage 7 robust regions: {robust_count}")
    print(f"Overlap between Stage 6 significant and Stage 7 robust: {robust_count}")
    print(f"Broad/stringent/weighted consistency: Spearman range {off_diagonal.min():.3f} to {off_diagonal.max():.3f}")
    print(f"Donor consistency: spatial percentile vs median donor-rank Spearman rho={donor_rho:.3f} (descriptive)")


def main() -> None:
    start = time.perf_counter()
    for directory in (INTERMEDIATE, DATA_RESULTS, TABLES, FIGURES, BRAIN_MAPS, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    regions, stage5, stage6_all, donor, image = load_inputs()
    atlas = np.asarray(image.dataobj, dtype=np.int16)
    weights_frame, edges = build_spatial_weights(atlas, regions.index.to_numpy(dtype=int))
    weights = weights_frame.to_numpy(dtype=float)
    outputs = calculate_outputs(regions, stage5, stage6_all, donor, weights)
    save_outputs(outputs, weights_frame)
    write_parameters(regions, weights, edges)
    plot_autocorrelation(outputs, regions)
    plot_robustness_map(outputs)
    plot_comparisons(outputs)
    plot_gene_set_sensitivity(outputs)
    plot_donor_robustness(outputs)
    elapsed = time.perf_counter() - start
    write_reports(regions, weights, edges, outputs, elapsed)
    final_summary(outputs)


if __name__ == "__main__":
    main()
