"""Stage 8: external Parkinson neuroanatomical phenotype validation.

The Stage 6/7 discovery outputs are immutable inputs.  ENIGMA-PD summary
statistics are mapped to AAL3 with a hand-audited anatomical crosswalk, and
inference is performed once per ENIGMA parcel (never once per propagated AAL3
parcel).  Propagated AAL3 values are produced only for atlas visualization.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import urllib.request
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mpl_colors
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import nibabel as nib
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import null_space
from scipy.stats import pearsonr, rankdata, spearmanr, t
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "data" / "validation"
RAW = VALIDATION / "raw"
PROCESSED = VALIDATION / "processed"
DATA_RESULTS = ROOT / "data" / "results"
FIGURES = ROOT / "results" / "figures"
BRAIN_MAPS = ROOT / "results" / "brain_maps"
REPORTS = ROOT / "reports"
REGIONS_FILE = ROOT / "data" / "processed" / "region_metadata.csv"
EXPRESSION_FILE = ROOT / "data" / "processed" / "brain_region_gene_expression.csv"
ATLAS_FILE = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.nii.gz"
GEOMETRY_FILE = ROOT / "web" / "public" / "data" / "aal3_regions.json"
STAGE6_FILE = DATA_RESULTS / "parkinson_regional_enrichment_all_methods.csv"
STAGE6_PRIMARY_FILE = DATA_RESULTS / "parkinson_regional_enrichment.csv"
STAGE7_FILE = DATA_RESULTS / "parkinson_spatial_robustness.csv"
STAGE5_FILE = DATA_RESULTS / "parkinson_regional_raw_scores.csv"
WEIGHTS_FILE = ROOT / "data" / "intermediate" / "brain_region_spatial_weights.csv"

ENIGMA_COMMIT = "b08974b55243060cbc1fad12c87048037446e8f7"
ENIGMA_BASE = (
    "https://raw.githubusercontent.com/MICA-MNI/ENIGMA/"
    f"{ENIGMA_COMMIT}/enigmatoolbox/datasets/summary_statistics"
)
RAW_FILES = {
    "cortical_thickness": "parkinsons_case-controls_CortThick_PDvsCN.csv",
    "subcortical_volume": "parkinsons_case-controls_Subvol_PDvsCN.csv",
}
N_BOOTSTRAP = 10_000
N_SPATIAL = 10_000
RANDOM_SEED = 20260914


# One dictionary entry is one independent ENIGMA validation unit. Values are
# AAL3 labels averaged to represent it. Empty/missing entries are excluded, not
# forced. Exact names are deliberately used instead of fuzzy string matching.
MAPPING: dict[str, tuple[list[int], str, str]] = {
    "L_caudalanteriorcingulate": ([155], "high", "DK caudal anterior cingulate -> AAL3 superior ACC, left"),
    "R_caudalanteriorcingulate": ([156], "high", "DK caudal anterior cingulate -> AAL3 superior ACC, right"),
    "L_cuneus": ([49], "high", "same named cortical structure and hemisphere"),
    "L_fusiform": ([59], "high", "same named cortical structure and hemisphere"),
    "R_fusiform": ([60], "high", "same named cortical structure and hemisphere"),
    "L_inferiorparietal": ([65, 69], "medium", "DK inferior parietal represented by AAL3 inferior parietal plus angular parcels"),
    "R_inferiorparietal": ([66, 70], "medium", "DK inferior parietal represented by AAL3 inferior parietal plus angular parcels"),
    "L_inferiortemporal": ([93], "high", "same named cortical structure and hemisphere"),
    "R_inferiortemporal": ([94], "high", "same named cortical structure and hemisphere"),
    "L_lateraloccipital": ([53, 55, 57], "medium", "DK lateral occipital represented by retained AAL3 lateral occipital gyri"),
    "R_lateraloccipital": ([56, 58], "medium", "retained AAL3 middle and inferior occipital gyri; superior right was not retained in discovery"),
    "L_lateralorbitofrontal": ([31], "high", "DK lateral orbitofrontal -> AAL3 lateral orbitofrontal"),
    "R_lateralorbitofrontal": ([32], "high", "DK lateral orbitofrontal -> AAL3 lateral orbitofrontal"),
    "L_lingual": ([51], "high", "same named cortical structure and hemisphere"),
    "R_lingual": ([52], "high", "same named cortical structure and hemisphere"),
    "L_medialorbitofrontal": ([21, 23, 25], "medium", "DK medial orbitofrontal represented by medial orbital frontal, rectus and medial OFC"),
    "R_medialorbitofrontal": ([22, 24, 26], "medium", "DK medial orbitofrontal represented by medial orbital frontal, rectus and medial OFC"),
    "L_middletemporal": ([89], "high", "same named cortical structure and hemisphere"),
    "R_middletemporal": ([90], "high", "same named cortical structure and hemisphere"),
    "L_parahippocampal": ([43], "high", "same named cortical structure and hemisphere"),
    "R_parahippocampal": ([44], "high", "same named cortical structure and hemisphere"),
    "L_paracentral": ([73], "high", "DK paracentral -> AAL3 paracentral lobule"),
    "R_paracentral": ([74], "high", "DK paracentral -> AAL3 paracentral lobule"),
    "L_parsopercularis": ([7], "high", "DK pars opercularis -> AAL3 inferior frontal opercular"),
    "R_parsopercularis": ([8], "high", "DK pars opercularis -> AAL3 inferior frontal opercular"),
    "L_parsorbitalis": ([11], "high", "DK pars orbitalis -> AAL3 inferior frontal orbital"),
    "R_parsorbitalis": ([12], "high", "DK pars orbitalis -> AAL3 inferior frontal orbital"),
    "L_parstriangularis": ([9], "high", "DK pars triangularis -> AAL3 inferior frontal triangular"),
    "R_parstriangularis": ([10], "high", "DK pars triangularis -> AAL3 inferior frontal triangular"),
    "L_pericalcarine": ([47], "high", "DK pericalcarine -> AAL3 calcarine"),
    "R_pericalcarine": ([48], "high", "DK pericalcarine -> AAL3 calcarine"),
    "L_postcentral": ([61], "high", "same named cortical structure and hemisphere"),
    "R_postcentral": ([62], "high", "same named cortical structure and hemisphere"),
    "L_posteriorcingulate": ([39], "high", "same named cortical structure and hemisphere"),
    "L_precentral": ([1], "high", "same named cortical structure and hemisphere"),
    "R_precentral": ([2], "high", "same named cortical structure and hemisphere"),
    "L_precuneus": ([71], "high", "same named cortical structure and hemisphere"),
    "R_precuneus": ([72], "high", "same named cortical structure and hemisphere"),
    "L_rostralanteriorcingulate": ([153], "high", "DK rostral anterior cingulate -> AAL3 pregenual ACC"),
    "R_rostralanteriorcingulate": ([154], "high", "DK rostral anterior cingulate -> AAL3 pregenual ACC"),
    "L_rostralmiddlefrontal": ([5], "medium", "DK rostral middle frontal -> AAL3 middle frontal"),
    "R_rostralmiddlefrontal": ([6], "medium", "DK rostral middle frontal -> AAL3 middle frontal"),
    "L_superiorfrontal": ([3, 19], "medium", "DK superior frontal represented by lateral and medial superior frontal AAL3 parcels"),
    "R_superiorfrontal": ([4, 20], "medium", "DK superior frontal represented by lateral and medial superior frontal AAL3 parcels"),
    "L_superiorparietal": ([63], "high", "same named cortical structure and hemisphere"),
    "R_superiorparietal": ([64], "high", "same named cortical structure and hemisphere"),
    "L_superiortemporal": ([85], "high", "same named cortical structure and hemisphere"),
    "R_superiortemporal": ([86], "high", "same named cortical structure and hemisphere"),
    "L_supramarginal": ([67], "high", "same named cortical structure and hemisphere"),
    "R_supramarginal": ([68], "high", "same named cortical structure and hemisphere"),
    "L_temporalpole": ([87, 91], "medium", "DK temporal pole represented by retained superior and middle temporal-pole AAL3 parcels"),
    "R_temporalpole": ([88], "medium", "retained superior temporal-pole parcel; middle right was not retained in discovery"),
    "L_transversetemporal": ([83], "high", "DK transverse temporal -> AAL3 Heschl gyrus"),
    "R_transversetemporal": ([84], "high", "DK transverse temporal -> AAL3 Heschl gyrus"),
    "L_insula": ([33], "high", "same named cortical structure and hemisphere"),
    "R_insula": ([34], "high", "same named cortical structure and hemisphere"),
    "Laccumb": ([157], "high", "FreeSurfer accumbens -> AAL3 nucleus accumbens"),
    "Raccumb": ([158], "high", "FreeSurfer accumbens -> AAL3 nucleus accumbens"),
    "Lamyg": ([45], "high", "FreeSurfer amygdala -> AAL3 amygdala"),
    "Ramyg": ([46], "high", "FreeSurfer amygdala -> AAL3 amygdala"),
    "Lcaud": ([75], "high", "FreeSurfer caudate -> AAL3 caudate"),
    "Rcaud": ([76], "high", "FreeSurfer caudate -> AAL3 caudate"),
    "Lhippo": ([41], "high", "FreeSurfer hippocampus -> AAL3 hippocampus"),
    "Rhippo": ([42], "high", "FreeSurfer hippocampus -> AAL3 hippocampus"),
    "Lpal": ([79], "high", "FreeSurfer pallidum -> AAL3 pallidum"),
    "Rpal": ([80], "high", "FreeSurfer pallidum -> AAL3 pallidum"),
    "Lput": ([77], "high", "FreeSurfer putamen -> AAL3 putamen"),
    "Rput": ([78], "high", "FreeSurfer putamen -> AAL3 putamen"),
    "Lthal": ([121, 125, 127, 129, 131, 135, 137, 139, 141, 143, 145, 147], "medium", "FreeSurfer whole thalamus represented by all retained left AAL3 thalamic nuclei"),
    "Rthal": ([122, 128, 130, 136, 140, 144, 146, 148], "medium", "FreeSurfer whole thalamus represented by all retained right AAL3 thalamic nuclei"),
}


def ensure_directories() -> None:
    for directory in (RAW, PROCESSED, DATA_RESULTS, FIGURES, BRAIN_MAPS, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)


def download_raw() -> dict[str, dict[str, object]]:
    provenance: dict[str, dict[str, object]] = {}
    for role, filename in RAW_FILES.items():
        url = f"{ENIGMA_BASE}/{filename}"
        destination = RAW / filename
        if not destination.exists():
            urllib.request.urlretrieve(url, destination)
        content = destination.read_bytes()
        provenance[role] = {
            "filename": filename,
            "url": url,
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
    return provenance


def read_validation() -> pd.DataFrame:
    cortical = pd.read_csv(RAW / RAW_FILES["cortical_thickness"]).rename(
        columns={"Structure": "validation_region"}
    )
    cortical["modality"] = "cortical_thickness"
    subcortical = pd.read_csv(RAW / RAW_FILES["subcortical_volume"]).rename(
        columns={"Structure": "validation_region"}
    )
    subcortical["modality"] = "subcortical_volume"
    frame = pd.concat([cortical, subcortical], ignore_index=True)
    frame["d_icv"] = pd.to_numeric(frame["d_icv"], errors="raise")
    frame["n_controls"] = pd.to_numeric(frame["n_controls"], errors="raise").astype(int)
    frame["n_patients"] = pd.to_numeric(frame["n_patients"], errors="raise").astype(int)
    frame["validation_score"] = -frame["d_icv"]
    frame["original_unit"] = "Cohen's d (PD minus control); sign-reversed so higher means thinner/smaller in PD"
    frame["sample_size_if_available"] = frame["n_controls"] + frame["n_patients"]
    if frame["validation_region"].duplicated().any() or len(frame) != 84:
        raise RuntimeError("Unexpected ENIGMA regional schema")
    return frame


def load_discovery() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    regions = pd.read_csv(REGIONS_FILE)
    regions = regions.loc[regions["retained_in_main_matrix"]].set_index("region_id").sort_index()
    stage6 = pd.read_csv(STAGE6_FILE)
    stage6_primary = pd.read_csv(STAGE6_PRIMARY_FILE).set_index("region_id").reindex(regions.index)
    stage7 = pd.read_csv(STAGE7_FILE).set_index("region_id").reindex(regions.index)
    stage5 = pd.read_csv(STAGE5_FILE).set_index("region_id").reindex(regions.index)
    if len(regions) != 138 or stage6_primary["z_score"].isna().any():
        raise RuntimeError("Frozen discovery input audit failed")
    return regions, stage6, stage6_primary, stage7, stage5


def make_mapping(validation: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    all_ids: list[int] = []
    for row in validation.itertuples(index=False):
        mapped = MAPPING.get(row.validation_region)
        if mapped is None:
            reason = "No non-overlapping, anatomically defensible AAL3 equivalent in the retained discovery matrix"
            if "LatVent" in row.validation_region:
                reason = "Lateral ventricle is not a tissue parcel and was excluded a priori"
            rows.append({
                "validation_region": row.validation_region,
                "gene2brain_region": "",
                "mapping_confidence": "excluded",
                "mapping_method": "manual anatomical crosswalk",
                "notes": reason,
                "gene2brain_region_ids": "",
                "included_in_inference": False,
            })
            continue
        ids, confidence, notes = mapped
        missing = sorted(set(ids) - set(regions.index))
        if missing:
            raise RuntimeError(f"Mapping for {row.validation_region} uses unavailable labels: {missing}")
        all_ids.extend(ids)
        names = regions.loc[ids, "region_name"].tolist()
        rows.append({
            "validation_region": row.validation_region,
            "gene2brain_region": " | ".join(names),
            "mapping_confidence": confidence,
            "mapping_method": "documented Desikan-Killiany/FreeSurfer-to-AAL3 anatomical crosswalk",
            "notes": notes,
            "gene2brain_region_ids": "|".join(str(value) for value in ids),
            "included_in_inference": True,
        })
    duplicates = pd.Series(all_ids).value_counts()
    if (duplicates > 1).any():
        raise RuntimeError(f"AAL3 parcels reused across validation units: {duplicates[duplicates > 1].to_dict()}")
    mapping = pd.DataFrame(rows)
    mapping.to_csv(VALIDATION / "parkinson_region_mapping.csv", index=False)
    return mapping


def aggregate(values: pd.Series, ids: list[int]) -> float:
    return float(values.reindex(ids).mean())


def build_analysis_frame(
    validation: pd.DataFrame,
    mapping: pd.DataFrame,
    regions: pd.DataFrame,
    stage6: pd.DataFrame,
    stage6_primary: pd.DataFrame,
    stage7: pd.DataFrame,
    stage5: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    included = mapping.loc[mapping["included_in_inference"]].set_index("validation_region")
    source = validation.set_index("validation_region").loc[included.index]
    profiles = {
        method: stage6.loc[stage6["gene_set_version"] == method].set_index("region_id")["z_score"]
        for method in ("broad", "stringent", "weighted")
    }
    fdr_rank = stage6_primary["fdr_p"].rank(method="average", ascending=False, pct=True)
    expression = pd.read_csv(EXPRESSION_FILE).set_index("region_id").mean(axis=1)
    atlas = np.asarray(nib.load(ATLAS_FILE).dataobj, dtype=np.int16)
    voxel_counts = pd.Series({region_id: int(np.count_nonzero(atlas == region_id)) for region_id in regions.index})
    rows: list[dict[str, object]] = []
    propagated: list[dict[str, object]] = []
    for validation_region, mapping_row in included.iterrows():
        ids = [int(value) for value in mapping_row["gene2brain_region_ids"].split("|")]
        source_row = source.loc[validation_region]
        discovery_z = aggregate(profiles["weighted"], ids)
        record = {
            "region": validation_region,
            "gene2brain_region": mapping_row["gene2brain_region"],
            "gene2brain_region_ids": mapping_row["gene2brain_region_ids"],
            "validation_score": float(source_row["validation_score"]),
            "original_effect_d": float(source_row["d_icv"]),
            "original_unit": source_row["original_unit"],
            "sample_size_if_available": int(source_row["sample_size_if_available"]),
            "n_patients": int(source_row["n_patients"]),
            "n_controls": int(source_row["n_controls"]),
            "modality": source_row["modality"],
            "mapping_confidence": mapping_row["mapping_confidence"],
            "weighted_z": discovery_z,
            "broad_z": aggregate(profiles["broad"], ids),
            "stringent_z": aggregate(profiles["stringent"], ids),
            "fdr_rank_enrichment": aggregate(fdr_rank, ids),
            "spatial_robustness": aggregate(stage7["spatial_robustness"], ids),
            "region_size_voxels": int(voxel_counts.reindex(ids).sum()),
            "log_region_size": math.log1p(voxel_counts.reindex(ids).sum()),
            "baseline_ahba_expression": aggregate(expression, ids),
            "anatomical_class": "cortex" if source_row["modality"] == "cortical_thickness" else "subcortex",
            "n_genes_represented": int(round(aggregate(stage5["n_weighted_genes_present"], ids))),
        }
        rows.append(record)
        for region_id in ids:
            propagated.append({
                "region_id": region_id,
                "region": regions.loc[region_id, "region_name"],
                "validation_region": validation_region,
                "validation_score": record["validation_score"],
                "weighted_z": float(profiles["weighted"].loc[region_id]),
                "mapping_confidence": mapping_row["mapping_confidence"],
            })
    analysis = pd.DataFrame(rows)
    if analysis["region"].duplicated().any() or len(analysis) < 50:
        raise RuntimeError("Too few unique mapped validation units")
    x_median = analysis["weighted_z"].median()
    y_median = analysis["validation_score"].median()
    analysis["agreement_status"] = np.select(
        [
            (analysis["weighted_z"] >= x_median) & (analysis["validation_score"] >= y_median),
            (analysis["weighted_z"] >= x_median) & (analysis["validation_score"] < y_median),
            (analysis["weighted_z"] < x_median) & (analysis["validation_score"] >= y_median),
        ],
        ["high_prediction_high_validation", "high_prediction_low_validation", "low_prediction_high_validation"],
        default="low_prediction_low_validation",
    )
    propagated_frame = pd.DataFrame(propagated).merge(
        analysis[["region", "agreement_status"]].rename(columns={"region": "validation_region"}),
        on="validation_region", how="left", validate="many_to_one",
    )
    analysis.to_csv(PROCESSED / "parkinson_validation_regional_scores.csv", index=False)
    propagated_frame.to_csv(PROCESSED / "parkinson_validation_aal3_scores.csv", index=False)
    return analysis, propagated_frame


def bootstrap_ci(x: np.ndarray, y: np.ndarray, method: str, rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    indices = rng.integers(0, n, size=(N_BOOTSTRAP, n))
    values: list[float] = []
    for sample in indices:
        if np.std(x[sample]) == 0 or np.std(y[sample]) == 0:
            continue
        result = pearsonr(x[sample], y[sample]).statistic if method == "pearson" else spearmanr(x[sample], y[sample]).statistic
        values.append(float(result))
    return tuple(float(value) for value in np.quantile(values, [0.025, 0.975]))


def correlation_row(name: str, x: pd.Series, y: pd.Series, rng: np.random.Generator, *, analysis_group: str = "prespecified_metric") -> dict[str, object]:
    x_values, y_values = x.to_numpy(dtype=float), y.to_numpy(dtype=float)
    pearson = pearsonr(x_values, y_values)
    spearman = spearmanr(x_values, y_values)
    pearson_ci = bootstrap_ci(x_values, y_values, "pearson", rng)
    spearman_ci = bootstrap_ci(x_values, y_values, "spearman", rng)
    return {
        "analysis_group": analysis_group,
        "discovery_metric": name,
        "validation_metric": "ENIGMA-PD sign-reversed Cohen's d",
        "pearson_r": pearson.statistic,
        "pearson_ci_low": pearson_ci[0],
        "pearson_ci_high": pearson_ci[1],
        "pearson_p": pearson.pvalue,
        "spearman_rho": spearman.statistic,
        "spearman_ci_low": spearman_ci[0],
        "spearman_ci_high": spearman_ci[1],
        "spearman_p": spearman.pvalue,
        "n_regions": len(x_values),
        "ci_method": f"{N_BOOTSTRAP} fixed-seed nonparametric region bootstraps",
    }


def residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(values)), covariates])
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def partial_row(frame: pd.DataFrame, covariates: list[str], name: str) -> dict[str, object]:
    x = frame["weighted_z"].to_numpy(dtype=float)
    y = frame["validation_score"].to_numpy(dtype=float)
    controls = frame[covariates].to_numpy(dtype=float)
    rx, ry = residualize(x, controls), residualize(y, controls)
    result = pearsonr(rx, ry)
    degrees = len(frame) - len(covariates) - 2
    statistic = result.statistic * math.sqrt(degrees / max(1e-15, 1 - result.statistic**2))
    return {
        "analysis_group": "confound_control",
        "discovery_metric": name,
        "validation_metric": "ENIGMA-PD sign-reversed Cohen's d",
        "pearson_r": result.statistic,
        "pearson_ci_low": np.nan,
        "pearson_ci_high": np.nan,
        "pearson_p": 2 * t.sf(abs(statistic), degrees),
        "spearman_rho": np.nan,
        "spearman_ci_low": np.nan,
        "spearman_ci_high": np.nan,
        "spearman_p": np.nan,
        "n_regions": len(frame),
        "ci_method": f"partial Pearson residuals; covariates={','.join(covariates)}",
    }


def weighted_correlation(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    weights = weights / weights.sum()
    x_centered = x - np.sum(weights * x)
    y_centered = y - np.sum(weights * y)
    return float(np.sum(weights * x_centered * y_centered) / math.sqrt(np.sum(weights * x_centered**2) * np.sum(weights * y_centered**2)))


def moran_surrogates(values: np.ndarray, weights: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = len(values)
    basis = null_space(np.ones((1, n)))
    reduced = basis.T @ weights @ basis
    eigenvalues, rotation = np.linalg.eigh(reduced)
    eigenvectors = basis @ rotation[:, np.argsort(eigenvalues)[::-1]]
    centered = values - values.mean()
    coefficients = eigenvectors.T @ centered / np.linalg.norm(centered)
    signs = rng.choice((-1.0, 1.0), size=(N_SPATIAL, eigenvectors.shape[1]))
    randomized = signs * coefficients
    return values.mean() + values.std(ddof=1) * math.sqrt(n - 1) * randomized @ eigenvectors.T


def spatial_null_test(frame: pd.DataFrame, stage6_primary: pd.DataFrame, rng: np.random.Generator) -> dict[str, object]:
    weights_frame = pd.read_csv(WEIGHTS_FILE, index_col=0)
    region_ids = stage6_primary.index.to_numpy(dtype=int)
    weights = weights_frame.reindex(index=region_ids, columns=[str(value) for value in region_ids]).to_numpy(dtype=float)
    surrogates = moran_surrogates(stage6_primary["z_score"].to_numpy(dtype=float), weights, rng)
    positions = {region_id: position for position, region_id in enumerate(region_ids)}
    aggregate_indices = [np.array([positions[int(value)] for value in token.split("|")]) for token in frame["gene2brain_region_ids"]]
    aggregated = np.column_stack([surrogates[:, indices].mean(axis=1) for indices in aggregate_indices])
    y = frame["validation_score"].to_numpy(dtype=float)
    y_centered = y - y.mean()
    null_r = ((aggregated - aggregated.mean(axis=1, keepdims=True)) @ y_centered) / (
        np.linalg.norm(aggregated - aggregated.mean(axis=1, keepdims=True), axis=1) * np.linalg.norm(y_centered)
    )
    observed = pearsonr(frame["weighted_z"], frame["validation_score"]).statistic
    p_value = (np.sum(np.abs(null_r) >= abs(observed)) + 1) / (N_SPATIAL + 1)
    pd.DataFrame({"spatial_null_pearson_r": null_r}).describe(percentiles=[0.025, 0.5, 0.975]).to_csv(
        PROCESSED / "parkinson_validation_spatial_null_summary.csv"
    )
    return {
        "analysis_group": "confound_control",
        "discovery_metric": "weighted_z_spatial_moran_null",
        "validation_metric": "ENIGMA-PD sign-reversed Cohen's d",
        "pearson_r": observed,
        "pearson_ci_low": float(np.quantile(null_r, 0.025)),
        "pearson_ci_high": float(np.quantile(null_r, 0.975)),
        "pearson_p": p_value,
        "spearman_rho": np.nan,
        "spearman_ci_low": np.nan,
        "spearman_ci_high": np.nan,
        "spearman_p": np.nan,
        "n_regions": len(frame),
        "ci_method": f"two-sided empirical test against {N_SPATIAL} Stage-7-compatible Moran singleton surrogates; interval is null quantile range",
    }


def calculate_statistics(frame: pd.DataFrame, stage6_primary: pd.DataFrame) -> pd.DataFrame:
    seeds = iter(np.random.SeedSequence(RANDOM_SEED).spawn(20))
    rows = [
        correlation_row("weighted_z_primary", frame["weighted_z"], frame["validation_score"], np.random.default_rng(next(seeds))),
        correlation_row("weighted_fdr_rank", frame["fdr_rank_enrichment"], frame["validation_score"], np.random.default_rng(next(seeds))),
        correlation_row("stage7_spatial_robustness", frame["spatial_robustness"], frame["validation_score"], np.random.default_rng(next(seeds))),
        correlation_row("broad_z", frame["broad_z"], frame["validation_score"], np.random.default_rng(next(seeds))),
        correlation_row("stringent_z", frame["stringent_z"], frame["validation_score"], np.random.default_rng(next(seeds))),
        correlation_row("weighted_z_sensitivity", frame["weighted_z"], frame["validation_score"], np.random.default_rng(next(seeds))),
    ]
    rows.extend([
        partial_row(frame, ["log_region_size"], "weighted_z_partial_log_region_size"),
        partial_row(frame, ["baseline_ahba_expression"], "weighted_z_partial_baseline_ahba_expression"),
        partial_row(frame, ["is_subcortex"], "weighted_z_partial_cortex_subcortex"),
        partial_row(frame, ["log_region_size", "baseline_ahba_expression", "is_subcortex"], "weighted_z_partial_joint_available_covariates"),
    ])
    for anatomical_class, subset in frame.groupby("anatomical_class"):
        rows.append(correlation_row(
            f"weighted_z_{anatomical_class}_only", subset["weighted_z"], subset["validation_score"],
            np.random.default_rng(next(seeds)), analysis_group="anatomical_sensitivity",
        ))
    reliability_r = weighted_correlation(
        frame["weighted_z"].to_numpy(), frame["validation_score"].to_numpy(), frame["sample_size_if_available"].to_numpy()
    )
    rows.append({
        "analysis_group": "confound_control", "discovery_metric": "weighted_z_sample_size_weighted",
        "validation_metric": "ENIGMA-PD sign-reversed Cohen's d", "pearson_r": reliability_r,
        "pearson_ci_low": np.nan, "pearson_ci_high": np.nan, "pearson_p": np.nan,
        "spearman_rho": np.nan, "spearman_ci_low": np.nan, "spearman_ci_high": np.nan,
        "spearman_p": np.nan, "n_regions": len(frame),
        "ci_method": "descriptive weighted Pearson using region-specific analyzed sample size as reliability weight",
    })
    rows.append(spatial_null_test(frame, stage6_primary, np.random.default_rng(next(seeds))))
    results = pd.DataFrame(rows)
    results.to_csv(DATA_RESULTS / "parkinson_independent_validation.csv", index=False)
    return results


def sensitivity_outputs(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    loo_rows = []
    for index, row in frame.iterrows():
        subset = frame.drop(index=index)
        loo_rows.append({
            "left_out_region": row["region"],
            "left_out_gene2brain_region": row["gene2brain_region"],
            "pearson_r": pearsonr(subset["weighted_z"], subset["validation_score"]).statistic,
            "pearson_p": pearsonr(subset["weighted_z"], subset["validation_score"]).pvalue,
            "spearman_rho": spearmanr(subset["weighted_z"], subset["validation_score"]).statistic,
            "spearman_p": spearmanr(subset["weighted_z"], subset["validation_score"]).pvalue,
            "n_regions": len(subset),
        })
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(DATA_RESULTS / "parkinson_leave_one_region_out.csv", index=False)
    primary = pearsonr(frame["weighted_z"], frame["validation_score"])
    sn = pd.DataFrame([
        {
            "analysis": "including_substantia_nigra",
            "status": "not_estimable",
            "pearson_r": np.nan,
            "spearman_rho": np.nan,
            "n_regions": len(frame),
            "notes": "ENIGMA-PD 2021 public FreeSurfer summary statistics do not contain substantia nigra; SN was not silently imputed or removed",
        },
        {
            "analysis": "primary_dataset_without_substantia_nigra_measure",
            "status": "reported_primary",
            "pearson_r": primary.statistic,
            "spearman_rho": spearmanr(frame["weighted_z"], frame["validation_score"]).statistic,
            "n_regions": len(frame),
            "notes": "All reliably mapped ENIGMA units included; this is not an SN-exclusion analysis because SN is absent from the phenotype",
        },
    ])
    sn.to_csv(DATA_RESULTS / "parkinson_substantia_nigra_sensitivity.csv", index=False)
    return loo, sn


def geometry_triangles() -> tuple[dict[int, list[np.ndarray]], np.ndarray]:
    payload = json.loads(GEOMETRY_FILE.read_text(encoding="utf-8"))
    parcels: dict[int, list[np.ndarray]] = {}
    all_vertices: list[np.ndarray] = []
    for parcel in payload["regions"]:
        positions = np.asarray(parcel["positions"], dtype=float).reshape(-1, 3)[:, [0, 2, 1]]
        triangles = positions[np.asarray(parcel["indices"], dtype=int).reshape(-1, 3)]
        parcels[int(parcel["region_id"])] = list(triangles)
        all_vertices.append(positions)
    return parcels, np.concatenate(all_vertices)


def draw_brain(axis: plt.Axes, parcels: dict[int, list[np.ndarray]], vertices: np.ndarray, values: dict[int, float], cmap: object, norm: object) -> None:
    triangles: list[np.ndarray] = []
    colors: list[object] = []
    for region_id, parcel_triangles in parcels.items():
        triangles.extend(parcel_triangles)
        color = "#dedede" if region_id not in values or not np.isfinite(values[region_id]) else cmap(norm(values[region_id]))
        colors.extend([color] * len(parcel_triangles))
    axis.add_collection3d(Poly3DCollection(triangles, facecolors=colors, linewidths=0, alpha=0.98))
    axis.auto_scale_xyz(vertices[:, 0], vertices[:, 1], vertices[:, 2])
    axis.set_box_aspect((1.0, 1.15, 1.0), zoom=1.25)
    axis.view_init(elev=18, azim=-48)
    axis.set_axis_off()


def plot_scatter(frame: pd.DataFrame, primary: pd.Series) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    figure, axis = plt.subplots(figsize=(9.2, 7.2))
    palette = {"cortex": "#346b73", "subcortex": "#c46a3b"}
    sns.scatterplot(data=frame, x="weighted_z", y="validation_score", hue="anatomical_class", palette=palette, s=62, ax=axis)
    sns.regplot(data=frame, x="weighted_z", y="validation_score", scatter=False, color="#202020", ci=95, line_kws={"linewidth": 1.8}, ax=axis)
    xz = np.abs((frame["weighted_z"] - frame["weighted_z"].mean()) / frame["weighted_z"].std(ddof=1))
    yz = np.abs((frame["validation_score"] - frame["validation_score"].mean()) / frame["validation_score"].std(ddof=1))
    labels = set(xz.nlargest(3).index) | set(yz.nlargest(3).index)
    for count, index in enumerate(sorted(labels)):
        row = frame.loc[index]
        axis.annotate(row["region"], (row["weighted_z"], row["validation_score"]), xytext=(6, 5 + 8 * (count % 3)), textcoords="offset points", fontsize=7.5)
    axis.axhline(0, color="#999999", linewidth=0.7)
    axis.axvline(0, color="#999999", linewidth=0.7)
    axis.set(
        xlabel="GENE2BRAIN Stage 6 L2G-weighted enrichment Z-score",
        ylabel="Independent Parkinson vulnerability (−Cohen's d)",
        title="GENE2BRAIN vs Independent Parkinson Vulnerability",
    )
    axis.text(
        0.02, 0.98,
        f"Pearson r = {primary['pearson_r']:.3f} (95% bootstrap CI {primary['pearson_ci_low']:.3f}, {primary['pearson_ci_high']:.3f})\n"
        f"Spearman ρ = {primary['spearman_rho']:.3f}; n = {int(primary['n_regions'])} ENIGMA parcels",
        transform=axis.transAxes, va="top", fontsize=10,
    )
    axis.legend(title="ENIGMA measurement", frameon=False)
    figure.tight_layout()
    figure.savefig(FIGURES / "stage_08_gene2brain_vs_independent_validation.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_brain_maps(propagated: pd.DataFrame) -> None:
    parcels, vertices = geometry_triangles()
    discovery = dict(zip(propagated["region_id"], propagated["weighted_z"]))
    validation = dict(zip(propagated["region_id"], propagated["validation_score"]))
    z_bound = max(abs(value) for value in discovery.values())
    v_bound = max(abs(value) for value in validation.values())
    cmap = plt.get_cmap("RdBu_r")
    figure = plt.figure(figsize=(15, 7.2))
    axes = [figure.add_subplot(1, 2, index + 1, projection="3d") for index in range(2)]
    draw_brain(axes[0], parcels, vertices, discovery, cmap, mpl_colors.TwoSlopeNorm(0, -z_bound, z_bound))
    draw_brain(axes[1], parcels, vertices, validation, cmap, mpl_colors.TwoSlopeNorm(0, -v_bound, v_bound))
    axes[0].set_title("GENE2BRAIN enrichment", fontsize=13, pad=2)
    axes[1].set_title("Independent ENIGMA-PD vulnerability", fontsize=13, pad=2)
    figure.colorbar(plt.cm.ScalarMappable(norm=mpl_colors.TwoSlopeNorm(0, -z_bound, z_bound), cmap=cmap), ax=axes[0], shrink=0.54, pad=0.01, label="Stage 6 weighted Z")
    figure.colorbar(plt.cm.ScalarMappable(norm=mpl_colors.TwoSlopeNorm(0, -v_bound, v_bound), cmap=cmap), ax=axes[1], shrink=0.54, pad=0.01, label="−Cohen's d")
    figure.suptitle("Predicted vs Independently Observed Parkinson Regional Vulnerability", fontsize=18, fontweight="bold")
    figure.text(0.5, 0.025, "Gray parcels were not represented by the ENIGMA cortical/subcortical phenotype; identical AAL3 geometry and orientation", ha="center")
    figure.savefig(BRAIN_MAPS / "stage_08_prediction_vs_validation.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    statuses = dict(zip(propagated["region_id"], propagated["agreement_status"]))
    colors = {
        "high_prediction_high_validation": "#7b3294",
        "high_prediction_low_validation": "#008837",
        "low_prediction_high_validation": "#e66101",
        "low_prediction_low_validation": "#5e5e5e",
    }
    figure = plt.figure(figsize=(10.8, 7.8))
    axis = figure.add_subplot(111, projection="3d")
    categorical_values = {region_id: float(list(colors).index(status)) for region_id, status in statuses.items()}
    categorical_cmap = mpl_colors.ListedColormap(list(colors.values()))
    draw_brain(axis, parcels, vertices, categorical_values, categorical_cmap, mpl_colors.BoundaryNorm(np.arange(-0.5, 4.5), 4))
    axis.set_title("Regional Agreement Map", fontsize=18, fontweight="bold", pad=2)
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=10, color=color, label=label.replace("_", " ")) for label, color in colors.items()]
    handles.append(Line2D([0], [0], marker="s", linestyle="", markersize=10, color="#dedede", label="not measured"))
    axis.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02), frameon=False, ncol=2)
    figure.text(0.5, 0.025, "Pre-specified median splits; equality assigned to high; categories are descriptive, not significance calls", ha="center")
    figure.savefig(BRAIN_MAPS / "stage_08_regional_agreement.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_provenance(provenance: dict[str, dict[str, object]], frame: pd.DataFrame, mapping: pd.DataFrame) -> None:
    payload = {
        "dataset": "ENIGMA-Parkinson's 2021 PD versus control summary statistics",
        "citation": "Laansma et al., Movement Disorders 2021, doi:10.1002/mds.28706, PMID:34288137",
        "official_repository": "https://github.com/MICA-MNI/ENIGMA",
        "repository_commit": ENIGMA_COMMIT,
        "access_date": date.today().isoformat(),
        "raw_files": provenance,
        "raw_regions": len(frame),
        "mapped_validation_units": int(mapping["included_in_inference"].sum()),
        "orientation": "validation_score = -d_icv; positive indicates PD-related thinning or smaller subcortical volume",
        "packages": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "nibabel": nib.__version__, "matplotlib": matplotlib.__version__},
    }
    (PROCESSED / "stage_08_validation_provenance.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_reports(frame: pd.DataFrame, results: pd.DataFrame, loo: pd.DataFrame) -> None:
    primary = results.loc[results["discovery_metric"] == "weighted_z_primary"].iloc[0]
    broad = results.loc[results["discovery_metric"] == "broad_z"].iloc[0]
    stringent = results.loc[results["discovery_metric"] == "stringent_z"].iloc[0]
    spatial = results.loc[results["discovery_metric"] == "weighted_z_spatial_moran_null"].iloc[0]
    cortex = results.loc[results["discovery_metric"] == "weighted_z_cortex_only"].iloc[0]
    subcortex = results.loc[results["discovery_metric"] == "weighted_z_subcortex_only"].iloc[0]
    if primary["pearson_r"] > 0 and primary["spearman_rho"] > 0 and (primary["pearson_p"] < 0.05 or primary["spearman_p"] < 0.05):
        label = "PARTIALLY SUPPORTED"
    else:
        label = "NOT SUPPORTED"
    source_review = """# Stage 08 validation-source review

## Decision

ENIGMA-Parkinson's 2021 was selected before any GENE2BRAIN correlation was
calculated. It supplies official machine-readable, hemisphere-specific standardized
PD-versus-control effects for cortical thickness and subcortical volume across 19
sites (2,357 people with PD and 1,182 controls).[^1] The official ENIGMA Toolbox
distributes the exact summary tables used here.[^2]

The selection prioritizes direct quantitative phenotype measurement, broad bilateral
coverage, sample size, standardized effect sizes, explicit FreeSurfer anatomy, and
reproducible numerical access. It does not use the sign or significance of the
GENE2BRAIN association.

## Alternatives and exclusions

Zeighami et al. provide a public whole-brain PPMI map and strong anatomical coverage,
including substantia nigra.[^3] However, the public voxel values are normalized-mixture
ICA component Z-scores, not regional PD-control effects; the paper explicitly makes
that distinction. Its direct regional atrophy t vector is not in the released source
tables. The Frigerio et al. postmortem study is biologically independent and directly
quantifies synaptic loss and pathology in eight cortical regions, but its supplements
do not release reusable region-level numerical observations; digitizing figures was
excluded.[^4] The ENIGMA cerebellar study is large but anatomically restricted,[^5]
and ENIGMA-DTI measures white-matter tracts that cannot be mapped to AAL3 gray-matter
parcels without an additional connectivity model.[^6]

## Independence qualification

The ENIGMA phenotype is independent of AHBA, the Stage 6 permutation model, Stage 7,
and the GENE2BRAIN analysis. It was measured with structural MRI and generated without
GENE2BRAIN outputs. Exact participant-disjointness from the upstream PD GWAS cannot be
guaranteed, because PPMI is one ENIGMA imaging site and PD GWAS consortia have used
PPMI genetic data. This analysis therefore supports an **external measurement and
analysis validation**, not an unqualified participant-independent replication. The
limitation is retained in every interpretation and prevents a `SUPPORTED` label.

## Sources

[^1]: Laansma et al. “[An International Multicenter Analysis of Brain Structure Across Clinical Stages of Parkinson's Disease](https://doi.org/10.1002/mds.28706).” *Movement Disorders* (2021).
[^2]: MICA-MNI. “[ENIGMA Toolbox: load summary statistics](https://enigma-toolbox.readthedocs.io/en/latest/pages/04.loadsumstats/).”
[^3]: Zeighami et al. “[Network structure of brain atrophy in de novo Parkinson's disease](https://doi.org/10.7554/eLife.08440).” *eLife* (2015); [NeuroVault collection 860](https://neurovault.org/collections/860/).
[^4]: Frigerio et al. “[Regional differences in synaptic degeneration](https://doi.org/10.1186/s40478-023-01711-w).” *Acta Neuropathologica Communications* (2024).
[^5]: Kerestes et al. “[Cerebellar Volume and Disease Staging in Parkinson's Disease](https://doi.org/10.1002/mds.29611).” *Movement Disorders* (2023).
[^6]: Owens-Walton et al. “[A worldwide study of white matter microstructural alterations](https://pubmed.ncbi.nlm.nih.gov/39128907/).” (2024).
"""
    (REPORTS / "stage_08_validation_source_review.md").write_text(source_review, encoding="utf-8")

    interpretation = f"""# Stage 08 interpretation

## Discovery and validation are separate

**Discovery** is the frozen Stage 6 L2G-weighted Parkinson gene-set enrichment
Z-score derived from GWAS prioritization and healthy-brain AHBA expression. **Validation**
is the sign-reversed ENIGMA-PD standardized PD-control difference in cortical thickness
or subcortical volume. ENIGMA values never altered the gene set, discovery scores,
permutations, atlas retention, or Stage 7 statistics.

## Primary result

Across {len(frame)} unique mapped ENIGMA parcels, Pearson r was
{primary['pearson_r']:.3f} (95% region-bootstrap CI {primary['pearson_ci_low']:.3f} to
{primary['pearson_ci_high']:.3f}, two-sided p={primary['pearson_p']:.4g}); Spearman rho
was {primary['spearman_rho']:.3f} (95% bootstrap CI {primary['spearman_ci_low']:.3f} to
{primary['spearman_ci_high']:.3f}, p={primary['spearman_p']:.4g}). The spatially
constrained Moran-surrogate p-value was {spatial['pearson_p']:.4g}.

The broad and stringent gene-set Pearson correlations were {broad['pearson_r']:.3f}
and {stringent['pearson_r']:.3f}. Cortex-only and subcortex-only estimates were
{cortex['pearson_r']:.3f} (n={int(cortex['n_regions'])}) and {subcortex['pearson_r']:.3f}
(n={int(subcortex['n_regions'])}), respectively. Leave-one-region-out Pearson r ranged
from {loo['pearson_r'].min():.3f} to {loo['pearson_r'].max():.3f}, with median
{loo['pearson_r'].median():.3f}.

## Substantia nigra

The ENIGMA-PD 2021 public FreeSurfer tables do not measure substantia nigra. An
including-versus-excluding-SN correlation is therefore not estimable for this
phenotype. SN was not imputed, substituted from another outcome, or used to change
the primary analysis.

## Controls and limitations

Region size, baseline AHBA expression, cortical/subcortical class, spatial dependence,
represented-gene count, and validation reliability were evaluated as pre-specified.
The represented-gene count is constant at {int(frame['n_genes_represented'].iloc[0])}
genes in every discovery parcel
and cannot explain between-region covariance. Region-specific analyzed sample size is
used only as a descriptive reliability weight because it is not a measurement-error
variance. AAL3-to-Desikan-Killiany mapping aggregates discovery parcels once per
ENIGMA unit and avoids pseudo-replication, but medium-confidence composite mappings
remain approximate.

The result is classified **{label}** using the complete evidence, not a single
p-value. The validation is independent in measurement and analysis, but exact
participant-disjointness from the upstream GWAS cannot be guaranteed. A positive
association would support convergence between healthy-brain genetic-expression
enrichment and cross-sectional PD morphometric vulnerability. It would not establish
causality, disease origin, temporal direction, cell-type mechanism, or individual
clinical prediction. A null association is retained as informative and does not
trigger discovery-model revision.

## References

1. Laansma et al. “[An International Multicenter Analysis of Brain Structure Across Clinical Stages of Parkinson's Disease](https://doi.org/10.1002/mds.28706).” *Movement Disorders* (2021).
2. MICA-MNI. “[ENIGMA Toolbox summary-statistics documentation](https://enigma-toolbox.readthedocs.io/en/latest/pages/04.loadsumstats/).”
3. Kim et al. “[Multi-ancestry genome-wide association meta-analysis of Parkinson's disease](https://doi.org/10.1038/s41588-023-01584-8).” *Nature Genetics* (2024).
"""
    (REPORTS / "stage_08_interpretation.md").write_text(interpretation, encoding="utf-8")


def final_summary(frame: pd.DataFrame, results: pd.DataFrame, loo: pd.DataFrame) -> None:
    get = lambda metric: results.loc[results["discovery_metric"] == metric].iloc[0]
    primary, broad, stringent, weighted = (get(name) for name in ("weighted_z_primary", "broad_z", "stringent_z", "weighted_z_sensitivity"))
    label = "PARTIALLY SUPPORTED" if primary["pearson_r"] > 0 and primary["spearman_rho"] > 0 and (primary["pearson_p"] < 0.05 or primary["spearman_p"] < 0.05) else "NOT SUPPORTED"
    print("\nINDEPENDENT PARKINSON VALIDATION")
    print("Validation dataset: ENIGMA-Parkinson's 2021 PD-versus-control MRI summary statistics")
    print("Validation cohort: 2,357 PD; 1,182 controls; 19 sites")
    print("Validation phenotype: sign-reversed Cohen's d for cortical thinning/subcortical volume loss")
    print(f"Number of matched regions: {len(frame)} unique ENIGMA parcels")
    print(f"Primary correlation: Pearson r={primary['pearson_r']:.4f}, p={primary['pearson_p']:.4g}; Spearman rho={primary['spearman_rho']:.4f}, p={primary['spearman_p']:.4g}")
    print(f"Broad gene-set correlation: r={broad['pearson_r']:.4f}")
    print(f"Stringent gene-set correlation: r={stringent['pearson_r']:.4f}")
    print(f"Weighted gene-set correlation: r={weighted['pearson_r']:.4f}")
    print("With substantia nigra: not estimable (SN absent from public phenotype)")
    print("Without substantia nigra: not an exclusion analysis; primary dataset contains no SN measure")
    print(f"Leave-one-region-out range: {loo['pearson_r'].min():.4f} to {loo['pearson_r'].max():.4f}; median={loo['pearson_r'].median():.4f}")
    print(f"Interpretation: {label}")
    print("Qualification: external measurement/analysis validation; exact participant-disjointness from the discovery GWAS is not guaranteed")


def main() -> None:
    ensure_directories()
    provenance = download_raw()
    validation = read_validation()
    regions, stage6, stage6_primary, stage7, stage5 = load_discovery()
    mapping = make_mapping(validation, regions)
    frame, propagated = build_analysis_frame(validation, mapping, regions, stage6, stage6_primary, stage7, stage5)
    frame["is_subcortex"] = (frame["anatomical_class"] == "subcortex").astype(int)
    # Save this derived field too; the four required columns and raw values remain present.
    frame.to_csv(PROCESSED / "parkinson_validation_regional_scores.csv", index=False)
    results = calculate_statistics(frame, stage6_primary)
    loo, _ = sensitivity_outputs(frame)
    plot_scatter(frame, results.loc[results["discovery_metric"] == "weighted_z_primary"].iloc[0])
    plot_brain_maps(propagated)
    write_provenance(provenance, validation, mapping)
    write_reports(frame, results, loo)
    final_summary(frame, results, loo)


if __name__ == "__main__":
    main()
