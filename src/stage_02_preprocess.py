"""Stage 2: build a documented AAL3 regional AHBA expression reference.

Scope is deliberately restricted to healthy AHBA preprocessing. No GWAS,
disease-gene, enrichment, permutation, or vulnerability analysis is performed.
"""

from __future__ import annotations

import inspect
import platform
import warnings
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import abagen
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
AHBA_DIR = ROOT / "data" / "ahba"
ATLAS_DIR = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3"
ATLAS_PATH = ATLAS_DIR / "AAL3v1.nii.gz"
ATLAS_XML = ATLAS_DIR / "AAL3v1.xml"
PROCESSED = ROOT / "data" / "processed"
INTERMEDIATE = ROOT / "data" / "intermediate"
FIGURES = ROOT / "results" / "figures"
REPORTS = ROOT / "reports"

PARAMETERS = {
    "ibf_threshold": 0.5,
    "probe_selection": "diff_stability",
    "donor_probes": "aggregate",
    "sim_threshold": None,
    "lr_mirror": None,
    "missing": None,
    "tolerance": 2,
    "sample_norm": "srs",
    "gene_norm": "srs",
    "norm_matched": True,
    "norm_structures": True,
    "region_agg": "donors",
    "agg_metric": "mean",
    "corrected_mni": True,
    "reannotated": True,
}
MIN_TOTAL_SAMPLES = 2
MIN_DONORS = 2
NEAR_ZERO_VARIANCE = 1e-12
RANDOM_SEED = 42


def ensure_dirs() -> None:
    for path in (PROCESSED, INTERMEDIATE, FIGURES, REPORTS):
        path.mkdir(parents=True, exist_ok=True)


def discover_ahba() -> tuple[dict[str, dict[str, Path]], dict]:
    """Inspect Stage 1 files instead of assuming a donor list or layout."""
    required = {
        "expression": "MicroarrayExpression.csv",
        "probes": "Probes.csv",
        "samples": "SampleAnnot.csv",
        "pacall": "PACall.csv",
        "ontology": "Ontology.csv",
    }
    donor_dirs = sorted(AHBA_DIR.glob("**/normalized_microarray_donor*"))
    files: dict[str, dict[str, Path]] = {}
    for donor_dir in donor_dirs:
        donor = donor_dir.name.removeprefix("normalized_microarray_donor")
        located = {kind: donor_dir / filename for kind, filename in required.items()}
        missing = [str(path) for path in located.values() if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"Donor {donor} is incomplete: {missing}")
        files[donor] = located
    if not files:
        raise FileNotFoundError(f"No AHBA donor directories found under {AHBA_DIR}")

    first = next(iter(files.values()))
    samples = [pd.read_csv(paths["samples"]) for paths in files.values()]
    probes = pd.read_csv(first["probes"])
    ontology = pd.read_csv(first["ontology"])
    audit = {
        "donors": sorted(files, key=int),
        "n_donors": len(files),
        "n_samples": int(sum(len(frame) for frame in samples)),
        "samples_by_donor": {donor: len(frame) for donor, frame in zip(files, samples)},
        "n_probes": len(probes),
        "n_raw_gene_symbols": int(probes["gene_symbol"].dropna().nunique()),
        "n_ontology_structures": len(ontology),
        "sample_fields": list(samples[0].columns),
        "probe_fields": list(probes.columns),
        "ontology_fields": list(ontology.columns),
        "coordinate_fields": [c for c in samples[0].columns if c.startswith(("mni_", "mri_"))],
    }
    return files, audit


def _hemisphere(label: str) -> str:
    if label.endswith("_L"):
        return "L"
    if label.endswith("_R"):
        return "R"
    return "B"


def _structure(label: str) -> str:
    if label.startswith(("Cerebellum", "Vermis")):
        return "cerebellum"
    subcortical = (
        "Caudate", "Putamen", "Pallidum", "Thal_", "N_Acc", "VTA_",
        "SN_", "Red_N", "LC_", "Raphe", "Hippocampus", "Amygdala",
    )
    if label.startswith(subcortical):
        return "subcortex/brainstem"
    return "cortex"


def _system(label: str) -> str:
    if label.startswith(("Cerebellum", "Vermis")):
        return "cerebellum"
    if label.startswith(("VTA_", "SN_", "Red_N", "LC_", "Raphe")):
        return "brainstem"
    if label.startswith(("Hippocampus", "ParaHippocampal")):
        return "hippocampal formation"
    if label.startswith(("Caudate", "Putamen", "Pallidum", "Thal_", "N_Acc", "Amygdala")):
        return "subcortex"
    return "cortex"


def build_atlas_metadata() -> pd.DataFrame:
    """Read actual AAL3 XML labels and derive abagen matching metadata."""
    if not ATLAS_PATH.is_file() or not ATLAS_XML.is_file():
        raise FileNotFoundError("AAL3v2 files are absent; see atlas_selection.md")
    image = nib.load(ATLAS_PATH)
    atlas_values = set(np.unique(np.asanyarray(image.dataobj)).astype(int)) - {0}
    labels = []
    for element in ET.parse(ATLAS_XML).findall("./data/label"):
        region_id = int(element.findtext("index"))
        if region_id in atlas_values:
            labels.append((region_id, element.findtext("name")))
    if {region_id for region_id, _ in labels} != atlas_values:
        raise ValueError("AAL3 XML labels and NIfTI integer values do not match")

    data = []
    atlas_data = np.asanyarray(image.dataobj)
    for region_id, label in labels:
        voxels = np.column_stack(np.where(atlas_data == region_id))
        center_voxel = voxels.mean(axis=0)
        center_mni = nib.affines.apply_affine(image.affine, center_voxel)
        structure = _structure(label)
        system = _system(label)
        data.append({
            "region_id": region_id,
            "region_name": label.replace("_", " "),
            "atlas_label": label,
            "hemisphere": _hemisphere(label),
            "structure": structure,
            "broad_system": system,
            "anatomical_hierarchy": f"brain/{structure}/{system}",
            "centroid_mni_x": center_mni[0],
            "centroid_mni_y": center_mni[1],
            "centroid_mni_z": center_mni[2],
        })
    metadata = pd.DataFrame(data).sort_values("region_id").reset_index(drop=True)
    abagen_info = metadata.rename(columns={"region_id": "id", "atlas_label": "label"})[
        ["id", "label", "hemisphere", "structure"]
    ]
    abagen_info.to_csv(INTERMEDIATE / "aal3v1_abagen_info.csv", index=False)
    abagen.images.check_atlas(ATLAS_PATH, abagen_info)
    return metadata


def _patch_abagen_pandas_compatibility() -> None:
    """Bridge two pandas APIs removed after abagen 0.1.3 was released."""
    if not hasattr(pd.DataFrame, "append"):
        def append(self, other, ignore_index=False, verify_integrity=False, sort=False):
            return pd.concat(
                [self, other], ignore_index=ignore_index,
                verify_integrity=verify_integrity, sort=sort,
            )
        pd.DataFrame.append = append  # type: ignore[attr-defined]
    if "inplace" not in inspect.signature(pd.DataFrame.set_axis).parameters:
        original = pd.DataFrame.set_axis
        def set_axis(self, labels, *, axis=0, inplace=None, copy=None):
            if inplace:
                raise ValueError("In-place set_axis is not supported by this compatibility layer")
            return original(self, labels, axis=axis, copy=copy)
        pd.DataFrame.set_axis = set_axis  # type: ignore[method-assign]


def process_expression(atlas_metadata: pd.DataFrame):
    """Run the inspected abagen 0.1.3 preprocessing workflow."""
    _patch_abagen_pandas_compatibility()
    atlas_info = atlas_metadata.rename(
        columns={"region_id": "id", "atlas_label": "label"}
    )[["id", "label", "hemisphere", "structure"]]
    print("abagen.get_expression_data", inspect.signature(abagen.get_expression_data))
    donor_expression, counts, abagen_report = abagen.get_expression_data(
        ATLAS_PATH,
        atlas_info=atlas_info,
        **PARAMETERS,
        return_counts=True,
        return_donors=True,
        return_report=True,
        donors="all",
        data_dir=AHBA_DIR,
        verbose=1,
        n_proc=1,
    )
    return donor_expression, counts, abagen_report


def filter_and_save(donor_expression, counts, atlas_metadata):
    """Apply explicit coverage rules, save donors, and form donor-equal mean."""
    counts.columns = counts.columns.astype(str)
    counts.index = counts.index.astype(int)
    total = counts.sum(axis=1)
    donor_coverage = counts.gt(0).sum(axis=1)
    keep = total.ge(MIN_TOTAL_SAMPLES) & donor_coverage.ge(MIN_DONORS)
    kept_ids = counts.index[keep]

    donor_paths = []
    aligned = {}
    for donor, matrix in donor_expression.items():
        donor = str(donor)
        matrix.index = matrix.index.astype(int)
        matrix = matrix.loc[kept_ids].sort_index().sort_index(axis=1)
        path = INTERMEDIATE / f"donor_{donor}_expression.csv"
        donor_paths.append(path)
        aligned[donor] = matrix

    stack = np.stack([matrix.to_numpy(dtype=float) for matrix in aligned.values()])
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pooled_values = np.nanmean(stack, axis=0)
    first = next(iter(aligned.values()))
    n_abagen_genes = first.shape[1]
    pooled = pd.DataFrame(pooled_values, index=first.index, columns=first.columns)
    # return_donors=True intentionally preserves normalization failures that
    # abagen would drop from its pooled output. Exclude these genes explicitly
    # from the complete reference rather than imputing them.
    invalid_genes = pooled.columns[~np.isfinite(pooled).all(axis=0)]
    pooled = pooled.drop(columns=invalid_genes)
    aligned = {donor: matrix.loc[:, pooled.columns] for donor, matrix in aligned.items()}
    pooled.index.name = "region_id"
    pooled.to_csv(PROCESSED / "brain_region_gene_expression.csv")
    regional_summary = pd.DataFrame({
        "mean_expression": pooled.mean(axis=1),
        "median_expression": pooled.median(axis=1),
        "std_expression": pooled.std(axis=1, ddof=1),
        "min_expression": pooled.min(axis=1),
        "max_expression": pooled.max(axis=1),
    })
    regional_summary.to_csv(PROCESSED / "regional_expression_summary.csv")

    # Save only after the common reliable gene basis is fixed.
    for donor, matrix in aligned.items():
        matrix.to_csv(INTERMEDIATE / f"donor_{donor}_expression.csv", index_label="region_id")

    region_metadata = atlas_metadata.set_index("region_id").copy()
    region_metadata["total_assigned_samples"] = total
    region_metadata["donors_with_samples"] = donor_coverage
    region_metadata["retained_in_main_matrix"] = keep
    region_metadata["coverage_status"] = np.select(
        [total.eq(0), ~keep],
        ["no_assigned_samples", "insufficient_coverage"],
        default="retained",
    )
    region_metadata.reset_index().to_csv(PROCESSED / "region_metadata.csv", index=False)
    processing_stats = {
        "abagen_genes": n_abagen_genes,
        "invalid_normalized_genes_removed": len(invalid_genes),
    }
    return pooled, aligned, region_metadata.reset_index(), donor_paths, processing_stats


def build_gene_metadata(genes: pd.Index, files) -> pd.DataFrame:
    """Document reannotated mappings and candidate AHBA probes per output gene."""
    from abagen import probes_  # current installed implementation
    raw = pd.read_csv(next(iter(files.values()))["probes"])
    reannotated = probes_.reannotate_probes(raw)
    reannotated = reannotated.reset_index()
    reannotated = reannotated[reannotated["gene_symbol"].isin(genes)]

    def joined(series):
        return ";".join(map(str, pd.unique(series.dropna())))

    grouped = reannotated.groupby("gene_symbol", sort=True).agg(
        entrez_id=("entrez_id", joined),
        ahba_probe_ids=("probe_id", joined),
        ahba_probe_names=("probe_name", joined),
        n_reannotated_candidate_probes=("probe_id", "nunique"),
    ).reset_index()
    grouped["mapping_source"] = "Arnatkeviciute et al. 2019 reannotation distributed with abagen"
    grouped["probe_selection"] = "differential stability after 0.5 intensity filtering"
    grouped["probe_note"] = (
        "Probe lists are reannotated candidates for the retained gene; abagen's public "
        "result exposes gene symbols rather than the selected probe ID."
    )
    grouped.to_csv(PROCESSED / "gene_metadata.csv", index=False)
    return grouped


def donor_correlations(aligned: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Pearson correlations of flattened common region-by-gene values."""
    donors = list(aligned)
    corr = pd.DataFrame(np.eye(len(donors)), index=donors, columns=donors)
    for i, first in enumerate(donors):
        for second in donors[i + 1:]:
            x = aligned[first].to_numpy(dtype=float).ravel()
            y = aligned[second].to_numpy(dtype=float).ravel()
            valid = np.isfinite(x) & np.isfinite(y)
            value = pearsonr(x[valid], y[valid]).statistic if valid.sum() > 2 else np.nan
            corr.loc[first, second] = corr.loc[second, first] = value
    corr.to_csv(PROCESSED / "donor_correlation_matrix.csv", index_label="donor")
    return corr


def key_region_coverage(region_metadata: pd.DataFrame) -> pd.DataFrame:
    definitions = {
        "substantia nigra": ["SN_pc", "SN_pr"],
        "caudate": ["Caudate"],
        "putamen": ["Putamen"],
        "globus pallidus": ["Pallidum"],
        "nucleus accumbens": ["N_Acc"],
        "hippocampus": ["Hippocampus"],
        "amygdala": ["Amygdala"],
        "thalamus": ["Thal_"],
        "cerebellum": ["Cerebellum", "Vermis"],
        "brainstem": ["VTA_", "SN_", "Red_N", "LC_", "Raphe"],
        "frontal cortex": ["Frontal", "OFC", "Rectus", "Precentral", "Supp_Motor"],
        "temporal cortex": ["Temporal", "Heschl", "Fusiform"],
        "parietal cortex": ["Parietal", "Postcentral", "SupraMarginal", "Angular", "Precuneus"],
        "occipital cortex": ["Occipital", "Calcarine", "Cuneus", "Lingual"],
    }
    rows = []
    for structure, patterns in definitions.items():
        matched = region_metadata[
            region_metadata["atlas_label"].apply(lambda value: any(p in value for p in patterns))
        ]
        retained = matched[matched["retained_in_main_matrix"]]
        rows.append({
            "major_structure": structure,
            "present_in_atlas": not matched.empty,
            "atlas_region_count": len(matched),
            "atlas_region_ids": ";".join(map(str, matched["region_id"])),
            "atlas_labels": ";".join(matched["atlas_label"]),
            "assigned_samples": int(matched["total_assigned_samples"].sum()),
            "max_donors_represented": int(matched["donors_with_samples"].max()) if len(matched) else 0,
            "retained_region_count": len(retained),
            "represented_in_main_matrix": not retained.empty,
            "retained_labels": ";".join(retained["atlas_label"]),
        })
    coverage = pd.DataFrame(rows)
    coverage.to_csv(PROCESSED / "key_region_coverage.csv", index=False)
    return coverage


def make_figures(pooled, aligned, region_metadata, correlations):
    """Generate the three prespecified Stage 2 visualizations."""
    metadata = region_metadata.set_index("region_id").loc[pooled.index]
    variances = pooled.var(axis=0, ddof=1).sort_values(ascending=False)
    top_genes = variances.head(50).index
    heat = pooled[top_genes]
    heat_z = (heat - heat.mean()) / heat.std(ddof=1)
    row_colors_map = dict(zip(
        sorted(metadata["broad_system"].unique()),
        sns.color_palette("colorblind", metadata["broad_system"].nunique()),
    ))
    row_colors = metadata["broad_system"].map(row_colors_map)
    row_colors.name = None
    grid = sns.clustermap(
        heat_z, method="average", metric="correlation", cmap="vlag", center=0,
        row_colors=row_colors, xticklabels=True, yticklabels=False,
        figsize=(17, 11), cbar_kws={"label": "Regional z-score"},
    )
    grid.fig.suptitle("Normal Human Brain Transcriptomic Landscape", y=1.02, fontsize=17)
    grid.ax_heatmap.set(xlabel="Top 50 genes by variance across retained regions", ylabel="AAL3 regions (clustered)")
    handles = [plt.Line2D([0], [0], marker="s", color=color, linestyle="", label=name)
               for name, color in row_colors_map.items()]
    grid.ax_heatmap.legend(handles=handles, title="Anatomical system", loc="upper left", bbox_to_anchor=(1.02, 1))
    grid.savefig(FIGURES / "stage_02_normal_transcriptomic_landscape.png", dpi=300, bbox_inches="tight")
    plt.close(grid.fig)

    usable = pooled.loc[:, pooled.var(axis=0, ddof=1) > NEAR_ZERO_VARIANCE]
    scaled = StandardScaler().fit_transform(usable)
    pca = PCA(n_components=2, random_state=RANDOM_SEED)
    scores = pca.fit_transform(scaled)
    plot_data = metadata.reset_index()[["region_id", "atlas_label", "broad_system"]].copy()
    plot_data[["PC1", "PC2"]] = scores
    fig, ax = plt.subplots(figsize=(11, 8))
    sns.scatterplot(
        data=plot_data, x="PC1", y="PC2", hue="broad_system",
        palette=row_colors_map, s=70, alpha=0.85, edgecolor="white", ax=ax,
    )
    ax.axhline(0, color="#dddddd", lw=0.8); ax.axvline(0, color="#dddddd", lw=0.8)
    ax.set_title("Regional Gene Expression Structure")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.legend(title="Anatomical system", bbox_to_anchor=(1.02, 1), loc="upper left")
    sns.despine(); fig.tight_layout()
    fig.savefig(FIGURES / "stage_02_regional_transcriptomic_structure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(correlations, annot=True, fmt=".3f", cmap="viridis", vmin=0, vmax=1,
                square=True, linewidths=0.5, cbar_kws={"label": "Pearson r"}, ax=ax)
    ax.set_title("Donor Concordance")
    ax.set(xlabel="Donor", ylabel="Donor")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_02_donor_concordance.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return list(top_genes), pca.explained_variance_ratio_


def write_atlas_selection(atlas_metadata, coverage):
    image = nib.load(ATLAS_PATH)
    key_status = coverage[["major_structure", "present_in_atlas"]].to_records(index=False)
    key_text = ", ".join(f"{name}: {'yes' if present else 'no'}" for name, present in key_status)
    content = f"""# Atlas selection: AAL3v1 (AAL3v2 distribution)

## Selection

- **Atlas:** Automated Anatomical Labeling atlas 3, AAL3v1 image from the April 2024 AAL3v2 SPM12 distribution
- **Regions:** {len(atlas_metadata)} nonzero labelled regions
- **Source:** GIN / CNRS AAL project; Rolls et al. (2020), *NeuroImage* 206:116189
- **Official archive:** https://www.gin.cnrs.fr/wp-content/uploads/AAL3v2_for_SPM12.tar.gz
- **Article:** https://doi.org/10.1016/j.neuroimage.2019.116189
- **Resolution used:** {image.header.get_zooms()[0]:g} mm isotropic
- **Image dimensions:** {image.shape}
- **Space:** MNI
- **Atlas image:** `data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz`
- **Label definition:** `data/atlases/aal_3v2/AAL3/AAL3v1.xml`
- **AHBA compatibility:** Deterministic integer-labelled volumetric NIfTI in MNI space, validated with `abagen.images.check_atlas`; hemisphere and tissue-class metadata are supplied to constrain sample matching.

## Why AAL3 was selected

AAL3 provides named bilateral cortical regions, basal ganglia, hippocampus,
amygdala, 15 thalamic subdivisions per side, cerebellar lobules/vermis, and
small disease-relevant nuclei including substantia nigra, nucleus accumbens,
VTA, red nucleus, locus coeruleus, and raphe. This coverage is materially
better aligned with a cross-neurological/psychiatric project than a
cortex-dominant parcellation.

The bundled 83-region Desikan–Killiany option was considered because it is
directly distributed by `abagen` and has robust cortical/subcortical regions,
but it lacks explicit cerebellar parcels and substantia nigra. AAL3 was not
selected merely for convenience; its anatomical scope is the deciding factor.

## Important limitation

AAL3's smallest brainstem and thalamic nuclei may be undersampled by the AHBA.
Atlas presence therefore does not imply usable transcriptomic coverage. No
missing parcel is imputed. Regions enter the main reference only when at least
{MIN_TOTAL_SAMPLES} samples from at least {MIN_DONORS} donors are assigned;
the full result is recorded in `key_region_coverage.csv`.

Atlas-level requested-structure presence: {key_text}.
"""
    (PROCESSED / "atlas_selection.md").write_text(content, encoding="utf-8")


def write_reports(audit, atlas_metadata, pooled, aligned, counts, correlations,
                  gene_metadata, coverage, abagen_report, top_genes, pca_ratio,
                  processing_stats):
    # The report returned by abagen 0.1.3 contains trailing spaces and, in one
    # sentence, an embedded backspace character. Sanitize formatting only;
    # preserve the generated scientific description verbatim otherwise.
    abagen_report = "\n".join(
        line.rstrip() for line in abagen_report.replace("\x08", "").splitlines()
    )
    all_genes_after_abagen = pooled.shape[1]
    near_zero = int((pooled.var(axis=0, ddof=1) <= NEAR_ZERO_VARIANCE).sum())
    missing_pct = float(pooled.isna().to_numpy().mean() * 100)
    regional_means = pooled.mean(axis=1)
    removed_regions = int((~atlas_metadata["retained_in_main_matrix"]).sum())
    raw_to_final_removed = audit["n_raw_gene_symbols"] - all_genes_after_abagen
    donor_missing = {d: float(m.isna().to_numpy().mean() * 100) for d, m in aligned.items()}
    offdiag = correlations.to_numpy()[np.triu_indices_from(correlations, 1)]

    qc = f"""# Stage 02 quality-control report

Generated: {date.today().isoformat()}

## Input audit

- AHBA donors: {audit['n_donors']} ({', '.join(audit['donors'])})
- AHBA samples: {audit['n_samples']:,}
- Raw probes: {audit['n_probes']:,}
- Raw annotated gene symbols: {audit['n_raw_gene_symbols']:,}
- AHBA ontology structures: {audit['n_ontology_structures']:,}
- Sample metadata fields: {', '.join(audit['sample_fields'])}
- Probe metadata fields: {', '.join(audit['probe_fields'])}
- Ontology metadata fields: {', '.join(audit['ontology_fields'])}
- Coordinate fields: {', '.join(audit['coordinate_fields'])}

## Output matrix

- Atlas parcels: {len(atlas_metadata)}
- Retained regions: {pooled.shape[0]}
- Regions removed for insufficient coverage: {removed_regions}
- Final complete-case genes in the main matrix: {pooled.shape[1]:,}
- Genes returned by abagen before complete-case QC: {processing_stats['abagen_genes']:,}
- Genes removed because normalized pooled values were non-finite: {processing_stats['invalid_normalized_genes_removed']:,}
- Raw-symbol-to-final gene difference: {raw_to_final_removed:,} (not a one-step filter count; raw annotations and reannotated symbols differ)
- Missing values in main matrix: {missing_pct:.6f}%
- Genes with variance <= {NEAR_ZERO_VARIANCE:g}: {near_zero}
- Coverage rule: >= {MIN_TOTAL_SAMPLES} assigned samples and >= {MIN_DONORS} donors

## Donor robustness

- Samples assigned by donor: {', '.join(f'{d}: {int(counts[d].sum())}' for d in counts)}
- Missing values in retained donor matrices: {', '.join(f'{d}: {v:.2f}%' for d, v in donor_missing.items())}
- Pairwise flattened-expression Pearson r: mean {np.nanmean(offdiag):.3f}, range {np.nanmin(offdiag):.3f}–{np.nanmax(offdiag):.3f}
- Correlations use gene values across regions observed in both donors; missing donor-region combinations are excluded pairwise.

## Distribution and visualization checks

- Main expression range: {np.nanmin(pooled.to_numpy()):.4f} to {np.nanmax(pooled.to_numpy()):.4f}
- Main expression median: {np.nanmedian(pooled.to_numpy()):.4f}
- Mean expression across genes by region: median {regional_means.median():.4f}, IQR {regional_means.quantile(.25):.4f}–{regional_means.quantile(.75):.4f}, range {regional_means.min():.4f}–{regional_means.max():.4f}
- Heatmap genes: top 50 by sample variance across retained regions (not hand-picked): {', '.join(top_genes)}
- PCA used all {pooled.shape[1] - near_zero:,} non-near-zero-variance genes after per-gene standardization; PC1/PC2 explain {pca_ratio[0]*100:.2f}%/{pca_ratio[1]*100:.2f}%.

## Key anatomical coverage

```text
{coverage.to_string(index=False)}
```

No disease, GWAS, enrichment, permutation, or vulnerability analysis was performed.
"""
    (REPORTS / "stage_02_qc.md").write_text(qc, encoding="utf-8")

    methods = f"""# Stage 02 methods

## Inputs and software

- AHBA source: Allen Human Brain Atlas normalized microarray dataset, downloaded in Stage 1 with `abagen.fetch_microarray(donors='all')`.
- Atlas: AAL3v1 from the April 2024 AAL3v2 distribution, 2 mm MNI NIfTI, {len(atlas_metadata)} labelled parcels.
- Python: {platform.python_version()}
- abagen: {abagen.__version__}
- Current inspected API: `abagen.get_expression_data{inspect.signature(abagen.get_expression_data)}`
- Random seed: {RANDOM_SEED} (PCA is deterministic here; recorded for reproducibility).

## Exact preprocessing parameters

```python
{PARAMETERS!r}
```

## Decisions, rationale, and alternatives

1. **Gene reannotation (`reannotated=True`).** Arnatkevičiūtė et al. mappings distributed with `abagen` replace outdated/ambiguous Allen probe annotations and discard probes without reliable mappings. Alternative: raw Allen annotations, which preserve more probes but increase mapping error.
2. **Poor-signal probes (`ibf_threshold=0.5`).** A probe must exceed background in at least 50% of samples across donors. This is the documented `abagen` default and avoids expression dominated by noise. Alternatives include more permissive/strict thresholds or no filtering.
3. **Multiple probes (`probe_selection='diff_stability'`, `donor_probes='aggregate'`).** One probe per gene is selected for the most reproducible regional pattern across donor pairs. Aggregate selection guarantees the same gene/probe basis across donors. Alternatives include maximum intensity/variance, RNA-seq concordance, or averaging probes.
4. **Sample quality (`sim_threshold=None`).** No inter-areal-similarity outlier filter is applied because a defensible threshold was not prespecified; discarded samples would otherwise be difficult to interpret. Alternative: a preregistered similarity threshold.
5. **Coordinates (`corrected_mni=True`).** Corrected coordinates distributed with `abagen` are used. Alternative: original Allen MNI coordinates.
6. **Sample assignment (`tolerance=2`, `exact=None`).** Samples inside a parcel or within 2 mm are assigned to the closest parcel centroid, constrained by atlas hemisphere and broad tissue class. Unmatched samples remain unassigned. Alternatives include exact-only assignment or a larger tolerance, the latter increasing questionable assignments.
7. **Hemisphere (`lr_mirror=None`).** Samples are not mirrored. Four donors have only left-hemisphere sampling, so right-hemisphere evidence is limited to the two bilaterally sampled donors. Mirroring can improve coverage but creates synthetic observations and was rejected for the reference matrix.
8. **Normalization (`sample_norm='srs'`, `gene_norm='srs'`).** Scaled robust sigmoid normalization is performed first across genes within each sample, then for each gene across matched samples separately within donor. `norm_matched=True` avoids unmatched samples affecting scaling. `norm_structures=True` performs gene normalization separately within cortex, subcortex/brainstem, and cerebellum so gross tissue-class differences do not dominate. Alternative global normalization (`norm_structures=False`) preserves class-level offsets but can overwhelm regional effects.
9. **Regional and donor aggregation (`region_agg='donors'`, `agg_metric='mean'`).** Samples are averaged within region separately for each donor. The main matrix is then the equal-weight mean across donors with data, preventing donors with more tissue samples from dominating. Median aggregation is a robust alternative.
10. **Missing data (`missing=None`).** No nearest-centroid or interpolation imputation is used. Donor-region absences remain `NaN`. The main matrix retains parcels with at least {MIN_TOTAL_SAMPLES} assigned samples across at least {MIN_DONORS} donors and averages available donors; this removes {removed_regions} of {len(atlas_metadata)} parcels. Alternatives include retaining single-donor regions or spatial imputation, both rejected to protect robustness.
11. **Gene filtering after abagen.** `abagen` returned {processing_stats['abagen_genes']:,} genes after reannotation, intensity filtering, and probe selection. We removed {processing_stats['invalid_normalized_genes_removed']:,} genes with at least one non-finite pooled regional value rather than imputing them, leaving {pooled.shape[1]:,}. The final matrix includes {near_zero} genes at or below the near-zero variance threshold ({NEAR_ZERO_VARIANCE:g}); these are reported, not silently deleted, and PCA excludes them only to avoid zero-scale features.

## Gene and region metadata

`gene_metadata.csv` records Entrez IDs and reannotated AHBA candidate probe IDs/names for each output symbol. Because the public `abagen` return object exposes genes rather than selected probe IDs, the candidate list is reported transparently rather than inferring an exact selected ID. `region_metadata.csv` retains every AAL3 parcel and marks its sample/donor coverage and inclusion status. `regional_expression_summary.csv` records the across-gene mean, median, standard deviation, minimum, and maximum for every retained region.

## Generated abagen report

{abagen_report}

No downstream disease analysis is part of this stage.
"""
    (REPORTS / "stage_02_methods.md").write_text(methods, encoding="utf-8")


def final_check(audit, atlas_metadata, pooled, gene_metadata, processing_stats):
    required = [
        PROCESSED / "brain_region_gene_expression.csv",
        PROCESSED / "region_metadata.csv",
        PROCESSED / "gene_metadata.csv",
        PROCESSED / "key_region_coverage.csv",
        PROCESSED / "atlas_selection.md",
        REPORTS / "stage_02_qc.md",
        REPORTS / "stage_02_methods.md",
        FIGURES / "stage_02_normal_transcriptomic_landscape.png",
        FIGURES / "stage_02_regional_transcriptomic_structure.png",
        FIGURES / "stage_02_donor_concordance.png",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing Stage 2 outputs: {missing}")
    print("\nSTAGE 02 FINAL SUMMARY")
    print(f"AHBA donors:     {audit['n_donors']}")
    print(f"AHBA samples:    {audit['n_samples']:,}")
    print(f"Brain regions:   {pooled.shape[0]}")
    print(f"Genes:           {pooled.shape[1]:,}")
    print(f"Missing values:  {pooled.isna().to_numpy().mean() * 100:.6f}%")
    print(f"Genes removed:   {audit['n_raw_gene_symbols'] - pooled.shape[1]:,} raw-to-final difference")
    print(f"  post-abagen:    {processing_stats['invalid_normalized_genes_removed']:,} non-finite genes")
    print(f"Regions removed: {len(atlas_metadata) - pooled.shape[0]}")
    print("Atlas:           AAL3v1 (AAL3v2 April 2024 distribution)")
    print(f"Atlas regions:   {len(atlas_metadata)}")
    print("Required files:  all present")
    print("Scope:           healthy reference only; no disease analysis")


def run():
    ensure_dirs()
    files, audit = discover_ahba()
    print("Input audit:", audit)
    atlas_metadata = build_atlas_metadata()
    donor_expression, counts, abagen_report = process_expression(atlas_metadata)
    pooled, aligned, region_metadata, _, processing_stats = filter_and_save(
        donor_expression, counts, atlas_metadata
    )
    gene_metadata = build_gene_metadata(pooled.columns, files)
    correlations = donor_correlations(aligned)
    coverage = key_region_coverage(region_metadata)
    top_genes, pca_ratio = make_figures(pooled, aligned, region_metadata, correlations)
    write_atlas_selection(atlas_metadata, coverage)
    write_reports(
        audit, region_metadata, pooled, aligned, counts, correlations,
        gene_metadata, coverage, abagen_report, top_genes, pca_ratio,
        processing_stats,
    )
    final_check(audit, region_metadata, pooled, gene_metadata, processing_stats)


if __name__ == "__main__":
    run()
