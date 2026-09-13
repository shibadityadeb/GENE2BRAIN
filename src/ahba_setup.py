"""Download, inventory, validate, and visualize the raw AHBA microarray data.

This module intentionally stops before probe reannotation, expression filtering,
normalization, sample-to-region matching, or regional aggregation.
"""

from __future__ import annotations

import csv
import inspect
from datetime import date
from importlib.metadata import version
from pathlib import Path

import abagen
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns
from nilearn import plotting


PROJECT_ROOT = Path(__file__).resolve().parents[1]
AHBA_ROOT = PROJECT_ROOT / "data" / "ahba"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
EXPECTED_DONORS = ("9861", "10021", "12876", "14380", "15496", "15697")
EXPECTED_TYPES = ("microarray", "ontology", "pacall", "probes", "annotation")
ALLEN_DOWNLOAD_BASE = (
    "https://human.brain-map.org/api/v2/well_known_file_download/"
)
PACKAGE_NAMES = (
    "abagen", "pandas", "numpy", "scipy", "nibabel", "nilearn",
    "matplotlib", "seaborn", "statsmodels", "scikit-learn", "requests",
    "jupyter", "setuptools",
)


def show_supported_api() -> str:
    """Print and return the installed fetcher signature."""
    signature = str(inspect.signature(abagen.fetch_microarray))
    print(f"Installed abagen version: {abagen.__version__}")
    print(f"Supported fetch API: abagen.fetch_microarray{signature}")
    return signature


def fetch_ahba(verbose: int = 1) -> dict[str, dict[str, str]]:
    """Fetch all AHBA microarray donors with the installed, inspected API."""
    AHBA_ROOT.mkdir(parents=True, exist_ok=True)
    return abagen.fetch_microarray(
        donors="all",
        data_dir=AHBA_ROOT,
        resume=True,
        verbose=verbose,
        convert=False,
        n_proc=1,
    )


def _count_lines(path: Path, chunk_size: int = 16 * 1024 * 1024) -> int:
    """Count lines without loading a large expression matrix into memory."""
    count = 0
    last = b""
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            count += chunk.count(b"\n")
            last = chunk[-1:]
    return count + (1 if last and last != b"\n" else 0)


def csv_shape(path: Path, has_header: bool = True) -> tuple[int, int]:
    """Return data rows and columns memory-efficiently."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        columns = len(next(csv.reader(stream)))
    return max(_count_lines(path) - int(has_header), 0), columns


def inventory(files: dict[str, dict[str, str]]) -> pd.DataFrame:
    """Inventory every downloaded file, including byte size and dimensions."""
    rows = []
    for donor, donor_files in files.items():
        for file_type in EXPECTED_TYPES:
            path = Path(donor_files[file_type]).resolve()
            exists = path.is_file()
            # Allen's expression and PA-call matrices are intentionally
            # headerless; their first column contains probe IDs.
            has_header = file_type not in {"microarray", "pacall"}
            n_rows, n_columns = csv_shape(path, has_header) if exists else (None, None)
            rows.append(
                {
                    "donor": str(donor),
                    "file_type": file_type,
                    "filename": path.name,
                    "relative_path": path.relative_to(PROJECT_ROOT).as_posix(),
                    "size_bytes": path.stat().st_size if exists else 0,
                    "size_mib": round(path.stat().st_size / 1024**2, 2) if exists else 0,
                    "n_rows": n_rows,
                    "n_columns": n_columns,
                    "present": exists,
                }
            )
    summary = pd.DataFrame(rows).sort_values(["donor", "file_type"])
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(PROCESSED_DIR / "ahba_metadata_summary.csv", index=False)
    return summary


def load_metadata(files: dict[str, dict[str, str]]):
    """Load sample annotations and one canonical copy of probe metadata."""
    annotations = []
    for donor, donor_files in files.items():
        frame = pd.read_csv(donor_files["annotation"])
        frame.insert(0, "donor", str(donor))
        annotations.append(frame)
    samples = pd.concat(annotations, ignore_index=True)
    first_donor = next(iter(files))
    probes = pd.read_csv(files[first_donor]["probes"])
    ontology = pd.read_csv(files[first_donor]["ontology"])
    return samples, probes, ontology


def validate(files, summary, samples, probes) -> list[str]:
    """Perform acquisition-level checks without preprocessing expression."""
    issues: list[str] = []
    donors = tuple(sorted(map(str, files), key=int))
    if donors != tuple(sorted(EXPECTED_DONORS, key=int)):
        issues.append(f"Unexpected donor set: {donors}")
    if set(summary["file_type"]) != set(EXPECTED_TYPES):
        issues.append("One or more expected file types are absent from the inventory")
    if not summary["present"].all() or (summary["size_bytes"] <= 0).any():
        issues.append("One or more expected files are missing or empty")
    coordinate_columns = {"mni_x", "mni_y", "mni_z"}
    if not coordinate_columns.issubset(samples.columns):
        issues.append("MNI coordinate columns are missing from sample annotations")
    elif samples[list(coordinate_columns)].isna().any().any():
        issues.append("Some sample MNI coordinates are missing")
    if "gene_symbol" not in probes.columns:
        issues.append("gene_symbol is missing from Probes.csv")

    for donor, donor_files in files.items():
        donor_samples = samples.loc[samples["donor"] == str(donor)]
        expression_shape = summary.query(
            "donor == @donor and file_type == 'microarray'"
        ).iloc[0]
        if int(expression_shape["n_rows"]) != len(probes):
            issues.append(f"Donor {donor}: expression rows do not match probe rows")
        if int(expression_shape["n_columns"]) - 1 != len(donor_samples):
            issues.append(f"Donor {donor}: expression columns do not match sample rows")
    return issues


def gene_count(probes: pd.DataFrame) -> int | None:
    """Count unique nonblank symbols directly supplied in raw Probes.csv."""
    if "gene_symbol" not in probes:
        return None
    symbols = probes["gene_symbol"].dropna().astype(str).str.strip()
    return int(symbols[symbols.ne("")].nunique())


def make_figures(samples: pd.DataFrame) -> None:
    """Create anatomical sampling coverage and per-donor QC figures."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    donor_order = sorted(samples["donor"].astype(str).unique(), key=int)
    palette = sns.color_palette("colorblind", n_colors=len(donor_order))
    donor_index = {donor: index for index, donor in enumerate(donor_order)}
    values = samples["donor"].astype(str).map(donor_index).to_numpy()
    coordinates = samples[["mni_x", "mni_y", "mni_z"]].to_numpy()
    cmap = ListedColormap(palette)

    figure = plt.figure(figsize=(15, 5.8), facecolor="white")
    display = plotting.plot_markers(
        values,
        coordinates,
        node_size=7,
        node_cmap=cmap,
        node_vmin=-0.5,
        node_vmax=len(donor_order) - 0.5,
        alpha=0.58,
        display_mode="lyrz",
        figure=figure,
        title="AHBA Sampling Coverage Map",
        annotate=True,
        black_bg=False,
        colorbar=False,
    )
    legend = [Patch(color=palette[i], label=f"Donor {d}") for i, d in enumerate(donor_order)]
    figure.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.055),
        ncol=6,
        frameon=False,
        fontsize=9,
    )
    figure.text(
        0.5, 0.012, "Tissue sample locations in MNI space; glass-brain anatomical projections",
        ha="center", fontsize=9, color="#444444",
    )
    figure.savefig(FIGURE_DIR / "ahba_sampling_coverage.png", dpi=300, bbox_inches="tight")
    display.close()

    counts = samples.groupby("donor", sort=False).size().reindex(donor_order)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(counts.index, counts.values, color=palette, edgecolor="white")
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set(title="AHBA Samples by Donor", xlabel="Donor ID", ylabel="Number of tissue samples")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, counts.max() * 1.13)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "ahba_samples_by_donor.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_provenance(files, n_genes: int | None) -> Path:
    """Write a human-readable provenance record with exact package versions."""
    donors = sorted(map(str, files), key=int)
    package_lines = "\n".join(f"  - `{name}=={version(name)}`" for name in PACKAGE_NAMES)
    donor_urls = []
    from abagen.datasets import WELL_KNOWN_IDS
    for donor in donors:
        donor_urls.append(f"  - Donor {donor}: {ALLEN_DOWNLOAD_BASE}{WELL_KNOWN_IDS.url[donor]}")
    content = f"""# AHBA data provenance

- **Source:** Allen Institute for Brain Science, Allen Human Brain Atlas
- **Project page:** https://human.brain-map.org/
- **Dataset name:** Allen Human Brain Atlas normalized microarray expression dataset
- **Access date:** {date.today().isoformat()}
- **Number of donors:** {len(donors)}
- **Donors:** {", ".join(donors)}
- **Data type:** Postmortem adult human brain microarray expression intensities, present/absent calls, probe annotations, sample annotations with MNI coordinates, and anatomical ontology
- **Genes directly represented in `Probes.csv`:** {n_genes if n_genes is not None else "not available"} unique nonblank gene symbols (no reannotation applied)
- **Preprocessing status:** Allen-distributed normalized microarray files only. `convert=False` prevents parquet conversion. No probe reannotation/filtering, sample filtering, normalization by this project, sample-to-region assignment, donor aggregation, or regional aggregation has been performed.
- **Download method:** `abagen.fetch_microarray(donors="all", resume=True, convert=False, n_proc=1)`
- **Download location:** `data/ahba/microarray/`
- **Assumptions:** Donor identifiers and download IDs are those shipped with the installed `abagen` release. `mni_x`, `mni_y`, and `mni_z` in `SampleAnnot.csv` are treated as MNI-space millimetre coordinates. The headerless expression and PA-call matrix columns follow the sample-row order supplied in the same donor bundle; their dimensions are checked against `SampleAnnot.csv`. Gene count describes raw annotation symbols, not a curated analysis-ready gene set.

## Allen well-known-file URLs

{chr(10).join(donor_urls)}

## Exact direct package versions

{package_lines}
"""
    path = PROCESSED_DIR / "ahba_provenance.md"
    path.write_text(content, encoding="utf-8")
    return path


def print_qc_report(files, summary, samples, probes, issues) -> None:
    """Print the concise terminal report requested for this acquisition stage."""
    probes_count = int(summary.loc[summary["file_type"] == "microarray", "n_rows"].max())
    print("\n" + "=" * 68)
    print("AHBA ACQUISITION QUALITY-CONTROL REPORT")
    print("=" * 68)
    print(f"Number of donors:       {len(files)}")
    print(f"Total tissue samples:   {len(samples):,}")
    print(f"Number of probes:       {probes_count:,}")
    genes = gene_count(probes)
    print(f"Genes (raw symbols):    {genes:,}" if genes is not None else "Genes: not directly available")
    print("Major metadata files:   SampleAnnot.csv, Probes.csv, Ontology.csv")
    print(f"Download location:      {AHBA_ROOT / 'microarray'}")
    print(f"Missing/suspicious:     {'; '.join(issues) if issues else 'None detected'}")
    print("Scope check:            No regional aggregation or disease analysis performed")


def run(verbose: int = 1):
    """Run the complete acquisition/QC stage and return its core objects."""
    show_supported_api()
    files = fetch_ahba(verbose=verbose)
    summary = inventory(files)
    samples, probes, ontology = load_metadata(files)
    issues = validate(files, summary, samples, probes)
    make_figures(samples)
    provenance = write_provenance(files, gene_count(probes))
    print(summary.to_string(index=False))
    print(f"\nSample annotations: {samples.shape}; probes: {probes.shape}; ontology: {ontology.shape}")
    print(f"Saved metadata summary: {PROCESSED_DIR / 'ahba_metadata_summary.csv'}")
    print(f"Saved provenance: {provenance}")
    print_qc_report(files, summary, samples, probes, issues)
    if issues:
        raise RuntimeError("AHBA acquisition QC failed: " + "; ".join(issues))
    return files, summary, samples, probes, ontology


if __name__ == "__main__":
    run()
