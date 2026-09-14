"""Stage 9: biological interpretation of the frozen Parkinson signal.

This module never modifies or reruns Stages 3--8. It uses the frozen Stage 4
gene sets, Stage 6 enrichment statistics, and Stage 7 spatial ranking as input.
Online functional annotations are cached verbatim for reproducible reruns.
"""

from __future__ import annotations

import hashlib
import gzip
import json
import math
import platform
import sys
import textwrap
import zipfile
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
import nibabel as nib
import numpy as np
import pandas as pd
import requests
import scipy
from scipy.stats import hypergeom, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
GENES = ROOT / "data" / "genes"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "results"
BIOLOGY = ROOT / "data" / "biology"
RAW = BIOLOGY / "raw"
FIGURES = ROOT / "results" / "figures"
REPORTS = ROOT / "reports"
WEB_DATA = ROOT / "web" / "public" / "data"

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
HPA_URL = "https://www.proteinatlas.org/download/tsv/rna_single_nuclei_cluster_type.tsv.zip"
HPA_CLASSES_URL = "https://www.proteinatlas.org/download/tsv/rna_single_nuclei_cluster_type_cluster_types.tsv.zip"
ACCESS_DATE = date.today().isoformat()
TOP_REGIONS = 10
TOP_DRIVERS_PER_REGION = 20
HEATMAP_GENES = 15
HPA_MARKERS_PER_TYPE = 200
MIN_TERM_SIZE = 5
MAX_TERM_SIZE = 1000
MIN_OVERLAP = 1
FDR_THRESHOLD = 0.05


def ensure_directories() -> None:
    for directory in (RAW, RESULTS, FIGURES, REPORTS, WEB_DATA):
        directory.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_inputs() -> dict[str, object]:
    expression = pd.read_csv(PROCESSED / "brain_region_gene_expression.csv", index_col="region_id")
    regions = pd.read_csv(PROCESSED / "region_metadata.csv")
    regions = regions.loc[regions["retained_in_main_matrix"]].set_index("region_id")
    stage6 = pd.read_csv(RESULTS / "parkinson_regional_enrichment.csv").set_index("region_id")
    stage7 = pd.read_csv(RESULTS / "parkinson_spatial_robustness.csv").set_index("region_id")
    validation = pd.read_csv(ROOT / "data" / "validation" / "processed" / "parkinson_validation_aal3_scores.csv")
    gene_sets = {
        name: pd.read_csv(GENES / f"parkinson_genes_{name}.csv")
        for name in ("broad", "stringent", "weighted")
    }
    if not expression.index.equals(regions.index):
        expression = expression.reindex(regions.index)
    if expression.isna().any().any() or len(expression) != 138:
        raise RuntimeError("Frozen AHBA matrix failed Stage 9 input audit")
    if not gene_sets["broad"]["gene"].equals(gene_sets["weighted"]["gene"]):
        raise RuntimeError("Stage 4 broad and weighted membership unexpectedly differ")
    return {
        "expression": expression,
        "regions": regions,
        "stage6": stage6,
        "stage7": stage7,
        "validation": validation,
        "gene_sets": gene_sets,
    }


def write_region_rule(stage7: pd.DataFrame) -> pd.DataFrame:
    selected = stage7.sort_values("robustness_rank").head(TOP_REGIONS).copy()
    table_frame = selected.reset_index()[[
        "robustness_rank", "region_id", "region_name", "stage6_z", "stage6_fdr",
        "spatial_null_fdr", "spatial_robustness_label",
    ]]
    header = "| " + " | ".join(table_frame.columns) + " |"
    separator = "| " + " | ".join("---" for _ in table_frame.columns) + " |"
    body = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in table_frame.itertuples(index=False, name=None)
    ]
    markdown_table = "\n".join([header, separator, *body])
    text = f"""# Stage 9 region-selection rule

## Frozen rule

Before any Stage 9 pathway or cell-type result was generated, the regional
interpretation set was fixed as the **top {TOP_REGIONS} AAL3 parcels by the
Stage 7 `robustness_rank` column**. Ties would be resolved by ascending
`region_id`. No anatomical name, Parkinson relevance, Stage 9 pathway result,
or Stage 8 validation value enters selection.

The Stage 7 result itself is not re-estimated. A high rank does not imply that
a parcel passed joint FDR; the original `spatial_robustness_label` is retained
and reported. This is therefore a ranked exploratory interpretation of frozen
regions, not a new discovery threshold.

Within each selected parcel, the top {TOP_DRIVERS_PER_REGION} genes are defined
by descending additive contribution `expression × L2G / sum(L2G)`. The heatmap
uses the {HEATMAP_GENES} genes with the highest mean contribution across these
parcels. These constants are set in `src/stage_09_biological_interpretation.py`.

## Selected regions

{markdown_table}
"""
    (REPORTS / "stage_09_region_selection.md").write_text(text)
    return selected


def gene_id_mapping(gene_sets: dict[str, pd.DataFrame], background: set[str]) -> pd.DataFrame:
    memberships: dict[str, list[str]] = {}
    for set_name, frame in gene_sets.items():
        for gene in frame["gene"]:
            memberships.setdefault(gene, []).append(set_name)
    broad = gene_sets["broad"].drop_duplicates("gene").copy()
    rows = []
    for row in broad.itertuples(index=False):
        ensembl = str(row.ensembl_id).split(".")[0] if pd.notna(row.ensembl_id) else ""
        if not ensembl:
            status = "unmapped_identifier"
        elif row.gene not in background:
            status = "mapped_not_represented_in_ahba"
        else:
            status = "mapped_and_represented_in_ahba"
        rows.append({
            "original_id": row.gene,
            "gene_symbol": row.gene,
            "ensembl_id": ensembl,
            "mapping_status": status,
            "source": f"Stage 4 Open Targets L2G; sets={';'.join(memberships[row.gene])}",
        })
    mapping = pd.DataFrame(rows)
    mapping.to_csv(RESULTS / "stage_09_gene_id_mapping.csv", index=False)
    return mapping


def request_gprofiler(
    cache_name: str,
    query: list[str],
    background: list[str],
    sources: list[str],
    ordered: bool = False,
) -> dict[str, object]:
    cache = RAW / cache_name
    if cache.exists():
        with gzip.open(cache, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    payload = {
        "organism": "hsapiens",
        "query": query,
        "sources": sources,
        "ordered": ordered,
        "all_results": True,
        "user_threshold": 1.0,
        "domain_scope": "custom",
        "background": background,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
    }
    response = requests.post(
        GPROFILER_URL,
        json=payload,
        headers={"User-Agent": "GENE2BRAIN-stage09/1.0"},
        timeout=180,
    )
    response.raise_for_status()
    data = response.json()
    with gzip.open(cache, "wt", encoding="utf-8") as stream:
        json.dump({"request": payload, "response": data}, stream, separators=(",", ":"))
    return {"request": payload, "response": data}


def extract_intersection(
    wrapped: dict[str, object], item: dict[str, object], query_labels: dict[str, str]
) -> list[str]:
    intersections = item.get("intersections", [])
    query_meta = wrapped["response"].get("meta", {}).get("genes_metadata", {}).get("query", {}).get("query_1", {})
    ensgs = query_meta.get("ensgs", [])
    mapping = query_meta.get("mapping", {})
    if len(intersections) != len(ensgs):
        return []
    ensg_to_original: dict[str, str] = {}
    for original, identifiers in mapping.items():
        for identifier in identifiers:
            ensg_to_original.setdefault(identifier, original)
    genes = []
    for ensg, evidence in zip(ensgs, intersections):
        if evidence:
            original = ensg_to_original.get(ensg, ensg)
            genes.append(query_labels.get(original, query_labels.get(ensg, original)))
    return list(dict.fromkeys(genes))


def parse_ora(
    wrapped: dict[str, object],
    query: list[str],
    gene_set: str,
    query_labels: dict[str, str] | None = None,
) -> pd.DataFrame:
    query_labels = query_labels or {value: value for value in query}
    records = []
    for item in wrapped["response"]["result"]:
        n = int(item["query_size"])
        k = int(item["intersection_size"])
        N = int(item["effective_domain_size"])
        K = int(item["term_size"])
        if K < MIN_TERM_SIZE or K > MAX_TERM_SIZE or k < MIN_OVERLAP or not all((n, N, K)):
            continue
        raw_p = float(hypergeom.sf(k - 1, N, K, n))
        fold = (k / n) / (K / N)
        records.append({
            "term_id": str(item["native"]).replace("REAC:", "").replace("WP:", ""),
            "term_name": item["name"],
            "ontology": item["source"],
            "overlap": k,
            "gene_count": n,
            "background_count": K,
            "effect_size": fold,
            "p_value": raw_p,
            "genes": ";".join(extract_intersection(wrapped, item, query_labels)),
            "gene_set": gene_set,
            "query_size": n,
            "effective_domain_size": N,
            "service_adjusted_p": float(item["p_value"]),
        })
    frame = pd.DataFrame(records)
    if frame.empty:
        return frame
    frame["fdr"] = np.nan
    for _, index in frame.groupby("ontology").groups.items():
        frame.loc[index, "fdr"] = multipletests(frame.loc[index, "p_value"], method="fdr_bh")[1]
    return frame.sort_values(["ontology", "fdr", "p_value", "effect_size"], ascending=[True, True, True, False])


def functional_enrichment(
    gene_sets: dict[str, pd.DataFrame], background: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    outputs: dict[str, pd.DataFrame] = {}
    versions: dict[str, object] = {}
    for set_name in ("weighted", "stringent"):
        frame = gene_sets[set_name]
        represented = frame.loc[frame["gene"].isin(background)].drop_duplicates("gene")
        query = represented["ensembl_id"].str.split(".").str[0].tolist()
        query_labels = dict(zip(query, represented["gene"]))
        wrapped = request_gprofiler(
            f"gprofiler_{set_name}_ora_ensembl.json.gz", query, background,
            ["GO:BP", "GO:MF", "GO:CC", "REAC", "WP"],
        )
        outputs[set_name] = parse_ora(wrapped, query, set_name, query_labels)
        versions[set_name] = wrapped["response"].get("meta", {})
    outputs["broad"] = outputs["weighted"].assign(gene_set="broad")

    combined = pd.concat(outputs.values(), ignore_index=True)
    go = combined.loc[combined["ontology"].isin(["GO:BP", "GO:MF", "GO:CC"])].copy()
    reactome = combined.loc[combined["ontology"] == "REAC"].copy()
    go = go.loc[go["gene_set"] == "weighted", [
        "term_id", "term_name", "ontology", "overlap", "gene_count",
        "background_count", "effect_size", "p_value", "fdr", "genes",
    ]]
    reactome_primary = reactome.loc[reactome["gene_set"] == "weighted", [
        "term_id", "term_name", "gene_count", "overlap", "p_value", "fdr",
        "genes", "effect_size",
    ]]
    go.to_csv(RESULTS / "stage_09_go_enrichment.csv", index=False)
    reactome_primary.to_csv(RESULTS / "stage_09_reactome_enrichment.csv", index=False)

    pathway_sensitivity = combined.loc[combined["ontology"].isin(["REAC", "WP"]), [
        "gene_set", "ontology", "term_id", "term_name", "overlap", "gene_count",
        "background_count", "effect_size", "p_value", "fdr", "genes",
    ]].sort_values(["gene_set", "ontology", "fdr", "p_value"])
    pathway_sensitivity.to_csv(RESULTS / "stage_09_pathway_sensitivity.csv", index=False)
    return go, reactome_primary, pathway_sensitivity, versions


def ranked_pathway_analysis(
    gene_sets: dict[str, pd.DataFrame], pathway_sensitivity: pd.DataFrame
) -> pd.DataFrame:
    weighted = gene_sets["weighted"].drop_duplicates("gene").set_index("gene")
    weights = weighted["gene_weight"].astype(float)
    rows = []
    primary = pathway_sensitivity.loc[pathway_sensitivity["gene_set"] == "weighted"]
    for item in primary.itertuples(index=False):
        members = [gene for gene in str(item.genes).split(";") if gene in weights.index]
        inside = weights.loc[members].to_numpy()
        outside = weights.loc[~weights.index.isin(members)].to_numpy()
        if len(inside) < 2 or len(outside) == 0:
            continue
        statistic, p_value = mannwhitneyu(inside, outside, alternative="greater")
        auc = float(statistic / (len(inside) * len(outside)))
        rows.append({
            "pathway_id": item.term_id,
            "pathway_name": item.term_name,
            "database": item.ontology,
            "overlap": len(inside),
            "rank_effect_auc": auc,
            "median_l2g_in_pathway": float(np.median(inside)),
            "median_l2g_outside_pathway": float(np.median(outside)),
            "p_value": float(p_value),
            "genes": ";".join(members),
            "ranking_interpretation": "higher values mean stronger Open Targets L2G prioritization evidence, not expression effect",
        })
    ranked = pd.DataFrame(rows)
    if not ranked.empty:
        ranked["fdr"] = np.nan
        for _, index in ranked.groupby("database").groups.items():
            ranked.loc[index, "fdr"] = multipletests(ranked.loc[index, "p_value"], method="fdr_bh")[1]
        ranked = ranked.sort_values(["database", "fdr", "p_value"])
    ranked.to_csv(RESULTS / "stage_09_ranked_pathway_analysis.csv", index=False)
    return ranked


def hpa_marker_sets(background: set[str]) -> tuple[dict[str, set[str]], pd.DataFrame]:
    expression_file = RAW / "rna_single_nuclei_cluster_type.tsv.zip"
    class_file = RAW / "rna_single_nuclei_cluster_type_cluster_types.tsv.zip"
    if not expression_file.exists() or not class_file.exists():
        raise FileNotFoundError("Download the official HPA single-nucleus brain files before Stage 9")
    hpa = pd.read_csv(expression_file, sep="\t", compression="zip")
    classes = pd.read_csv(class_file, sep="\t", compression="zip")
    # HPA is keyed by stable Ensembl ID; a small number of current records share
    # a display symbol. AHBA is symbol-indexed, so collapse those deterministically
    # with the maximum observed nCPM rather than duplicating a background gene.
    matrix = hpa.pivot_table(
        index="Gene name", columns="Cluster type", values="nCPM", aggfunc="max"
    )
    mean_other = (matrix.sum(axis=1).to_numpy()[:, None] - matrix.to_numpy()) / (matrix.shape[1] - 1)
    specificity = np.log2((matrix.to_numpy() + 1.0) / (mean_other + 1.0))
    marker_sets: dict[str, set[str]] = {}
    audit_rows = []
    class_map = classes.set_index("Cluster type")["Cell type class"].to_dict()
    for column_index, cell_type in enumerate(matrix.columns):
        candidates = pd.DataFrame({
            "gene": matrix.index,
            "expression": matrix.iloc[:, column_index].to_numpy(),
            "specificity": specificity[:, column_index],
        })
        candidates = candidates.loc[
            (candidates["gene"].isin(background))
            & (candidates["expression"] >= 1.0)
            & (candidates["specificity"] > 1.0)
        ].sort_values(["specificity", "expression", "gene"], ascending=[False, False, True])
        selected = candidates.head(HPA_MARKERS_PER_TYPE)
        marker_sets[cell_type] = set(selected["gene"])
        for row in selected.itertuples(index=False):
            audit_rows.append({
                "cell_type": cell_type,
                "cell_type_class": class_map.get(cell_type, ""),
                "gene": row.gene,
                "nCPM": row.expression,
                "specificity_log2_ratio": row.specificity,
                "marker_rank": len([x for x in audit_rows if x["cell_type"] == cell_type]) + 1,
            })
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(BIOLOGY / "hpa_brain_cell_type_markers.csv", index=False)
    return marker_sets, classes


def cell_type_enrichment(
    gene_sets: dict[str, pd.DataFrame], background: set[str], marker_sets: dict[str, set[str]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_rows = []
    for set_name in ("weighted", "broad", "stringent"):
        query = set(gene_sets[set_name]["gene"]) & background
        for cell_type, markers in marker_sets.items():
            overlap = sorted(query & markers)
            p_value = hypergeom.sf(len(overlap) - 1, len(background), len(markers), len(query))
            effect = (len(overlap) / len(query)) / (len(markers) / len(background))
            all_rows.append({
                "gene_set": set_name,
                "cell_type": cell_type,
                "marker_source": f"Human Protein Atlas v25.1 single-nucleus brain; top {HPA_MARKERS_PER_TYPE} specificity-ranked AHBA genes",
                "overlap": len(overlap),
                "gene_count": len(query),
                "background_count": len(markers),
                "effect_size": effect,
                "p_value": p_value,
                "genes": ";".join(overlap),
            })
    sensitivity = pd.DataFrame(all_rows)
    sensitivity["fdr"] = np.nan
    for _, index in sensitivity.groupby("gene_set").groups.items():
        sensitivity.loc[index, "fdr"] = multipletests(sensitivity.loc[index, "p_value"], method="fdr_bh")[1]
    sensitivity = sensitivity.sort_values(["gene_set", "fdr", "p_value"])
    primary = sensitivity.loc[sensitivity["gene_set"] == "weighted", [
        "cell_type", "marker_source", "overlap", "gene_count", "effect_size",
        "p_value", "fdr", "genes",
    ]]
    primary.to_csv(RESULTS / "stage_09_cell_type_enrichment.csv", index=False)
    sensitivity.to_csv(RESULTS / "stage_09_cell_type_sensitivity.csv", index=False)
    return primary, sensitivity


def regional_drivers(
    expression: pd.DataFrame,
    weighted_genes: pd.DataFrame,
    selected_regions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    weights = weighted_genes.drop_duplicates("gene").set_index("gene")["gene_weight"]
    represented = weights.index.intersection(expression.columns)
    weights = weights.loc[represented].astype(float)
    rows = []
    for region_id, region in selected_regions.sort_values("robustness_rank").iterrows():
        values = expression.loc[region_id, represented].astype(float)
        additive = values * weights / weights.sum()
        ranks = additive.rank(method="min", ascending=False).astype(int)
        for gene in represented:
            rows.append({
                "region_id": region_id,
                "region": region.region_name,
                "robustness_rank": int(region.robustness_rank),
                "spatial_robustness_label": region.spatial_robustness_label,
                "gene": gene,
                "gene_expression": values[gene],
                "l2g_score": weights[gene],
                "weighted_contribution": additive[gene],
                "regional_rank": ranks[gene],
                "top_driver": bool(ranks[gene] <= TOP_DRIVERS_PER_REGION),
            })
    drivers = pd.DataFrame(rows).sort_values(["robustness_rank", "regional_rank", "gene"])
    drivers.to_csv(RESULTS / "stage_09_regional_gene_drivers.csv", index=False)
    top = drivers.loc[drivers["top_driver"]]
    summary = top.groupby("gene").agg(
        l2g_score=("l2g_score", "first"),
        number_of_high_rank_regions=("region", "nunique"),
        mean_expression=("gene_expression", "mean"),
        mean_weighted_contribution=("weighted_contribution", "mean"),
        maximum_contribution=("weighted_contribution", "max"),
        regions=("region", lambda values: ";".join(sorted(set(values)))),
    ).reset_index().sort_values(
        ["number_of_high_rank_regions", "mean_weighted_contribution", "l2g_score"],
        ascending=[False, False, False],
    )
    summary.to_csv(RESULTS / "stage_09_driver_genes.csv", index=False)
    return drivers, summary


def region_pathways(
    drivers: pd.DataFrame, background: list[str], selected_regions: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for region_id, group in drivers.loc[drivers["top_driver"]].groupby("region_id", sort=False):
        query = group.sort_values("regional_rank")["gene"].tolist()
        if len(query) < 10:
            rows.append({
                "region_id": region_id, "region": group["region"].iloc[0],
                "status": "insufficient genes for reliable pathway enrichment",
            })
            continue
        wrapped = request_gprofiler(
            f"gprofiler_region_{int(region_id)}_reactome.json.gz", query, background, ["REAC"]
        )
        parsed = parse_ora(wrapped, query, str(region_id))
        if parsed.empty:
            rows.append({
                "region_id": region_id, "region": group["region"].iloc[0],
                "status": "no eligible Reactome terms",
            })
            continue
        for item in parsed.itertuples(index=False):
            rows.append({
                "region_id": region_id,
                "region": group["region"].iloc[0],
                "pathway_id": item.term_id,
                "pathway": item.term_name,
                "overlap": item.overlap,
                "effect_size": item.effect_size,
                "p_value": item.p_value,
                "within_region_fdr": item.fdr,
                "genes": item.genes,
                "status": "tested",
            })
    regional = pd.DataFrame(rows)
    tested = regional["status"].eq("tested")
    regional["regional_fdr"] = np.nan
    if tested.any():
        regional.loc[tested, "regional_fdr"] = multipletests(
            regional.loc[tested, "p_value"], method="fdr_bh"
        )[1]
    regional.to_csv(RESULTS / "stage_09_region_pathway_enrichment.csv", index=False)
    return regional


def biological_evidence_summary(
    selected: pd.DataFrame,
    stage6: pd.DataFrame,
    regional_pathways: pd.DataFrame,
    cell_types: pd.DataFrame,
    validation: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        "region", "region_enrichment", "fdr", "gene_set", "pathway",
        "pathway_fdr", "driver_genes", "cell_type", "cell_type_fdr",
        "validation_status",
    ]
    supported = regional_pathways.loc[
        regional_pathways.get("regional_fdr", pd.Series(dtype=float)).fillna(1) < FDR_THRESHOLD
    ]
    significant_cells = cell_types.loc[cell_types["fdr"] < FDR_THRESHOLD]
    validation_map = {}
    if "region_id" in validation.columns:
        for row in validation.itertuples(index=False):
            score = getattr(row, "validation_score", np.nan)
            validation_map[int(row.region_id)] = "measured; Stage 8 validation not supported" if pd.notna(score) else "not measured"
    rows = []
    for item in supported.itertuples(index=False):
        driver_set = set(str(item.genes).split(";"))
        relevant = significant_cells.loc[
            significant_cells["genes"].fillna("").map(lambda value: bool(driver_set & set(value.split(";"))))
        ]
        cell = relevant.iloc[0] if not relevant.empty else None
        rows.append({
            "region": item.region,
            "region_enrichment": stage6.loc[item.region_id, "z_score"],
            "fdr": stage6.loc[item.region_id, "fdr_p"],
            "gene_set": "weighted",
            "pathway": item.pathway,
            "pathway_fdr": item.regional_fdr,
            "driver_genes": item.genes,
            "cell_type": cell.cell_type if cell is not None else np.nan,
            "cell_type_fdr": cell.fdr if cell is not None else np.nan,
            "validation_status": validation_map.get(int(item.region_id), "not measured"),
        })
    summary = pd.DataFrame(rows, columns=columns)
    summary.to_csv(RESULTS / "stage_09_biological_evidence_summary.csv", index=False)
    return summary


def plot_pathways(reactome: pd.DataFrame) -> None:
    top = reactome.sort_values(["fdr", "p_value", "effect_size"], ascending=[True, True, False]).head(15).copy()
    top = top.sort_values("effect_size")
    fig, ax = plt.subplots(figsize=(10, 7))
    if top.empty:
        ax.text(0.5, 0.5, "No eligible Reactome terms", ha="center", va="center")
    else:
        size = 35 + top["overlap"] * 24
        color_values = -np.log10(top["fdr"].clip(lower=1e-300))
        scatter = ax.scatter(top["effect_size"], range(len(top)), s=size, c=color_values, cmap="viridis", edgecolor="black", linewidth=.4)
        labels = ["\n".join(textwrap.wrap(value, width=48)) for value in top["term_name"]]
        ax.set_yticks(range(len(top)), labels=labels)
        ax.axvline(1, color="#888", linestyle="--", linewidth=1)
        ax.set_xlabel("Fold enrichment")
        ax.set_ylabel("Reactome pathway")
        fig.colorbar(scatter, ax=ax, label="−log10(BH FDR)")
    status = "No primary Reactome pathway passed BH FDR < 0.05"
    ax.set_title(f"Biological Programs Behind Parkinson Genetic Risk\nTop tested terms; {status}", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_09_parkinson_pathway_enrichment.png", dpi=250)
    plt.close(fig)


def plot_cell_types(cell_types: pd.DataFrame) -> None:
    top = cell_types.sort_values(["fdr", "p_value", "effect_size"], ascending=[True, True, False]).head(20).sort_values("effect_size")
    fig, ax = plt.subplots(figsize=(10, 8))
    colors_ = np.where(top["fdr"] < FDR_THRESHOLD, "#b2182b", "#6b7280")
    bars = ax.barh(top["cell_type"], top["effect_size"], color=colors_)
    ax.axvline(1, color="#222", linestyle="--", linewidth=1)
    for bar, row in zip(bars, top.itertuples(index=False)):
        ax.text(bar.get_width() + .04, bar.get_y() + bar.get_height()/2, f"q={row.fdr:.3g} · n={row.overlap}", va="center", fontsize=8)
    ax.set_xlabel("Fold enrichment")
    ax.set_title("Parkinson Cell-Type Enrichment\nHPA v25.1 human-brain single-nucleus markers; none passed BH FDR < 0.05")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_09_parkinson_cell_type_enrichment.png", dpi=250)
    plt.close(fig)


def plot_drivers(drivers: pd.DataFrame) -> None:
    gene_order = (
        drivers.groupby("gene")["weighted_contribution"].mean().nlargest(HEATMAP_GENES).index.tolist()
    )
    matrix = drivers.loc[drivers["gene"].isin(gene_order)].pivot(
        index="region", columns="gene", values="weighted_contribution"
    )
    region_order = drivers[["region", "robustness_rank"]].drop_duplicates().sort_values("robustness_rank")["region"]
    matrix = matrix.reindex(index=region_order, columns=gene_order)
    standardized = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(12, 7))
    image = ax.imshow(standardized, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(gene_order)), labels=gene_order, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix)), labels=matrix.index)
    ax.set_title("Regional Gene Drivers\nRow-standardized weighted contribution; frozen Stage 7 top 10 regions")
    fig.colorbar(image, ax=ax, label="Within-region z-score")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_09_regional_gene_drivers.png", dpi=250)
    plt.close(fig)


def plot_sensitivity(pathways: pd.DataFrame) -> None:
    reactome = pathways.loc[pathways["ontology"] == "REAC"].copy()
    terms = reactome.groupby("term_name")["fdr"].min().nsmallest(15).index
    matrix = reactome.loc[reactome["term_name"].isin(terms)].pivot_table(index="term_name", columns="gene_set", values="fdr", aggfunc="min")
    matrix = matrix.reindex(columns=["broad", "stringent", "weighted"])
    matrix = matrix.loc[matrix.min(axis=1).sort_values().index]
    values = -np.log10(matrix.clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(9, 7))
    image = ax.imshow(values, aspect="auto", cmap="magma")
    ax.set_xticks(range(3), labels=matrix.columns)
    ax.set_yticks(range(len(matrix)), labels=matrix.index)
    ax.set_title("Pathway Sensitivity Across Frozen Gene-Set Definitions")
    fig.colorbar(image, ax=ax, label="−log10(BH FDR)")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_09_pathway_sensitivity.png", dpi=250)
    plt.close(fig)


def plot_genetic_region_biology(
    drivers: pd.DataFrame, driver_summary: pd.DataFrame, regional: pd.DataFrame,
    weighted_genes: pd.DataFrame,
) -> None:
    top_genes = driver_summary.head(8)["gene"].tolist()
    loci = weighted_genes.drop_duplicates("gene").set_index("gene").reindex(top_genes)["top_locus"].fillna("unmapped locus")
    top_regions = drivers[["region", "robustness_rank"]].drop_duplicates().sort_values("robustness_rank").head(5)["region"].tolist()
    supported = regional.loc[regional.get("regional_fdr", pd.Series(dtype=float)).fillna(1) < FDR_THRESHOLD]
    pathways = supported["pathway"].drop_duplicates().head(6).tolist()
    columns = [list(dict.fromkeys(loci.tolist())), top_genes, top_regions, pathways]
    x_positions = [0, 1, 2, 3]
    positions = {}
    fig, ax = plt.subplots(figsize=(15, 9))
    for x, labels in zip(x_positions, columns):
        if not labels:
            continue
        ys = np.linspace(.9, .1, len(labels))
        for y, label in zip(ys, labels):
            positions[(x, label)] = (x, y)
            ax.text(x, y, label, ha="center", va="center", fontsize=8,
                    bbox={"boxstyle": "round,pad=.3", "fc": "white", "ec": "#555"})
    for gene in top_genes:
        locus = loci.get(gene)
        if (0, locus) in positions:
            ax.plot([positions[(0, locus)][0], positions[(1, gene)][0]], [positions[(0, locus)][1], positions[(1, gene)][1]], color="#9ca3af", lw=.8, zorder=-1)
        gene_regions = drivers.loc[(drivers["gene"] == gene) & drivers["top_driver"] & drivers["region"].isin(top_regions), "region"]
        for region in gene_regions:
            ax.plot([1, 2], [positions[(1, gene)][1], positions[(2, region)][1]], color="#4c78a8", alpha=.55, lw=1)
    for row in supported.loc[supported["region"].isin(top_regions) & supported["pathway"].isin(pathways)].itertuples(index=False):
        ax.plot([2, 3], [positions[(2, row.region)][1], positions[(3, row.pathway)][1]], color="#b2182b", alpha=.65, lw=1)
    if not pathways:
        ax.text(3, .5, "No region-specific Reactome\nrelationship survived BH FDR", ha="center", va="center", color="#555")
    for x, heading in zip(x_positions, ["Genetic loci", "Prioritized genes", "Frozen top regions", "Supported pathways"]):
        ax.text(x, 1.02, heading, ha="center", weight="bold")
    ax.set_xlim(-.45, 3.45); ax.set_ylim(0, 1.08); ax.axis("off")
    ax.set_title("Genetic Risk → Region → Biology\nOnly observed locus–gene, top-contributor, and FDR-supported region–pathway links", pad=20)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_09_genetic_region_biology.png", dpi=250)
    plt.close(fig)


def plot_summary(
    regions: pd.DataFrame, stage6: pd.DataFrame, reactome: pd.DataFrame, cell_types: pd.DataFrame
) -> None:
    fig = plt.figure(figsize=(18, 6))
    ax1 = fig.add_subplot(131, projection="3d")
    coords = regions[["centroid_mni_x", "centroid_mni_y", "centroid_mni_z"]].to_numpy(float)
    scores = stage6.reindex(regions.index)["z_score"].to_numpy()
    norm = colors.TwoSlopeNorm(vmin=min(scores.min(), -1), vcenter=0, vmax=max(scores.max(), 1))
    points = ax1.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=scores, cmap="RdBu_r", norm=norm, s=28, alpha=.9)
    ax1.set_title("A · GENE2BRAIN enrichment")
    ax1.set_xlabel("MNI x"); ax1.set_ylabel("MNI y"); ax1.set_zlabel("MNI z")
    fig.colorbar(points, ax=ax1, fraction=.03, pad=.08, label="Stage 6 Z")

    ax2 = fig.add_subplot(132)
    pathways = reactome.sort_values(["fdr", "p_value"]).head(10).sort_values("effect_size")
    ax2.barh(pathways["term_name"], pathways["effect_size"], color="#4c78a8")
    ax2.axvline(1, color="#333", ls="--", lw=1); ax2.set_xlabel("Fold enrichment")
    ax2.set_title("B · Top tested Reactome\nNone passed BH FDR < 0.05")

    ax3 = fig.add_subplot(133)
    cells = cell_types.sort_values(["fdr", "p_value"]).head(10).sort_values("effect_size")
    ax3.barh(cells["cell_type"], cells["effect_size"], color="#b279a2")
    ax3.axvline(1, color="#333", ls="--", lw=1); ax3.set_xlabel("Fold enrichment")
    ax3.set_title("C · Top tested brain cell types\nNone passed BH FDR < 0.05")
    fig.suptitle("Parkinson Spatial Vulnerability Biology", fontsize=17, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_09_parkinson_biological_summary.png", dpi=250)
    plt.close(fig)


def web_payload(
    drivers: pd.DataFrame, regional: pd.DataFrame, cell_types: pd.DataFrame,
    selected: pd.DataFrame, expression: pd.DataFrame,
    weighted_genes: pd.DataFrame, regions: pd.DataFrame,
) -> dict[str, object]:
    significant_cells = cell_types.loc[cell_types["fdr"] < FDR_THRESHOLD]
    weights = weighted_genes.drop_duplicates("gene").set_index("gene")["gene_weight"].astype(float)
    represented = weights.index.intersection(expression.columns)
    weights = weights.loc[represented]
    region_records = []
    for region_id, region_meta in regions.iterrows():
        selected_for_interpretation = region_id in selected.index
        if selected_for_interpretation:
            selection = selected.loc[region_id]
            top = drivers.loc[(drivers["region_id"] == region_id) & drivers["top_driver"]].sort_values("regional_rank").head(5)
        else:
            values = expression.loc[region_id, represented].astype(float)
            contribution = values * weights / weights.sum()
            top_genes = contribution.nlargest(5)
            top = pd.DataFrame({
                "gene": top_genes.index,
                "gene_expression": values.loc[top_genes.index].to_numpy(),
                "l2g_score": weights.loc[top_genes.index].to_numpy(),
                "weighted_contribution": top_genes.to_numpy(),
                "regional_rank": range(1, len(top_genes) + 1),
            })
        pathways = regional.loc[
            (regional["region_id"] == region_id)
            & (regional.get("regional_fdr", pd.Series(index=regional.index, dtype=float)).fillna(1) < FDR_THRESHOLD)
        ].sort_values("regional_fdr").head(5)
        driver_set = set(top["gene"])
        cells = significant_cells.loc[
            significant_cells["genes"].fillna("").map(lambda value: bool(driver_set & set(value.split(";"))))
        ].head(3)
        region_records.append({
            "region_id": int(region_id),
            "selected_for_regional_interpretation": selected_for_interpretation,
            "selection_rule": f"top {TOP_REGIONS} by frozen Stage 7 robustness rank" if selected_for_interpretation else "not in pre-specified Stage 9 regional set",
            "robustness_rank": int(selection.robustness_rank) if selected_for_interpretation else None,
            "top_genes": [
                {"gene": row.gene, "expression": row.gene_expression, "l2g_score": row.l2g_score,
                 "weighted_contribution": row.weighted_contribution, "regional_rank": int(row.regional_rank)}
                for row in top.itertuples(index=False)
            ],
            "pathways": [
                {"id": row.pathway_id, "name": row.pathway, "fdr": row.regional_fdr, "genes": str(row.genes).split(";")}
                for row in pathways.itertuples(index=False)
            ],
            "cell_types": [
                {"name": row.cell_type, "fdr": row.fdr, "genes": sorted(driver_set & set(str(row.genes).split(";")))}
                for row in cells.itertuples(index=False)
            ],
        })
    payload = {
        "schema_version": "1.0.0",
        "generated_on": ACCESS_DATE,
        "evidence_note": "Genes and enrichment statistics are evidence; pathway and cell-type labels are biological interpretation and do not establish causality.",
        "regions": region_records,
    }
    serialized = json.dumps(payload, indent=2) + "\n"
    (ROOT / "data" / "web" / "parkinson_biological_interpretation.json").write_text(serialized)
    (WEB_DATA / "parkinson_biological_interpretation.json").write_text(serialized)
    return payload


def update_web_metadata(
    mapping: pd.DataFrame, go: pd.DataFrame, reactome: pd.DataFrame,
    ranked: pd.DataFrame, cells: pd.DataFrame, evidence: pd.DataFrame,
) -> None:
    canonical = ROOT / "data" / "web" / "project_metadata.json"
    metadata = json.loads(canonical.read_text())
    metadata["analysis"]["biological_interpretation"] = {
        "status": "post-hoc interpretation of frozen discovery outputs",
        "primary_gene_set": "Stage 4 L2G-weighted",
        "genes_analyzed": int(len(mapping)),
        "genes_represented_in_ahba": int(mapping["mapping_status"].eq("mapped_and_represented_in_ahba").sum()),
        "background": "15,632 genes retained in the Stage 2 AHBA matrix",
        "go_significant_terms": int((go["fdr"] < FDR_THRESHOLD).sum()),
        "reactome_significant_pathways": int((reactome["fdr"] < FDR_THRESHOLD).sum()),
        "ranked_significant_pathways": int((ranked["fdr"] < FDR_THRESHOLD).sum()) if not ranked.empty else 0,
        "cell_types_significant": int((cells["fdr"] < FDR_THRESHOLD).sum()),
        "supported_region_pathway_rows": int(len(evidence)),
        "region_rule": f"top {TOP_REGIONS} by frozen Stage 7 robustness rank",
        "interpretation_note": "Enrichment and cell-type associations do not establish causality, pathology, pathway activation, or clinical utility",
    }
    additions = [
        {"name": "g:Profiler functional enrichment", "url": "https://biit.cs.ut.ee/gprofiler/"},
        {"name": "Gene Ontology", "url": "https://geneontology.org/"},
        {"name": "Reactome release 97", "url": "https://reactome.org/"},
        {"name": "WikiPathways September 2026", "url": "https://data.wikipathways.org/current/gmt/"},
        {"name": "Human Protein Atlas v25.1 single-nucleus brain", "url": "https://www.proteinatlas.org/humanproteome/single+cell/single+nuclei+brain/data"},
    ]
    existing = {source["name"] for source in metadata["sources"]}
    metadata["sources"].extend(source for source in additions if source["name"] not in existing)
    serialized = json.dumps(metadata, indent=2) + "\n"
    canonical.write_text(serialized)
    (WEB_DATA / "project_metadata.json").write_text(serialized)


def reports(
    mapping: pd.DataFrame, gene_sets: dict[str, pd.DataFrame], go: pd.DataFrame,
    reactome: pd.DataFrame, ranked: pd.DataFrame, cell_types: pd.DataFrame,
    driver_summary: pd.DataFrame, selected: pd.DataFrame, pathways: pd.DataFrame,
    versions: dict[str, object], evidence: pd.DataFrame,
) -> None:
    represented = mapping["mapping_status"].eq("mapped_and_represented_in_ahba").sum()
    go_sig = go.loc[go["fdr"] < FDR_THRESHOLD]
    reac_sig = reactome.loc[reactome["fdr"] < FDR_THRESHOLD]
    cell_sig = cell_types.loc[cell_types["fdr"] < FDR_THRESHOLD]
    version = versions.get("weighted", {}).get("version", "recorded in raw g:Profiler response")
    methods = f"""# Stage 9 methods: biological interpretation

## Analysis boundary

Stage 9 is post-discovery interpretation. Stages 3--8 were frozen and were not
rerun or tuned. The pre-specified primary gene set is the Stage 4 L2G-weighted
set; broad and stringent sets are sensitivity analyses. L2G is genetic
prioritization evidence and is not interpreted as an expression effect or
causal probability.

## Identifiers and background

Original Stage 4 symbols and Ensembl identifiers were preserved. The custom
background is the {len(pd.read_csv(PROCESSED / 'brain_region_gene_expression.csv', nrows=0).columns) - 1:,}
genes retained in the Stage 2 AHBA expression matrix: the genes that could
realistically enter spatial scoring. Unrepresented genes remain in the mapping
audit but are excluded by the custom universe rather than silently discarded.
Stage 4 queries use stable Ensembl identifiers; the symbol-indexed AHBA
background is normalized by g:Profiler. HPA records sharing a display symbol
are collapsed by maximum nCPM within cluster type before symbol-based joining.

## GO and pathway over-representation

g:Profiler `{version}` was accessed programmatically on {ACCESS_DATE}. Sources
were GO Biological Process, Molecular Function, Cellular Component, Reactome,
and WikiPathways. Queries used the custom AHBA background. Raw one-sided
hypergeometric probabilities were recomputed from returned query, overlap,
term and domain counts. Terms required {MIN_TERM_SIZE}--{MAX_TERM_SIZE}
background genes and at least {MIN_OVERLAP} overlapping query gene. Benjamini--Hochberg
FDR was applied separately to GO:BP, GO:MF, GO:CC, Reactome and WikiPathways.
GO terms are overlapping and hierarchical; FDR does not make them independent.
Reactome release { (RAW / 'reactome_version.txt').read_text().strip() } was current at access.

## Ranked pathway analysis

For each eligible Reactome and WikiPathways term, a one-sided Mann--Whitney
rank-sum test asks whether member genes have larger Stage 4 L2G scores than
other prioritized genes. AUC is the rank effect (0.5 under no shift). FDR is
separate by pathway database. This ranks genetic prioritization evidence; it
is not GSEA of differential expression.

## Cell types

Human Protein Atlas v25.1 single-nucleus brain data (34 cluster types across
11 brain regions; Siletti et al. source data) were downloaded on {ACCESS_DATE}.
Within the AHBA background, markers require nCPM >= 1 and log2(expression + 1
over mean-other-types + 1) > 1; the top {HPA_MARKERS_PER_TYPE} per cell type by
specificity form marker sets. Hypergeometric ORA and BH FDR were run as one
34-test family per gene-set definition. Categories are retained under HPA
names; no dopaminergic label was manufactured where the resource does not
provide one.

## Regions and drivers

The top {TOP_REGIONS} parcels by frozen Stage 7 robustness rank were specified
before interpretation. A contribution is `AHBA expression × L2G / sum(L2G)`.
The top {TOP_DRIVERS_PER_REGION} per parcel define regional contributors.
"Spatially contributing prioritized genes" are ranked by the number of these
parcels in which they enter that top-{TOP_DRIVERS_PER_REGION} list, then mean
contribution. They are not called causal drivers.

Reactome ORA was run on each parcel's top contributors. BH correction across
all eligible region--term pairs is the regional family used for support;
within-region FDR is retained descriptively. Fewer than 10 genes would be
reported as insufficient rather than analyzed.

## Network decision

A STRING network was not added. Network centrality would answer a different,
optional question and could invite causal overinterpretation without changing
the requested pathway, cell-type or regional conclusions.

## Multiple testing families

- GO: separate GO:BP, GO:MF and GO:CC families.
- Pathways: separate Reactome and WikiPathways families, per frozen gene set.
- Ranked pathways: separate Reactome and WikiPathways families.
- Cell types: 34 HPA cluster types, separately per frozen gene set.
- Regional pathways: all eligible top-region × Reactome-term tests together.

## Reproducibility and provenance

Raw API requests/responses and official HPA downloads are stored under
`data/biology/raw`. HPA SHA-256: `{sha256(RAW / 'rna_single_nuclei_cluster_type.tsv.zip')}`.
Python {platform.python_version()}, pandas {pd.__version__}, NumPy {np.__version__},
SciPy {scipy.__version__}, requests {requests.__version__}, matplotlib {matplotlib.__version__}.
"""
    (REPORTS / "stage_09_methods.md").write_text(methods)

    def names(frame: pd.DataFrame, column: str, n: int = 8) -> str:
        if frame.empty:
            return "No robust enrichment detected."
        return "; ".join(frame.sort_values(["fdr", "p_value"])[column].head(n))

    interpretation = f"""# Stage 9 interpretation: Parkinson spatial vulnerability biology

## Scope

This post-hoc interpretation used {len(gene_sets['weighted'])} frozen weighted
Stage 4 genes, of which {represented} were represented in AHBA. It did not
alter discovery, null matching, thresholds, spatial robustness or validation.

## Gene Ontology

{names(go_sig.loc[go_sig['ontology'] == 'GO:BP'], 'term_name')}

Molecular function: {names(go_sig.loc[go_sig['ontology'] == 'GO:MF'], 'term_name')}

Cellular component: {names(go_sig.loc[go_sig['ontology'] == 'GO:CC'], 'term_name')}

These results suggest enrichment only where BH FDR is below 0.05. Related GO
terms share genes and ancestry and should be read as programs, not independent
replications.

## Pathways and ranked genetic evidence

Reactome: {names(reac_sig, 'term_name')}

The rank-based analysis asks whether genes assigned to a pathway tend to carry
higher L2G prioritization scores; it does not estimate pathway activation.
{len(ranked.loc[ranked['fdr'] < FDR_THRESHOLD]) if not ranked.empty else 0}
ranked pathway results survived their database-specific BH correction.

## Human-brain cell types

{names(cell_sig, 'cell_type')}

No neuronal, glial or immune category was assumed in advance. HPA cluster-type
labels and data-driven marker specificity determine the ranking.

## Spatially contributing prioritized genes

{'; '.join(driver_summary.head(12)['gene'])}.

These genes repeatedly enter the top-{TOP_DRIVERS_PER_REGION} additive
contributors across frozen high-ranking parcels. This combines healthy-brain
expression with Stage 4 prioritization; it does not make them causal drivers.

## Regional interpretation

{len(evidence)} region--pathway relationships survived the single regional BH
family. The surviving rows all concern FCGR activation in Frontal Sup Medial R,
Hippocampus L, and Amygdala L, and each is supported by only the same two genes,
FYN and FCGR2A. This narrow, repeated two-gene result is exploratory and should
not be generalized into a broad regional mechanism. If no rows had survived,
the correct conclusion would have been “no robust enrichment detected,” not an
expanded gene list. The selected regions retain their original Stage 7 labels;
none is upgraded by Stage 9 biology.

## Sensitivity and uncertainty

Broad and weighted sets have identical membership and therefore identical ORA
but differ conceptually because weighted analysis preserves L2G in spatial and
rank analyses. The stringent set tests dependence on the Stage 4 cutoff.
Only the stringent set yielded an FDR-significant Reactome result (FCGR
activation, FYN/FCGR2A); it was absent in primary broad/weighted ORA, so pathway
robustness across definitions is not supported.
Results remain limited by GWAS-to-gene uncertainty, AHBA donor and sampling
coverage, incomplete pathway/cell annotations, correlated ontologies, and the
negative Stage 8 ENIGMA validation. Biological enrichment does not establish
causality, pathology, clinical utility, or a disease-origin region.
"""
    (REPORTS / "stage_09_interpretation.md").write_text(interpretation)


def write_provenance(versions: dict[str, object]) -> None:
    payload = {
        "stage": 9,
        "access_date": ACCESS_DATE,
        "analysis_boundary": "post-hoc interpretation; Stages 3-8 frozen",
        "sources": {
            "gprofiler": {
                "url": GPROFILER_URL,
                "version": versions.get("weighted", {}).get("version"),
                "cached_responses": sorted(path.name for path in RAW.glob("gprofiler_*.json.gz")),
            },
            "gene_ontology": {"via": "g:Profiler", "url": "https://geneontology.org/"},
            "reactome": {
                "release": (RAW / "reactome_version.txt").read_text().strip(),
                "version_url": "https://reactome.org/ContentService/data/database/version",
            },
            "wikipathways": {
                "release": "2026-09-10",
                "via": "g:Profiler",
                "url": "https://data.wikipathways.org/current/gmt/",
            },
            "human_protein_atlas": {
                "version": "25.1",
                "url": HPA_URL,
                "class_url": HPA_CLASSES_URL,
                "expression_sha256": sha256(RAW / "rna_single_nuclei_cluster_type.tsv.zip"),
                "classes_sha256": sha256(RAW / "rna_single_nuclei_cluster_type_cluster_types.tsv.zip"),
            },
        },
        "frozen_inputs": {
            path.name: sha256(path) for path in [
                GENES / "parkinson_genes_broad.csv",
                GENES / "parkinson_genes_stringent.csv",
                GENES / "parkinson_genes_weighted.csv",
                RESULTS / "parkinson_regional_enrichment.csv",
                RESULTS / "parkinson_spatial_robustness.csv",
                PROCESSED / "brain_region_gene_expression.csv",
            ]
        },
        "packages": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "requests": requests.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    (BIOLOGY / "stage_09_provenance.json").write_text(json.dumps(payload, indent=2) + "\n")


def print_summary(
    mapping: pd.DataFrame, go: pd.DataFrame, reactome: pd.DataFrame,
    cells: pd.DataFrame, drivers: pd.DataFrame, selected: pd.DataFrame,
    sensitivity: pd.DataFrame,
) -> None:
    def top(frame: pd.DataFrame, name: str, n: int = 5) -> str:
        significant = frame.loc[frame["fdr"] < FDR_THRESHOLD]
        return ", ".join(significant.sort_values(["fdr", "p_value"])[name].head(n)) or "No robust enrichment detected."
    sets = {
        key: set(group.loc[group["fdr"] < FDR_THRESHOLD, "term_id"])
        for key, group in sensitivity.loc[sensitivity["ontology"] == "REAC"].groupby("gene_set")
    }
    common = set.intersection(*sets.values()) if sets else set()
    print("\nPARKINSON BIOLOGICAL INTERPRETATION")
    print(f"Genes analyzed: {len(mapping)}")
    print(f"Genes represented in AHBA: {(mapping.mapping_status == 'mapped_and_represented_in_ahba').sum()}")
    print(f"Top GO biological processes: {top(go.loc[go.ontology == 'GO:BP'], 'term_name')}")
    print(f"Top Reactome pathways: {top(reactome, 'term_name')}")
    print(f"Top cell types: {top(cells, 'cell_type')}")
    print(f"Top spatially contributing genes: {', '.join(drivers.head(8).gene)}")
    print(f"Top robust regions: {', '.join(selected.sort_values('robustness_rank').head(5).region_name)}")
    print(f"Pathway robustness across gene-set definitions: {len(common)} Reactome terms significant in all sets")
    print(f"Cell-type robustness: {len(cells.loc[cells.fdr < FDR_THRESHOLD])} weighted-set cell types at FDR < 0.05")
    print("Independent validation status: NOT SUPPORTED (frozen Stage 8 result)")


def main() -> None:
    ensure_directories()
    inputs = read_inputs()
    expression = inputs["expression"]
    regions = inputs["regions"]
    stage6 = inputs["stage6"]
    stage7 = inputs["stage7"]
    validation = inputs["validation"]
    gene_sets = inputs["gene_sets"]
    background = expression.columns.tolist()
    selected = write_region_rule(stage7)
    mapping = gene_id_mapping(gene_sets, set(background))
    go, reactome, sensitivity, versions = functional_enrichment(gene_sets, background)
    ranked = ranked_pathway_analysis(gene_sets, sensitivity)
    marker_sets, _ = hpa_marker_sets(set(background))
    cell_types, cell_sensitivity = cell_type_enrichment(gene_sets, set(background), marker_sets)
    drivers, driver_summary = regional_drivers(expression, gene_sets["weighted"], selected)
    regional = region_pathways(drivers, background, selected)
    evidence = biological_evidence_summary(selected, stage6, regional, cell_types, validation)
    go.to_csv(WEB_DATA / "stage_09_go_enrichment.csv", index=False)
    reactome.to_csv(WEB_DATA / "stage_09_reactome_enrichment.csv", index=False)
    cell_types.to_csv(WEB_DATA / "stage_09_cell_type_enrichment.csv", index=False)
    evidence.to_csv(WEB_DATA / "stage_09_biological_evidence_summary.csv", index=False)
    plot_pathways(reactome)
    plot_cell_types(cell_types)
    plot_drivers(drivers)
    plot_sensitivity(sensitivity)
    plot_genetic_region_biology(drivers, driver_summary, regional, gene_sets["weighted"])
    plot_summary(regions, stage6, reactome, cell_types)
    web_payload(drivers, regional, cell_types, selected, expression, gene_sets["weighted"], regions)
    update_web_metadata(mapping, go, reactome, ranked, cell_types, evidence)
    reports(mapping, gene_sets, go, reactome, ranked, cell_types, driver_summary,
            selected, regional, versions, evidence)
    write_provenance(versions)
    print_summary(mapping, go, reactome, cell_types, driver_summary, selected, sensitivity)


if __name__ == "__main__":
    main()
