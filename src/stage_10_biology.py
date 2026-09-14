"""Stage 10 standardized functional, cell-type, and regional-driver analysis."""

from __future__ import annotations

from datetime import date
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

import stage_09_biological_interpretation as reference_stage9


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
GENES = ROOT / "data" / "genes" / "multidisease"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "results"
BIOLOGY = ROOT / "data" / "biology"
RAW = BIOLOGY / "raw"
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
DISEASE_REPORTS = REPORTS / "diseases"

METHODS = ("weighted", "broad", "stringent")
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "WP"]
TOP_REGIONS = 10
TOP_DRIVERS = 20
FDR = 0.05


def load_json_yaml(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def analyzed_diseases() -> list[dict]:
    panel = load_json_yaml(CONFIG / "disease_panel.yaml")["diseases"]
    qc = pd.read_csv(RESULTS / "multidisease_qc_summary.csv").set_index("disease")
    return [item for item in panel if item["disease_name"] in qc.index and qc.loc[item["disease_name"], "analysis_status"] == "analyzable"]


def functional_enrichment(diseases: list[dict], background: list[str]) -> pd.DataFrame:
    frames = []
    for disease in diseases:
        disease_id = disease["disease_id"]
        for method in ("weighted", "stringent"):
            genes = pd.read_csv(GENES / f"{disease_id}_{method}.csv")
            genes = genes.loc[genes["gene"].isin(background)].drop_duplicates("gene")
            query = genes["ensembl_id"].astype(str).str.split(".").str[0].tolist()
            labels = dict(zip(query, genes["gene"]))
            wrapped = reference_stage9.request_gprofiler(
                f"stage10_{disease_id}_{method}_ora.json.gz",
                query,
                background,
                SOURCES,
            )
            parsed = reference_stage9.parse_ora(wrapped, query, method, labels)
            if parsed.empty:
                continue
            parsed.insert(0, "disease", disease["disease_name"])
            parsed.insert(1, "disease_id", disease_id)
            frames.append(parsed)
        weighted = next((frame for frame in reversed(frames) if not frame.empty and frame.iloc[0]["disease_id"] == disease_id and frame.iloc[0]["gene_set"] == "weighted"), None)
        if weighted is not None:
            broad = weighted.copy()
            broad["gene_set"] = "broad"
            frames.append(broad)
    columns = [
        "disease", "disease_id", "gene_set", "ontology", "term_id", "term_name",
        "overlap", "gene_count", "background_count", "effect_size", "p_value", "fdr", "genes",
        "effective_domain_size", "service_adjusted_p",
    ]
    output = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=columns)
    output = output[[column for column in columns if column in output.columns]]
    output.to_csv(RESULTS / "multidisease_pathway_enrichment.csv", index=False)
    return output


def cell_type_enrichment(diseases: list[dict], background: set[str]) -> pd.DataFrame:
    markers = pd.read_csv(BIOLOGY / "hpa_brain_cell_type_markers.csv")
    marker_sets = {name: set(group["gene"]) & background for name, group in markers.groupby("cell_type")}
    marker_classes = markers.drop_duplicates("cell_type").set_index("cell_type")["cell_type_class"]
    rows = []
    for disease in diseases:
        for method in METHODS:
            genes = pd.read_csv(GENES / f"{disease['disease_id']}_{method}.csv")
            query = set(genes["gene"]) & background
            for cell_type, marker_set in marker_sets.items():
                overlap = sorted(query & marker_set)
                p_value = float(hypergeom.sf(len(overlap) - 1, len(background), len(marker_set), len(query)))
                effect = (len(overlap) / len(query)) / (len(marker_set) / len(background))
                rows.append({
                    "disease": disease["disease_name"], "disease_id": disease["disease_id"],
                    "gene_set": method, "cell_type": cell_type,
                    "cell_type_class": marker_classes.get(cell_type, ""),
                    "marker_source": "Human Protein Atlas v25.1 single-nucleus brain; top 200 specificity-ranked AHBA genes per cluster type",
                    "overlap": len(overlap), "gene_count": len(query),
                    "background_count": len(marker_set), "effect_size": effect,
                    "p_value": p_value, "genes": ";".join(overlap),
                })
    output = pd.DataFrame(rows)
    output["fdr"] = np.nan
    for _, index in output.groupby(["disease", "gene_set"]).groups.items():
        output.loc[index, "fdr"] = multipletests(output.loc[index, "p_value"], method="fdr_bh")[1]
    output = output.sort_values(["disease", "gene_set", "fdr", "p_value"])
    output.to_csv(RESULTS / "multidisease_cell_type_enrichment.csv", index=False)
    return output


def driver_genes(diseases: list[dict], expression: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    spatial = pd.read_csv(RESULTS / "multidisease_spatial_robustness.csv")
    rows = []
    for disease in diseases:
        selected = spatial.loc[
            (spatial["disease"] == disease["disease_name"]) & (spatial["gene_set"] == "weighted")
        ].nsmallest(TOP_REGIONS, "robustness_rank")
        weighted = pd.read_csv(GENES / f"{disease['disease_id']}_weighted.csv").drop_duplicates("gene").set_index("gene")
        genes = weighted.index.intersection(expression.columns)
        weights = weighted.loc[genes, "gene_weight"].astype(float)
        for region in selected.itertuples(index=False):
            values = expression.loc[int(region.region_id), genes].astype(float)
            contributions = values * weights / weights.sum()
            ranks = contributions.rank(method="min", ascending=False).astype(int)
            for gene in genes:
                rows.append({
                    "disease": disease["disease_name"], "disease_id": disease["disease_id"],
                    "region": region.region, "region_id": region.region_id,
                    "region_robustness_rank": region.robustness_rank,
                    "region_spatially_robust": region.spatially_robust,
                    "gene": gene, "ensembl_id": weighted.loc[gene, "ensembl_id"],
                    "expression": values[gene], "l2g_score": weights[gene],
                    "weighted_contribution": contributions[gene], "regional_rank": ranks[gene],
                    "top_driver": bool(ranks[gene] <= TOP_DRIVERS),
                })
    regional = pd.DataFrame(rows).sort_values(["disease", "region_robustness_rank", "regional_rank"])
    top = regional.loc[regional["top_driver"]]
    summary = top.groupby(["disease", "disease_id", "gene", "ensembl_id"], as_index=False).agg(
        l2g_score=("l2g_score", "first"),
        number_of_top_regions=("region", "nunique"),
        mean_expression=("expression", "mean"),
        mean_weighted_contribution=("weighted_contribution", "mean"),
        maximum_weighted_contribution=("weighted_contribution", "max"),
        regions=("region", lambda value: ";".join(sorted(set(value)))),
    ).sort_values(["disease", "number_of_top_regions", "mean_weighted_contribution"], ascending=[True, False, False])
    summary.to_csv(RESULTS / "multidisease_driver_genes.csv", index=False)
    regional.to_csv(RESULTS / "multidisease_regional_gene_drivers.csv", index=False)
    return regional, summary


def biological_similarity(pathways: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    spatial = pd.read_csv(RESULTS / "disease_spatial_similarity.csv")
    spatial = spatial.loc[(spatial["metric"] == "pearson") & (spatial["disease_1"] != spatial["disease_2"])].copy()
    spatial["pair"] = spatial.apply(lambda row: "|".join(sorted((row.disease_1, row.disease_2))), axis=1)
    spatial = spatial.drop_duplicates("pair")
    weighted_pathways = pathways.loc[pathways["gene_set"] == "weighted"]
    weighted_cells = cells.loc[cells["gene_set"] == "weighted"]
    rows = []
    for item in spatial.itertuples(index=False):
        path_a = set(weighted_pathways.loc[(weighted_pathways.disease == item.disease_1) & (weighted_pathways.fdr < FDR), "term_id"])
        path_b = set(weighted_pathways.loc[(weighted_pathways.disease == item.disease_2) & (weighted_pathways.fdr < FDR), "term_id"])
        cell_a = set(weighted_cells.loc[(weighted_cells.disease == item.disease_1) & (weighted_cells.fdr < FDR), "cell_type"])
        cell_b = set(weighted_cells.loc[(weighted_cells.disease == item.disease_2) & (weighted_cells.fdr < FDR), "cell_type"])
        path_union = path_a | path_b
        cell_union = cell_a | cell_b
        rows.append({
            "disease_1": item.disease_1, "disease_2": item.disease_2,
            "spatial_pearson_r": item.correlation,
            "significant_pathway_jaccard": len(path_a & path_b) / len(path_union) if path_union else np.nan,
            "significant_cell_type_jaccard": len(cell_a & cell_b) / len(cell_union) if cell_union else np.nan,
            "n_pathways_disease_1": len(path_a), "n_pathways_disease_2": len(path_b),
            "n_cell_types_disease_1": len(cell_a), "n_cell_types_disease_2": len(cell_b),
            "interpretation": "pairwise descriptive comparison; NA when both FDR-significant sets are empty; no inferential biological-similarity claim",
        })
    output = pd.DataFrame(rows)
    output.to_csv(RESULTS / "disease_spatial_vs_biological_similarity.csv", index=False)
    return output


def biology_figure(pathways: pd.DataFrame) -> None:
    primary = pathways.loc[(pathways["gene_set"] == "weighted") & (pathways["ontology"].isin(["GO:BP", "REAC", "WP"]))].copy()
    if primary.empty:
        return
    primary["strength"] = -np.log10(primary["fdr"].clip(lower=1e-300))
    top_terms = primary.groupby(["term_id", "term_name"], as_index=False)["strength"].max().nlargest(24, "strength")
    chosen = set(top_terms["term_id"])
    matrix = primary.loc[primary["term_id"].isin(chosen)].pivot_table(index="disease", columns="term_name", values="strength", aggfunc="max", fill_value=0)
    matrix = matrix.reindex(columns=top_terms.sort_values("strength", ascending=False)["term_name"].drop_duplicates())
    fig, ax = plt.subplots(figsize=(18, 8))
    sns.heatmap(matrix, cmap="mako", ax=ax, cbar_kws={"label": "−log10(BH FDR); 0 = term not returned"})
    ax.set(title="Biological Programs Across Disease Spatial Signatures", xlabel="GO biological process / pathway", ylabel="Disease")
    ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_10_biological_programs.png", dpi=220)
    plt.close(fig)


def update_master_and_reports(pathways: pd.DataFrame, cells: pd.DataFrame, drivers: pd.DataFrame) -> None:
    master_path = TABLES / "gene2brain_master_disease_results.csv"
    master = pd.read_csv(master_path)
    for index, row in master.iterrows():
        disease = row["disease"]
        if pd.isna(row["top_region"]):
            continue
        go = pathways.loc[(pathways.disease == disease) & (pathways.gene_set == "weighted") & (pathways.ontology == "GO:BP") & (pathways.fdr < FDR)]
        cell = cells.loc[(cells.disease == disease) & (cells.gene_set == "weighted") & (cells.fdr < FDR)]
        master.loc[index, "top_biological_process"] = go.sort_values(["fdr", "p_value"]).iloc[0].term_name if not go.empty else "No BH-FDR-significant GO:BP term"
        master.loc[index, "top_cell_type"] = cell.sort_values(["fdr", "p_value"]).iloc[0].cell_type if not cell.empty else "No BH-FDR-significant HPA brain cell type"
        report_path = DISEASE_REPORTS / f"{next(item['disease_id'] for item in analyzed_diseases() if item['disease_name'] == disease)}.md"
        text = report_path.read_text(encoding="utf-8")
        biology = (
            f"Weighted-set GO biological processes at BH-FDR < 0.05: {len(go)}. "
            f"Weighted-set HPA brain cell types at BH-FDR < 0.05: {len(cell)}. "
            f"Top repeated spatial driver genes: {', '.join(drivers.loc[drivers.disease == disease].head(5).gene)}."
        )
        text = text.replace(
            "Stage 10 standardized pathway, cell-type, and driver-gene tables are populated\nonly for diseases that pass the spatial-analysis gate. Empty/significance-null\nresults are preserved.",
            biology,
        )
        report_path.write_text(text, encoding="utf-8")
    master.to_csv(master_path, index=False)


def provenance_report(pathways: pd.DataFrame, cells: pd.DataFrame, drivers: pd.DataFrame, similarity: pd.DataFrame) -> None:
    significant_pathways = int(((pathways.gene_set == "weighted") & (pathways.fdr < FDR)).sum())
    significant_cells = int(((cells.gene_set == "weighted") & (cells.fdr < FDR)).sum())
    report = f"""# Stage 10 biological interpretation

- Access date: {date.today().isoformat()}
- Diseases analyzed: {len(analyzed_diseases())}
- Pathway/ontology rows: {len(pathways):,}
- Weighted-set ontology/pathway rows at BH-FDR < 0.05: {significant_pathways}
- Cell-type rows: {len(cells):,}
- Weighted-set HPA cell-type rows at BH-FDR < 0.05: {significant_cells}
- Driver-gene summary rows: {len(drivers):,}
- Disease-pair spatial/biological comparisons: {len(similarity)}

Functional enrichment uses g:Profiler with GO biological process, molecular
function, cellular component, Reactome, and WikiPathways sources and the frozen
15,632-gene AHBA background. BH correction is performed separately by ontology
inside each disease and gene-set request, matching Stage 9. Cell-type enrichment
uses the frozen HPA v25.1 human-brain single-nucleus marker table and a
hypergeometric test with BH correction within disease and gene-set version.

Regional drivers are additive `AHBA expression × L2G / sum(L2G)` contributions
in the ten regions selected by the already-computed spatial robustness rank.
Because no disease has a jointly spatially robust parcel, these driver results
are explicitly ranked exploratory interpretation, not regional discovery.

The spatial-versus-biological table is descriptive. It does not test or claim
that spatially similar diseases share causal biology, particularly when both
diseases have no FDR-significant biological terms.
"""
    (REPORTS / "stage_10_biological_interpretation.md").write_text(report, encoding="utf-8")


def main() -> None:
    for directory in (RESULTS, RAW, FIGURES, REPORTS, DISEASE_REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    diseases = analyzed_diseases()
    expression = pd.read_csv(PROCESSED / "brain_region_gene_expression.csv", index_col="region_id")
    background = expression.columns.tolist()
    print(f"Running standardized biology for {len(diseases)} diseases", flush=True)
    pathways = functional_enrichment(diseases, background)
    cells = cell_type_enrichment(diseases, set(background))
    _, drivers = driver_genes(diseases, expression)
    similarity = biological_similarity(pathways, cells)
    biology_figure(pathways)
    update_master_and_reports(pathways, cells, drivers)
    provenance_report(pathways, cells, drivers, similarity)
    print(f"Pathway rows: {len(pathways):,}; weighted FDR-significant: {((pathways.gene_set == 'weighted') & (pathways.fdr < FDR)).sum()}")
    print(f"Cell-type rows: {len(cells):,}; weighted FDR-significant: {((cells.gene_set == 'weighted') & (cells.fdr < FDR)).sum()}")
    print(f"Driver-gene rows: {len(drivers):,}")


if __name__ == "__main__":
    main()
