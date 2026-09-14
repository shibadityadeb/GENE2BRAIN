"""Finalize the Stage 10 quantitative summary and result manifest."""

from __future__ import annotations

from datetime import date
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results"
REPORTS = ROOT / "reports"
TABLES = ROOT / "results" / "tables"
CONFIG = ROOT / "config"

SCIENTIFIC_PACKAGES = (
    "abagen", "pandas", "numpy", "scipy", "nibabel", "nilearn",
    "matplotlib", "seaborn", "statsmodels", "scikit-learn", "requests", "jupyter",
)


def software_versions() -> dict[str, str]:
    versions = {"Python": platform.python_version()}
    for package in SCIENTIFIC_PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not installed in this execution environment"
    return versions


def write_environment_report(versions: dict[str, str]) -> None:
    rows = "\n".join(f"| {name} | {value} |" for name, value in versions.items())
    (REPORTS / "stage_10_software_environment.md").write_text(
        "# Stage 10 software environment\n\n"
        "Exact versions captured from the Python environment that finalized the "
        "committed Stage 10 artifacts. The declared clean-environment dependencies "
        "remain in `requirements.txt`.\n\n"
        "| Component | Version |\n| --- | --- |\n"
        f"{rows}\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def statistical_summary() -> dict[str, object]:
    similarity = pd.read_csv(RESULTS / "disease_spatial_similarity.csv")
    pearson_pairs = similarity.loc[(similarity.metric == "pearson") & (similarity.disease_1 != similarity.disease_2)].copy()
    pearson_pairs["pair"] = pearson_pairs.apply(lambda row: "|".join(sorted((row.disease_1, row.disease_2))), axis=1)
    pearson_pairs = pearson_pairs.drop_duplicates("pair")
    similar = pearson_pairs.nlargest(1, "correlation").iloc[0]
    distinct = pearson_pairs.nsmallest(1, "correlation").iloc[0]
    shared = pd.read_csv(RESULTS / "shared_brain_region_enrichment.csv")
    specificity = pd.read_csv(RESULTS / "disease_specific_regional_signatures.csv")
    power = pd.read_csv(RESULTS / "disease_power_relationships.csv")
    robustness = pd.read_csv(RESULTS / "multidisease_gene_set_robustness.csv")
    biology = pd.read_csv(RESULTS / "disease_spatial_vs_biological_similarity.csv")
    valid_biology = biology.dropna(subset=["significant_pathway_jaccard"])
    spatial_biology = (
        spearmanr(valid_biology.spatial_pearson_r, valid_biology.significant_pathway_jaccard)
        if len(valid_biology) >= 3 else None
    )
    power_tests = {}
    for predictor in ("gwas_sample_size", "gene_set_size"):
        for outcome in ("number_of_enriched_regions", "maximum_z_score", "number_of_spatially_robust_regions"):
            if power[outcome].nunique() > 1:
                result = spearmanr(power[predictor], power[outcome])
                power_tests[f"{predictor}_vs_{outcome}"] = {"rho": float(result.statistic), "p": float(result.pvalue)}
            else:
                power_tests[f"{predictor}_vs_{outcome}"] = {"rho": None, "p": None, "note": "outcome constant"}
    return {
        "most_similar": (similar.disease_1, similar.disease_2, float(similar.correlation), float(similar.p_value)),
        "most_distinct": (distinct.disease_1, distinct.disease_2, float(distinct.correlation), float(distinct.p_value)),
        "correlation_range": (float(pearson_pairs.correlation.min()), float(pearson_pairs.correlation.max())),
        "median_abs_correlation": float(pearson_pairs.correlation.abs().median()),
        "most_recurrent": shared.iloc[0].to_dict(),
        "specific_count": int((specificity.fdr_p < 0.05).sum()),
        "top_specific": specificity.iloc[0].to_dict(),
        "power_tests": power_tests,
        "gene_set_robustness_range": (float(robustness.spearman_rho.min()), float(robustness.spearman_rho.max())),
        "spatial_biology": None if spatial_biology is None else {
            "rho": float(spatial_biology.statistic), "p": float(spatial_biology.pvalue), "n_pairs": len(valid_biology)
        },
    }


def write_final_report(summary: dict[str, object]) -> None:
    master = pd.read_csv(TABLES / "gene2brain_master_disease_results.csv")
    analyzed = master.loc[master.top_region.notna()]
    not_analyzed = master.loc[master.top_region.isna()]
    weighted = pd.read_csv(RESULTS / "multidisease_pathway_enrichment.csv")
    cells = pd.read_csv(RESULTS / "multidisease_cell_type_enrichment.csv")
    p1, p2, pr, pp = summary["most_similar"]
    d1, d2, dr, dp = summary["most_distinct"]
    recurrent = summary["most_recurrent"]
    top_specific = summary["top_specific"]
    power_lines = []
    for name, result in summary["power_tests"].items():
        if result.get("rho") is None:
            power_lines.append(f"- {name}: not estimable ({result['note']})")
        else:
            power_lines.append(f"- {name}: Spearman rho={result['rho']:.3f}, p={result['p']:.4g}")
    biology_result = summary["spatial_biology"]
    biology_line = (
        f"Across {biology_result['n_pairs']} disease pairs with at least one significant pathway set, spatial correlation versus pathway Jaccard similarity had Spearman rho={biology_result['rho']:.3f}, p={biology_result['p']:.4g}."
        if biology_result else
        "Too few disease pairs had non-empty FDR-significant pathway sets for a correlation test."
    )
    disease_rows = "\n".join(
        f"| {row.disease} | {row.gwas_sample_size if pd.notna(row.gwas_sample_size) else 'NA'} | {row.n_gwas_loci} | {row.n_prioritized_genes} | {row.ahba_gene_coverage if pd.notna(row.ahba_gene_coverage) else 'NA'} | {row.n_significant_regions} | {row.n_spatially_robust_regions} | {row.validation_status} | {row.top_region if pd.notna(row.top_region) else 'not analyzed'} |"
        for row in master.itertuples()
    )
    text = f"""# GENE2BRAIN Stage 10 final multi-disease summary

## Scope

- Panel screened: {len(master)} diseases
- Diseases analyzed: {len(analyzed)}
- Excluded or held at needs-review: {len(not_analyzed)}
- Brain regions: 138 AAL3 parcels
- AHBA genes: 15,632
- Primary comparison statistic: regional matched-gene-set permutation Z score
- Weighted-set regional BH-FDR discoveries: {int(analyzed.n_significant_regions.sum())}
- Jointly spatially robust discoveries: {int(analyzed.n_spatially_robust_regions.sum())}

## Disease results

| Disease | GWAS N | Loci | Prioritized genes | AHBA coverage % | Significant regions | Spatially robust | Validation | Top Z region |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
{disease_rows}

## Quantitative cross-disease questions

1. **Distinct regional signatures.** Pairwise Pearson correlations span {summary['correlation_range'][0]:.3f} to {summary['correlation_range'][1]:.3f}; median absolute correlation is {summary['median_abs_correlation']:.3f}. This quantifies heterogeneity but is not a disease-classification test.
2. **Most similar pair.** {p1} and {p2}: Pearson r={pr:.3f}, nominal p={pp:.4g} across 138 regions.
3. **Most distinct pair.** {d1} and {d2}: Pearson r={dr:.3f}, nominal p={dp:.4g} across 138 regions.
4. **Repeated region.** {recurrent['region']} has {int(recurrent['number_of_diseases_with_significant_enrichment'])} diseases with gene-set FDR significance and {int(recurrent['number_of_diseases_with_spatially_robust_enrichment'])} with joint spatial robustness. Recurrence is not universal vulnerability.
5. **Disease-specific regions.** {summary['specific_count']} disease-region cells pass the exploratory global BH-FDR test based on leave-one-disease-out specificity. The strongest is {top_specific['disease']} in {top_specific['region']} (specificity Z={top_specific['cross_disease_specificity_z']:.3f}, FDR={top_specific['fdr_p']:.4g}).
6. **Spatial versus biological similarity.** {biology_line} This is exploratory and does not establish shared causal biology.
7. **Gene-set robustness.** Broad/stringent/weighted regional Spearman correlations range from {summary['gene_set_robustness_range'][0]:.3f} to {summary['gene_set_robustness_range'][1]:.3f} across diseases. Epilepsy's two-gene stringent set is explicitly unstable.

## GWAS power and gene-set size

{chr(10).join(power_lines)}

These ten-disease tests have limited power and correlated regional observations.
They diagnose possible confounding; they do not prove that GWAS power causes the
spatial results.

## Biology

The standardized tables contain {len(weighted):,} GO/pathway rows and {len(cells):,}
HPA brain cell-type rows. Weighted-set BH-FDR discoveries total
{int(((weighted.gene_set == 'weighted') & (weighted.fdr < 0.05)).sum())} ontology/pathway
rows and {int(((cells.gene_set == 'weighted') & (cells.fdr < 0.05)).sum())} cell-type
rows. Correlated and nested ontology terms are not counted as independent
mechanisms.

## Negative-results statement

No disease produced a region passing both gene-set FDR and spatial-null FDR.
This null spatial-robustness result is retained. No threshold, GWAS, gene-set
size, matching variable, null model, or region list was changed to improve it.
Parkinson independent validation remains NOT SUPPORTED; other analyzed diseases
remain pending rather than receiving artificial validation data.

## Final visual package

1. `stage_10_gene2brain_pipeline.png`
2. existing Parkinson genetic enrichment figure from Stage 6
3. existing Parkinson independent-validation figure from Stage 8
4. `stage_10_multidisease_signature_matrix.png`
5. `stage_10_disease_spatial_similarity.png`
6. `stage_10_shared_vs_specific.png`
7. `stage_10_disease_pca.png`
8. `stage_10_biological_programs.png`

The atlas is a research visualization of healthy-brain expression enrichment.
It is not a clinical classifier, diagnostic tool, patient-level prediction, or
causal map.
"""
    (REPORTS / "stage_10_final_summary.md").write_text(text, encoding="utf-8")


def manifest(versions: dict[str, str]) -> pd.DataFrame:
    patterns = [
        "config/disease_panel.yaml", "config/statistical_parameters.yaml",
        "data/gwas/multidisease_*.csv", "data/genes/multidisease/*.csv",
        "data/genes/multidisease_gene_locus_evidence.csv",
        "data/results/multidisease_*.csv", "data/results/disease_*.csv",
        "data/results/shared_brain_region_enrichment.csv",
        "data/validation/multidisease_validation_metadata.csv",
        "results/tables/gene2brain_master_disease_results.csv",
        "results/figures/stage_10_*.png", "reports/stage_10_*.md",
        "reports/diseases/*.md", "data/web/multidisease_atlas.json",
    ]
    paths: set[Path] = set()
    for pattern in patterns:
        paths.update(ROOT.glob(pattern))
    rows = []
    version_text = "; ".join(f"{name} {value}" for name, value in versions.items())
    for path in sorted(paths):
        relative = path.relative_to(ROOT)
        disease = "all"
        if relative.parts[:3] == ("data", "genes", "multidisease"):
            disease = path.stem.rsplit("_", 1)[0]
        elif relative.parts[:2] == ("reports", "diseases"):
            disease = path.stem
        rows.append({
            "disease": disease, "file": str(relative), "analysis_stage": 10,
            "source": "generated by GENE2BRAIN Stage 10 from frozen GWAS/Open Targets/AHBA inputs",
            "generated_date": date.today().isoformat(),
            "software_version": f"GENE2BRAIN stage10; {version_text}",
            "status": "complete", "bytes": path.stat().st_size, "sha256": sha256(path),
        })
    output = pd.DataFrame(rows)
    output.to_csv(RESULTS / "gene2brain_data_manifest.csv", index=False)
    return output


def main() -> None:
    versions = software_versions()
    write_environment_report(versions)
    summary = statistical_summary()
    write_final_report(summary)
    output = manifest(versions)
    print(f"Finalized Stage 10: {len(output)} manifest entries")


if __name__ == "__main__":
    main()
