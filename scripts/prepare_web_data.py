#!/usr/bin/env python3
"""Prepare read-only web artifacts from validated GENE2BRAIN outputs.

This script never recomputes scientific statistics. It serializes Stage 6 values
and extracts rectilinear parcel surfaces from the exact AAL3v1 NIfTI labels used
in Stage 2. Adjacent coplanar voxel faces are greedily merged to reduce transfer
size without changing parcel boundaries or region identifiers.
"""

from __future__ import annotations

import csv
import json
import math
import re
import shutil
from collections.abc import Iterable
from datetime import date
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results" / "parkinson_regional_enrichment.csv"
REGIONS = ROOT / "data" / "processed" / "region_metadata.csv"
PARAMETERS = ROOT / "data" / "results" / "stage_06_permutation_parameters.json"
SPATIAL_RESULTS = ROOT / "data" / "results" / "parkinson_spatial_robustness.csv"
SPATIAL_PARAMETERS = ROOT / "data" / "results" / "stage_07_spatial_null_parameters.json"
GLOBAL_SPATIAL = ROOT / "data" / "results" / "parkinson_global_spatial_test.csv"
VALIDATION_AAL3 = ROOT / "data" / "validation" / "processed" / "parkinson_validation_aal3_scores.csv"
VALIDATION_RESULTS = ROOT / "data" / "results" / "parkinson_independent_validation.csv"
VALIDATION_REGIONAL = ROOT / "data" / "validation" / "processed" / "parkinson_validation_regional_scores.csv"
BIOLOGY_JSON = ROOT / "data" / "web" / "parkinson_biological_interpretation.json"
STAGE9_GO = ROOT / "data" / "results" / "stage_09_go_enrichment.csv"
STAGE9_REACTOME = ROOT / "data" / "results" / "stage_09_reactome_enrichment.csv"
STAGE9_RANKED = ROOT / "data" / "results" / "stage_09_ranked_pathway_analysis.csv"
STAGE9_CELLS = ROOT / "data" / "results" / "stage_09_cell_type_enrichment.csv"
STAGE9_EVIDENCE = ROOT / "data" / "results" / "stage_09_biological_evidence_summary.csv"
ATLAS = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.nii.gz"
ATLAS_XML = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.xml"
AHBA_SUMMARY = ROOT / "data" / "processed" / "ahba_metadata_summary.csv"
EXPRESSION = ROOT / "data" / "processed" / "brain_region_gene_expression.csv"
STUDY_SELECTION = ROOT / "data" / "gwas" / "parkinson_study_selection.csv"
CREDIBLE_SETS = ROOT / "data" / "gwas" / "parkinson_credible_sets.csv"
WEB_DATA = ROOT / "data" / "web"
PUBLIC_DATA = ROOT / "web" / "public" / "data"

NUMERIC_FIELDS = (
    "observed_score", "random_mean", "random_std", "z_score", "empirical_p",
    "fdr_p", "effect_size",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def json_number(value: str) -> int | float:
    """Parse a CSV token without rounding beyond standard JSON number precision."""
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite numeric token: {value!r}")
    return parsed


def greedy_rectangles(mask: np.ndarray) -> Iterable[tuple[int, int, int, int]]:
    """Yield maximal non-overlapping [row, row_end) x [col, col_end) rectangles."""
    remaining = mask.copy()
    rows, cols = remaining.shape
    for row in range(rows):
        col = 0
        while col < cols:
            if not remaining[row, col]:
                col += 1
                continue
            col_end = col + 1
            while col_end < cols and remaining[row, col_end]:
                col_end += 1
            row_end = row + 1
            while row_end < rows and remaining[row_end, col:col_end].all():
                row_end += 1
            remaining[row:row_end, col:col_end] = False
            yield row, row_end, col, col_end
            col = col_end


def exposed_slice(mask: np.ndarray, axis: int, index: int, direction: int) -> np.ndarray:
    current = np.take(mask, index, axis=axis)
    neighbor_index = index + direction
    if neighbor_index < 0 or neighbor_index >= mask.shape[axis]:
        return current
    return current & ~np.take(mask, neighbor_index, axis=axis)


def make_quad(
    axis: int,
    fixed: float,
    row_axis: int,
    col_axis: int,
    row_start: float,
    row_end: float,
    col_start: float,
    col_end: float,
) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    for row, col in (
        (row_start, col_start), (row_end, col_start),
        (row_end, col_end), (row_start, col_end),
    ):
        voxel = [0.0, 0.0, 0.0]
        voxel[axis] = fixed
        voxel[row_axis] = row
        voxel[col_axis] = col
        points.append(tuple(voxel))
    return points


def parcel_mesh(mask: np.ndarray, affine: np.ndarray) -> tuple[list[float], list[int]]:
    vertices: list[float] = []
    indices: list[int] = []
    vertex_lookup: dict[tuple[float, float, float], int] = {}

    def add_vertex(voxel: tuple[float, float, float]) -> int:
        mni = nib.affines.apply_affine(affine, voxel)
        # Web scene convention: MNI left/right -> X, superior -> Y, anterior -> Z.
        web = (float(mni[0]), float(mni[2]), float(mni[1]))
        key = tuple(round(value, 5) for value in web)
        if key not in vertex_lookup:
            vertex_lookup[key] = len(vertices) // 3
            vertices.extend(web)
        return vertex_lookup[key]

    for axis in range(3):
        other = [candidate for candidate in range(3) if candidate != axis]
        row_axis, col_axis = other
        occupied = np.where(mask.any(axis=tuple(other)))[0]
        for index in occupied:
            for direction in (-1, 1):
                face_mask = exposed_slice(mask, axis, int(index), direction)
                if not face_mask.any():
                    continue
                for row0, row1, col0, col1 in greedy_rectangles(face_mask):
                    quad = make_quad(
                        axis, index + direction * 0.5, row_axis, col_axis,
                        row0 - 0.5, row1 - 0.5, col0 - 0.5, col1 - 0.5,
                    )
                    face = [add_vertex(point) for point in quad]
                    if direction > 0:
                        indices.extend((face[0], face[1], face[2], face[0], face[2], face[3]))
                    else:
                        indices.extend((face[0], face[2], face[1], face[0], face[3], face[2]))
    return vertices, indices


def expression_dimensions() -> tuple[int, int]:
    with EXPRESSION.open(encoding="utf-8") as handle:
        header = next(csv.reader(handle))
        row_count = sum(1 for _ in handle)
    return row_count, len(header) - 1


def ahba_counts() -> tuple[int, int]:
    rows = read_csv(AHBA_SUMMARY)
    annotations = [row for row in rows if row["file_type"] == "annotation"]
    return len({row["donor"] for row in annotations}), sum(int(row["n_rows"]) for row in annotations)


def write_json(path: Path, payload: Any, *, compact: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=None if compact else 2, separators=(",", ":") if compact else None)
        handle.write("\n")


def build() -> None:
    results = read_csv(RESULTS)
    metadata = {int(row["region_id"]): row for row in read_csv(REGIONS)}
    parameters = json.loads(PARAMETERS.read_text(encoding="utf-8"))
    spatial_parameters = json.loads(SPATIAL_PARAMETERS.read_text(encoding="utf-8"))
    spatial_results = {int(row["region_id"]): row for row in read_csv(SPATIAL_RESULTS)}
    validation_results = {int(row["region_id"]): row for row in read_csv(VALIDATION_AAL3)}
    atlas_image = nib.load(ATLAS)
    atlas = np.asarray(atlas_image.dataobj, dtype=np.int16)

    records: list[dict[str, Any]] = []
    meshes: list[dict[str, Any]] = []
    for result in results:
        region_id = int(result["region_id"])
        region = metadata[region_id]
        record: dict[str, Any] = {
            "region_id": region_id,
            "region_name": result["region_name"],
            "atlas_id": f"AAL3v1:{region_id}",
            "atlas_label": region["atlas_label"],
            "hemisphere": region["hemisphere"],
            "broad_system": region["broad_system"],
            "centroid_mni": [
                json_number(region["centroid_mni_x"]),
                json_number(region["centroid_mni_y"]),
                json_number(region["centroid_mni_z"]),
            ],
            **{field: json_number(result[field]) for field in NUMERIC_FIELDS},
            "number_of_genes": int(parameters["genes_per_set"][result["gene_set_version"]]),
            "spatial_null_p": json_number(spatial_results[region_id]["spatial_null_p"]),
            "spatial_null_fdr": json_number(spatial_results[region_id]["spatial_null_fdr"]),
            "spatial_robustness": json_number(spatial_results[region_id]["spatial_robustness"]),
            "spatial_robustness_label": spatial_results[region_id]["spatial_robustness_label"],
            "robustness_rank": int(spatial_results[region_id]["robustness_rank"]),
            "spatial_isolate": spatial_results[region_id]["spatial_isolate"] == "True",
            "validation_score": json_number(validation_results[region_id]["validation_score"]) if region_id in validation_results else None,
            "validation_region": validation_results[region_id]["validation_region"] if region_id in validation_results else None,
            "validation_mapping_confidence": validation_results[region_id]["mapping_confidence"] if region_id in validation_results else None,
            "agreement_status": validation_results[region_id]["agreement_status"] if region_id in validation_results else "not_measured",
        }
        records.append(record)
        positions, indices = parcel_mesh(atlas == region_id, atlas_image.affine)
        if not indices:
            raise ValueError(f"Atlas label {region_id} ({result['region_name']}) has no geometry")
        meshes.append({
            "region_id": region_id,
            "positions": positions,
            "indices": indices,
            "voxel_count": int(np.count_nonzero(atlas == region_id)),
        })

    enrichment = {
        "schema_version": "1.0.0",
        "disease_id": "parkinson-disease",
        "disease_name": "Parkinson disease",
        "gene_set_version": parameters["primary_gene_set_version"],
        "source_file": "data/results/parkinson_regional_enrichment.csv",
        "spatial_source_file": "data/results/parkinson_spatial_robustness.csv",
        "generated_on": date.today().isoformat(),
        "regions": records,
    }
    region_count, gene_count = expression_dimensions()
    donor_count, sample_count = ahba_counts()
    primary_studies = [row for row in read_csv(STUDY_SELECTION) if row["primary_selected"] == "True"]
    if len(primary_studies) != 1:
        raise ValueError("Expected exactly one primary GWAS in Stage 3 output")
    primary_study = primary_studies[0]
    credible_set_count = len(read_csv(CREDIBLE_SETS))
    z_bound = max(abs(record["z_score"]) for record in records)
    observed_values = [record["observed_score"] for record in records]
    global_spatial = [row for row in read_csv(GLOBAL_SPATIAL) if row["gene_set_version"] == "weighted"]
    if len(global_spatial) != 1:
        raise ValueError("Expected one weighted Stage 7 global result")
    validation_statistics = read_csv(VALIDATION_RESULTS)
    primary_validation = [row for row in validation_statistics if row["discovery_metric"] == "weighted_z_primary"]
    spatial_validation = [row for row in validation_statistics if row["discovery_metric"] == "weighted_z_spatial_moran_null"]
    if len(primary_validation) != 1 or len(spatial_validation) != 1:
        raise ValueError("Expected one primary and one spatial Stage 8 validation result")
    primary_validation = primary_validation[0]
    spatial_validation = spatial_validation[0]
    validation_values = [record["validation_score"] for record in records if record["validation_score"] is not None]
    validation_bound = max(abs(value) for value in validation_values)
    metadata_payload = {
        "schema_version": "1.0.0",
        "project": "GENE2BRAIN",
        "title": "From Genetic Risk to Spatial Brain Vulnerability",
        "generated_on": date.today().isoformat(),
        "available_diseases": [{"id": "parkinson-disease", "name": "Parkinson disease", "status": "available"}],
        "future_diseases": [
            "Alzheimer disease", "Huntington disease", "ALS", "Schizophrenia", "Bipolar disorder",
        ],
        "counts": {
            "ahba_donors": donor_count,
            "ahba_samples": sample_count,
            "atlas_parcels": len(metadata),
            "regions_analyzed": region_count,
            "genes_analyzed": gene_count,
            "parkinson_genes_represented": parameters["genes_per_set"]["weighted"],
            "random_gene_sets": parameters["n_permutations_per_gene_set"],
            "credible_sets": credible_set_count,
            "gwas_sample_size": int(primary_study["sample_size"]),
        },
        "analysis": {
            "primary_gene_set": parameters["primary_gene_set_version"],
            "gwas_accession": primary_study["study_accession"],
            "gwas_title": "Multi-ancestry genome-wide association meta-analysis of Parkinson's disease",
            "gene_prioritization": "Open Targets Platform Locus-to-Gene (L2G)",
            "null_model": parameters["matching_method"],
            "matching_variables": parameters["matching_variables"],
            "multiple_testing": parameters["multiple_testing"],
            "fdr_threshold": parameters["fdr_threshold"],
            "significant_regions": sum(record["fdr_p"] < parameters["fdr_threshold"] for record in records),
            "null_draws_available": False,
            "null_summary_available": True,
            "null_summary_note": "Per-region random mean and standard deviation are retained; individual permutation draws were not exported by Stage 6.",
            "spatial_sensitivity": {
                "method": spatial_parameters["spatial_null"],
                "neighbor_definition": spatial_parameters["neighbor_definition"],
                "permutations": spatial_parameters["n_spatial_permutations"],
                "global_statistic": spatial_parameters["global_test"],
                "global_p": json_number(global_spatial[0]["spatial_null_p"]),
                "robust_rule": spatial_parameters["robust_rule"],
                "robust_regions": sum(record["spatial_robustness_label"] == "robust" for record in records),
                "isolated_regions": len(spatial_parameters["isolated_region_ids"]),
            },
            "independent_validation": {
                "dataset": "ENIGMA-Parkinson's 2021 PD-versus-control MRI summary statistics",
                "phenotype": "Sign-reversed Cohen's d for cortical thickness or subcortical volume; higher means thinner/smaller in PD",
                "matched_validation_units": int(primary_validation["n_regions"]),
                "mapped_aal3_regions": len(validation_values),
                "pearson_r": json_number(primary_validation["pearson_r"]),
                "pearson_p": json_number(primary_validation["pearson_p"]),
                "pearson_ci": [json_number(primary_validation["pearson_ci_low"]), json_number(primary_validation["pearson_ci_high"])],
                "spearman_rho": json_number(primary_validation["spearman_rho"]),
                "spearman_p": json_number(primary_validation["spearman_p"]),
                "spatial_null_p": json_number(spatial_validation["pearson_p"]),
                "interpretation": "NOT SUPPORTED",
                "independence_note": "Independent measurement and analysis; exact participant-disjointness from the upstream GWAS cannot be guaranteed",
                "substantia_nigra": "Not measured in the public ENIGMA-PD FreeSurfer tables",
            },
        },
        "metrics": {
            "z_score": {"label": "Z-score", "domain": [-z_bound, z_bound], "scale": "diverging, symmetric about zero"},
            "observed_score": {"label": "Observed expression", "domain": [min(observed_values), max(observed_values)], "scale": "sequential"},
            "fdr_p": {"label": "FDR significance", "threshold": parameters["fdr_threshold"], "scale": "binary significance status"},
            "spatial_robustness": {"label": "Spatial robustness", "domain": [0, 1], "scale": "sequential spatial-null percentile"},
            "validation_score": {"label": "Independent validation", "domain": [-validation_bound, validation_bound], "scale": "diverging, symmetric about zero; gray is not measured"},
            "agreement": {"label": "Agreement", "scale": "four pre-specified median-split categories; gray is not measured"},
        },
        "atlas": {
            "name": "Automated Anatomical Labeling atlas 3 (AAL3v1)",
            "distribution": "AAL3v2 for SPM12, April 2024",
            "space": "MNI",
            "resolution_mm": [2, 2, 2],
            "dimensions": list(atlas.shape),
            "source_image": "data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz",
            "mapping": "region_id equals the nonzero integer voxel label in AAL3v1.nii.gz",
            "geometry_method": "Exposed voxel faces extracted from retained labels and greedily merged only when coplanar and contiguous",
            "render_coordinate_mapping": "MNI (x, y, z) is rendered as WebGL (x, z, y) so superior is up",
        },
        "sources": [
            {"name": "Allen Human Brain Atlas", "url": "https://human.brain-map.org/"},
            {"name": "abagen", "url": "https://abagen.readthedocs.io/"},
            {"name": "NHGRI-EBI GWAS Catalog", "url": "https://www.ebi.ac.uk/gwas/"},
            {"name": "Selected Parkinson GWAS GCST90308590", "url": "https://www.ebi.ac.uk/gwas/studies/GCST90308590"},
            {"name": "Open Targets Platform", "url": "https://platform.opentargets.org/"},
            {"name": "AAL3 atlas publication", "url": "https://doi.org/10.1016/j.neuroimage.2019.116189"},
            {"name": "ENIGMA-Parkinson's structural MRI validation", "url": "https://doi.org/10.1002/mds.28706"},
            {"name": "ENIGMA Toolbox summary statistics", "url": "https://enigma-toolbox.readthedocs.io/en/latest/pages/04.loadsumstats/"},
        ],
    }
    stage9_files = [STAGE9_GO, STAGE9_REACTOME, STAGE9_RANKED, STAGE9_CELLS, STAGE9_EVIDENCE]
    if all(path.is_file() for path in stage9_files):
        go_rows = read_csv(STAGE9_GO)
        reactome_rows = read_csv(STAGE9_REACTOME)
        ranked_rows = read_csv(STAGE9_RANKED)
        cell_rows = read_csv(STAGE9_CELLS)
        evidence_rows = read_csv(STAGE9_EVIDENCE)
        metadata_payload["analysis"]["biological_interpretation"] = {
            "status": "post-hoc interpretation of frozen discovery outputs",
            "primary_gene_set": "Stage 4 L2G-weighted",
            "genes_analyzed": 149,
            "genes_represented_in_ahba": parameters["genes_per_set"]["weighted"],
            "background": f"{gene_count:,} genes retained in the Stage 2 AHBA matrix",
            "go_significant_terms": sum(float(row["fdr"]) < 0.05 for row in go_rows),
            "reactome_significant_pathways": sum(float(row["fdr"]) < 0.05 for row in reactome_rows),
            "ranked_significant_pathways": sum(float(row["fdr"]) < 0.05 for row in ranked_rows),
            "cell_types_significant": sum(float(row["fdr"]) < 0.05 for row in cell_rows),
            "supported_region_pathway_rows": len(evidence_rows),
            "region_rule": "top 10 by frozen Stage 7 robustness rank",
            "interpretation_note": "Enrichment and cell-type associations do not establish causality, pathology, pathway activation, or clinical utility",
        }
        metadata_payload["sources"].extend([
            {"name": "g:Profiler functional enrichment", "url": "https://biit.cs.ut.ee/gprofiler/"},
            {"name": "Gene Ontology", "url": "https://geneontology.org/"},
            {"name": "Reactome release 97", "url": "https://reactome.org/"},
            {"name": "WikiPathways September 2026", "url": "https://data.wikipathways.org/current/gmt/"},
            {"name": "Human Protein Atlas v25.1 single-nucleus brain", "url": "https://www.proteinatlas.org/humanproteome/single+cell/single+nuclei+brain/data"},
        ])
    geometry = {
        "schema_version": "1.0.0",
        "atlas": "AAL3v1",
        "coordinate_system": "WebGL (MNI x, MNI z, MNI y), millimetres",
        "region_count": len(meshes),
        "regions": meshes,
    }

    write_json(WEB_DATA / "parkinson_brain_enrichment.json", enrichment)
    write_json(WEB_DATA / "project_metadata.json", metadata_payload)
    write_json(PUBLIC_DATA / "parkinson_brain_enrichment.json", enrichment)
    write_json(PUBLIC_DATA / "project_metadata.json", metadata_payload)
    write_json(PUBLIC_DATA / "aal3_regions.json", geometry, compact=True)
    shutil.copy2(RESULTS, PUBLIC_DATA / "parkinson_regional_enrichment.csv")
    shutil.copy2(SPATIAL_RESULTS, PUBLIC_DATA / "parkinson_spatial_robustness.csv")
    shutil.copy2(VALIDATION_REGIONAL, PUBLIC_DATA / "parkinson_independent_validation_regional_scores.csv")
    shutil.copy2(VALIDATION_RESULTS, PUBLIC_DATA / "parkinson_independent_validation_statistics.csv")
    if all(path.is_file() for path in stage9_files) and BIOLOGY_JSON.is_file():
        shutil.copy2(BIOLOGY_JSON, PUBLIC_DATA / BIOLOGY_JSON.name)
        for path in (STAGE9_GO, STAGE9_REACTOME, STAGE9_CELLS, STAGE9_EVIDENCE):
            shutil.copy2(path, PUBLIC_DATA / path.name)
    print(f"Prepared {len(records)} research records and {len(meshes)} atlas meshes")
    print(f"Geometry: {(PUBLIC_DATA / 'aal3_regions.json').stat().st_size / 1024 / 1024:.2f} MiB")


if __name__ == "__main__":
    build()
