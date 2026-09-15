#!/usr/bin/env python3
"""Validate the anatomical surface → AAL3 ID → real research-record join."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "web/public/data/anatomical_brain.json"
CANONICAL = ROOT / "data/web/anatomical_brain.json"
PARKINSON = ROOT / "data/web/parkinson_brain_enrichment.json"
MULTIDISEASE = ROOT / "data/web/multidisease_atlas.json"
REGIONAL_CSV = ROOT / "data/results/parkinson_regional_enrichment.csv"


def fail(message: str) -> None:
    raise SystemExit(f"Brain region mapping validation failed: {message}")


def validate_mesh(mesh: dict, label: str) -> None:
    vertices = mesh.get("positions")
    triangles = mesh.get("indices")
    if not isinstance(vertices, list) or len(vertices) < 9 or len(vertices) % 3:
        fail(f"{label} has an invalid vertex array")
    if not isinstance(triangles, list) or len(triangles) < 3 or len(triangles) % 3:
        fail(f"{label} has an invalid triangle array")
    if any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in vertices):
        fail(f"{label} has a non-finite coordinate")
    if any(not isinstance(value, int) or value < 0 or value >= len(vertices) // 3 for value in triangles):
        fail(f"{label} has an invalid vertex reference")


def main() -> None:
    for path in (GEOMETRY, CANONICAL, PARKINSON, MULTIDISEASE, REGIONAL_CSV):
        if not path.is_file():
            fail(f"missing input {path.relative_to(ROOT)}")
    if GEOMETRY.read_bytes() != CANONICAL.read_bytes():
        fail("deployed anatomical geometry differs from its canonical copy")
    geometry = json.loads(GEOMETRY.read_text(encoding="utf-8"))
    parkinson = json.loads(PARKINSON.read_text(encoding="utf-8"))
    multi = json.loads(MULTIDISEASE.read_text(encoding="utf-8"))
    if geometry.get("schema_version") != "2.0.0" or "voxel exterior" in geometry.get("geometry_source", ""):
        fail("active geometry is not the anatomical surface release")
    validate_mesh(geometry.get("anatomy", {}), "continuous anatomical surface")
    mesh_rows = geometry.get("regions", [])
    mesh_ids = [row.get("region_id") for row in mesh_rows]
    if len(mesh_ids) != len(set(mesh_ids)) or len(mesh_ids) != geometry.get("region_count"):
        fail("geometry IDs are duplicated or the count is inconsistent")
    if any(row.get("source") not in (
        "fsaverage6 pial surface; nearest AAL3 registration-fusion label",
        "AAL3 label isosurface; 1 mm occupancy interpolation",
    ) for row in mesh_rows):
        fail("geometry contains an undocumented label mapping method")
    for row in mesh_rows:
        validate_mesh(row, f"AAL3 label {row['region_id']}")

    pd_rows = parkinson.get("regions", [])
    pd_ids = [row.get("region_id") for row in pd_rows]
    if len(pd_ids) != len(set(pd_ids)) or set(pd_ids) != set(mesh_ids):
        fail("Parkinson research records do not join one-to-one to anatomical regions")
    for disease in multi.get("diseases", []):
        ids = [row.get("region_id") for row in disease.get("regions", [])]
        if len(ids) != len(set(ids)) or set(ids) != set(mesh_ids):
            fail(f"{disease.get('disease_name')} does not join one-to-one to anatomical regions")
    with REGIONAL_CSV.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row.get("gene_set_version") == "weighted"]
    actual = {int(row["region_id"]): row for row in rows}
    if set(actual) != set(mesh_ids):
        fail("weighted Parkinson CSV IDs differ from anatomical regions")
    for record in pd_rows:
        source = actual[record["region_id"]]
        if record.get("atlas_id") != f"AAL3v1:{record['region_id']}":
            fail(f"invalid atlas ID for region {record['region_id']}")
        for field in ("z_score", "fdr_p", "observed_score", "random_mean", "effect_size"):
            if not math.isclose(float(source[field]), float(record[field]), rel_tol=1e-12, abs_tol=1e-12):
                fail(f"published Parkinson {field} differs from CSV for region {record['region_id']}")
    print(f"Anatomical mapping passed: {len(mesh_ids)} AAL3 regions; {len(multi['diseases'])} real disease datasets")


if __name__ == "__main__":
    main()
