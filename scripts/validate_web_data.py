#!/usr/bin/env python3
"""Fail loudly when public GENE2BRAIN records and AAL3 geometry diverge."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required web artifact is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def validate(data_dir: Path) -> list[str]:
    errors: list[str] = []
    enrichment = load(data_dir / "parkinson_brain_enrichment.json")
    project = load(data_dir / "project_metadata.json")
    geometry = load(data_dir / "aal3_regions.json")
    records = enrichment.get("regions", [])
    meshes = geometry.get("regions", [])
    record_ids = [row.get("region_id") for row in records]
    mesh_ids = [row.get("region_id") for row in meshes]

    if len(record_ids) != len(set(record_ids)):
        errors.append("duplicate region_id values in research records")
    if len(mesh_ids) != len(set(mesh_ids)):
        errors.append("duplicate region_id values in atlas geometry")
    if set(record_ids) != set(mesh_ids):
        errors.append(
            f"research/geometry ID mismatch: data-only={sorted(set(record_ids) - set(mesh_ids))}, "
            f"geometry-only={sorted(set(mesh_ids) - set(record_ids))}"
        )
    if len(records) != project.get("counts", {}).get("regions_analyzed"):
        errors.append("research record count disagrees with project metadata")

    required_text = ("region_name", "atlas_id", "hemisphere")
    required_numbers = (
        "observed_score", "random_mean", "random_std", "z_score", "empirical_p",
        "fdr_p", "effect_size", "number_of_genes",
    )
    for index, row in enumerate(records):
        prefix = f"record[{index}] region_id={row.get('region_id')}:"
        for field in required_text:
            if not isinstance(row.get(field), str) or not row[field].strip():
                errors.append(f"{prefix} missing/invalid {field}")
        for field in required_numbers:
            value = row.get(field)
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
                errors.append(f"{prefix} missing/non-finite {field}")
        for field in ("empirical_p", "fdr_p"):
            value = row.get(field)
            if isinstance(value, (int, float)) and not 0 <= value <= 1:
                errors.append(f"{prefix} {field} is outside [0, 1]")
        if isinstance(row.get("random_std"), (int, float)) and row["random_std"] < 0:
            errors.append(f"{prefix} random_std is negative")
        if isinstance(row.get("number_of_genes"), int) and row["number_of_genes"] <= 0:
            errors.append(f"{prefix} number_of_genes is not positive")

    for index, mesh in enumerate(meshes):
        positions = mesh.get("positions")
        indices = mesh.get("indices")
        prefix = f"geometry[{index}] region_id={mesh.get('region_id')}:"
        if not isinstance(positions, list) or len(positions) < 9 or len(positions) % 3:
            errors.append(f"{prefix} invalid positions")
            continue
        if not all(isinstance(value, (int, float)) and math.isfinite(value) for value in positions):
            errors.append(f"{prefix} non-finite position")
        if not isinstance(indices, list) or len(indices) < 3 or len(indices) % 3:
            errors.append(f"{prefix} invalid triangle indices")
        elif min(indices) < 0 or max(indices) >= len(positions) // 3:
            errors.append(f"{prefix} triangle index outside vertex array")

    threshold = project.get("analysis", {}).get("fdr_threshold")
    if not isinstance(threshold, (int, float)) or not 0 < threshold < 1:
        errors.append("invalid configured FDR threshold")
    expected_significant = sum(row["fdr_p"] < threshold for row in records) if isinstance(threshold, (int, float)) else -1
    if expected_significant != project.get("analysis", {}).get("significant_regions"):
        errors.append("configured significant-region count is inconsistent")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir", type=Path, default=ROOT / "web" / "public" / "data",
        help="Directory containing enrichment, metadata, and geometry JSON files",
    )
    args = parser.parse_args()
    errors = validate(args.data_dir)
    if errors:
        raise SystemExit("Web data validation failed:\n- " + "\n- ".join(errors))
    if args.data_dir.resolve() == (ROOT / "web" / "public" / "data").resolve():
        for filename in ("parkinson_brain_enrichment.json", "project_metadata.json"):
            canonical = ROOT / "data" / "web" / filename
            deployed = args.data_dir / filename
            if canonical.read_bytes() != deployed.read_bytes():
                raise SystemExit(f"Web data validation failed:\n- deployed {filename} differs from data/web source")
    print(f"Web data validation passed: {args.data_dir}")


if __name__ == "__main__":
    main()
