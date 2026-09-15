#!/usr/bin/env python3
"""Build public multi-disease JSON directly from Stage 10 result tables."""

from __future__ import annotations

from datetime import date
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results"
PROCESSED = ROOT / "data" / "processed"
DATA_WEB = ROOT / "data" / "web"
PUBLIC = ROOT / "web" / "public" / "data"


def finite(value: object) -> float | int | None:
    if value is None or pd.isna(value):
        return None
    number = float(value)
    return int(number) if number.is_integer() else number


def main() -> None:
    DATA_WEB.mkdir(parents=True, exist_ok=True)
    PUBLIC.mkdir(parents=True, exist_ok=True)
    enrichment = pd.read_csv(RESULTS / "multidisease_regional_enrichment.csv")
    enrichment = enrichment.loc[enrichment["gene_set"] == "weighted"]
    spatial = pd.read_csv(RESULTS / "multidisease_spatial_robustness.csv")
    spatial = spatial.loc[spatial["gene_set"] == "weighted"]
    raw = pd.read_csv(RESULTS / "multidisease_regional_raw_scores.csv")
    raw = raw.loc[raw["gene_set"] == "weighted"]
    regions = pd.read_csv(PROCESSED / "region_metadata.csv")
    regions = regions.loc[regions["retained_in_main_matrix"].astype(bool)].set_index("region_id")
    drivers = pd.read_csv(RESULTS / "multidisease_regional_gene_drivers.csv")
    spatial_weights = pd.read_csv(ROOT / "data" / "intermediate" / "brain_region_spatial_weights.csv", index_col="region_id")
    isolated = set(spatial_weights.index[spatial_weights.sum(axis=1) == 0].astype(int))
    existing_pd = json.loads((DATA_WEB / "parkinson_brain_enrichment.json").read_text(encoding="utf-8"))
    pd_validation = {row["region_id"]: row for row in existing_pd["regions"]}
    datasets = []
    for disease, disease_enrichment in enrichment.groupby("disease", sort=False):
        disease_id = disease_enrichment["disease_id"].iloc[0]
        e = disease_enrichment.set_index("region_id")
        s = spatial.loc[spatial["disease"] == disease].set_index("region_id")
        r = raw.loc[raw["disease"] == disease].set_index("region_id")
        disease_regions = []
        for region_id in regions.index:
            meta = regions.loc[region_id]
            validation = pd_validation.get(int(region_id), {}) if disease_id == "parkinson" else {}
            top = drivers.loc[
                (drivers["disease"] == disease) & (drivers["region_id"] == region_id) & drivers["top_driver"].astype(bool)
            ].nsmallest(10, "regional_rank")
            biology = {
                "region_id": int(region_id),
                "selected_for_regional_interpretation": not top.empty,
                "selection_rule": "top 10 regions by frozen spatial robustness rank; top 20 additive weighted-expression drivers",
                "robustness_rank": int(s.loc[region_id, "robustness_rank"]),
                "top_genes": [
                    {
                        "gene": row.gene, "expression": finite(row.expression),
                        "l2g_score": finite(row.l2g_score),
                        "weighted_contribution": finite(row.weighted_contribution),
                        "regional_rank": int(row.regional_rank),
                    }
                    for row in top.itertuples()
                ],
                "pathways": [], "cell_types": [],
            }
            # Preserve the independently validated Stage 9 Parkinson regional
            # interpretation where it exists. The Stage 10 table supplies the
            # same fields for newly analysed diseases, but intentionally only
            # for regions selected by the frozen regional-interpretation rule.
            if disease_id == "parkinson" and validation.get("biology"):
                biology = validation["biology"]
            disease_regions.append({
                "region_id": int(region_id), "region_name": meta["region_name"],
                "atlas_id": f"AAL3v1:{int(region_id)}", "atlas_label": str(meta["atlas_label"]),
                "hemisphere": meta["hemisphere"], "broad_system": meta["broad_system"],
                "centroid_mni": [finite(meta["centroid_mni_x"]), finite(meta["centroid_mni_y"]), finite(meta["centroid_mni_z"])],
                "observed_score": finite(e.loc[region_id, "observed_score"]),
                "random_mean": finite(e.loc[region_id, "random_mean"]),
                "random_std": finite(e.loc[region_id, "random_std"]),
                "z_score": finite(e.loc[region_id, "z_score"]),
                "empirical_p": finite(e.loc[region_id, "empirical_p"]),
                "fdr_p": finite(e.loc[region_id, "fdr_p"]),
                "effect_size": finite(e.loc[region_id, "effect_size"]),
                "number_of_genes": int(r.loc[region_id, "number_of_genes"]),
                "spatial_null_p": finite(s.loc[region_id, "spatial_null_p"]),
                "spatial_null_fdr": finite(s.loc[region_id, "spatial_null_fdr"]),
                "spatial_robustness": finite(s.loc[region_id, "spatial_robustness"]),
                "spatial_robustness_label": "robust" if bool(s.loc[region_id, "spatially_robust"]) else "not_robust",
                "robustness_rank": int(s.loc[region_id, "robustness_rank"]),
                "spatial_isolate": int(region_id) in isolated,
                "validation_score": validation.get("validation_score"),
                "validation_region": validation.get("validation_region"),
                "validation_mapping_confidence": validation.get("validation_mapping_confidence"),
                "agreement_status": validation.get("agreement_status", "not_measured"),
                "biology": biology,
            })
        datasets.append({
            "schema_version": "2.0", "disease_id": disease_id,
            "disease_name": disease, "gene_set_version": "weighted",
            "source_file": "data/results/multidisease_regional_enrichment.csv",
            "generated_on": date.today().isoformat(), "regions": disease_regions,
        })
    all_z = enrichment["z_score"].to_numpy(float)
    limit = float(np.quantile(np.abs(all_z), 0.99))
    config = json.loads((ROOT / "config" / "disease_panel.yaml").read_text(encoding="utf-8"))
    similarity = pd.read_csv(RESULTS / "disease_spatial_similarity.csv")
    similarity = similarity.loc[similarity["metric"] == "pearson"]
    payload = {
        "schema_version": "2.0", "generated_on": date.today().isoformat(),
        "primary_metric": "matched gene-set permutation Z score",
        "comparison_note": "All diseases share one global color domain and the frozen Parkinson matched-null method. A-B means Z_A minus Z_B at the same AAL3 parcel; it is descriptive, not a test of the difference.",
        "z_domain": [-limit, limit],
        "pearson_similarity": [
            {"disease_1": row.disease_1, "disease_2": row.disease_2,
             "correlation": finite(row.correlation), "p_value": finite(row.p_value),
             "n_regions": int(row.n_regions)}
            for row in similarity.itertuples()
        ],
        "diseases": datasets,
        "excluded_or_needs_review": [
            {"disease_id": item["disease_id"], "disease_name": item["disease_name"], "status": item["status"], "reason": item["selection_note"]}
            for item in config["diseases"] if item["disease_name"] not in {dataset["disease_name"] for dataset in datasets}
        ],
    }
    canonical = DATA_WEB / "multidisease_atlas.json"
    canonical.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")
    shutil.copy2(canonical, PUBLIC / canonical.name)
    for filename in (
        "disease_region_enrichment_matrix.csv", "disease_region_fdr_matrix.csv",
        "disease_region_robustness_matrix.csv", "disease_spatial_similarity.csv",
        "shared_brain_region_enrichment.csv", "disease_specific_regional_signatures.csv",
    ):
        shutil.copy2(RESULTS / filename, PUBLIC / filename)
    print(f"Built web atlas: {len(datasets)} diseases × {len(regions)} regions")


if __name__ == "__main__":
    main()
