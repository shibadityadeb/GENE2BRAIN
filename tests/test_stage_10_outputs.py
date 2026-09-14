from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class Stage10Outputs(unittest.TestCase):
    def test_panel_and_frozen_parameters(self) -> None:
        panel = json.loads((ROOT / "config" / "disease_panel.yaml").read_text())
        params = json.loads((ROOT / "config" / "statistical_parameters.yaml").read_text())
        self.assertEqual(len(panel["diseases"]), 15)
        self.assertEqual(sum(item["status"] == "ready" for item in panel["diseases"]), 10)
        self.assertEqual(params["gene_set_permutations"], 10_000)
        self.assertEqual(params["spatial_permutations"], 10_000)
        self.assertEqual(params["fdr_threshold"], 0.05)

    def test_common_matrices(self) -> None:
        z = pd.read_csv(ROOT / "data" / "results" / "disease_region_enrichment_matrix.csv", index_col=0)
        fdr = pd.read_csv(ROOT / "data" / "results" / "disease_region_fdr_matrix.csv", index_col=0)
        robust = pd.read_csv(ROOT / "data" / "results" / "disease_region_robustness_matrix.csv", index_col=0)
        self.assertEqual(z.shape, (10, 138))
        self.assertTrue(z.index.equals(fdr.index) and z.columns.equals(fdr.columns))
        self.assertTrue(z.index.equals(robust.index) and z.columns.equals(robust.columns))
        self.assertTrue(np.isfinite(z.to_numpy()).all())
        self.assertTrue(((fdr >= 0) & (fdr <= 1)).all().all())
        self.assertTrue(robust.isin([0, 1]).all().all())

    def test_enrichment_and_spatial_contract(self) -> None:
        enrichment = pd.read_csv(ROOT / "data" / "results" / "multidisease_regional_enrichment.csv")
        spatial = pd.read_csv(ROOT / "data" / "results" / "multidisease_spatial_robustness.csv")
        self.assertEqual(len(enrichment), 10 * 3 * 138)
        self.assertEqual(len(spatial), 10 * 3 * 138)
        self.assertEqual(set(enrichment["n_permutations"]), {10_000})
        self.assertEqual(set(spatial["n_spatial_permutations"]), {10_000})
        self.assertEqual(int(spatial["spatially_robust"].sum()), 0)

    def test_gene_sets_and_qc_exclusions(self) -> None:
        qc = pd.read_csv(ROOT / "data" / "results" / "multidisease_qc_summary.csv")
        self.assertEqual((qc["analysis_status"] == "analyzable").sum(), 10)
        self.assertEqual(qc.loc[qc["disease"] == "Autism spectrum disorder", "analysis_status"].iloc[0], "needs_review_not_analyzed")
        self.assertEqual(qc.loc[qc["disease"] == "Progressive supranuclear palsy", "analysis_status"].iloc[0], "needs_review_not_analyzed")
        for disease in qc.loc[qc["analysis_status"] == "analyzable", "disease"]:
            disease_id = next(item["disease_id"] for item in json.loads((ROOT / "config" / "disease_panel.yaml").read_text())["diseases"] if item["disease_name"] == disease)
            for method in ("broad", "stringent", "weighted"):
                self.assertTrue((ROOT / "data" / "genes" / "multidisease" / f"{disease_id}_{method}.csv").is_file())

    def test_biology_and_manifest(self) -> None:
        pathways = pd.read_csv(ROOT / "data" / "results" / "multidisease_pathway_enrichment.csv")
        cells = pd.read_csv(ROOT / "data" / "results" / "multidisease_cell_type_enrichment.csv")
        manifest = pd.read_csv(ROOT / "data" / "results" / "gene2brain_data_manifest.csv")
        self.assertEqual(pathways["disease"].nunique(), 10)
        self.assertEqual(cells["disease"].nunique(), 10)
        self.assertGreater(len(manifest), 80)
        self.assertTrue((manifest["status"] == "complete").all())

    def test_web_atlas_contains_only_completed_real_results(self) -> None:
        atlas = json.loads((ROOT / "data" / "web" / "multidisease_atlas.json").read_text())
        self.assertEqual(len(atlas["diseases"]), 10)
        self.assertEqual(len(atlas["excluded_or_needs_review"]), 5)
        self.assertTrue(all(len(disease["regions"]) == 138 for disease in atlas["diseases"]))


if __name__ == "__main__":
    unittest.main()
