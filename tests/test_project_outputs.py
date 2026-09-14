"""Fast integrity checks for committed AHBA and Stage 2 outputs.

The raw AHBA archive is intentionally not downloaded in CI because it is about
4 GB. These tests verify the reproducible code, notebook, atlas, principal
region-by-gene matrix, metadata, and quality-control artifacts committed to the
repository.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "results" / "figures"
ATLAS = ROOT / "data" / "atlases" / "aal_3v2" / "AAL3" / "AAL3v1.nii.gz"


class ProjectStructureTests(unittest.TestCase):
    def test_required_project_files_exist(self) -> None:
        required = [
            ROOT / "README.md",
            ROOT / "requirements.txt",
            ROOT / "notebooks" / "01_download_ahba.ipynb",
            ROOT / "src" / "stage_02_preprocess.py",
            PROCESSED / "ahba_metadata_summary.csv",
            PROCESSED / "ahba_provenance.md",
            PROCESSED / "atlas_selection.md",
            ROOT / "reports" / "stage_02_methods.md",
            ROOT / "reports" / "stage_02_qc.md",
            FIGURES / "ahba_sampling_coverage.png",
            FIGURES / "ahba_samples_by_donor.png",
            FIGURES / "stage_02_donor_concordance.png",
            FIGURES / "stage_02_normal_transcriptomic_landscape.png",
            FIGURES / "stage_02_regional_transcriptomic_structure.png",
        ]
        missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
        self.assertEqual(missing, [], f"Missing required files: {missing}")

    def test_download_notebook_is_valid_json(self) -> None:
        path = ROOT / "notebooks" / "01_download_ahba.ipynb"
        with path.open(encoding="utf-8") as stream:
            notebook = json.load(stream)
        self.assertEqual(notebook.get("nbformat"), 4)
        self.assertGreater(len(notebook.get("cells", [])), 0)


class ScientificArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.expression = pd.read_csv(
            PROCESSED / "brain_region_gene_expression.csv", index_col="region_id"
        )
        cls.regions = pd.read_csv(PROCESSED / "region_metadata.csv")
        cls.genes = pd.read_csv(PROCESSED / "gene_metadata.csv")

    def test_aal3_atlas_has_expected_parcels(self) -> None:
        image = nib.load(ATLAS)
        labels = np.unique(np.asanyarray(image.dataobj).astype(np.int16))
        non_background = labels[labels != 0]
        self.assertEqual(len(non_background), 166)
        self.assertEqual(image.shape, (91, 109, 91))

    def test_main_matrix_dimensions_and_values(self) -> None:
        self.assertEqual(self.expression.shape, (138, 15_632))
        self.assertTrue(self.expression.index.is_unique)
        self.assertTrue(self.expression.columns.is_unique)
        self.assertTrue(np.isfinite(self.expression.to_numpy()).all())

    def test_gene_metadata_matches_matrix(self) -> None:
        self.assertEqual(len(self.genes), 15_632)
        self.assertTrue(self.genes["gene_symbol"].is_unique)
        self.assertEqual(set(self.genes["gene_symbol"]), set(self.expression.columns))

    def test_region_metadata_matches_matrix(self) -> None:
        self.assertEqual(len(self.regions), 166)
        retained = self.regions.loc[self.regions["retained_in_main_matrix"]]
        self.assertEqual(len(retained), 138)
        self.assertEqual(set(retained["region_id"]), set(self.expression.index))

    def test_donor_correlation_matrix(self) -> None:
        correlations = pd.read_csv(
            PROCESSED / "donor_correlation_matrix.csv", index_col="donor"
        )
        self.assertEqual(correlations.shape, (6, 6))
        self.assertTrue(np.allclose(correlations, correlations.T))
        self.assertTrue(np.allclose(np.diag(correlations), 1.0))

    def test_key_structures_are_represented(self) -> None:
        coverage = pd.read_csv(PROCESSED / "key_region_coverage.csv")
        self.assertEqual(len(coverage), 14)
        self.assertTrue(coverage["present_in_atlas"].all())
        self.assertTrue(coverage["represented_in_main_matrix"].all())


class ParkinsonGwasArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gwas = ROOT / "data" / "gwas"
        cls.candidates = pd.read_csv(cls.gwas / "parkinson_candidate_studies.csv")
        cls.selection = pd.read_csv(cls.gwas / "parkinson_study_selection.csv")
        cls.associations = pd.read_csv(cls.gwas / "parkinson_gwas_associations_raw.csv")
        cls.loci = pd.read_csv(cls.gwas / "parkinson_initial_loci.csv")
        cls.mapped = pd.read_csv(cls.gwas / "parkinson_gwas_mapped_genes.csv")

    def test_required_stage_03_outputs_exist(self) -> None:
        required = [
            self.gwas / "parkinson_candidate_studies.csv",
            self.gwas / "parkinson_study_selection.csv",
            self.gwas / "parkinson_gwas_associations_raw.csv",
            self.gwas / "parkinson_initial_loci.csv",
            self.gwas / "parkinson_gwas_mapped_genes.csv",
            self.gwas / "parkinson_summary_statistics_metadata.csv",
            ROOT / "reports" / "stage_03_gwas_api.md",
            ROOT / "reports" / "stage_03_gwas_qc.md",
            ROOT / "reports" / "stage_03_methods.md",
            FIGURES / "stage_03_parkinson_gwas_manhattan.png",
            FIGURES / "stage_03_parkinson_locus_overview.png",
            FIGURES / "stage_03_gwas_to_gene_concept.png",
        ]
        missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
        self.assertEqual(missing, [], f"Missing Stage 3 files: {missing}")

    def test_ontology_candidate_set_and_primary_selection(self) -> None:
        self.assertEqual(len(self.candidates), 113)
        self.assertTrue(self.candidates["study_accession"].is_unique)
        self.assertEqual(set(self.candidates["ontology_id"]), {"MONDO_0005180"})
        selected = self.selection.loc[self.selection["primary_selected"]]
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected.iloc[0]["study_accession"], "GCST90308590")
        self.assertEqual(int(selected.iloc[0]["study_rank"]), 1)

    def test_significant_associations_are_valid(self) -> None:
        self.assertEqual(len(self.associations), 109)
        self.assertTrue(self.associations["variant"].is_unique)
        self.assertTrue(self.associations["rsid"].str.fullmatch(r"rs\d+").all())
        self.assertTrue(self.associations["chromosome"].between(1, 22).all())
        self.assertTrue(self.associations["position"].gt(0).all())
        self.assertTrue(self.associations["p_value"].between(0, 5e-8).all())
        self.assertEqual(set(self.associations["study_accession"]), {"GCST90308590"})

    def test_provisional_loci_cover_every_association(self) -> None:
        self.assertEqual(len(self.loci), 78)
        self.assertEqual(int(self.loci["number_of_associations"].sum()), 109)
        self.assertFalse(self.associations["locus_id"].isna().any())
        self.assertEqual(set(self.associations["locus_id"]), set(self.loci["locus_id"]))

    def test_mapped_genes_remain_candidate_annotations(self) -> None:
        self.assertEqual(self.mapped["gene"].nunique(), 763)
        self.assertTrue(set(self.mapped["locus_id"]).issubset(set(self.loci["locus_id"])))
        self.assertTrue(self.mapped["mapping_source"].str.contains("not causal").all())

    def test_open_companion_summary_statistics_are_documented(self) -> None:
        metadata = pd.read_csv(
            self.gwas / "parkinson_summary_statistics_metadata.csv"
        )
        primary = metadata.loc[metadata["study_accession"] == "GCST90308590"].iloc[0]
        self.assertTrue(bool(primary["download_available"]))
        self.assertEqual(primary["source_accession"], "GCST90275127")
        self.assertEqual(primary["genome_build"], "GRCh37")
        self.assertEqual(int(primary["sample_size"]), 611_485)
        validation = metadata.loc[metadata["study_accession"] == "GCST90480008"].iloc[0]
        self.assertEqual(validation["genome_build"], "GRCh38")
        self.assertEqual(validation["file_size"], "827 MB")


class ParkinsonLocusToGeneArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gwas = ROOT / "data" / "gwas"
        cls.gene_dir = ROOT / "data" / "genes"
        cls.result_data = ROOT / "data" / "results"
        cls.credible = pd.read_csv(cls.gwas / "parkinson_credible_sets.csv")
        cls.l2g = pd.read_csv(cls.gene_dir / "parkinson_l2g_predictions.csv")
        cls.broad = pd.read_csv(cls.gene_dir / "parkinson_genes_broad.csv")
        cls.stringent = pd.read_csv(cls.gene_dir / "parkinson_genes_stringent.csv")
        cls.weighted = pd.read_csv(cls.gene_dir / "parkinson_genes_weighted.csv")

    def test_required_stage_04_outputs_exist(self) -> None:
        required = [
            self.gwas / "parkinson_opentargets_study_match.csv",
            self.gwas / "parkinson_credible_sets.csv",
            self.gwas / "opentargets_query_log.txt",
            self.gene_dir / "parkinson_l2g_predictions.csv",
            self.gene_dir / "parkinson_genes_broad.csv",
            self.gene_dir / "parkinson_genes_stringent.csv",
            self.gene_dir / "parkinson_genes_weighted.csv",
            self.gene_dir / "parkinson_gene_locus_evidence.csv",
            self.result_data / "parkinson_gwas_vs_l2g_gene_comparison.csv",
            self.result_data / "parkinson_gene_ahba_coverage.csv",
            ROOT / "results" / "tables" / "parkinson_top_prioritized_genes.csv",
            ROOT / "reports" / "stage_04_opentargets_method.md",
            ROOT / "reports" / "stage_04_methods.md",
            FIGURES / "stage_04_parkinson_gwas_to_gene.png",
            FIGURES / "stage_04_parkinson_l2g_distribution.png",
            FIGURES / "stage_04_parkinson_gene_mapping_comparison.png",
        ]
        missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
        self.assertEqual(missing, [], f"Missing Stage 4 files: {missing}")

    def test_study_and_credible_sets_are_release_consistent(self) -> None:
        match = pd.read_csv(self.gwas / "parkinson_opentargets_study_match.csv")
        self.assertEqual(len(match), 1)
        self.assertEqual(match.iloc[0]["opentargets_study"], "GCST90308590")
        self.assertEqual(match.iloc[0]["match_status"], "matched")
        self.assertEqual(len(self.credible), 67)
        self.assertTrue(self.credible["credible_set_id"].is_unique)
        self.assertEqual(set(self.credible["fine_mapping_method"]), {"PICS"})
        self.assertTrue(self.credible["stage_03_locus_id"].notna().all())
        self.assertTrue(self.credible["lead_posterior_inclusion_probability"].between(0, 1).all())

    def test_l2g_scores_and_gene_sets_preserve_evidence(self) -> None:
        self.assertEqual(len(self.l2g), 149)
        self.assertFalse(self.l2g.duplicated(["credible_set_id", "gene"]).any())
        self.assertTrue(self.l2g["l2g_score"].between(0, 1).all())
        self.assertTrue(self.l2g["l2g_score"].gt(0.05).all())
        self.assertEqual(set(self.stringent["gene"]) - set(self.broad["gene"]), set())
        self.assertEqual(set(self.weighted["gene"]), set(self.broad["gene"]))
        self.assertTrue(np.allclose(self.weighted["gene_weight"], self.weighted["max_l2g_score"]))

    def test_master_and_top_tables_have_required_semantics(self) -> None:
        master = pd.read_csv(self.gene_dir / "parkinson_gene_locus_evidence.csv")
        required = {
            "disease", "gene", "ensembl_id", "locus_id", "credible_set_id",
            "lead_variant", "chromosome", "position", "p_value_if_available",
            "l2g_score", "fine_mapping_method", "credible_set_confidence",
            "study", "evidence_source",
        }
        self.assertTrue(required.issubset(master.columns))
        self.assertEqual(len(master), len(self.l2g))
        top = pd.read_csv(ROOT / "results" / "tables" / "parkinson_top_prioritized_genes.csv")
        self.assertEqual(top["rank"].tolist(), list(range(1, len(top) + 1)))
        self.assertTrue(top["l2g_score"].is_monotonic_decreasing)
        self.assertTrue(top["evidence"].str.contains("not proof of causality").all())

    def test_ahba_check_is_coverage_only(self) -> None:
        coverage = pd.read_csv(self.result_data / "parkinson_gene_ahba_coverage.csv")
        self.assertEqual(set(coverage["gene_set"]), {"broad", "stringent", "weighted"})
        self.assertTrue(coverage["coverage_only_no_expression_scoring"].all())
        self.assertTrue(coverage["coverage_percentage"].between(0, 100).all())


class ParkinsonSpatialSignalArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result_data = ROOT / "data" / "results"
        cls.tables = ROOT / "results" / "tables"
        cls.brain_maps = ROOT / "results" / "brain_maps"
        cls.raw = pd.read_csv(cls.result_data / "parkinson_regional_raw_scores.csv")
        cls.expression = pd.read_csv(PROCESSED / "brain_region_gene_expression.csv", index_col="region_id")
        cls.weighted_genes = pd.read_csv(ROOT / "data" / "genes" / "parkinson_genes_weighted.csv")

    def test_required_stage_05_outputs_exist(self) -> None:
        required = [
            self.result_data / "parkinson_regional_raw_scores.csv",
            self.result_data / "parkinson_gene_match_audit.csv",
            self.result_data / "parkinson_top_region_gene_contributions.csv",
            self.result_data / "parkinson_donor_regional_scores.csv",
            self.result_data / "parkinson_donor_pattern_correlations.csv",
            self.result_data / "parkinson_score_method_correlations.csv",
            self.result_data / "parkinson_expression_scale_sensitivity.csv",
            ROOT / "reports" / "stage_05_methods.md",
            self.tables / "parkinson_top_regions_broad.csv",
            self.tables / "parkinson_top_regions_stringent.csv",
            self.tables / "parkinson_top_regions_weighted.csv",
            self.brain_maps / "stage_05_parkinson_broad_expression.png",
            self.brain_maps / "stage_05_parkinson_stringent_expression.png",
            self.brain_maps / "stage_05_parkinson_weighted_expression.png",
            FIGURES / "stage_05_parkinson_top_regions.png",
            FIGURES / "stage_05_parkinson_donor_concordance.png",
            FIGURES / "stage_05_parkinson_gene_contributions.png",
            FIGURES / "stage_05_parkinson_method_comparison.png",
        ]
        missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
        self.assertEqual(missing, [], f"Missing Stage 5 files: {missing}")

    def test_primary_scores_have_expected_dimensions_and_scale(self) -> None:
        self.assertEqual(len(self.raw), 138)
        self.assertTrue(self.raw["region_id"].is_unique)
        self.assertEqual(set(self.raw["analysis_scale"]), {"raw_scale_analysis"})
        score_columns = [
            "broad_mean_expression", "stringent_mean_expression", "weighted_mean_expression"
        ]
        self.assertTrue(np.isfinite(self.raw[score_columns].to_numpy()).all())
        self.assertTrue(self.raw[score_columns].apply(lambda x: x.between(0, 1).all()).all())
        self.assertEqual(set(self.raw["n_broad_genes_present"]), {123})
        self.assertEqual(set(self.raw["n_stringent_genes_present"]), {34})
        self.assertEqual(set(self.raw["n_weighted_genes_present"]), {123})

    def test_weighted_formula_is_reproducible(self) -> None:
        weights = self.weighted_genes.set_index("gene")["gene_weight"]
        genes = [gene for gene in weights.index if gene in self.expression.columns]
        expected = self.expression[genes].mul(weights.loc[genes], axis=1).sum(axis=1) / weights.loc[genes].sum()
        observed = self.raw.set_index("region_id")["weighted_mean_expression"]
        self.assertTrue(np.allclose(expected.sort_index(), observed.sort_index()))

    def test_rankings_are_descending_and_match_primary_scores(self) -> None:
        mapping = {
            "broad": "broad_mean_expression",
            "stringent": "stringent_mean_expression",
            "weighted": "weighted_mean_expression",
        }
        raw = self.raw.set_index("region_id")
        for method, column in mapping.items():
            ranking = pd.read_csv(self.tables / f"parkinson_top_regions_{method}.csv")
            self.assertEqual(len(ranking), 25)
            self.assertEqual(ranking["rank"].tolist(), list(range(1, 26)))
            self.assertTrue(ranking["score"].is_monotonic_decreasing)
            self.assertTrue(np.allclose(ranking["score"], raw.loc[ranking["region_id"], column]))

    def test_contributions_sum_to_weighted_scores(self) -> None:
        contributions = pd.read_csv(
            self.result_data / "parkinson_top_region_gene_contributions.csv"
        )
        grouped = contributions.groupby("region_id")
        self.assertEqual(grouped.ngroups, 10)
        self.assertTrue(np.allclose(grouped["weighted_contribution"].sum(), grouped["regional_weighted_score"].first()))
        self.assertTrue(np.allclose(grouped["normalized_contribution"].sum(), 1.0))

    def test_donor_and_sensitivity_outputs_are_complete(self) -> None:
        donor = pd.read_csv(self.result_data / "parkinson_donor_regional_scores.csv")
        correlations = pd.read_csv(self.result_data / "parkinson_donor_pattern_correlations.csv")
        self.assertEqual(donor["donor"].nunique(), 6)
        self.assertEqual(len(donor), 6 * 138)
        self.assertEqual(len(correlations), 3 * 15)
        self.assertTrue(correlations["pearson_r"].between(-1, 1).all())
        sensitivity = pd.read_csv(self.result_data / "parkinson_expression_scale_sensitivity.csv")
        self.assertEqual(len(sensitivity), 3 * 138)
        self.assertEqual(set(sensitivity["method"]), {"broad", "stringent", "weighted"})
        self.assertEqual(set(sensitivity["raw_scale_label"]), {"raw_scale_analysis"})
        self.assertEqual(set(sensitivity["standardized_scale_label"]), {"gene_standardized_analysis"})

    def test_method_correlations_are_symmetric(self) -> None:
        correlations = pd.read_csv(
            self.result_data / "parkinson_score_method_correlations.csv", index_col="method"
        )
        self.assertEqual(correlations.shape, (3, 3))
        self.assertTrue(np.allclose(correlations, correlations.T))
        self.assertTrue(np.allclose(np.diag(correlations), 1.0))


class ParkinsonEnrichmentArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result_data = ROOT / "data" / "results"
        cls.brain_maps = ROOT / "results" / "brain_maps"
        cls.primary = pd.read_csv(cls.result_data / "parkinson_regional_enrichment.csv")
        cls.all_methods = pd.read_csv(
            cls.result_data / "parkinson_regional_enrichment_all_methods.csv"
        )

    def test_required_stage_06_outputs_exist(self) -> None:
        required = [
            self.result_data / "parkinson_regional_enrichment.csv",
            self.result_data / "parkinson_regional_enrichment_all_methods.csv",
            self.result_data / "parkinson_null_model_method_comparison.csv",
            self.result_data / "stage_06_gene_matching_balance.csv",
            self.result_data / "stage_06_permutation_parameters.json",
            ROOT / "results" / "tables" / "parkinson_top_enriched_regions.csv",
            ROOT / "reports" / "stage_06_gene_set_bias_assessment.md",
            ROOT / "reports" / "stage_06_interpretation.md",
            ROOT / "reports" / "stage_06_performance.md",
            self.brain_maps / "stage_06_parkinson_zscore_map.png",
            self.brain_maps / "stage_06_parkinson_fdr_map.png",
            FIGURES / "stage_06_gene_matching_expression.png",
            FIGURES / "stage_06_gene_matching_coverage.png",
            FIGURES / "stage_06_parkinson_null_distributions.png",
            FIGURES / "stage_06_observed_vs_random.png",
            FIGURES / "stage_06_parkinson_enrichment_ranking.png",
            FIGURES / "stage_06_method_robustness_heatmap.png",
            FIGURES / "stage_06_donor_robustness.png",
        ]
        missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
        self.assertEqual(missing, [], f"Missing Stage 6 files: {missing}")

    def test_primary_empirical_statistics_are_reproducible(self) -> None:
        self.assertEqual(len(self.primary), 138)
        self.assertEqual(set(self.primary["gene_set_version"]), {"weighted"})
        self.assertEqual(set(self.primary["n_permutations"]), {10_000})
        expected_p = (self.primary["number_random_equal_or_greater"] + 1) / 10_001
        self.assertTrue(np.allclose(self.primary["empirical_p"], expected_p))
        self.assertTrue(
            np.allclose(
                self.primary["effect_size"],
                self.primary["observed_score"] - self.primary["random_mean"],
            )
        )
        expected_z = self.primary["effect_size"] / self.primary["random_std"]
        self.assertTrue(np.allclose(self.primary["z_score"], expected_z))
        expected_fdr = multipletests(self.primary["empirical_p"], method="fdr_bh")[1]
        self.assertTrue(np.allclose(self.primary["fdr_p"], expected_fdr))

    def test_all_gene_set_versions_have_complete_region_tests(self) -> None:
        self.assertEqual(len(self.all_methods), 3 * 138)
        self.assertEqual(
            self.all_methods.groupby("gene_set_version")["region_id"].nunique().to_dict(),
            {"broad": 138, "stringent": 138, "weighted": 138},
        )
        self.assertTrue(self.all_methods["empirical_p"].between(1 / 10_001, 1).all())
        self.assertTrue(self.all_methods["fdr_p"].between(0, 1).all())
        self.assertEqual(int((self.primary["fdr_p"] < 0.05).sum()), 0)

    def test_matching_parameters_and_balance_are_valid(self) -> None:
        with (self.result_data / "stage_06_permutation_parameters.json").open() as stream:
            parameters = json.load(stream)
        self.assertEqual(parameters["random_seed"], 20260913)
        self.assertEqual(parameters["n_permutations_per_gene_set"], 10_000)
        self.assertEqual(parameters["genes_per_set"], {"broad": 123, "stringent": 34, "weighted": 123})
        self.assertTrue(all(parameters["uniqueness_verified"].values()))
        balance = pd.read_csv(self.result_data / "stage_06_gene_matching_balance.csv")
        self.assertEqual(set(balance["gene_set_version"]), {"broad", "stringent", "weighted"})
        matched = balance.loc[balance["used_for_matching"]]
        self.assertTrue(matched["standardized_mean_difference_after"].abs().lt(0.1).all())

    def test_primary_and_comparison_tables_are_well_formed(self) -> None:
        ordered = self.primary.sort_values(
            ["fdr_p", "z_score"], ascending=[True, False]
        ).reset_index(drop=True)
        pd.testing.assert_frame_equal(self.primary, ordered)
        top = pd.read_csv(ROOT / "results" / "tables" / "parkinson_top_enriched_regions.csv")
        pd.testing.assert_frame_equal(top, self.primary.head(25))
        comparison = pd.read_csv(
            self.result_data / "parkinson_null_model_method_comparison.csv"
        )
        self.assertEqual(len(comparison), 3)
        self.assertTrue(comparison["pearson_r"].between(-1, 1).all())
        self.assertTrue(comparison["spearman_rho"].between(-1, 1).all())

    def test_interpretation_records_spatial_null_limitation(self) -> None:
        text = (ROOT / "reports" / "stage_06_interpretation.md").read_text()
        self.assertIn("gene-set permutation significance", text)
        self.assertIn("do not remove spatial autocorrelation", text)


if __name__ == "__main__":
    unittest.main()
