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


if __name__ == "__main__":
    unittest.main()
