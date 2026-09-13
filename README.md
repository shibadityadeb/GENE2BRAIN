# From Genetic Risk to Spatial Brain Vulnerability

This repository supports a research pipeline testing whether disease-associated
genes are expressed unusually strongly in specific regions of the healthy human
brain. Stages 1–2 acquire the Allen Human Brain Atlas (AHBA) adult human
microarray data and build a documented healthy-brain region × gene reference.

## Current contents

- `notebooks/01_download_ahba.ipynb`: documented, executable AHBA download and QC
- `data/ahba/`: raw AHBA files downloaded by `abagen` (not tracked by Git)
- `data/processed/ahba_metadata_summary.csv`: concise donor/file inventory
- `data/processed/ahba_provenance.md`: source, assumptions, and package versions
- `results/figures/`: sampling coverage and donor sample-count plots
- `src/stage_02_preprocess.py`: reproducible AAL3 regional-expression workflow
- `data/processed/brain_region_gene_expression.csv`: final region × gene matrix
- `data/processed/region_metadata.csv`: all atlas parcels and coverage status
- `data/processed/gene_metadata.csv`: output-gene/probe mapping information
- `data/intermediate/donor_*_expression.csv`: generated donor-level regional
  matrices (kept local because they total about 151 MB)
- `reports/stage_02_methods.md` and `stage_02_qc.md`: exact methods and QC
- `.github/workflows/ci.yml`: reproducibility and artifact-integrity checks
- `src/stage_03_parkinson_gwas.py`: ontology-backed Parkinson GWAS Catalog
  retrieval, study ranking, association QC, and provisional locus grouping
- `data/gwas/parkinson_*.csv`: Stage 3 candidate studies, selected-study
  associations, provisional loci, mapped candidate genes, and summary-statistics
  metadata
- `reports/stage_03_*.md`: exact API audit, methods, and GWAS QC findings

Work currently stops at the Stage 3 Parkinson GWAS evidence layer. Catalog gene
mappings are labeled candidate/mapped genes, never causal genes. No rigorous
locus-to-gene prioritization, AHBA disease integration, enrichment, permutation,
or vulnerability mapping has been performed.

## Reproduce from a clean environment

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate       # Windows: .venv\\Scripts\\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=-1 notebooks/01_download_ahba.ipynb
python src/stage_02_preprocess.py
python src/stage_03_parkinson_gwas.py
python -m unittest discover -s tests -v
```

The notebook resolves all paths relative to the repository and calls the
currently documented API:

```python
abagen.fetch_microarray(
    donors="all",
    data_dir=PROJECT_ROOT / "data" / "ahba",
    resume=True,
    convert=False,
)
```

The complete six-donor archive requires roughly 4 GB plus temporary download
space. Downloads can be resumed.

The CI suite intentionally does not redownload the approximately 4 GB AHBA
archive. It validates the notebook, source code, atlas definition, committed
regional matrix, metadata, and QC outputs. A clean end-to-end data rebuild is
performed with the two commands above before scientific releases.

## Data source

The AHBA was created by the Allen Institute for Brain Science. Raw data are
retrieved from the Allen Brain Map well-known-file API by `abagen`; see the
generated provenance record for URLs, versions, and the access date.
