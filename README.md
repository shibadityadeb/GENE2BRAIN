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
- `src/stage_04_parkinson_l2g.py`: current Open Targets Platform credible-set
  and Locus-to-Gene retrieval, thresholding, comparison, and AHBA symbol coverage
- `data/genes/parkinson_*.csv`: full L2G evidence plus broad, stringent, and
  continuous-weight gene sets
- `reports/stage_04_*.md`: Open Targets release/query details, methods, caveats,
  and Stage 4 quality control

Work currently stops at the Stage 4 Parkinson genetic evidence layer. L2G scores
are explicitly treated as prioritisation evidence, never proof of causality. No
AHBA disease scoring, enrichment, permutation, pathway analysis, or vulnerability
mapping has been performed.

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
python src/stage_04_parkinson_l2g.py
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
regional matrix, metadata, and Stage 2–4 scientific artifacts. Stage 4 makes
targeted calls to the current Open Targets GraphQL API when explicitly rebuilt;
CI validates the release-stamped committed results without relying on network
availability. A clean end-to-end data rebuild is performed with the commands
above before scientific releases.

## Data source

The AHBA was created by the Allen Institute for Brain Science. Raw data are
retrieved from the Allen Brain Map well-known-file API by `abagen`; see the
generated provenance record for URLs, versions, and the access date.
