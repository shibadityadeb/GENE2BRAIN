# AHBA data provenance

- **Source:** Allen Institute for Brain Science, Allen Human Brain Atlas
- **Project page:** https://human.brain-map.org/
- **Dataset name:** Allen Human Brain Atlas normalized microarray expression dataset
- **Access date:** 2026-09-13
- **Number of donors:** 6
- **Donors:** 9861, 10021, 12876, 14380, 15496, 15697
- **Data type:** Postmortem adult human brain microarray expression intensities, present/absent calls, probe annotations, sample annotations with MNI coordinates, and anatomical ontology
- **Genes directly represented in `Probes.csv`:** 29131 unique nonblank gene symbols (no reannotation applied)
- **Preprocessing status:** Allen-distributed normalized microarray files only. `convert=False` prevents parquet conversion. No probe reannotation/filtering, sample filtering, normalization by this project, sample-to-region assignment, donor aggregation, or regional aggregation has been performed.
- **Download method:** `abagen.fetch_microarray(donors="all", resume=True, convert=False, n_proc=1)`
- **Download location:** `data/ahba/microarray/`
- **Assumptions:** Donor identifiers and download IDs are those shipped with the installed `abagen` release. `mni_x`, `mni_y`, and `mni_z` in `SampleAnnot.csv` are treated as MNI-space millimetre coordinates. The headerless expression and PA-call matrix columns follow the sample-row order supplied in the same donor bundle; their dimensions are checked against `SampleAnnot.csv`. Gene count describes raw annotation symbols, not a curated analysis-ready gene set.

## Allen well-known-file URLs

  - Donor 9861: https://human.brain-map.org/api/v2/well_known_file_download/178238387
  - Donor 10021: https://human.brain-map.org/api/v2/well_known_file_download/178238373
  - Donor 12876: https://human.brain-map.org/api/v2/well_known_file_download/178238359
  - Donor 14380: https://human.brain-map.org/api/v2/well_known_file_download/178238316
  - Donor 15496: https://human.brain-map.org/api/v2/well_known_file_download/178238266
  - Donor 15697: https://human.brain-map.org/api/v2/well_known_file_download/178236545

## Exact direct package versions

  - `abagen==0.1.3`
  - `pandas==2.2.3`
  - `numpy==1.26.4`
  - `scipy==1.13.1`
  - `nibabel==5.3.2`
  - `nilearn==0.11.1`
  - `matplotlib==3.9.4`
  - `seaborn==0.13.2`
  - `statsmodels==0.14.4`
  - `scikit-learn==1.6.1`
  - `requests==2.32.3`
  - `jupyter==1.1.1`
  - `setuptools==75.8.2`
