# Stage 08 validation-source review

## Decision

ENIGMA-Parkinson's 2021 was selected before any GENE2BRAIN correlation was
calculated. It supplies official machine-readable, hemisphere-specific standardized
PD-versus-control effects for cortical thickness and subcortical volume across 19
sites (2,357 people with PD and 1,182 controls).[^1] The official ENIGMA Toolbox
distributes the exact summary tables used here.[^2]

The selection prioritizes direct quantitative phenotype measurement, broad bilateral
coverage, sample size, standardized effect sizes, explicit FreeSurfer anatomy, and
reproducible numerical access. It does not use the sign or significance of the
GENE2BRAIN association.

## Alternatives and exclusions

Zeighami et al. provide a public whole-brain PPMI map and strong anatomical coverage,
including substantia nigra.[^3] However, the public voxel values are normalized-mixture
ICA component Z-scores, not regional PD-control effects; the paper explicitly makes
that distinction. Its direct regional atrophy t vector is not in the released source
tables. The Frigerio et al. postmortem study is biologically independent and directly
quantifies synaptic loss and pathology in eight cortical regions, but its supplements
do not release reusable region-level numerical observations; digitizing figures was
excluded.[^4] The ENIGMA cerebellar study is large but anatomically restricted,[^5]
and ENIGMA-DTI measures white-matter tracts that cannot be mapped to AAL3 gray-matter
parcels without an additional connectivity model.[^6]

## Independence qualification

The ENIGMA phenotype is independent of AHBA, the Stage 6 permutation model, Stage 7,
and the GENE2BRAIN analysis. It was measured with structural MRI and generated without
GENE2BRAIN outputs. Exact participant-disjointness from the upstream PD GWAS cannot be
guaranteed, because PPMI is one ENIGMA imaging site and PD GWAS consortia have used
PPMI genetic data. This analysis therefore supports an **external measurement and
analysis validation**, not an unqualified participant-independent replication. The
limitation is retained in every interpretation and prevents a `SUPPORTED` label.

## Sources

[^1]: Laansma et al. “[An International Multicenter Analysis of Brain Structure Across Clinical Stages of Parkinson's Disease](https://doi.org/10.1002/mds.28706).” *Movement Disorders* (2021).
[^2]: MICA-MNI. “[ENIGMA Toolbox: load summary statistics](https://enigma-toolbox.readthedocs.io/en/latest/pages/04.loadsumstats/).”
[^3]: Zeighami et al. “[Network structure of brain atrophy in de novo Parkinson's disease](https://doi.org/10.7554/eLife.08440).” *eLife* (2015); [NeuroVault collection 860](https://neurovault.org/collections/860/).
[^4]: Frigerio et al. “[Regional differences in synaptic degeneration](https://doi.org/10.1186/s40478-023-01711-w).” *Acta Neuropathologica Communications* (2024).
[^5]: Kerestes et al. “[Cerebellar Volume and Disease Staging in Parkinson's Disease](https://doi.org/10.1002/mds.29611).” *Movement Disorders* (2023).
[^6]: Owens-Walton et al. “[A worldwide study of white matter microstructural alterations](https://pubmed.ncbi.nlm.nih.gov/39128907/).” (2024).
