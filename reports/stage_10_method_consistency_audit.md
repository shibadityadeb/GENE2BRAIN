# Stage 10 method consistency audit

All 10 analyzed diseases used the same frozen AHBA matrix (138 AAL3
regions × 15,632 genes), region definitions, raw-scale Stage 5 equations,
exact-symbol matching, Stage 6 matching covariates (mean expression, log10
expression variance, log1p reannotated probe count), 200-nearest-neighbor unique
sampling, 10,000 gene-set permutations, one-sided +1 empirical p-values,
Benjamini-Hochberg correction, binary 26-neighbor AAL3 graph, Moran singleton
spectral randomization, 10,000 spatial permutations, and FDR 0.05.

The weighted matched null reuses broad matched replacements and carries the
fixed disease L2G weights, exactly as for Parkinson. Cross-disease comparisons
use the primary matched-null regional Z score; raw expression is never compared
between diseases.

## Exceptions

- Parkinson uses its frozen Stage 3–9 source artifacts without modification.
- Major depression uses the same-publication European component because the
  larger bi-ancestry Catalog record has zero Open Targets credible sets.
- Autism and progressive supranuclear palsy have real L2G outputs but fail the
  prespecified five-gene QC gate and are not spatially analyzed.
- Epilepsy passes the primary weighted-set gate with seven AHBA genes, but its
  two-gene stringent sensitivity result is explicitly flagged as unstable.
- Huntington, multiple system atrophy, and Tourette syndrome fail genetic/L2G
  inclusion and remain documented exclusions.
- Independent validation type is allowed to differ by disease. Only Parkinson
  has a completed Stage 8 analysis; its outcome remains NOT SUPPORTED. Other
  analyzed diseases remain pending rather than receiving proxy data.
