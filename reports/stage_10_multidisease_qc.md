# Stage 10 multi-disease quality control

This report was generated before spatial inference using the frozen panel and
thresholds. 10 diseases passed and 5 did not.

## Analyzed

- Parkinson disease
- Alzheimer disease
- Amyotrophic lateral sclerosis
- Schizophrenia
- Bipolar disorder
- Major depressive disorder
- Attention deficit-hyperactivity disorder
- Epilepsy
- Migraine
- Multiple sclerosis

## Not analyzed

- Huntington disease: Catalog ontology has only GCST004691 for Huntington progression, not common-variant disease risk.
- Multiple system atrophy: Eligible risk GWAS exists, but current Open Targets has zero credible sets; older eligible studies also have zero.
- Progressive supranuclear palsy: Largest disease-specific genome-wide array study with Open Targets credible sets, but only four prioritized genes were returned; it fails the frozen minimum-five-gene AHBA analysis gate.
- Autism spectrum disorder: The larger SPARK record has zero Open Targets credible sets; this disease-specific alternative returns only one prioritized gene and fails the frozen minimum-five-gene AHBA analysis gate.
- Tourette syndrome: Eligible risk GWAS exists, but current Open Targets has zero credible sets; older eligible studies also have zero.

Checks include gene-set size, exact AHBA symbol overlap, duplicate identifiers,
missing Ensembl identifiers, L2G distribution, ancestry, GWAS sample size, and
the frozen minimum-five-gene / 50%-coverage gate. Full values are in
`data/results/multidisease_qc_summary.csv` and
`data/results/multidisease_ahba_gene_coverage.csv`.
