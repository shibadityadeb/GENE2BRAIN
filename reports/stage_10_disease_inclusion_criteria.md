# Stage 10 disease inclusion criteria

## Frozen decision point

These criteria were fixed on 2026-09-14 before any multi-disease AHBA regional score, enrichment statistic, spatial-null result, pathway result, or website map was calculated. Parkinson remains the reference pipeline. Disease-specific thresholds will not be tuned to obtain positive results.

## Required criteria

A disease is `ready` only when all of the following are satisfied:

1. The NHGRI-EBI GWAS Catalog contains a reproducible, variant-level, genome-wide study of disease risk. Progression, age-at-onset, treatment, imaging, cross-trait, gene-burden-only, chromosome-only, and proxy-only phenotypes are not primary-eligible.
2. The phenotype is sufficiently specific to represent the configured disease and the study reports interpretable sample and ancestry metadata.
3. The study is available in the current Open Targets Platform and has at least one fine-mapped credible set with displayed Locus-to-Gene predictions above the Platform's release filter of 0.05.
4. Variant coordinates and study metadata are compatible with the Open Targets GRCh38 representation. No silent liftover or assembly mixing is permitted.
5. A broad, stringent, and L2G-weighted gene set can be produced by the unchanged Stage 4 rules.
6. At least five prioritized genes and at least 50% of the selected genes must overlap the frozen AHBA matrix before spatial inference. Failure of this downstream check changes the analysis status to `needs_review`; it does not trigger a new GWAS or relaxed gene threshold.

Independent regional validation is desirable but is not required for genetic inclusion. Its availability and independence are classified separately. A missing validation dataset is preserved as a negative availability result, not replaced with a non-independent correlate.

## Consistent primary-study selection

Within each ontology-linked candidate pool, selection prioritizes a direct disease-risk phenotype, disease-focused consortium/meta-analysis design, adequate sample size, reproducible accession, and usable Open Targets credible-set/L2G evidence. Larger generic biobank, cross-trait, phenotype-proxy, chromosome-specific, or method-development records do not automatically supersede a dedicated disease GWAS. When the largest otherwise eligible record has no Open Targets credible sets, a same-phenotype alternative may be selected only if the exception is explicit in `config/disease_panel.yaml`.

## Pre-analysis panel decision

- **Ready after gene-set QC (10):** Parkinson disease, Alzheimer disease, amyotrophic lateral sclerosis, schizophrenia, bipolar disorder, major depressive disorder, attention deficit-hyperactivity disorder, epilepsy, migraine, and multiple sclerosis.
- **Needs review after gene-set QC (2):** progressive supranuclear palsy (four prioritized genes) and autism spectrum disorder (one prioritized gene). Their real source results are retained, but they do not enter spatial inference.
- **Excluded (3):** Huntington disease (only a progression GWAS is ontology-linked), multiple system atrophy (eligible risk GWAS but no Open Targets credible sets), and Tourette syndrome (eligible risk GWAS but no Open Targets credible sets).

Exclusion is about compatibility with this specific common-variant GWAS-to-L2G framework. It is not evidence that a disease lacks genetic architecture or biological importance.

## Frozen analysis parameters

All diseases load the thresholds, matching variables, permutation counts, FDR rule, and robustness rule from `config/statistical_parameters.yaml`. No disease-specific override is permitted. The cross-disease primary statistic is the matched-gene-set permutation Z score, not raw expression, because it standardizes each observed regional score against gene sets matched using the Parkinson reference procedure.

## Sources audited

- NHGRI-EBI GWAS Catalog REST API: ontology-linked studies, accessions, publication metadata, sample descriptions, ancestry, association availability, and full-p-value-set flag.
- Open Targets Platform GraphQL API, Platform 26.06 / API 26.6.3 at the audit date: study identity, GRCh38 credible sets, fine-mapping metadata, and Locus-to-Gene predictions.

The generated candidate and primary tables retain the live metadata snapshot used for the final decision. API evolution can change availability on a future rerun; such changes require a versioned panel amendment, not silent reselection.
