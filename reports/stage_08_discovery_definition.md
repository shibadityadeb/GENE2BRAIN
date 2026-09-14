# Stage 08 frozen discovery definition

This definition was recorded before candidate validation outcomes were obtained or
tested. Stage 8 does not rerun, tune, filter, or otherwise modify Stages 1–7.

## Primary discovery statistic

The primary GENE2BRAIN discovery statistic is the **Stage 6 L2G-weighted regional
enrichment Z-score** (`z_score`, `gene_set_version=weighted`) for every one of the
138 retained AAL3 regions in
`data/results/parkinson_regional_enrichment.csv`. Higher values represent stronger
expression of the frozen Parkinson-prioritized gene set relative to technical-
feature-matched random gene sets. All validation regions that can be mapped with
the pre-specified anatomical rules will be included, irrespective of direction or
influence on the validation association.

The primary validation test is the association between this weighted Z-score and
the independent regional phenotype. Pearson correlation is the primary effect-size
summary; Spearman rank correlation is co-primary for robustness to scale and
outliers. Two-sided p-values, 95% confidence intervals, and the number of matched
regions will be reported.

## Pre-specified secondary discovery metrics

1. Stage 6 weighted FDR rank, oriented so larger values mean stronger evidence.
2. Stage 7 spatial robustness (`1 - spatial_null_p`).
3. Stage 6 broad-gene-set Z-score.
4. Stage 6 stringent-gene-set Z-score.
5. Stage 6 L2G-weighted Z-score (repeated in the sensitivity table for direct
   comparability across gene-set definitions).

The weighted Z-score remains primary even if another metric correlates more strongly.
Stage 7 produced no regions satisfying its predeclared joint FDR rule, so validation
will use its continuous spatial-robustness measure as a sensitivity statistic and
will not relabel any region as robust.

## Pre-specified source-selection criteria

Candidate validation datasets will be ranked without examining their correlation
with GENE2BRAIN. The preferred source must provide: (1) quantitative Parkinson-versus-
control regional measurements; (2) a cohort independent of the GWAS and AHBA donors;
(3) numerical data released by the study or an official repository; (4) documented
anatomical definitions that permit rule-based mapping to AAL3; (5) broad bilateral
coverage where possible; (6) adequate sample size; and (7) reproducible provenance.

When several sources qualify, priority is: anatomical coverage, cohort size,
machine-readable numerical availability, mapping clarity, and measurement
reliability. The sign of every phenotype will be oriented before testing so larger
values consistently mean greater Parkinson-related vulnerability. A source will not
be selected because it yields a favorable association.

## Pre-specified analyses

- Primary Pearson and Spearman correlations use all reliably mapped regions.
- Correlation confidence intervals use a fixed-seed nonparametric region bootstrap.
- Leave-one-region-out analysis uses every matched region without selective removal.
- Substantia nigra is retained in the primary analysis and excluded only in a named
  sensitivity analysis if the selected dataset measures it.
- Potential confounds will be assessed only when the required covariate is available
  for the matched regions. Region size, baseline AHBA expression, broad anatomical
  class, represented-gene count, spatial dependence, and phenotype reliability will
  be reported separately rather than used for outcome-driven model selection.
- Agreement categories use median splits fixed here: high/high, high/low, low/high,
  and low/low after each vector is oriented so higher means greater predicted or
  observed vulnerability. Values equal to the median enter the high group.

No discovery statistic, source-selection criterion, anatomical mapping rule, or
agreement threshold may be changed in response to the validation result.
