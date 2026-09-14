# Stage 6 gene-set bias assessment

## Candidate confounders

The Stage 5 matrix contains scaled-robust-sigmoid normalized expression for 15,632
genes across 138 retained AAL3 regions. Before random sampling, Parkinson broad-set
genes were compared with the eligible AHBA background on global mean expression,
regional expression variance, reannotated candidate probe count, and the number of
regions with a finite measurement.

The matching model uses mean expression, log10 variance, and log1p probe count.
Each variable is standardized over the AHBA universe and given equal unit-variance
weight in Euclidean distance. These variables address baseline abundance, spatial
variability, and microarray representation. They summarize measurement properties
without matching the region-by-region pattern being tested.

Coverage was inspected but not used because Stage 2 complete-case filtering makes
it exactly 138 regions for every eligible gene. Gene length is not available in the
existing Stage 2 metadata and was not inferred or silently substituted. Candidate
probe count is an imperfect technical proxy: `abagen` exposes the retained gene,
while Stage 2 records reannotated candidate probes rather than claiming the exact
selected probe identity.

## Background and matching

All 123 AHBA-represented broad Parkinson genes were excluded from the common
background, leaving 15,509 genes. For each disease gene, the algorithm computes its
200 nearest background neighbors, randomizes disease-gene processing order within
each permutation, and draws an unused candidate. It expands beyond 200 only if a
pool is exhausted. This yields unique genes within every set while keeping the
procedure transparent. The same broad matched sets are reused for the weighted null.

## Balance results

| gene_set_version | variable | parkinson_mean | all_eligible_background_mean | matched_draw_mean | standardized_mean_difference_before | standardized_mean_difference_after | used_for_matching |
|---|---|---|---|---|---|---|---|
| broad | mean_expression | 0.5077 | 0.5051 | 0.5077 | 0.1178 | -3.011e-05 | True |
| broad | log10_expression_variance | -1.764 | -1.783 | -1.768 | 0.1419 | 0.02696 | True |
| broad | log1p_probe_count | 1.208 | 1.17 | 1.2 | 0.1471 | 0.03333 | True |
| broad | measurable_region_count | 138 | 138 | 138 | nan | nan | False |
| stringent | mean_expression | 0.5084 | 0.5051 | 0.5086 | 0.1494 | -0.009215 | True |
| stringent | log10_expression_variance | -1.752 | -1.783 | -1.754 | 0.2292 | 0.01232 | True |
| stringent | log1p_probe_count | 1.224 | 1.17 | 1.216 | 0.2072 | 0.02947 | True |
| stringent | measurable_region_count | 138 | 138 | 138 | nan | nan | False |
| weighted | mean_expression | 0.5077 | 0.5051 | 0.5077 | 0.1178 | -3.011e-05 | True |
| weighted | log10_expression_variance | -1.764 | -1.783 | -1.768 | 0.1419 | 0.02696 | True |
| weighted | log1p_probe_count | 1.208 | 1.17 | 1.2 | 0.1471 | 0.03333 | True |
| weighted | measurable_region_count | 138 | 138 | 138 | nan | nan | False |

Standardized mean differences after matching are diagnostics rather than outcome
tests. The figures show complete empirical distributions and a fixed 100,000-draw
subsample of the matched distribution for plotting efficiency.
