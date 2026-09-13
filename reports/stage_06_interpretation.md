# Stage 6 interpretation

## What this analysis means

For each region, the observed Stage 5 L2G-weighted Parkinson score was compared with
10,000 unique-gene random sets matched on global expression mean,
expression variance, and candidate probe count. A positive effect size or Z score
means the Parkinson-prioritised genes have higher expression than those matched
sets in that region. The one-sided empirical p-value includes the prespecified +1
correction, and Benjamini-Hochberg correction is applied across 138 regions at
q < 0.05. 0 regions pass that threshold in the primary weighted analysis.

This result is labelled **gene-set permutation significance**. It does not show that
a region causes Parkinson disease, initiates pathology, is specific to Parkinson,
or will be affected in an individual. L2G scores remain prioritisation evidence,
and healthy adult post-mortem expression is not a direct measurement of disease.

## Assumptions and limitations

- The selected matching summaries adequately control major measurable gene-level
  abundance, variance, and probe-representation differences.
- Matching does not control unrecorded properties such as gene length, GC content,
  cell-type specificity, network degree, or all probe-design effects.
- Excluding Parkinson genes from the background avoids direct contamination but
  slightly changes the eligible-gene universe.
- Fixed Parkinson L2G weights are carried to their gene-specific matched replacements,
  preserving the exact weight vector without inventing random weights.
- The six AHBA donors are few and unevenly sampled; four are predominantly left-sided.
- Benjamini-Hochberg correction addresses multiple regional tests under its standard
  assumptions, but the gene-set permutations do not remove spatial autocorrelation
  among neighboring atlas parcels. A spatially informed null remains a separate
  sensitivity analysis for later work.
- Results depend on the AAL3 parcellation, Stage 2 structural normalization, Stage 4
  gene thresholds, and Open Targets release.
