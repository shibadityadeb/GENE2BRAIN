# Stage 07 robustness definition

This rule was fixed before Stage 7 results were generated.

A region is called **robust** only if all three conditions hold:

1. its Stage 6 weighted enrichment Z-score is positive;
2. its Stage 6 gene-set permutation BH-FDR is below 0.05; and
3. its one-sided Stage 7 MSR regional p-value passes BH-FDR across 138 regions at q < 0.05.

The Stage 7 regional p-value is `(1 + number of surrogate values >= observed) /
(10000 + 1)` at the same fixed atlas parcel. `spatial_robustness` is the
descriptive percentile-like quantity `1 - spatial_null_p`; larger values indicate
that the observed positive score lies higher in that parcel's spatial-null
distribution. It is not a posterior probability. `robustness_rank` sorts this
measure descending with deterministic minimum ranks.

The joint rule prevents a spatial sensitivity result from reviving a region that
did not pass the discovery-stage multiplicity criterion. It also means a negative
Stage 6 result must remain negative even if its descriptive spatial rank is high.
