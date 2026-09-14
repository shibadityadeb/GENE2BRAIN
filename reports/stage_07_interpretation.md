# Stage 07 interpretation

## What each stage measures

- **Stage 5 — raw expression:** regional averages of Parkinson-prioritized genes on
  the preserved AHBA normalization scale; no null-model inference.
- **Stage 6 — gene-set enrichment:** whether those genes score higher than technical-
  feature-matched random gene sets. This is the discovery analysis.
- **Stage 7 — spatial robustness:** whether the location and peak structure of the
  Stage 6 Z profile remain unusual against MSR maps with comparable atlas-graph
  spatial dependence. This is a sensitivity analysis, not a replacement result.

## Findings

The weighted Stage 6 Z profile had global Moran's I =
0.2131
(two-sided unrestricted-label permutation p =
9.999e-05). The spatially constrained
global peak test gave p = 0.5301.

Stage 6 contained 0 weighted regions at BH-FDR q < 0.05;
therefore 0 regions met the pre-specified joint Stage 6 + Stage 7 rule.
This negative joint result is retained explicitly. Continuous spatial-null
percentiles are provided for sensitivity and ranking, but they must not be presented
as corrected truth.

Broad/stringent/weighted spatial-robustness rank correlations ranged from
0.584
to 0.869.
The correlation between weighted spatial robustness and the median donor rank
percentile was 0.557. Donor summaries are descriptive because six donors,
uneven sampling, and donor missingness do not support donor-level significance.

## Interpretation boundaries

Parkinson-associated genes show a regional spatial pattern, but no region passed the
current joint robustness rule. This does not show where Parkinson starts, identify a
causal region, map pathology, prove selective vulnerability, or provide clinical
prediction. MSR depends on the chosen AAL3 contact graph, assumes the graph is an
adequate proxy for spatial dependence. Singleton MSR exactly preserves the graph
power spectrum but may be conservative because randomized maps can remain correlated
with the observed map. Four isolated brainstem parcels have weaker spatial constraints.
Independent pathological, longitudinal, and replication evidence would
be required for stronger claims.

Runtime for this build: 7.3 seconds using Python 3.12.7,
NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.3, and
nibabel 5.3.2.
