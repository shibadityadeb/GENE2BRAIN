# Stage 10 biological interpretation

- Access date: 2026-09-14
- Diseases analyzed: 10
- Pathway/ontology rows: 88,086
- Weighted-set ontology/pathway rows at BH-FDR < 0.05: 1026
- Cell-type rows: 1,020
- Weighted-set HPA cell-type rows at BH-FDR < 0.05: 12
- Driver-gene summary rows: 344
- Disease-pair spatial/biological comparisons: 45

Functional enrichment uses g:Profiler with GO biological process, molecular
function, cellular component, Reactome, and WikiPathways sources and the frozen
15,632-gene AHBA background. BH correction is performed separately by ontology
inside each disease and gene-set request, matching Stage 9. Cell-type enrichment
uses the frozen HPA v25.1 human-brain single-nucleus marker table and a
hypergeometric test with BH correction within disease and gene-set version.

Regional drivers are additive `AHBA expression × L2G / sum(L2G)` contributions
in the ten regions selected by the already-computed spatial robustness rank.
Because no disease has a jointly spatially robust parcel, these driver results
are explicitly ranked exploratory interpretation, not regional discovery.

The spatial-versus-biological table is descriptive. It does not test or claim
that spatially similar diseases share causal biology, particularly when both
diseases have no FDR-significant biological terms.
