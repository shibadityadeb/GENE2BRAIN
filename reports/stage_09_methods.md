# Stage 9 methods: biological interpretation

## Analysis boundary

Stage 9 is post-discovery interpretation. Stages 3--8 were frozen and were not
rerun or tuned. The pre-specified primary gene set is the Stage 4 L2G-weighted
set; broad and stringent sets are sensitivity analyses. L2G is genetic
prioritization evidence and is not interpreted as an expression effect or
causal probability.

## Identifiers and background

Original Stage 4 symbols and Ensembl identifiers were preserved. The custom
background is the 15,632
genes retained in the Stage 2 AHBA expression matrix: the genes that could
realistically enter spatial scoring. Unrepresented genes remain in the mapping
audit but are excluded by the custom universe rather than silently discarded.
Stage 4 queries use stable Ensembl identifiers; the symbol-indexed AHBA
background is normalized by g:Profiler. HPA records sharing a display symbol
are collapsed by maximum nCPM within cluster type before symbol-based joining.

## GO and pathway over-representation

g:Profiler `e114_eg62_p19_27110d83` was accessed programmatically on 2026-09-14. Sources
were GO Biological Process, Molecular Function, Cellular Component, Reactome,
and WikiPathways. Queries used the custom AHBA background. Raw one-sided
hypergeometric probabilities were recomputed from returned query, overlap,
term and domain counts. Terms required 5--1000
background genes and at least 1 overlapping query gene. Benjamini--Hochberg
FDR was applied separately to GO:BP, GO:MF, GO:CC, Reactome and WikiPathways.
GO terms are overlapping and hierarchical; FDR does not make them independent.
Reactome release 97 was current at access.

## Ranked pathway analysis

For each eligible Reactome and WikiPathways term, a one-sided Mann--Whitney
rank-sum test asks whether member genes have larger Stage 4 L2G scores than
other prioritized genes. AUC is the rank effect (0.5 under no shift). FDR is
separate by pathway database. This ranks genetic prioritization evidence; it
is not GSEA of differential expression.

## Cell types

Human Protein Atlas v25.1 single-nucleus brain data (34 cluster types across
11 brain regions; Siletti et al. source data) were downloaded on 2026-09-14.
Within the AHBA background, markers require nCPM >= 1 and log2(expression + 1
over mean-other-types + 1) > 1; the top 200 per cell type by
specificity form marker sets. Hypergeometric ORA and BH FDR were run as one
34-test family per gene-set definition. Categories are retained under HPA
names; no dopaminergic label was manufactured where the resource does not
provide one.

## Regions and drivers

The top 10 parcels by frozen Stage 7 robustness rank were specified
before interpretation. A contribution is `AHBA expression × L2G / sum(L2G)`.
The top 20 per parcel define regional contributors.
"Spatially contributing prioritized genes" are ranked by the number of these
parcels in which they enter that top-20 list, then mean
contribution. They are not called causal drivers.

Reactome ORA was run on each parcel's top contributors. BH correction across
all eligible region--term pairs is the regional family used for support;
within-region FDR is retained descriptively. Fewer than 10 genes would be
reported as insufficient rather than analyzed.

## Network decision

A STRING network was not added. Network centrality would answer a different,
optional question and could invite causal overinterpretation without changing
the requested pathway, cell-type or regional conclusions.

## Multiple testing families

- GO: separate GO:BP, GO:MF and GO:CC families.
- Pathways: separate Reactome and WikiPathways families, per frozen gene set.
- Ranked pathways: separate Reactome and WikiPathways families.
- Cell types: 34 HPA cluster types, separately per frozen gene set.
- Regional pathways: all eligible top-region × Reactome-term tests together.

## Reproducibility and provenance

Raw API requests/responses and official HPA downloads are stored under
`data/biology/raw`. HPA SHA-256: `a287c1abab593c9d6802349a82f118efc0d4aea411be7d897b8feb0b31b02f2a`.
Python 3.12.7, pandas 2.2.3, NumPy 1.26.4,
SciPy 1.13.1, requests 2.32.3, matplotlib 3.9.4.
