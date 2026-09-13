# Stage 4 methods and quality control

## Scope

Stage 4 matched the primary Parkinson GWAS `GCST90308590` to Open Targets, retrieved
fine-mapped credible sets and all displayed L2G predictions, and generated broad,
stringent, and weighted gene sets. No AHBA regional scores, enrichment tests,
permutations, pathways, or disease brain maps were computed.

## Why locus-to-gene prioritisation is required

A GWAS association identifies a statistical signal at a variant, not a gene. Lead
variants can merely tag nearby correlated variants through linkage disequilibrium,
and regulatory variants can act on genes that are not the nearest gene. Fine mapping
uses association statistics and LD to distribute posterior probability across a set
of plausible variants; the resulting credible set is intended to contain the causal
variant at a stated probability, rather than declaring its lead variant causal.

Open Targets L2G applies a gradient-boosting model to rank protein-coding genes near
each credible set. Its inputs include posterior-weighted gene/TSS distance,
colocalisation with eQTL, sQTL and pQTL signals, predicted variant consequences,
enhancer-to-gene evidence, and neighbourhood/context features. Multiple genes are
retained because the model supplies relative evidence, not a definitive causal call.

## Stage 3 consistency gate

The audit confirmed 109 genome-wide significant Stage 3 associations in 78
provisional loci on GRCh38.p14, with each reported lead variant equal to the minimum
p-value association in its locus. All 887 mapped-gene rows refer to valid loci (763
unique symbols). The selected study is `GCST90308590`, sample size 2,525,730, with
European, East Asian, Hispanic or Latin American, and African-unspecified discovery
ancestries. No inconsistency triggered the stop condition.

## Open Targets study match and retrieval

The Platform record has the exact accession, Parkinson's disease trait, sample size
2,525,730, and matching ancestry structure. Open Targets reports GRCh38 coordinates;
Stage 3 records GRCh38.p14, a patch-level label for the same major assembly.
67 of 67 credible-set lead
variants matched Stage 3 associations by rsID or coordinate. The Platform returned
67 PICS credible sets and 149 gene–credible-set L2G predictions
for 149 unique genes.

## Gene-set construction

- Broad: 149 genes with a returned score > 0.05.
- Stringent: 38 genes with maximum score ≥ 0.615029, the eligible-score Q3.
- Weighted: 149 genes; weight is the gene's maximum score across credible sets.
- Multiple candidates were retained per locus. The full feature values and SHAP
  contributions are serialized in the master table, with separate distance,
  colocalisation, variant-effect, and other-evidence fields.

## AHBA identifier coverage (lookup only)

| Gene set | Genes | In AHBA | Absent | Coverage |
|---|---:|---:|---:|---:|
| broad | 149 | 123 | 26 | 82.6% |
| stringent | 38 | 34 | 4 | 89.5% |
| weighted | 149 | 123 | 26 | 82.6% |

This is only an identifier-coverage check against the Stage 2 AHBA gene metadata;
expression values were not scored or aggregated by disease.

Genes were matched by exact approved gene symbol against the columns of
`data/processed/brain_region_gene_expression.csv`. This deliberately avoids alias
substitution at this stage; absent symbols may reflect microarray coverage,
reannotation, or symbol-version differences rather than absence of brain expression.

Three sets were retained so later work can test sensitivity to inclusion: broad
maximises coverage under the Platform's displayed-evidence rule, stringent focuses
on the upper quartile for this release, and weighted retains broad coverage without
treating a 0.90 and 0.08 score as equivalent.

## Caveats

All study credible sets have PICS/top-hit confidence and quality-control flags noting
the absence of in-sample LD. The Platform reports that harmonised summary statistics
for this study are unavailable or empty. `locusStart` and `locusEnd` may consequently
be absent; variant-level coordinates and posterior probabilities are preserved in
JSON. Stage 3 distance-grouped provisional loci and Open Targets credible sets are
different analytical objects and are linked only when their lead variant matches.
