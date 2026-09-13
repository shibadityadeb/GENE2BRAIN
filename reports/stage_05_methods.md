# Stage 5 methods: Parkinson × AHBA observed spatial signal

## Scope and interpretation

This stage describes the observed spatial expression of genes prioritised for
Parkinson disease in the healthy adult AHBA reference. A high score means only
that the selected Parkinson-associated genes have higher expression on the Stage 2
scale in that region. No random gene sets, null distributions, p-values, FDR,
spatial permutations, pathways, or other diseases were analysed.

## Inputs and identifier matching

- AHBA matrix: 138 retained AAL3 regions × 15,632 genes.
- Broad set: 149 total, 123 exact-symbol matches.
- Stringent set: 38 total, 34 exact-symbol matches.
- Weighted set: 149 total, 123 exact-symbol matches.
- Matching used the approved gene-symbol column from Stage 4 against the columns of
  `brain_region_gene_expression.csv`. `gene_metadata.csv` was checked to have the
  same 15,632-symbol universe. No identifier conversion was needed.
- Every included and missing identifier is recorded in
  `data/results/parkinson_gene_match_audit.csv`; nothing was silently discarded.

Missing broad/weighted genes (26): BRIP1, BST1, CCDC178, CCDC18, CRHR1, DEPDC1B, FCGR2B, FGF20, FOXA1, GCH1, GPR65, IGSF9B, IL1R2, IL1RL2, IQGAP2, KCNMB3, LCORL, MIX23, MMRN1, MYOC, NOD2, NUP42, OLAH, TTC6, WDR5B, WNT9B.

Missing stringent genes (4): BST1, IGSF9B, LCORL, NOD2.

## Expression scale

Stage 2 already applied scaled robust sigmoid normalization across genes within each
sample and then across matched samples for each gene within donor, separately for
cortex, subcortex/brainstem, and cerebellum. Samples were averaged within region and
donor, then available donors were equally averaged. Values therefore lie in [0, 1]
and are neither raw microarray intensities nor z-scores. The primary analysis,
labelled `raw_scale_analysis`, preserves this Stage 2 scale without another
transformation.

The separate `gene_standardized_analysis` sensitivity analysis transforms each gene
across the 138 pooled regions as

$$Z_{rg} = (E_{rg} - \bar E_g) / \sigma_g,$$

using the population standard deviation (`ddof=0`). It is not mixed with or used to
replace the primary score.

## Regional scores

For region $r$, broad genes $G_B$, and stringent genes $G_S$:

$$S_B(r) = \frac{1}{|G_B|} \sum_{g \in G_B} E_{rg},$$

$$S_S(r) = \frac{1}{|G_S|} \sum_{g \in G_S} E_{rg}.$$

For weighted genes $G_W$ and their unmodified Stage 4 maximum L2G score $w_g$:

$$S_W(r) = \frac{\sum_{g \in G_W} E_{rg} w_g}{\sum_{g \in G_W} w_g}.$$

Weights were not rescaled, normalized, thresholded again, or exponentiated. Missing
genes are excluded before both numerator and denominator are formed. Regions are
ranked independently for each score; ranking is descriptive.

## Gene contributions

For the ten highest weighted-score regions, the additive contribution of gene $g$ is
$E_{rg}w_g / \sum_g w_g$. These contributions sum to the region's weighted score.
The heatmap additionally divides each additive contribution by the regional score,
so its displayed fractions sum to one within each region.

## Donor robustness

The six Stage 2 donor matrices were scored with the same matched genes and weights.
Unobserved donor-region combinations remain missing; they were not imputed.
Pairwise Pearson correlations use only regions observed in both donors. Mean
pairwise correlations were broad=0.168,
stringent=0.321, and weighted=0.251.

## Sensitivity analysis

Raw-scale versus gene-standardized regional Pearson correlations were
broad=0.995, stringent=0.995, and
weighted=0.997. The long-form sensitivity file preserves
both profiles for every region and method.

## Limitations

The six post-mortem donors are few and unevenly sampled, four predominantly in the
left hemisphere. The pooled matrix retains only regions meeting the Stage 2 coverage
rule and is an adult healthy-brain reference. Exact-symbol matching leaves genes
without reliable AHBA representation out of the score. Broad anatomical-structure
normalization limits interpretation of absolute offsets between cortex,
subcortex/brainstem, and cerebellum. L2G is prioritisation evidence rather than a
definitive gene assignment. These observed rankings require an appropriate matched
null model before inferential interpretation.
