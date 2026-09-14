# GENE2BRAIN Stage 10 final multi-disease summary

## Scope

- Panel screened: 15 diseases
- Diseases analyzed: 10
- Excluded or held at needs-review: 5
- Brain regions: 138 AAL3 parcels
- AHBA genes: 15,632
- Primary comparison statistic: regional matched-gene-set permutation Z score
- Weighted-set regional BH-FDR discoveries: 9
- Jointly spatially robust discoveries: 0

## Disease results

| Disease | GWAS N | Loci | Prioritized genes | AHBA coverage % | Significant regions | Spatially robust | Validation | Top Z region |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| Parkinson disease | 2525730.0 | 67 | 149 | 82.5503355704698 | 0 | 0 | partially_validated | OFCant R |
| Alzheimer disease | 487511.0 | 83 | 113 | 82.30088495575221 | 0 | 0 | pending | Thal PuM L |
| Huntington disease | NA | 0 | 0 | NA | 0 | 0 | not_available | not analyzed |
| Amyotrophic lateral sclerosis | 152268.0 | 17 | 25 | 88.0 | 0 | 0 | pending | Occipital Mid R |
| Multiple system atrophy | 8016.0 | 0 | 0 | NA | 0 | 0 | not_available | not analyzed |
| Progressive supranuclear palsy | 12308.0 | 2 | 4 | 100.0 | 0 | 0 | not_available | not analyzed |
| Schizophrenia | 175799.0 | 222 | 451 | 84.4789356984479 | 7 | 0 | pending | Putamen R |
| Bipolar disorder | 413466.0 | 37 | 100 | 80.0 | 0 | 0 | pending | Temporal Mid L |
| Major depressive disorder | 1154267.0 | 94 | 178 | 81.46067415730337 | 0 | 0 | pending | Hippocampus L |
| Autism spectrum disorder | 46350.0 | 1 | 1 | 100.0 | 0 | 0 | not_available | not analyzed |
| Attention deficit-hyperactivity disorder | 225534.0 | 22 | 35 | 88.57142857142857 | 2 | 0 | pending | OFCmed R |
| Epilepsy | 82482.0 | 4 | 8 | 87.5 | 0 | 0 | pending | Occipital Inf R |
| Migraine | 873341.0 | 98 | 205 | 81.46341463414635 | 0 | 0 | pending | Cerebellum 8 R |
| Multiple sclerosis | 41505.0 | 140 | 297 | 73.06397306397307 | 0 | 0 | pending | Vermis 9 |
| Tourette syndrome | 19698.0 | 0 | 0 | NA | 0 | 0 | not_available | not analyzed |

## Quantitative cross-disease questions

1. **Distinct regional signatures.** Pairwise Pearson correlations span -0.527 to 0.593; median absolute correlation is 0.182. This quantifies heterogeneity but is not a disease-classification test.
2. **Most similar pair.** Amyotrophic lateral sclerosis and Epilepsy: Pearson r=0.593, nominal p=1.741e-14 across 138 regions.
3. **Most distinct pair.** Alzheimer disease and Schizophrenia: Pearson r=-0.527, nominal p=3.077e-11 across 138 regions.
4. **Repeated region.** OFCant R has 1 diseases with gene-set FDR significance and 0 with joint spatial robustness. Recurrence is not universal vulnerability.
5. **Disease-specific regions.** 11 disease-region cells pass the exploratory global BH-FDR test based on leave-one-disease-out specificity. The strongest is Schizophrenia in Hippocampus R (specificity Z=5.145, FDR=0.000126).
6. **Spatial versus biological similarity.** Across 45 disease pairs with at least one significant pathway set, spatial correlation versus pathway Jaccard similarity had Spearman rho=0.145, p=0.3417. This is exploratory and does not establish shared causal biology.
7. **Gene-set robustness.** Broad/stringent/weighted regional Spearman correlations range from 0.427 to 0.965 across diseases. Epilepsy's two-gene stringent set is explicitly unstable.

## GWAS power and gene-set size

- gwas_sample_size_vs_number_of_enriched_regions: Spearman rho=-0.182, p=0.6155
- gwas_sample_size_vs_maximum_z_score: Spearman rho=-0.200, p=0.5796
- gwas_sample_size_vs_number_of_spatially_robust_regions: not estimable (outcome constant)
- gene_set_size_vs_number_of_enriched_regions: Spearman rho=0.234, p=0.5161
- gene_set_size_vs_maximum_z_score: Spearman rho=0.527, p=0.1173
- gene_set_size_vs_number_of_spatially_robust_regions: not estimable (outcome constant)

These ten-disease tests have limited power and correlated regional observations.
They diagnose possible confounding; they do not prove that GWAS power causes the
spatial results.

## Biology

The standardized tables contain 88,086 GO/pathway rows and 1,020
HPA brain cell-type rows. Weighted-set BH-FDR discoveries total
1026 ontology/pathway
rows and 12 cell-type
rows. Correlated and nested ontology terms are not counted as independent
mechanisms.

## Negative-results statement

No disease produced a region passing both gene-set FDR and spatial-null FDR.
This null spatial-robustness result is retained. No threshold, GWAS, gene-set
size, matching variable, null model, or region list was changed to improve it.
Parkinson independent validation remains NOT SUPPORTED; other analyzed diseases
remain pending rather than receiving artificial validation data.

## Final visual package

1. `stage_10_gene2brain_pipeline.png`
2. existing Parkinson genetic enrichment figure from Stage 6
3. existing Parkinson independent-validation figure from Stage 8
4. `stage_10_multidisease_signature_matrix.png`
5. `stage_10_disease_spatial_similarity.png`
6. `stage_10_shared_vs_specific.png`
7. `stage_10_disease_pca.png`
8. `stage_10_biological_programs.png`

The atlas is a research visualization of healthy-brain expression enrichment.
It is not a clinical classifier, diagnostic tool, patient-level prediction, or
causal map.
