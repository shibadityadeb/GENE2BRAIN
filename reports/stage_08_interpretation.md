# Stage 08 interpretation

## Discovery and validation are separate

**Discovery** is the frozen Stage 6 L2G-weighted Parkinson gene-set enrichment
Z-score derived from GWAS prioritization and healthy-brain AHBA expression. **Validation**
is the sign-reversed ENIGMA-PD standardized PD-control difference in cortical thickness
or subcortical volume. ENIGMA values never altered the gene set, discovery scores,
permutations, atlas retention, or Stage 7 statistics.

## Primary result

Across 70 unique mapped ENIGMA parcels, Pearson r was
-0.028 (95% region-bootstrap CI -0.217 to
0.176, two-sided p=0.8186); Spearman rho
was -0.015 (95% bootstrap CI -0.247 to
0.219, p=0.903). The spatially
constrained Moran-surrogate p-value was 0.8313.

The broad and stringent gene-set Pearson correlations were 0.139
and -0.094. Cortex-only and subcortex-only estimates were
-0.022 (n=56) and -0.283
(n=14), respectively. Leave-one-region-out Pearson r ranged
from -0.072 to 0.003, with median
-0.028.

## Substantia nigra

The ENIGMA-PD 2021 public FreeSurfer tables do not measure substantia nigra. An
including-versus-excluding-SN correlation is therefore not estimable for this
phenotype. SN was not imputed, substituted from another outcome, or used to change
the primary analysis.

## Controls and limitations

Region size, baseline AHBA expression, cortical/subcortical class, spatial dependence,
represented-gene count, and validation reliability were evaluated as pre-specified.
The represented-gene count is constant at 123 genes in every discovery parcel
and cannot explain between-region covariance. Region-specific analyzed sample size is
used only as a descriptive reliability weight because it is not a measurement-error
variance. AAL3-to-Desikan-Killiany mapping aggregates discovery parcels once per
ENIGMA unit and avoids pseudo-replication, but medium-confidence composite mappings
remain approximate.

The result is classified **NOT SUPPORTED** using the complete evidence, not a single
p-value. The validation is independent in measurement and analysis, but exact
participant-disjointness from the upstream GWAS cannot be guaranteed. A positive
association would support convergence between healthy-brain genetic-expression
enrichment and cross-sectional PD morphometric vulnerability. It would not establish
causality, disease origin, temporal direction, cell-type mechanism, or individual
clinical prediction. A null association is retained as informative and does not
trigger discovery-model revision.

## References

1. Laansma et al. “[An International Multicenter Analysis of Brain Structure Across Clinical Stages of Parkinson's Disease](https://doi.org/10.1002/mds.28706).” *Movement Disorders* (2021).
2. MICA-MNI. “[ENIGMA Toolbox summary-statistics documentation](https://enigma-toolbox.readthedocs.io/en/latest/pages/04.loadsumstats/).”
3. Kim et al. “[Multi-ancestry genome-wide association meta-analysis of Parkinson's disease](https://doi.org/10.1038/s41588-023-01584-8).” *Nature Genetics* (2024).
