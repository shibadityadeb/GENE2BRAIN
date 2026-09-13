# Stage 02 methods

## Inputs and software

- AHBA source: Allen Human Brain Atlas normalized microarray dataset, downloaded in Stage 1 with `abagen.fetch_microarray(donors='all')`.
- Atlas: AAL3v1 from the April 2024 AAL3v2 distribution, 2 mm MNI NIfTI, 166 labelled parcels.
- Python: 3.12.7
- abagen: 0.1.3
- Current inspected API: `abagen.get_expression_data(atlas, atlas_info=None, *, ibf_threshold=0.5, probe_selection='diff_stability', donor_probes='aggregate', sim_threshold=None, lr_mirror=None, exact=None, missing=None, tolerance=2, sample_norm='srs', gene_norm='srs', norm_matched=True, norm_structures=False, region_agg='donors', agg_metric='mean', corrected_mni=True, reannotated=True, return_counts=False, return_donors=False, return_report=False, donors='all', data_dir=None, verbose=0, n_proc=1)`
- Random seed: 42 (PCA is deterministic here; recorded for reproducibility).

## Exact preprocessing parameters

```python
{'ibf_threshold': 0.5, 'probe_selection': 'diff_stability', 'donor_probes': 'aggregate', 'sim_threshold': None, 'lr_mirror': None, 'missing': None, 'tolerance': 2, 'sample_norm': 'srs', 'gene_norm': 'srs', 'norm_matched': True, 'norm_structures': True, 'region_agg': 'donors', 'agg_metric': 'mean', 'corrected_mni': True, 'reannotated': True}
```

## Decisions, rationale, and alternatives

1. **Gene reannotation (`reannotated=True`).** Arnatkevičiūtė et al. mappings distributed with `abagen` replace outdated/ambiguous Allen probe annotations and discard probes without reliable mappings. Alternative: raw Allen annotations, which preserve more probes but increase mapping error.
2. **Poor-signal probes (`ibf_threshold=0.5`).** A probe must exceed background in at least 50% of samples across donors. This is the documented `abagen` default and avoids expression dominated by noise. Alternatives include more permissive/strict thresholds or no filtering.
3. **Multiple probes (`probe_selection='diff_stability'`, `donor_probes='aggregate'`).** One probe per gene is selected for the most reproducible regional pattern across donor pairs. Aggregate selection guarantees the same gene/probe basis across donors. Alternatives include maximum intensity/variance, RNA-seq concordance, or averaging probes.
4. **Sample quality (`sim_threshold=None`).** No inter-areal-similarity outlier filter is applied because a defensible threshold was not prespecified; discarded samples would otherwise be difficult to interpret. Alternative: a preregistered similarity threshold.
5. **Coordinates (`corrected_mni=True`).** Corrected coordinates distributed with `abagen` are used. Alternative: original Allen MNI coordinates.
6. **Sample assignment (`tolerance=2`, `exact=None`).** Samples inside a parcel or within 2 mm are assigned to the closest parcel centroid, constrained by atlas hemisphere and broad tissue class. Unmatched samples remain unassigned. Alternatives include exact-only assignment or a larger tolerance, the latter increasing questionable assignments.
7. **Hemisphere (`lr_mirror=None`).** Samples are not mirrored. Four donors have only left-hemisphere sampling, so right-hemisphere evidence is limited to the two bilaterally sampled donors. Mirroring can improve coverage but creates synthetic observations and was rejected for the reference matrix.
8. **Normalization (`sample_norm='srs'`, `gene_norm='srs'`).** Scaled robust sigmoid normalization is performed first across genes within each sample, then for each gene across matched samples separately within donor. `norm_matched=True` avoids unmatched samples affecting scaling. `norm_structures=True` performs gene normalization separately within cortex, subcortex/brainstem, and cerebellum so gross tissue-class differences do not dominate. Alternative global normalization (`norm_structures=False`) preserves class-level offsets but can overwhelm regional effects.
9. **Regional and donor aggregation (`region_agg='donors'`, `agg_metric='mean'`).** Samples are averaged within region separately for each donor. The main matrix is then the equal-weight mean across donors with data, preventing donors with more tissue samples from dominating. Median aggregation is a robust alternative.
10. **Missing data (`missing=None`).** No nearest-centroid or interpolation imputation is used. Donor-region absences remain `NaN`. The main matrix retains parcels with at least 2 assigned samples across at least 2 donors and averages available donors; this removes 28 of 166 parcels. Alternatives include retaining single-donor regions or spatial imputation, both rejected to protect robustness.
11. **Gene filtering after abagen.** `abagen` returned 15,633 genes after reannotation, intensity filtering, and probe selection. We removed 1 genes with at least one non-finite pooled regional value rather than imputing them, leaving 15,632. The final matrix includes 0 genes at or below the near-zero variance threshold (1e-12); these are reported, not silently deleted, and PCA excludes them only to avoid zero-scale features.

## Gene and region metadata

`gene_metadata.csv` records Entrez IDs and reannotated AHBA candidate probe IDs/names for each output symbol. Because the public `abagen` return object exposes genes rather than selected probe IDs, the candidate list is reported transparently rather than inferring an exact selected ID. `region_metadata.csv` retains every AAL3 parcel and marks its sample/donor coverage and inclusion status. `regional_expression_summary.csv` records the across-gene mean, median, standard deviation, minimum, and maximum for every retained region.

## Generated abagen report

Regional microarry expression data were obtained from 6 post-mortem brains (1 female, ages 24.0--57.0, 42.50 +/- 13.38) provided by the Allen Human Brain Atlas (AHBA, https://human.brain-map.org; [H2012N]). Data were processed with the abagen toolbox (version 0.1.3; https://github.com/rmarkello/abagen) using a 166-region volumetric atlas in MNI space.

First, microarray probes were reannotated using data provided by [A2019N]; probes not matched to a valid Entrez ID were discarded. Next, probes were filtered based on their expression intensity relative to background noise [Q2002N], such that probes with intensity less than the background in >=50.00% of samples across donors were discarded , yielding 31,569 probes . When multiple probes indexed the expression of the same gene, we selected and used the probe with the most consistent pattern of regional variation across donors (i.e., differential stability; [H2015N]), calculated with:

$$ \Delta_{{S}}(p) = \frac{{1}}{{\binom{{N}}{{2}}}} \, \sum_{{i=1}}^{{N-1}} \sum_{{j=i+1}}^{{N}} \rho[B_{{i}}(p), B_{{j}}(p)] $$

where $ \rho $ is Spearman's rank correlation of the expression of a single probe, p, across regions in two donors $B_{{i}}$ and $B_{{j}}$, and N is the total number of donors. Here, regions correspond to the structural designations provided in the ontology from the AHBA.

The MNI coordinates of tissue samples were updated to those generated via non-linear registration using the Advanced Normalization Tools (ANTs; https://github.com/chrisfilo/alleninf). Samples were assigned to brain regions in the provided atlas if their MNI coordinates were within 2 mm of a given parcel. To reduce the potential for misassignment, sample-to-region matching was constrained by hemisphere and gross structural divisions (i.e., cortex, subcortex/brainstem, and cerebellum, such that e.g., a sample in the left cortex could only be assigned to an atlas parcel in the left cortex; [A2019N]). All tissue samples not assigned to a brain region in the provided atlas were discarded.

Inter-subject variation was addressed by normalizing tissue sample expression values across genes using a robust sigmoid function [F2013J]:

$$ x_{{norm}} = \frac{{1}}{{1 + \exp(-\frac{{(x-\langle x \rangle)}} {{\text{{IQR}}_{{x}}}})}} $$

where $\langle x \rangle$ is the median and $\text{{IQR}}_{{x}}$ is the normalized interquartile range of the expression of a single tissue sample across genes. Normalized expression values were then rescaled to the unit interval:

$$ x_{{scaled}} = \frac{{x_{{norm}} - \min(x_{{norm}})}} {{\max(x_{{norm}}) - \min(x_{{norm}})}} $$

Gene expression values were then normalized across tissue samples using an identical procedure. Normalization was performed separately for samples in distinct structural classes (i.e., cortex, subcortex/brainstem, cerebellum). Samples assigned to the same brain region were averaged separately for each donor, yielding a regional expression matrix for each donor with 166 rows, corresponding to brain regions, and 15,633 columns, corresponding to the retained genes .

REFERENCES
----------
[A2019N]: Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). A practical guide to linking brain-wide gene expression and neuroimaging data. Neuroimage, 189, 353-367.
[F2013J]: Fulcher, B. D., Little, M. A., & Jones, N. S. (2013). Highly comparative time-series analysis: the empirical structure of time series and their methods. Journal of the Royal Society Interface, 10(83), 20130048.
[H2012N]: Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng, L., Miller, J. A., ... & Jones, A. R. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416), 391-399.
[H2015N]: Hawrylycz, M., Miller, J. A., Menon, V., Feng, D., Dolbeare, T., Guillozet-Bongaarts, A. L., ... & Lein, E. (2015). Canonical genetic signatures of the adult human brain. Nature Neuroscience, 18(12), 1832.
[Q2002N]: Quackenbush, J. (2002). Microarray data normalization and transformation. Nature Genetics, 32(4), 496-501.


No downstream disease analysis is part of this stage.
