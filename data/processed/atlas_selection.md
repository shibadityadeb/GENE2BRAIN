# Atlas selection: AAL3v1 (AAL3v2 distribution)

## Selection

- **Atlas:** Automated Anatomical Labeling atlas 3, AAL3v1 image from the April 2024 AAL3v2 SPM12 distribution
- **Regions:** 166 nonzero labelled regions
- **Source:** GIN / CNRS AAL project; Rolls et al. (2020), *NeuroImage* 206:116189
- **Official archive:** https://www.gin.cnrs.fr/wp-content/uploads/AAL3v2_for_SPM12.tar.gz
- **Article:** https://doi.org/10.1016/j.neuroimage.2019.116189
- **Resolution used:** 2 mm isotropic
- **Image dimensions:** (91, 109, 91)
- **Space:** MNI
- **Atlas image:** `data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz`
- **Label definition:** `data/atlases/aal_3v2/AAL3/AAL3v1.xml`
- **AHBA compatibility:** Deterministic integer-labelled volumetric NIfTI in MNI space, validated with `abagen.images.check_atlas`; hemisphere and tissue-class metadata are supplied to constrain sample matching.

## Why AAL3 was selected

AAL3 provides named bilateral cortical regions, basal ganglia, hippocampus,
amygdala, 15 thalamic subdivisions per side, cerebellar lobules/vermis, and
small disease-relevant nuclei including substantia nigra, nucleus accumbens,
VTA, red nucleus, locus coeruleus, and raphe. This coverage is materially
better aligned with a cross-neurological/psychiatric project than a
cortex-dominant parcellation.

The bundled 83-region Desikan–Killiany option was considered because it is
directly distributed by `abagen` and has robust cortical/subcortical regions,
but it lacks explicit cerebellar parcels and substantia nigra. AAL3 was not
selected merely for convenience; its anatomical scope is the deciding factor.

## Important limitation

AAL3's smallest brainstem and thalamic nuclei may be undersampled by the AHBA.
Atlas presence therefore does not imply usable transcriptomic coverage. No
missing parcel is imputed. Regions enter the main reference only when at least
2 samples from at least 2 donors are assigned;
the full result is recorded in `key_region_coverage.csv`.

Atlas-level requested-structure presence: substantia nigra: yes, caudate: yes, putamen: yes, globus pallidus: yes, nucleus accumbens: yes, hippocampus: yes, amygdala: yes, thalamus: yes, cerebellum: yes, brainstem: yes, frontal cortex: yes, temporal cortex: yes, parietal cortex: yes, occipital cortex: yes.
