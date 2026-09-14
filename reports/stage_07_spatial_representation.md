# Stage 07 spatial representation

## Selected representation

Stage 7 uses a **binary, symmetric volumetric AAL3 adjacency graph** derived from
`data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz`, the exact 2 mm MNI image used in Stage 2
and the web geometry. Two retained parcels are neighbors when any of their voxels
touch by a face, edge, or corner (26-connectivity). This criterion is deterministic
and respects atlas contact; it does not claim that Euclidean proximity is biological
connectivity.

- Regions: 138 retained Stage 2 parcels
- Undirected edges: 512
- Degree: minimum 0, median 7.0, maximum 15
- Connected-component sizes: 134, 1, 1, 1, 1
- Isolated zero-weight parcels (4): SN pr R, LC L, LC R, Raphe M
- Weights: binary (1=voxel contact, 0=no contact), diagonal zero
- Normalization: not row-normalized; symmetry is preserved for Moran eigenvectors

Isolated parcels remain in every matrix and output. Their spatial-weight rows are
zero and `spatial_isolate=True`; no missing edge is invented. Their regional
spectral-null results are consequently less anatomically constrained and must be
interpreted cautiously. The matrix is saved with numeric region IDs as both rows
and columns.

## Why this representation and null

A cortical spin requires spherical cortical coordinates and cannot validly rotate
irregular subcortical, cerebellar, or brainstem parcels. AAL3 is volumetric and the
analysis deliberately retains those structures. Moran spectral randomization (MSR)
instead accepts an explicit symmetric neighbor matrix for irregularly arranged
observations. The singleton procedure independently flips each observed Moran-basis
coefficient's sign, preserving the mean, variance, complete Moran power spectrum,
and global Moran's I exactly while changing hotspot locations. This strong constraint
is intentionally conservative and can leave surrogates correlated with the input.

Primary method reference: Wagner HH, Dray S (2015), *Methods in Ecology and
Evolution* 6:1169–1178, https://doi.org/10.1111/2041-210X.12407. BrainSpace's
independent implementation documents the same singleton procedure:
https://brainspace.readthedocs.io/en/latest/generated/brainspace.null_models.moran.moran_randomization.html.
