# Brain visualization problem audit

## Current approach

The web application uses React, TypeScript, React Three Fiber, and Three.js. `scripts/prepare_web_data.py` reads the project’s AAL3v1 NIfTI, a 2 mm MNI-space label volume, and creates one mesh for each of the 138 retained research parcels by enumerating every exposed voxel face. Adjacent coplanar faces are merged into rectangles. `web/public/data/aal3_regions.json` contains positions, triangle indices, and AAL3 region IDs. `BrainScene.tsx` combines these into one GPU mesh and resolves hover/click from the triangle’s region ID. Scientific values come from generated Stage 6–10 JSON, not from geometry.

## Problem

The brain is visually an external surface of labeled 2 mm voxels. Greedy face merging reduces file size, but does not add anatomical curvature, gyri, or sulci. Cortical and subcortical parcels are extracted identically. A lighting change or normal recomputation cannot remove the rectilinear steps. The current atlas-validation screenshot visibly confirms the block artifacts.

## Root cause

The rendered shape is derived from the low-resolution *parcellation volume*, not a high-resolution anatomical brain surface. Each visible polygon is axis-aligned to the atlas grid. Region colors dominate because parcel surfaces are the entire visible anatomy. This is `atlas label → voxel exterior → colored rectangle`, not `anatomical surface → atlas label → research overlay`.

## Proposed replacement

Use a continuous 1 mm MNI152 tissue-derived cortical surface for recognizable anatomy, with AAL3v1 labels sampled only to assign the research overlay and hit-test region IDs. Preserve the original AAL3v1 integer IDs as the sole statistical join key. Generate smooth subcortical surfaces from the same AAL3 label volume with marching-cubes extraction and explicitly validate each mesh’s region ID, triangle indices, and research-data join. Keep anatomy and numerical research records as separate deployable assets. Document mismatches, interpolation, source resolution, and any structures not represented rather than inventing them.

The first replacement will be accepted only after anatomy-only visual inspection, region-ID parity validation, and browser hover/click checks for known cortical and subcortical parcels. No GWAS, AHBA, L2G, null-model, FDR, spatial, or validation statistics will be recalculated.
