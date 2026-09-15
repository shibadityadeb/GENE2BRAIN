# Public website and folded brain QC · 2026-09-15

## Scope

This is a display-only research website revision on `feature/web-anatomical-research-release`. It audits and replaces the old block-like AAL3 voxel-face geometry, expands the public site navigation and documentation, adds shareable URLs and a development atlas-ID view, and keeps the ten completed Stage 10 disease results. Parkinson remains the default worked example. No GWAS, gene-prioritization, AHBA expression, permutation, enrichment, validation, or biology source result was recomputed for the redesign.

## Browser and data verification

- Folded source geometry: 81,924 vertices and 163,840 pial triangles from two fsaverage6 hemispheres; 78 cortical parcels projected by nearest-neighbor neuromaps MNI152 registration fusion, 60 visual AAL3 isosurfaces for noncortical/fallback parcels.
- Region integrity: 138 nonempty mesh IDs exactly match the 138 Parkinson records and each of the ten completed disease datasets. The browser uses the same integer `region_id` join for click, hover, and metrics. The cross-disease Pearson value is serialized from the frozen Stage 10 CSV, not estimated in the client.
- Real Parkinson values: the validator compares observed score, random mean, effect size, Z-score, and FDR q directly against the frozen Stage 6 weighted regional CSV without rounding or substituting values.
- Anatomy-only mode removes research coloring. Diagnostic atlas mode uses deterministic ID-only colors and exposes all 138 region names at `/dev/atlas-validation`. Research mode overlays the stored values on the folded anatomy.
- Left, right, and superior screenshots visibly show cortical gyri/sulci, an anatomically plausible hemispheric divide, and no exposed voxel cubes. Subcortical/fallback parcel meshes use smooth isosurfaces; they are independent from the pial cortex.
- URL state is preserved for disease, metric, hemisphere, camera, and selected AAL3 region. The production preview path fallback and root-relative data assets were browser-tested.
- The paired comparison brains load only when requested, avoiding three automatically active WebGL canvases on the main page. Download links and released CSV/JSON remain static.

Checks performed locally: `python -m unittest discover -s tests -q` (53 pass), `python scripts/validate_web_data.py` (pass), `python scripts/validate_brain_region_mapping.py` (pass), `npm test` (7 pass), `npm run build` (pass), and the full desktop/mobile Playwright suite (13 pass, 7 intentional project-specific skips). The left/right/superior anatomical-orientation browser case passed. The production build still reports a large lazy Three.js chunk (~889 KB uncompressed) and the anatomical JSON is ~21 MiB uncompressed; hosts should serve compressed assets. No formal static-host deployment was performed.

Temporary tunnel check: the existing `unruly-footprint-clump.ngrok-free.dev` endpoint returned Vite HTTP 403 (“This host is not allowed”) before Vite was restarted with `GENE2BRAIN_TUNNEL_HOST` set to that exact host. After restart, the local page and public tunnel both returned HTTP 200 `text/html`, while an unrelated Host header still returned HTTP 403. Starting another ngrok process against the occupied endpoint produced `ERR_NGROK_334`; the existing tunnel was left running. This temporary dev link is not a formal deployment or an asset-license clearance.

## Visual evidence

- [Old voxel-face validation image](../results/figures/web_brain_before.png) — copied from the pre-redesign Git commit for an apples-to-apples diagnostic comparison.
- [Folded anatomy only](../results/figures/web_brain_anatomical.png), [left](../results/figures/web_brain_left_anatomical.png), [right](../results/figures/web_brain_right_anatomical.png), [superior](../results/figures/web_brain_superior_anatomical.png).
- [Real Parkinson overlay](../results/figures/web_brain_gene2brain.png), [AAL3 ID validation](../results/figures/web_atlas_region_validation.png).
- Website sections: [home](../results/figures/web/web_home.png), [region detail](../results/figures/web/web_region_detail.png), [methods](../results/figures/web/web_methods.png), [results](../results/figures/web/web_results.png), [biology](../results/figures/web/web_biology.png).

## Cautions before public deployment

The projected cortical atlas border is an approximation at 2 mm AAL3 resolution, and some edges remain visibly angular despite real folded anatomy. fsaverage cortex and MNI AAL3 subcortex have ID correspondence, not guaranteed point-for-point spatial registration. Visualization geometry is not the analytical atlas mask. The site explicitly warns that genetic-expression enrichment is not a pathology, disease-origin, causal, diagnostic, or patient-prediction map.

Researcher identity, affiliation, contact, and publication status remain unconfirmed placeholders. The neuromaps transformation package identifies CC BY-NC-SA 4.0 terms; Nilearn lists the fsaverage6 dataset license as unknown. Review surface/atlas redistribution and derivative-asset terms before deploying. See [geometry method](web_brain_geometry_method.md) and [deployment guide](../web/DEPLOYMENT.md).
