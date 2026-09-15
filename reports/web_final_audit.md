# GENE2BRAIN public website final audit · 2026-09-15

## Architecture

The website is a Vite/React/TypeScript static single-page research atlas in `web/`. The hero explorer has a lazily imported React Three Fiber/Three.js brain, a disease/metric/hemisphere selector, region search, camera presets, hover tooltip, click detail, an observed-versus-retained-null summary, and research caveats. Comparison brains are opt-in to avoid running three WebGL canvases by default. Research, before/after geometry, Methods, Data, Results, Biology, About/Researcher, and Citations are addressable sections. Query parameters preserve disease, metric, hemisphere, view, and region; `/dev/atlas-validation` shows diagnostic region-ID colors. `web/public/data/` contains read-only JSON/CSV and a digest-stamped `release.json`. No score, permutation, or statistical significance is produced by the browser. The displayed cross-disease Pearson coefficients are copied from the frozen Stage 10 similarity table.

## Brain geometry source and atlas mapping

`scripts/build_anatomical_brain.py` downloads the official Nilearn FreeSurfer **fsaverage6 pial** surfaces and uses the neuromaps 0.0.7 41k MNI152-to-fsaverage registration-fusion coordinates with nearest-neighbor sampling of the same `AAL3v1.nii.gz` 2 mm integer labels used in Stage 2. Pial triangles with majority-agreeing retained cortical labels form 78 region meshes; 60 noncortical/fallback visual meshes use AAL3 label isosurfaces. The active web asset has 81,924 pial vertices, 163,840 anatomy triangles, and 138 nonempty regional meshes. Research CSV, atlas image, regional mesh, and web record are keyed to the **same original integer AAL3 region ID**. The former exposed-voxel-face JSON was removed from the public bundle, while the before screenshot and Git history retain it for audit. The surface, mapping, and coordinate caveats are detailed in [web brain geometry method](web_brain_geometry_method.md).

## Data source and current disease

The healthy expression reference is the Allen Human Brain Atlas microarray resource processed by abagen into the frozen AAL3 regional gene matrix. GWAS and Open Targets L2G prioritization, matched-gene permutation scores, spatial sensitivity, independent Parkinson ENIGMA structural-MRI phenotype comparison, and downstream biology are read from committed Stage 3–10 result tables. The initial hero is Parkinson disease, the documented reference analysis. The existing ten completed eligible diseases remain available; five screened diseases excluded or under review are explicitly listed rather than given placeholder maps.

## Current metrics

Anatomy only (no research colors), matched-null Z-score, observed expression, FDR significance, spatial robustness, independent validation (nullable), and prediction/validation agreement (nullable). Diagnostic atlas validation colors are ID colors only. Missing validation is never imputed. Significance and robustness calls retain the frozen source FDR rules.

## Validation status

Local Python project suite: 53 passed. Web unit suite: 7 passed. `scripts/validate_web_data.py` and `scripts/validate_brain_region_mapping.py` pass, including canonical/deployed data identity, release SHA-256 digests, 138-region/ten-disease one-to-one mapping, geometry indices, finite fields, and real weighted Parkinson Stage 6 value equality. Production `npm run build` passes. The browser suite passed 14 tests with 8 intentional platform-specific skips after the long exploration test was split and the redundant atlas screenshot test was removed, including left/right/superior anatomy validation. Browser checks cover WebGL hover/click, search/detail, disease/metric controls, downloads, responsive mobile gestures, deep-link preservation, and developer atlas route. See [visual/release QC](web_research_release_qc.md) and its screenshot evidence.

## Deployment status

The production static bundle is locally buildable and previewable but has **not been formally deployed**. A temporary existing ngrok dev tunnel was tested after Vite's exact-host allowlist was configured: the public HTML returned HTTP 200 and unrelated Host headers remained blocked. See [deployment instructions](../web/DEPLOYMENT.md). The main anatomical JSON is approximately 21 MiB uncompressed; serve gzip/Brotli and cache immutable assets. The lazy Three.js chunk is approximately 889 KB uncompressed. Review formal asset redistribution/attribution before a public host is launched: the neuromaps distribution identifies CC BY-NC-SA 4.0 terms, and Nilearn's fsaverage6 dataset page lists the surface data license as unknown.

## Known limitations and manual information

The nearest-neighbor AAL3 atlas border can remain angular at 2 mm resolution, even though the cortex is continuously folded. Registration fusion maps labels to fsaverage vertices but does not make fsaverage pial and MNI subcortical geometry perfectly co-registered; displayed parcel outlines are approximate. The six AHBA donors are unevenly sampled. Parkinson external structural-MRI validation has an overall **NOT SUPPORTED** interpretation and is not a pathology map. None of the visual patterns establish causality, origin, diagnosis, or patient-level vulnerability.

Researcher **name, affiliation, contact**, and **publication status** remain pending confirmation in the website and project metadata. No names, portraits, institutional relationships, or endorsement claims have been fabricated. Those details and upstream redistribution review require manual resolution before formal public release.
