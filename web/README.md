# GENE2BRAIN interactive brain

This React + TypeScript application visualizes the validated Stage 10 atlas for
the 10 diseases that passed the frozen eligibility gates. It provides a primary
disease explorer, shared-scale side-by-side brains, a region-first cross-disease
ranking, Stage 7 spatial sensitivity, available external validation, and biological
interpretation on the same AAL3v1 parcels used by the scientific pipeline. Missing
validation stays missing. The browser displays authoritative values; it does not
calculate scores, permutations, p-values, FDR values, or Z-scores.

## Run locally

From `web/`:

```bash
npm install
npm run dev
npm test
npm run build
```

The development server runs at `http://127.0.0.1:8000` (also reachable via
`localhost` on most machines).

To share a temporary development link through ngrok, start `ngrok http 8000`,
copy the hostname from its HTTPS URL, and restart Vite with that exact hostname:

```bash
GENE2BRAIN_TUNNEL_HOST=your-assigned-host.ngrok-free.dev npm run dev
```

For example, set the variable to `unruly-footprint-clump.ngrok-free.dev` if
your URL is `https://unruly-footprint-clump.ngrok-free.dev/`. Vite otherwise
returns HTTP 403 with “This host is not allowed.” Do not use
`server.allowedHosts: true`: it would admit arbitrary Host headers. If ngrok
reports `ERR_NGROK_334`, the endpoint is already online in another session;
use that existing session or stop it before opening another. A free ngrok
browser-warning page is separate from Vite's host rejection.

The production bundle is written to `web/dist/`. `npm run preview` serves that
bundle locally. Atlas validation mode is available at `/dev/atlas-validation`
or the older `?mode=atlas`. The disease, metric, hemisphere, region and camera
view are encoded in the URL; use **Share this view** to copy a reproducible link.

Before a build, the repository-level data checks can be run with:

```bash
python scripts/validate_web_data.py
```

Regenerate checked-in web artifacts after an authoritative pipeline output changes:

```bash
python scripts/prepare_web_data.py
python scripts/build_multidisease_web_data.py
python scripts/build_anatomical_brain.py
python scripts/validate_brain_region_mapping.py
```

## Research data

`data/web/multidisease_atlas.json` is serialized from the Stage 10 regional,
spatial, validation, and biological result tables and joined to
`data/processed/region_metadata.csv`. Each completed disease contains the exact
stored observed score, null mean and standard deviation, matched-null Z-score,
empirical p-value, FDR q-value, effect size, gene-set size, and the separately
identified spatial-null values and joint robustness call. A single global Z domain
is used for cross-disease visual comparison. Paired Pearson correlations are
serialized from the frozen Stage 10 similarity table, not recomputed in the
browser. Canonical and deployed copies are
required to be byte-identical.

Stage 8 adds nullable ENIGMA-PD validation scores, the source ENIGMA parcel,
crosswalk confidence, and the pre-specified median-split agreement status. The
scientific correlation uses 70 unique ENIGMA parcels; visualization-only propagation
to AAL3 never increases the inferential sample size.

Stage 9–10 add biological-interpretation payloads. Only parcels selected in advance
by the frozen Stage 7 robustness-rank rule receive detailed regional drivers. GO,
Reactome, WikiPathways and Human Protein Atlas human-brain cell-type results are
read from committed analysis outputs; the browser does not perform enrichment.
The “Why is this region highlighted?” disclosure separates observed evidence
from post-hoc biological interpretation and reports corrected negative results.

Stage 6 retained per-region null mean and standard deviation but not the 10,000
individual scores. The detail panel therefore shows a clearly labeled observed vs.
null mean ± 1 SD summary. It does not synthesize a histogram or assume normality.

Stage 7 uses a binary 26-neighbor AAL3 voxel-contact graph and 10,000 singleton
Moran spectral randomizations. The frontend reads those pipeline outputs; it never
generates surrogates or changes the pre-specified joint robustness rule.

## Atlas and region-ID mapping

The active mesh is a real folded fsaverage6 pial cortical surface, not a generic
brain image or a blocky region-colored mask. `scripts/build_anatomical_brain.py`
projects labels from the exact 2 mm AAL3 MNI-space image used in Stage 2:

`data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz`

For each of the 138 analyzed regions, `region_id` equals the integer AAL3 voxel
label. Official neuromaps registration-fusion coordinates map those labels to
fsaverage6 vertices with nearest-neighbor sampling; a majority rule assigns
retained cortical IDs to pial triangles. Other regions use separate interpolated
AAL3 label isosurfaces. The original integer atlas and research tables are never
resampled for statistical analysis. The browser combines indexed overlay meshes
and retains a triangle-to-region lookup for hover/click selection. RAS `(x, y, z)`
is displayed as WebGL `(x, z, y)`, keeping superior up. fsaverage cortex and MNI
subcortex are an approximate visual combination, not perfectly registered
boundaries. See [geometry method](../reports/web_brain_geometry_method.md).

## Adding a future disease

To add a disease after it passes the frozen source, genetic-evidence, and AHBA
coverage gates, regenerate the Stage 10 tables and run
`scripts/build_multidisease_web_data.py`. Never add placeholders or copy values
from another disease.

## Tests

- `scripts/validate_web_data.py` and `scripts/validate_brain_region_mapping.py`
  check ID parity, duplicates, required fields,
  finite values, probabilities, geometry indices, and threshold consistency.
- Vitest checks client-side rejection of invalid/missing data and color semantics.
- Playwright exercises desktop and mobile production builds, selectors, search,
  region details, null-summary disclosure, reset, downloads, and atlas validation.
