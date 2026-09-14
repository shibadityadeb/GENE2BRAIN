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

The development server runs at `http://localhost:8000`.

The production bundle is written to `web/dist/`. `npm run preview` serves that
bundle locally. Atlas validation mode is available at `?mode=atlas`.

Before a build, the repository-level data checks can be run with:

```bash
python scripts/validate_web_data.py
```

Regenerate checked-in web artifacts after an authoritative pipeline output changes:

```bash
python scripts/prepare_web_data.py
python scripts/build_multidisease_web_data.py
```

## Research data

`data/web/multidisease_atlas.json` is serialized from the Stage 10 regional,
spatial, validation, and biological result tables and joined to
`data/processed/region_metadata.csv`. Each completed disease contains the exact
stored observed score, null mean and standard deviation, matched-null Z-score,
empirical p-value, FDR q-value, effect size, gene-set size, and the separately
identified spatial-null values and joint robustness call. A single global Z domain
is used for cross-disease visual comparison. Canonical and deployed copies are
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

The mesh is not a generic brain. `scripts/prepare_web_data.py` reads the exact 2 mm
MNI-space image used in Stage 2:

`data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz`

For each of the 138 analyzed regions, `region_id` equals the integer AAL3 voxel
label. Exposed voxel faces are extracted and neighboring coplanar faces are merged
into rectangles. This reduces triangles without moving or smoothing boundaries.
The resulting `aal3_regions.json` maintains a separate indexed geometry for each
region ID. At load time the app combines these into one GPU mesh and retains a
triangle-to-region lookup for hover/click selection. MNI `(x, y, z)` is displayed as
WebGL `(x, z, y)`, keeping superior up. The validation script requires exact equality
between geometry and research record IDs.

## Adding a future disease

To add a disease after it passes the frozen source, genetic-evidence, and AHBA
coverage gates, regenerate the Stage 10 tables and run
`scripts/build_multidisease_web_data.py`. Never add placeholders or copy values
from another disease.

## Tests

- `scripts/validate_web_data.py` checks ID parity, duplicates, required fields,
  finite values, probabilities, geometry indices, and threshold consistency.
- Vitest checks client-side rejection of invalid/missing data and color semantics.
- Playwright exercises desktop and mobile production builds, selectors, search,
  region details, null-summary disclosure, reset, downloads, and atlas validation.
