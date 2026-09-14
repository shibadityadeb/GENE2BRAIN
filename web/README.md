# GENE2BRAIN interactive brain

This React + TypeScript application visualizes the validated Stage 6 Parkinson
regional-enrichment results, Stage 7 spatial sensitivity, and Stage 8 external
phenotype and agreement layers on the same AAL3v1 parcels used by the scientific
pipeline. The Stage 8 layers are explicitly separate from discovery and unmeasured
parcels remain missing. The browser displays authoritative values; it does not calculate scores,
permutations, p-values, FDR values, or Z-scores.

## Run locally

From `web/`:

```bash
npm install
npm run dev
npm test
npm run build
```

The production bundle is written to `web/dist/`. `npm run preview` serves that
bundle locally. Atlas validation mode is available at `?mode=atlas`.

Before a build, the repository-level data checks can be run with:

```bash
python scripts/validate_web_data.py
```

Regenerate checked-in web artifacts after an authoritative pipeline output changes:

```bash
python scripts/prepare_web_data.py
```

## Research data

`data/web/parkinson_brain_enrichment.json` is serialized from
`data/results/parkinson_regional_enrichment.csv`, joined to
`data/processed/region_metadata.csv`. Its `regions` array contains the exact stored
Stage 6 observed score, null mean and standard deviation, Z-score, empirical p-value,
FDR q-value, effect size, and gene-set size, plus the separately identified Stage 7
spatial-null p-value, FDR, percentile, rank, and joint robustness call. The deployed copy is under
`web/public/data/`. `project_metadata.json` supplies counts, data sources, scale
domains, and the pipeline-defined FDR threshold (q < 0.05).

Stage 8 adds nullable ENIGMA-PD validation scores, the source ENIGMA parcel,
crosswalk confidence, and the pre-specified median-split agreement status. The
scientific correlation uses 70 unique ENIGMA parcels; visualization-only propagation
to AAL3 never increases the inferential sample size.

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

Only Parkinson disease is registered because it is the only completed analysis.
To add a disease after its pipeline outputs are validated:

1. Create a web JSON file using the same region schema without recalculating values.
2. Validate every record against the atlas geometry.
3. Add the dataset to `diseaseRegistry` in `src/data.ts`.
4. Add disease-specific metric domains and analysis metadata.

Never add placeholders or copy Parkinson values to another disease.

## Tests

- `scripts/validate_web_data.py` checks ID parity, duplicates, required fields,
  finite values, probabilities, geometry indices, and threshold consistency.
- Vitest checks client-side rejection of invalid/missing data and color semantics.
- Playwright exercises desktop and mobile production builds, selectors, search,
  region details, null-summary disclosure, reset, downloads, and atlas validation.
