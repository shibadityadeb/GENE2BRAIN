import { describe, expect, it } from 'vitest'
import { colorForRegion } from './colors'
import type { RegionRecord } from './types'

const base = {
  region_id: 1, region_name: 'Region', atlas_id: 'AAL3v1:1', atlas_label: 'Region_L', hemisphere: 'L',
  broad_system: 'cortex', centroid_mni: [0, 0, 0], observed_score: 0.5, random_mean: 0.5,
  random_std: 0.1, empirical_p: 0.5, effect_size: 0, number_of_genes: 123,
  spatial_null_p: 0.2, spatial_null_fdr: 0.8, spatial_robustness: 0.8,
  spatial_robustness_label: 'not_robust', robustness_rank: 4, spatial_isolate: false,
} as RegionRecord
const config = { zDomain: [-4, 4] as [number, number], observedDomain: [0.4, 0.6] as [number, number], fdrThreshold: 0.05 }

describe('scientific color mapping', () => {
  it('uses opposite diverging colors for negative and positive Z scores', () => {
    const negative = colorForRegion({ ...base, z_score: -4, fdr_p: 0.5 }, 'z_score', config)
    const positive = colorForRegion({ ...base, z_score: 4, fdr_p: 0.5 }, 'z_score', config)
    expect(negative.getHexString()).not.toBe(positive.getHexString())
  })
  it('uses the configured FDR threshold', () => {
    const passing = colorForRegion({ ...base, z_score: 0, fdr_p: 0.049 }, 'fdr_p', config)
    const failing = colorForRegion({ ...base, z_score: 0, fdr_p: 0.05 }, 'fdr_p', config)
    expect(passing.getHexString()).not.toBe(failing.getHexString())
  })
  it('maps the Stage 7 spatial percentile monotonically', () => {
    const low = colorForRegion({ ...base, z_score: 0, fdr_p: 0.5, spatial_robustness: 0 }, 'spatial_robustness', config)
    const high = colorForRegion({ ...base, z_score: 0, fdr_p: 0.5, spatial_robustness: 1 }, 'spatial_robustness', config)
    expect(low.getHexString()).not.toBe(high.getHexString())
  })
})
