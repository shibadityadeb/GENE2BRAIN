import { describe, expect, it } from 'vitest'
import { validateClientData } from './data'
import type { AtlasGeometry, EnrichmentData, ProjectMetadata, RegionRecord } from './types'

const region: RegionRecord = {
  region_id: 77, region_name: 'Putamen L', atlas_id: 'AAL3v1:77', atlas_label: 'Putamen_L',
  hemisphere: 'L', broad_system: 'subcortex', centroid_mni: [-24, 4, 2], observed_score: 0.5,
  random_mean: 0.49, random_std: 0.01, z_score: 1, empirical_p: 0.1, fdr_p: 0.2,
  effect_size: 0.01, number_of_genes: 123,
  spatial_null_p: 0.2, spatial_null_fdr: 0.8, spatial_robustness: 0.8,
  spatial_robustness_label: 'not_robust', robustness_rank: 4, spatial_isolate: false,
  validation_score: 0.136, validation_region: 'Lput', validation_mapping_confidence: 'high',
  agreement_status: 'high_prediction_high_validation',
  biology: null,
}
const enrichment = { regions: [region] } as EnrichmentData
const geometry = { regions: [{ region_id: 77, positions: [0, 0, 0], indices: [0], voxel_count: 1 }] } as AtlasGeometry
const metadata = { analysis: { fdr_threshold: 0.05 } } as ProjectMetadata

describe('client data validation', () => {
  it('accepts a one-to-one finite record mapping', () => {
    expect(() => validateClientData(enrichment, metadata, geometry)).not.toThrow()
  })
  it('rejects missing atlas geometry', () => {
    expect(() => validateClientData(enrichment, metadata, { ...geometry, regions: [] })).toThrow(/one-to-one/)
  })
  it('rejects impossible probabilities and NaN', () => {
    expect(() => validateClientData({ ...enrichment, regions: [{ ...region, fdr_p: 1.2 }] }, metadata, geometry)).toThrow(/probability/)
    expect(() => validateClientData({ ...enrichment, regions: [{ ...region, z_score: Number.NaN }] }, metadata, geometry)).toThrow(/Invalid z_score/)
  })
})
