import type { AtlasGeometry, EnrichmentData, ProjectMetadata, RegionRecord } from './types'

const DATA_ROOT = './data'

async function fetchJson<T>(name: string): Promise<T> {
  const response = await fetch(`${DATA_ROOT}/${name}`)
  if (!response.ok) throw new Error(`Could not load ${name} (${response.status})`)
  return response.json() as Promise<T>
}

export async function loadResearchData() {
  const [enrichment, metadata, geometry] = await Promise.all([
    fetchJson<EnrichmentData>('parkinson_brain_enrichment.json'),
    fetchJson<ProjectMetadata>('project_metadata.json'),
    fetchJson<AtlasGeometry>('aal3_regions.json'),
  ])
  validateClientData(enrichment, metadata, geometry)
  return { enrichment, metadata, geometry }
}

export function validateClientData(
  enrichment: EnrichmentData,
  metadata: ProjectMetadata,
  geometry: AtlasGeometry,
) {
  const recordIds = enrichment.regions.map((region) => region.region_id)
  const geometryIds = geometry.regions.map((region) => region.region_id)
  if (new Set(recordIds).size !== recordIds.length) throw new Error('Duplicate research region IDs')
  if (new Set(geometryIds).size !== geometryIds.length) throw new Error('Duplicate atlas region IDs')
  if (recordIds.length !== geometryIds.length || recordIds.some((id) => !geometryIds.includes(id))) {
    throw new Error('Research records and atlas geometry do not map one-to-one')
  }
  const numeric: Array<keyof RegionRecord> = [
    'observed_score', 'random_mean', 'random_std', 'z_score', 'empirical_p', 'fdr_p', 'effect_size',
    'spatial_null_p', 'spatial_null_fdr', 'spatial_robustness', 'robustness_rank',
  ]
  enrichment.regions.forEach((region) => {
    numeric.forEach((field) => {
      if (!Number.isFinite(region[field] as number)) throw new Error(`Invalid ${field} for region ${region.region_id}`)
    })
    if ([region.fdr_p, region.empirical_p, region.spatial_null_p, region.spatial_null_fdr, region.spatial_robustness].some((value) => value < 0 || value > 1)) {
      throw new Error(`Invalid probability for region ${region.region_id}`)
    }
    if (region.validation_score !== null && !Number.isFinite(region.validation_score)) {
      throw new Error(`Invalid validation_score for region ${region.region_id}`)
    }
    if (!['high_prediction_high_validation', 'high_prediction_low_validation', 'low_prediction_high_validation', 'low_prediction_low_validation', 'not_measured'].includes(region.agreement_status)) {
      throw new Error(`Invalid agreement_status for region ${region.region_id}`)
    }
  })
  if (!(metadata.analysis.fdr_threshold > 0 && metadata.analysis.fdr_threshold < 1)) {
    throw new Error('Invalid configured FDR threshold')
  }
}

export const diseaseRegistry = [{ id: 'parkinson-disease', label: 'Parkinson disease', available: true }]
