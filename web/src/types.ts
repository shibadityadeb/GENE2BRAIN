export type Metric = 'z_score' | 'observed_score' | 'fdr_p' | 'spatial_robustness' | 'validation_score' | 'agreement'
export type Hemisphere = 'whole' | 'L' | 'R'
export type ViewPreset = 'reset' | 'anterior' | 'posterior' | 'superior' | 'inferior'

export interface RegionRecord {
  region_id: number
  region_name: string
  atlas_id: string
  atlas_label: string
  hemisphere: 'L' | 'R' | 'B'
  broad_system: string
  centroid_mni: [number, number, number]
  observed_score: number
  random_mean: number
  random_std: number
  z_score: number
  empirical_p: number
  fdr_p: number
  effect_size: number
  number_of_genes: number
  spatial_null_p: number
  spatial_null_fdr: number
  spatial_robustness: number
  spatial_robustness_label: 'robust' | 'not_robust'
  robustness_rank: number
  spatial_isolate: boolean
  validation_score: number | null
  validation_region: string | null
  validation_mapping_confidence: 'high' | 'medium' | null
  agreement_status: 'high_prediction_high_validation' | 'high_prediction_low_validation' | 'low_prediction_high_validation' | 'low_prediction_low_validation' | 'not_measured'
}

export interface EnrichmentData {
  schema_version: string
  disease_id: string
  disease_name: string
  gene_set_version: string
  source_file: string
  generated_on: string
  regions: RegionRecord[]
}

export interface GeometryRegion {
  region_id: number
  positions: number[]
  indices: number[]
  voxel_count: number
}

export interface AtlasGeometry {
  schema_version: string
  atlas: string
  coordinate_system: string
  region_count: number
  regions: GeometryRegion[]
}

export interface ProjectMetadata {
  project: string
  title: string
  generated_on: string
  counts: Record<string, number>
  analysis: {
    fdr_threshold: number
    significant_regions: number
    primary_gene_set: string
    gwas_accession: string
    gwas_title: string
    gene_prioritization: string
    null_model: string
    null_draws_available: boolean
    null_summary_available: boolean
    null_summary_note: string
    spatial_sensitivity: {
      method: string
      neighbor_definition: string
      permutations: number
      global_statistic: string
      global_p: number
      robust_rule: string
      robust_regions: number
      isolated_regions: number
    }
    independent_validation: {
      dataset: string
      phenotype: string
      matched_validation_units: number
      mapped_aal3_regions: number
      pearson_r: number
      pearson_p: number
      pearson_ci: [number, number]
      spearman_rho: number
      spearman_p: number
      spatial_null_p: number
      interpretation: 'SUPPORTED' | 'PARTIALLY SUPPORTED' | 'NOT SUPPORTED'
      independence_note: string
      substantia_nigra: string
    }
  }
  metrics: Record<string, { label: string; domain?: number[]; threshold?: number; scale: string }>
  atlas: Record<string, unknown>
  sources: Array<{ name: string; url: string }>
}
