export type Metric = 'z_score' | 'observed_score' | 'fdr_p'
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
  }
  metrics: Record<string, { label: string; domain?: number[]; threshold?: number; scale: string }>
  atlas: Record<string, unknown>
  sources: Array<{ name: string; url: string }>
}
