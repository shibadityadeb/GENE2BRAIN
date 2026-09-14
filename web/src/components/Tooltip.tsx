import { formatValue } from '../format'
import type { RegionRecord } from '../types'

export function Tooltip({ region, point }: { region: RegionRecord; point: { x: number; y: number } }) {
  return (
    <div className="tooltip" style={{ left: point.x + 14, top: point.y + 14 }} role="status">
      <strong>{region.region_name}</strong>
      <span>Z-score <b>{formatValue(region.z_score)}</b></span>
      <span>FDR <b>{formatValue(region.fdr_p, 4)}</b></span>
      <span>Observed <b>{formatValue(region.observed_score, 4)}</b></span>
      <span>Null mean <b>{formatValue(region.random_mean, 4)}</b></span>
      <span>Effect size <b>{formatValue(region.effect_size, 4)}</b></span>
      <span>Spatial percentile <b>{formatValue(region.spatial_robustness, 4)}</b></span>
      <span>Spatial result <b>{region.spatial_robustness_label === 'robust' ? 'Robust' : 'Not robust'}</b></span>
      <span>Genes <b>{region.number_of_genes}</b></span>
    </div>
  )
}
