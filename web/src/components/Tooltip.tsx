import { formatValue } from '../format'
import type { RegionRecord } from '../types'

export function Tooltip({ region, point }: { region: RegionRecord; point: { x: number; y: number } }) {
  return (
    <div className="tooltip" style={{ left: point.x + 14, top: point.y + 14 }} role="status">
      <strong>{region.region_name}</strong>
      <span>{region.z_score > 0 ? 'Disease-linked genes show a relatively strong signal here.' : 'Disease-linked genes show a relatively lower signal here.'}</span>
      <span>Click to view the evidence behind this signal →</span>
    </div>
  )
}
