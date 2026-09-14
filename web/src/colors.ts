import { Color } from 'three'
import type { Metric, RegionRecord } from './types'

const NEGATIVE = new Color('#2166ac')
const NEUTRAL = new Color('#f2efe5')
const POSITIVE = new Color('#b2182b')
const OBSERVED_LOW = new Color('#d9eee9')
const OBSERVED_HIGH = new Color('#006d77')
const NOT_SIGNIFICANT = new Color('#526068')
const SIGNIFICANT = new Color('#ffb000')
const SPATIAL_LOW = new Color('#e8e7ee')
const SPATIAL_HIGH = new Color('#54278f')

export function colorForRegion(
  region: RegionRecord,
  metric: Metric,
  metadata: { zDomain: [number, number]; observedDomain: [number, number]; fdrThreshold: number },
  validationMode = false,
) {
  if (validationMode) {
    return new Color().setHSL((region.region_id * 0.61803398875) % 1, 0.48, 0.61)
  }
  if (metric === 'fdr_p') return (region.fdr_p < metadata.fdrThreshold ? SIGNIFICANT : NOT_SIGNIFICANT).clone()
  if (metric === 'spatial_robustness') return SPATIAL_LOW.clone().lerp(SPATIAL_HIGH, region.spatial_robustness)
  if (metric === 'observed_score') {
    const [low, high] = metadata.observedDomain
    const t = Math.max(0, Math.min(1, (region.observed_score - low) / (high - low)))
    return OBSERVED_LOW.clone().lerp(OBSERVED_HIGH, t)
  }
  const bound = Math.max(Math.abs(metadata.zDomain[0]), Math.abs(metadata.zDomain[1]))
  const t = Math.min(1, Math.abs(region.z_score) / bound)
  return NEUTRAL.clone().lerp(region.z_score < 0 ? NEGATIVE : POSITIVE, t)
}
