import type { Metric } from '../types'
import { formatValue } from '../format'

export function Legend({ metric, zDomain, observedDomain, threshold, significant, validationMode }: {
  metric: Metric
  zDomain: [number, number]
  observedDomain: [number, number]
  threshold: number
  significant: number
  validationMode: boolean
}) {
  if (validationMode) return (
    <div className="legend" aria-label="Atlas validation legend">
      <strong>Atlas Region Validation</strong><span>Deterministic colors distinguish region IDs; they do not encode disease values.</span>
    </div>
  )
  if (metric === 'fdr_p') return (
    <div className="legend" aria-label="FDR significance legend">
      <strong>FDR significance</strong>
      <div className="binary-legend"><i className="sig" /> q &lt; {threshold} <i className="nonsig" /> not significant</div>
      <span>{significant} of 138 regions meet the predefined threshold.</span>
    </div>
  )
  const observed = metric === 'observed_score'
  const domain = observed ? observedDomain : zDomain
  return (
    <div className="legend" aria-label={`${observed ? 'Observed expression' : 'Z-score'} color legend`}>
      <strong>{observed ? 'Observed expression' : 'Z-score'}</strong>
      <div className={`gradient-key ${observed ? 'sequential' : 'diverging'}`} />
      <div className="legend-values"><span>{formatValue(domain[0], 2)}</span>{!observed && <span>0</span>}<span>{formatValue(domain[1], 2)}</span></div>
      <span>{observed ? 'Lower to higher normalized expression score' : 'Blue: below null · Red: above null'}</span>
    </div>
  )
}
