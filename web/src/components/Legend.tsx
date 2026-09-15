import type { Metric } from '../types'
import { formatValue } from '../format'

export function Legend({ metric, zDomain, observedDomain, validationDomain, threshold, significant, spatialRobustRegions, mappedValidationRegions, validationMode }: {
  metric: Metric
  zDomain: [number, number]
  observedDomain: [number, number]
  validationDomain: [number, number]
  threshold: number
  significant: number
  spatialRobustRegions: number
  mappedValidationRegions: number
  validationMode: boolean
}) {
  if (validationMode) return (
    <div className="legend" aria-label="Atlas validation legend">
      <strong>Atlas Region Validation</strong><span>Deterministic colors distinguish region IDs; they do not encode disease values.</span>
    </div>
  )
  if (metric === 'anatomy') return (
    <div className="legend" aria-label="Anatomy only legend">
      <strong>Anatomy only</strong><span>Folded fsaverage6 pial cortical surface. No research values or region colors are shown.</span>
    </div>
  )
  if (metric === 'fdr_p') return (
    <div className="legend" aria-label="FDR significance legend">
      <strong>FDR significance</strong>
      <div className="binary-legend"><i className="sig" /> q &lt; {threshold} <i className="nonsig" /> not significant</div>
      <span>{significant} of 138 regions meet the predefined threshold.</span>
    </div>
  )
  if (metric === 'spatial_robustness') return (
    <div className="legend" aria-label="Spatial robustness legend">
      <strong>Spatial robustness</strong>
      <div className="gradient-key spatial" />
      <div className="legend-values"><span>0</span><span>Spatial-null percentile</span><span>1</span></div>
      <span>{spatialRobustRegions} regions meet the joint Stage 6 + Stage 7 rule. A high percentile alone is not a robust call.</span>
    </div>
  )
  if (metric === 'validation_score') return (
    <div className="legend" aria-label="Independent validation legend">
      <strong>Independent validation data</strong>
      <div className="gradient-key diverging" />
      <div className="legend-values"><span>{formatValue(validationDomain[0], 2)}</span><span>0</span><span>{formatValue(validationDomain[1], 2)}</span></div>
      <span>Higher values mean thinner cortex or smaller subcortical volume in PD. Gray: not measured ({mappedValidationRegions} AAL3 parcels mapped).</span>
    </div>
  )
  if (metric === 'agreement') return (
    <div className="legend" aria-label="Regional agreement legend">
      <strong>Agreement · median split</strong>
      <div className="agreement-legend"><i className="hh" /> high / high <i className="hl" /> high / low <i className="lh" /> low / high <i className="ll" /> low / low</div>
      <span>Descriptive categories only. Gray parcels were not measured.</span>
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
