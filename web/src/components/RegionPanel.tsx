import { useState } from 'react'
import { formatValue } from '../format'
import type { RegionRecord } from '../types'

function NullSummary({ region }: { region: RegionRecord }) {
  const low = Math.min(region.observed_score, region.random_mean - 3 * region.random_std)
  const high = Math.max(region.observed_score, region.random_mean + 3 * region.random_std)
  const scale = (value: number) => 8 + 84 * (value - low) / (high - low || 1)
  return (
    <div className="null-summary">
      <svg viewBox="0 0 100 34" role="img" aria-label="Observed score compared with null mean plus or minus one standard deviation">
        <line x1={scale(region.random_mean - region.random_std)} x2={scale(region.random_mean + region.random_std)} y1="18" y2="18" />
        <circle className="null-dot" cx={scale(region.random_mean)} cy="18" r="3" />
        <line className="observed-mark" x1={scale(region.observed_score)} x2={scale(region.observed_score)} y1="7" y2="29" />
      </svg>
      <div className="null-key"><span>● Null mean ± 1 SD</span><span>│ Observed</span></div>
      <p>Individual permutation draws were not retained by Stage 6, so this shows the real stored summary—not a reconstructed distribution.</p>
    </div>
  )
}

export function RegionPanel({ region, threshold, onClose }: { region: RegionRecord; threshold: number; onClose: () => void }) {
  const [showNull, setShowNull] = useState(false)
  const direction = region.z_score > 0 ? 'above' : region.z_score < 0 ? 'below' : 'at'
  return (
    <aside className="region-panel" aria-live="polite" aria-label={`Details for ${region.region_name}`}>
      <button className="close-panel" onClick={onClose} aria-label="Close region details">×</button>
      <span className="eyebrow">AAL3 · {region.atlas_id}</span>
      <h2>{region.region_name}</h2>
      <p className="disease-label">Parkinson disease · weighted gene set</p>
      <h3 className="panel-section-title">Stage 6 · Gene-set enrichment</h3>
      <div className="stat-grid">
        <div><span>Z-score</span><strong>{formatValue(region.z_score)}</strong></div>
        <div><span>FDR q-value</span><strong>{formatValue(region.fdr_p, 4)}</strong></div>
        <div><span>Observed score</span><strong>{formatValue(region.observed_score, 4)}</strong></div>
        <div><span>Null mean</span><strong>{formatValue(region.random_mean, 4)}</strong></div>
        <div><span>Null SD</span><strong>{formatValue(region.random_std, 4)}</strong></div>
        <div><span>Effect size</span><strong>{formatValue(region.effect_size, 4)}</strong></div>
        <div><span>Genes in score</span><strong>{region.number_of_genes}</strong></div>
        <div><span>FDR status</span><strong>{region.fdr_p < threshold ? 'Significant' : 'Not significant'}</strong></div>
      </div>
      <h3 className="panel-section-title">Stage 7 · Spatial sensitivity</h3>
      <div className="stat-grid">
        <div><span>Spatial percentile</span><strong>{formatValue(region.spatial_robustness, 4)}</strong></div>
        <div><span>Spatial p-value</span><strong>{formatValue(region.spatial_null_p, 4)}</strong></div>
        <div><span>Spatial FDR</span><strong>{formatValue(region.spatial_null_fdr, 4)}</strong></div>
        <div><span>Joint result</span><strong>{region.spatial_robustness_label === 'robust' ? 'Robust' : 'Not robust'}</strong></div>
        <div><span>Robustness rank</span><strong>{region.robustness_rank}</strong></div>
        <div><span>Graph status</span><strong>{region.spatial_isolate ? 'Isolated parcel' : 'Connected parcel'}</strong></div>
      </div>
      <h3 className="panel-section-title">Stage 8 · Independent validation data</h3>
      <div className="stat-grid">
        <div><span>Validation score</span><strong>{formatValue(region.validation_score)}</strong></div>
        <div><span>Agreement</span><strong>{region.agreement_status === 'not_measured' ? 'Not measured' : region.agreement_status.replaceAll('_', ' ')}</strong></div>
        <div><span>ENIGMA parcel</span><strong>{region.validation_region ?? 'N/A'}</strong></div>
        <div><span>Mapping confidence</span><strong>{region.validation_mapping_confidence ?? 'N/A'}</strong></div>
      </div>
      <p className="caution">The ENIGMA-PD score is external to the discovery model. Higher values mean thinner cortex or smaller subcortical volume in PD; N/A means the phenotype did not measure this AAL3 parcel.</p>
      <p className="interpretation">The observed Parkinson-associated gene-expression score is {formatValue(Math.abs(region.z_score))} standard deviations {direction} the matched gene-set expectation. It {region.fdr_p < threshold ? 'meets' : 'does not meet'} the project’s FDR threshold.</p>
      <p className="caution">This does not indicate where Parkinson disease occurs or establish a causal brain region.</p>
      <button className="outline-button" onClick={() => setShowNull((visible) => !visible)} aria-expanded={showNull}>
        {showNull ? 'Hide null summary' : 'Show null distribution'}
      </button>
      {showNull && <NullSummary region={region} />}
    </aside>
  )
}
