import { useEffect, useState } from 'react'
import type { Hemisphere, Metric, RegionRecord, ViewPreset } from '../types'

export function Controls({
  metric, hemisphere, regions, selected, onMetric, onHemisphere, onSelect, onView,
}: {
  metric: Metric
  hemisphere: Hemisphere
  regions: RegionRecord[]
  selected: RegionRecord | null
  onMetric: (metric: Metric) => void
  onHemisphere: (hemisphere: Hemisphere) => void
  onSelect: (region: RegionRecord) => void
  onView: (view: ViewPreset) => void
}) {
  const [query, setQuery] = useState('')
  useEffect(() => { setQuery(selected?.region_name ?? '') }, [selected])
  return (
    <div className="control-stack">
      <div className="control-row primary-controls">
        <label>Disease
          <select value="parkinson-disease" aria-label="Disease">
            <option value="parkinson-disease">Parkinson disease</option>
          </select>
        </label>
        <label>Metric
          <select value={metric} onChange={(event) => onMetric(event.target.value as Metric)} aria-label="Metric">
            <option value="z_score">Z-score</option>
            <option value="observed_score">Observed expression</option>
            <option value="fdr_p">FDR significance</option>
          </select>
        </label>
        <label className="search-control">Search brain region
          <input
            list="region-options"
            placeholder="Putamen, caudate…"
            value={query}
            onChange={(event) => {
              setQuery(event.target.value)
              const match = regions.find((region) => region.region_name.toLowerCase() === event.target.value.toLowerCase())
              if (match) onSelect(match)
            }}
            aria-label="Search brain region"
          />
          <datalist id="region-options">
            {regions.map((region) => <option key={region.region_id} value={region.region_name} />)}
          </datalist>
        </label>
      </div>
      <div className="control-row compact-controls" aria-label="Anatomical navigation">
        <div className="segmented">
          {([['whole', 'Whole'], ['L', 'Left'], ['R', 'Right']] as const).map(([value, label]) => (
            <button key={value} className={hemisphere === value ? 'active' : ''} onClick={() => onHemisphere(value)}>{label}</button>
          ))}
        </div>
        <div className="view-buttons">
          {(['anterior', 'posterior', 'superior', 'inferior'] as ViewPreset[]).map((view) => (
            <button key={view} onClick={() => onView(view)}>{view[0].toUpperCase() + view.slice(1)}</button>
          ))}
          <button className="reset-button" onClick={() => onView('reset')} aria-label="Reset camera">Reset view</button>
        </div>
      </div>
    </div>
  )
}
