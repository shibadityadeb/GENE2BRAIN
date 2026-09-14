import { lazy, Suspense, useEffect, useMemo, useState } from 'react'
import { Controls } from './components/Controls'
import { Legend } from './components/Legend'
import { RegionPanel } from './components/RegionPanel'
import { ResearchSections } from './components/ResearchSections'
import { Tooltip } from './components/Tooltip'
import { loadResearchData } from './data'
import type { AtlasGeometry, EnrichmentData, Hemisphere, Metric, ProjectMetadata, RegionRecord, ViewPreset } from './types'

type Loaded = { enrichment: EnrichmentData; metadata: ProjectMetadata; geometry: AtlasGeometry }
const BrainScene = lazy(() => import('./components/BrainScene').then((module) => ({ default: module.BrainScene })))

export default function App() {
  const [loaded, setLoaded] = useState<Loaded | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [metric, setMetric] = useState<Metric>('z_score')
  const [hemisphere, setHemisphere] = useState<Hemisphere>('whole')
  const [selected, setSelected] = useState<RegionRecord | null>(null)
  const [hovered, setHovered] = useState<RegionRecord | null>(null)
  const [tooltipPoint, setTooltipPoint] = useState({ x: 0, y: 0 })
  const [focusNonce, setFocusNonce] = useState(0)
  const [viewPreset, setViewPreset] = useState<ViewPreset>('reset')
  const validationMode = useMemo(() => new URLSearchParams(window.location.search).get('mode') === 'atlas', [])

  useEffect(() => {
    loadResearchData().then(setLoaded).catch((reason: unknown) => setError(reason instanceof Error ? reason.message : 'Unknown data error'))
  }, [])

  if (error) return <main className="load-state error-state"><h1>Visualization unavailable</h1><p>{error}</p><p>No scientific values have been displayed.</p></main>
  if (!loaded) return <main className="load-state"><div className="loader" /><p>Loading verified AAL3 atlas geometry…</p></main>

  const { enrichment, metadata, geometry } = loaded
  const zDomain = metadata.metrics.z_score.domain as [number, number]
  const observedDomain = metadata.metrics.observed_score.domain as [number, number]
  const selectRegion = (region: RegionRecord) => {
    setHemisphere('whole')
    setSelected(region)
    setFocusNonce((value) => value + 1)
  }
  const setView = (view: ViewPreset) => {
    setViewPreset(view)
    setSelected(null)
    setFocusNonce(0)
  }

  return (
    <div className="app-shell">
      <header className="topbar">
        <a className="brand" href="#top" aria-label="GENE2BRAIN home"><strong>GENE<span>2</span>BRAIN</strong><small>Spatial genetic enrichment</small></a>
        <nav aria-label="Primary navigation"><a href="#how">How it works</a><a href="#data">Data</a><a href="#methods">Methods</a></nav>
      </header>

      <section className="hero" id="top">
        <div className="hero-heading">
          <div><p className="eyebrow">Interactive research atlas · Parkinson disease</p><h1>From Genetic Risk to<br /><em>Spatial Brain Vulnerability</em></h1></div>
          <p>Explore where prioritized disease-associated genes show unusually high or low expression across healthy human brain regions.</p>
        </div>
        <Controls
          metric={metric}
          hemisphere={hemisphere}
          regions={enrichment.regions}
          selected={selected}
          onMetric={setMetric}
          onHemisphere={setHemisphere}
          onSelect={selectRegion}
          onView={setView}
        />
        <div className={`brain-stage ${selected ? 'has-panel' : ''}`}>
          <Suspense fallback={<div className="scene-loading">Preparing WebGL brain…</div>}><BrainScene
            geometryData={geometry}
            regions={enrichment.regions}
            metric={metric}
            hemisphere={hemisphere}
            selected={selected}
            focusNonce={focusNonce}
            viewPreset={viewPreset}
            validationMode={validationMode}
            zDomain={zDomain}
            observedDomain={observedDomain}
            fdrThreshold={metadata.analysis.fdr_threshold}
            onHover={(region, point) => { setHovered(region); if (point) setTooltipPoint(point) }}
            onSelect={selectRegion}
            onReset={() => setView('reset')}
          /></Suspense>
          <div className="stage-label"><span>{validationMode ? 'ATLAS REGION VALIDATION' : 'AAL3 · 138 ANALYZED REGIONS'}</span><span>Drag to rotate · scroll/pinch to zoom · right-drag to pan</span></div>
          <Legend
            metric={metric}
            zDomain={zDomain}
            observedDomain={observedDomain}
            threshold={metadata.analysis.fdr_threshold}
            significant={metadata.analysis.significant_regions}
            spatialRobustRegions={metadata.analysis.spatial_sensitivity.robust_regions}
            validationMode={validationMode}
          />
          {validationMode && (
            <div className="validation-list" aria-label="Atlas region list"><strong>Region ID map</strong>{[...enrichment.regions].sort((a, b) => a.region_id - b.region_id).map((region) => (
              <button key={region.region_id} onClick={() => selectRegion(region)}><span>{region.region_id}</span>{region.region_name}</button>
            ))}</div>
          )}
          {selected && <RegionPanel region={selected} threshold={metadata.analysis.fdr_threshold} onClose={() => setSelected(null)} />}
        </div>
      </section>
      {hovered && <Tooltip region={hovered} point={tooltipPoint} />}
      <ResearchSections metadata={metadata} />
      <footer><strong>GENE2BRAIN</strong><span>Research visualization · values generated by the validated scientific pipeline</span></footer>
    </div>
  )
}
