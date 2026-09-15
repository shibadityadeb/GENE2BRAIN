import { lazy, Suspense, useMemo, useState } from 'react'
import type { AtlasGeometry, EnrichmentData, MultidiseaseAtlas } from '../types'

const BrainScene = lazy(() => import('./BrainScene').then((module) => ({ default: module.BrainScene })))

function MiniBrain({ disease, geometry, zDomain }: { disease: EnrichmentData; geometry: AtlasGeometry; zDomain: [number, number] }) {
  const observed = disease.regions.map((region) => region.observed_score)
  return <div className="compare-brain"><Suspense fallback={<div className="scene-loading">Preparing comparison brain…</div>}><BrainScene
    geometryData={geometry} regions={disease.regions} metric="z_score" hemisphere="whole"
    selected={null} focusNonce={0} viewPreset="reset" validationMode={false}
    zDomain={zDomain} observedDomain={[Math.min(...observed), Math.max(...observed)]}
    validationDomain={[-1, 1]} fdrThreshold={0.05}
    onHover={() => undefined} onSelect={() => undefined} onReset={() => undefined}
  /></Suspense></div>
}

export function MultiDiseaseSections({ atlas, geometry }: { atlas: MultidiseaseAtlas; geometry: AtlasGeometry }) {
  const [firstId, setFirstId] = useState(atlas.diseases[0].disease_id)
  const [secondId, setSecondId] = useState(atlas.diseases[1].disease_id)
  const [regionId, setRegionId] = useState(atlas.diseases[0].regions[0].region_id)
  const [showBrains, setShowBrains] = useState(false)
  const first = atlas.diseases.find((disease) => disease.disease_id === firstId)!
  const second = atlas.diseases.find((disease) => disease.disease_id === secondId)!
  const comparison = useMemo(() => {
    const secondById = new Map(second.regions.map((region) => [region.region_id, region]))
    const rows = first.regions.map((region) => ({
      region: region.region_name,
      first: region.z_score,
      second: secondById.get(region.region_id)!.z_score,
      difference: region.z_score - secondById.get(region.region_id)!.z_score,
    }))
    return {
      rows: rows.sort((a, b) => Math.abs(b.difference) - Math.abs(a.difference)),
    }
  }, [first, second])
  const similarity = atlas.pearson_similarity.find((row) => row.disease_1 === first.disease_name && row.disease_2 === second.disease_name)
    ?? atlas.pearson_similarity.find((row) => row.disease_1 === second.disease_name && row.disease_2 === first.disease_name)
  if (!similarity) throw new Error(`Missing frozen Pearson similarity for ${first.disease_name} and ${second.disease_name}`)
  const regionName = first.regions.find((region) => region.region_id === regionId)!.region_name
  const topRegions = [...first.regions].sort((a, b) => b.z_score - a.z_score).slice(0, 5)
  const lowRegions = [...first.regions].sort((a, b) => a.z_score - b.z_score).slice(0, 5)
  const ranked = atlas.diseases.map((disease) => ({
    disease: disease.disease_name,
    z: disease.regions.find((region) => region.region_id === regionId)!.z_score,
    fdr: disease.regions.find((region) => region.region_id === regionId)!.fdr_p,
  })).sort((a, b) => b.z - a.z)

  return <div className="multidisease-content">
    <section id="results" className="compare-section">
      <p className="section-kicker">What we found</p><h2>What does the current map show?</h2><p className="lead">Parkinson disease is the default worked example. The map summarizes where Parkinson-linked genes show stronger or lower combined activity across 138 healthy brain regions. Explore a region to see the evidence and its limits.</p>
      <p className="section-kicker">Compare diseases</p>
      <h2>Compare Diseases</h2>
      <div id="compare" />
      <div className="compare-controls">
        <label>Disease A<select value={firstId} onChange={(event) => setFirstId(event.target.value)}>{atlas.diseases.map((disease) => <option key={disease.disease_id} value={disease.disease_id}>{disease.disease_name}</option>)}</select></label>
        <label>Disease B<select value={secondId} onChange={(event) => setSecondId(event.target.value)}>{atlas.diseases.map((disease) => <option key={disease.disease_id} value={disease.disease_id}>{disease.disease_name}</option>)}</select></label>
        <div className="correlation-card"><span>How similar are the brain patterns?</span><strong>{similarity.correlation.toFixed(3)}</strong><small>Correlation across {similarity.n_regions} shared brain regions.</small></div>
      </div>
      <button className="compare-toggle" onClick={() => setShowBrains((value) => !value)}>{showBrains ? 'Hide 3D comparison brains' : 'Load 3D comparison brains'}</button>
      {showBrains ? <div className="brain-comparison">
        <article><h3>{first.disease_name}</h3><MiniBrain disease={first} geometry={geometry} zDomain={atlas.z_domain} /></article>
        <article><h3>{second.disease_name}</h3><MiniBrain disease={second} geometry={geometry} zDomain={atlas.z_domain} /></article>
      </div> : <p className="comparison-note">The paired 3D views load on demand to keep the main atlas responsive. The correlation and ranked differences above are always available from the frozen regional tables.</p>}
      <p className="comparison-note">Both brains use the same color scale. The table shows where the two stored regional patterns differ most. This is a descriptive comparison, not proof that one disease affects a region more than another.</p>
      <div className="difference-table"><strong>Largest absolute standardized differences</strong>{comparison.rows.slice(0, 8).map((row) => <div key={row.region}><span>{row.region}</span><b>{row.difference.toFixed(2)}</b></div>)}</div>
      <div className="regional-highlights">
        <div><h3>{first.disease_name} · strongest signals</h3>{topRegions.map((region) => <p key={region.region_id}><span>{region.region_name}</span><b>Higher</b><small>{region.z_score.toFixed(2)} standardized score</small></p>)}</div>
        <div><h3>{first.disease_name} · lower signals</h3>{lowRegions.map((region) => <p key={region.region_id}><span>{region.region_name}</span><b>Lower</b><small>{region.z_score.toFixed(2)} standardized score</small></p>)}</div>
      </div>
      <p className="comparison-note">These rankings summarize the project’s stored statistical evidence. They are not maps of pathology or individual disease risk.</p>
    </section>

    <section id="atlas" className="master-atlas-section">
      <p className="section-kicker">Explore a region across diseases</p>
      <h2>Which diseases show a signal here?</h2>
      <label>Brain region<select value={regionId} onChange={(event) => setRegionId(Number(event.target.value))}>{first.regions.map((region) => <option key={region.region_id} value={region.region_id}>{region.region_name}</option>)}</select></label>
      <h3>{regionName}</h3>
      <div className="region-ranking">{ranked.map((row, index) => <div key={row.disease}><span>{index + 1}</span><strong>{row.disease}</strong><b>{row.z > 0 ? 'Higher signal' : 'Lower signal'}</b><small>{row.z.toFixed(2)} standardized score</small></div>)}</div>
      <div className="downloads">
        <a className="download-button" href="/data/disease_region_enrichment_matrix.csv" download>Download enrichment matrix</a>
        <a className="outline-button" href="/data/disease_spatial_similarity.csv" download>Download similarity table</a>
        <a className="outline-button" href="/data/shared_brain_region_enrichment.csv" download>Download shared regions</a>
      </div>
      <details className="excluded-list"><summary>Diseases not shown ({atlas.excluded_or_needs_review.length})</summary>{atlas.excluded_or_needs_review.map((item) => <p key={item.disease_id}><strong>{item.disease_name} · {item.status}</strong><br />{item.reason}</p>)}</details>
    </section>

    <section id="limitations" className="story-stage10"><p className="section-kicker">Important context</p><h2>What this map does not mean</h2><p className="comparison-note">This is not a diagnostic tool, a prediction of an individual’s disease, or proof that a highlighted region causes disease. Gene activity in healthy tissue is not the same as measuring disease damage. Some diseases and regions have stronger evidence than others.</p></section>
  </div>
}
