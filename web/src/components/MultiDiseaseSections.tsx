import { useMemo, useState } from 'react'
import { BrainScene } from './BrainScene'
import type { AtlasGeometry, EnrichmentData, MultidiseaseAtlas } from '../types'

function correlation(first: number[], second: number[]) {
  const meanA = first.reduce((sum, value) => sum + value, 0) / first.length
  const meanB = second.reduce((sum, value) => sum + value, 0) / second.length
  const numerator = first.reduce((sum, value, index) => sum + (value - meanA) * (second[index] - meanB), 0)
  const denominator = Math.sqrt(
    first.reduce((sum, value) => sum + (value - meanA) ** 2, 0)
    * second.reduce((sum, value) => sum + (value - meanB) ** 2, 0),
  )
  return numerator / denominator
}

function MiniBrain({ disease, geometry, zDomain }: { disease: EnrichmentData; geometry: AtlasGeometry; zDomain: [number, number] }) {
  const observed = disease.regions.map((region) => region.observed_score)
  return <div className="compare-brain"><BrainScene
    geometryData={geometry} regions={disease.regions} metric="z_score" hemisphere="whole"
    selected={null} focusNonce={0} viewPreset="reset" validationMode={false}
    zDomain={zDomain} observedDomain={[Math.min(...observed), Math.max(...observed)]}
    validationDomain={[-1, 1]} fdrThreshold={0.05}
    onHover={() => undefined} onSelect={() => undefined} onReset={() => undefined}
  /></div>
}

export function MultiDiseaseSections({ atlas, geometry }: { atlas: MultidiseaseAtlas; geometry: AtlasGeometry }) {
  const [firstId, setFirstId] = useState(atlas.diseases[0].disease_id)
  const [secondId, setSecondId] = useState(atlas.diseases[1].disease_id)
  const [regionId, setRegionId] = useState(atlas.diseases[0].regions[0].region_id)
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
      r: correlation(rows.map((row) => row.first), rows.map((row) => row.second)),
      rows: rows.sort((a, b) => Math.abs(b.difference) - Math.abs(a.difference)),
    }
  }, [first, second])
  const regionName = first.regions.find((region) => region.region_id === regionId)!.region_name
  const ranked = atlas.diseases.map((disease) => ({
    disease: disease.disease_name,
    z: disease.regions.find((region) => region.region_id === regionId)!.z_score,
    fdr: disease.regions.find((region) => region.region_id === regionId)!.fdr_p,
  })).sort((a, b) => b.z - a.z)

  return <div className="multidisease-content">
    <section id="compare" className="compare-section">
      <p className="section-kicker">Cross-disease analysis</p>
      <h2>Compare Diseases</h2>
      <div className="compare-controls">
        <label>Disease A<select value={firstId} onChange={(event) => setFirstId(event.target.value)}>{atlas.diseases.map((disease) => <option key={disease.disease_id} value={disease.disease_id}>{disease.disease_name}</option>)}</select></label>
        <label>Disease B<select value={secondId} onChange={(event) => setSecondId(event.target.value)}>{atlas.diseases.map((disease) => <option key={disease.disease_id} value={disease.disease_id}>{disease.disease_name}</option>)}</select></label>
        <div className="correlation-card"><span>Regional Pearson correlation</span><strong>{comparison.r.toFixed(3)}</strong><small>138 matched-null AAL3 Z scores</small></div>
      </div>
      <div className="brain-comparison">
        <article><h3>{first.disease_name}</h3><MiniBrain disease={first} geometry={geometry} zDomain={atlas.z_domain} /></article>
        <article><h3>{second.disease_name}</h3><MiniBrain disease={second} geometry={geometry} zDomain={atlas.z_domain} /></article>
      </div>
      <p className="comparison-note">Both brains use the same global color scale. The tabulated difference is Z<sub>A</sub> − Z<sub>B</sub> in the same parcel. It is a descriptive difference between standardized matched-null scores, not a formal test that diseases differ.</p>
      <div className="difference-table"><strong>Largest absolute standardized differences</strong>{comparison.rows.slice(0, 8).map((row) => <div key={row.region}><span>{row.region}</span><b>{row.difference.toFixed(2)}</b></div>)}</div>
    </section>

    <section id="atlas" className="master-atlas-section">
      <p className="section-kicker">Region-first view</p>
      <h2>GENE2BRAIN Atlas</h2>
      <label>Brain region<select value={regionId} onChange={(event) => setRegionId(Number(event.target.value))}>{first.regions.map((region) => <option key={region.region_id} value={region.region_id}>{region.region_name}</option>)}</select></label>
      <h3>{regionName}</h3>
      <div className="region-ranking">{ranked.map((row, index) => <div key={row.disease}><span>{index + 1}</span><strong>{row.disease}</strong><b>Z {row.z.toFixed(2)}</b><small>FDR {row.fdr.toPrecision(3)}</small></div>)}</div>
      <div className="downloads">
        <a className="download-button" href="./data/disease_region_enrichment_matrix.csv" download>Download enrichment matrix</a>
        <a className="outline-button" href="./data/disease_spatial_similarity.csv" download>Download similarity table</a>
        <a className="outline-button" href="./data/shared_brain_region_enrichment.csv" download>Download shared regions</a>
      </div>
      <details className="excluded-list"><summary>Diseases not shown ({atlas.excluded_or_needs_review.length})</summary>{atlas.excluded_or_needs_review.map((item) => <p key={item.disease_id}><strong>{item.disease_name} · {item.status}</strong><br />{item.reason}</p>)}</details>
    </section>

    <section id="story" className="story-stage10">
      <p className="section-kicker">Research story</p><h2>From Genetic Risk to Brain</h2>
      <div className="stage10-flow">{['GWAS', 'Loci', 'Genes', 'Healthy brain', 'Regional expression', 'Enrichment', 'Spatial robustness', 'Independent validation', 'Biological interpretation', 'Cross-disease atlas'].map((step, index) => <div key={step}><span>{String(index + 1).padStart(2, '0')}</span><strong>{step}</strong></div>)}</div>
      <p className="comparison-note">15 diseases were screened; 10 passed all genetic and AHBA gates. Every displayed map contains 138 regions, uses 10,000 matched gene-set permutations and 10,000 spatial permutations, and retains null results. No map is a diagnosis, patient prediction, or causal anatomical claim.</p>
    </section>
  </div>
}
