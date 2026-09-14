import { formatCount } from '../format'
import type { ProjectMetadata } from '../types'

const pipeline = [
  'GWAS', 'Disease-associated loci', 'Gene prioritization', 'Allen Human Brain Atlas',
  'Regional gene expression', 'Matched random gene sets', 'Enrichment statistics', '3D brain map',
]

export function ResearchSections({ metadata }: { metadata: ProjectMetadata }) {
  const counts = metadata.counts
  return (
    <main className="research-content">
      <section id="how" className="method-section">
        <p className="section-kicker">Research context</p>
        <h2>How this map is generated</h2>
        <div className="pipeline" aria-label={pipeline.join(' then ')}>
          {pipeline.map((step, index) => <div key={step}><span>{String(index + 1).padStart(2, '0')}</span>{step}</div>)}
        </div>
        <p className="lead">The colors represent how strongly Parkinson-associated genes are expressed in each healthy brain region relative to matched random gene sets.</p>
        <p className="scientific-caution">This is a genetically informed spatial enrichment map, not a direct map of where disease pathology occurs.</p>
      </section>

      <section className="story-section">
        <p className="section-kicker">From Genetic Risk to Brain</p>
        <div className="story-flow">
          <article><span>01</span><h3>GWAS</h3><strong>{formatCount(counts.gwas_sample_size)}</strong><p>participants in selected multi-ancestry study {metadata.analysis.gwas_accession}</p></article>
          <article><span>02</span><h3>Prioritized genes</h3><strong>{formatCount(counts.parkinson_genes_represented)}</strong><p>weighted L2G genes represented in the AHBA matrix, from {formatCount(counts.credible_sets)} credible sets</p></article>
          <article><span>03</span><h3>Brain enrichment</h3><strong>{formatCount(counts.regions_analyzed)}</strong><p>AAL3 regions compared with {formatCount(counts.random_gene_sets)} matched random gene sets</p></article>
        </div>
      </section>

      <section id="data" className="data-methods-grid">
        <div>
          <p className="section-kicker">Data</p>
          <h2>Research inputs</h2>
          <dl className="data-list">
            <div><dt>AHBA donors</dt><dd>{formatCount(counts.ahba_donors)}</dd></div>
            <div><dt>AHBA tissue samples</dt><dd>{formatCount(counts.ahba_samples)}</dd></div>
            <div><dt>AAL3 parcels / analyzed</dt><dd>{formatCount(counts.atlas_parcels)} / {formatCount(counts.regions_analyzed)}</dd></div>
            <div><dt>AHBA genes analyzed</dt><dd>{formatCount(counts.genes_analyzed)}</dd></div>
            <div><dt>Parkinson genes represented</dt><dd>{formatCount(counts.parkinson_genes_represented)}</dd></div>
            <div><dt>Matched random gene sets</dt><dd>{formatCount(counts.random_gene_sets)}</dd></div>
            <div><dt>Significance threshold</dt><dd>FDR q &lt; {metadata.analysis.fdr_threshold}</dd></div>
          </dl>
          <div className="downloads">
            <a className="download-button" href="./data/parkinson_regional_enrichment.csv" download>Download regional results</a>
            <a className="outline-button" href="./data/project_metadata.json" download>Download analysis metadata</a>
          </div>
        </div>
        <div id="methods">
          <p className="section-kicker">Methods &amp; Sources</p>
          <h2>Traceable by design</h2>
          <h3>Data sources</h3>
          <ul className="source-list">
            {metadata.sources.map((source) => <li key={source.name}><a href={source.url} target="_blank" rel="noreferrer">{source.name}</a></li>)}
          </ul>
          <h3>Analysis performed by GENE2BRAIN</h3>
          <p>GENE2BRAIN selected the GWAS, prioritized genes with {metadata.analysis.gene_prioritization}, processed AHBA expression into AAL3 parcels, and tested regional scores against expression-, variance-, and probe-count-matched gene sets.</p>
          <p>Source organizations provide data and methods; their inclusion does not imply endorsement of GENE2BRAIN.</p>
        </div>
      </section>
    </main>
  )
}
