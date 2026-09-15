import { formatCount } from '../format'
import type { ProjectMetadata } from '../types'

const pipeline = [
  'GWAS', 'Disease-associated loci', 'Gene prioritization', 'Allen Human Brain Atlas',
  'Regional gene expression', 'Matched random gene sets', 'Enrichment statistics', '3D brain map',
  'Independent PD phenotype validation', 'Biological interpretation',
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
        <p className="lead">The selected disease map shows how strongly its prioritized genes are expressed in each healthy brain region relative to matched random gene sets. The detailed audit below uses Parkinson disease as the frozen reference analysis.</p>
        <p className="scientific-caution">This is a genetically informed spatial enrichment map, not a direct map of where disease pathology occurs.</p>
      </section>

      <section id="before-after" className="before-after-section">
        <p className="section-kicker">Visual method redesign</p>
        <h2>From atlas voxels to folded anatomy</h2>
        <p>The old diagnostic view displayed the outer faces of 2 mm AAL3 labels as block-like parcels. The new diagnostic view shows the same region IDs projected onto an actual pial cortical surface. Neither image is a patient brain or a disease pathology map.</p>
        <div className="before-after-grid">
          <figure><img src="/images/web_brain_before.png" alt="Before: block-like exposed AAL3 voxel-face visualization" loading="lazy" /><figcaption>Before · 2 mm atlas voxel faces</figcaption></figure>
          <figure><img src="/images/web_brain_after.png" alt="After: folded pial cortical anatomy with diagnostic atlas-ID colors" loading="lazy" /><figcaption>After · fsaverage6 pial cortex with AAL3 ID projection</figcaption></figure>
        </div>
        <p>Research values are overlaid separately in the interactive explorer; the pastel diagnostic colors above encode IDs only. Projected borders are approximate at the 2 mm atlas resolution.</p>
      </section>

      <section className="story-section">
        <p className="section-kicker">From Genetic Risk to Brain</p>
        <div className="story-flow">
          <article><span>01</span><h3>GWAS</h3><strong>{formatCount(counts.gwas_sample_size)}</strong><p>participants in selected multi-ancestry study {metadata.analysis.gwas_accession}</p></article>
          <article><span>02</span><h3>Prioritized genes</h3><strong>{formatCount(counts.parkinson_genes_represented)}</strong><p>weighted L2G genes represented in the AHBA matrix, from {formatCount(counts.credible_sets)} credible sets</p></article>
          <article><span>03</span><h3>Brain enrichment</h3><strong>{formatCount(counts.regions_analyzed)}</strong><p>AAL3 regions compared with {formatCount(counts.random_gene_sets)} matched random gene sets</p></article>
          <article><span>04</span><h3>External validation</h3><strong>{metadata.analysis.independent_validation.pearson_r.toFixed(3)}</strong><p>primary Pearson r across {formatCount(metadata.analysis.independent_validation.matched_validation_units)} unique ENIGMA parcels · {metadata.analysis.independent_validation.interpretation.toLowerCase()}</p></article>
        </div>
      </section>

      <section id="data" className="data-methods-grid">
        <div>
          <p className="section-kicker">Parkinson reference analysis</p>
          <h2>Worked-example inputs</h2>
          <dl className="data-list">
            <div><dt>AHBA donors</dt><dd>{formatCount(counts.ahba_donors)}</dd></div>
            <div><dt>AHBA tissue samples</dt><dd>{formatCount(counts.ahba_samples)}</dd></div>
            <div><dt>AAL3 parcels / analyzed</dt><dd>{formatCount(counts.atlas_parcels)} / {formatCount(counts.regions_analyzed)}</dd></div>
            <div><dt>AHBA genes analyzed</dt><dd>{formatCount(counts.genes_analyzed)}</dd></div>
            <div><dt>Parkinson genes represented</dt><dd>{formatCount(counts.parkinson_genes_represented)}</dd></div>
            <div><dt>Matched random gene sets</dt><dd>{formatCount(counts.random_gene_sets)}</dd></div>
            <div><dt>Spatial null realizations</dt><dd>{formatCount(metadata.analysis.spatial_sensitivity.permutations)}</dd></div>
            <div><dt>Jointly robust regions</dt><dd>{formatCount(metadata.analysis.spatial_sensitivity.robust_regions)}</dd></div>
            <div><dt>Validation units / mapped AAL3</dt><dd>{formatCount(metadata.analysis.independent_validation.matched_validation_units)} / {formatCount(metadata.analysis.independent_validation.mapped_aal3_regions)}</dd></div>
            <div><dt>Primary validation</dt><dd>r = {metadata.analysis.independent_validation.pearson_r.toFixed(3)}; p = {metadata.analysis.independent_validation.pearson_p.toFixed(3)}</dd></div>
            <div><dt>Significant GO terms</dt><dd>{formatCount(metadata.analysis.biological_interpretation.go_significant_terms)}</dd></div>
            <div><dt>Significant Reactome pathways</dt><dd>{formatCount(metadata.analysis.biological_interpretation.reactome_significant_pathways)}</dd></div>
            <div><dt>Significant brain cell types</dt><dd>{formatCount(metadata.analysis.biological_interpretation.cell_types_significant)}</dd></div>
            <div><dt>Significance threshold</dt><dd>FDR q &lt; {metadata.analysis.fdr_threshold}</dd></div>
          </dl>
          <div className="downloads">
            <a className="download-button" href="/data/parkinson_regional_enrichment.csv" download>Download regional results</a>
            <a className="outline-button" href="/data/parkinson_atlas.json" download>Download Parkinson atlas JSON</a>
            <a className="outline-button" href="/data/parkinson_spatial_robustness.csv" download>Download spatial robustness</a>
            <a className="outline-button" href="/data/parkinson_independent_validation_regional_scores.csv" download>Download validation scores</a>
            <a className="outline-button" href="/data/parkinson_independent_validation_statistics.csv" download>Download validation statistics</a>
            <a className="outline-button" href="/data/stage_09_go_enrichment.csv" download>Download GO enrichment</a>
            <a className="outline-button" href="/data/stage_09_reactome_enrichment.csv" download>Download pathway enrichment</a>
            <a className="outline-button" href="/data/stage_09_cell_type_enrichment.csv" download>Download cell-type enrichment</a>
            <a className="outline-button" href="/data/stage_09_biological_evidence_summary.csv" download>Download biological evidence</a>
            <a className="outline-button" href="/data/project_metadata.json" download>Download analysis metadata</a>
            <a className="outline-button" href="/data/release.json" download>Download release manifest</a>
          </div>
        </div>
        <div id="methods">
          <p className="section-kicker">Reference methods &amp; sources</p>
          <h2>Parkinson audit trail</h2>
          <h3>Data sources</h3>
          <ul className="source-list">
            {metadata.sources.map((source) => <li key={source.name}><a href={source.url} target="_blank" rel="noreferrer">{source.name}</a></li>)}
          </ul>
          <h3>Discovery analysis · Stage 6</h3>
          <p>The gene-set null asks whether Parkinson-prioritized genes score higher than random genes matched for mean expression, expression variance, and probe count. It produces the displayed Z-scores, empirical p-values, and BH-FDR values.</p>
          <h3>Spatial sensitivity analysis · Stage 7</h3>
          <p>Neighboring brain regions are not independent. Stage 7 uses {metadata.analysis.spatial_sensitivity.neighbor_definition.toLowerCase()} and {metadata.analysis.spatial_sensitivity.method} to move enrichment peaks while preserving the observed atlas-graph spatial spectrum. It uses {formatCount(metadata.analysis.spatial_sensitivity.permutations)} surrogate maps and does not replace Stage 6.</p>
          <p>The weighted global peak test has spatial-null p = {metadata.analysis.spatial_sensitivity.global_p.toFixed(4)}. A region is only called robust when it passes both stages’ pre-specified FDR rules; {metadata.analysis.spatial_sensitivity.robust_regions} regions do so here.</p>
          <h3>Independent phenotype validation · Stage 8</h3>
          <p>{metadata.analysis.independent_validation.dataset} supplies {metadata.analysis.independent_validation.phenotype.toLowerCase()}. The AAL3 discovery parcels were mapped by an explicit anatomical crosswalk and averaged once within each ENIGMA validation parcel, preventing duplicated phenotype values from inflating the correlation sample size.</p>
          <p>Pearson r = {metadata.analysis.independent_validation.pearson_r.toFixed(3)} (95% region-bootstrap CI {metadata.analysis.independent_validation.pearson_ci[0].toFixed(3)} to {metadata.analysis.independent_validation.pearson_ci[1].toFixed(3)}, p = {metadata.analysis.independent_validation.pearson_p.toFixed(3)}); Spearman ρ = {metadata.analysis.independent_validation.spearman_rho.toFixed(3)} (p = {metadata.analysis.independent_validation.spearman_p.toFixed(3)}). The Stage-7-compatible spatial-null p is {metadata.analysis.independent_validation.spatial_null_p.toFixed(3)}. Overall interpretation: <strong>{metadata.analysis.independent_validation.interpretation}</strong>.</p>
          <p>The validation phenotype is an independent measurement and was not used in discovery. {metadata.analysis.independent_validation.independence_note}. {metadata.analysis.independent_validation.substantia_nigra}.</p>
          <h3>Biological interpretation · Stage 9</h3>
          <p><strong>Genetic evidence:</strong> the frozen Stage 4 L2G-weighted set is primary; broad and stringent definitions are sensitivities. L2G ranks prioritization evidence and is not a biological expression effect.</p>
          <p><strong>Spatial transcriptomic evidence:</strong> {metadata.analysis.biological_interpretation.genes_represented_in_ahba} of {metadata.analysis.biological_interpretation.genes_analyzed} prioritized genes are represented against {metadata.analysis.biological_interpretation.background.toLowerCase()}. The regional rule is {metadata.analysis.biological_interpretation.region_rule.toLowerCase()}.</p>
          <p><strong>Statistical enrichment:</strong> custom-background hypergeometric tests use GO and Reactome, with WikiPathways as a pathway sensitivity. Human Protein Atlas v25.1 single-nucleus brain profiles define 34 cell-type marker sets. BH correction is applied separately to GO ontologies, pathway databases, cell types, ranked pathways, and regional pathway tests.</p>
          <p><strong>Biological interpretation:</strong> {metadata.analysis.biological_interpretation.go_significant_terms} GO terms, {metadata.analysis.biological_interpretation.reactome_significant_pathways} primary Reactome pathways, and {metadata.analysis.biological_interpretation.cell_types_significant} cell types survive their respective corrections. {metadata.analysis.biological_interpretation.interpretation_note}.</p>
          <h3>Limitations</h3>
          <p>The adjacency graph encodes voxel contact, not neural connectivity. Four isolated small brainstem parcels are retained but less spatially constrained. Singleton MSR is conservative, the AHBA has six unevenly sampled donors, the atlas crosswalk includes composite mappings, ontologies contain overlapping terms, and ENIGMA structural MRI is a cross-sectional anatomical phenotype rather than direct histopathology. The maps do not establish pathology, causality, disease origin, pathway activation, temporal direction, or individual prediction.</p>
          <h3>Analysis performed by GENE2BRAIN</h3>
          <p>GENE2BRAIN selected the GWAS, prioritized genes with {metadata.analysis.gene_prioritization}, processed AHBA expression into AAL3 parcels, and performed the gene-set and spatial sensitivity analyses.</p>
          <p>Source organizations provide data and methods; their inclusion does not imply endorsement of GENE2BRAIN.</p>
        </div>
      </section>
      <section id="biology" className="biology-section">
        <p className="section-kicker">Biological interpretation</p>
        <h2>Evidence, not a pathology claim</h2>
        <div className="biology-cards">
          <article><span>Gene ontology</span><strong>{formatCount(metadata.analysis.biological_interpretation.go_significant_terms)}</strong><p>Terms surviving the ontology-specific false-discovery correction for the frozen Parkinson gene set.</p></article>
          <article><span>Reactome</span><strong>{formatCount(metadata.analysis.biological_interpretation.reactome_significant_pathways)}</strong><p>Primary pathway results under a custom AHBA-measured gene background.</p></article>
          <article><span>Brain cell types</span><strong>{formatCount(metadata.analysis.biological_interpretation.cell_types_significant)}</strong><p>Marker-set associations using Human Protein Atlas single-nucleus profiles; not cellular localization of pathology.</p></article>
        </div>
        <p className="biology-note">{metadata.analysis.biological_interpretation.interpretation_note} Region-level biology is shown in the brain detail panel only where the predeclared spatial robustness rule supports that interpretation.</p>
      </section>
      <section id="about" className="about-section">
        <p className="section-kicker">About the project</p>
        <h2>GENE2BRAIN</h2>
        <p className="lead">From Genetic Risk to Spatial Brain Vulnerability is an open, reproducible research visualization of disease-gene expression patterns in the healthy human brain.</p>
        <p>The analysis and visual website are developed under the GENE2BRAIN project. Researcher name, affiliation, contact address, and publication status have not been confirmed for this public release; no invented attribution is shown. Source organizations provide data and methods but do not endorse these interpretations.</p>
        <div id="author" className="author-profile" aria-label="Researcher profile pending confirmation"><h3>Researcher profile</h3><dl><div><dt>Name</dt><dd>Pending confirmation</dd></div><div><dt>Affiliation</dt><dd>Pending confirmation</dd></div><div><dt>Contact</dt><dd>Pending confirmation</dd></div><div><dt>Publication</dt><dd>Not confirmed</dd></div></dl></div>
        <p className="scientific-caution">These atlas maps are research summaries. They cannot diagnose patients or establish where disease begins.</p>
      </section>
      <section id="citations" className="citations-section">
        <p className="section-kicker">Sources and citation trail</p>
        <h2>Citations</h2>
        <p>Follow the source links and downloadable metadata for accessions, provenance, analysis versions, and the exact numerical data displayed. Cite primary data providers alongside any GENE2BRAIN analysis you use.</p>
        <ol className="bibliography">
          <li>Hawrylycz et al. (2012). <a href="https://doi.org/10.1038/nature11405" target="_blank" rel="noreferrer">An anatomically comprehensive atlas of the adult human brain transcriptome</a>. <em>Nature</em>. AHBA primary resource.</li>
          <li>Kim et al. (2023). <a href="https://doi.org/10.1038/s41588-023-01584-8" target="_blank" rel="noreferrer">Multi-ancestry genome-wide association meta-analysis of Parkinson’s disease</a>. <em>Nature Genetics</em>. <a href="https://www.ebi.ac.uk/gwas/studies/GCST90308590" target="_blank" rel="noreferrer">GWAS Catalog GCST90308590</a>.</li>
          <li>Rolls et al. (2020). <a href="https://doi.org/10.1016/j.neuroimage.2019.116189" target="_blank" rel="noreferrer">Automated Anatomical Labeling atlas 3</a>. <em>NeuroImage</em>. Atlas identity.</li>
          <li>Laansma et al. (2021). <a href="https://doi.org/10.1002/mds.28706" target="_blank" rel="noreferrer">An international multicenter analysis of brain structure across clinical stages of Parkinson’s disease</a>. <em>Movement Disorders</em>. External phenotype.</li>
          <li><a href="https://doi.org/10.1038/s41592-022-01625-w" target="_blank" rel="noreferrer">neuromaps surface-transform method</a> and <a href="https://nilearn.github.io/stable/modules/generated/nilearn.datasets.fetch_surf_fsaverage.html" target="_blank" rel="noreferrer">Nilearn fsaverage pial dataset</a>. Visual geometry and label projection.</li>
        </ol>
        <p className="asset-credit">3D geometry adapts the Nilearn-distributed fsaverage6 FreeSurfer pial template and projects AAL3 labels with neuromaps registration-fusion coordinates. The registration-fusion distribution is identified as <a href="https://creativecommons.org/licenses/by-nc-sa/4.0/" target="_blank" rel="noreferrer">CC BY-NC-SA 4.0</a>; changes here comprise nearest-neighbor AAL3 label sampling, indexed regional mesh extraction, axis conversion, and rendering. The fsaverage6 dataset page lists the surface data license as unknown; formal redistribution review is pending.</p>
        <ul className="source-list">{metadata.sources.map((source) => <li key={source.name}><a href={source.url} target="_blank" rel="noreferrer">{source.name}</a></li>)}</ul>
      </section>
    </main>
  )
}
