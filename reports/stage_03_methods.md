# Stage 03 methods: Parkinson disease GWAS

Generated/accessed: 2026-09-13

## Source and search

The source was the current NHGRI-EBI GWAS Catalog HAL REST API at `https://www.ebi.ac.uk/gwas/rest/api`. The live version audit is recorded in `stage_03_gwas_api.md`. The exact Parkinson disease ontology resource was `MONDO_0005180` (http://purl.obolibrary.org/obo/MONDO_0005180); all 113 studies linked from that resource were retained as candidates. This avoids defining the candidate set by publication-title or trait-string matching.

## Reproducible study selection

Candidates were first classified as primary-eligible only when the reported trait was direct Parkinson disease, the design was variant-level genome-wide genotyping/sequencing, at least one ontology-linked association was present, and the trait did not represent progression, severity, subtype, familial disease, proxy-only disease, medication, interaction, pleiotropy, CNV, or gene-burden analysis.

Eligible studies received an ordinal score: discovery sample size (0–4 points), reported associations (1–3), full p-value set (0 or 2), publication from 2020 onward (0 or 1), ancestry breadth (0–2), and a replication stage (0 or 1). Ties were ordered by discovery sample size, association count, then accession. All component data and reasons are preserved in `parkinson_study_selection.csv`.

The selected primary study is **GCST90308590**, rank 1: Multi-ancestry genome-wide association meta-analysis of Parkinson's disease. (Kim JJ, 2023, PMID 38155330). It has a direct Parkinson phenotype, 2,525,730 discovery/meta-analysis participants across European; East Asian; Hispanic or Latin American; African unspecified, and the largest curated reported-association set (121) among candidates. Although the Catalog does not flag a downloadable full p-value set for this accession, its power, phenotype, ancestry breadth, recency, and curated associations outweighed that limitation.

Candidate sensitivity/replication datasets, not analyzed here, are `GCST90480008` (a publication-distinct multi-ancestry biobank analysis), `GCST90270939` (Chinese whole-genome sequencing GWAS), and `GCST90828116` (a later Taiwanese GWAS). These are publication/cohort-distinct from the component datasets listed for the open Kim companion meta-analysis. Participant-level overlap cannot be ruled out from Catalog metadata alone and must be audited before treating any result as statistically independent replication.

## Associations and genome assembly

The standard threshold **p ≤ 5×10⁻⁸** was applied without modification. Of 121 Catalog records for GCST90308590, 109 met the threshold and 12 were excluded. The output is a significant-association table, not a complete summary-statistics dataset. Variant locations are current GWAS Catalog SNP-resource coordinates on **GRCh38.p14**, retrieved from the API metadata endpoint; the output labels this explicitly. No coordinate liftover was performed.

## Initial locus grouping

Significant variants were sorted by chromosome and position. Within each chromosome, adjacent variants separated by **≤ 250,000 bp** were merged using single-link distance grouping; a gap > 250,000 bp began a new locus. Locus bounds are the minimum and maximum reported significant-variant positions, and the lowest-p variant is the lead. This FUMA-inspired merge distance is a transparent initial consolidation, yielding 78 loci from 109 associations. It is **not LD clumping**, because an ancestry-matched LD reference was not introduced in Stage 3, and it does not imply signal independence or causality. LD-aware refinement belongs in a later genetics stage.

## Mapped candidate genes

The selected study has no author-reported genes in its association objects. Therefore `parkinson_gwas_mapped_genes.csv` preserves genes from the GWAS Catalog SNP `genomicContexts` field, with Entrez IDs, Ensembl/NCBI source, mapping method, and positional relationship where available. These are explicitly **GWAS-mapped candidate genes, not causal genes or confirmed risk genes**. No functional prioritization was performed.

## Summary statistics and QC

The Catalog `fullPvalueSet` field is false for GCST90308590, but the publication's data-availability statement identifies an immediately accessible open companion meta-analysis under **GCST90275127**. Official FTP metadata describes it as a 1.2 GB, unharmonized GRCh37 GWAS-SSF v1.0 TSV excluding 23andMe; the full meta-analysis is available to qualified researchers under a 23andMe agreement. The exact open URL and scope are recorded, but the large file was not downloaded. Validation candidate GCST90480008 has a separate verified 827 MB, unharmonized GRCh38 GWAS-SSF v1.0 compressed TSV; its URL is recorded without downloading. The other candidate validation accessions are not flagged as full p-value sets.

QC tested duplicate variants, rsID/coordinate/chromosome validity, missing/invalid p-values, duplicate study accessions, genome-build labels, mapped-gene coverage, allele completeness, and effect-size availability/range. Findings are reported without silent repairs in `stage_03_gwas_qc.md`.

No Open Targets analysis, causal-gene inference, AHBA integration, regional enrichment, permutations, or brain mapping was performed.
