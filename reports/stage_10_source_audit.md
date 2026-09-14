# Stage 10 source and gene-prioritization audit

- Access date: 2026-09-14
- GWAS source: NHGRI-EBI GWAS Catalog REST API
- Gene prioritization: Open Targets Platform GraphQL API
- Open Targets version metadata: `{"api": {"suffix": null, "x": "26", "y": "6", "z": "3"}, "data": {"iteration": null, "month": "06", "year": "26"}}`
- Ready diseases: 10
- Needs-review diseases: 2
- Excluded diseases: 3
- Combined L2G evidence rows: 1,621
- Unique disease-gene pairs: 1,566

The source stage applies the frozen Stage 4 broad rule (`L2G > 0.05`), the same
upper-quartile stringent rule, and maximum-L2G gene weights. Parkinson source
files are copied into the multi-disease namespace without modifying the frozen
Stage 4 files. Open Targets responses are cached as gzip-compressed JSON under
`data/intermediate/stage_10/opentargets/`.

Candidate selection is ontology-derived and bounded to the twelve highest-sample
candidate records per disease plus the configured primary when necessary. The
candidate table records phenotype eligibility and the current Catalog association
count. The three excluded diseases remain in metadata and configuration; no
positional-gene or literature-gene substitution was made.

Generated at 2026-09-14T15:25:33.498569+00:00.
