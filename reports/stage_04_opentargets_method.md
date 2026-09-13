# Stage 4 Open Targets method

## Current resource and access

- Source: Open Targets Platform, not the retired Open Targets Genetics-only API.
- Access date: 2026-09-13
- GraphQL endpoint: `https://api.platform.opentargets.org/api/v4/graphql`
- Platform API version: 26.6.3
- Platform data release: 26.06 (`platform2606`)
- Documentation: [GraphQL API](https://platform-docs.opentargets.org/data-access/graphql-api), [credible sets](https://platform-docs.opentargets.org/credible-set), and [Locus-to-Gene](https://platform-docs.opentargets.org/gentropy/locus-to-gene-l2g).
- Exact queries, variables, request timestamps, and response sizes: `data/gwas/opentargets_query_log.txt`.

The current GraphQL API was selected because this is a targeted retrieval for one
study. Open Targets recommends its Parquet downloads for systematic bulk analysis;
that route was not necessary for 67 study credible sets. The study query was made
once, followed by a parameterised detail query for each returned credible-set ID.

## Interpretation and thresholds

Open Targets defines a credible set as variants whose posterior probabilities
collectively cover the top 95% likelihood of containing a causal variant. The lead
variant is therefore not assumed causal. For this study, all sets were produced by
PICS from curated reported top hits and carry the corresponding lower-confidence
annotation; no in-sample LD fine mapping or study summary statistics are available
in the Platform.

The current Platform displays protein-coding L2G predictions only when score >
0.05. That documented release filter defines the **broad** eligible set. The
continuous score is retained. The **stringent** operating threshold is the empirical
75th percentile of eligible gene–credible-set scores in this retrieval (0.615029);
it selects genes whose maximum score reaches that boundary. This relative threshold
is reproducible for the pinned release but has no intrinsic biological or causal
meaning. The weighted set retains every broad gene and uses its maximum observed L2G
score as its weight, avoiding inflation from genes appearing at several loci.

Runtime package versions: pandas 2.2.3, NumPy 1.26.4, requests
2.32.3, seaborn 0.13.2, and matplotlib
3.9.4.

L2G is a machine-learning prioritisation score, not proof of causality. Scores can
change between releases as the model and source evidence are updated.
