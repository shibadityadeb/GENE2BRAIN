# Stage 03 GWAS Catalog API record

- **Source:** NHGRI-EBI GWAS Catalog
- **Official documentation:** https://www.ebi.ac.uk/gwas/rest/docs/api
- **Access date:** 2026-09-13
- **Response format:** JSON using HAL (`_embedded` and `_links`)
- **Live API root used:** `https://www.ebi.ac.uk/gwas/rest/api`
- **API version finding:** the current official root and documentation expose an unnumbered REST API. Direct checks of `/gwas/rest/api/v2`, `/gwas/api/v2`, and their trailing-slash variants returned HTTP 404 on 2026-09-13. The retired summary-statistics API root returned HTTP 410. Therefore this workflow uses the current live API advertised by the official root; it does not mislabel it as v2.

## Endpoints and parameters

1. `https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findByEfoTrait`
   - `trait=Parkinson disease`, `size=100`
   - resolved uniquely to `MONDO_0005180` (http://purl.obolibrary.org/obo/MONDO_0005180).
2. `https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005180/studies`
   - `size=500`
   - returned 113 ontology-linked candidate studies.
3. `https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByEfoTrait`
   - `efoTrait=Parkinson disease`, `projection=associationByStudy`, `size=5000`
   - returned 814 associations used only to count records per candidate study.
4. `https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90308590/associations`
   - `projection=associationByStudy`, `size=500`
   - returned 121 reported primary-study associations before the fixed p-value filter.
5. `https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90308590/snps`
   - `projection=snp`, `size=500`
   - supplied rsIDs, current Catalog coordinates, functional class, cytogenetic region, and genomic-context gene mappings.
6. `https://www.ebi.ac.uk/gwas/rest/api/metadata`
   - no query parameters
   - supplied the Catalog's current mapping assembly: `GRCh38.p14`.

The ontology relation, rather than a free-text study query, defines the candidate set. A label query is used only to resolve the ontology resource itself.
