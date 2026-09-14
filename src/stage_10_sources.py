"""Stage 10 source audit and standardized multi-disease L2G retrieval.

This module freezes candidate/primary GWAS metadata and creates broad,
stringent, and weighted gene sets with the unchanged Stage 4 Open Targets
rules. It deliberately performs no AHBA scoring or regional inference.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timezone
import gzip
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
GWAS = ROOT / "data" / "gwas"
GENES = ROOT / "data" / "genes"
MULTI_GENES = GENES / "multidisease"
RESULTS = ROOT / "data" / "results"
CACHE = ROOT / "data" / "intermediate" / "stage_10" / "opentargets"
REPORTS = ROOT / "reports"

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
ACCESS_DATE = date.today().isoformat()

STUDY_QUERY = """query Stage10Study($id: String!) {
  meta { name product apiVersion { x y z suffix }
    dataVersion { year month iteration } dataPrefix downloads }
  study(studyId: $id) {
    id traitFromSource traitFromSourceMappedIds projectId studyType
    publicationTitle publicationFirstAuthor publicationDate pubmedId publicationJournal
    nSamples nCases nControls initialSampleSize hasSumstats summarystatsLocation
    discoverySamples { ancestry sampleSize }
    replicationSamples { ancestry sampleSize }
    credibleSets(page: {index: 0, size: 500}) { count rows {
      studyLocusId studyId chromosome position region locusStart locusEnd
      finemappingMethod confidence credibleSetlog10BF credibleSetIndex
      purityMeanR2 purityMinR2 pValueMantissa pValueExponent sampleSize
      qualityControls studyType
      variant { id rsIds chromosome position referenceAllele alternateAllele }
    } }
  }
}"""

DETAIL_QUERY = """query Stage10CredibleSet($id: String!) {
  credibleSet(studyLocusId: $id) {
    studyLocusId studyId chromosome position region locusStart locusEnd
    finemappingMethod confidence credibleSetlog10BF credibleSetIndex
    purityMeanR2 purityMinR2 pValueMantissa pValueExponent sampleSize
    qualityControls studyType
    variant { id rsIds chromosome position referenceAllele alternateAllele }
    locus(page: {index: 0, size: 3000}) { count rows {
      is95CredibleSet is99CredibleSet posteriorProbability logBF beta
      standardError r2Overall pValueMantissa pValueExponent
      variant { id rsIds chromosome position referenceAllele alternateAllele }
    } }
    l2GPredictions(page: {index: 0, size: 3000}) { count rows {
      score shapBaseValue target { id approvedSymbol approvedName }
      features { name value shapValue }
    } }
  }
}"""


def load_json_yaml(path: Path) -> dict[str, Any]:
    """Load JSON-compatible YAML without adding a parser dependency."""
    return json.loads(path.read_text(encoding="utf-8"))


def panel() -> list[dict[str, Any]]:
    return load_json_yaml(CONFIG / "disease_panel.yaml")["diseases"]


def parameters() -> dict[str, Any]:
    return load_json_yaml(CONFIG / "statistical_parameters.yaml")


def retry_session(methods: tuple[str, ...]) -> requests.Session:
    client = requests.Session()
    retry = Retry(
        total=6,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=methods,
    )
    client.mount("https://", HTTPAdapter(max_retries=retry))
    client.headers.update(
        {"Accept": "application/json", "User-Agent": "GENE2BRAIN-stage-10/1.0 (research)"}
    )
    return client


def get_json(client: requests.Session, url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    response = client.get(url, params=params, timeout=240)
    response.raise_for_status()
    return response.json()


def embedded(payload: dict[str, Any], key: str) -> list[dict[str, Any]]:
    return payload.get("_embedded", {}).get(key, [])


def paged(client: requests.Session, url: str, key: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    page = 0
    while True:
        payload = get_json(client, url, {"size": 500, "page": page})
        rows.extend(embedded(payload, key))
        metadata = payload.get("page") or {}
        if page + 1 >= int(metadata.get("totalPages", 1)):
            return rows
        page += 1


def sample_total(study: dict[str, Any], sample_type: str = "initial") -> int:
    return int(
        sum(
            group.get("numberOfIndividuals") or 0
            for group in study.get("ancestries") or []
            if group.get("type") == sample_type
        )
    )


def ancestry(study: dict[str, Any]) -> str:
    values: list[str] = []
    for cohort in study.get("ancestries") or []:
        for group in cohort.get("ancestralGroups") or []:
            value = group.get("ancestralGroup")
            if value and value not in values:
                values.append(value)
    return "; ".join(values) if values else "Not reported"


def normalize_trait(value: str) -> str:
    value = value.lower().replace("parkinson's", "parkinson").replace("alzheimer's", "alzheimer")
    return re.sub(r"[^a-z0-9]+", "", value)


def direct_trait(disease: dict[str, Any], reported: str) -> bool:
    expected = disease["disease_name"]
    aliases = {
        "parkinson": ["Parkinson disease", "Parkinson's disease"],
        "alzheimer": ["Alzheimer disease", "Alzheimer's disease"],
        "als": ["Amyotrophic lateral sclerosis"],
        "migraine": ["Migraine"],
        "adhd": ["Attention deficit-hyperactivity disorder", "Attention deficit hyperactivity disorder"],
    }
    allowed = aliases.get(disease["disease_id"], [expected])
    return normalize_trait(reported) in {normalize_trait(value) for value in allowed}


def candidate_eligible(disease: dict[str, Any], study: dict[str, Any]) -> tuple[bool, str]:
    reported = ((study.get("diseaseTrait") or {}).get("trait") or "").strip()
    technology = "; ".join(
        item.get("genotypingTechnology") or ""
        for item in study.get("genotypingTechnologies") or []
    ).lower()
    excluded = (
        "progression", "age at onset", "age of onset", "pleiotropy", "cross-trait",
        "mtag", "medication", "gene-based burden", "cnv", "chromosome x",
    )
    if not direct_trait(disease, reported):
        return False, "phenotype is not the configured direct disease-risk trait"
    if any(term in reported.lower() for term in excluded):
        return False, "progression/proxy/cross-trait/specialized phenotype"
    if "genome-wide" not in technology:
        return False, "not a variant-level genome-wide design"
    return True, "direct disease-risk phenotype and variant-level genome-wide design"


def association_count(client: requests.Session, accession: str) -> int:
    payload = get_json(
        client,
        f"{GWAS_API}/studies/{accession}/associations",
        {"projection": "associationByStudy", "size": 1},
    )
    return int((payload.get("page") or {}).get("totalElements", len(embedded(payload, "associations"))))


def independent_association_count(accession: str) -> int:
    """Fetch one count with an independent session for safe concurrency."""
    with retry_session(("GET",)) as client:
        return association_count(client, accession)


def audit_catalog() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create a bounded ontology-derived candidate snapshot and primary table."""
    client = retry_session(("GET",))
    records: list[dict[str, Any]] = []
    primary_records: list[dict[str, Any]] = []
    for disease in panel():
        studies = paged(
            client,
            f"{GWAS_API}/efoTraits/{disease['ontology_id']}/studies",
            "studies",
        )
        ranked: list[tuple[bool, int, str, dict[str, Any], str]] = []
        for study in studies:
            eligible, reason = candidate_eligible(disease, study)
            publication = study.get("publicationInfo") or {}
            year = (publication.get("publicationDate") or "")[:4]
            ranked.append((eligible, sample_total(study), year, study, reason))
        ranked.sort(key=lambda row: (row[0], row[1], row[2]), reverse=True)
        selected_accession = disease.get("gwas_accession")
        short = ranked[:12]
        if selected_accession and all(row[3].get("accessionId") != selected_accession for row in short):
            selected = next((row for row in ranked if row[3].get("accessionId") == selected_accession), None)
            if selected:
                short.append(selected)
        with ThreadPoolExecutor(max_workers=8) as pool:
            count_futures = {
                study.get("accessionId"): pool.submit(
                    independent_association_count, study.get("accessionId")
                )
                for _, _, _, study, _ in short
            }
            counts = {accession: future.result() for accession, future in count_futures.items()}
        for eligible, _, _, study, reason in short:
            publication = study.get("publicationInfo") or {}
            accession = study.get("accessionId")
            count = counts[accession]
            row = {
                "disease": disease["disease_name"],
                "disease_id": disease["disease_id"],
                "study_accession": accession,
                "trait": (study.get("diseaseTrait") or {}).get("trait"),
                "publication": publication.get("title"),
                "pmid": publication.get("pubmedId"),
                "publication_year": (publication.get("publicationDate") or "")[:4],
                "discovery_sample_size": sample_total(study),
                "replication_sample_size": sample_total(study, "replication"),
                "ancestry": ancestry(study),
                "number_of_associations": count,
                "summary_statistics_available": bool(study.get("fullPvalueSet")),
                "genome_build": "GRCh38.p14 (current GWAS Catalog mapping)",
                "candidate_eligible": eligible,
                "primary_selected": accession == selected_accession,
                "selection_assessment": reason,
                "source": f"{GWAS_API}/studies/{accession}",
                "access_date": ACCESS_DATE,
            }
            records.append(row)
        chosen = next((row for row in records if row["disease_id"] == disease["disease_id"] and row["primary_selected"]), None)
        primary_records.append(
            {
                "disease": disease["disease_name"],
                "disease_id": disease["disease_id"],
                "study_accession": selected_accession,
                "study": disease.get("primary_gwas_study"),
                "status": disease["status"],
                "sample_size": disease.get("sample_size"),
                "ancestry": disease.get("ancestry"),
                "genome_build": disease.get("genome_build"),
                "summary_statistics_available": chosen.get("summary_statistics_available") if chosen else pd.NA,
                "number_of_associations": chosen.get("number_of_associations") if chosen else pd.NA,
                "selection_note": disease["selection_note"],
                "source": chosen.get("source") if chosen else f"{GWAS_API}/efoTraits/{disease['ontology_id']}/studies",
            }
        )
    candidate = pd.DataFrame(records)
    primary = pd.DataFrame(primary_records)
    candidate.to_csv(GWAS / "multidisease_candidate_studies.csv", index=False)
    primary.to_csv(GWAS / "multidisease_primary_gwas.csv", index=False)
    return candidate, primary


def graphql(query: str, variables: dict[str, Any]) -> dict[str, Any]:
    with retry_session(("POST",)) as client:
        response = client.post(OT_API, json={"query": query, "variables": variables}, timeout=300)
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(json.dumps(payload["errors"], indent=2))
    return payload["data"]


def cache_json(path: Path, value: Any) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, separators=(",", ":"))


def get_study_l2g(accession: str, disease_id: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    study_cache = CACHE / f"{disease_id}_study.json.gz"
    details_cache = CACHE / f"{disease_id}_credible_sets.json.gz"
    if study_cache.exists() and details_cache.exists():
        with gzip.open(study_cache, "rt", encoding="utf-8") as stream:
            study = json.load(stream)
        with gzip.open(details_cache, "rt", encoding="utf-8") as stream:
            details = json.load(stream)
        if study.get("id") == accession:
            return study, details
    data = graphql(STUDY_QUERY, {"id": accession})
    study = data.get("study")
    if not study or study.get("id") != accession:
        raise RuntimeError(f"Open Targets study not found: {accession}")
    summaries = (study.get("credibleSets") or {}).get("rows") or []
    count = int((study.get("credibleSets") or {}).get("count") or 0)
    if count != len(summaries):
        raise RuntimeError(f"Credible-set pagination incomplete for {accession}: {len(summaries)}/{count}")
    details: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {
            pool.submit(graphql, DETAIL_QUERY, {"id": row["studyLocusId"]}): row["studyLocusId"]
            for row in summaries
        }
        for future in as_completed(futures):
            identifier = futures[future]
            detail = future.result().get("credibleSet")
            if not detail or detail.get("studyLocusId") != identifier:
                raise RuntimeError(f"Incomplete credible-set response: {identifier}")
            details.append(detail)
    details.sort(key=lambda row: (str(row.get("chromosome")), int(row.get("position") or 0)))
    study["_meta"] = data.get("meta")
    cache_json(CACHE / f"{disease_id}_study.json.gz", study)
    cache_json(CACHE / f"{disease_id}_credible_sets.json.gz", details)
    return study, details


def p_value(mantissa: Any, exponent: Any) -> float | None:
    if mantissa is None or exponent is None:
        return None
    return float(mantissa) * 10.0 ** int(exponent)


def l2g_table(disease: dict[str, Any], details: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for credible_set in details:
        variant = credible_set.get("variant") or {}
        lead_variant = ((variant.get("rsIds") or []) + [variant.get("id")])[0]
        for prediction in (credible_set.get("l2GPredictions") or {}).get("rows") or []:
            target = prediction.get("target") or {}
            rows.append(
                {
                    "disease": disease["disease_name"],
                    "disease_id": disease["disease_id"],
                    "gene": target.get("approvedSymbol"),
                    "ensembl_id": target.get("id"),
                    "locus_id": credible_set.get("region") or credible_set.get("studyLocusId"),
                    "credible_set_id": credible_set.get("studyLocusId"),
                    "lead_variant": lead_variant,
                    "chromosome": credible_set.get("chromosome"),
                    "position": credible_set.get("position"),
                    "p_value_if_available": p_value(credible_set.get("pValueMantissa"), credible_set.get("pValueExponent")),
                    "l2g_score": float(prediction["score"]),
                    "study": disease["gwas_accession"],
                    "fine_mapping_method": credible_set.get("finemappingMethod"),
                    "credible_set_confidence": credible_set.get("confidence"),
                    "evidence_source": "Open Targets Platform Locus-to-Gene",
                    "genome_build": "GRCh38",
                    "score_interpretation": "gene prioritization evidence; not proof of causality",
                }
            )
    if not rows:
        raise RuntimeError(f"No L2G predictions returned for {disease['disease_id']}")
    return pd.DataFrame(rows).sort_values(["l2g_score", "gene"], ascending=[False, True])


def write_gene_sets(disease: dict[str, Any], evidence: pd.DataFrame) -> dict[str, pd.DataFrame]:
    broad_threshold = float(parameters()["l2g_broad_threshold"])
    eligible = evidence.loc[evidence["l2g_score"] > broad_threshold].copy()
    if eligible.empty:
        raise RuntimeError(f"No predictions above the frozen L2G threshold for {disease['disease_id']}")
    threshold = float(eligible["l2g_score"].quantile(0.75))
    genes = (
        eligible.sort_values(["l2g_score", "gene"], ascending=[False, True])
        .groupby(["gene", "ensembl_id"], as_index=False, dropna=False)
        .agg(
            gene_weight=("l2g_score", "max"),
            max_l2g_score=("l2g_score", "max"),
            mean_l2g_score=("l2g_score", "mean"),
            n_credible_sets=("credible_set_id", "nunique"),
            study=("study", "first"),
        )
        .sort_values(["gene_weight", "gene"], ascending=[False, True])
    )
    genes["score_interpretation"] = "Open Targets L2G prioritization score; not proof of causality"
    broad = genes.copy()
    broad["selection_rule"] = f"maximum L2G score > {broad_threshold}"
    stringent = genes.loc[genes["max_l2g_score"] >= threshold].copy()
    stringent["selection_rule"] = f"maximum L2G score >= eligible prediction Q3 ({threshold:.6f})"
    weighted = genes.copy()
    weighted["selection_rule"] = "eligible under broad rule; gene_weight=max L2G score"
    for name, frame in {"broad": broad, "stringent": stringent, "weighted": weighted}.items():
        frame.to_csv(MULTI_GENES / f"{disease['disease_id']}_{name}.csv", index=False)
    return {"broad": broad, "stringent": stringent, "weighted": weighted}


def copy_frozen_parkinson() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    evidence = pd.read_csv(GENES / "parkinson_gene_locus_evidence.csv")
    if "disease_id" not in evidence:
        evidence.insert(1, "disease_id", "parkinson")
    sets: dict[str, pd.DataFrame] = {}
    for name in ("broad", "stringent", "weighted"):
        frame = pd.read_csv(GENES / f"parkinson_genes_{name}.csv")
        frame.to_csv(MULTI_GENES / f"parkinson_{name}.csv", index=False)
        sets[name] = frame
    return evidence, sets


def retrieve_l2g() -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence_frames: list[pd.DataFrame] = []
    metadata: list[dict[str, Any]] = []
    for disease in panel():
        if disease["status"] == "excluded":
            metadata.append(
                {
                    "disease": disease["disease_name"], "study": disease.get("gwas_accession"),
                    "trait": disease["disease_name"], "sample_size": disease.get("sample_size"),
                    "replication_sample_size": pd.NA, "ancestry": disease.get("ancestry"),
                    "genome_build": disease.get("genome_build"), "number_of_significant_associations": pd.NA,
                    "number_of_loci": 0, "summary_statistics_available": pd.NA,
                    "publication_year": pd.NA, "source": GWAS_API, "status": disease["status"],
                }
            )
            continue
        if disease["disease_id"] == "parkinson":
            evidence, _ = copy_frozen_parkinson()
            study_data = graphql(STUDY_QUERY, {"id": disease["gwas_accession"]})
            study = study_data["study"]
        else:
            study, details = get_study_l2g(disease["gwas_accession"], disease["disease_id"])
            evidence = l2g_table(disease, details)
            write_gene_sets(disease, evidence)
        evidence_frames.append(evidence)
        metadata.append(
            {
                "disease": disease["disease_name"],
                "study": disease["gwas_accession"],
                "trait": study.get("traitFromSource"),
                "sample_size": study.get("nSamples"),
                "replication_sample_size": sum(x.get("sampleSize") or 0 for x in study.get("replicationSamples") or []),
                "ancestry": "; ".join(sorted({x.get("ancestry") for x in study.get("discoverySamples") or [] if x.get("ancestry")})),
                "genome_build": "GRCh38",
                "number_of_significant_associations": pd.NA,
                "number_of_loci": int((study.get("credibleSets") or {}).get("count") or 0),
                "summary_statistics_available": bool(study.get("hasSumstats")),
                "publication_year": (study.get("publicationDate") or "")[:4],
                "source": f"{OT_API} study={disease['gwas_accession']}",
                "status": disease["status"],
            }
        )
    combined = pd.concat(evidence_frames, ignore_index=True, sort=False)
    combined.to_csv(GENES / "multidisease_gene_locus_evidence.csv", index=False)
    gwas_metadata = pd.DataFrame(metadata)
    candidate_path = GWAS / "multidisease_candidate_studies.csv"
    if candidate_path.exists():
        candidates = pd.read_csv(candidate_path)
        counts = candidates.loc[candidates["primary_selected"].astype(bool)].set_index("disease")["number_of_associations"]
        gwas_metadata["number_of_significant_associations"] = gwas_metadata["disease"].map(counts)
    gwas_metadata.to_csv(RESULTS / "multidisease_gwas_metadata.csv", index=False)
    return combined, gwas_metadata


def write_source_report(evidence: pd.DataFrame, metadata: pd.DataFrame) -> None:
    ready = [item for item in panel() if item["status"] == "ready"]
    needs_review = [item for item in panel() if item["status"] == "needs_review"]
    excluded = [item for item in panel() if item["status"] == "excluded"]
    version = "unknown"
    cached = next(CACHE.glob("*_study.json.gz"), None)
    if cached:
        with gzip.open(cached, "rt", encoding="utf-8") as stream:
            study = json.load(stream)
        meta = study.get("_meta") or {}
        version = json.dumps({"api": meta.get("apiVersion"), "data": meta.get("dataVersion")}, sort_keys=True)
    report = f"""# Stage 10 source and gene-prioritization audit

- Access date: {ACCESS_DATE}
- GWAS source: NHGRI-EBI GWAS Catalog REST API
- Gene prioritization: Open Targets Platform GraphQL API
- Open Targets version metadata: `{version}`
- Ready diseases: {len(ready)}
- Needs-review diseases: {len(needs_review)}
- Excluded diseases: {len(excluded)}
- Combined L2G evidence rows: {len(evidence):,}
- Unique disease-gene pairs: {evidence.groupby(['disease', 'gene']).ngroups:,}

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

Generated at {datetime.now(timezone.utc).isoformat()}.
"""
    (REPORTS / "stage_10_source_audit.md").write_text(report, encoding="utf-8")


def validate_outputs() -> None:
    required = [
        GWAS / "multidisease_candidate_studies.csv",
        GWAS / "multidisease_primary_gwas.csv",
        RESULTS / "multidisease_gwas_metadata.csv",
        GENES / "multidisease_gene_locus_evidence.csv",
    ]
    required.extend(
        MULTI_GENES / f"{disease['disease_id']}_{kind}.csv"
        for disease in panel() if disease["status"] in {"ready", "needs_review"}
        for kind in ("broad", "stringent", "weighted")
    )
    missing = [path for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError("Missing Stage 10 source outputs: " + ", ".join(map(str, missing)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-catalog", action="store_true", help="reuse existing candidate/primary tables")
    args = parser.parse_args()
    for directory in (GWAS, MULTI_GENES, RESULTS, CACHE, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    if not args.skip_catalog:
        audit_catalog()
    evidence, metadata = retrieve_l2g()
    write_source_report(evidence, metadata)
    validate_outputs()
    print("GENE2BRAIN STAGE 10 SOURCE SUMMARY")
    print(f"Panel diseases: {len(panel())}")
    print(f"Ready diseases: {sum(x['status'] == 'ready' for x in panel())}")
    print(f"Needs-review diseases: {sum(x['status'] == 'needs_review' for x in panel())}")
    print(f"Excluded diseases: {sum(x['status'] == 'excluded' for x in panel())}")
    print(f"L2G evidence rows: {len(evidence):,}")
    print("STOP: no AHBA scoring or regional inference performed by this module.")


if __name__ == "__main__":
    main()
