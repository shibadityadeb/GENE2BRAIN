"""Stage 4: Open Targets credible sets and Locus-to-Gene prioritisation for PD.

This stage converts the Stage 3 Parkinson GWAS evidence into broad, stringent,
and weighted candidate-gene sets. It deliberately does not calculate regional
AHBA scores, enrichment, permutations, pathways, or disease brain maps.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
GENES = ROOT / "data" / "genes"
RESULT_DATA = ROOT / "data" / "results"
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"

STUDY_ID = "GCST90308590"
ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
GRAPHQL_DOCS = "https://platform-docs.opentargets.org/data-access/graphql-api"
L2G_DOCS = "https://platform-docs.opentargets.org/gentropy/locus-to-gene-l2g"
CREDIBLE_SET_DOCS = "https://platform-docs.opentargets.org/credible-set"
BROAD_THRESHOLD = 0.05
ACCESS_DATE = date.today().isoformat()

STUDY_QUERY = """query Stage4Study($id: String!) {
  meta { name product apiVersion { x y z suffix }
    dataVersion { year month iteration } dataPrefix downloads }
  study(studyId: $id) {
    id traitFromSource traitFromSourceMappedIds projectId studyType
    publicationTitle publicationFirstAuthor publicationDate pubmedId publicationJournal
    nSamples nCases nControls initialSampleSize hasSumstats summarystatsLocation
    discoverySamples { ancestry sampleSize }
    replicationSamples { ancestry sampleSize }
    ldPopulationStructure { ldPopulation relativeSampleSize }
    qualityControls analysisFlags cohorts diseases { id name }
    credibleSets(page: {index: 0, size: 500}) { count rows {
      studyLocusId studyId chromosome position region locusStart locusEnd
      finemappingMethod confidence credibleSetlog10BF credibleSetIndex
      purityMeanR2 purityMinR2 pValueMantissa pValueExponent sampleSize
      qualityControls studyType
      variant { id rsIds chromosome position referenceAllele alternateAllele }
    } }
  }
}"""

DETAIL_QUERY = """query Stage4CredibleSet($id: String!) {
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


def make_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=6,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("POST",),
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.headers.update(
        {"User-Agent": "GENE2BRAIN-stage-04/1.0 (research)", "Accept": "application/json"}
    )
    return session


def graphql(query: str, variables: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    started = datetime.now(timezone.utc)
    with make_session() as session:
        response = session.post(
            ENDPOINT,
            json={"query": query, "variables": variables},
            timeout=240,
        )
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(json.dumps(payload["errors"], indent=2))
    metadata = {
        "timestamp_utc": started.isoformat(),
        "variables": variables,
        "http_status": response.status_code,
        "elapsed_seconds": round(response.elapsed.total_seconds(), 3),
        "response_bytes": len(response.content),
    }
    return payload["data"], metadata


def p_value(mantissa: Any, exponent: Any) -> float | None:
    if mantissa is None or exponent is None:
        return None
    return float(mantissa) * (10.0 ** int(exponent))


def json_compact(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def audit_stage_03() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assoc = pd.read_csv(GWAS / "parkinson_gwas_associations_raw.csv")
    loci = pd.read_csv(GWAS / "parkinson_initial_loci.csv")
    mapped = pd.read_csv(GWAS / "parkinson_gwas_mapped_genes.csv")
    selection = pd.read_csv(GWAS / "parkinson_study_selection.csv")
    selected = selection.loc[selection["primary_selected"]]
    errors: list[str] = []
    if len(selected) != 1 or selected.iloc[0]["study_accession"] != STUDY_ID:
        errors.append("Stage 3 primary study is not uniquely GCST90308590")
    if set(assoc["study_accession"]) != {STUDY_ID}:
        errors.append("Stage 3 association study accession is inconsistent")
    if set(assoc["genomic_assembly"]) != {"GRCh38.p14"}:
        errors.append("Stage 3 association build is not GRCh38.p14")
    if set(assoc["locus_id"]) != set(loci["locus_id"]):
        errors.append("Association and provisional-locus identifiers differ")
    if int(loci["number_of_associations"].sum()) != len(assoc):
        errors.append("Provisional-locus association counts do not sum to associations")
    for _, locus in loci.iterrows():
        members = assoc.loc[assoc["locus_id"] == locus["locus_id"]]
        lead = members.loc[members["p_value"].idxmin()]
        if lead["rsid"] != locus["lead_variant"] or not np.isclose(
            float(lead["p_value"]), float(locus["lead_p_value"]), rtol=1e-12
        ):
            errors.append(f"Lead variant mismatch in {locus['locus_id']}")
    if not set(mapped["locus_id"]).issubset(set(loci["locus_id"])):
        errors.append("Mapped-gene table contains an unknown locus")
    if errors:
        raise RuntimeError("Stage 3 consistency audit failed:\n- " + "\n- ".join(errors))
    return assoc, loci, mapped


def get_open_targets() -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    data, study_log = graphql(STUDY_QUERY, {"id": STUDY_ID})
    study = data.get("study")
    if not study or study.get("id") != STUDY_ID:
        raise RuntimeError(f"Open Targets study {STUDY_ID} was not found")
    summaries = study["credibleSets"]["rows"]
    if study["credibleSets"]["count"] != len(summaries):
        raise RuntimeError("Credible-set pagination was incomplete")

    calls: list[dict[str, Any]] = [dict(study_log, query_name="Stage4Study")]
    details: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        jobs = {
            executor.submit(graphql, DETAIL_QUERY, {"id": row["studyLocusId"]}): row["studyLocusId"]
            for row in summaries
        }
        for job in as_completed(jobs):
            csid = jobs[job]
            payload, call = job.result()
            detail = payload.get("credibleSet")
            if not detail or detail.get("studyLocusId") != csid:
                raise RuntimeError(f"Incomplete detail response for credible set {csid}")
            details.append(detail)
            calls.append(dict(call, query_name="Stage4CredibleSet"))
    details.sort(key=lambda x: (int(x["chromosome"]), int(x["position"]), x["studyLocusId"]))
    study["_meta"] = data["meta"]
    return study, details, calls


def stage3_locus_for_variant(
    variant: dict[str, Any], assoc: pd.DataFrame
) -> tuple[str | None, str | None]:
    rsids = variant.get("rsIds") or []
    match = assoc.loc[assoc["rsid"].isin(rsids)] if rsids else assoc.iloc[0:0]
    if match.empty:
        match = assoc.loc[
            (assoc["chromosome"].astype(str) == str(variant.get("chromosome")))
            & (assoc["position"] == variant.get("position"))
        ]
    if match.empty:
        return None, None
    return str(match.iloc[0]["locus_id"]), str(match.iloc[0]["rsid"])


def build_study_match(study: dict[str, Any], details: list[dict[str, Any]], assoc: pd.DataFrame) -> pd.DataFrame:
    ot_rsids = {
        rsid for cs in details for rsid in (cs.get("variant") or {}).get("rsIds", [])
    }
    shared = sorted(ot_rsids & set(assoc["rsid"]))
    discovery = study.get("discoverySamples") or []
    ancestries = sorted({x["ancestry"] for x in discovery if x.get("ancestry")})
    evidence = (
        f"Exact accession; trait={study.get('traitFromSource')}; "
        f"mapped disease IDs={';'.join(study.get('traitFromSourceMappedIds') or [])}; "
        f"sample size={study.get('nSamples')}; {len(shared)}/{len(details)} Open Targets "
        "credible-set lead rsIDs occur among the 109 Stage 3 significant variants; "
        "Open Targets variant coordinates are GRCh38 and Stage 3 coordinates are "
        "GRCh38.p14 (same major assembly, patch label differs)."
    )
    return pd.DataFrame(
        [{
            "gwas_catalog_study": STUDY_ID,
            "opentargets_study": study["id"],
            "trait": study.get("traitFromSource"),
            "genome_build": "GRCh38 (Open Targets); GRCh38.p14 (Stage 3)",
            "ancestry": "; ".join(ancestries),
            "match_status": "matched",
            "matching_evidence": evidence,
        }]
    )


def build_credible_sets(details: list[dict[str, Any]], assoc: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for cs in details:
        variant = cs["variant"]
        locus_rows = cs["locus"]["rows"]
        stage3_locus, lead_rsid = stage3_locus_for_variant(variant, assoc)
        lead = next(
            (x for x in locus_rows if (x.get("variant") or {}).get("id") == variant.get("id")),
            None,
        )
        records.append(
            {
                "credible_set_id": cs["studyLocusId"],
                "study": cs["studyId"],
                "lead_variant": lead_rsid or ((variant.get("rsIds") or [variant["id"]])[0]),
                "lead_variant_ot_id": variant["id"],
                "chromosome": str(cs["chromosome"]),
                "position": int(cs["position"]),
                "number_of_variants": int(cs["locus"]["count"]),
                "fine_mapping_method": cs.get("finemappingMethod"),
                "lead_posterior_inclusion_probability": (lead or {}).get("posteriorProbability"),
                "max_posterior_inclusion_probability": max(
                    (float(x["posteriorProbability"]) for x in locus_rows), default=np.nan
                ),
                "credible_set_confidence": cs.get("confidence"),
                "credible_set_index": cs.get("credibleSetIndex"),
                "credible_set_log10_bf": cs.get("credibleSetlog10BF"),
                "purity_min_r2": cs.get("purityMinR2"),
                "purity_mean_r2": cs.get("purityMeanR2"),
                "p_value": p_value(cs.get("pValueMantissa"), cs.get("pValueExponent")),
                "locus_start": cs.get("locusStart"),
                "locus_end": cs.get("locusEnd"),
                "region": cs.get("region"),
                "stage_03_locus_id": stage3_locus,
                "genome_build": "GRCh38",
                "quality_controls": "; ".join(cs.get("qualityControls") or []),
                "variant_level_information_json": json_compact(locus_rows),
            }
        )
    frame = pd.DataFrame(records)
    frame["_chromosome_order"] = pd.to_numeric(frame["chromosome"])
    return frame.sort_values(["_chromosome_order", "position"]).drop(columns="_chromosome_order")


def grouped_features(features: list[dict[str, Any]], category: str) -> dict[str, float | None]:
    def belongs(name: str) -> bool:
        lower = name.lower()
        if category == "distance":
            return "distance" in lower
        if category == "colocalisation":
            return "coloc" in lower
        if category == "variant_effect":
            return lower.startswith("vep")
        return not ("distance" in lower or "coloc" in lower or lower.startswith("vep"))
    return {x["name"]: x.get("value") for x in features if belongs(x["name"])}


def build_l2g(details: list[dict[str, Any]], credible: pd.DataFrame) -> pd.DataFrame:
    cs_lookup = credible.set_index("credible_set_id").to_dict("index")
    rows: list[dict[str, Any]] = []
    for cs in details:
        info = cs_lookup[cs["studyLocusId"]]
        for pred in cs["l2GPredictions"]["rows"]:
            target = pred.get("target") or {}
            features = pred.get("features") or []
            values = {x["name"]: x.get("value") for x in features}
            shap = {x["name"]: x.get("shapValue") for x in features}
            rows.append(
                {
                    "credible_set_id": cs["studyLocusId"],
                    "gene": target.get("approvedSymbol"),
                    "ensembl_id": target.get("id"),
                    "gene_name": target.get("approvedName"),
                    "l2g_score": float(pred["score"]),
                    "lead_variant": info["lead_variant"],
                    "lead_variant_ot_id": info["lead_variant_ot_id"],
                    "locus": info["stage_03_locus_id"],
                    "chromosome": info["chromosome"],
                    "position": info["position"],
                    "p_value_if_available": info["p_value"],
                    "distance_features_json": json_compact(grouped_features(features, "distance")),
                    "colocalisation_features_json": json_compact(grouped_features(features, "colocalisation")),
                    "variant_effect_features_json": json_compact(grouped_features(features, "variant_effect")),
                    "other_evidence_json": json_compact(grouped_features(features, "other")),
                    "all_feature_values_json": json_compact(values),
                    "all_feature_shap_values_json": json_compact(shap),
                    "shap_base_value": pred.get("shapBaseValue"),
                    "fine_mapping_method": info["fine_mapping_method"],
                    "credible_set_confidence": info["credible_set_confidence"],
                    "study": cs["studyId"],
                    "genome_build": "GRCh38",
                }
            )
    if not rows:
        raise RuntimeError("Open Targets returned no L2G predictions")
    frame = pd.DataFrame(rows).sort_values("l2g_score", ascending=False)
    stringent = float(frame.loc[frame["l2g_score"] > BROAD_THRESHOLD, "l2g_score"].quantile(0.75))
    frame["evidence_class"] = np.select(
        [frame["l2g_score"] >= stringent, frame["l2g_score"] > BROAD_THRESHOLD],
        ["high_prioritization", "moderate_prioritization"],
        default="low_prioritization",
    )
    frame["broad_threshold"] = BROAD_THRESHOLD
    frame["stringent_threshold"] = stringent
    return frame.reset_index(drop=True)


def aggregate_genes(l2g: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (gene, ensembl), group in l2g.loc[l2g["l2g_score"] > BROAD_THRESHOLD].groupby(
        ["gene", "ensembl_id"], dropna=False
    ):
        top = group.sort_values("l2g_score", ascending=False).iloc[0]
        records.append(
            {
                "gene": gene,
                "ensembl_id": ensembl,
                "gene_name": top["gene_name"],
                "gene_weight": float(group["l2g_score"].max()),
                "max_l2g_score": float(group["l2g_score"].max()),
                "mean_l2g_score": float(group["l2g_score"].mean()),
                "n_credible_sets": int(group["credible_set_id"].nunique()),
                "n_stage_03_loci": int(group["locus"].dropna().nunique()),
                "credible_set_ids": ";".join(sorted(group["credible_set_id"].unique())),
                "stage_03_loci": ";".join(sorted(group["locus"].dropna().astype(str).unique())),
                "lead_variants": ";".join(sorted(group["lead_variant"].dropna().astype(str).unique())),
                "top_credible_set_id": top["credible_set_id"],
                "top_lead_variant": top["lead_variant"],
                "top_locus": top["locus"],
                "study": STUDY_ID,
                "score_interpretation": "Open Targets L2G prioritisation score; not proof of causality",
            }
        )
    return pd.DataFrame(records).sort_values(["gene_weight", "gene"], ascending=[False, True])


def comparison_table(mapped: pd.DataFrame, broad: pd.DataFrame, l2g: pd.DataFrame) -> pd.DataFrame:
    gwas = set(mapped["gene"].dropna().astype(str))
    prioritised = set(broad["gene"].dropna().astype(str))
    both = gwas & prioritised
    union = gwas | prioritised
    multi = int((l2g.loc[l2g["l2g_score"] > BROAD_THRESHOLD].groupby("credible_set_id")["gene"].nunique() > 1).sum())
    rows = [
        ("gwas_mapped_unique_genes", len(gwas), len(gwas), 100.0, gwas),
        ("l2g_broad_unique_genes", len(prioritised), len(prioritised), 100.0, prioritised),
        ("genes_only_in_gwas_mapped", len(gwas - prioritised), len(union), 100 * len(gwas - prioritised) / len(union), gwas - prioritised),
        ("genes_only_in_l2g_broad", len(prioritised - gwas), len(union), 100 * len(prioritised - gwas) / len(union), prioritised - gwas),
        ("genes_in_both", len(both), len(union), 100 * len(both) / len(union), both),
        ("overlap_jaccard_percentage", len(both), len(union), 100 * len(both) / len(union), both),
        ("credible_sets_with_multiple_l2g_candidates", multi, l2g["credible_set_id"].nunique(), 100 * multi / l2g["credible_set_id"].nunique(), set()),
    ]
    return pd.DataFrame(
        [{"metric": m, "value": v, "denominator": d, "percentage": p, "genes": ";".join(sorted(gs))} for m, v, d, p, gs in rows]
    )


def master_evidence_table(l2g: pd.DataFrame) -> pd.DataFrame:
    """Retain one row per credible-set/gene prediction with required provenance."""
    master = l2g.copy()
    master.insert(0, "disease", "Parkinson disease")
    master = master.rename(columns={"locus": "locus_id"})
    master["evidence_source"] = "Open Targets Platform Locus-to-Gene"
    required_first = [
        "disease", "gene", "ensembl_id", "locus_id", "credible_set_id",
        "lead_variant", "chromosome", "position", "p_value_if_available",
        "l2g_score", "fine_mapping_method", "credible_set_confidence", "study",
        "evidence_source",
    ]
    return master[required_first + [c for c in master.columns if c not in required_first]]


def top_gene_table(l2g: pd.DataFrame, limit: int = 50) -> pd.DataFrame:
    """Return genes ranked by their highest-scoring credible-set prediction."""
    top = (
        l2g.sort_values("l2g_score", ascending=False)
        .drop_duplicates("gene")
        .head(limit)
        .copy()
    )
    top.insert(0, "rank", np.arange(1, len(top) + 1))
    top["evidence"] = (
        "highest-prioritized candidate gene; Open Targets L2G; "
        + top["evidence_class"].astype(str)
        + "; score is not proof of causality"
    )
    top = top.rename(
        columns={"locus": "locus", "credible_set_id": "credible_set"}
    )
    return top[
        ["rank", "gene", "ensembl_id", "locus", "credible_set", "lead_variant",
         "l2g_score", "fine_mapping_method", "evidence"]
    ]


def ahba_coverage(broad: pd.DataFrame, stringent: pd.DataFrame, weighted: pd.DataFrame) -> pd.DataFrame:
    matrix_path = ROOT / "data" / "processed" / "brain_region_gene_expression.csv"
    available = set(pd.read_csv(matrix_path, nrows=0).columns) - {"region_id"}
    rows = []
    for name, frame in (("broad", broad), ("stringent", stringent), ("weighted", weighted)):
        genes = set(frame["gene"].dropna().astype(str))
        present, absent = genes & available, genes - available
        rows.append(
            {
                "gene_set": name,
                "n_genes": len(genes),
                "n_present_in_ahba": len(present),
                "n_absent_from_ahba": len(absent),
                "coverage_percentage": 100 * len(present) / len(genes) if genes else np.nan,
                "genes_absent_from_ahba": ";".join(sorted(absent)),
                "coverage_only_no_expression_scoring": True,
            }
        )
    return pd.DataFrame(rows)


def make_figures(l2g: pd.DataFrame, comparison: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    # Show all candidates in the 24 loci with the strongest maximum prediction,
    # rather than only the top gene, so multi-candidate and weak tails remain visible.
    order = list(
        l2g.groupby("lead_variant")["l2g_score"].max().nlargest(24).index
    )
    plot = l2g.loc[l2g["lead_variant"].isin(order)].copy()
    gene_order = list(
        plot.assign(_x=plot["lead_variant"].map({x: i for i, x in enumerate(order)}))
        .sort_values(["_x", "l2g_score"], ascending=[False, True])["gene"]
        .drop_duplicates()
    )
    xmap, ymap = {x: i for i, x in enumerate(order)}, {x: i for i, x in enumerate(gene_order)}
    fig, ax = plt.subplots(figsize=(15, max(11, len(gene_order) * 0.28)))
    scatter = ax.scatter(
        plot["lead_variant"].map(xmap), plot["gene"].map(ymap),
        s=30 + 320 * plot["l2g_score"], c=plot["l2g_score"],
        cmap="viridis", vmin=0, vmax=1, edgecolor="white", linewidth=0.5,
    )
    ax.set_xticks(range(len(order)), order, rotation=75, ha="right", fontsize=8)
    ax.set_yticks(range(len(gene_order)), gene_order, fontsize=9)
    ax.set_xlabel("Credible-set lead variant")
    ax.set_ylabel("Top prioritised genes")
    ax.set_title("Parkinson GWAS-to-Gene Evidence Landscape")
    fig.colorbar(scatter, ax=ax, label="L2G score", shrink=0.75)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_04_parkinson_gwas_to_gene.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    strict = float(l2g["stringent_threshold"].iloc[0])
    fig, ax = plt.subplots(figsize=(11, 7))
    sns.histplot(l2g["l2g_score"], bins=25, kde=True, color="#35618f", ax=ax)
    ax.axvline(BROAD_THRESHOLD, color="#d95f02", linestyle="--", label=f"Broad > {BROAD_THRESHOLD:.2f}")
    ax.axvline(strict, color="#1b9e77", linestyle="--", label=f"Stringent ≥ Q3 = {strict:.3f}")
    ax.set(xlabel="Open Targets L2G score", ylabel="Gene–credible-set predictions", title="Distribution of Parkinson L2G scores")
    summary = f"n = {len(l2g)}\nmedian = {l2g['l2g_score'].median():.3f}"
    ax.text(0.04, 0.96, summary, transform=ax.transAxes, ha="left", va="top",
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9})
    strongest = l2g.nlargest(5, "l2g_score")
    ymax = ax.get_ylim()[1]
    for index, (_, row) in enumerate(strongest.iterrows()):
        ax.annotate(
            f"{row['gene']} ({row['l2g_score']:.2f})",
            xy=(row["l2g_score"], ymax * (0.68 - index * 0.085)),
            xytext=(4, 0), textcoords="offset points", fontsize=9,
        )
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_04_parkinson_l2g_distribution.png", dpi=200)
    plt.close(fig)

    vals = comparison.set_index("metric")["value"]
    labels = ["GWAS mapped only", "Both", "L2G broad only"]
    counts = [vals["genes_only_in_gwas_mapped"], vals["genes_in_both"], vals["genes_only_in_l2g_broad"]]
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(labels, counts, color=["#8da0cb", "#66c2a5", "#fc8d62"])
    ax.bar_label(bars, padding=5)
    ax.set(xlabel="Number of unique genes", title="Stage 3 mapped genes versus Open Targets L2G genes")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_04_parkinson_gene_mapping_comparison.png", dpi=200)
    plt.close(fig)


def write_query_log(meta: dict[str, Any], calls: list[dict[str, Any]], details: list[dict[str, Any]]) -> None:
    api = meta["apiVersion"]
    dv = meta["dataVersion"]
    lines = [
        "GENE2BRAIN Stage 4 Open Targets query log",
        f"Access date: {ACCESS_DATE}",
        f"Endpoint: {ENDPOINT}",
        f"API version: {api['x']}.{api['y']}.{api['z']}{api.get('suffix') or ''}",
        f"Data release: {dv['year']}.{int(dv['month']):02d}" + (f".{dv['iteration']}" if dv.get("iteration") else ""),
        f"Data prefix: {meta.get('dataPrefix')}",
        f"Study: {STUDY_ID}",
        f"Credible-set detail responses: {len(details)}",
        "",
        "STUDY QUERY",
        STUDY_QUERY,
        "Variables: " + json_compact({"id": STUDY_ID}),
        "",
        "CREDIBLE-SET DETAIL QUERY TEMPLATE",
        DETAIL_QUERY,
        "Credible-set IDs: " + ";".join(x["studyLocusId"] for x in details),
        "",
        "RESPONSE METADATA (one JSON object per request)",
    ]
    lines.extend(json_compact(x) for x in calls)
    (GWAS / "opentargets_query_log.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_reports(study: dict[str, Any], credible: pd.DataFrame, l2g: pd.DataFrame, broad: pd.DataFrame, stringent: pd.DataFrame, coverage: pd.DataFrame) -> None:
    api = study["_meta"]["apiVersion"]
    dv = study["_meta"]["dataVersion"]
    api_version = f"{api['x']}.{api['y']}.{api['z']}{api.get('suffix') or ''}"
    release = f"{dv['year']}.{int(dv['month']):02d}" + (f".{dv['iteration']}" if dv.get("iteration") else "")
    strict = float(l2g["stringent_threshold"].iloc[0])
    coverage_lines = [
        "| Gene set | Genes | In AHBA | Absent | Coverage |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, row in coverage.iterrows():
        coverage_lines.append(
            f"| {row['gene_set']} | {int(row['n_genes'])} | "
            f"{int(row['n_present_in_ahba'])} | {int(row['n_absent_from_ahba'])} | "
            f"{float(row['coverage_percentage']):.1f}% |"
        )
    coverage_markdown = "\n".join(coverage_lines)
    methods = f"""# Stage 4 Open Targets method

## Current resource and access

- Source: Open Targets Platform, not the retired Open Targets Genetics-only API.
- Access date: {ACCESS_DATE}
- GraphQL endpoint: `{ENDPOINT}`
- Platform API version: {api_version}
- Platform data release: {release} (`{study['_meta'].get('dataPrefix')}`)
- Documentation: [GraphQL API]({GRAPHQL_DOCS}), [credible sets]({CREDIBLE_SET_DOCS}), and [Locus-to-Gene]({L2G_DOCS}).
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
75th percentile of eligible gene–credible-set scores in this retrieval ({strict:.6f});
it selects genes whose maximum score reaches that boundary. This relative threshold
is reproducible for the pinned release but has no intrinsic biological or causal
meaning. The weighted set retains every broad gene and uses its maximum observed L2G
score as its weight, avoiding inflation from genes appearing at several loci.

Runtime package versions: pandas {pd.__version__}, NumPy {np.__version__}, requests
{requests.__version__}, seaborn {sns.__version__}, and matplotlib
{plt.matplotlib.__version__}.

L2G is a machine-learning prioritisation score, not proof of causality. Scores can
change between releases as the model and source evidence are updated.
"""
    (REPORTS / "stage_04_opentargets_method.md").write_text(methods, encoding="utf-8")

    report = f"""# Stage 4 methods and quality control

## Scope

Stage 4 matched the primary Parkinson GWAS `{STUDY_ID}` to Open Targets, retrieved
fine-mapped credible sets and all displayed L2G predictions, and generated broad,
stringent, and weighted gene sets. No AHBA regional scores, enrichment tests,
permutations, pathways, or disease brain maps were computed.

## Why locus-to-gene prioritisation is required

A GWAS association identifies a statistical signal at a variant, not a gene. Lead
variants can merely tag nearby correlated variants through linkage disequilibrium,
and regulatory variants can act on genes that are not the nearest gene. Fine mapping
uses association statistics and LD to distribute posterior probability across a set
of plausible variants; the resulting credible set is intended to contain the causal
variant at a stated probability, rather than declaring its lead variant causal.

Open Targets L2G applies a gradient-boosting model to rank protein-coding genes near
each credible set. Its inputs include posterior-weighted gene/TSS distance,
colocalisation with eQTL, sQTL and pQTL signals, predicted variant consequences,
enhancer-to-gene evidence, and neighbourhood/context features. Multiple genes are
retained because the model supplies relative evidence, not a definitive causal call.

## Stage 3 consistency gate

The audit confirmed 109 genome-wide significant Stage 3 associations in 78
provisional loci on GRCh38.p14, with each reported lead variant equal to the minimum
p-value association in its locus. All 887 mapped-gene rows refer to valid loci (763
unique symbols). The selected study is `{STUDY_ID}`, sample size 2,525,730, with
European, East Asian, Hispanic or Latin American, and African-unspecified discovery
ancestries. No inconsistency triggered the stop condition.

## Open Targets study match and retrieval

The Platform record has the exact accession, Parkinson's disease trait, sample size
2,525,730, and matching ancestry structure. Open Targets reports GRCh38 coordinates;
Stage 3 records GRCh38.p14, a patch-level label for the same major assembly.
{credible['stage_03_locus_id'].notna().sum()} of {len(credible)} credible-set lead
variants matched Stage 3 associations by rsID or coordinate. The Platform returned
{len(credible)} PICS credible sets and {len(l2g)} gene–credible-set L2G predictions
for {l2g['gene'].nunique()} unique genes.

## Gene-set construction

- Broad: {len(broad)} genes with a returned score > {BROAD_THRESHOLD:.2f}.
- Stringent: {len(stringent)} genes with maximum score ≥ {strict:.6f}, the eligible-score Q3.
- Weighted: {len(broad)} genes; weight is the gene's maximum score across credible sets.
- Multiple candidates were retained per locus. The full feature values and SHAP
  contributions are serialized in the master table, with separate distance,
  colocalisation, variant-effect, and other-evidence fields.

## AHBA identifier coverage (lookup only)

{coverage_markdown}

This is only an identifier-coverage check against the Stage 2 AHBA gene metadata;
expression values were not scored or aggregated by disease.

Genes were matched by exact approved gene symbol against the columns of
`data/processed/brain_region_gene_expression.csv`. This deliberately avoids alias
substitution at this stage; absent symbols may reflect microarray coverage,
reannotation, or symbol-version differences rather than absence of brain expression.

Three sets were retained so later work can test sensitivity to inclusion: broad
maximises coverage under the Platform's displayed-evidence rule, stringent focuses
on the upper quartile for this release, and weighted retains broad coverage without
treating a 0.90 and 0.08 score as equivalent.

## Caveats

All study credible sets have PICS/top-hit confidence and quality-control flags noting
the absence of in-sample LD. The Platform reports that harmonised summary statistics
for this study are unavailable or empty. `locusStart` and `locusEnd` may consequently
be absent; variant-level coordinates and posterior probabilities are preserved in
JSON. Stage 3 distance-grouped provisional loci and Open Targets credible sets are
different analytical objects and are linked only when their lead variant matches.
"""
    (REPORTS / "stage_04_methods.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    for directory in (GWAS, GENES, RESULT_DATA, FIGURES, TABLES, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    assoc, _, mapped = audit_stage_03()
    study, details, calls = get_open_targets()
    study_match = build_study_match(study, details, assoc)
    credible = build_credible_sets(details, assoc)
    l2g = build_l2g(details, credible)
    genes = aggregate_genes(l2g)
    strict_threshold = float(l2g["stringent_threshold"].iloc[0])
    broad = genes.copy()
    broad["selection_rule"] = f"maximum L2G score > {BROAD_THRESHOLD}"
    stringent = genes.loc[genes["max_l2g_score"] >= strict_threshold].copy()
    stringent["selection_rule"] = f"maximum L2G score >= eligible prediction Q3 ({strict_threshold:.6f})"
    weighted = genes.copy()
    weighted["selection_rule"] = "eligible under broad rule; gene_weight=max L2G score"
    comparison = comparison_table(mapped, broad, l2g)
    coverage = ahba_coverage(broad, stringent, weighted)
    master = master_evidence_table(l2g)
    top = top_gene_table(l2g)

    study_match.to_csv(GWAS / "parkinson_opentargets_study_match.csv", index=False)
    credible.to_csv(GWAS / "parkinson_credible_sets.csv", index=False)
    l2g.to_csv(GENES / "parkinson_l2g_predictions.csv", index=False)
    broad.to_csv(GENES / "parkinson_genes_broad.csv", index=False)
    stringent.to_csv(GENES / "parkinson_genes_stringent.csv", index=False)
    weighted.to_csv(GENES / "parkinson_genes_weighted.csv", index=False)
    master.to_csv(GENES / "parkinson_gene_locus_evidence.csv", index=False)
    comparison.to_csv(RESULT_DATA / "parkinson_gwas_vs_l2g_gene_comparison.csv", index=False)
    coverage.to_csv(RESULT_DATA / "parkinson_gene_ahba_coverage.csv", index=False)
    top.to_csv(TABLES / "parkinson_top_prioritized_genes.csv", index=False)
    make_figures(l2g, comparison)
    write_query_log(study["_meta"], calls, details)
    write_reports(study, credible, l2g, broad, stringent, coverage)

    missing = [
        str(path.relative_to(ROOT)) for path in (
            GWAS / "parkinson_opentargets_study_match.csv",
            GWAS / "parkinson_credible_sets.csv",
            GENES / "parkinson_l2g_predictions.csv",
            GENES / "parkinson_genes_broad.csv",
            GENES / "parkinson_genes_stringent.csv",
            GENES / "parkinson_genes_weighted.csv",
            FIGURES / "stage_04_parkinson_gwas_to_gene.png",
            FIGURES / "stage_04_parkinson_l2g_distribution.png",
            FIGURES / "stage_04_parkinson_gene_mapping_comparison.png",
        ) if not path.is_file()
    ]
    print("\nPARKINSON GENETIC EVIDENCE SUMMARY")
    print(f"GWAS study: {STUDY_ID}")
    print(f"Independent loci (Stage 3 provisional): {assoc['locus_id'].nunique()}")
    print(f"Credible sets: {len(credible)}")
    print(f"Total L2G gene predictions: {len(l2g)}")
    print(f"Broad gene set: {len(broad)}")
    print(f"Stringent gene set: {len(stringent)}")
    print(f"Weighted gene set: {len(weighted)}")
    print(f"Unique prioritized genes: {l2g['gene'].nunique()}")
    print(f"AHBA-represented genes (broad): {int(coverage.iloc[0]['n_present_in_ahba'])}")
    print(f"AHBA coverage (broad): {float(coverage.iloc[0]['coverage_percentage']):.1f}%")
    print(f"Highest L2G score: {l2g['l2g_score'].max():.6f}")
    print(f"Median L2G score: {l2g['l2g_score'].median():.6f}")
    print("Top 20 prioritized candidate genes:")
    for _, row in top.head(20).iterrows():
        print(f"  {int(row['rank']):2d}. {row['gene']} ({row['l2g_score']:.6f})")
    print(f"Missing outputs: {missing or 'none'}")
    print("  STOP: no regional scoring, enrichment, permutations, pathways, or brain maps performed.")


if __name__ == "__main__":
    main()
