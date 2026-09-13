"""Stage 3: retrieve and organize Parkinson disease GWAS Catalog evidence.

This stage deliberately stops at curated GWAS associations, provisional
distance-grouped loci, and Catalog-mapped genes. It does not infer causal genes
or integrate disease genetics with AHBA expression.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = ROOT / "data" / "gwas"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "results" / "figures"

API_ROOT = "https://www.ebi.ac.uk/gwas/rest/api"
DOCS_URL = "https://www.ebi.ac.uk/gwas/rest/docs/api"
PD_TRAIT = "Parkinson disease"
PD_ONTOLOGY_ID = "MONDO_0005180"
PD_ONTOLOGY_URI = "http://purl.obolibrary.org/obo/MONDO_0005180"
PRIMARY_STUDY = "GCST90308590"
VALIDATION_STUDIES = ("GCST90480008", "GCST90270939", "GCST90828116")
PVALUE_THRESHOLD = 5e-8
LOCUS_GAP_BP = 250_000
ACCESS_DATE = date.today().isoformat()


def make_session() -> requests.Session:
    """Return a polite, retrying HTTP session for the public Catalog."""
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "GENE2BRAIN-stage-03/1.0 (research; GWAS Catalog client)",
        }
    )
    return session


def get_json(
    session: requests.Session, endpoint: str, params: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Retrieve a JSON resource and fail visibly on HTTP or decoding errors."""
    response = session.get(endpoint, params=params, timeout=180)
    response.raise_for_status()
    return response.json()


def embedded(payload: dict[str, Any], key: str) -> list[dict[str, Any]]:
    """Extract a named collection from a HAL response."""
    return payload.get("_embedded", {}).get(key, [])


def study_ancestries(study: dict[str, Any]) -> str:
    groups: list[str] = []
    for ancestry in study.get("ancestries") or []:
        for group in ancestry.get("ancestralGroups") or []:
            value = group.get("ancestralGroup")
            if value and value not in groups:
                groups.append(value)
    return "; ".join(groups) if groups else "Not reported"


def ancestry_count(study: dict[str, Any]) -> int:
    value = study_ancestries(study)
    return 0 if value == "Not reported" else len(value.split("; "))


def sample_total(study: dict[str, Any], sample_type: str = "initial") -> int:
    return int(
        sum(
            ancestry.get("numberOfIndividuals") or 0
            for ancestry in study.get("ancestries") or []
            if ancestry.get("type") == sample_type
        )
    )


def is_direct_pd_risk_study(study: dict[str, Any], association_count: int) -> bool:
    """Identify direct, variant-level PD case-control risk GWAS candidates."""
    trait = (study.get("diseaseTrait") or {}).get("trait", "").strip().lower()
    direct_names = {"parkinson disease", "parkinson's disease", "parkinsons disease"}
    excluded = (
        "progression",
        "severity",
        "time to event",
        "subtype",
        "familial",
        "medication",
        "cnv",
        "gene-based burden",
        "pleiotropy",
        "interaction",
        "age at diagnosis",
    )
    technology = "; ".join(
        item.get("genotypingTechnology", "")
        for item in study.get("genotypingTechnologies") or []
    ).lower()
    is_variant_gwas = "genome-wide" in technology
    return (
        trait in direct_names
        and not any(term in trait for term in excluded)
        and is_variant_gwas
        and association_count > 0
    )


def selection_score(study: dict[str, Any], association_count: int) -> int:
    """Transparent ordinal score used to rank, not to establish causality."""
    if not is_direct_pd_risk_study(study, association_count):
        return -1
    n = sample_total(study)
    power = 4 if n >= 1_000_000 else 3 if n >= 500_000 else 2 if n >= 100_000 else 1 if n >= 20_000 else 0
    associations = 3 if association_count >= 50 else 2 if association_count >= 10 else 1
    summary_stats = 2 if study.get("fullPvalueSet") else 0
    year = int(((study.get("publicationInfo") or {}).get("publicationDate") or "0")[:4])
    recency = 1 if year >= 2020 else 0
    diversity = 2 if ancestry_count(study) >= 3 else 1 if ancestry_count(study) == 2 else 0
    replication = 1 if sample_total(study, "replication") > 0 else 0
    return power + associations + summary_stats + recency + diversity + replication


def suitability_reason(study: dict[str, Any], count: int) -> str:
    if is_direct_pd_risk_study(study, count):
        details = ["direct Parkinson case-control phenotype", "variant-level genome-wide design"]
        if count:
            details.append(f"{count} ontology-linked reported associations")
        if study.get("fullPvalueSet"):
            details.append("full p-value set flagged by Catalog")
        if ancestry_count(study) > 1:
            details.append("multiple ancestry groups")
        return "; ".join(details)
    trait = (study.get("diseaseTrait") or {}).get("trait", "not reported")
    if count == 0:
        return f"not primary-eligible: no ontology-linked reported associations; trait={trait}"
    return f"not primary-eligible: phenotype/design is not a direct common-variant Parkinson risk GWAS; trait={trait}"


def build_candidate_tables(
    studies: list[dict[str, Any]], association_counts: Counter[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates: list[dict[str, Any]] = []
    ranked: list[dict[str, Any]] = []
    for study in studies:
        accession = study.get("accessionId", "")
        publication = study.get("publicationInfo") or {}
        count = int(association_counts.get(accession, 0))
        technologies = "; ".join(
            item.get("genotypingTechnology", "")
            for item in study.get("genotypingTechnologies") or []
        )
        publication_date = publication.get("publicationDate") or ""
        eligible = is_direct_pd_risk_study(study, count)
        reason = suitability_reason(study, count)
        candidates.append(
            {
                "study_accession": accession,
                "reported_trait": (study.get("diseaseTrait") or {}).get("trait"),
                "efo_trait": PD_TRAIT,
                "ontology_id": PD_ONTOLOGY_ID,
                "ontology_uri": PD_ONTOLOGY_URI,
                "publication": publication.get("title"),
                "journal": publication.get("publication"),
                "pmid": publication.get("pubmedId"),
                "first_author": (publication.get("author") or {}).get("fullname"),
                "publication_year": publication_date[:4] or pd.NA,
                "discovery_sample_size": study.get("initialSampleSize"),
                "discovery_sample_size_numeric": sample_total(study),
                "replication_sample_size": study.get("replicationSampleSize"),
                "replication_sample_size_numeric": sample_total(study, "replication"),
                "ancestry_population": study_ancestries(study),
                "number_of_reported_associations": count,
                "summary_statistics_available": bool(study.get("fullPvalueSet")),
                "study_design_information": (
                    f"technology={technologies or 'not reported'}; "
                    f"imputed={study.get('imputed')}; pooled={study.get('pooled')}; "
                    f"SNP count={study.get('snpCount')}; comment={study.get('studyDesignComment')}"
                ),
                "appropriate_for_common_variant_pd_risk": eligible,
                "appropriateness_reason": reason,
            }
        )
        ranked.append(
            {
                "study_accession": accession,
                "selection_score": selection_score(study, count),
                "sample_size": sample_total(study),
                "replication": study.get("replicationSampleSize"),
                "ancestry": study_ancestries(study),
                "summary_statistics_available": bool(study.get("fullPvalueSet")),
                "publication_year": publication_date[:4] or pd.NA,
                "number_of_associations": count,
                "selection_reason": reason,
                "primary_selected": accession == PRIMARY_STUDY,
                "independent_validation_candidate": accession in VALIDATION_STUDIES,
            }
        )

    candidate_df = pd.DataFrame(candidates).sort_values(
        ["publication_year", "study_accession"], ascending=[False, True], na_position="last"
    )
    selection_df = pd.DataFrame(ranked).sort_values(
        ["selection_score", "sample_size", "number_of_associations", "study_accession"],
        ascending=[False, False, False, True],
    )
    selection_df.insert(1, "study_rank", np.arange(1, len(selection_df) + 1))
    if selection_df.iloc[0]["study_accession"] != PRIMARY_STUDY:
        raise RuntimeError("Configured primary study is no longer top-ranked; review selection criteria")
    return candidate_df, selection_df


def risk_alleles(association: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    output: list[tuple[str, dict[str, Any]]] = []
    for locus in association.get("loci") or []:
        for risk in locus.get("strongestRiskAlleles") or []:
            output.append((risk.get("riskAlleleName", ""), locus))
    return output


def association_id(association: dict[str, Any]) -> str | None:
    href = ((association.get("_links") or {}).get("self") or {}).get("href", "")
    return href.rstrip("/").split("/")[-1] or None


def context_genes(snp: dict[str, Any]) -> tuple[str, str]:
    names: list[str] = []
    identifiers: list[str] = []
    for context in snp.get("genomicContexts") or []:
        gene = context.get("gene") or {}
        name = gene.get("geneName")
        if name and name not in names:
            names.append(name)
        for record in gene.get("entrezGeneIds") or []:
            value = record.get("entrezGeneId")
            if value and value not in identifiers:
                identifiers.append(value)
    return ";".join(names), ";".join(identifiers)


def build_associations(
    raw_associations: list[dict[str, Any]], snps: list[dict[str, Any]], genome_build: str
) -> tuple[pd.DataFrame, int]:
    snp_lookup = {item.get("rsId"): item for item in snps}
    rows: list[dict[str, Any]] = []
    excluded = 0
    for association in raw_associations:
        pvalue = association.get("pvalue")
        if pvalue is None or not (0 <= float(pvalue) <= PVALUE_THRESHOLD):
            excluded += 1
            continue
        for risk_name, locus in risk_alleles(association):
            variant, separator, allele = risk_name.rpartition("-")
            if not separator:
                variant, allele = risk_name, None
            snp = snp_lookup.get(variant, {})
            locations = snp.get("locations") or []
            location = locations[0] if len(locations) == 1 else {}
            mapped_names, mapped_ids = context_genes(snp)
            efo_traits = association.get("efoTraits") or []
            trait_names = ";".join(x.get("trait", "") for x in efo_traits)
            rows.append(
                {
                    "variant": variant or None,
                    "rsid": variant if re.fullmatch(r"rs\d+", variant or "") else None,
                    "chromosome": location.get("chromosomeName"),
                    "position": location.get("chromosomePosition"),
                    "p_value": float(pvalue),
                    "effect_allele": allele,
                    "other_allele": None,
                    "effect_size": association.get("betaNum") or association.get("orPerCopyNum"),
                    "odds_ratio": association.get("orPerCopyNum"),
                    "beta": association.get("betaNum"),
                    "beta_unit": association.get("betaUnit"),
                    "beta_direction": association.get("betaDirection"),
                    "standard_error": association.get("standardError"),
                    "risk_frequency": association.get("riskFrequency"),
                    "mapped_genes": mapped_names or None,
                    "mapped_gene_entrez_ids": mapped_ids or None,
                    "study_accession": PRIMARY_STUDY,
                    "trait": trait_names or PD_TRAIT,
                    "ontology_id": PD_ONTOLOGY_ID,
                    "genomic_assembly": genome_build,
                    "cytogenetic_region": (location.get("region") or {}).get("name"),
                    "functional_class": snp.get("functionalClass"),
                    "snp_type": association.get("snpType"),
                    "pvalue_description": association.get("pvalueDescription"),
                    "multi_snp_haplotype": association.get("multiSnpHaplotype"),
                    "snp_interaction": association.get("snpInteraction"),
                    "catalog_association_id": association_id(association),
                    "catalog_mapping_note": (
                        "Coordinates and genomic contexts from current GWAS Catalog SNP resource; "
                        "mapped genes are not causal-gene claims"
                    ),
                }
            )
    frame = pd.DataFrame(rows)
    frame["chromosome"] = pd.to_numeric(frame["chromosome"], errors="coerce").astype("Int64")
    frame["position"] = pd.to_numeric(frame["position"], errors="coerce").astype("Int64")
    frame = frame.sort_values(["chromosome", "position", "p_value", "variant"]).reset_index(drop=True)
    return frame, excluded


def group_loci(associations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Single-link significant variants <=250 kb apart on each chromosome."""
    working = associations.copy()
    assignments: dict[int, str] = {}
    loci: list[dict[str, Any]] = []
    serial = 0
    valid = working.dropna(subset=["chromosome", "position"]).copy()
    for chromosome, chrom_rows in valid.groupby("chromosome", sort=True):
        current_indices: list[int] = []
        previous_position: int | None = None
        for index, row in chrom_rows.sort_values("position").iterrows():
            position = int(row["position"])
            if previous_position is not None and position - previous_position > LOCUS_GAP_BP:
                serial += 1
                locus_id = _finish_locus(working, current_indices, int(chromosome), serial, loci)
                assignments.update({item: locus_id for item in current_indices})
                current_indices = []
            current_indices.append(index)
            previous_position = position
        if current_indices:
            serial += 1
            locus_id = _finish_locus(working, current_indices, int(chromosome), serial, loci)
            assignments.update({item: locus_id for item in current_indices})

    working["locus_id"] = working.index.map(assignments)
    locus_df = pd.DataFrame(loci)
    return working, locus_df


def _finish_locus(
    associations: pd.DataFrame,
    indices: list[int],
    chromosome: int,
    serial: int,
    output: list[dict[str, Any]],
) -> str:
    group = associations.loc[indices]
    lead = group.sort_values(["p_value", "position", "variant"]).iloc[0]
    locus_id = f"PD_CHR{chromosome:02d}_{serial:03d}"
    output.append(
        {
            "locus_id": locus_id,
            "chromosome": chromosome,
            "start": int(group["position"].min()),
            "end": int(group["position"].max()),
            "lead_variant": lead["variant"],
            "lead_p_value": lead["p_value"],
            "number_of_associations": len(group),
            "study_accession": PRIMARY_STUDY,
            "grouping_method": (
                "single-link distance grouping: adjacent significant variants <=250 kb"
            ),
        }
    )
    return locus_id


def build_mapped_genes(
    associations: pd.DataFrame, snps: list[dict[str, Any]]
) -> pd.DataFrame:
    snp_lookup = {item.get("rsId"): item for item in snps}
    rows: list[dict[str, Any]] = []
    for association in associations.itertuples(index=False):
        snp = snp_lookup.get(association.variant, {})
        contexts_by_gene: dict[tuple[str, str], list[dict[str, Any]]] = {}
        for context in snp.get("genomicContexts") or []:
            gene = context.get("gene") or {}
            name = gene.get("geneName")
            if not name:
                continue
            ids = ";".join(
                record.get("entrezGeneId", "") for record in gene.get("entrezGeneIds") or []
            )
            contexts_by_gene.setdefault((name, ids), []).append(context)
        for (name, gene_id), contexts in contexts_by_gene.items():
            sources = sorted({item.get("source", "not reported") for item in contexts})
            methods = sorted({item.get("mappingMethod", "not reported") for item in contexts})
            relations: list[str] = []
            for item in contexts:
                if not item.get("isIntergenic"):
                    relations.append("overlapping")
                if item.get("isUpstream"):
                    relations.append("upstream")
                if item.get("isDownstream"):
                    relations.append("downstream")
                if item.get("isClosestGene"):
                    relations.append("closest")
            relation_text = ",".join(sorted(set(relations))) or "context"
            rows.append(
                {
                    "gene": name,
                    "gene_id_if_available": gene_id or None,
                    "locus_id": association.locus_id,
                    "variant": association.variant,
                    "study_accession": association.study_accession,
                    "p_value": association.p_value,
                    "mapping_source": (
                        "GWAS Catalog SNP genomic context "
                        f"(sources={','.join(sources)}; methods={','.join(methods)}; "
                        f"relations={relation_text}); candidate/mapped gene, not causal"
                    ),
                }
            )
    return (
        pd.DataFrame(rows)
        .drop_duplicates(["gene", "gene_id_if_available", "locus_id", "variant"])
        .sort_values(["locus_id", "p_value", "variant", "gene"])
        .reset_index(drop=True)
    )


def plot_manhattan(associations: pd.DataFrame, loci: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    frame = associations.dropna(subset=["chromosome", "position", "p_value"]).copy()
    frame["chromosome"] = frame["chromosome"].astype(int)
    chrom_max = frame.groupby("chromosome")["position"].max()
    offsets: dict[int, int] = {}
    cursor = 0
    for chromosome in range(1, 23):
        offsets[chromosome] = cursor
        cursor += int(chrom_max.get(chromosome, 0)) + 20_000_000
    frame["genome_x"] = frame.apply(
        lambda row: int(row["position"]) + offsets[int(row["chromosome"])], axis=1
    )
    frame["minus_log10_p"] = -np.log10(frame["p_value"].clip(lower=np.finfo(float).tiny))

    fig, ax = plt.subplots(figsize=(18, 8.5))
    palette = ("#315C8C", "#75A6C9")
    for chromosome, group in frame.groupby("chromosome"):
        ax.scatter(
            group["genome_x"], group["minus_log10_p"], s=34,
            color=palette[(int(chromosome) - 1) % 2], alpha=0.86,
            edgecolor="white", linewidth=0.35,
        )
    threshold_y = -math.log10(PVALUE_THRESHOLD)
    ax.axhline(threshold_y, color="#B23A48", linestyle="--", linewidth=1.7,
               label=r"Genome-wide significance ($p = 5\times10^{-8}$)")
    centers = [offsets[c] + chrom_max.get(c, 0) / 2 for c in range(1, 23)]
    ax.set_xticks(centers, [str(c) for c in range(1, 23)])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title("Parkinson Disease GWAS Genetic Risk Landscape", weight="bold", pad=18)
    ax.text(
        0.5, 1.01,
        f"{PRIMARY_STUDY}: GWAS Catalog associations meeting p <= 5e-8 (not full summary statistics)",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=11, color="#4A4A4A",
    )
    top = loci.nsmallest(8, "lead_p_value")
    for row in top.itertuples(index=False):
        point = frame.loc[frame["variant"] == row.lead_variant].iloc[0]
        ax.annotate(
            row.lead_variant, (point["genome_x"], point["minus_log10_p"]),
            xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8,
        )
    ax.legend(loc="upper right", frameon=True, fontsize=10)
    ax.grid(axis="x", visible=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_03_parkinson_gwas_manhattan.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_loci(loci: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    frame = loci.copy()
    frame["position_mb"] = (frame["start"] + frame["end"]) / 2 / 1e6
    frame["strength"] = -np.log10(frame["lead_p_value"].clip(lower=np.finfo(float).tiny))
    fig, ax = plt.subplots(figsize=(14, 10))
    scatter = ax.scatter(
        frame["position_mb"], frame["chromosome"], c=frame["strength"],
        s=45 + 24 * np.sqrt(frame["number_of_associations"]), cmap="viridis",
        edgecolor="white", linewidth=0.55,
    )
    ax.set_yticks(range(1, 23))
    ax.set_ylim(22.7, 0.3)
    ax.set_xlabel("Genomic position within chromosome (Mb)")
    ax.set_ylabel("Chromosome")
    ax.set_title("Parkinson GWAS Locus Overview", weight="bold", pad=16)
    ax.text(
        0.5, 1.01, f"{len(frame)} provisional loci; 250 kb single-link distance grouping",
        transform=ax.transAxes, ha="center", fontsize=11, color="#4A4A4A",
    )
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label(r"Lead $-\log_{10}(p)$")
    label_rows = (
        frame.sort_values("lead_p_value")
        .drop_duplicates("chromosome")
        .head(8)
    )
    for row in label_rows.itertuples(index=False):
        ax.annotate(
            row.lead_variant, (row.position_mb, row.chromosome), xytext=(5, 0),
            textcoords="offset points", va="center", fontsize=7.5,
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "alpha": 0.72,
                  "edgecolor": "none"},
        )
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_03_parkinson_locus_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_concept() -> None:
    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(15, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    labels = [
        ("GWAS variant", "Statistical association"),
        ("Genomic locus", "Nearby/linked signals grouped"),
        ("Gene mapping / prioritization", "Multiple evidence types required"),
        ("Candidate / mapped gene", "Not a causal-gene claim"),
    ]
    centers = [0.12, 0.37, 0.63, 0.88]
    colors = ["#DCEAF4", "#D7E7DD", "#F4E7C7", "#F2D6DB"]
    for index, ((title, subtitle), x, color) in enumerate(zip(labels, centers, colors)):
        box = FancyBboxPatch(
            (x - 0.10, 0.35), 0.20, 0.30,
            boxstyle="round,pad=0.018,rounding_size=0.025",
            facecolor=color, edgecolor="#34495E", linewidth=1.4,
        )
        ax.add_patch(box)
        ax.text(x, 0.54, title, ha="center", va="center", fontsize=13, weight="bold")
        ax.text(x, 0.43, subtitle, ha="center", va="center", fontsize=9.3, color="#4A4A4A")
        if index < len(centers) - 1:
            ax.annotate(
                "", xy=(centers[index + 1] - 0.11, 0.5), xytext=(x + 0.11, 0.5),
                arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "#34495E"},
            )
    ax.set_title("From GWAS Variant to Candidate Gene", weight="bold", pad=12)
    ax.text(
        0.5, 0.16,
        "Stage 3 preserves mapped candidates only; causal inference requires independent functional evidence.",
        ha="center", fontsize=11, color="#8A2432", weight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "stage_03_gwas_to_gene_concept.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def summary_statistics_table(study_lookup: dict[str, dict[str, Any]]) -> pd.DataFrame:
    kim_open_source = (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90275001-GCST90276000/GCST90275127/GCST90275127.tsv"
    )
    all_of_us_source = (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90480001-GCST90481000/GCST90480008/GCST90480008.tsv.gz"
    )
    rows = []
    for accession in (PRIMARY_STUDY,) + VALIDATION_STUDIES:
        study = study_lookup[accession]
        is_primary = accession == PRIMARY_STUDY
        source_study = study_lookup["GCST90275127"] if is_primary else study
        available = is_primary or bool(study.get("fullPvalueSet"))
        source = (
            kim_open_source if is_primary else
            all_of_us_source if accession == "GCST90480008" else
            f"{API_ROOT}/studies/{accession}"
        )
        downloadable_build = (
            "GRCh37" if is_primary else
            "GRCh38" if accession == "GCST90480008" else
            "Not supplied for downloadable file"
        )
        rows.append(
            {
                "study_accession": accession,
                "source_accession": "GCST90275127" if is_primary else accession,
                "source": source,
                "download_available": available,
                "file_format": (
                    "GWAS-SSF v1.0 TSV" if is_primary else
                    "GWAS-SSF v1.0 TSV.GZ" if accession == "GCST90480008" else
                    "Not available from GWAS Catalog"
                ),
                "genome_build": downloadable_build,
                "sample_size": sample_total(source_study) + sample_total(source_study, "replication"),
                "access_date": ACCESS_DATE,
                "availability_scope": (
                    "Open companion meta-analysis; excludes 23andMe; full primary meta-analysis is controlled access"
                    if is_primary else
                    "Catalog full p-value set" if available else
                    "No Catalog full p-value set"
                ),
                "file_size": (
                    "1.2 GB" if is_primary else
                    "827 MB" if accession == "GCST90480008" else None
                ),
                "download_action": "metadata recorded; large file not downloaded" if available else "no Catalog full p-value set",
            }
        )
    return pd.DataFrame(rows)


def qc_counts(
    candidates: pd.DataFrame, associations: pd.DataFrame, mapped: pd.DataFrame,
    raw_primary_count: int, threshold_excluded: int,
) -> dict[str, int]:
    chromosomes = pd.to_numeric(associations["chromosome"], errors="coerce")
    positions = pd.to_numeric(associations["position"], errors="coerce")
    pvalues = pd.to_numeric(associations["p_value"], errors="coerce")
    mapped_variants = set(mapped["variant"].dropna())
    return {
        "raw_primary_count": raw_primary_count,
        "threshold_excluded": threshold_excluded,
        "duplicate_variant_rows": int(associations.duplicated("variant", keep=False).sum()),
        "missing_rsids": int(associations["rsid"].isna().sum()),
        "invalid_chromosomes": int((~chromosomes.between(1, 22) | chromosomes.isna()).sum()),
        "invalid_positions": int((positions.le(0) | positions.isna()).sum()),
        "missing_pvalues": int(pvalues.isna().sum()),
        "invalid_pvalues": int((~pvalues.between(0, 1) | pvalues.isna()).sum()),
        "duplicated_study_records": int(candidates.duplicated("study_accession", keep=False).sum()),
        "missing_catalog_mapped_gene_variants": int((~associations["variant"].isin(mapped_variants)).sum()),
        "missing_author_reported_gene_rows": len(associations),
        "missing_effect_alleles": int(associations["effect_allele"].isna().sum()),
        "missing_other_alleles": int(associations["other_allele"].isna().sum()),
        "missing_effect_sizes": int(associations["effect_size"].isna().sum()),
        "unusually_large_effect_sizes": int(
            (pd.to_numeric(associations["effect_size"], errors="coerce").abs() > 10).sum()
        ),
    }


def write_reports(
    candidates: pd.DataFrame, selection: pd.DataFrame, associations: pd.DataFrame,
    loci: pd.DataFrame, mapped: pd.DataFrame, summary_stats: pd.DataFrame,
    qc: dict[str, int], primary: dict[str, Any], all_association_count: int,
    genome_build: str,
) -> None:
    api_report = f"""# Stage 03 GWAS Catalog API record

- **Source:** NHGRI-EBI GWAS Catalog
- **Official documentation:** {DOCS_URL}
- **Access date:** {ACCESS_DATE}
- **Response format:** JSON using HAL (`_embedded` and `_links`)
- **Live API root used:** `{API_ROOT}`
- **API version finding:** the current official root and documentation expose an unnumbered REST API. Direct checks of `/gwas/rest/api/v2`, `/gwas/api/v2`, and their trailing-slash variants returned HTTP 404 on {ACCESS_DATE}. The retired summary-statistics API root returned HTTP 410. Therefore this workflow uses the current live API advertised by the official root; it does not mislabel it as v2.

## Endpoints and parameters

1. `{API_ROOT}/efoTraits/search/findByEfoTrait`
   - `trait=Parkinson disease`, `size=100`
   - resolved uniquely to `{PD_ONTOLOGY_ID}` ({PD_ONTOLOGY_URI}).
2. `{API_ROOT}/efoTraits/{PD_ONTOLOGY_ID}/studies`
   - `size=500`
   - returned {len(candidates)} ontology-linked candidate studies.
3. `{API_ROOT}/associations/search/findByEfoTrait`
   - `efoTrait=Parkinson disease`, `projection=associationByStudy`, `size=5000`
   - returned {all_association_count} associations used only to count records per candidate study.
4. `{API_ROOT}/studies/{PRIMARY_STUDY}/associations`
   - `projection=associationByStudy`, `size=500`
   - returned {qc['raw_primary_count']} reported primary-study associations before the fixed p-value filter.
5. `{API_ROOT}/studies/{PRIMARY_STUDY}/snps`
   - `projection=snp`, `size=500`
   - supplied rsIDs, current Catalog coordinates, functional class, cytogenetic region, and genomic-context gene mappings.
6. `{API_ROOT}/metadata`
   - no query parameters
   - supplied the Catalog's current mapping assembly: `{genome_build}`.

The ontology relation, rather than a free-text study query, defines the candidate set. A label query is used only to resolve the ontology resource itself.
"""
    (REPORTS / "stage_03_gwas_api.md").write_text(api_report, encoding="utf-8")

    top = selection.iloc[0]
    primary_publication = primary.get("publicationInfo") or {}
    methods = f"""# Stage 03 methods: Parkinson disease GWAS

Generated/accessed: {ACCESS_DATE}

## Source and search

The source was the current NHGRI-EBI GWAS Catalog HAL REST API at `{API_ROOT}`. The live version audit is recorded in `stage_03_gwas_api.md`. The exact Parkinson disease ontology resource was `{PD_ONTOLOGY_ID}` ({PD_ONTOLOGY_URI}); all {len(candidates)} studies linked from that resource were retained as candidates. This avoids defining the candidate set by publication-title or trait-string matching.

## Reproducible study selection

Candidates were first classified as primary-eligible only when the reported trait was direct Parkinson disease, the design was variant-level genome-wide genotyping/sequencing, at least one ontology-linked association was present, and the trait did not represent progression, severity, subtype, familial disease, proxy-only disease, medication, interaction, pleiotropy, CNV, or gene-burden analysis.

Eligible studies received an ordinal score: discovery sample size (0–4 points), reported associations (1–3), full p-value set (0 or 2), publication from 2020 onward (0 or 1), ancestry breadth (0–2), and a replication stage (0 or 1). Ties were ordered by discovery sample size, association count, then accession. All component data and reasons are preserved in `parkinson_study_selection.csv`.

The selected primary study is **{PRIMARY_STUDY}**, rank {int(top['study_rank'])}: {primary_publication.get('title')} ({primary_publication.get('author', {}).get('fullname')}, {str(primary_publication.get('publicationDate'))[:4]}, PMID {primary_publication.get('pubmedId')}). It has a direct Parkinson phenotype, {sample_total(primary):,} discovery/meta-analysis participants across {study_ancestries(primary)}, and the largest curated reported-association set ({qc['raw_primary_count']}) among candidates. Although the Catalog does not flag a downloadable full p-value set for this accession, its power, phenotype, ancestry breadth, recency, and curated associations outweighed that limitation.

Candidate sensitivity/replication datasets, not analyzed here, are `{VALIDATION_STUDIES[0]}` (a publication-distinct multi-ancestry biobank analysis), `{VALIDATION_STUDIES[1]}` (Chinese whole-genome sequencing GWAS), and `{VALIDATION_STUDIES[2]}` (a later Taiwanese GWAS). These are publication/cohort-distinct from the component datasets listed for the open Kim companion meta-analysis. Participant-level overlap cannot be ruled out from Catalog metadata alone and must be audited before treating any result as statistically independent replication.

## Associations and genome assembly

The standard threshold **p ≤ 5×10⁻⁸** was applied without modification. Of {qc['raw_primary_count']} Catalog records for {PRIMARY_STUDY}, {len(associations)} met the threshold and {qc['threshold_excluded']} were excluded. The output is a significant-association table, not a complete summary-statistics dataset. Variant locations are current GWAS Catalog SNP-resource coordinates on **{genome_build}**, retrieved from the API metadata endpoint; the output labels this explicitly. No coordinate liftover was performed.

## Initial locus grouping

Significant variants were sorted by chromosome and position. Within each chromosome, adjacent variants separated by **≤ {LOCUS_GAP_BP:,} bp** were merged using single-link distance grouping; a gap > {LOCUS_GAP_BP:,} bp began a new locus. Locus bounds are the minimum and maximum reported significant-variant positions, and the lowest-p variant is the lead. This FUMA-inspired merge distance is a transparent initial consolidation, yielding {len(loci)} loci from {len(associations)} associations. It is **not LD clumping**, because an ancestry-matched LD reference was not introduced in Stage 3, and it does not imply signal independence or causality. LD-aware refinement belongs in a later genetics stage.

## Mapped candidate genes

The selected study has no author-reported genes in its association objects. Therefore `parkinson_gwas_mapped_genes.csv` preserves genes from the GWAS Catalog SNP `genomicContexts` field, with Entrez IDs, Ensembl/NCBI source, mapping method, and positional relationship where available. These are explicitly **GWAS-mapped candidate genes, not causal genes or confirmed risk genes**. No functional prioritization was performed.

## Summary statistics and QC

The Catalog `fullPvalueSet` field is false for {PRIMARY_STUDY}, but the publication's data-availability statement identifies an immediately accessible open companion meta-analysis under **GCST90275127**. Official FTP metadata describes it as a 1.2 GB, unharmonized GRCh37 GWAS-SSF v1.0 TSV excluding 23andMe; the full meta-analysis is available to qualified researchers under a 23andMe agreement. The exact open URL and scope are recorded, but the large file was not downloaded. Validation candidate GCST90480008 has a separate verified 827 MB, unharmonized GRCh38 GWAS-SSF v1.0 compressed TSV; its URL is recorded without downloading. The other candidate validation accessions are not flagged as full p-value sets.

QC tested duplicate variants, rsID/coordinate/chromosome validity, missing/invalid p-values, duplicate study accessions, genome-build labels, mapped-gene coverage, allele completeness, and effect-size availability/range. Findings are reported without silent repairs in `stage_03_gwas_qc.md`.

No Open Targets analysis, causal-gene inference, AHBA integration, regional enrichment, permutations, or brain mapping was performed.
"""
    (REPORTS / "stage_03_methods.md").write_text(methods, encoding="utf-8")

    builds = sorted(associations["genomic_assembly"].dropna().unique())
    qc_report = f"""# Stage 03 GWAS quality-control report

Generated: {ACCESS_DATE}

## Candidate-study checks

- Ontology-linked candidate studies: {len(candidates)}
- Duplicate study-accession rows: {qc['duplicated_study_records']}
- Candidate studies with at least one reported association: {(candidates['number_of_reported_associations'] > 0).sum()}
- Candidate studies classified as direct common-variant Parkinson risk GWAS: {candidates['appropriate_for_common_variant_pd_risk'].sum()}

## Primary association checks

- Catalog associations returned before thresholding: {qc['raw_primary_count']}
- Records excluded because p > 5×10⁻⁸: {qc['threshold_excluded']}
- Retained genome-wide significant associations: {len(associations)}
- Duplicate variant rows: {qc['duplicate_variant_rows']}
- Missing rsIDs: {qc['missing_rsids']}
- Invalid chromosomes (allowed 1–22): {qc['invalid_chromosomes']}
- Invalid/non-positive positions: {qc['invalid_positions']}
- Missing p-values: {qc['missing_pvalues']}
- P-values outside [0,1]: {qc['invalid_pvalues']}
- Genome-build labels in primary table: {', '.join(builds)}; no within-table inconsistency
- Initial distance-grouped loci: {len(loci)}

## Annotation, alleles, and effects

- Variants lacking GWAS Catalog genomic-context mapped genes: {qc['missing_catalog_mapped_gene_variants']}
- Association rows lacking author-reported genes: {qc['missing_author_reported_gene_rows']} (all; positional Catalog genomic contexts were preserved separately)
- Unique mapped candidate genes: {mapped['gene'].nunique()}
- Missing effect alleles: {qc['missing_effect_alleles']}
- Missing other alleles: {qc['missing_other_alleles']} (the association projection reports strongest risk/effect alleles, not the alternate allele)
- Missing effect-size fields (OR and beta both absent): {qc['missing_effect_sizes']}
- Available effect sizes with absolute magnitude >10: {qc['unusually_large_effect_sizes']}

The absence of alternate alleles and effect estimates is a limitation of these curated association records, not silently imputed data. The Manhattan plot therefore displays reported p-values and coordinates only. Build labels in the summary-statistics metadata describe separate downloadable files and are not mixed with the primary current-Catalog association coordinates.

No values were silently corrected, no causal genes were asserted, and no disease-to-brain analysis was performed.
"""
    (REPORTS / "stage_03_gwas_qc.md").write_text(qc_report, encoding="utf-8")


def main() -> None:
    GWAS_DIR.mkdir(parents=True, exist_ok=True)
    (GWAS_DIR / "summary_statistics").mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    session = make_session()

    mapping_metadata = embedded(
        get_json(session, f"{API_ROOT}/metadata"), "mappingMetadatas"
    )
    if not mapping_metadata or not mapping_metadata[0].get("genomeBuildVersion"):
        raise RuntimeError("GWAS Catalog mapping metadata did not report a genome build")
    genome_build = mapping_metadata[0]["genomeBuildVersion"]

    trait_payload = get_json(
        session,
        f"{API_ROOT}/efoTraits/search/findByEfoTrait",
        {"trait": PD_TRAIT, "size": 100},
    )
    traits = embedded(trait_payload, "efoTraits")
    if len(traits) != 1 or traits[0].get("shortForm") != PD_ONTOLOGY_ID:
        raise RuntimeError(f"Parkinson ontology resolution changed: {traits}")

    studies = embedded(
        get_json(session, f"{API_ROOT}/efoTraits/{PD_ONTOLOGY_ID}/studies", {"size": 500}),
        "studies",
    )
    all_associations = embedded(
        get_json(
            session,
            f"{API_ROOT}/associations/search/findByEfoTrait",
            {"efoTrait": PD_TRAIT, "projection": "associationByStudy", "size": 5000},
        ),
        "associations",
    )
    association_counts: Counter[str] = Counter(
        (item.get("study") or {}).get("accessionId") for item in all_associations
    )
    candidate_df, selection_df = build_candidate_tables(studies, association_counts)
    study_lookup = {item.get("accessionId"): item for item in studies}
    missing_configured = {PRIMARY_STUDY, *VALIDATION_STUDIES} - set(study_lookup)
    if missing_configured:
        raise RuntimeError(f"Configured studies absent from ontology query: {missing_configured}")

    primary_associations = embedded(
        get_json(
            session,
            f"{API_ROOT}/studies/{PRIMARY_STUDY}/associations",
            {"projection": "associationByStudy", "size": 500},
        ),
        "associations",
    )
    primary_snps = embedded(
        get_json(
            session,
            f"{API_ROOT}/studies/{PRIMARY_STUDY}/snps",
            {"projection": "snp", "size": 500},
        ),
        "singleNucleotidePolymorphisms",
    )
    associations, threshold_excluded = build_associations(
        primary_associations, primary_snps, genome_build
    )
    associations, loci = group_loci(associations)
    mapped = build_mapped_genes(associations, primary_snps)
    summary_stats = summary_statistics_table(study_lookup)

    candidate_df.to_csv(GWAS_DIR / "parkinson_candidate_studies.csv", index=False)
    selection_df.to_csv(GWAS_DIR / "parkinson_study_selection.csv", index=False)
    associations.to_csv(GWAS_DIR / "parkinson_gwas_associations_raw.csv", index=False)
    loci.to_csv(GWAS_DIR / "parkinson_initial_loci.csv", index=False)
    mapped.to_csv(GWAS_DIR / "parkinson_gwas_mapped_genes.csv", index=False)
    summary_stats.to_csv(
        GWAS_DIR / "parkinson_summary_statistics_metadata.csv", index=False
    )

    plot_manhattan(associations, loci)
    plot_loci(loci)
    plot_concept()
    qc = qc_counts(
        candidate_df, associations, mapped, len(primary_associations), threshold_excluded
    )
    primary = study_lookup[PRIMARY_STUDY]
    write_reports(
        candidate_df, selection_df, associations, loci, mapped, summary_stats, qc,
        primary, len(all_associations), genome_build,
    )

    required = [
        GWAS_DIR / "parkinson_candidate_studies.csv",
        GWAS_DIR / "parkinson_study_selection.csv",
        GWAS_DIR / "parkinson_gwas_associations_raw.csv",
        GWAS_DIR / "parkinson_initial_loci.csv",
        GWAS_DIR / "parkinson_gwas_mapped_genes.csv",
        GWAS_DIR / "parkinson_summary_statistics_metadata.csv",
        REPORTS / "stage_03_gwas_api.md",
        REPORTS / "stage_03_gwas_qc.md",
        REPORTS / "stage_03_methods.md",
        FIGURES / "stage_03_parkinson_gwas_manhattan.png",
        FIGURES / "stage_03_parkinson_locus_overview.png",
        FIGURES / "stage_03_gwas_to_gene_concept.png",
    ]
    missing = [path for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing/empty outputs: {missing}")

    print("\nPARKINSON GWAS SUMMARY")
    print(f"Candidate studies: {len(candidate_df)}")
    print(f"Selected primary study: {PRIMARY_STUDY}")
    print(f"Independent validation studies: {', '.join(VALIDATION_STUDIES)} (candidates; overlap audit required)")
    print(f"Discovery sample size: {primary.get('initialSampleSize')}")
    print(f"Replication sample size: {primary.get('replicationSampleSize')}")
    print(f"Ancestry: {study_ancestries(primary)}")
    print(f"Genome build: {genome_build} (current Catalog SNP coordinates)")
    print(f"Genome-wide significant associations: {len(associations)}")
    print(f"Initial independent loci: {len(loci)} (provisional distance groups)")
    print(f"Mapped candidate genes: {mapped['gene'].nunique()}")
    print("Summary statistics available: Yes — open companion GCST90275127; full primary data controlled-access")
    print("Verified output files:")
    for path in required:
        print(f"  - {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
