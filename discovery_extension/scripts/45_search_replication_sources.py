#!/usr/bin/env python3
"""Run and retain the result-free independent-source search for 44 phenotypes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import ssl
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api/studies/search"
FINNGEN_API = "https://r13.finngen.fi/api/phenos"

# Terms are exact GWAS Catalog disease-trait/EFO-label searches. Multiple terms are
# retained where UK spelling, coding labels, or a common scientific synonym differs.
SEARCH_TERMS: dict[str, list[str]] = {
    "panukbb_phecode__785__both_sexes__na__na": ["abdominal pain"],
    "panukbb_categorical__6159__both_sexes__7__na": ["pain", "chronic pain"],
    "panukbb_prescriptions__tramadol__both_sexes__na__na": ["tramadol use", "tramadol"],
    "panukbb_icd10__k21__both_sexes__na__na": ["gastroesophageal reflux disease"],
    "panukbb_categorical__1538__both_sexes__1__na": ["dietary change"],
    "panukbb_categorical__6145__both_sexes__6__na": ["stressful life events", "stress"],
    "panukbb_phecode__550_2__both_sexes__na__na": ["diaphragmatic hernia"],
    "panukbb_prescriptions__paracetamol__both_sexes__na__na": ["acetaminophen use", "paracetamol use"],
    "panukbb_phecode__558__both_sexes__na__na": ["noninfective gastroenteritis", "gastroenteritis"],
    "panukbb_icd10__z88__both_sexes__na__na": ["drug allergy", "adverse drug reaction"],
    "panukbb_phecode__496_2__both_sexes__na__na": ["chronic bronchitis"],
    "panukbb_prescriptions__orlistat__both_sexes__na__na": ["orlistat use", "orlistat"],
    "panukbb_phecode__318__both_sexes__na__na": ["tobacco use disorder", "smoking dependence"],
    "panukbb_continuous__20520__both_sexes__na__na": ["worry", "generalized anxiety disorder"],
    "panukbb_continuous__20515__both_sexes__na__na": ["trouble relaxing", "generalized anxiety disorder"],
    "panukbb_icd10__j44__both_sexes__na__na": ["chronic obstructive pulmonary disease"],
    "panukbb_prescriptions__gabapentin__both_sexes__na__na": ["gabapentin use", "gabapentin"],
    "panukbb_icd10__m17__both_sexes__na__na": ["knee osteoarthritis", "gonarthrosis"],
    "panukbb_phecode__740__both_sexes__na__na": ["osteoarthritis", "arthrosis"],
    "panukbb_prescriptions__citalopram__both_sexes__na__na": ["citalopram use", "citalopram"],
    "panukbb_phecode__960_2__both_sexes__na__na": ["penicillin allergy"],
    "panukbb_categorical__20117__both_sexes__2__na": ["alcohol drinking status", "alcohol consumption"],
    "panukbb_categorical__6138__both_sexes__5__na": ["educational attainment", "qualifications"],
    "panukbb_continuous__874__both_sexes__na__irnt": ["walking duration", "walking"],
    "panukbb_icd10__m51__both_sexes__na__na": ["intervertebral disc disorder", "intervertebral disc degeneration"],
    "panukbb_prescriptions__bendroflumethiazide__both_sexes__na__na": ["bendroflumethiazide use", "bendroflumethiazide"],
    "panukbb_categorical__6164__both_sexes__5__na": ["physical activity"],
    "panukbb_icd10__z95__both_sexes__na__na": ["cardiac implant", "vascular implant"],
    "panukbb_prescriptions__mebeverine__both_sexes__na__na": ["mebeverine use", "mebeverine"],
    "panukbb_continuous__20485__both_sexes__na__na": ["suicidal ideation", "self harm ideation"],
    "panukbb_continuous__398__both_sexes__na__na": ["cognitive performance", "matching task"],
    "panukbb_prescriptions__amlodipine__both_sexes__na__na": ["amlodipine use", "amlodipine"],
    "panukbb_prescriptions__clarithromycin__both_sexes__na__na": ["clarithromycin use", "clarithromycin"],
    "panukbb_categorical__6149__both_sexes__4__na": ["dental problems", "oral health"],
    "panukbb_phecode__594__both_sexes__na__na": ["urolithiasis", "kidney stone"],
    "panukbb_categorical__6142__both_sexes__6__na": ["employment status"],
    "panukbb_continuous__1588__both_sexes__na__na": ["beer consumption", "alcohol consumption"],
    "panukbb_icd10__k80__both_sexes__na__na": ["cholelithiasis", "gallstone disease"],
    "panukbb_categorical__41200__both_sexes__j183__na": ["operative procedure", "surgical procedure"],
    "panukbb_categorical__6151__both_sexes__7__na": ["bone fracture", "fracture"],
    "panukbb_categorical__1418__both_sexes__1__na": ["milk type", "milk consumption"],
    "panukbb_continuous__699__both_sexes__na__irnt": ["residential duration", "length of residence"],
    "panukbb_categorical__4294__both_sexes__0__na": ["cognitive performance"],
    "panukbb_continuous__20409__both_sexes__na__na": ["alcohol-related guilt", "alcohol consumption"],
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def one_line(value: object) -> str:
    """Keep remotely supplied metadata inside one auditable TSV cell."""
    return " ".join(str(value if value not in {None, ""} else "NA").split())


def ssl_context() -> ssl.SSLContext:
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def get_json(url: str) -> object:
    request = urllib.request.Request(url, headers={"User-Agent": "sleep-gwas-atlas-extension/1.0"})
    with urllib.request.urlopen(request, context=ssl_context(), timeout=120) as response:  # nosec B310
        return json.load(response)


def head(url: str) -> dict[str, str]:
    request = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": "sleep-gwas-atlas-extension/1.0"}
    )
    with urllib.request.urlopen(request, context=ssl_context(), timeout=120) as response:  # nosec B310
        return {key.lower(): value for key, value in response.headers.items()}


def catalog_query(endpoint: str, term: str) -> tuple[str, str, list[dict[str, object]]]:
    parameter = "efoTrait" if endpoint == "findByEfoTrait" else "diseaseTrait"
    url = f"{GWAS_API}/{endpoint}?{urllib.parse.urlencode({parameter: term, 'size': 500})}"
    payload = get_json(url)
    assert isinstance(payload, dict)
    studies = payload.get("_embedded", {}).get("studies", [])
    return endpoint, url, studies


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/replication/replication_source_queue.tsv"))
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/replication_finngen_r13_sources.tsv"))
    parser.add_argument("--mvp-sources", type=Path, default=Path("discovery_extension/config/replication_mvp_sources.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/replication/replication_unique_phenotype_search.tsv"))
    parser.add_argument("--catalog-candidates-out", type=Path, default=Path("discovery_extension/results/replication/replication_gwas_catalog_candidates.tsv"))
    parser.add_argument("--source-snapshot", type=Path, default=Path("discovery_extension/provenance/replication_sources/finngen_r13_source_snapshot.tsv"))
    parser.add_argument("--mvp-source-snapshot", type=Path, default=Path("discovery_extension/provenance/replication_sources/mvp_source_snapshot.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/replication_sources/search_provenance.json"))
    parser.add_argument("--search-date", default=datetime.now(timezone.utc).date().isoformat())
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    queue = read_tsv(args.queue)
    sources = read_tsv(args.sources)
    mvp_sources = read_tsv(args.mvp_sources)
    source_by_trait = {row["extension_trait_id"]: row for row in sources}
    mvp_by_trait = {row["extension_trait_id"]: row for row in mvp_sources}
    if set(source_by_trait) & set(mvp_by_trait):
        raise SystemExit("ERROR: a trait has more than one selected independent source")
    unique: dict[str, str] = {}
    for row in queue:
        unique.setdefault(row["extension_trait_id"], row["external_phenotype_name"])
    if len(queue) != 217 or len(unique) != 44:
        raise SystemExit(f"ERROR: expected frozen 217-pair/44-phenotype family, found {len(queue)}/{len(unique)}")
    if set(unique) != set(SEARCH_TERMS):
        raise SystemExit(f"ERROR: search-term coverage differs: missing={sorted(set(unique)-set(SEARCH_TERMS))} extra={sorted(set(SEARCH_TERMS)-set(unique))}")
    if args.validate_only:
        existing = read_tsv(args.out)
        if [row["extension_trait_id"] for row in existing] != list(unique) or any(row["search_status"] != "COMPLETE_BEFORE_RESULTS" for row in existing):
            raise SystemExit("ERROR: retained replication source search is incomplete or reordered")
        print(f"REPLICATION_SOURCE_SEARCH_VALID phenotypes={len(existing)} sha256={sha256(args.out)}")
        return

    phenos_payload = get_json(FINNGEN_API)
    if not isinstance(phenos_payload, list):
        raise SystemExit("ERROR: FinnGen R13 phenotype API did not return a list")
    phenos = {str(row["phenocode"]): row for row in phenos_payload}
    jobs = [(trait_id, endpoint, term) for trait_id in unique for term in SEARCH_TERMS[trait_id] for endpoint in ("findByEfoTrait", "findByDiseaseTrait")]
    catalog: dict[str, list[tuple[str, str, str, list[dict[str, object]]]]] = {trait_id: [] for trait_id in unique}
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(catalog_query, endpoint, term): (trait_id, endpoint, term) for trait_id, endpoint, term in jobs}
        for future in as_completed(futures):
            trait_id, _, term = futures[future]
            endpoint, url, studies = future.result()
            catalog[trait_id].append((endpoint, term, url, studies))

    snapshot_rows: list[dict[str, str]] = []
    for source in sources:
        code = source["finngen_phenocode"]
        if code not in phenos:
            raise SystemExit(f"ERROR: selected FinnGen endpoint absent from R13 API: {code}")
        pheno = phenos[code]
        if (str(pheno["phenostring"]) != source["replication_phenotype_definition"] or str(pheno["num_cases"]) != source["cases"] or str(pheno["num_controls"]) != source["controls"]):
            raise SystemExit(f"ERROR: selected FinnGen endpoint metadata changed: {code}")
        base_url = f"https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{code}.gz"
        locked_url = f"{base_url}?generation={source['source_generation']}"
        headers = head(locked_url)
        observed_etag = headers.get("etag", "").strip('"')
        if (headers.get("x-goog-generation") != source["source_generation"] or observed_etag != source["etag"] or headers.get("content-length") != source["content_length_bytes"]):
            raise SystemExit(f"ERROR: selected FinnGen remote identity changed: {code}")
        snapshot_rows.append({
            "replication_source_id": f"finngen_r13_{code}", "finngen_phenocode": code,
            "phenostring": str(pheno["phenostring"]), "cases": str(pheno["num_cases"]),
            "controls": str(pheno["num_controls"]), "source_url": locked_url,
            "source_generation": source["source_generation"], "etag": observed_etag,
            "md5_hex": source["md5_hex"], "content_length_bytes": headers["content-length"],
            "snapshot_status": "VERIFIED_BEFORE_RESULTS",
        })

    mvp_snapshot_rows: list[dict[str, str]] = []
    all_catalog_studies: dict[str, dict[str, object]] = {}
    for queried in catalog.values():
        for _, _, _, result in queried:
            for study in result:
                all_catalog_studies[str(study["accessionId"])] = study
    for source in mvp_sources:
        accession = source["study_accession"]
        study = all_catalog_studies.get(accession)
        if study is None or study.get("fullPvalueSet") is not True:
            raise SystemExit(f"ERROR: selected MVP catalog study missing or lacks full statistics: {accession}")
        if str((study.get("diseaseTrait") or {}).get("trait")) != source["replication_phenotype_definition"]:
            raise SystemExit(f"ERROR: selected MVP phenotype metadata changed: {accession}")
        headers = head(source["source_url"])
        observed_etag = headers.get("etag", "").strip('"')
        if observed_etag != source["etag"] or headers.get("content-length") != source["content_length_bytes"]:
            raise SystemExit(f"ERROR: selected MVP remote identity changed: {accession}")
        mvp_snapshot_rows.append({
            "replication_source_id": f"gwas_catalog_{accession}", "study_accession": accession,
            "phenotype": source["replication_phenotype_definition"], "cases": source["cases"],
            "controls": source["controls"], "source_url": source["source_url"],
            "etag": observed_etag, "md5_hex": source["md5_hex"],
            "content_length_bytes": headers["content-length"],
            "snapshot_status": "VERIFIED_BEFORE_RESULTS",
        })

    fields = [
        "extension_trait_id", "external_phenotype_name", "search_terms", "gwas_catalog_query_urls",
        "gwas_catalog_unique_study_count", "gwas_catalog_full_pvalue_study_count",
        "gwas_catalog_full_pvalue_accessions", "finngen_r13_selected_source_id",
        "selected_independent_source_id",
        "curation_decision", "search_date", "search_status",
    ]
    output: list[dict[str, str]] = []
    candidate_rows: list[dict[str, str]] = []
    for trait_id, phenotype_name in unique.items():
        queried = sorted(catalog[trait_id], key=lambda value: (value[0], value[1]))
        studies: dict[str, dict[str, object]] = {}
        urls: list[str] = []
        for _, _, url, result in queried:
            urls.append(url)
            for study in result:
                studies[str(study["accessionId"])] = study
        full = sorted(accession for accession, study in studies.items() if study.get("fullPvalueSet") is True)
        for accession, study in sorted(studies.items()):
            publication = study.get("publicationInfo") or {}
            ancestries = study.get("ancestries") or []
            ancestry_groups = sorted({
                str(group.get("ancestralGroup", ""))
                for ancestry in ancestries for group in (ancestry.get("ancestralGroups") or [])
                if group.get("ancestralGroup")
            })
            recruitment_countries = sorted({
                str(country.get("countryName", ""))
                for ancestry in ancestries for country in (ancestry.get("countryOfRecruitment") or [])
                if country.get("countryName")
            })
            candidate_rows.append({
                "extension_trait_id": trait_id, "external_phenotype_name": phenotype_name,
                "study_accession": accession,
                "catalog_disease_trait": one_line((study.get("diseaseTrait") or {}).get("trait", "NA")),
                "full_pvalue_set": "YES" if study.get("fullPvalueSet") is True else "NO",
                "initial_sample_size": one_line(study.get("initialSampleSize")),
                "ancestral_groups": ",".join(ancestry_groups) if ancestry_groups else "NA",
                "recruitment_countries": ",".join(recruitment_countries) if recruitment_countries else "NA",
                "pubmed_id": one_line(publication.get("pubmedId")),
                "publication_title": one_line(publication.get("title")),
                "study_url": str((study.get("_links") or {}).get("self", {}).get("href", "NA")),
            })
        selected = source_by_trait.get(trait_id)
        selected_mvp = mvp_by_trait.get(trait_id)
        selected_id = (
            f"finngen_r13_{selected['finngen_phenocode']}" if selected else
            f"gwas_catalog_{selected_mvp['study_accession']}" if selected_mvp else "NONE"
        )
        output.append({
            "extension_trait_id": trait_id, "external_phenotype_name": phenotype_name,
            "search_terms": " | ".join(SEARCH_TERMS[trait_id]),
            "gwas_catalog_query_urls": " | ".join(urls),
            "gwas_catalog_unique_study_count": str(len(studies)),
            "gwas_catalog_full_pvalue_study_count": str(len(full)),
            "gwas_catalog_full_pvalue_accessions": ",".join(full) if full else "NONE",
            "finngen_r13_selected_source_id": f"finngen_r13_{selected['finngen_phenocode']}" if selected else "NONE",
            "selected_independent_source_id": selected_id,
            "curation_decision": (
                "FINNGEN_R13_SELECTED_FOR_NONOVERLAPPING_SLEEP_TRAITS" if selected else
                "MVP_EUR_SELECTED" if selected_mvp else "NO_ELIGIBLE_INDEPENDENT_SOURCE_IDENTIFIED"
            ),
            "search_date": args.search_date, "search_status": "COMPLETE_BEFORE_RESULTS",
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    candidate_fields = [
        "extension_trait_id", "external_phenotype_name", "study_accession", "catalog_disease_trait",
        "full_pvalue_set", "initial_sample_size", "ancestral_groups", "recruitment_countries",
        "pubmed_id", "publication_title", "study_url",
    ]
    with args.catalog_candidates_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=candidate_fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(candidate_rows)
    args.source_snapshot.parent.mkdir(parents=True, exist_ok=True)
    with args.source_snapshot.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(snapshot_rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(snapshot_rows)
    with args.mvp_source_snapshot.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(mvp_snapshot_rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(mvp_snapshot_rows)
    provenance = {
        "schema_version": "1.0.0", "completed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "before_replication_result_access", "replication_results_accessed": False,
        "candidate_pair_count": len(queue), "unique_phenotype_count": len(unique),
        "gwas_catalog_query_count": len(jobs), "gwas_catalog_api": GWAS_API,
        "finngen_phenotype_api": FINNGEN_API, "finngen_release": "R13",
        "queue_sha256": sha256(args.queue), "source_config_sha256": sha256(args.sources),
        "mvp_source_config_sha256": sha256(args.mvp_sources),
        "search_output_sha256": sha256(args.out), "source_snapshot_sha256": sha256(args.source_snapshot),
        "mvp_source_snapshot_sha256": sha256(args.mvp_source_snapshot),
        "catalog_candidate_count": len(candidate_rows),
        "catalog_candidates_sha256": sha256(args.catalog_candidates_out),
        "search_scope_limit": "Exact GWAS Catalog disease-trait/EFO-label queries plus documented synonyms and FinnGen R13 endpoint review; absence is not proof that no controlled-access or unindexed source exists.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"REPLICATION_SOURCE_SEARCH_COMPLETE phenotypes={len(output)} catalog_queries={len(jobs)} selected_sources={len(sources) + len(mvp_sources)} sha256={sha256(args.out)}")


if __name__ == "__main__":
    main()
