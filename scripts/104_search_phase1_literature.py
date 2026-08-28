#!/usr/bin/env python3
"""Create a fixed-date Europe PMC search snapshot for all 153 Phase-1 discoveries."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import ssl
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PRIMARY = "PRIMARY_PHASE1"
SLEEP_TERMS = {
    "insomnia": '"insomnia"',
    "sleepdur": '"sleep duration"',
    "shortsleep": '"short sleep"',
    "longsleep": '"long sleep"',
    "chronotype": '("chronotype" OR "morningness")',
    "sleepiness": '"daytime sleepiness"',
    "napping": '"daytime napping"',
    "snoring": '"snoring"',
    "sleep_apnea": '("sleep apnea" OR "sleep apnoea" OR "obstructive sleep apnea")',
    "sleep_efficiency": '("sleep efficiency" AND (actigraphy OR accelerometer))',
    "accel_sleep_duration": '("sleep duration" AND (actigraphy OR accelerometer))',
    "sleep_timing": '(("sleep timing" OR "sleep midpoint") AND (actigraphy OR accelerometer))',
}
EXTERNAL_TERMS = {
    "mdd": '("major depressive disorder" OR "major depression")',
    "scz": '"schizophrenia"',
    "bipolar": '"bipolar disorder"',
    "adhd": '("attention deficit hyperactivity disorder" OR "ADHD")',
    "ibd": '"inflammatory bowel disease"',
    "crohn": '("Crohn disease" OR "Crohn\'s disease")',
    "ra": '"rheumatoid arthritis"',
    "asthma": '"asthma"',
    "bmi": '("body mass index" OR "BMI")',
    "ldl": '"LDL cholesterol"',
    "hdl": '"HDL cholesterol"',
    "triglycerides": '"triglycerides"',
    "cad": '("coronary artery disease" OR "coronary heart disease")',
    "stroke": '"stroke"',
    "atrial_fibrillation": '"atrial fibrillation"',
    "sbp": '("systolic blood pressure" OR "SBP")',
    "longevity": '"longevity"',
    "parental_lifespan": '("parental lifespan" OR "parents age at death")',
    "healthspan": '"healthspan"',
    "frailty": '("frailty" OR "frailty index")',
    "telomere_length": '("telomere length" OR "leukocyte telomere length")',
    "grip_strength": '("grip strength" OR "hand grip strength")',
    "breast_cancer": '"breast cancer"',
    "prostate_cancer": '"prostate cancer"',
    "colorectal_cancer": '("colorectal cancer" OR "bowel cancer")',
    "lung_cancer": '"lung cancer"',
    "ovarian_cancer": '"ovarian cancer"',
    "parkinson": '("Parkinson disease" OR "Parkinson\'s disease")',
}
LENSES = {
    "DIRECT_RG": '("genetic correlation" OR "genome-wide genetic correlation")',
    "LDSC": '("LD score regression" OR "LDSC")',
    "RELATED_GENETIC": '("GWAS" OR "genome-wide association" OR "pleiotropy" OR "colocalization" OR "shared loci")',
    "MR_OBSERVATIONAL": '("Mendelian randomization" OR "Mendelian randomisation" OR "observational")',
}


def tls_context() -> ssl.SSLContext:
    candidates = [os.environ.get("SSL_CERT_FILE"), "/etc/ssl/cert.pem"]
    try:
        import certifi
        candidates.append(certifi.where())
    except ImportError:
        pass
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return ssl.create_default_context(cafile=candidate)
    return ssl.create_default_context()


TLS_CONTEXT = tls_context()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "NA") for field in fields})
    return buffer.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def search(query: str, page_size: int, attempts: int = 4) -> dict[str, object]:
    url = API + "?" + urlencode({
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": page_size,
        "sort": "CITED desc",
    })
    last = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": "sleep-gwas-atlas-literature-audit/1.0"})
            with urlopen(request, timeout=30, context=TLS_CONTEXT) as response:
                payload = json.load(response)
            return {"url": url, "payload": payload}
        except Exception as error:  # network/API errors are retried and then fail closed
            last = error
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Europe PMC query failed after {attempts} attempts: {last}")


def normalize(value: object) -> str:
    if value in (None, ""):
        return "NA"
    return " ".join(str(value).split()).replace("\t", " ")


def execute_query(row: dict[str, str], page_size: int) -> tuple[dict[str, object], list[dict[str, object]]]:
    response = search(row["europe_pmc_query"], page_size)
    payload = response["payload"]
    result_list = payload.get("resultList", {}).get("result", [])
    query_row: dict[str, object] = dict(row)
    query_row.update({
        "europe_pmc_url": response["url"],
        "hit_count": int(payload.get("hitCount", 0)),
        "returned_result_count": len(result_list),
        "query_status": "COMPLETE",
    })
    candidates = []
    for rank, result in enumerate(result_list, 1):
        candidates.append({
            "query_id": row["query_id"],
            "sleep_trait": row["sleep_trait"],
            "external_trait": row["external_trait"],
            "evidence_lens": row["evidence_lens"],
            "result_rank": rank,
            "pmid": normalize(result.get("pmid")),
            "pmcid": normalize(result.get("pmcid")),
            "doi": normalize(result.get("doi")),
            "title": normalize(result.get("title")),
            "author_string": normalize(result.get("authorString")),
            "journal": normalize(result.get("journalTitle")),
            "publication_year": normalize(result.get("pubYear")),
            "cited_by_count": normalize(result.get("citedByCount")),
            "publication_type": normalize(result.get("pubType")),
            "is_open_access": normalize(result.get("isOpenAccess")),
            "europe_pmc_record_url": (
                f"https://europepmc.org/article/MED/{result['pmid']}"
                if result.get("pmid") else "NA"
            ),
        })
    return query_row, candidates


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--as-of-date", default="2026-08-28")
    parser.add_argument("--page-size", type=int, default=8)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    try:
        date.fromisoformat(args.as_of_date)
    except ValueError:
        fail("--as-of-date must use YYYY-MM-DD")
    if not 1 <= args.page_size <= 25 or not 1 <= args.workers <= 8:
        fail("page size must be 1..25 and workers 1..8")
    root = Path(args.root).resolve()
    master_path = root / "results/analysis/phase1_master_analysis.tsv"
    master = read_tsv(master_path)
    pairs = [
        row for row in master
        if row["primary_or_sensitivity"] == PRIMARY and float(row["fdr"]) <= 0.05
    ]
    if len(pairs) != 153:
        fail(f"expected 153 locked primary discoveries, observed {len(pairs)}")
    if len({(row["sleep_trait"], row["external_trait"]) for row in pairs}) != 153:
        fail("significant pair keys are not unique")
    queries: list[dict[str, str]] = []
    for pair in pairs:
        sleep = pair["sleep_trait"]
        external = pair["external_trait"]
        if sleep not in SLEEP_TERMS or external not in EXTERNAL_TERMS:
            fail(f"missing controlled literature term for {sleep}__{external}")
        for lens, evidence in LENSES.items():
            expression = (
                f"TITLE_ABS:({SLEEP_TERMS[sleep]}) AND "
                f"TITLE_ABS:({EXTERNAL_TERMS[external]}) AND "
                f"TITLE_ABS:({evidence}) AND "
                f"FIRST_PDATE:[1900-01-01 TO {args.as_of_date}]"
            )
            queries.append({
                "query_id": f"{sleep}__{external}__{lens}",
                "sleep_trait": sleep,
                "external_trait": external,
                "sleep_trait_label": pair["sleep_trait_label"],
                "external_trait_label": pair["external_trait_label"],
                "external_domain": pair["external_domain"],
                "evidence_lens": lens,
                "as_of_date": args.as_of_date,
                "europe_pmc_query": expression,
            })
    query_path = root / "results/analysis/literature_search_queries.tsv"
    candidate_path = root / "results/analysis/literature_search_candidates.tsv"
    provenance_path = root / "results/analysis/literature_search_provenance.json"
    if args.validate_only:
        observed_queries = read_tsv(query_path)
        observed_candidates = read_tsv(candidate_path)
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if len(observed_queries) != 612 or {row["query_id"] for row in observed_queries} != {row["query_id"] for row in queries}:
            fail("literature query snapshot does not cover exactly four queries for all 153 pairs")
        if provenance.get("query_count") != 612 or provenance.get("pair_count") != 153:
            fail("literature search provenance counts drifted")
        if any(row["query_status"] != "COMPLETE" for row in observed_queries):
            fail("literature query snapshot contains an incomplete query")
        print(f"PHASE1_LITERATURE_SEARCH_VALID queries={len(observed_queries)} candidates={len(observed_candidates)}")
        return 0
    completed: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(execute_query, row, args.page_size): row for row in queries}
        for future in as_completed(futures):
            try:
                query_row, result_rows = future.result()
                completed.append(query_row)
                candidates.extend(result_rows)
            except Exception as error:
                failures.append((futures[future]["query_id"], str(error)))
            if (len(completed) + len(failures)) % 50 == 0:
                print(
                    f"LITERATURE_SEARCH_PROGRESS complete={len(completed)} "
                    f"failed={len(failures)} total={len(queries)}",
                    flush=True,
                )
    if failures:
        fail("literature search failed closed: " + "; ".join(f"{key}: {value}" for key, value in failures[:5]))
    completed.sort(key=lambda row: row["query_id"])
    candidates.sort(key=lambda row: (row["query_id"], int(row["result_rank"])))
    query_fields = list(queries[0]) + ["europe_pmc_url", "hit_count", "returned_result_count", "query_status"]
    candidate_fields = [
        "query_id", "sleep_trait", "external_trait", "evidence_lens", "result_rank",
        "pmid", "pmcid", "doi", "title", "author_string", "journal",
        "publication_year", "cited_by_count", "publication_type", "is_open_access",
        "europe_pmc_record_url",
    ]
    query_payload = tsv_text(completed, query_fields)
    candidate_payload = tsv_text(candidates, candidate_fields)
    provenance = {
        "schema_version": "sleep-atlas-phase1-literature-search.1",
        "search_provider": "Europe PMC REST API",
        "api_endpoint": API,
        "as_of_date": args.as_of_date,
        "pair_count": 153,
        "evidence_lenses_per_pair": list(LENSES),
        "query_count": len(completed),
        "candidate_row_count": len(candidates),
        "page_size_per_query": args.page_size,
        "sorting": "CITED desc",
        "master_input_path": str(master_path.relative_to(root)),
        "master_input_sha256": hashlib.sha256(master_path.read_bytes()).hexdigest(),
        "query_snapshot_sha256": hashlib.sha256(query_payload.encode()).hexdigest(),
        "candidate_snapshot_sha256": hashlib.sha256(candidate_payload.encode()).hexdigest(),
        "interpretation_limit": "Search candidates require pair-level human review; hit counts do not establish novelty or absence of evidence.",
    }
    atomic_text(query_path, query_payload)
    atomic_text(candidate_path, candidate_payload)
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"PHASE1_LITERATURE_SEARCH_BUILT queries={len(completed)} candidates={len(candidates)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
