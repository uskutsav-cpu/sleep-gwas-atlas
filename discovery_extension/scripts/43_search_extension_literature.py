#!/usr/bin/env python3
"""Run a resumable four-lens Europe PMC search for every extension-FDR pair."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import ssl
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
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
LENSES = {
    "DIRECT_RG": '("genetic correlation" OR "genome-wide genetic correlation")',
    "LDSC": '("LD score regression" OR "LDSC")',
    "RELATED_GENETIC": '("GWAS" OR "genome-wide association" OR "pleiotropy" OR "colocalization" OR "shared loci")',
    "MR_OBSERVATIONAL": '("Mendelian randomization" OR "Mendelian randomisation" OR "observational")',
}
GENERIC_NAMES = {
    "operation code", "operative procedures - main opcs4",
    "non-cancer illness code, self-reported", "treatment/medication code",
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "NA") for field in fields})
    return buffer.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def normalize(value: object) -> str:
    if value in (None, ""):
        return "NA"
    return " ".join(str(value).split()).replace("\t", " ")


def search_label(panel_row: dict[str, str]) -> str:
    name = panel_row["phenotype_name"].strip()
    if name.lower() in GENERIC_NAMES:
        parts = [part.strip() for part in panel_row["phenotype_definition"].split("|") if part.strip()]
        if len(parts) > 1:
            name = parts[1]
    name = re.sub(r"^(?:FI\d+\s*:\s*synonym|Medication for pain relief, inflammation, or fever)\s*[:|]?\s*", "", name, flags=re.I)
    name = re.sub(r"\s+", " ", name).strip(" .;:")
    return name or panel_row["phenotype_name"]


def quoted_phrase(value: str) -> str:
    clean = value.replace('"', " ").replace("/", " ")
    clean = " ".join(clean.split())
    return f'"{clean}"'


def search(query: str, page_size: int, attempts: int = 5) -> tuple[str, dict[str, object]]:
    url = API + "?" + urlencode({
        "query": query, "format": "json", "resultType": "core",
        "pageSize": page_size, "sort": "CITED desc",
    })
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": "sleep-gwas-atlas-extension-audit/1.0"})
            with urlopen(request, timeout=45, context=TLS_CONTEXT) as response:
                return url, json.load(response)
        except Exception as error:
            last = error
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Europe PMC query failed after {attempts} attempts: {last}")


def execute_query(row: dict[str, str], page_size: int) -> tuple[dict[str, object], list[dict[str, object]]]:
    url, payload = search(row["europe_pmc_query"], page_size)
    result_list = payload.get("resultList", {}).get("result", [])
    query_row: dict[str, object] = dict(row)
    query_row.update({
        "europe_pmc_url": url,
        "hit_count": int(payload.get("hitCount", 0)),
        "returned_result_count": len(result_list),
        "query_status": "COMPLETE",
    })
    candidates = []
    for rank, result in enumerate(result_list, 1):
        pmid = normalize(result.get("pmid"))
        candidates.append({
            "query_id": row["query_id"], "pair_id": row["pair_id"],
            "sleep_trait": row["sleep_trait"], "extension_trait_id": row["extension_trait_id"],
            "phenotype_name": row["phenotype_name"], "evidence_lens": row["evidence_lens"],
            "result_rank": rank, "pmid": pmid, "pmcid": normalize(result.get("pmcid")),
            "doi": normalize(result.get("doi")), "title": normalize(result.get("title")),
            "abstract": normalize(result.get("abstractText")),
            "author_string": normalize(result.get("authorString")),
            "journal": normalize(result.get("journalTitle")),
            "publication_year": normalize(result.get("pubYear")),
            "cited_by_count": normalize(result.get("citedByCount")),
            "publication_type": normalize(result.get("pubType")),
            "is_open_access": normalize(result.get("isOpenAccess")),
            "europe_pmc_record_url": f"https://europepmc.org/article/MED/{pmid}" if pmid != "NA" else "NA",
        })
    return query_row, candidates


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--as-of-date", default="2026-08-29")
    parser.add_argument("--page-size", type=int, default=5)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    try:
        date.fromisoformat(args.as_of_date)
    except ValueError:
        fail("--as-of-date must use YYYY-MM-DD")
    if not 1 <= args.page_size <= 10 or not 1 <= args.workers <= 8:
        fail("page size must be 1..10 and workers 1..8")
    root = Path(args.root).resolve()
    audit_path = root / "discovery_extension/results/novelty/extension_novelty_audit.tsv"
    panel_path = root / "discovery_extension/config/candidate_traits.tsv"
    query_path = root / "discovery_extension/results/novelty/literature_search_queries.tsv"
    candidate_path = root / "discovery_extension/results/novelty/literature_search_candidates.tsv"
    provenance_path = root / "discovery_extension/provenance/extension_literature_search.json"
    audit = read_tsv(audit_path)
    if len(audit) != 603 or len({row["pair_id"] for row in audit}) != 603:
        fail(f"expected 603 unique extension-FDR pairs, observed {len(audit)}")
    panel = {row["extension_trait_id"]: row for row in read_tsv(panel_path)}
    queries: list[dict[str, str]] = []
    for pair in audit:
        sleep = pair["sleep_trait"]
        trait_id = pair["extension_trait_id"]
        if sleep not in SLEEP_TERMS or trait_id not in panel:
            fail(f"missing controlled search term or panel row: {pair['pair_id']}")
        external = search_label(panel[trait_id])
        for lens, evidence in LENSES.items():
            expression = (
                f"TITLE_ABS:({SLEEP_TERMS[sleep]}) AND "
                f"TITLE_ABS:({quoted_phrase(external)}) AND "
                f"TITLE_ABS:({evidence}) AND "
                f"FIRST_PDATE:[1900-01-01 TO {args.as_of_date}]"
            )
            queries.append({
                "query_id": f"{pair['pair_id']}__{lens}", "pair_id": pair["pair_id"],
                "sleep_trait": sleep, "extension_trait_id": trait_id,
                "phenotype_name": pair["phenotype_name"], "external_search_label": external,
                "evidence_lens": lens, "as_of_date": args.as_of_date,
                "europe_pmc_query": expression,
            })
    expected = {row["query_id"]: row for row in queries}
    if len(expected) != 603 * len(LENSES):
        fail("query IDs are not unique or the four-lens family is incomplete")
    existing_queries = read_tsv(query_path) if query_path.is_file() else []
    existing_candidates = read_tsv(candidate_path) if candidate_path.is_file() else []
    complete: dict[str, dict[str, object]] = {}
    candidates_by_query: dict[str, list[dict[str, object]]] = {}
    for row in existing_queries:
        query_id = row.get("query_id", "")
        if query_id in expected and row.get("query_status") == "COMPLETE" and row.get("europe_pmc_query") == expected[query_id]["europe_pmc_query"]:
            complete[query_id] = row
    for row in existing_candidates:
        if row.get("query_id") in complete:
            candidates_by_query.setdefault(row["query_id"], []).append(row)
    if args.validate_only:
        if set(complete) != set(expected):
            fail(f"literature search is incomplete: {len(complete)}/{len(expected)}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if provenance.get("query_count") != len(expected) or provenance.get("pair_count") != 603:
            fail("literature search provenance counts drifted")
        print(f"EXTENSION_LITERATURE_SEARCH_VALID queries={len(complete)} candidates={len(existing_candidates)}")
        return 0

    query_fields = list(queries[0]) + ["europe_pmc_url", "hit_count", "returned_result_count", "query_status"]
    candidate_fields = [
        "query_id", "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name",
        "evidence_lens", "result_rank", "pmid", "pmcid", "doi", "title", "abstract",
        "author_string", "journal", "publication_year", "cited_by_count", "publication_type",
        "is_open_access", "europe_pmc_record_url",
    ]

    def checkpoint() -> None:
        query_rows = sorted(complete.values(), key=lambda row: str(row["query_id"]))
        candidate_rows = sorted(
            (row for query_id in complete for row in candidates_by_query.get(query_id, [])),
            key=lambda row: (str(row["query_id"]), int(row["result_rank"])),
        )
        atomic_text(query_path, tsv_text(query_rows, query_fields))
        atomic_text(candidate_path, tsv_text(candidate_rows, candidate_fields))

    pending = [row for row in queries if row["query_id"] not in complete]
    failures: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(execute_query, row, args.page_size): row for row in pending}
        for future in as_completed(futures):
            source = futures[future]
            try:
                query_row, candidate_rows = future.result()
                complete[source["query_id"]] = query_row
                candidates_by_query[source["query_id"]] = candidate_rows
            except Exception as error:
                failures.append((source["query_id"], str(error)))
            observed = len(complete) + len(failures)
            if observed % 100 == 0 or observed == len(expected):
                checkpoint()
                print(f"EXTENSION_LITERATURE_SEARCH_PROGRESS complete={len(complete)} failed={len(failures)} total={len(expected)}", flush=True)
    checkpoint()
    if failures:
        fail("literature search failed closed after retaining resumable progress: " + "; ".join(f"{key}: {value}" for key, value in failures[:5]))
    completed_rows = sorted(complete.values(), key=lambda row: str(row["query_id"]))
    candidate_rows = sorted(
        (row for query_id in complete for row in candidates_by_query.get(query_id, [])),
        key=lambda row: (str(row["query_id"]), int(row["result_rank"])),
    )
    query_payload = tsv_text(completed_rows, query_fields)
    candidate_payload = tsv_text(candidate_rows, candidate_fields)
    provenance = {
        "schema_version": "sleep-atlas-extension-literature-search.1",
        "search_provider": "Europe PMC REST API", "api_endpoint": API,
        "as_of_date": args.as_of_date, "pair_count": 603,
        "evidence_lenses_per_pair": list(LENSES), "query_count": len(completed_rows),
        "candidate_row_count": len(candidate_rows), "page_size_per_query": args.page_size,
        "sorting": "CITED desc", "source_audit": str(audit_path.relative_to(root)),
        "source_audit_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        "query_snapshot_sha256": hashlib.sha256(query_payload.encode()).hexdigest(),
        "candidate_snapshot_sha256": hashlib.sha256(candidate_payload.encode()).hexdigest(),
        "interpretation_limit": "Search candidates require pair-level evidence review; hit counts do not establish novelty or absence of evidence.",
    }
    atomic_text(query_path, query_payload)
    atomic_text(candidate_path, candidate_payload)
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"EXTENSION_LITERATURE_SEARCH_BUILT queries={len(completed_rows)} candidates={len(candidate_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
