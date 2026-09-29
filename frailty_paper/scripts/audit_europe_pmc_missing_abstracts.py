#!/usr/bin/env python3
"""Audit public Europe PMC abstract availability for blank PubMed queue rows.

Abstract text is read transiently from the API response to calculate its
length. It is never written to disk, and this script does not modify screening
records or decisions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from html.parser import HTMLParser
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests


ENDPOINT = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
QUEUE = Path("frailty_paper/review/screening/title_abstract_queue.tsv")
OUT_TSV = Path("frailty_paper/review/europe_pmc_missing_abstract_audit.tsv")
OUT_JSON = Path("frailty_paper/analysis/europe_pmc_missing_abstract_audit_2026-09-24.json")
BATCH_SIZE = 50


class AbstractBodyParser(HTMLParser):
    """Collect abstract text outside structural heading elements."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_heading = False
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() in {"h1", "h2", "h3", "h4", "h5", "h6"}:
            self.in_heading = True

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() in {"h1", "h2", "h3", "h4", "h5", "h6"}:
            self.in_heading = False

    def handle_data(self, data: str) -> None:
        if not self.in_heading and data.strip():
            self.parts.append(data.strip())


def abstract_body(abstract: str) -> str:
    parser = AbstractBodyParser()
    parser.feed(abstract)
    return " ".join(parser.parts).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def get_queue_pmids(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if not {"pmid", "abstract"}.issubset(rows.fieldnames or []):
            raise ValueError("Queue lacks required pmid and abstract columns")
        ids = [row["pmid"].strip() for row in rows if not row["abstract"].strip()]
    if len(ids) != len(set(ids)):
        raise ValueError("Blank-abstract queue contains duplicate PMIDs")
    return ids


def fetch_batch(session: requests.Session, pmids: list[str]) -> list[dict[str, Any]]:
    query = " OR ".join(f"EXT_ID:{pmid}" for pmid in pmids)
    response = session.get(
        ENDPOINT,
        params={"query": query, "format": "json", "resultType": "core", "pageSize": len(pmids)},
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()
    if "resultList" not in payload or "hitCount" not in payload:
        raise ValueError("Europe PMC response lacks resultList or hitCount")
    return payload["resultList"].get("result", [])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=QUEUE)
    parser.add_argument("--out-tsv", type=Path, default=OUT_TSV)
    parser.add_argument("--out-json", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    pmids = get_queue_pmids(args.queue)
    by_pmid: dict[str, dict[str, Any]] = {}
    batches = [pmids[i : i + BATCH_SIZE] for i in range(0, len(pmids), BATCH_SIZE)]
    with requests.Session() as session:
        session.headers.update({"User-Agent": "frailty-review-metadata-audit/1.0 (Europe PMC public API)"})
        for i, batch in enumerate(batches, 1):
            result = None
            for attempt in range(4):
                try:
                    result = fetch_batch(session, batch)
                    break
                except (requests.RequestException, ValueError):
                    if attempt == 3:
                        raise
                    time.sleep(2**attempt)
            for record in result or []:
                pmid = str(record.get("pmid") or "").strip()
                if pmid not in batch:
                    continue
                abstract = record.get("abstractText") or ""
                body = abstract_body(abstract)
                by_pmid[pmid] = {
                    "pmid": pmid,
                    "europe_pmc_match": True,
                    "abstract_field_characters": len(abstract.strip()),
                    "abstract_body_available": bool(body),
                    "abstract_body_characters": len(body),
                    "record_source": record.get("source", ""),
                }
            print(f"batch={i}/{len(batches)} returned={len(result or [])} matched_total={len(by_pmid)}", flush=True)
            if i < len(batches):
                time.sleep(0.25)

    rows = [
        by_pmid.get(
            pmid,
            {
                "pmid": pmid,
                "europe_pmc_match": False,
                "abstract_field_characters": 0,
                "abstract_body_available": False,
                "abstract_body_characters": 0,
                "record_source": "",
            },
        )
        for pmid in pmids
    ]
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint": ENDPOINT,
        "query_field": "Europe PMC EXT_ID (PMID)",
        "result_type": "core",
        "batch_size": BATCH_SIZE,
        "queue_path": str(args.queue),
        "queue_sha256": sha256(args.queue),
        "blank_abstract_pmids": len(pmids),
        "europe_pmc_matched": sum(r["europe_pmc_match"] for r in rows),
        "abstracts_with_body": sum(r["abstract_body_available"] for r in rows),
        "matched_without_abstract_body": sum(r["europe_pmc_match"] and not r["abstract_body_available"] for r in rows),
        "heading_only_or_empty_abstract_fields": sum(r["europe_pmc_match"] and r["abstract_field_characters"] > 0 and not r["abstract_body_available"] for r in rows),
        "not_matched": sum(not r["europe_pmc_match"] for r in rows),
        "audit_tsv": str(args.out_tsv),
        "audit_tsv_sha256": sha256(args.out_tsv),
        "abstract_text_persisted": False,
        "screening_queue_modified": False,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
