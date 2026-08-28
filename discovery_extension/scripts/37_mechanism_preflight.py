#!/usr/bin/env python3
"""Validate the mechanistic source/search contract without treating landing pages as evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_SOURCES = {
    "GTEX_PORTAL", "EQTL_CATALOGUE", "OPEN_TARGETS_PLATFORM", "PSYCHENCODE_SYNAPSE",
    "ENCODE_DCC", "SCREEN_CCRE", "SINGLE_CELL_EXPRESSION_ATLAS", "GWAS_CATALOG",
    "REACTOME", "BIOSTUDIES_GEO",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/mechanistic_sources.tsv"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/mechanistic_annotation_contract.json"))
    parser.add_argument("--fine-mapping", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_colocalization.tsv"))
    parser.add_argument("--fine-mapping-lock", type=Path, default=Path("discovery_extension/config/fine_mapping_manifest.lock.json"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/mechanism/mechanism_readiness.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/mechanism_dependencies.json"))
    args = parser.parse_args()

    sources = read_tsv(args.sources)
    if {row["source_id"] for row in sources} != EXPECTED_SOURCES or len(sources) != len(EXPECTED_SOURCES):
        raise SystemExit("ERROR: mechanistic source registry is incomplete or duplicated")
    for row in sources:
        if not row["official_url"].startswith("https://"):
            raise SystemExit(f"ERROR: mechanistic source is not an HTTPS official landing page: {row['source_id']}")
        if row["exact_release_status"] != "RESULT_DEPENDENT_NOT_LOCKED" or not row["claim_limit"].strip():
            raise SystemExit(f"ERROR: mechanistic source improperly claims an exact release or lacks a claim limit: {row['source_id']}")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    tasks = contract["required_search_tasks_per_signal"]
    if len(tasks) != 7 or len(set(tasks)) != len(tasks):
        raise SystemExit("ERROR: mechanistic search task family must contain seven unique tasks")
    if contract["source_policy"]["registry"] != "config/mechanistic_sources.tsv":
        raise SystemExit("ERROR: mechanistic contract points to the wrong source registry")
    fine_mapping_ready = args.fine_mapping.is_file() and args.fine_mapping_lock.is_file()
    readiness = "READY_FOR_RESULT_FREE_SOURCE_PLAN" if fine_mapping_ready else "PROTOCOL_READY_UPSTREAM_FINE_MAPPING_BLOCKED"
    output = [{
        "source_id": row["source_id"], "domain": row["domain"], "official_url": row["official_url"],
        "landing_page_status": "VERIFIED_REGISTRY_ENTRY", "exact_release_status": row["exact_release_status"],
        "protocol_status": "PASS", "upstream_fine_mapping_status": "PRESENT" if fine_mapping_ready else "BLOCKED_UPSTREAM",
        "readiness": readiness,
    } for row in sources]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(output[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall_status": readiness,
        "source_count": len(sources),
        "search_task_count": len(tasks),
        "sources_sha256": sha256(args.sources),
        "contract_sha256": sha256(args.contract),
        "fine_mapping_present": args.fine_mapping.is_file(),
        "fine_mapping_lock_present": args.fine_mapping_lock.is_file(),
        "evidence_search_started": False,
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "warning": "Verified landing-page registry entries are search routes only. Exact releases, accessions, files, snapshots, citations, and results remain unlocked until a real signal family exists.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"MECHANISM_PREFLIGHT_{readiness} sources={len(sources)} fine_mapping_ready={str(fine_mapping_ready).lower()}")


if __name__ == "__main__":
    main()
