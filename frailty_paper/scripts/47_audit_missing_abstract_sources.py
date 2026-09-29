#!/usr/bin/env python3
"""Check whether empty screening abstracts are also absent from retained PubMed XML."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def audit(review_dir: Path, script_path: Path, external_storage_root: Path | None = None) -> dict:
    review_root = review_dir.resolve()
    queue_path = review_root / "screening/title_abstract_queue.tsv"
    source_manifest = review_root / "pubmed_source_files.tsv"
    queue_fields, queue_rows = read_tsv(queue_path)
    required_queue = {"pmid", "abstract"}
    if not required_queue.issubset(queue_fields):
        raise ValueError(f"queue missing columns: {sorted(required_queue - set(queue_fields))}")
    empty_rows = [row for row in queue_rows if not (row.get("abstract") or "").strip()]
    missing_ids = {row["pmid"].strip() for row in empty_rows if row.get("pmid", "").strip()}
    if len(missing_ids) != len([row for row in empty_rows if row.get("pmid", "").strip()]):
        raise ValueError("duplicate PMID among empty-abstract screening rows")

    manifest_fields, source_rows = read_tsv(source_manifest)
    required_manifest = {"source_file", "bytes", "sha256"}
    if not required_manifest.issubset(manifest_fields):
        raise ValueError(f"source manifest missing columns: {sorted(required_manifest - set(manifest_fields))}")

    occurrences: Counter[str] = Counter()
    allowed_roots = [review_root]
    if external_storage_root is not None:
        allowed_roots.append(external_storage_root.resolve())
    abstract_ids: set[str] = set()
    records_seen = 0
    for row in source_rows:
        relative = Path(row["source_file"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"unsafe source path in manifest: {relative}")
        path = (review_root / relative).resolve()
        if not any(path == allowed or allowed in path.parents for allowed in allowed_roots):
            raise ValueError(f"source XML resolves outside review and explicitly allowed external roots: {relative}")
        if not path.is_file():
            raise ValueError(f"missing source XML: {relative}")
        if path.stat().st_size != int(row["bytes"]):
            raise ValueError(f"source XML size mismatch: {relative}")
        if sha256(path) != row["sha256"].lower():
            raise ValueError(f"source XML SHA-256 mismatch: {relative}")
        for _, node in ET.iterparse(path, events=("end",)):
            if node.tag not in {"PubmedArticle", "PubmedBookArticle"}:
                continue
            records_seen += 1
            pmid_node = node.find(".//PMID")
            pmid = (pmid_node.text or "").strip() if pmid_node is not None else ""
            if pmid in missing_ids:
                occurrences[pmid] += 1
                if any("".join(item.itertext()).strip() for item in node.findall(".//Abstract/AbstractText")):
                    abstract_ids.add(pmid)
            node.clear()

    unmatched = sorted(missing_ids - set(occurrences))
    duplicate_ids = sorted(pmid for pmid, count in occurrences.items() if count > 1)
    status = "PASS" if not unmatched and not abstract_ids else "FAIL"
    return {
        "schema_version": 1,
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "screening_queue": "review/screening/title_abstract_queue.tsv",
        "screening_queue_sha256": sha256(queue_path),
        "source_manifest": "review/pubmed_source_files.tsv",
        "source_manifest_sha256": sha256(source_manifest),
        "allowed_external_storage_root": str(external_storage_root.resolve()) if external_storage_root else None,
        "source_xml_files_verified": len(source_rows),
        "source_xml_records_parsed": records_seen,
        "empty_abstract_queue_rows": len(empty_rows),
        "empty_abstract_pmids": len(missing_ids),
        "empty_abstract_pmids_found_in_source_xml": len(occurrences),
        "empty_abstract_pmids_absent_from_source_xml": unmatched,
        "pmids_with_duplicate_source_occurrences": len(duplicate_ids),
        "duplicate_source_pmids": duplicate_ids,
        "pmids_with_abstract_in_any_source_occurrence": sorted(abstract_ids),
        "script_sha256": sha256(script_path),
        "interpretation": (
            "Every empty-abstract PMID is present in verified retained XML, and no source occurrence "
            "contains an AbstractText element. This rules out abstract loss during retained-record "
            "priority selection; missing abstracts remain missing information, not exclusion evidence."
            if status == "PASS"
            else "Investigate source/query identity or builder field loss before using this audit to describe missing abstracts."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=".")
    parser.add_argument("--review-dir", default="frailty_paper/review")
    parser.add_argument("--external-storage-root", help="Allow source XML files resolving under this mounted data root")
    parser.add_argument("--output", default="frailty_paper/analysis/missing_abstract_source_audit_2026-09-24.json")
    args = parser.parse_args()
    root = Path(args.repo).resolve()
    review_dir = Path(args.review_dir)
    if not review_dir.is_absolute():
        review_dir = root / review_dir
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    external_root = Path(args.external_storage_root) if args.external_storage_root else None
    result = audit(review_dir, Path(__file__).resolve(), external_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(
        f"{result['status']} empty_abstract_pmids={result['empty_abstract_pmids']} "
        f"source_xml_found={result['empty_abstract_pmids_found_in_source_xml']} "
        f"with_abstract={len(result['pmids_with_abstract_in_any_source_occurrence'])} "
        f"duplicate_pmids={result['pmids_with_duplicate_source_occurrences']}"
    )
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
