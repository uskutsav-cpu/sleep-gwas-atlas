#!/usr/bin/env python3
"""Build a separate, unscreened index for a reconciled post-cutoff PubMed window.

This deliberately does not edit the frozen PubMed record corpus, screening
queue, or PRISMA counts. All occurrences remain available alongside a
window-deduplicated, blank-decision title/abstract review list.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BUILDER_PATH = ROOT / "frailty_paper/scripts/08_build_review_records.py"
SPEC = importlib.util.spec_from_file_location("review_record_builder", BUILDER_PATH)
assert SPEC and SPEC.loader
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)

OCCURRENCE_FIELDS = [
    "window_record_id", "record_type", "pmid", "doi", "normalized_title", "title",
    "abstract", "year", "journal", "query_id", "source_file", "source_database",
    "source_record_id", "authors", "source_xml_bytes", "source_xml_sha256",
    "window_duplicate_group", "window_retained_record_id", "window_duplicate_role",
    "window_duplicate_count", "window_duplicate_sources", "frozen_snapshot_match",
    "frozen_snapshot_match_rule", "screening_status",
]
QUEUE_FIELDS = [
    "screening_id", "pmid", "doi", "year", "title", "abstract", "query_id",
    "duplicate_count", "duplicate_sources", "frozen_snapshot_match",
    "frozen_snapshot_match_rule", "screening_status", "title_abstract_decision_r1",
    "title_abstract_reason_r1", "reviewer_1", "title_abstract_decision_r2",
    "title_abstract_reason_r2", "reviewer_2", "adjudication", "adjudicator", "notes",
]
QUERY_PRIORITY = {"primary_sleep_frailty": 0, "genetic_sleep_frailty": 1,
                  "frailty_neighborhood": 2}
csv.field_size_limit(100_000_000)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.partial")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            record = io.StringIO(newline="")
            csv.DictWriter(record, fieldnames=fields, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n").writerow(row)
            serialized = record.getvalue()
            # Keep empty final TSV cells parseable without ending the line in tabs.
            if serialized.endswith("\t\n"):
                serialized = serialized[:-1] + '""\n'
            stream.write(serialized)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def guard_unique_queue(path: Path) -> None:
    if not path.is_file():
        return
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not set(QUEUE_FIELDS).issubset(reader.fieldnames):
            raise RuntimeError(f"existing supplemental queue is malformed: {path}")
        decision_columns = {
            "title_abstract_decision_r1", "title_abstract_reason_r1", "reviewer_1",
            "title_abstract_decision_r2", "title_abstract_reason_r2", "reviewer_2",
            "adjudication", "adjudicator", "notes",
        }
        for row in reader:
            if any((row.get(key) or "").strip() for key in decision_columns):
                raise RuntimeError(
                    f"refusing to overwrite supplemental reviewer data in {path}"
                )


def normalized_doi(value: str) -> str:
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", (value or "").strip().lower())


def snapshot_indexes(path: Path) -> tuple[set[str], dict[str, set[str]], dict[str, set[str]]]:
    pmids: set[str] = set()
    dois: dict[str, set[str]] = defaultdict(set)
    titles: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"pmid", "doi", "normalized_title"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RuntimeError(f"frozen all-records table lacks columns {sorted(required)}")
        for row in reader:
            pmid = (row.get("pmid") or "").strip()
            doi = normalized_doi(row.get("doi") or "")
            title = (row.get("normalized_title") or "").strip()
            if pmid:
                pmids.add(pmid)
            if doi:
                dois[doi].add(title)
            if title:
                titles[title].add(pmid or doi)
    return pmids, dois, titles


def snapshot_match(record: dict, indexes: tuple[set[str], dict[str, set[str]], dict[str, set[str]]]) -> str:
    pmids, dois, titles = indexes
    pmid = (record.get("pmid") or "").strip()
    doi = normalized_doi(record.get("doi") or "")
    title = (record.get("normalized_title") or "").strip()
    if pmid:
        return "PMID" if pmid in pmids else ""
    if doi and doi in dois:
        return "DOI_AND_TITLE" if title in dois[doi] else "DOI_TITLE_CONFLICT"
    if title and title in titles:
        candidates = titles[title]
        return "TITLE_ONLY" if len(candidates) == 1 and not next(iter(candidates)) else "TITLE_CANDIDATE"
    return ""


def reconcile_inputs(increment_dir: Path, reconciliation_path: Path) -> tuple[dict, dict, list[Path]]:
    window_path = increment_dir / "window_manifest.json"
    if not window_path.is_file():
        raise RuntimeError(f"missing bounded-window manifest: {window_path}")
    if not reconciliation_path.is_file():
        raise RuntimeError(f"missing ESearch reconciliation report: {reconciliation_path}")
    window = json.loads(window_path.read_text(encoding="utf-8"))
    reconciliation = json.loads(reconciliation_path.read_text(encoding="utf-8"))
    if window.get("integration_status") != "SEPARATE_CACHE_NOT_MERGED":
        raise RuntimeError("expected a separate, unmerged bounded PubMed window")
    if reconciliation.get("query_errors") or reconciliation.get("missing_segment_id_occurrences") or reconciliation.get("extra_segment_id_occurrences"):
        raise RuntimeError("bounded-window ESearch reconciliation is not exact")
    if not window.get("queries"):
        raise RuntimeError("bounded-window manifest contains no query summaries")
    files = sorted(path for path in increment_dir.rglob("*.xml")
                   if path.is_file() and not path.name.startswith("._"))
    expected = {row["query_id"]: int(row["fetched"]) for row in window["queries"]}
    actual: Counter = Counter()
    for path in files:
        actual[path.relative_to(increment_dir).parts[0]] += sum(
            1 for _, element in ET.iterparse(path, events=("end",))
            if element.tag in {"PubmedArticle", "PubmedBookArticle"}
        )
    if dict(actual) != {key: value for key, value in expected.items() if value}:
        raise RuntimeError(f"bounded-window XML counts differ from its search log: expected {expected}, found {dict(actual)}")
    if sum(actual.values()) != int(window.get("records_returned_across_queries", -1)):
        raise RuntimeError("bounded-window XML total differs from window manifest")
    reconciled_nodes = int(reconciliation.get("pubmed_article_nodes", 0)) + int(reconciliation.get("pubmed_book_article_nodes", 0))
    if reconciled_nodes != sum(actual.values()):
        raise RuntimeError("ESearch reconciliation node count differs from XML files")
    return window, reconciliation, files


def build(increment_dir: Path, review_dir: Path, outdir: Path,
          reconciliation_path: Path, source_prefix: str) -> dict:
    increment_dir = increment_dir.resolve()
    review_dir = review_dir.resolve()
    outdir = outdir.resolve()
    occurrence_path = outdir / "pubmed_window_occurrences.tsv"
    queue_path = outdir / "pubmed_window_unique_review_queue.tsv"
    guard_unique_queue(queue_path)

    window, reconciliation, xml_files = reconcile_inputs(increment_dir, reconciliation_path)
    snapshot_output_manifest = review_dir / "record_build_outputs.tsv"
    snapshot_manifest_path = review_dir / "record_build_manifest.json"
    snapshot_records_path = review_dir / "all_records.tsv"
    if not all(path.is_file() for path in (snapshot_output_manifest, snapshot_manifest_path, snapshot_records_path)):
        raise RuntimeError("frozen PubMed snapshot records/manifests are incomplete")
    with snapshot_output_manifest.open(encoding="utf-8", newline="") as stream:
        output_rows = {row["file"]: row for row in csv.DictReader(stream, delimiter="\t")}
    expected_snapshot = output_rows.get("all_records.tsv")
    if not expected_snapshot or expected_snapshot["sha256"] != sha256(snapshot_records_path):
        raise RuntimeError("frozen PubMed all_records.tsv no longer matches its build output manifest")
    snapshot_build = json.loads(snapshot_manifest_path.read_text(encoding="utf-8"))
    with snapshot_records_path.open(encoding="utf-8", newline="") as stream:
        snapshot_row_count = sum(1 for _ in csv.DictReader(stream, delimiter="\t"))
    if snapshot_build.get("total_source_records") != snapshot_row_count:
        raise RuntimeError("frozen PubMed all-record count differs from its build manifest")
    indexes = snapshot_indexes(snapshot_records_path)

    occurrences: list[dict] = []
    file_rows: list[dict] = []
    for path in xml_files:
        source_rel = (Path(source_prefix) / path.relative_to(increment_dir)).as_posix()
        file_hash = sha256(path)
        n = 0
        try:
            for _, element in ET.iterparse(path, events=("end",)):
                if element.tag in {"PubmedArticle", "PubmedBookArticle"}:
                    n += 1
                    record = BUILDER.parse_article(
                        element, path, increment_dir, n, Path(source_prefix)
                    )
                    record["window_record_id"] = record.pop("record_id")
                    record["source_xml_bytes"] = path.stat().st_size
                    record["source_xml_sha256"] = file_hash
                    record["_snapshot_match"] = snapshot_match(record, indexes)
                    occurrences.append(record)
                    element.clear()
        except Exception as exc:
            raise RuntimeError(f"failed to parse supplemental PubMed XML {path}: {exc}") from exc
        if n and any(row["source_file"] != source_rel for row in occurrences[-n:]):
            raise RuntimeError(f"source locator mismatch for {path}")
        file_rows.append({"source_file": source_rel, "bytes": path.stat().st_size,
                          "sha256": file_hash, "record_count": n})

    occurrences.sort(key=lambda row: (QUERY_PRIORITY.get(row["query_id"], 9),
                                      row["source_file"], row["window_record_id"]))
    groups = BUILDER.group_duplicate_records(occurrences)
    group_by_record: dict[str, dict] = {}
    for group_index, (rule, members) in enumerate(
        sorted(groups, key=lambda group: min(BUILDER.dedup_keys(member) for member in group[1])), 1
    ):
        members.sort(key=lambda row: (QUERY_PRIORITY.get(row["query_id"], 9),
                                       row["source_file"], row["window_record_id"]))
        group_id = f"W20260925-D{group_index:04d}"
        retained = members[0]["window_record_id"]
        sources = ";".join(sorted({row["source_file"] for row in members}))
        for member in members:
            group_by_record[member["window_record_id"]] = {
                "duplicate_group": group_id,
                "retained": retained,
                "duplicate_count": len(members),
                "duplicate_sources": sources,
                "dedup_rule": rule,
            }

    snapshot_match_count = 0
    occurrence_rows = []
    unique_rows = []
    seen_groups = set()
    for row in occurrences:
        group = group_by_record[row["window_record_id"]]
        match_rule = row["_snapshot_match"]
        frozen_match = bool(match_rule and match_rule not in {"DOI_TITLE_CONFLICT", "TITLE_CANDIDATE"})
        snapshot_match_count += bool(frozen_match and row["window_record_id"] == group["retained"])
        status = "DUPLICATE_OF_FROZEN_SNAPSHOT_UNSCREENED" if frozen_match else "SUPPLEMENTAL_OUT_OF_WINDOW_UNSCREENED"
        occurrence_rows.append({
            **{key: row.get(key, "") for key in BUILDER.FIELDS if key != "record_id"},
            "window_record_id": row["window_record_id"],
            "source_xml_bytes": row["source_xml_bytes"],
            "source_xml_sha256": row["source_xml_sha256"],
            "window_duplicate_group": group["duplicate_group"],
            "window_retained_record_id": group["retained"],
            "window_duplicate_role": "RETAINED" if row["window_record_id"] == group["retained"] else "DUPLICATE_OCCURRENCE",
            "window_duplicate_count": group["duplicate_count"],
            "window_duplicate_sources": group["duplicate_sources"],
            "frozen_snapshot_match": "YES" if frozen_match else "NO",
            "frozen_snapshot_match_rule": match_rule,
            "screening_status": status,
        })
        if group["duplicate_group"] in seen_groups:
            continue
        seen_groups.add(group["duplicate_group"])
        unique_rows.append({
            "screening_id": group["retained"],
            "pmid": row.get("pmid", ""), "doi": row.get("doi", ""),
            "year": row.get("year", ""), "title": row.get("title", ""),
            "abstract": row.get("abstract", ""), "query_id": row.get("query_id", ""),
            "duplicate_count": group["duplicate_count"],
            "duplicate_sources": group["duplicate_sources"],
            "frozen_snapshot_match": "YES" if frozen_match else "NO",
            "frozen_snapshot_match_rule": match_rule,
            "screening_status": status,
            "title_abstract_decision_r1": "", "title_abstract_reason_r1": "", "reviewer_1": "",
            "title_abstract_decision_r2": "", "title_abstract_reason_r2": "", "reviewer_2": "",
            "adjudication": "", "adjudicator": "", "notes": "",
        })

    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / "pubmed_window_source_manifest.tsv",
              ["source_file", "bytes", "sha256", "record_count"], file_rows)
    write_tsv(occurrence_path, OCCURRENCE_FIELDS, occurrence_rows)
    write_tsv(queue_path, QUEUE_FIELDS, unique_rows)
    outputs = []
    for path in (outdir / "pubmed_window_source_manifest.tsv", occurrence_path, queue_path):
        outputs.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    output_manifest_path = outdir / "pubmed_window_review_outputs.tsv"
    write_tsv(output_manifest_path, ["file", "bytes", "sha256"], outputs)
    provenance = {
        "schema_version": 1,
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": "python frailty_paper/scripts/55_build_pubmed_window_review_records.py",
        "script_sha256": sha256(Path(__file__).resolve()),
        "window_manifest": str(increment_dir / "window_manifest.json"),
        "window_manifest_sha256": sha256(increment_dir / "window_manifest.json"),
        "reconciliation_report": str(reconciliation_path.resolve()),
        "reconciliation_sha256": sha256(reconciliation_path.resolve()),
        "window_start_date": window["start_date"], "window_end_date": window["end_date"],
        "frozen_protocol_cutoff": "2026-09-22",
        "xml_occurrences": len(occurrences), "window_unique_records": len(unique_rows),
        "duplicate_occurrences_within_window": len(occurrences) - len(unique_rows),
        "prior_snapshot_matches": snapshot_match_count,
        "snapshot_record_build_manifest_sha256": sha256(snapshot_manifest_path),
        "snapshot_all_records_sha256": sha256(snapshot_records_path),
        "snapshot_unique_records": snapshot_build["unique_records"],
        "snapshot_xml_records": snapshot_build["xml_records"],
        "main_screening_queue_mutated": False,
        "prisma_counts_changed": False,
        "screening_decisions_created": 0,
        "integration_status": "SEPARATE_SUPPLEMENTAL_WINDOW_UNSCREENED",
        "source_manifest": "pubmed_window_source_manifest.tsv",
        "output_manifest": output_manifest_path.name,
        "output_manifest_sha256": sha256(output_manifest_path),
        "outputs": outputs,
        "esearch_reconciliation": {
            "checked_at_utc": reconciliation.get("checked_at_utc"),
            "query_errors": reconciliation.get("query_errors"),
            "missing_segment_id_occurrences": reconciliation.get("missing_segment_id_occurrences"),
            "extra_segment_id_occurrences": reconciliation.get("extra_segment_id_occurrences"),
            "pubmed_article_nodes": reconciliation.get("pubmed_article_nodes"),
            "pubmed_book_article_nodes": reconciliation.get("pubmed_book_article_nodes"),
        },
    }
    provenance_path = outdir / "pubmed_window_review_provenance.json"
    temporary = provenance_path.with_name(provenance_path.name + f".{os.getpid()}.partial")
    temporary.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, provenance_path)
    print(json.dumps(provenance, indent=2, sort_keys=True))
    return provenance


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--increment-dir", type=Path, required=True)
    parser.add_argument("--review-dir", type=Path, default=ROOT / "frailty_paper/review")
    parser.add_argument("--outdir", type=Path, default=ROOT / "frailty_paper/review/reconciliation_2026-09-25_window")
    parser.add_argument("--reconciliation", type=Path, default=ROOT / "frailty_paper/review/reconciliation_2026-09-25_window/pubmed_retrieval_reconciliation.json")
    parser.add_argument("--source-prefix", default="pubmed_cache_increments/2026-09-23_2026-09-25")
    args = parser.parse_args(argv)
    try:
        build(args.increment_dir, args.review_dir, args.outdir, args.reconciliation, args.source_prefix)
    except (OSError, ValueError, RuntimeError, KeyError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
