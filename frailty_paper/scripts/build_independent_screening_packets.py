#!/usr/bin/env python3
"""Build separate, auditable title/abstract screening packets for two reviewers.

The exporter never makes screening decisions and never edits the canonical queue.
Packets contain only their own reviewer's blank decision fields, avoiding exposure
to the other reviewer's in-progress judgments.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_QUEUE = ROOT / "frailty_paper/review/screening/title_abstract_queue.tsv"
DEFAULT_PROTOCOL = ROOT / "frailty_paper/review/protocol.md"
SOURCE_FIELDS = (
    "screening_id", "pmid", "doi", "year", "title", "abstract",
    "query_id", "duplicate_count", "duplicate_sources",
)
RESPONSE_FIELDS = ("decision", "reason", "reviewer", "notes")
DECISION_COLUMNS = (
    "title_abstract_decision_r1", "title_abstract_reason_r1", "reviewer_1",
    "title_abstract_decision_r2", "title_abstract_reason_r2", "reviewer_2",
    "adjudication", "adjudicator",
)
VALID_DECISIONS = {"include", "exclude", "unclear"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_rows_sha256(rows: Iterable[dict[str, str]], fields: tuple[str, ...]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        encoded = json.dumps([row.get(field, "") for field in fields],
                             ensure_ascii=False, separators=(",", ":"))
        digest.update(encoded.encode("utf-8") + b"\n")
    return digest.hexdigest()


def read_queue(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        missing = set(SOURCE_FIELDS + DECISION_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"queue missing required columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    ids = [row["screening_id"].strip() for row in rows]
    if any(not item for item in ids):
        raise ValueError("queue contains an empty screening_id")
    if len(ids) != len(set(ids)):
        raise ValueError("queue contains duplicate screening_id values")
    decided = [row["screening_id"] for row in rows
               if any(row.get(column, "").strip() for column in DECISION_COLUMNS)]
    if decided:
        raise ValueError(f"refusing to export over {len(decided)} records with existing review data")
    return rows, list(reader.fieldnames or ())


def reviewer_row(source: dict[str, str], slot: str) -> dict[str, str]:
    return {**{field: source.get(field, "") for field in SOURCE_FIELDS},
            "packet_slot": slot, **{field: "" for field in RESPONSE_FIELDS}}


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())


def build_packets(queue: Path, protocol: Path, output_dir: Path,
                  batch_size: int = 500) -> dict[str, object]:
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    rows, _ = read_queue(queue)
    response_fields = ("screening_id", *SOURCE_FIELDS[1:], "packet_slot", *RESPONSE_FIELDS)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        manifest: dict[str, object] = {
            "schema_version": 1,
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "queue_path": str(queue.resolve()),
            "queue_sha256": sha256_file(queue),
            "protocol_path": str(protocol.resolve()),
            "protocol_sha256": sha256_file(protocol),
            "record_count": len(rows),
            "batch_size": batch_size,
            "source_fields": list(SOURCE_FIELDS),
            "response_fields": list(RESPONSE_FIELDS),
            "reviewers": {},
        }
        for slot in ("reviewer_1", "reviewer_2"):
            packets = []
            for start in range(0, len(rows), batch_size):
                batch = rows[start:start + batch_size]
                packet_rows = [reviewer_row(row, slot) for row in batch]
                relative = Path(slot) / f"batch_{start // batch_size + 1:04d}.tsv"
                packet_path = temporary / relative
                packet_path.parent.mkdir(parents=True, exist_ok=True)
                write_tsv(packet_path, response_fields, packet_rows)
                packets.append({
                    "path": relative.as_posix(),
                    "row_count": len(packet_rows),
                    "screening_ids_sha256": stable_rows_sha256(packet_rows, ("screening_id",)),
                    "source_rows_sha256": stable_rows_sha256(packet_rows, SOURCE_FIELDS),
                    "template_sha256": sha256_file(packet_path),
                })
            manifest["reviewers"][slot] = {"packet_count": len(packets), "packets": packets}
        (temporary / "packet_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        os.replace(temporary, output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=500)
    args = parser.parse_args()
    try:
        manifest = build_packets(args.queue, args.protocol, args.output_dir, args.batch_size)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir.resolve()),
                      "records_per_reviewer": manifest["record_count"],
                      "reviewer_1_packets": manifest["reviewers"]["reviewer_1"]["packet_count"],
                      "reviewer_2_packets": manifest["reviewers"]["reviewer_2"]["packet_count"],
                      "source_queue_sha256": manifest["queue_sha256"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
