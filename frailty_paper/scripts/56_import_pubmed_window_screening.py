#!/usr/bin/env python3
"""Validate and optionally import one reviewer slot for a supplemental PubMed window."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WINDOW_ROOT = ROOT / "frailty_paper/review/incremental_windows/2026-09-26_2026-09-27"
DEFAULT_QUEUE = WINDOW_ROOT / "review/pubmed_window_unique_review_queue.tsv"
DEFAULT_PROTOCOL = ROOT / "frailty_paper/review/protocol.md"
DEFAULT_PACKET_DIR = WINDOW_ROOT / "reviewer_packets"
EXPECTED_PACKET_MANIFEST_SHA256 = "83989c4c13b4354940dd2b125a0afbfde886c89659b0c393cb4604efd859f8a4"
EXPECTED_QUEUE_SHA256 = "7461d548b509d5a035c2798e211d9b90d464018861a9193be4f23b56bbe03ba7"
EXPECTED_PROTOCOL_SHA256 = "1013526fa40de6708ed9e547e09c1cb0ae3de0c27d8cbe2d5d04a99ffa796abb"
VALID_DECISIONS = {"include", "exclude", "unclear"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_rows_sha256(rows: list[dict[str, str]], fields: tuple[str, ...]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        encoded = json.dumps([row.get(field, "") for field in fields],
                             ensure_ascii=False, separators=(",", ":"))
        digest.update(encoded.encode("utf-8") + b"\n")
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return list(reader), fields


def reviewer_columns(slot: str) -> tuple[str, str, str, str]:
    suffix = "r1" if slot == "reviewer_1" else "r2"
    return (f"title_abstract_decision_{suffix}", f"title_abstract_reason_{suffix}",
            slot, "reviewer_2" if slot == "reviewer_1" else "reviewer_1")


def _note_for_slot(existing: str, slot: str, incoming: str) -> str:
    incoming = incoming.strip()
    if not incoming:
        return existing
    marker = f"{slot}: {incoming}"
    current = [part for part in existing.split(" | ") if part]
    if marker not in current:
        current.append(marker)
    return " | ".join(current)


def import_packets(packet_dir: Path, queue_path: Path, protocol_path: Path,
                   slot: str, apply: bool = False) -> dict:
    packet_dir, queue_path, protocol_path = packet_dir.resolve(), queue_path.resolve(), protocol_path.resolve()
    manifest_path = packet_dir / "packet_manifest.json"
    if sha256_file(manifest_path) != EXPECTED_PACKET_MANIFEST_SHA256:
        raise ValueError("supplemental packet manifest hash differs from the frozen packet set")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported supplemental packet manifest schema")
    source_fields = tuple(manifest.get("source_fields", []))
    response_fields = tuple(manifest.get("response_fields", []))
    if not source_fields or response_fields != ("decision", "reason", "reviewer", "notes"):
        raise ValueError("unexpected packet source or response field definition")
    if (manifest.get("queue_sha256") != EXPECTED_QUEUE_SHA256
            or manifest.get("protocol_sha256") != EXPECTED_PROTOCOL_SHA256):
        raise ValueError("packet manifest queue/protocol hashes differ from the frozen supplemental set")
    if sha256_file(protocol_path) != EXPECTED_PROTOCOL_SHA256:
        raise ValueError("frozen review protocol hash differs from packet manifest")
    if sha256_file(queue_path) != manifest.get("queue_sha256"):
        backup = queue_path.with_name(
            f"{queue_path.stem}.before_reviewers_{manifest['queue_sha256'][:12]}{queue_path.suffix}")
        if not backup.is_file() or sha256_file(backup) != manifest.get("queue_sha256"):
            raise ValueError("supplemental queue differs from its packet baseline; verified backup missing")

    queue_rows, queue_fields = read_tsv(queue_path)
    required_queue_fields = set(source_fields) | {
        "title_abstract_decision_r1", "title_abstract_reason_r1", "reviewer_1",
        "title_abstract_decision_r2", "title_abstract_reason_r2", "reviewer_2",
        "adjudication", "adjudicator", "notes",
    }
    if not required_queue_fields.issubset(queue_fields):
        raise ValueError("supplemental queue is missing required source/review columns")
    queue_by_id = {row["screening_id"]: row for row in queue_rows}
    if len(queue_by_id) != len(queue_rows) or len(queue_rows) != manifest.get("record_count"):
        raise ValueError("supplemental queue identity/count differs from packet manifest")

    reviewer_meta = manifest.get("reviewers", {}).get(slot)
    if not reviewer_meta or len(reviewer_meta.get("packets", [])) != reviewer_meta.get("packet_count"):
        raise ValueError(f"invalid packet manifest for reviewer slot {slot}")
    decision_col, reason_col, reviewer_col, other_reviewer_col = reviewer_columns(slot)
    seen: set[str] = set()
    slot_reviewer: str | None = None
    imported = already_present = 0
    for packet_info in reviewer_meta["packets"]:
        relative = Path(packet_info["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("unsafe packet path in manifest")
        path = packet_dir / relative
        rows, fields = read_tsv(path)
        expected_fields = [*source_fields, "packet_slot", *response_fields]
        if fields != expected_fields:
            raise ValueError(f"unexpected packet columns: {path}")
        if len(rows) != packet_info.get("row_count"):
            raise ValueError(f"packet row count changed: {path}")
        if stable_rows_sha256(rows, ("screening_id",)) != packet_info.get("screening_ids_sha256"):
            raise ValueError(f"packet screening IDs changed: {path}")
        if stable_rows_sha256(rows, source_fields) != packet_info.get("source_rows_sha256"):
            raise ValueError(f"packet source rows changed: {path}")
        for row in rows:
            sid = row.get("screening_id", "")
            if sid in seen:
                raise ValueError(f"duplicate screening ID across packets: {sid}")
            seen.add(sid)
            current = queue_by_id.get(sid)
            if current is None or any(current.get(key, "") != row.get(key, "") for key in source_fields):
                raise ValueError(f"packet source metadata no longer matches supplemental queue: {sid}")
            if row.get("packet_slot") != slot:
                raise ValueError(f"packet assigned to wrong reviewer slot: {path}")
            decision = row.get("decision", "").strip().lower()
            reason = row.get("reason", "").strip()
            reviewer = row.get("reviewer", "").strip()
            note = row.get("notes", "").strip()
            if not decision:
                if reason or reviewer or note:
                    raise ValueError(f"response metadata without a decision for {sid}")
                if apply:
                    raise ValueError(f"cannot apply incomplete reviewer packet; no decision for {sid}")
                continue
            if decision not in VALID_DECISIONS:
                raise ValueError(f"invalid decision {decision!r} for {sid}")
            if not reviewer:
                raise ValueError(f"reviewer identity required for {sid}")
            if slot_reviewer is not None and slot_reviewer != reviewer:
                raise ValueError(f"reviewer identity changes within {slot}")
            slot_reviewer = reviewer
            if current[other_reviewer_col].strip() and current[other_reviewer_col].strip() == reviewer:
                raise ValueError(f"independent decisions name the same reviewer for {sid}")
            previous = (current[decision_col].strip().lower(), current[reason_col].strip(),
                        current[reviewer_col].strip())
            incoming = (decision, reason, reviewer)
            if any(previous) and previous != incoming:
                raise ValueError(f"refusing to overwrite conflicting {slot} decision for {sid}")
            if previous == incoming:
                already_present += 1
            else:
                current[decision_col], current[reason_col], current[reviewer_col] = incoming
                imported += 1
            current["notes"] = _note_for_slot(current["notes"], slot, note)

    if seen != set(queue_by_id):
        raise ValueError(f"packet coverage mismatch (missing={len(set(queue_by_id)-seen)}, extra={len(seen-set(queue_by_id))})")
    complete_pairs = sum(bool(r["title_abstract_decision_r1"] and r["title_abstract_decision_r2"])
                         for r in queue_rows)
    disagreements = [r["screening_id"] for r in queue_rows
                     if r["title_abstract_decision_r1"] and r["title_abstract_decision_r2"]
                     and r["title_abstract_decision_r1"].lower() != r["title_abstract_decision_r2"].lower()]
    backup = None
    if apply and imported:
        baseline_hash = manifest["queue_sha256"]
        backup = queue_path.with_name(
            f"{queue_path.stem}.before_reviewers_{baseline_hash[:12]}{queue_path.suffix}")
        if not backup.exists():
            if sha256_file(queue_path) != baseline_hash:
                raise ValueError("refusing to create baseline backup from an already modified queue")
            shutil.copy2(queue_path, backup)
        elif sha256_file(backup) != baseline_hash:
            raise ValueError("existing supplemental queue backup does not match packet baseline")
        fd, temp_name = tempfile.mkstemp(prefix=f".{queue_path.name}.", suffix=".tmp", dir=queue_path.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8-sig", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=queue_fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(queue_rows)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temp_name, queue_path)
        except BaseException:
            try:
                os.unlink(temp_name)
            except FileNotFoundError:
                pass
            raise
    return {
        "reviewer_slot": slot, "queue_records": len(queue_rows),
        "packet_files": len(reviewer_meta["packets"]), "imported_decisions": imported,
        "already_present": already_present, "complete_independent_pairs": complete_pairs,
        "unresolved_disagreement_ids": disagreements, "applied": bool(apply and imported),
        "backup_path": str(backup) if backup else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-dir", type=Path, default=DEFAULT_PACKET_DIR)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--slot", choices=("reviewer_1", "reviewer_2"), required=True)
    parser.add_argument("--apply", action="store_true", help="atomically import validated decisions")
    args = parser.parse_args()
    try:
        result = import_packets(args.packet_dir, args.queue, args.protocol, args.slot, args.apply)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
