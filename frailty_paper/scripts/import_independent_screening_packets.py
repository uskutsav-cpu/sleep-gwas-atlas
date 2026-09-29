#!/usr/bin/env python3
"""Validate and merge one reviewer's packet responses into the canonical queue."""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import tempfile
from pathlib import Path

from build_independent_screening_packets import (
    DEFAULT_PROTOCOL, DEFAULT_QUEUE, RESPONSE_FIELDS, SOURCE_FIELDS,
    VALID_DECISIONS, sha256_file, stable_rows_sha256,
)


def load_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return list(reader), list(reader.fieldnames or ())


def _row_note(existing: str, slot: str, note: str) -> str:
    if not note:
        return existing
    prefix = f"{slot}_note: "
    lines = [line for line in existing.splitlines() if not line.startswith(prefix)]
    lines.append(prefix + note)
    return "\n".join(lines)


def import_packets(packet_dir: Path, queue_path: Path, protocol_path: Path,
                   slot: str, apply: bool = False) -> dict[str, object]:
    if slot not in {"reviewer_1", "reviewer_2"}:
        raise ValueError("slot must be reviewer_1 or reviewer_2")
    manifest_path = packet_dir / "packet_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported packet manifest schema")
    if sha256_file(protocol_path) != manifest.get("protocol_sha256"):
        raise ValueError("protocol changed since packet creation; rebuild packets under the amended protocol")
    packet_meta = manifest.get("reviewers", {}).get(slot)
    if not isinstance(packet_meta, dict):
        raise ValueError(f"manifest has no packets for {slot}")
    queue_rows, queue_fields = load_tsv(queue_path)
    required_queue = set(SOURCE_FIELDS) | {
        "title_abstract_decision_r1", "title_abstract_reason_r1", "reviewer_1",
        "title_abstract_decision_r2", "title_abstract_reason_r2", "reviewer_2",
        "adjudication", "adjudicator", "notes",
    }
    missing = required_queue - set(queue_fields)
    if missing:
        raise ValueError(f"canonical queue missing columns: {', '.join(sorted(missing))}")
    queue_by_id = {row["screening_id"]: row for row in queue_rows}
    if len(queue_by_id) != len(queue_rows):
        raise ValueError("canonical queue contains duplicate screening IDs")
    packet_fields = [*SOURCE_FIELDS, "packet_slot", *RESPONSE_FIELDS]
    expected_packets = packet_meta.get("packets", [])
    if len(expected_packets) != packet_meta.get("packet_count"):
        raise ValueError(f"manifest packet count mismatch for {slot}")
    seen: set[str] = set()
    imported = already_present = changed_rows = 0
    decision_col = "title_abstract_decision_r1" if slot == "reviewer_1" else "title_abstract_decision_r2"
    reason_col = "title_abstract_reason_r1" if slot == "reviewer_1" else "title_abstract_reason_r2"
    reviewer_col = "reviewer_1" if slot == "reviewer_1" else "reviewer_2"
    other_reviewer_col = "reviewer_2" if slot == "reviewer_1" else "reviewer_1"
    existing_slot_reviewers = {row[reviewer_col].strip() for row in queue_rows if row[reviewer_col].strip()}
    existing_other_reviewers = {row[other_reviewer_col].strip() for row in queue_rows if row[other_reviewer_col].strip()}
    if len(existing_slot_reviewers) > 1:
        raise ValueError(f"canonical queue has inconsistent reviewer identities in {slot}")
    slot_reviewer: str | None = next(iter(existing_slot_reviewers), None)

    for packet in expected_packets:
        relative = Path(packet["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("unsafe packet path in manifest")
        path = packet_dir / relative
        rows, fields = load_tsv(path)
        if fields != packet_fields:
            raise ValueError(f"unexpected packet columns: {path}")
        if len(rows) != packet["row_count"]:
            raise ValueError(f"packet row count changed: {path}")
        if stable_rows_sha256(rows, ("screening_id",)) != packet["screening_ids_sha256"]:
            raise ValueError(f"packet screening IDs changed: {path}")
        if stable_rows_sha256(rows, tuple(SOURCE_FIELDS)) != packet["source_rows_sha256"]:
            raise ValueError(f"packet source content changed: {path}")
        for row in rows:
            sid = row["screening_id"]
            if sid in seen:
                raise ValueError(f"duplicate screening ID across packets: {sid}")
            seen.add(sid)
            current = queue_by_id.get(sid)
            if current is None:
                raise ValueError(f"packet record is absent from canonical queue: {sid}")
            if any(current.get(field, "") != row.get(field, "") for field in SOURCE_FIELDS):
                raise ValueError(f"packet source metadata no longer matches canonical queue: {sid}")
            if row["packet_slot"] != slot:
                raise ValueError(f"packet assigned to wrong reviewer slot: {path}")
            decision = row["decision"].strip().lower()
            reason, reviewer, note = row["reason"].strip(), row["reviewer"].strip(), row["notes"].strip()
            if not decision:
                if reason or reviewer or note:
                    raise ValueError(f"response metadata without a decision for {sid}")
                continue
            if decision not in VALID_DECISIONS:
                raise ValueError(f"invalid decision {decision!r} for {sid}")
            if not reviewer:
                raise ValueError(f"reviewer identity required for {sid}")
            if slot_reviewer is None:
                slot_reviewer = reviewer
            elif reviewer != slot_reviewer:
                raise ValueError(f"reviewer identity changed within {slot}: {slot_reviewer!r} vs {reviewer!r}")
            if reviewer in existing_other_reviewers:
                raise ValueError("the two independent decisions name the same reviewer")
            if current[other_reviewer_col].strip() and current[other_reviewer_col].strip() == reviewer:
                raise ValueError(f"the two independent decisions name the same reviewer for {sid}")
            previous = (current[decision_col].strip().lower(), current[reason_col].strip(), current[reviewer_col].strip())
            incoming = (decision, reason, reviewer)
            if any(previous) and previous != incoming:
                raise ValueError(f"refusing to overwrite conflicting {slot} decision for {sid}")
            if previous == incoming:
                already_present += 1
            else:
                current[decision_col], current[reason_col], current[reviewer_col] = decision, reason, reviewer
                imported += 1
                changed_rows += 1
            updated_note = _row_note(current["notes"], slot, note)
            if updated_note != current["notes"]:
                current["notes"] = updated_note
                changed_rows += 1

    if seen != set(queue_by_id):
        missing_ids = len(set(queue_by_id) - seen)
        extra_ids = len(seen - set(queue_by_id))
        raise ValueError(f"packet coverage mismatch (missing={missing_ids}, extra={extra_ids})")

    if apply and changed_rows:
        old_hash = sha256_file(queue_path)
        backup = queue_path.with_name(f"{queue_path.stem}.before_packets_{old_hash[:12]}{queue_path.suffix}")
        if not backup.exists():
            shutil.copy2(queue_path, backup)
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
    return {"reviewer_slot": slot, "queue_records": len(queue_rows),
            "packet_files": len(expected_packets), "imported_decisions": imported,
            "already_present": already_present, "changed_rows": changed_rows,
            "applied": bool(apply and changed_rows),
            "backup_path": str(backup) if apply and changed_rows else None}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-dir", type=Path, required=True)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--slot", choices=("reviewer_1", "reviewer_2"), required=True)
    parser.add_argument("--apply", action="store_true", help="atomically merge validated decisions into the queue")
    args = parser.parse_args()
    try:
        result = import_packets(args.packet_dir, args.queue, args.protocol, args.slot, args.apply)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
