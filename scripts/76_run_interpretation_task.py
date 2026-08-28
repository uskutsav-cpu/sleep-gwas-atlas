#!/usr/bin/env python3
"""Dispatch one locked interpretation task from an explicit curator import record."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--import-dir", default="work/interpretation_imports")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest_path = root / args.manifest
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        tasks = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    import_path = root / args.import_dir / f"{args.task_id}.json"
    if not import_path.is_file():
        fail(
            f"interpretation task requires an explicit curator record: {import_path.relative_to(root)}; "
            "record terminal status/reason and a real normalized input path for COMPLETED"
        )
    record = json.loads(import_path.read_text(encoding="utf-8"))
    if set(record) - {"status", "reason", "input"} or not record.get("status") or not record.get("reason"):
        fail("interpretation curator record has unexpected or missing fields")
    command = [
        sys.executable, str(root / "scripts/76_record_interpretation_task.py"), args.task_id,
        "--root", str(root), "--status", str(record["status"]), "--reason", str(record["reason"]),
    ]
    if record.get("input"):
        command.extend(["--input", str(record["input"])])
    if not args.execute:
        print("INTERPRETATION_TASK_READY " + " ".join(command))
        return 0
    result = subprocess.run(command, cwd=root, check=False)
    if result.returncode:
        fail(f"interpretation task recorder failed: {args.task_id}")
    print(
        f"INTERPRETATION_TASK_DISPATCHED task={args.task_id} "
        f"family={task['analysis_family']} source={task['source_id']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
