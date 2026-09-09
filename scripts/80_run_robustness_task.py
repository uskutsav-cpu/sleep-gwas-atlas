#!/usr/bin/env python3
"""Dispatch one robustness task; prespecified non-applicability is automatic."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import downstream_contract

def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/robustness_task_manifest.tsv")
    parser.add_argument("--import-dir", default="work/robustness_imports")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    tasks, _ = downstream_contract.validate_robustness_manifest(
        root, root / args.manifest,
        root / "results/tables/robustness_task_manifest.lock.json",
        root / "config/interpretation_analysis_policy.json",
        root / "config/downstream_analysis_policy.json",
    )
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate robustness task: {args.task_id}")
    task = selected[0]
    command = [
        sys.executable, str(root / "scripts/80_record_robustness_task.py"), args.task_id,
        "--root", str(root),
    ]
    if task["applicable"] == "TRUE":
        import_path = root / args.import_dir / f"{args.task_id}.json"
        if not import_path.is_file():
            fail(
                f"applicable robustness task requires an explicit curator record: {import_path.relative_to(root)}"
            )
        record = json.loads(import_path.read_text(encoding="utf-8"))
        if set(record) != {"input"} or not record["input"]:
            fail("robustness curator record must contain only a real sensitivity input path")
        command.extend(["--input", str(record["input"])])
    if not args.execute:
        print("ROBUSTNESS_TASK_READY " + " ".join(command))
        return 0
    result = subprocess.run(command, cwd=root, check=False)
    if result.returncode:
        fail(f"robustness task recorder failed: {args.task_id}")
    print(
        f"ROBUSTNESS_TASK_DISPATCHED task={args.task_id} "
        f"family={task['robustness_family']} applicable={task['applicable']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
