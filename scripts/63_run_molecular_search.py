#!/usr/bin/env python3
"""Dispatch one locked molecular search, failing closed for curator-mediated sources."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import molecular_contract


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("search_task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--plan", default="results/tables/molecular_search_plan.tsv")
    parser.add_argument("--plan-lock", default="results/tables/molecular_search_plan.lock.json")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        fail("source search requires explicit --execute after reviewing the locked plan")
    root = Path(args.root).resolve()
    plan_path, lock_path = root / args.plan, root / args.plan_lock
    policy_path = root / "config/molecular_analysis_policy.json"
    preflight_path = root / "results/tables/molecular_preflight.json"
    plan, lock = molecular_contract.validate_search_plan(
        root, plan_path, lock_path, policy_path, preflight_path,
    )
    selected = [row for row in plan if row["search_task_id"] == args.search_task_id]
    if len(selected) != 1:
        fail("search_task_id must identify exactly one locked task")
    task = selected[0]
    if task["query_mode"] == "TABIX_GRCH38_INTERVAL":
        command = [
            sys.executable, str(root / "scripts/63_query_eqtl_catalogue.py"), args.search_task_id,
            "--root", str(root), "--plan", args.plan, "--plan-lock", args.plan_lock, "--execute",
        ]
        return subprocess.run(command, cwd=root).returncode
    fail(
        f"{args.search_task_id} uses curator-mediated mode {task['query_mode']}; obtain the exact "
        "source snapshot, normalize it if accessible, then use scripts/63_record_molecular_search.py "
        "to record DATA_READY or a documented terminal outcome"
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
