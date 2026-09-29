#!/usr/bin/env python3
"""Validate the frozen frailty plan and every file identity in its lock."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(root: Path) -> list[str]:
    root = root.resolve()
    plan_rel = "frailty_paper/config/analysis_plan_v1.yaml"
    lock_rel = "frailty_paper/config/analysis_plan_v1.lock.json"
    plan_path, lock_path = root / plan_rel, root / lock_rel
    errors: list[str] = []
    try:
        plan: Any = yaml.safe_load(plan_path.read_text(encoding="utf-8"))
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, yaml.YAMLError) as exc:
        return [f"cannot load frozen plan or lock: {exc}"]
    if not isinstance(plan, dict):
        return ["plan YAML root must be a mapping"]

    if lock.get("schema_version") != "frailty-analysis-plan-lock.1":
        errors.append("unexpected plan lock schema")
    if lock.get("plan_id") != "frailty_sleep_genetics_v1" or plan.get("plan_id") != lock.get("plan_id"):
        errors.append("plan ID does not match the frozen lock")
    if lock.get("plan_path") != plan_rel:
        errors.append("locked plan path is unexpected")
    elif sha256(plan_path) != lock.get("plan_sha256"):
        errors.append("frozen plan SHA-256 mismatch")
    if plan.get("version") != 1 or plan.get("status") != "FROZEN_PRE_NEW_FRAILTY_ANALYSES":
        errors.append("plan version/status differs from the frozen v1 contract")

    basis = lock.get("basis_file_sha256")
    if not isinstance(basis, dict) or not basis:
        errors.append("plan lock has no basis-file checksum map")
    else:
        for rel, expected in sorted(basis.items()):
            path = root / rel
            if not path.is_file():
                errors.append(f"locked basis file is missing: {rel}")
            elif sha256(path) != expected:
                errors.append(f"locked basis file SHA-256 mismatch: {rel}")

    endpoint = plan.get("estimand_and_hierarchy", {}).get("primary_endpoint", {})
    if endpoint.get("source_id") != "GCST90020053" or endpoint.get("trait_id") != "frailty":
        errors.append("primary endpoint differs from the registered FI")
    if plan.get("locked_inputs", {}).get("sleep_trait_count") != 12:
        errors.append("locked sleep-trait count must remain 12")
    if plan.get("multiplicity", {}).get("original_atlas_scope") != 396:
        errors.append("original all-396 multiplicity scope must remain locked")

    panel_path = root / "config/analysis_panel.tsv"
    try:
        import csv

        with panel_path.open(encoding="utf-8", newline="") as handle:
            panel = list(csv.DictReader(handle, delimiter="\t"))
        if len(panel) != 45 or sum(row.get("domain") == "sleep" for row in panel) != 12:
            errors.append("locked atlas panel must contain 45 traits and 12 sleep traits")
    except OSError as exc:
        errors.append(f"cannot read locked atlas panel: {exc}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=".")
    args = parser.parse_args()
    root = Path(args.repo).resolve()
    errors = validate(root)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("ANALYSIS_PLAN_LOCK_OK version=1 sleep_traits=12 multiplicity=396")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
