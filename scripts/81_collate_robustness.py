#!/usr/bin/env python3
"""Collate the complete locked robustness matrix with deterministic provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import downstream_contract

def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, str]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--downstream-policy", default="config/downstream_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/robustness_task_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/robustness_task_manifest.lock.json")
    parser.add_argument("--out", default="results/tables/robustness_summary.tsv")
    parser.add_argument("--provenance-out", default="results/tables/robustness.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, downstream_path = root / args.policy, root / args.downstream_policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    downstream = json.loads(downstream_path.read_text(encoding="utf-8"))
    manifest_path, lock_path = root / args.manifest, root / args.manifest_lock
    tasks, lock = downstream_contract.validate_robustness_manifest(
        root, manifest_path, lock_path, policy_path, downstream_path,
    )
    output_fields = downstream["robustness"]["table_fields"]
    rows: list[dict[str, str]] = []
    provenance_hashes: dict[str, str] = {}
    for task in tasks:
        result_path, provenance_path = root / task["normalized_result_path"], root / task["provenance_path"]
        if not result_path.is_file() or not provenance_path.is_file():
            fail(f"robustness task is incomplete: {task['task_id']}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if (
            provenance.get("task_id") != task["task_id"]
            or provenance.get("policy_sha256") != sha256(policy_path)
            or provenance.get("downstream_policy_sha256") != sha256(downstream_path)
            or provenance.get("task_manifest_sha256") != sha256(manifest_path)
            or provenance.get("task_manifest_lock_sha256") != sha256(lock_path)
            or provenance.get("result_sha256") != sha256(result_path)
            or provenance.get("script_sha256")
            != downstream_contract.script_hashes(root, "robustness")
        ):
            fail(f"robustness task provenance drifted: {task['task_id']}")
        result_fields, values = read_tsv(result_path)
        if result_fields != output_fields or len(values) != 1:
            fail(f"robustness task result schema/count drifted: {task['task_id']}")
        row = values[0]
        if (
            row["conclusion_id"] != task["conclusion_id"]
            or row["conclusion"] != task["conclusion"]
            or row["robustness_family"] != task["robustness_family"]
            or row["applicable"] != task["applicable"]
            or row["primary_result"] != task["primary_result"]
        ):
            fail(f"robustness task identity drifted: {task['task_id']}")
        normalized = dict(row)
        normalized["provenance_id"] = "ROB:" + sha256(provenance_path)
        rows.append(normalized)
        provenance_hashes[task["task_id"]] = sha256(provenance_path)
    expected_rows = lock.get("conclusion_count", 0) * lock.get("family_count", 0)
    if len(rows) != expected_rows or len(rows) != len(tasks):
        fail("robustness matrix is not exact conclusion-by-family coverage")
    payload = table_text(output_fields, rows)
    out_path = root / args.out
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "downstream_policy_sha256": sha256(downstream_path),
        "task_manifest_sha256": sha256(manifest_path), "task_manifest_lock_sha256": sha256(lock_path),
        "task_provenance_sha256": provenance_hashes, "conclusion_count": lock["conclusion_count"],
        "family_count": lock["family_count"], "row_count": len(rows),
        "output": args.out, "output_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "script_sha256": downstream_contract.script_hashes(root, "robustness"),
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        if not out_path.is_file() or out_path.read_text(encoding="utf-8") != payload:
            fail("robustness summary differs from deterministic recomputation")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            fail("robustness provenance drifted")
    else:
        if out_path.exists() or provenance_path.exists():
            fail("robustness aggregate outputs already exist; refusing overwrite")
        atomic_text(out_path, payload)
        atomic_text(provenance_path, provenance_text)
    if not args.quiet:
        print(f"ROBUSTNESS_COLLATION_OK conclusions={lock['conclusion_count']} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
