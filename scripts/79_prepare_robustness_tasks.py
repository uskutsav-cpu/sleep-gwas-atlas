#!/usr/bin/env python3
"""Freeze the exact major-conclusion by nine-family robustness matrix."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


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


def task_id(conclusion_id: str, family: str) -> str:
    digest = hashlib.sha256(f"{conclusion_id}\x1f{family}".encode()).hexdigest()[:20]
    return f"ROB__{digest}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--downstream-policy", default="config/downstream_analysis_policy.json")
    parser.add_argument("--conclusions", default="results/tables/major_conclusions.tsv")
    parser.add_argument("--out", default="results/tables/robustness_task_manifest.tsv")
    parser.add_argument("--lock-out", default="results/tables/robustness_task_manifest.lock.json")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, downstream_path = root / args.policy, root / args.downstream_policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    downstream = json.loads(downstream_path.read_text(encoding="utf-8"))
    families = policy["robustness"]["required_families"]
    if families != downstream["robustness"]["required_families"]:
        fail("robustness family differs between locked policies")
    conclusion_path = root / args.conclusions
    fields, conclusions = read_tsv(conclusion_path)
    if fields != policy["major_conclusion_fields"] or not conclusions:
        fail("major conclusion family is empty or schema-drifted")
    conclusion_ids = [row["conclusion_id"] for row in conclusions]
    if len(conclusion_ids) != len(set(conclusion_ids)):
        fail("major conclusion IDs are duplicated")
    allowed_types = set(policy["robustness"]["major_conclusion_edge_levels"])
    applicability = policy["robustness"]["applicability_by_conclusion_type"]
    if set(applicability) != allowed_types:
        fail("robustness applicability matrix differs from major-conclusion types")
    tasks: list[dict[str, str]] = []
    for conclusion in conclusions:
        conclusion_type = conclusion["conclusion_type"]
        if conclusion_type not in allowed_types:
            fail(f"unknown major-conclusion type: {conclusion_type}")
        matrix = applicability[conclusion_type]
        if list(matrix) != families:
            fail(f"applicability matrix is not in locked family order: {conclusion_type}")
        evidence = root / conclusion["evidence_path"]
        if not evidence.is_file() or sha256(evidence) != conclusion["evidence_sha256"]:
            fail(f"major-conclusion evidence differs: {conclusion['conclusion_id']}")
        for family in families:
            reason = matrix[family]
            applicable = reason == "APPLICABLE"
            identity = task_id(conclusion["conclusion_id"], family)
            base = f"results/robustness/runs/{identity}"
            tasks.append({
                "task_id": identity, "conclusion_id": conclusion["conclusion_id"],
                "conclusion": conclusion["conclusion"], "conclusion_type": conclusion_type,
                "entity_id": conclusion["entity_id"], "robustness_family": family,
                "applicable": "TRUE" if applicable else "FALSE",
                "applicability_reason": reason, "primary_result": conclusion["primary_result"],
                "primary_evidence_path": conclusion["evidence_path"],
                "primary_evidence_sha256": conclusion["evidence_sha256"],
                "normalized_result_path": f"{base}/normalized.tsv",
                "provenance_path": f"{base}/provenance.json",
            })
    task_fields = policy["robustness"]["task_manifest_fields"]
    if any(list(row) != task_fields for row in tasks):
        fail("robustness task manifest schema differs from policy")
    task_ids = [row["task_id"] for row in tasks]
    if len(task_ids) != len(set(task_ids)) or len(tasks) != len(conclusions) * len(families):
        fail("robustness task matrix is duplicated or incomplete")
    run_root = root / "results/robustness/runs"
    if run_root.exists() and any(path.is_file() for path in run_root.rglob("*")):
        fail("robustness results exist before the task matrix is locked")
    payload = table_text(task_fields, tasks)
    lock = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "downstream_policy_sha256": sha256(downstream_path),
        "major_conclusions_sha256": sha256(conclusion_path),
        "task_manifest": args.out, "task_manifest_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "task_ids_in_locked_order": task_ids, "conclusion_ids_in_locked_order": conclusion_ids,
        "conclusion_count": len(conclusions), "task_count": len(tasks), "family_count": len(families),
        "applicable_task_count": sum(row["applicable"] == "TRUE" for row in tasks),
        "robustness_results_accessed_before_task_lock": False,
    }
    lock_text = json.dumps(lock, indent=2, sort_keys=True) + "\n"
    out_path, lock_path = root / args.out, root / args.lock_out
    if args.validate_only:
        if not out_path.is_file() or out_path.read_text(encoding="utf-8") != payload:
            fail("robustness task manifest differs from deterministic recomputation")
        if not lock_path.is_file() or lock_path.read_text(encoding="utf-8") != lock_text:
            fail("robustness task lock differs from deterministic recomputation")
    else:
        if out_path.exists() or lock_path.exists():
            fail("robustness task lock already exists; refusing overwrite")
        atomic_text(out_path, payload)
        atomic_text(lock_path, lock_text)
    print(
        f"ROBUSTNESS_TASKS_LOCKED conclusions={len(conclusions)} tasks={len(tasks)} "
        f"applicable={lock['applicable_task_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
