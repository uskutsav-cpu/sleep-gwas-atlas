#!/usr/bin/env python3
"""Record one locked robustness sensitivity result or prespecified non-applicability."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


SENSITIVITY_FIELDS = [
    "sensitivity_result", "direction_concordant", "significance_concordant",
    "material_contradiction", "resolution", "evidence_path", "evidence_sha256",
]
FORBIDDEN = ("SYNTHETIC", "PLACEHOLDER", "FAKE_RESULT", "SMOKE_TEST")


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
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--downstream-policy", default="config/downstream_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/robustness_task_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/robustness_task_manifest.lock.json")
    parser.add_argument("--input")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, downstream_path = root / args.policy, root / args.downstream_policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    downstream = json.loads(downstream_path.read_text(encoding="utf-8"))
    manifest_path, lock_path = root / args.manifest, root / args.manifest_lock
    fields, tasks = read_tsv(manifest_path)
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if (
        fields != policy["robustness"]["task_manifest_fields"]
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("downstream_policy_sha256") != sha256(downstream_path)
        or lock.get("task_manifest_sha256") != sha256(manifest_path)
        or lock.get("task_ids_in_locked_order") != [row["task_id"] for row in tasks]
    ):
        fail("robustness task manifest differs from its pre-result lock")
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate robustness task: {args.task_id}")
    task = selected[0]
    evidence = root / task["primary_evidence_path"]
    if not evidence.is_file() or sha256(evidence) != task["primary_evidence_sha256"]:
        fail("primary conclusion evidence differs from the robustness task lock")
    applicable = task["applicable"] == "TRUE"
    source_input: dict[str, object] | None = None
    if applicable:
        if not args.input:
            fail("applicable robustness task requires a real sensitivity result")
        input_path = Path(args.input).resolve()
        input_fields, values = read_tsv(input_path)
        if input_fields != SENSITIVITY_FIELDS or len(values) != 1:
            fail("robustness sensitivity input must contain one exact-schema row")
        sensitivity = values[0]
        for field in ("direction_concordant", "significance_concordant", "material_contradiction"):
            if sensitivity[field] not in {"TRUE", "FALSE"}:
                fail(f"invalid {field} for {args.task_id}")
        if sensitivity["sensitivity_result"] in {"", "NA"}:
            fail("applicable robustness task lacks a sensitivity result")
        if any(marker in "\t".join(sensitivity.values()).upper() for marker in FORBIDDEN):
            fail("robustness sensitivity input contains a forbidden non-real marker")
        relative = Path(sensitivity["evidence_path"])
        sensitivity_evidence = root / relative
        if (
            relative.is_absolute() or ".." in relative.parts or not relative.parts
            or relative.parts[0] != "results" or not sensitivity_evidence.is_file()
            or sha256(sensitivity_evidence) != sensitivity["evidence_sha256"]
        ):
            fail("robustness sensitivity evidence path/checksum is invalid")
        if sensitivity["material_contradiction"] == "TRUE" and sensitivity["resolution"] in {"", "NA", "UNRESOLVED"}:
            fail("material contradiction lacks a recorded resolution")
        source_input = {"path": str(input_path), "bytes": input_path.stat().st_size, "sha256": sha256(input_path)}
    else:
        if args.input:
            fail("prespecified non-applicable robustness task must not have a result input")
        sensitivity = {
            "sensitivity_result": "NA", "direction_concordant": "NA",
            "significance_concordant": "NA", "material_contradiction": "FALSE",
            "resolution": task["applicability_reason"], "evidence_path": "NA", "evidence_sha256": "NA",
        }
    row = {
        "conclusion_id": task["conclusion_id"], "conclusion": task["conclusion"],
        "robustness_family": task["robustness_family"], "applicable": task["applicable"],
        "primary_result": task["primary_result"], **sensitivity,
        "provenance_id": "ROBTASK:" + sha256(manifest_path) + ":" + args.task_id,
    }
    output_fields = downstream["robustness"]["table_fields"]
    if list(row) != output_fields:
        fail("robustness output schema differs from downstream policy")
    payload = table_text(output_fields, [row])
    result_path, provenance_path = root / task["normalized_result_path"], root / task["provenance_path"]
    if result_path.exists() or provenance_path.exists():
        fail("robustness task output already exists; refusing overwrite")
    atomic_text(result_path, payload)
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": args.task_id, "conclusion_id": task["conclusion_id"],
        "robustness_family": task["robustness_family"], "applicable": applicable,
        "policy_sha256": sha256(policy_path), "downstream_policy_sha256": sha256(downstream_path),
        "task_manifest_sha256": sha256(manifest_path), "source_input": source_input,
        "result_path": task["normalized_result_path"], "result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"ROBUSTNESS_TASK_RECORDED task={args.task_id} applicable={task['applicable']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
