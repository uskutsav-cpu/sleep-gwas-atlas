#!/usr/bin/env python3
"""Normalize one locked LDSC-SEG domain task from its all-tissue trait run."""
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
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--provenance-out", required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["ldsc_seg_gtex"]
    manifest_path = root / args.manifest
    _, tasks = read_tsv(manifest_path)
    selected_task = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected_task) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected_task[0]
    if (
        task["analysis_family"] != "cell_type" or task["method"] != "stratified_LDSC"
        or task["source_id"] != spec["source_id"]
        or task["domain"] not in policy["cell_types"]["required_domains"]
    ):
        fail("task is not supported by the LDSC-SEG adapter")
    selection_path = root / spec["selection_config"]
    if sha256(selection_path) != spec["selection_config_sha256"]:
        fail("LDSC-SEG tissue selection differs from policy")
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    trait_path = root / spec["trait_result_path_template"].format(trait_id=task["trait_id"])
    trait_provenance_path = root / spec["trait_provenance_path_template"].format(
        trait_id=task["trait_id"]
    )
    fields, rows = read_tsv(trait_path)
    expected_fields = ["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"]
    if fields != expected_fields or not trait_provenance_path.is_file():
        fail("LDSC-SEG trait result or provenance is incomplete")
    trait_provenance = json.loads(trait_provenance_path.read_text(encoding="utf-8"))
    if (
        trait_provenance.get("trait_id") != task["trait_id"]
        or trait_provenance.get("policy_sha256") != sha256(policy_path)
        or trait_provenance.get("task_manifest_sha256") != sha256(manifest_path)
        or trait_provenance.get("result_sha256") != sha256(trait_path)
        or trait_provenance.get("task_input_scope_sha256", {}).get(task["domain"])
        != task["input_scope_sha256"]
    ):
        fail("LDSC-SEG trait result differs from its pre-result task lock")
    selected_tissues = [
        row for row in selection["selected_tissues"] if row["domain"] == task["domain"]
    ]
    source_rows = {row["Name"]: row for row in rows}
    normalized: list[dict[str, str]] = []
    for tissue in selected_tissues:
        source_label = tissue["source_label"]
        if source_label not in source_rows:
            fail(f"LDSC-SEG trait result omits selected tissue: {source_label}")
        row = source_rows[source_label]
        normalized.append({
            "cell_type_id": f"LDSC_SEG_GTEX::{source_label}",
            "cell_type": tissue["tissue"],
            "tissue": tissue["tissue"],
            "domain": task["domain"],
            "method": task["method"],
            "trait_or_locus_id": task["trait_id"],
            "effect": row["Coefficient"],
            "se": row["Coefficient_std_error"],
            "p_value": row["Coefficient_P_value"],
            "fdr": "NA",
            "evidence_level": "UNADJUSTED_LDSC_SEG_COEFFICIENT",
            "dataset": task["source_id"],
            "version": task["source_release"],
            "provenance_id": "LDSC_SEG_TASK_PENDING",
        })
    canonical_fields = policy["cell_types"]["canonical_fields"]
    if len(normalized) != spec["selected_tissues_by_domain"][task["domain"]]:
        fail("LDSC-SEG normalized domain family is incomplete")
    if any(list(row) != canonical_fields for row in normalized):
        fail("LDSC-SEG normalized row differs from the locked cell-type schema")
    payload = table_text(canonical_fields, normalized)
    out_path = root / args.out
    atomic_text(out_path, payload)
    provenance_path = root / args.provenance_out
    provenance = {
        "schema_version": policy["schema_version"],
        "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"],
        "terminal_status": "COMPLETED",
        "terminal_reason": (
            f"LDSC-SEG h2-cts completed for all {len(normalized)} prespecified "
            f"{task['domain']} GTEx tissue annotations"
        ),
        "policy_sha256": sha256(policy_path),
        "task_manifest_sha256": sha256(manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"],
        "source_id": task["source_id"],
        "source_release": task["source_release"],
        "trait_result_provenance_sha256": sha256(trait_provenance_path),
        "result_rows": len(normalized),
        "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "multiple_testing_family": spec["multiple_testing_family"],
        "claim_limit": spec["claim_limit"],
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"LDSC_SEG_TASK_OK task={task['task_id']} rows={len(normalized)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
