#!/usr/bin/env python3
"""Materialize one locked ABC enhancer-gene task from the global overlap cache."""
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--cache")
    parser.add_argument("--cache-provenance")
    parser.add_argument("--out", required=True)
    parser.add_argument("--provenance-out", required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    manifest_path = root / args.manifest
    _, tasks = read_tsv(manifest_path)
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    if (
        task["analysis_family"] != "regulatory" or task["method"] != "activity_by_contact"
        or task["source_id"] != policy["abc_2021"]["source_id"]
    ):
        fail("task is not supported by the ABC task adapter")
    cache_path = root / (args.cache or policy["abc_2021"]["cache_path"])
    cache_provenance_path = root / (
        args.cache_provenance or policy["abc_2021"]["cache_provenance_path"]
    )
    if not cache_path.is_file() or not cache_provenance_path.is_file():
        fail("ABC overlap cache and provenance must be complete before a task runs")
    cache_provenance = json.loads(cache_provenance_path.read_text(encoding="utf-8"))
    variants_path, genes_path = root / "results/atlas/variants.tsv", root / "results/atlas/genes.tsv"
    if (
        cache_provenance.get("policy_sha256") != sha256(policy_path)
        or cache_provenance.get("source_id") != task["source_id"]
        or cache_provenance.get("source_release") != task["source_release"]
        or cache_provenance.get("variants_sha256") != sha256(variants_path)
        or cache_provenance.get("genes_sha256") != sha256(genes_path)
        or cache_provenance.get("cache_sha256") != sha256(cache_path)
    ):
        fail("ABC overlap cache differs from its checksum-bound inputs or task source")
    fields, cached = read_tsv(cache_path)
    if fields != policy["regulatory_mapping"]["canonical_fields"]:
        fail("ABC overlap cache differs from the locked regulatory schema")
    rows = [
        row for row in cached
        if row["locus_id"] == task["locus_id"] and row["context_domain"] == task["domain"]
        and row["source_dataset"] == task["source_id"] and row["source_version"] == task["source_release"]
    ]
    rows.sort(key=lambda row: (
        row["variant_id"], row["regulatory_element_id"], row["biosample"],
        row["target_gene_id"], row["regulatory_evidence_id"],
    ))
    payload = table_text(fields, rows)
    out_path = root / args.out
    atomic_text(out_path, payload)
    status = "COMPLETED" if rows else "NO_EVIDENCE_FOUND"
    reason = (
        f"Pinned ABC cache produced {len(rows)} qualified enhancer-gene evidence rows"
        if rows else
        "No eligible fine-mapped variant overlapped a non-self-promoter ABC link to an exactly matched supported same-locus gene in the locked domain"
    )
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"], "terminal_status": status, "terminal_reason": reason,
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"], "source_id": task["source_id"],
        "source_release": task["source_release"], "cache_path": str(cache_path),
        "cache_sha256": sha256(cache_path), "cache_provenance_path": str(cache_provenance_path),
        "cache_provenance_sha256": sha256(cache_provenance_path), "result_rows": len(rows),
        "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "claim_limit": policy["abc_2021"]["claim_limit"],
    }
    provenance_path = root / args.provenance_out
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"ABC_TASK_OK task={task['task_id']} status={status} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
