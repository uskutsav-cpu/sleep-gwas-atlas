#!/usr/bin/env python3
"""Validate and record one locked interpretation task with terminal provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path


SHA256 = re.compile(r"^[0-9a-f]{64}$")
MISSING = {"", "NA"}
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


def probability(value: str, field: str, identity: str, *, missing_ok: bool = True) -> None:
    if value in MISSING and missing_ok:
        return
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"nonnumeric {field} for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        fail(f"invalid {field} for {identity}")


def finite(value: str, field: str, identity: str, *, missing_ok: bool = True) -> None:
    if value in MISSING and missing_ok:
        return
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"nonnumeric {field} for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed):
        fail(f"nonfinite {field} for {identity}")


def validate_rows(
    family: str, fields: list[str], rows: list[dict[str, str]], task: dict[str, str],
    policy: dict[str, object], root: Path,
) -> None:
    section_name = {
        "regulatory": "regulatory_mapping", "cell_type": "cell_types",
        "pathway": "pathways", "causal": "causal_inference",
    }[family]
    expected_fields = policy[section_name]["canonical_fields"]
    if fields != expected_fields:
        fail(f"{family} normalized input header differs from the locked schema")
    keys = {
        "regulatory": ["regulatory_evidence_id"],
        "cell_type": ["cell_type_id"], "pathway": ["pathway_id", "trait_or_locus_id"],
        "causal": ["causal_test_id"],
    }[family]
    observed: set[tuple[str, ...]] = set()
    traits = {row["trait_id"] for row in read_tsv(root / "results/atlas/traits.tsv")[1]}
    loci = {row["locus_id"] for row in read_tsv(root / "results/atlas/loci.tsv")[1]}
    variants = {
        (row["locus_id"], row["variant_id"])
        for row in read_tsv(root / "results/atlas/variants.tsv")[1]
    }
    genes = {
        (row["locus_id"], row["gene_id"])
        for row in read_tsv(root / "results/atlas/genes.tsv")[1]
    }
    for line, row in enumerate(rows, start=2):
        identity = "/".join(row[field] for field in keys)
        key = tuple(row[field] for field in keys)
        if any(value in MISSING for value in key) or key in observed:
            fail(f"missing or duplicate {family} result key at line {line}")
        observed.add(key)
        if any(marker in "\t".join(row.values()).upper() for marker in FORBIDDEN):
            fail(f"forbidden non-real marker in {family} result at line {line}")
        if row.get("provenance_id", "") in MISSING:
            fail(f"missing source provenance for {identity}")
        if family == "regulatory":
            if row["regulatory_element_id"] in MISSING or row["variant_id"] in MISSING:
                fail(f"regulatory overlap lacks an element or variant for {identity}")
            if (
                row["locus_id"] != task["locus_id"]
                or row["element_type"] != task_source_layer(task, root, policy)
            ):
                fail(f"regulatory task identity differs for {identity}")
            if row["context_domain"] != task["domain"]:
                fail(f"regulatory result is outside the task domain for {identity}")
            if (row["locus_id"], row["variant_id"]) not in variants:
                fail(f"regulatory result references an absent variant for {identity}")
            target = row["target_gene_id"]
            if target not in MISSING and (row["locus_id"], target) not in genes:
                fail(f"regulatory result references an absent gene for {identity}")
            linked_layers = {"promoter", "enhancer_promoter_link", "3D_contact_where_available"}
            if row["element_type"] in linked_layers and target in MISSING:
                fail(f"linked regulatory layer lacks a supported target gene for {identity}")
            finite(row["effect"], "effect", identity)
            probability(row["p_value"], "p_value", identity)
            if row["source_dataset"] != task["source_id"] or row["source_version"] != task["source_release"]:
                fail(f"regulatory source identity differs for {identity}")
        elif family == "cell_type":
            if (
                row["method"] != task["method"] or row["domain"] != task["domain"]
                or row["trait_or_locus_id"] != task["trait_id"]
                or row["dataset"] != task["source_id"] or row["version"] != task["source_release"]
            ):
                fail(f"cell-type task identity differs for {identity}")
            finite(row["effect"], "effect", identity)
            finite(row["se"], "se", identity)
            probability(row["p_value"], "p_value", identity, missing_ok=False)
            probability(row["fdr"], "fdr", identity)
        elif family == "pathway":
            if (
                row["resource"] != task_source_layer(task, root, policy)
                or row["dataset_version"] != task["source_release"]
            ):
                fail(f"pathway resource differs for {identity}")
            if row["trait_or_locus_id"] not in traits | loci:
                fail(f"pathway result references an absent trait/locus for {identity}")
            try:
                size = int(row["gene_set_size"])
            except ValueError as exc:
                fail(f"invalid gene_set_size for {identity}")
                raise AssertionError from exc
            if not policy["pathways"]["minimum_gene_set_size"] <= size <= policy["pathways"]["maximum_gene_set_size"]:
                fail(f"out-of-policy gene set size for {identity}")
            finite(row["effect"], "effect", identity)
            probability(row["p_value"], "p_value", identity, missing_ok=False)
            probability(row["fdr"], "fdr", identity)
        else:
            allowed_method = (
                row["method"] == task["method"]
                or task["method"] == "CAUSE_or_LHC_MR_where_appropriate" and row["method"] in {"CAUSE", "LHC_MR"}
            )
            if (
                not allowed_method or row["exposure"] != task["exposure"]
                or row["outcome"] != task["outcome"] or row["direction"] != task["direction"]
            ):
                fail(f"causal task identity differs for {identity}")
            try:
                instruments = int(row["n_instruments"])
            except ValueError as exc:
                fail(f"invalid n_instruments for {identity}")
                raise AssertionError from exc
            if instruments < policy["causal_inference"]["minimum_instruments"]:
                fail(f"too few instruments for canonical causal result {identity}")
            finite(row["effect"], "effect", identity, missing_ok=False)
            finite(row["se"], "se", identity, missing_ok=False)
            finite(row["f_statistic_min"], "f_statistic_min", identity, missing_ok=False)
            if float(row["f_statistic_min"]) < policy["causal_inference"]["minimum_F_statistic"]:
                fail(f"weak instrument in canonical causal result {identity}")
            for field in ("p_value", "heterogeneity_p", "pleiotropy_p"):
                probability(row[field], field, identity, missing_ok=field != "p_value")


def task_source_layer(task: dict[str, str], root: Path, policy: dict[str, object]) -> str:
    registry = read_tsv(root / policy["source_registry"])[1]
    source = next((row for row in registry if row["source_id"] == task["source_id"]), None)
    if source is None:
        fail(f"task source is absent from registry: {task['source_id']}")
    return source["layer_or_resource"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/interpretation_task_manifest.lock.json")
    parser.add_argument("--status", required=True)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--input")
    parser.add_argument("--adapter-provenance")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    manifest_path, lock_path = root / args.manifest, root / args.manifest_lock
    manifest_fields, tasks = read_tsv(manifest_path)
    if manifest_fields != policy["task_manifest_fields"]:
        fail("interpretation task manifest schema drifted")
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if (
        lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("task_manifest_sha256") != sha256(manifest_path)
        or lock.get("task_ids_in_locked_order") != [row["task_id"] for row in tasks]
    ):
        fail("interpretation task manifest differs from its pre-result lock")
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    status = args.status.upper()
    if status not in policy["allowed_terminal_statuses"]:
        fail(f"invalid terminal status: {status}")
    reason = args.reason.strip()
    if not reason or any(marker in reason.upper() for marker in FORBIDDEN):
        fail("a real, non-placeholder terminal reason is required")
    family = task["analysis_family"]
    section = {
        "regulatory": "regulatory_mapping", "cell_type": "cell_types",
        "pathway": "pathways", "causal": "causal_inference",
    }[family]
    fields = policy[section]["canonical_fields"]
    rows: list[dict[str, str]] = []
    external: dict[str, object] | None = None
    adapter_provenance: dict[str, object] | None = None
    adapter_result_path: Path | None = None
    if args.adapter_provenance:
        adapter_path = Path(args.adapter_provenance).resolve()
        if not adapter_path.is_file() or adapter_path.stat().st_size == 0:
            fail("automatic adapter provenance is absent")
        adapter = json.loads(adapter_path.read_text(encoding="utf-8"))
        if adapter.get("task_id") != args.task_id or adapter.get("terminal_status") != status:
            fail("automatic adapter provenance differs from the recorded task/status")
        adapter_result_path = Path(str(adapter.get("normalized_result_path", ""))).resolve()
        if (
            not adapter_result_path.is_file()
            or sha256(adapter_result_path) != adapter.get("normalized_result_sha256")
        ):
            fail("automatic adapter normalized result differs from its provenance")
        adapter_provenance = {
            "path": str(adapter_path), "bytes": adapter_path.stat().st_size,
            "sha256": sha256(adapter_path),
        }
    if status == "COMPLETED":
        if not args.input:
            fail("COMPLETED requires a normalized real result input")
        input_path = Path(args.input).resolve()
        if adapter_result_path is not None and input_path != adapter_result_path:
            fail("automatic adapter input is not the provenance-bound normalized result")
        input_fields, rows = read_tsv(input_path)
        if not rows:
            fail("COMPLETED result contains no rows; use NO_EVIDENCE_FOUND")
        validate_rows(family, input_fields, rows, task, policy, root)
        external = {"path": str(input_path), "bytes": input_path.stat().st_size, "sha256": sha256(input_path)}
    elif args.input:
        fail(f"{status} must not include a result file")
    elif adapter_result_path is not None:
        adapter_fields, adapter_rows = read_tsv(adapter_result_path)
        if adapter_fields != fields or adapter_rows:
            fail("non-completed automatic adapter result must be a locked-schema header-only table")
    payload = table_text(fields, rows)
    result_path, provenance_path = root / task["normalized_result_path"], root / task["provenance_path"]
    if result_path.exists() or provenance_path.exists():
        fail("interpretation task output already exists; refusing overwrite")
    atomic_text(result_path, payload)
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": args.task_id, "analysis_family": family, "terminal_status": status,
        "terminal_reason": reason, "policy_sha256": sha256(policy_path),
        "task_manifest_sha256": sha256(manifest_path), "task_input_scope_sha256": task["input_scope_sha256"],
        "source_id": task["source_id"], "source_release": task["source_release"],
        "external_input": external, "adapter_provenance": adapter_provenance,
        "result_path": task["normalized_result_path"],
        "result_rows": len(rows), "result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"INTERPRETATION_TASK_RECORDED task={args.task_id} status={status} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
