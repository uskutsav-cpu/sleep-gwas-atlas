#!/usr/bin/env python3
"""Validate all locked interpretation tasks and build four canonical evidence tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
from collections import defaultdict
from pathlib import Path

import downstream_contract

MISSING = {"", "NA"}
ROOT = Path(__file__).resolve().parents[1]
RECORDER_SPEC = importlib.util.spec_from_file_location(
    "interpretation_recorder", ROOT / "scripts/76_record_interpretation_task.py"
)
if RECORDER_SPEC is None or RECORDER_SPEC.loader is None:
    raise RuntimeError("could not load interpretation task validation helpers")
recorder = importlib.util.module_from_spec(RECORDER_SPEC)
RECORDER_SPEC.loader.exec_module(recorder)


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


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
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


def bh(values: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted values with stable ties/order."""
    if not values:
        return []
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in values):
        fail("BH input contains an invalid probability")
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    adjusted = [1.0] * len(values)
    running = 1.0
    total = len(values)
    for rank_index in range(total - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[index] * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def fmt(value: float) -> str:
    return format(value, ".12g")


def canonicalize_identities(family: str, row: dict[str, str], task: dict[str, str]) -> None:
    if family == "regulatory":
        native = row["regulatory_element_id"]
        context = "\x1f".join((
            task["source_id"], native, row["variant_id"], row["locus_id"],
            row["biosample"], row["tissue"], row["cell_type"],
            row["context_domain"], row["target_gene_id"], row["link_method"],
            row["annotation"], row["effect"], row["p_value"],
        ))
        row["regulatory_element_id"] = f"REG::{task['source_id']}::{native}"
        row["regulatory_evidence_id"] = f"REGEV::{hashlib.sha256(context.encode()).hexdigest()[:24]}"
    elif family == "cell_type":
        native = row["cell_type_id"]
        context = "\x1f".join((
            task["trait_id"], task["source_id"], row["method"], row["domain"],
            row["tissue"], row["cell_type"], native,
        ))
        row["cell_type_id"] = f"CT::{task['source_id']}::{hashlib.sha256(context.encode()).hexdigest()[:20]}"
    elif family == "pathway":
        row["pathway_id"] = f"PATH::{row['resource']}::{row['pathway_id']}"
    else:
        native = row["causal_test_id"]
        context = "\x1f".join((task["task_id"], row["method"], native))
        row["causal_test_id"] = f"MR::{hashlib.sha256(context.encode()).hexdigest()[:24]}"


def corrected_rows(
    rows: list[dict[str, str]], group_fields: tuple[str, ...], threshold: float,
) -> list[dict[str, str]]:
    groups: dict[tuple[str, ...], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        groups[tuple(row[field] for field in group_fields)].append(index)
    result = [dict(row) for row in rows]
    for indices in groups.values():
        values = [float(rows[index]["p_value"]) for index in indices]
        for index, adjusted in zip(indices, bh(values)):
            result[index]["fdr"] = fmt(adjusted)
    return [row for row in result if float(row["fdr"]) <= threshold]


def classify_cells(rows: list[dict[str, str]], minimum_strategies: int) -> None:
    methods: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for row in rows:
        key = (row["trait_or_locus_id"], row["cell_type"], row["tissue"], row["domain"])
        methods[key].add(row["method"])
    for row in rows:
        key = (row["trait_or_locus_id"], row["cell_type"], row["tissue"], row["domain"])
        row["evidence_level"] = (
            "MULTI_STRATEGY_CELL_TYPE"
            if len(methods[key]) >= minimum_strategies else "SINGLE_STRATEGY_CELL_TYPE"
        )


def classify_causal(rows: list[dict[str, str]]) -> dict[str, float]:
    fdr_by_test: dict[str, float] = {}
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        groups[(row["direction"], row["method"])].append(index)
    for indices in groups.values():
        for index, adjusted in zip(indices, bh([float(rows[index]["p_value"]) for index in indices])):
            fdr_by_test[rows[index]["causal_test_id"]] = adjusted
    by_direction: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_direction[(row["exposure"], row["outcome"], row["direction"])].append(row)
    for values in by_direction.values():
        significant = [row for row in values if fdr_by_test[row["causal_test_id"]] <= 0.05]
        signs = {math.copysign(1, float(row["effect"])) for row in significant if float(row["effect"]) != 0}
        methods = {row["method"] for row in significant}
        robust_estimators = bool(methods & {"weighted_median", "CAUSE", "LHC_MR"})
        pleiotropy_clear = all(
            row["pleiotropy_p"] in MISSING or float(row["pleiotropy_p"]) >= 0.05 for row in values
        )
        diagnostics_clear = all(
            row["steiger_status"] == "CORRECT_DIRECTION"
            and row["sample_overlap_status"] in {"CHARACTERIZED", "ADJUSTED"}
            and row["colocalization_consistency"] == "COMPATIBLE"
            for row in values
        )
        robust = (
            "IVW" in methods and robust_estimators and len(methods) >= 2 and len(signs) == 1
            and pleiotropy_clear and diagnostics_clear
        )
        for row in values:
            row["evidence_level"] = (
                "ROBUST_CAUSAL_EVIDENCE" if robust else "ESTIMATOR_RESULT_NOT_ROBUST_CAUSAL"
            )
    return fdr_by_test


def require_unique(rows: list[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    keys = [tuple(row[field] for field in fields) for row in rows]
    if len(keys) != len(set(keys)):
        fail(f"canonical {label} primary keys are duplicated")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/interpretation_task_manifest.lock.json")
    parser.add_argument("--regulatory-out", default="results/atlas/regulatory_elements.tsv")
    parser.add_argument("--cell-out", default="results/atlas/cell_types.tsv")
    parser.add_argument("--pathway-out", default="results/atlas/pathways.tsv")
    parser.add_argument("--causal-out", default="results/atlas/causal_tests.tsv")
    parser.add_argument("--coverage-out", default="results/tables/interpretation_coverage.tsv")
    parser.add_argument("--provenance-out", default="results/atlas/interpretation.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    manifest_path, lock_path = root / args.manifest, root / args.manifest_lock
    tasks, lock = downstream_contract.validate_interpretation_manifest(
        root, manifest_path, lock_path, policy_path,
    )
    by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    coverage: list[dict[str, object]] = []
    task_hashes: dict[str, str] = {}
    for task in tasks:
        result_path, provenance_path = root / task["normalized_result_path"], root / task["provenance_path"]
        if not result_path.is_file() or not provenance_path.is_file():
            fail(f"interpretation task is incomplete: {task['task_id']}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        status = provenance.get("terminal_status")
        if (
            provenance.get("task_id") != task["task_id"]
            or provenance.get("policy_sha256") != sha256(policy_path)
            or provenance.get("task_manifest_sha256") != sha256(manifest_path)
            or provenance.get("task_manifest_lock_sha256") != sha256(lock_path)
            or provenance.get("task_input_scope_sha256") != task["input_scope_sha256"]
            or status not in policy["allowed_terminal_statuses"]
            or provenance.get("result_sha256") != sha256(result_path)
            or provenance.get("script_sha256")
            != downstream_contract.script_hashes(root, "interpretation")
        ):
            fail(f"interpretation task provenance drifted: {task['task_id']}")
        family = task["analysis_family"]
        section = {
            "regulatory": "regulatory_mapping", "cell_type": "cell_types",
            "pathway": "pathways", "causal": "causal_inference",
        }[family]
        result_fields, result_rows = read_tsv(result_path)
        if result_fields != policy[section]["canonical_fields"]:
            fail(f"interpretation result schema drifted: {task['task_id']}")
        if len(result_rows) != provenance.get("result_rows"):
            fail(f"interpretation result row count drifted: {task['task_id']}")
        if (status == "COMPLETED") != bool(result_rows):
            fail(f"terminal status/result presence mismatch: {task['task_id']}")
        recorder.validate_rows(family, result_fields, result_rows, task, policy, root)
        provenance_id = "INT:" + sha256(provenance_path)
        for row in result_rows:
            normalized = dict(row)
            canonicalize_identities(family, normalized, task)
            normalized["provenance_id"] = provenance_id
            by_family[family].append(normalized)
        counts = {
            "analyzed_unit_count": int(status == "COMPLETED"),
            "no_evidence_unit_count": int(status == "NO_EVIDENCE_FOUND"),
            "access_blocked_unit_count": int(status == "ACCESS_BLOCKED"),
            "not_applicable_unit_count": int(status == "NOT_APPLICABLE"),
        }
        coverage.append({
            "coverage_id": "ICOV__" + task["task_id"], "analysis_family": family,
            "unit_id": task["unit_id"], "method": task["method"], "domain": task["domain"],
            "source_id": task["source_id"], "required_unit_count": 1,
            "completed_unit_count": 1, **counts, "coverage_status": "COMPLETE",
            "evidence_path": task["normalized_result_path"], "evidence_sha256": sha256(result_path),
            "provenance_id": provenance_id,
        })
        task_hashes[task["task_id"]] = sha256(provenance_path)

    cells = corrected_rows(
        by_family["cell_type"], ("trait_or_locus_id", "method", "dataset"),
        float(policy["cell_types"]["support_fdr"]),
    )
    classify_cells(cells, int(policy["cell_types"]["minimum_independent_strategies_for_supported_claim"]))
    pathways = corrected_rows(
        by_family["pathway"], ("trait_or_locus_id", "resource"),
        float(policy["pathways"]["support_fdr"]),
    )
    for row in pathways:
        row["evidence_level"] = "FDR_SUPPORTED_PATHWAY"
    causal = [dict(row) for row in by_family["causal"]]
    causal_fdr = classify_causal(causal)
    regulatory = [dict(row) for row in by_family["regulatory"]]
    require_unique(regulatory, ("regulatory_evidence_id",), "regulatory")
    require_unique(cells, ("cell_type_id",), "cell-type")
    require_unique(pathways, ("pathway_id", "trait_or_locus_id"), "pathway")
    require_unique(causal, ("causal_test_id",), "causal")

    regulatory.sort(key=lambda row: (row["locus_id"], row["variant_id"], row["regulatory_element_id"], row["target_gene_id"]))
    cells.sort(key=lambda row: (row["trait_or_locus_id"], row["domain"], row["cell_type"], row["method"], row["cell_type_id"]))
    pathways.sort(key=lambda row: (row["trait_or_locus_id"], row["resource"], row["fdr"], row["pathway_id"]))
    causal.sort(key=lambda row: (row["exposure"], row["outcome"], row["direction"], row["method"], row["causal_test_id"]))
    coverage.sort(key=lambda row: (str(row["analysis_family"]), str(row["coverage_id"])))
    payloads = {
        root / args.regulatory_out: table_text(policy["regulatory_mapping"]["canonical_fields"], regulatory),
        root / args.cell_out: table_text(policy["cell_types"]["canonical_fields"], cells),
        root / args.pathway_out: table_text(policy["pathways"]["canonical_fields"], pathways),
        root / args.causal_out: table_text(policy["causal_inference"]["canonical_fields"], causal),
        root / args.coverage_out: table_text(policy["coverage_fields"], coverage),
    }
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(manifest_path),
        "task_manifest_lock_sha256": sha256(lock_path), "task_provenance_sha256": task_hashes,
        "task_count": len(tasks), "coverage_rows": len(coverage),
        "regulatory_rows": len(regulatory), "cell_type_rows": len(cells),
        "pathway_rows": len(pathways), "causal_rows": len(causal),
        "causal_bh_fdr_by_test_id": {key: fmt(value) for key, value in sorted(causal_fdr.items())},
        "outputs": {
            str(path.relative_to(root)): hashlib.sha256(text.encode()).hexdigest()
            for path, text in payloads.items()
        },
        "script_sha256": downstream_contract.script_hashes(root, "interpretation"),
        "claim_limit": policy["claim_limit"],
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        for path, text in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                fail(f"interpretation aggregate differs from deterministic recomputation: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            fail("interpretation aggregate provenance drifted")
    else:
        if any(path.exists() for path in [*payloads, provenance_path]):
            fail("interpretation aggregate outputs already exist; refusing overwrite")
        for path, text in payloads.items():
            atomic_text(path, text)
        atomic_text(provenance_path, provenance_text)
    if not args.quiet:
        print(
            f"INTERPRETATION_OK tasks={len(tasks)} regulatory={len(regulatory)} "
            f"cells={len(cells)} pathways={len(pathways)} causal={len(causal)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
