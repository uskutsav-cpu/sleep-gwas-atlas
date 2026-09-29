"""Create an immutable canonical-family decision after LAVA v3 is complete.

This command is deliberately separate from the production runner. It refuses
to write any decision while the canonical 7 x 2,495 family is incomplete,
re-audits all receipts, verifies the runner's final aggregate against the
receipt-bound cell files, and records whether frozen QC permits promotion.
It never launches or changes production analyses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
sys.path.insert(0, str(ROOT / "brain6/scripts"))
import audit_lava_canonical_v3_partial as partial
import run_lava_canonical_v3 as canonical

AGGREGATE_FIELDS = [
    "phen", "locus_id", "chromosome", "status", "n_snps", "n_components",
    "h2.obs", "h2.latent", "p", "reason",
]
CELL_FIELDS = ("status", "n_snps", "n_components", "h2.obs", "h2.latent", "p", "reason")


def json_bytes(value: dict[str, Any]) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def csv_value(value: Any) -> str:
    return "" if value is None else str(value)


def verify_runner_identity(latest: dict[str, Any], identity_doc: dict[str, Any],
                           audit: dict[str, Any]) -> None:
    """Verify the runner's hashed execution identity against its receipt audit.

    ``latest_audit.json`` stores the full run-identity inputs used to derive the
    run ID, while ``run_identity.json`` stores a smaller canonical projection.
    Compare their shared bindings and verify the full identity self-hashes to
    the immutable run ID instead of comparing those different schemas as equal.
    """
    run_id = audit["run_id"]
    runtime_identity = latest.get("run_identity")
    if not isinstance(runtime_identity, dict):
        raise TypeError("Runner's final audit is missing its checksum-bound run identity")
    if (latest.get("run_id") != run_id or identity_doc.get("run_id") != run_id or
            hashlib.sha256(canonical.canonical_json(runtime_identity)).hexdigest() != run_id):
        raise ValueError("Runner's final audit is missing or bound to a different run identity")

    expected = {
        "analysis_id": identity_doc["analysis_id"],
        "family_lock_sha256": audit["family_lock_sha256"],
        "execution_lock_sha256": audit["execution_lock_sha256"],
        "input_provenance_sha256": audit["input_provenance_sha256"],
        "reference_provenance_sha256": audit["reference_provenance_sha256"],
        "batch_worker_sha256": identity_doc["script_sha256"],
        "runtime_validator_sha256": identity_doc["runtime_validator_sha256"],
        "runtime_validation_output_sha256": identity_doc["runtime_validation_output_sha256"],
        "R": identity_doc["runtime"]["R"],
        "LAVA": identity_doc["runtime"]["LAVA"],
    }
    if any(runtime_identity.get(key) != value for key, value in expected.items()):
        raise ValueError("Runner's final audit identity disagrees with verified run provenance")
    if latest.get("manifest_sha256") != identity_doc.get("canonical_manifest_sha256"):
        raise ValueError("Runner's final audit manifest differs from the frozen run identity")
    if latest.get("analysis_id") != audit["analysis_id"] or latest.get("workers") != audit["worker_count"]:
        raise ValueError("Runner's final audit analysis or worker count differs from the receipt audit")


def verify_aggregate(path: Path, run_dir: Path, loci: list[dict[str, str]],
                     traits: tuple[str, ...]) -> int:
    """Verify every aggregate row against its receipt-bound canonical cell."""
    if not path.is_file():
        raise FileNotFoundError(f"Canonical aggregate is missing: {path}")
    expected_rows: list[dict[str, str]] = []
    for locus in loci:
        cells = canonical.load_cells(canonical.result_path(run_dir, locus["LOC"]))
        by_trait = {row["phen"]: row for row in cells}
        if set(by_trait) != set(traits) or len(by_trait) != len(cells):
            raise ValueError(f"Receipt-bound cells are not unique and complete at {locus['LOC']}")
        for trait in traits:
            cell = by_trait[trait]
            expected_rows.append({
                "phen": trait,
                "locus_id": locus["LOC"],
                "chromosome": locus["CHR"],
                **{field: csv_value(cell.get(field, "")) for field in CELL_FIELDS},
            })

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != AGGREGATE_FIELDS:
            raise ValueError("Canonical aggregate schema differs from the frozen v3 schema")
        observed_rows = list(reader)
    if len(observed_rows) != len(expected_rows):
        raise ValueError(f"Canonical aggregate has {len(observed_rows)} rows; expected {len(expected_rows)}")
    seen: set[tuple[str, str]] = set()
    for index, (observed, expected) in enumerate(zip(observed_rows, expected_rows, strict=True)):
        key = (observed["phen"], observed["locus_id"])
        if key in seen:
            raise ValueError(f"Duplicate canonical aggregate identity: {key}")
        seen.add(key)
        if observed != expected:
            raise ValueError(f"Canonical aggregate row {index + 2} differs from receipt-bound cells")
    return len(observed_rows)


def build_decision(audit: dict[str, Any], family: dict[str, Any], run_identity: dict[str, Any],
                   aggregate: dict[str, Any], aggregate_rows: int) -> dict[str, Any]:
    """Derive the family decision from frozen locks and a complete audit."""
    canonical_lock = family["canonical_univariate"]
    planned_cells = int(canonical_lock["n_tests"])
    planned_loci = int(canonical_lock["loci"])
    if audit["state"] != "COMPLETE":
        raise ValueError("Refusing to finalize an incomplete canonical family")
    if audit["verified_loci"] != planned_loci or audit["planned_loci"] != planned_loci:
        raise ValueError("Completed audit does not cover the exact frozen locus family")
    if audit["verified_cells"] != planned_cells or audit["planned_cells"] != planned_cells:
        raise ValueError("Completed audit does not cover the exact frozen cell family")
    if audit["invalid_receipt_count"] != 0:
        raise ValueError("Cannot decide a family with invalid or unknown receipts")
    if aggregate_rows != planned_cells or aggregate.get("rows") != planned_cells:
        raise ValueError("Aggregate row count differs from the frozen canonical family")
    if run_identity.get("run_id") != audit["run_id"]:
        raise ValueError("Run identity differs from the receipt audit")

    gate = family["univariate_gate"]
    alpha = float(gate["familywise_alpha"])
    denominator = int(gate["n_tests"])
    threshold = float(gate["p_threshold_strictly_less_than"])
    if denominator != planned_cells or threshold != alpha / denominator:
        raise ValueError("Frozen univariate multiplicity gate is internally inconsistent")

    qc_pass = bool(audit["family_qc_pass"])
    untested = int(audit["currently_uncompleted_or_untested_cells"])
    max_untested = int(audit["maximum_allowed_untested_cells"])
    status_counts = audit["statuses"]
    return {
        "schema_version": 1,
        "decision_type": "CANONICAL_UNIVARIATE_FAMILY",
        "analysis_id": audit["analysis_id"],
        "run_id": audit["run_id"],
        "run_identity_sha256": run_identity["run_identity_sha256"],
        "family_lock_sha256": audit["family_lock_sha256"],
        "execution_lock_sha256": audit["execution_lock_sha256"],
        "canonical_manifest_sha256": run_identity["canonical_manifest_sha256"],
        "canonical_script_sha256": audit["canonical_script_sha256"],
        "aggregate": aggregate,
        "aggregate_rows_verified": aggregate_rows,
        "planned_loci": planned_loci,
        "planned_cells": planned_cells,
        "verified_cells": audit["verified_cells"],
        "status_counts": status_counts,
        "untested_cells_including_missing": untested,
        "maximum_allowed_untested_cells": max_untested,
        "untested_fraction": untested / planned_cells,
        "qc_pass": qc_pass,
        "overall_status": "PASS" if qc_pass else "FAILED_QC_NOT_PROMOTED",
        "promotion_permitted": qc_pass,
        "pairwise_stage_authorized": qc_pass,
        "pairwise_stage_executed": False,
        "multiple_testing": {
            "method": "Bonferroni fixed canonical univariate gate",
            "familywise_alpha": alpha,
            "denominator": denominator,
            "strict_p_threshold": threshold,
            "inference_promoted": qc_pass,
        },
        "failure_reason": None if qc_pass else "FROZEN_MAXIMUM_UNTESTED_FRACTION_EXCEEDED",
    }


def write_immutable(path: Path, payload: bytes) -> str:
    """Atomically create a decision; identical reruns are idempotent."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise ValueError(f"Refusing to overwrite an existing different decision: {path}")
        return canonical.sha256(path)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    temp = Path(temp_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temp, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise ValueError(f"Refusing to overwrite an existing different decision: {path}")
        return canonical.sha256(path)
    finally:
        temp.unlink(missing_ok=True)


def finalize(run_dir: Path, input_root: Path, reference_provenance: Path,
             family_path: Path = partial.DEFAULT_FAMILY,
             execution_path: Path = partial.DEFAULT_EXECUTION) -> tuple[dict[str, Any], str]:
    family = canonical.read_json(family_path)
    audit = partial.audit(run_dir, input_root, reference_provenance, family_path, execution_path)
    if audit["state"] != "COMPLETE":
        raise ValueError("Refusing to write a family decision before all frozen cells are audited")

    run_identity_path = run_dir / "run_identity.json"
    latest_audit_path = run_dir / "latest_audit.json"
    run_identity = canonical.read_json(run_identity_path)
    latest = canonical.read_json(latest_audit_path)
    verify_runner_identity(latest, run_identity, audit)
    if latest.get("audit", {}).get("qc_pass") != audit["family_qc_pass"]:
        raise ValueError("Runner audit QC decision differs from independent receipt audit")
    for key in ("intended_cells", "valid_cells", "statuses", "untested_fraction"):
        if latest.get("audit", {}).get(key) != {
            "intended_cells": audit["planned_cells"],
            "valid_cells": audit["verified_cells"],
            "statuses": audit["statuses"],
            "untested_fraction": audit["currently_uncompleted_or_untested_cells"] / audit["planned_cells"],
        }[key]:
            raise ValueError(f"Runner audit field {key} differs from independent receipt audit")

    loci = canonical.read_loci(
        ROOT / family["locus_definition"]["path"],
        family["locus_definition"]["sha256"],
        int(family["canonical_univariate"]["loci"]),
    )
    aggregate_path = run_dir / "results" / "canonical_family_results.tsv"
    recorded_aggregate = latest.get("aggregate", {})
    if Path(recorded_aggregate.get("path", "")).resolve() != aggregate_path.resolve():
        raise ValueError("Runner audit aggregate path differs from the canonical output")
    if recorded_aggregate.get("sha256") != canonical.sha256(aggregate_path):
        raise ValueError("Runner audit aggregate checksum does not match its output")
    aggregate_rows = verify_aggregate(aggregate_path, run_dir, loci, tuple(family["trait_ids"]))
    decision = build_decision(audit, family, run_identity, recorded_aggregate, aggregate_rows)
    decision_path = run_dir / "canonical_family_decision.json"
    digest = write_immutable(decision_path, json_bytes(decision))
    return decision, digest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference-provenance", type=Path, required=True)
    parser.add_argument("--family-lock", type=Path, default=partial.DEFAULT_FAMILY)
    parser.add_argument("--execution-lock", type=Path, default=partial.DEFAULT_EXECUTION)
    args = parser.parse_args()
    decision, digest = finalize(args.run_dir, args.input_root, args.reference_provenance,
                                args.family_lock, args.execution_lock)
    print(json.dumps({"decision": decision, "decision_sha256": digest}, indent=2, sort_keys=True))
    return 0 if decision["promotion_permitted"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
