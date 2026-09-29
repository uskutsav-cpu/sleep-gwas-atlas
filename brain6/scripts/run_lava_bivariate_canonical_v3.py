#!/usr/bin/env python3
"""Run the frozen five-pair LAVA v3 stage after canonical univariate QC passes.

The coordinator refuses to start before all 17,465 canonical tests have valid
receipts and pass their frozen missingness gate. Four locus-scoped R workers
consume those verified canonical rows; pairwise inputs and the LD reference
are loaded only for pairs that pass the canonical eligibility gate.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
sys.path.insert(0, str(ROOT / "brain6/scripts"))
import run_lava_canonical_v3 as canonical
from brain6.io import read_json, sha256
from brain6.stats import bh

LOCK = ROOT / "brain6/config/lava_bivariate_canonical_v3.json"
WORKER = ROOT / "brain6/scripts/run_lava_family_locus_v3.R"
TRAITS = canonical.TRAITS
PAIRS = ("insomnia__adhd", "insomnia__mdd", "longsleep__scz",
         "longsleep__bipolar", "longsleep__parkinson")
PAIR_STATUSES = {"TESTED", "UNIVARIATE_UNDERPOWERED", "NO_OVERLAP", "FAILED", "NOT_RUN"}
RESULT_FIELDS = ("pair_id", "locus_id", "status", "p", "local_rg", "reason")
CANONICAL_AGGREGATE_FIELDS = (
    "phen", "locus_id", "chromosome", "status", "n_snps", "n_components",
    "h2.obs", "h2.latent", "p", "reason",
)
PROTECTED_RUN_IDS = canonical.PROTECTED_RUN_IDS


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def atomic_json(path: Path, value: Any) -> None:
    payload = json.dumps(value, sort_keys=True, indent=2).encode() + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try: os.unlink(temporary)
        except FileNotFoundError: pass
        raise


def immutable_bytes(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Refusing to replace immutable bivariate input/output: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try: os.unlink(temporary)
        except FileNotFoundError: pass
        raise


def tsv_bytes(fields: list[str], records: list[dict[str, Any]]) -> bytes:
    import io
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(records)
    return stream.getvalue().encode()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def file_identity(path: Path) -> dict[str, Any]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def bivariate_summary(rows: list[dict[str, str]], maximum_failure_fraction: float) -> dict[str, Any]:
    if len(rows) != len(PAIRS) * 2495:
        raise ValueError("Bivariate family must preserve all 12,475 pair-locus slots")
    if any(row.get("status") not in PAIR_STATUSES for row in rows):
        raise ValueError("Bivariate family contains an unknown status")
    counts = {status: sum(row["status"] == status for row in rows) for status in sorted(PAIR_STATUSES)}
    # The frozen 1% limit is explicitly an execution-failure ceiling.
    # NO_OVERLAP is a valid scientific ineligibility status, like a gate miss;
    # it remains in the 12,475-slot BH denominator but must not consume the
    # execution-failure allowance.
    execution_failures = counts["FAILED"] + counts["NOT_RUN"]
    evidence_gaps = counts["NO_OVERLAP"] + execution_failures
    tested = [row for row in rows if row["status"] == "TESTED"]
    p_values: list[float | None] = []
    for row in rows:
        if row["status"] == "TESTED":
            p = float(row["p"]); rg = float(row["local_rg"])
            if not (0 <= p <= 1 and -1 <= rg <= 1):
                raise ValueError("Invalid tested bivariate estimate")
            p_values.append(p)
        else:
            p_values.append(None)
    q_values = bh(p_values, family_size=len(rows))
    result = []
    for row, q in zip(rows, q_values):
        result.append({**{field: row.get(field, "") for field in RESULT_FIELDS},
                       "family_fdr": q if row["status"] == "TESTED" else "NA"})
    fraction = execution_failures / len(rows)
    status = ("FAILED_QC_NOT_CONSUMED" if fraction > maximum_failure_fraction else
              "INSUFFICIENT_EVIDENCE" if evidence_gaps or not tested else "PASS")
    return {"rows": result, "status": status, "counts": counts,
            "tested_slots": len(tested), "execution_failures_or_missing": execution_failures,
            "execution_failure_fraction": fraction,
            "maximum_execution_failure_fraction": maximum_failure_fraction,
            "planned_slots": len(rows), "bh_family_size": len(rows)}


def canonical_rows_for_locus(path: Path, locus_id: str) -> list[dict[str, str]]:
    rows = [row for row in read_tsv(path) if row.get("locus_id") == locus_id]
    if len(rows) != len(TRAITS) or {row.get("phen") for row in rows} != set(TRAITS):
        raise ValueError(f"Canonical aggregate does not contain exactly seven cells at locus {locus_id}")
    if any(row["status"] not in canonical.VALID_STATUSES for row in rows):
        raise ValueError(f"Unknown canonical status at locus {locus_id}")
    return rows


def verify_canonical_aggregate(aggregate_path: Path, canonical_root: Path,
                               loci: list[dict[str, str]]) -> None:
    """Require every aggregate row and value to match its verified locus cells."""
    observed_rows = read_tsv(aggregate_path)
    expected_keys = [(trait, locus["LOC"]) for locus in loci for trait in TRAITS]
    observed_keys = [(row.get("phen", ""), row.get("locus_id", "")) for row in observed_rows]
    if len(observed_rows) != len(expected_keys) or observed_keys != expected_keys:
        raise ValueError("Canonical aggregate does not preserve the exact ordered 17,465-cell family")
    if any(tuple(row.keys()) != CANONICAL_AGGREGATE_FIELDS for row in observed_rows):
        raise ValueError("Canonical aggregate columns differ from the frozen output schema")
    index = 0
    for locus in loci:
        cells_path = canonical.result_path(canonical_root, locus["LOC"])
        cells = canonical.load_cells(cells_path)
        if (len(cells) != len(TRAITS) or
            {row.get("phen") for row in cells} != set(TRAITS) or
            any(row.get("locus_id") != locus["LOC"] for row in cells)):
            raise ValueError(f"Canonical receipt cells are incomplete at locus {locus['LOC']}")
        by_trait = {row["phen"]: row for row in cells}
        for trait in TRAITS:
            cell = by_trait[trait]
            expected = {
                "phen": trait, "locus_id": locus["LOC"], "chromosome": locus["CHR"],
                "status": cell["status"], "n_snps": cell.get("n_snps", ""),
                "n_components": cell.get("n_components", ""), "h2.obs": cell.get("h2.obs", ""),
                "h2.latent": cell.get("h2.latent", ""), "p": cell.get("p", ""),
                "reason": cell.get("reason", ""),
            }
            if any(observed_rows[index].get(key, "") != value for key, value in expected.items()):
                raise ValueError(f"Canonical aggregate differs from verified cells at {trait}/{locus['LOC']}")
            index += 1


def validate_pair_coverage(rows: list[dict[str, str]], expected_keys: set[tuple[str, str]],
                           minimum_shared_variants: int) -> dict[tuple[str, str], int]:
    if minimum_shared_variants < 1:
        raise ValueError("Minimum shared-reference-variant threshold must be positive")
    coverage: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row.get("pair_id", ""), row.get("locus_id", ""))
        if key in coverage:
            raise ValueError(f"Duplicate pair-locus shared-variant coverage row: {key}")
        try:
            count = int(row["shared_reference_variants"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Invalid shared-reference-variant count at {key}") from error
        if count < 0:
            raise ValueError(f"Negative shared-reference-variant count at {key}")
        expected_status = "READY" if count >= minimum_shared_variants else "NO_OVERLAP_LT_MIN_K"
        if row.get("status") != expected_status:
            raise ValueError(f"Coverage status disagrees with variant count at {key}")
        coverage[key] = count
    if set(coverage) != expected_keys:
        raise ValueError("Coverage table does not exactly cover the locked pair-locus family")
    return coverage


def pair_gate_outcomes(pairs: list[dict[str, Any]], canonical_rows: list[dict[str, str]],
                       minimum_shared_variants: int, p_threshold: float
                       ) -> tuple[list[dict[str, str]], list[str]]:
    if minimum_shared_variants < 1 or not math.isfinite(p_threshold) or not 0 < p_threshold <= 1:
        raise ValueError("Pair gate thresholds must be finite and positive")
    if len(pairs) != len(PAIRS) or {str(pair.get("pair_id")) for pair in pairs} != set(PAIRS):
        raise ValueError("Pair gate requires the exact five unique locked phenotype pairs")
    by_trait = {row["phen"]: row for row in canonical_rows}
    if len(by_trait) != len(TRAITS) or set(by_trait) != set(TRAITS):
        raise ValueError("Pair gate requires the exact seven canonical trait rows")
    if any(row.get("status") not in canonical.VALID_STATUSES for row in canonical_rows):
        raise ValueError("Pair gate received an unknown canonical status")
    locus_ids = {str(row.get("locus_id", "")) for row in canonical_rows}
    if len(locus_ids) != 1 or not next(iter(locus_ids)):
        raise ValueError("Canonical pair-gate rows must identify one common locus")
    locus_id = next(iter(locus_ids))
    results: list[dict[str, str]] = []
    eligible: list[str] = []
    for pair in pairs:
        pair_id = str(pair["pair_id"])
        locus_id = str(pair["locus_id"])
        if locus_id not in locus_ids:
            raise ValueError(f"Pair/canonical locus mismatch: {pair_id}/{locus_id}")
        a, b = pair["phenotypes"]
        if a not in by_trait or b not in by_trait or a == b:
            raise ValueError(f"Invalid locked phenotype pair: {pair_id}")
        if (a, b) != tuple(pair_id.split("__", 1)):
            raise ValueError(f"Pair identifier does not match its locked phenotype tuple: {pair_id}")
        left, right = by_trait[a], by_trait[b]
        status, reason = "", ""
        if int(pair["shared_reference_variants"]) < minimum_shared_variants:
            status, reason = "NO_OVERLAP", "FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"
        elif left["status"] != "TESTED" or right["status"] != "TESTED":
            status, reason = "UNIVARIATE_UNDERPOWERED", "ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_NOT_RUN"
        else:
            try:
                p_left, p_right = float(left["p"]), float(right["p"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Invalid tested canonical P value for {pair_id}/{locus_id}") from error
            if not (0 <= p_left <= 1 and 0 <= p_right <= 1):
                raise ValueError(f"Out-of-range tested canonical P value for {pair_id}/{locus_id}")
            if not (p_left < p_threshold and p_right < p_threshold):
                status, reason = "UNIVARIATE_UNDERPOWERED", "ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_FAILED_FROZEN_GATE"
            else:
                eligible.append(pair_id)
                continue
        results.append({"pair_id": pair_id, "locus_id": locus_id, "status": status,
                        "p": "NA", "local_rg": "NA", "reason": reason})
    return results, eligible


def write_gate_only_unit(target: Path, locus_id: str, bivariate_run_id: str,
                         lock_sha: str, canonical_run_id: str, canonical_receipt_sha: str,
                         canonical_digest: str, canonical_rows: list[dict[str, str]],
                         pair_results: list[dict[str, str]], pair_gate_inputs: list[dict[str, Any]],
                         minimum_shared_variants: int, p_threshold: float,
                         worker_script_sha: str, maximum_attempts: int) -> dict[str, Any]:
    if len(pair_results) != len(PAIRS) or {row["pair_id"] for row in pair_results} != set(PAIRS):
        raise ValueError("Gate-only execution must preserve all five pair-locus slots")
    if any(str(row.get("locus_id")) != locus_id for row in pair_results):
        raise ValueError("Gate-only pair results must all belong to their receipt locus")
    expected_results, eligible_pairs = pair_gate_outcomes(
        pair_gate_inputs, canonical_rows, minimum_shared_variants, p_threshold)
    if eligible_pairs or pair_results != expected_results:
        raise ValueError("Gate-only receipt is permitted only when every frozen pair fails its computed gate")
    attempts_root = target.parent / ".attempts"
    attempts_root.mkdir(parents=True, exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix=f".{locus_id}.gate-only-", dir=attempts_root))
    output = attempt / "worker_output"
    output.mkdir()
    univariate_fields = ["phen", "locus_id", "status", "h2.obs", "h2.latent", "p", "reason"]
    univariate_rows = [{field: row.get(field, "NA") if row.get(field, "") != "" else "NA"
                        for field in univariate_fields} for row in canonical_rows]
    immutable_bytes(output / "pair_results.tsv", tsv_bytes(list(RESULT_FIELDS), pair_results))
    immutable_bytes(output / "univariate.tsv", tsv_bytes(univariate_fields, univariate_rows))
    n_underpowered = sum(row["status"] == "UNIVARIATE_UNDERPOWERED" for row in pair_results)
    status = {"status": "PASS", "analysis_id": "brain6-lava-local-rg-v3",
        "bivariate_run_id": bivariate_run_id, "canonical_run_id": canonical_run_id,
        "canonical_cell_sha256": canonical_digest,
        "canonical_receipt_sha256": canonical_receipt_sha,
        "locus_id": locus_id, "n_pairs": len(pair_results), "n_tested": 0,
        "n_underpowered": n_underpowered, "n_failed": 0,
        "unique_univariate_tested": sum(row["status"] == "TESTED" for row in canonical_rows),
        "unique_univariate_not_run": sum(row["status"] != "TESTED" for row in canonical_rows),
        "participant_level_overlap_verified": False, "lava_version": "0.1.5",
        "execution_mode": "PYTHON_GATE_ONLY"}
    atomic_json(output / "status.json", status)
    config = {"schema_version": 1, "analysis_id": "brain6-lava-local-rg-v3",
        "family_lock_sha256": lock_sha, "execution_lock_sha256": lock_sha,
        "locus_id": locus_id, "run_id": bivariate_run_id,
        "bivariate_run_id": bivariate_run_id,
        "canonical_run_id": canonical_run_id, "canonical_receipt_sha256": canonical_receipt_sha,
        "canonical_cell_sha256": canonical_digest, "canonical_univariate_rows": canonical_rows,
        "pair_gate_inputs": pair_gate_inputs,
        "minimum_shared_reference_variants": minimum_shared_variants,
        "strict_canonical_p_threshold": p_threshold,
        "execution_mode": "PYTHON_GATE_ONLY", "r_process_attempts": 0,
        "maximum_r_process_attempts_if_eligible": maximum_attempts}
    config_path = attempt / "worker_config.json"
    atomic_json(config_path, config)
    outputs = [{"path": str(path.relative_to(attempt)), "bytes": path.stat().st_size,
                "sha256": sha256(path)} for path in sorted(output.iterdir()) if path.is_file()]
    receipt = {"schema_version": 1, "locus_id": locus_id, "run_id": bivariate_run_id,
        "canonical_run_id": canonical_run_id, "family_lock_sha256": lock_sha,
        "execution_lock_sha256": lock_sha, "canonical_receipt_sha256": canonical_receipt_sha,
        "canonical_cell_sha256": canonical_digest, "worker_script_sha256": worker_script_sha,
        "worker_config_sha256": sha256(config_path), "execution_mode": "PYTHON_GATE_ONLY",
        "outputs": outputs, "status": "COMPLETE"}
    atomic_json(attempt / "receipt.json", receipt)
    if target.exists():
        raise FileExistsError(f"Bivariate target already exists without a valid receipt: {target}")
    os.rename(attempt, target)
    return {"locus_id": locus_id, "resumed": False, "attempts": 0,
            "execution_mode": "PYTHON_GATE_ONLY", "worker_log_sha256": ""}


def validate_lock(path: Path) -> tuple[dict[str, Any], str]:
    lock_sha = sha256(path)
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if not sidecar.is_file() or sidecar.read_text().split()[0] != lock_sha:
        raise ValueError("Bivariate v3 lock is missing or its checksum sidecar disagrees")
    lock = read_json(path)
    if lock.get("analysis_id") != "brain6-lava-local-rg-v3" or lock.get("pair_family", {}).get("n_slots") != 12475:
        raise ValueError("Bivariate lock identity or denominator changed")
    if lock.get("pair_family", {}).get("pair_ids") != list(PAIRS):
        raise ValueError("Bivariate pair family order differs from the reviewed lock")
    if lock.get("execution", {}).get("worker_count") != 4:
        raise ValueError("Bivariate execution lock must use the supported four-worker setting")
    expected_protocol = (
        "Python coordinator writes gate-only receipts without R when no pair passes the frozen "
        "univariate and coverage gates; up to four one-locus R workers run only for eligible loci")
    if lock.get("execution", {}).get("worker_protocol") != expected_protocol:
        raise ValueError("Bivariate execution lock differs from the validated gate-first scheduler")
    if (lock.get("pair_family", {}).get("n_loci") != 2495 or
        lock.get("pair_family", {}).get("maximum_univariate_untested_fraction") != 0.05 or
        lock.get("pair_family", {}).get("maximum_execution_failure_fraction") != 0.01):
        raise ValueError("Bivariate family size or frozen QC limits changed")
    if lock.get("provenance", {}).get("worker_script_sha256") != sha256(WORKER):
        raise ValueError("Bivariate worker source changed after the family lock was written")
    source_exec = ROOT / lock["provenance"]["source_execution_policy_path"]
    source_family = ROOT / lock["provenance"]["source_family_policy_path"]
    if sha256(source_exec) != lock["provenance"]["source_execution_policy_sha256"]:
        raise ValueError("Inherited bivariate parameter source checksum changed")
    if sha256(source_family) != lock["provenance"]["source_family_policy_sha256"]:
        raise ValueError("Inherited family QC source checksum changed")
    if read_json(source_exec)["bivariate"] != lock["bivariate_parameters"]:
        raise ValueError("Bivariate method parameters differ from their frozen pre-v3 source")
    if read_json(source_family)["bivariate_family"]["n_slots"] != lock["pair_family"]["n_slots"]:
        raise ValueError("Bivariate family denominator differs from the frozen source policy")
    return lock, lock_sha


def verify_canonical_source(lock: dict[str, Any], canonical_root: Path,
                            family_path: Path, execution_path: Path,
                            input_root: Path, reference_root: Path) -> tuple[dict, list[dict[str, str]], dict]:
    family = read_json(family_path)
    expected = lock["canonical_source"]
    if family["analysis_id"] != expected["analysis_id"] or sha256(family_path) != expected["family_lock_sha256"]:
        raise ValueError("Canonical family lock differs from the frozen bivariate source")
    if sha256(execution_path) != expected["execution_lock_sha256"]:
        raise ValueError("Canonical execution lock differs from the frozen bivariate source")
    identity = read_json(canonical_root / "run_identity.json")
    if identity.get("run_id") != expected["run_id"] or identity.get("run_identity_sha256") != expected["run_id"]:
        raise ValueError("Canonical output root is not the bound production run")
    if identity.get("family_lock_sha256") != expected["family_lock_sha256"] or identity.get("execution_lock_sha256") != expected["execution_lock_sha256"]:
        raise ValueError("Canonical run identity does not bind the locked source files")
    latest_path = canonical_root / "latest_audit.json"
    if not latest_path.is_file():
        raise RuntimeError("Canonical v3 run has no completed full-family audit yet")
    latest = read_json(latest_path)
    if latest.get("run_id") != expected["run_id"] or latest.get("audit", {}).get("qc_pass") is not True:
        raise RuntimeError("Canonical v3 family audit has not passed its frozen QC gate")
    loci_path = ROOT / family["locus_definition"]["path"]
    loci = canonical.read_loci(loci_path, family["locus_definition"]["sha256"], 2495)
    reference_prov = reference_root / "reference.provenance.json"
    input_prov_path = input_root / "provenance.json"
    if sha256(reference_prov) != lock["provenance"]["official_reference_provenance_sha256"]:
        raise ValueError("Official reference provenance differs from the bivariate lock")
    if sha256(input_prov_path) != lock["provenance"]["input_provenance_sha256"]:
        raise ValueError("Materialized input provenance differs from the bivariate lock")
    if lock["pair_family"]["pair_lock_sha256"] != family["pair_lock_sha256"]:
        raise ValueError("Bivariate pair selection lock differs from the canonical family")
    if lock["provenance"]["overlap_matrices"] != family["overlap_gate"]["pairwise_sample_overlap_file_manifest"]["matrices"]:
        raise ValueError("Bivariate overlap covariance matrix manifest differs from the canonical family")
    for matrix in lock["provenance"]["overlap_matrices"]:
        matrix_path = ROOT / matrix["path"]
        if not matrix_path.is_file() or sha256(matrix_path) != matrix["sha256"]:
            raise ValueError(f"Pairwise covariance matrix checksum mismatch: {matrix['pair_id']}")
    audit = canonical.audit_family(canonical_root, loci, tuple(TRAITS), family["analysis_id"],
        expected["family_lock_sha256"], expected["execution_lock_sha256"], input_root,
        reference_prov, identity["script_sha256"])
    if audit["qc_pass"] is not True or audit["valid_cells"] != 17465 or audit["invalid_loci"]:
        raise RuntimeError("Recomputed canonical audit fails coverage or frozen 5% untested QC")
    aggregate = canonical_root / "results/canonical_family_results.tsv"
    if not aggregate.is_file() or latest.get("aggregate", {}).get("sha256") != sha256(aggregate):
        raise ValueError("Canonical family aggregate is missing or differs from its completed audit")
    verify_canonical_aggregate(aggregate, canonical_root, loci)
    input_prov = read_json(input_prov_path)
    coverage_meta = input_prov["pair_coverage"]
    coverage_path = Path(coverage_meta["path"])
    if not coverage_path.is_file() or sha256(coverage_path) != coverage_meta["sha256"]:
        raise ValueError("Pair-locus shared-variant coverage table fails its input provenance hash")
    coverage_rows = read_tsv(coverage_path)
    minimum_shared = int(lock["eligibility"]["minimum_shared_reference_variants"])
    if minimum_shared != int(lock["locus_processing"]["min_K"]):
        raise ValueError("Pair eligibility and LAVA min_K thresholds disagree")
    coverage = validate_pair_coverage(coverage_rows,
        {(pair, locus["LOC"]) for pair in PAIRS for locus in loci}, minimum_shared)
    return family, loci, {"identity": identity, "run_id": identity["run_id"], "audit": audit, "aggregate": aggregate,
                          "aggregate_sha256": sha256(aggregate), "input_provenance": input_prov,
                          "coverage": coverage, "reference_provenance": reference_prov}


def verify_unit(target: Path, locus_id: str, run_id: str, family_sha: str,
                execution_sha: str, canonical_run_id: str,
                canonical_receipt_sha: str) -> dict[str, Any] | None:
    receipt_path = target / "receipt.json"
    if not receipt_path.is_file():
        return None
    receipt = read_json(receipt_path)
    if (receipt.get("locus_id") != locus_id or receipt.get("run_id") != run_id or
        receipt.get("family_lock_sha256") != family_sha or
        receipt.get("execution_lock_sha256") != execution_sha or
        receipt.get("canonical_run_id") != canonical_run_id or
        receipt.get("canonical_receipt_sha256") != canonical_receipt_sha):
        raise ValueError(f"Bivariate checkpoint identity differs: {locus_id}")
    config_path = target / "worker_config.json"
    if not config_path.is_file() or sha256(config_path) != receipt.get("worker_config_sha256"):
        raise ValueError(f"Bivariate worker config checksum mismatch: {locus_id}")
    config = read_json(config_path)
    mode = receipt.get("execution_mode")
    if (mode not in {"R_BIVARIATE", "PYTHON_GATE_ONLY"} or
        config.get("execution_mode") != mode or config.get("locus_id") != locus_id or
        config.get("analysis_id") != "brain6-lava-local-rg-v3" or
        config.get("family_lock_sha256") != family_sha or
        config.get("execution_lock_sha256") != execution_sha or
        config.get("canonical_run_id") != canonical_run_id or
        config.get("canonical_receipt_sha256") != canonical_receipt_sha or
        config.get("canonical_cell_sha256") != receipt.get("canonical_cell_sha256") or
        config.get("bivariate_run_id", config.get("run_id")) != run_id):
        raise ValueError(f"Bivariate execution mode is missing or differs from its worker config: {locus_id}")
    output_map = {item["path"]: item for item in receipt.get("outputs", [])}
    if set(output_map) != {"worker_output/pair_results.tsv", "worker_output/univariate.tsv", "worker_output/status.json"}:
        raise ValueError(f"Bivariate checkpoint output inventory differs: {locus_id}")
    for rel, item in output_map.items():
        path = target / rel
        if not path.is_file() or path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Bivariate checkpoint output checksum mismatch: {path}")
    results = read_tsv(target / "worker_output/pair_results.tsv")
    if (len(results) != 5 or {r["pair_id"] for r in results} != set(PAIRS) or
        any(r["status"] not in PAIR_STATUSES or str(r.get("locus_id")) != locus_id for r in results)):
        raise ValueError(f"Bivariate checkpoint does not contain the exact five pair slots: {locus_id}")
    if mode == "PYTHON_GATE_ONLY":
        try:
            expected_results, eligible_pairs = pair_gate_outcomes(
                config["pair_gate_inputs"], config["canonical_univariate_rows"],
                int(config["minimum_shared_reference_variants"]),
                float(config["strict_canonical_p_threshold"]))
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Gate-only worker config does not reproduce its frozen gate: {locus_id}") from error
        if eligible_pairs or config.get("r_process_attempts") != 0 or results != expected_results:
            raise ValueError(f"Gate-only receipt contains an eligible or inconsistent pair: {locus_id}")
    status = read_json(target / "worker_output/status.json")
    if (status.get("status") != "PASS" or
        status.get("analysis_id") != "brain6-lava-local-rg-v3" or
        str(status.get("locus_id")) != locus_id or
        status.get("canonical_cell_sha256") != receipt.get("canonical_cell_sha256") or
        status.get("bivariate_run_id") != run_id or
        status.get("canonical_run_id") != canonical_run_id or
        status.get("canonical_receipt_sha256") != canonical_receipt_sha):
        raise ValueError(f"Bivariate worker status identity differs: {locus_id}")
    if mode == "PYTHON_GATE_ONLY" and status.get("execution_mode") != mode:
        raise ValueError(f"Gate-only worker status omits its execution mode: {locus_id}")
    return receipt


def write_worker_unit(target: Path, locus: dict[str, str], family: dict,
                      lock: dict, lock_sha: str, bivariate_run_id: str, canonical_root: Path,
                      source: dict, input_root: Path, reference_root: Path,
                      pair_lock: dict, rscript: Path, runtime_sha: str) -> dict[str, Any]:
    locus_id, chrom = locus["LOC"], int(locus["CHR"])
    canonical_receipt = canonical.receipt_path(canonical_root, locus_id)
    if not canonical_receipt.is_file():
        raise ValueError(f"Canonical locus receipt is missing: {locus_id}")
    existing = verify_unit(target, locus_id, bivariate_run_id, lock_sha,
                           lock_sha, source["run_id"], sha256(canonical_receipt))
    if existing:
        return {"locus_id": locus_id, "resumed": True, "attempts": 0,
                "execution_mode": existing["execution_mode"], "worker_log_sha256": ""}
    if (target / "receipt.json").exists():
        raise ValueError(f"Existing bivariate receipt is invalid; refusing overwrite: {target}")
    canonical_file = canonical.result_path(canonical_root, locus_id)
    if sha256(canonical_file) != read_json(canonical_receipt).get("output", {}).get("sha256"):
        raise ValueError(f"Canonical locus output changed after source audit: {locus_id}")
    canonical_rows = canonical_rows_for_locus(canonical_file, locus_id)
    pair_metadata = []
    for pair_id in PAIRS:
        locked_pair = pair_lock[pair_id]
        pair_metadata.append({"pair_id": pair_id,
            "phenotypes": [locked_pair["sleep_trait"], locked_pair["disease_trait"]],
            "locus_id": locus_id,
            "shared_reference_variants": source["coverage"][(pair_id, locus_id)]})
    policy = {"locus_processing": lock["locus_processing"],
        "univariate": {"cap_estimates": True, "gate_p_strictly_less_than": lock["eligibility"]["strict_p_threshold"]},
        "bivariate": lock["bivariate_parameters"]}
    canonical_digest = sha256(canonical_file)
    canonical_receipt_sha = sha256(canonical_receipt)
    gate_results, eligible_pairs = pair_gate_outcomes(pair_metadata, canonical_rows,
        int(lock["eligibility"]["minimum_shared_reference_variants"]),
        float(lock["eligibility"]["strict_p_threshold"]))
    if not eligible_pairs:
        return write_gate_only_unit(target, locus_id, bivariate_run_id, lock_sha,
            source["run_id"], canonical_receipt_sha, canonical_digest, canonical_rows,
            gate_results, pair_metadata,
            int(lock["eligibility"]["minimum_shared_reference_variants"]),
            float(lock["eligibility"]["strict_p_threshold"]),
            lock["provenance"]["worker_script_sha256"],
            int(lock["execution"]["max_process_attempts_per_locus"]))
    pair_matrix = {row["pair_id"]: row for row in family["overlap_gate"]["pairwise_sample_overlap_file_manifest"]["matrices"]}
    pairs = []
    for pair in pair_metadata:
        pair_id = pair["pair_id"]
        locked_pair = pair_lock[pair_id]
        matrix = pair_matrix[pair_id]
        matrix_path = ROOT / matrix["path"]
        pairs.append({**pair,
            "phenotypes": [locked_pair["sleep_trait"], locked_pair["disease_trait"]],
            "input_info": str(input_root / f"locus_{locus_id}/input_info.tsv"),
            "sample_overlap_file": str(matrix_path),
            "reference_prefix": str(reference_root / f"lava-ukb-v1.1_chr{chrom}")})
    attempts_root = target.parent / ".attempts"
    attempts_root.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt_no in range(1, int(lock["execution"]["max_process_attempts_per_locus"]) + 1):
        attempt = Path(tempfile.mkdtemp(prefix=f".{locus_id}.attempt-{attempt_no}-", dir=attempts_root))
        out = attempt / "worker_output"; out.mkdir()
        cfg = {"schema_version": 1, "analysis_id": lock["analysis_id"],
            "execution_mode": "R_BIVARIATE",
            "family_lock_sha256": lock_sha, "execution_lock_sha256": lock_sha,
            "run_id": bivariate_run_id, "bivariate_run_id": bivariate_run_id,
            "canonical_run_id": source["run_id"], "canonical_cell_sha256": canonical_digest,
            "canonical_receipt_sha256": canonical_receipt_sha,
            "canonical_univariate_rows": canonical_rows,
            "runtime_validator_passed": True, "runtime_validation_output_sha256": runtime_sha,
            "locus_id": locus_id, "loci_file": str(ROOT / family["locus_definition"]["path"]),
            "random_seed": int(lock["execution"].get("random_seed", 20260922)),
            "execution_policy": policy, "output_dir": str(out),
            "trait_ids": list(TRAITS), "pairs": pairs}
        cfg_path = attempt / "worker_config.json"
        atomic_json(cfg_path, cfg)
        proc = subprocess.run([str(rscript), str(WORKER), str(cfg_path)], cwd=ROOT,
            capture_output=True, text=True, check=False)
        log_path = attempt / "worker.log"
        log_path.write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode == 0:
            results = read_tsv(out / "pair_results.tsv") if (out / "pair_results.tsv").is_file() else []
            univariate = read_tsv(out / "univariate.tsv") if (out / "univariate.tsv").is_file() else []
            status = read_json(out / "status.json") if (out / "status.json").is_file() else {}
            actual_univariate = {row.get("phen"): row for row in univariate}
            expected_univariate = {row["phen"]: row for row in canonical_rows}
            univariate_matches = len(actual_univariate) == 7 and set(actual_univariate) == set(TRAITS)
            if univariate_matches:
                for trait in TRAITS:
                    actual, expected = actual_univariate[trait], expected_univariate[trait]
                    if actual.get("status") != expected.get("status") or actual.get("reason", "") != expected.get("reason", ""):
                        univariate_matches = False; break
                    for field in ("h2.obs", "h2.latent", "p"):
                        av, ev = actual.get(field, "NA"), expected.get(field, "NA")
                        if av in {"", "NA", "NaN"} or ev in {"", "NA", "NaN"}:
                            if av not in {"", "NA", "NaN"} or ev not in {"", "NA", "NaN"}:
                                univariate_matches = False; break
                        elif abs(float(av)-float(ev)) > 1e-12:
                            univariate_matches = False; break
            univariate_matches = univariate_matches and status.get("canonical_receipt_sha256") == canonical_receipt_sha
            if (len(results) == 5 and {r.get("pair_id") for r in results} == set(PAIRS) and
                all(r.get("status") in PAIR_STATUSES for r in results) and
                univariate_matches and status.get("bivariate_run_id") == bivariate_run_id and
                status.get("canonical_run_id") == source["run_id"] and
                status.get("canonical_cell_sha256") == canonical_digest):
                outputs = []
                for item in sorted(out.iterdir()):
                    if item.is_file():
                        outputs.append({"path": str(item.relative_to(attempt)), "bytes": item.stat().st_size,
                                        "sha256": sha256(item)})
                receipt = {"schema_version": 1, "locus_id": locus_id, "run_id": bivariate_run_id,
                    "canonical_run_id": source["run_id"],
                    "family_lock_sha256": lock_sha, "execution_lock_sha256": lock_sha,
                    "canonical_receipt_sha256": canonical_receipt_sha,
                    "canonical_cell_sha256": canonical_digest,
                    "worker_script_sha256": lock["provenance"]["worker_script_sha256"],
                    "worker_config_sha256": sha256(cfg_path), "execution_mode": "R_BIVARIATE",
                    "outputs": outputs, "status": "COMPLETE"}
                atomic_json(attempt / "receipt.json", receipt)
                if target.exists():
                    raise FileExistsError(f"Bivariate target already exists without a receipt: {target}")
                os.rename(attempt, target)
                return {"locus_id": locus_id, "resumed": False, "attempts": attempt_no,
                        "execution_mode": "R_BIVARIATE",
                        "worker_log_sha256": sha256(target / "worker.log")}
            proc = subprocess.CompletedProcess(proc.args, 1, proc.stdout,
                "Worker output failed exact five-pair/seven-canonical-row/source identity validation")
        last_error = {"returncode": proc.returncode, "stdout_tail": proc.stdout[-3000:],
                      "stderr_tail": proc.stderr[-3000:]}
        atomic_json(attempt / "failure.json", {"locus_id": locus_id, "attempt": attempt_no,
            **last_error, "worker_log_sha256": sha256(log_path)})
    raise RuntimeError(f"LAVA bivariate worker failed at locus {locus_id}: {last_error}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-run-dir", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--rscript", type=Path, required=True)
    parser.add_argument("--lock", type=Path, default=LOCK)
    parser.add_argument("--family-lock", type=Path, default=ROOT / "brain6/config/lava_family_canonical_v3.json")
    parser.add_argument("--canonical-execution-lock", type=Path, default=ROOT / "brain6/config/lava_execution_canonical_v3.json")
    parser.add_argument("--pair-lock", type=Path, default=ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json")
    args = parser.parse_args()
    lock, lock_sha = validate_lock(args.lock)
    if not args.rscript.is_file(): raise FileNotFoundError(args.rscript)
    if not args.input_root.is_dir() or not args.reference_root.is_dir():
        raise FileNotFoundError("Materialized input/reference root is missing")
    output_root = ROOT / lock["execution"]["output_root"]
    if not canonical.isolated_output_root(output_root) or output_root.name != "lava-local-rg-v3":
        raise ValueError("Bivariate output root is not isolated from preserved v1/v2/roundoff runs")
    reference_prov = args.reference_root / "reference.provenance.json"
    family, loci, source = verify_canonical_source(lock, args.canonical_run_dir,
        args.family_lock, args.canonical_execution_lock, args.input_root, args.reference_root)
    if len(loci) != 2495 or len(PAIRS) != 5:
        raise ValueError("Bivariate family is not the exact frozen five-pair by 2,495-locus design")

    validator = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
    runtime = subprocess.run([str(args.rscript), str(validator)], cwd=ROOT,
        capture_output=True, text=True, check=False)
    if runtime.returncode != 0 or runtime.stdout.count("TRACK_B_LAVA_RUNTIME\tstatus=PASS\t") != 1:
        raise RuntimeError("Pinned R/LAVA runtime preflight failed: " + (runtime.stdout + runtime.stderr)[-3000:])
    reference_check = subprocess.run([sys.executable, "scripts/lava_contract.py", "--verify-reference"],
        cwd=ROOT, capture_output=True, text=True, check=False)
    if reference_check.returncode != 0 or "LAVA_REFERENCE_VALIDATED" not in reference_check.stdout:
        raise RuntimeError("Official LAVA reference preflight failed: " + (reference_check.stdout + reference_check.stderr)[-3000:])
    runtime_sha = hashlib.sha256(runtime.stdout.encode()).hexdigest()
    identity = {"analysis_id": lock["analysis_id"], "bivariate_lock_sha256": lock_sha,
        "canonical_run_id": source["identity"]["run_id"],
        "canonical_aggregate_sha256": source["aggregate_sha256"],
        "canonical_audit": source["audit"], "input_provenance_sha256": sha256(args.input_root / "provenance.json"),
        "official_reference_provenance_sha256": sha256(reference_prov),
        "locus_definition_sha256": sha256(ROOT / family["locus_definition"]["path"]),
        "coordinator_sha256": sha256(Path(__file__)), "worker_sha256": sha256(WORKER),
        "runtime_validator_sha256": sha256(validator), "runtime_output_sha256": runtime_sha,
        "reference_validation_output_sha256": hashlib.sha256(reference_check.stdout.encode()).hexdigest(),
        "rscript_path": str(args.rscript.resolve()), "rscript_sha256": sha256(args.rscript.resolve()),
        "R": "4.3.3", "LAVA": "0.1.5", "reference_build": "GRCh37/hg19"}
    run_id = hashlib.sha256(canonical_json(identity)).hexdigest()
    run_dir = output_root / run_id
    output_root.mkdir(parents=True, exist_ok=True)
    lockfile = output_root / f".{run_id}.lock"
    try: fd = os.open(lockfile, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error: raise RuntimeError(f"Bivariate run is already locked: {lockfile}") from error
    try:
        with os.fdopen(fd, "w") as handle: json.dump({"pid": os.getpid(), "run_id": run_id}, handle)
        run_dir.mkdir(parents=True, exist_ok=True)
        identity_doc = {"analysis_id": lock["analysis_id"], "run_id": run_id,
                        "run_identity_sha256": run_id, "identity": identity}
        identity_path = run_dir / "run_identity.json"
        if identity_path.exists() and read_json(identity_path) != identity_doc:
            raise ValueError("Existing bivariate run identity differs; refusing resume")
        if not identity_path.exists(): atomic_json(identity_path, identity_doc)
        immutable_bytes(run_dir / "bivariate_lock.json", args.lock.read_bytes())
        planned = [{"pair_id": pair, "locus_id": locus["LOC"]} for locus in loci for pair in PAIRS]
        immutable_bytes(run_dir / "planned_units.tsv", tsv_bytes(["pair_id", "locus_id"], planned))
        pair_lock_doc = read_json(args.pair_lock)
        pair_lock_body = {key: value for key, value in pair_lock_doc.items() if key != "lock_sha256"}
        pair_lock_sha = hashlib.sha256(canonical_json(pair_lock_body)).hexdigest()
        if pair_lock_sha != lock["pair_family"]["pair_lock_sha256"]:
            raise ValueError("Reviewed pair selection lock differs from the bivariate family binding")
        pair_lock = {row["pair_id"]: row for row in pair_lock_doc["pairs"] if row.get("pair_id") in PAIRS}
        if set(pair_lock) != set(PAIRS): raise ValueError("Reviewed pair lock differs from the frozen five pairs")
        units = run_dir / "loci"; units.mkdir(exist_ok=True)
        started = time.monotonic(); before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(write_worker_unit, units / locus["LOC"], locus, family,
                lock, lock_sha, run_id, args.canonical_run_dir, source, args.input_root,
                args.reference_root, pair_lock, args.rscript, runtime_sha) for locus in loci]
            outcomes = []
            for future in concurrent.futures.as_completed(futures):
                outcomes.append(future.result())
                completed = len(outcomes)
                if completed % 100 == 0 or completed == len(loci):
                    print(f"BRAIN6_LAVA_BIVARIATE_PROGRESS loci={completed}/2495", flush=True)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        all_rows = []
        receipt_rows = []
        for locus in loci:
            loc_id = locus["LOC"]
            receipt = canonical.receipt_path(args.canonical_run_dir, loc_id)
            verified = verify_unit(units / loc_id, loc_id, run_id,
                lock_sha, lock_sha, source["identity"]["run_id"], sha256(receipt))
            if verified is None: raise ValueError(f"Completed bivariate locus lacks receipt: {loc_id}")
            receipt_rows.append({"locus_id": loc_id, "receipt_path": str(units / loc_id / "receipt.json"),
                                 "receipt_sha256": sha256(units / loc_id / "receipt.json")})
            all_rows.extend(read_tsv(units / loc_id / "worker_output/pair_results.tsv"))
        if len(all_rows) != 12475: raise ValueError("Bivariate output denominator differs from 12,475")
        ordered = {(r["pair_id"], r["locus_id"]): r for r in all_rows}
        if len(ordered) != len(all_rows): raise ValueError("Duplicate bivariate pair-locus output")
        rows = [ordered[(pair, locus["LOC"])] for locus in loci for pair in PAIRS]
        summary = bivariate_summary(rows, lock["pair_family"]["maximum_execution_failure_fraction"])
        results_path = run_dir / "results/bivariate_family_results.tsv"
        fields = [*RESULT_FIELDS, "family_fdr"]
        immutable_bytes(results_path, tsv_bytes(fields, summary["rows"]))
        receipts_path = run_dir / "worker_receipts.tsv"
        immutable_bytes(receipts_path, tsv_bytes(["locus_id", "receipt_path", "receipt_sha256"], receipt_rows))
        decision = {k: v for k, v in summary.items() if k != "rows"}
        decision.update({"analysis_id": lock["analysis_id"], "run_id": run_id,
            "bivariate_lock_sha256": lock_sha, "canonical_run_id": source["identity"]["run_id"],
            "canonical_aggregate_sha256": source["aggregate_sha256"],
            "participant_level_overlap_verified": False,
            "results_sha256": sha256(results_path), "receipt_index_sha256": sha256(receipts_path)})
        decision_path = run_dir / "family_decision.json"
        immutable_bytes(decision_path, (json.dumps(decision, sort_keys=True, indent=2) + "\n").encode())
        report = {"run_id": run_id, "workers": 4, "loci": 2495,
            "wall_seconds": time.monotonic()-started,
            "child_cpu_seconds": (after.ru_utime-before.ru_utime)+(after.ru_stime-before.ru_stime),
            "execution_modes": dict(Counter(row["execution_mode"] for row in outcomes)),
            "r_loci_with_eligible_pairs": sum(row["execution_mode"] == "R_BIVARIATE" for row in outcomes),
            "r_process_attempts": sum(row["attempts"] for row in outcomes if row["execution_mode"] == "R_BIVARIATE"),
            "gate_only_loci": sum(row["execution_mode"] == "PYTHON_GATE_ONLY" for row in outcomes),
            "decision": decision}
        atomic_json(run_dir / "latest_audit.json", report)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if summary["status"] != "FAILED_QC_NOT_CONSUMED" else 2
    finally:
        lockfile.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
