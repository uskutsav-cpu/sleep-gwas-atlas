#!/usr/bin/env python3
"""Export the complete failed-QC Track B LAVA family as diagnostics only.

This is deliberately not a second canonical-publication route.  It accepts
only a deeply validated V2 ``TERMINAL_FAILED_QC`` READY bundle, retains every
predeclared pair-locus state, and refuses to operate if the success-only LAVA
provenance file exists.  The requested tables are useful for audit and for
independent downstream methods, but they never attest that LAVA passed its
frozen scientific QC gates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "scripts/134_track_b_lava_continuation_contract.py"
SUPERVISOR_PATH = ROOT / "scripts/136_run_track_b_lava_continuation.py"
RESULT_VALIDATOR_PATH = ROOT / "scripts/138_validate_track_b_lava_results_v2.py"

OUTPUT_DIR = ROOT / "results/track_b/lava"
SUMMARY_PATH = ROOT / "results/track_b/LAVA_DISCOVERY_SUMMARY.md"
FORBIDDEN_CANONICAL_PROVENANCE = ROOT / "results/track_b/local/lava_results.provenance.json"

EXPORT_SCHEMA = "track-b-lava-terminal-diagnostic-export.1"
SCIENTIFIC_STATE = "TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY"
FINGERPRINT_RE = re.compile(r"^[0-9a-f]{64}$")
NA_VALUES = {"", "NA", "NaN", "nan"}

UNIV_EXTRA_FIELDS = [
    "correction_only_family_p", "correction_only_role", "terminal_scientific_state",
]
BIVARIATE_GRID_FIELDS = [
    "LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2",
    "locus_status", "trait1_univariate_status", "trait2_univariate_status",
    "trait1_local_h2_obs", "trait2_local_h2_obs", "both_traits_local_h2_eligible",
    "analysis_status", "eligibility_reason", "local_covariance", "rho", "rho.lower",
    "rho.upper", "r2", "r2.lower", "r2.upper", "p",
    "correction_only_family_p", "correction_only_role", "bivariate_test_family_n",
    "predeclared_pair_locus_family_n", "p_fdr", "fdr_significant", "n_snps", "K",
    "error", "numerical_stability", "terminal_scientific_state",
]
BIVARIATE_FDR_EXTRA_FIELDS = [
    "correction_only_family_p", "correction_only_role", "multiple_testing_rule",
    "terminal_scientific_state",
]
FAILURE_FIELDS = [
    "record_class", "stage", "pair_id", "LOC", "CHR", "START", "STOP",
    "phenotype", "conditional_model_id", "analysis_status", "error",
    "correction_treatment", "terminal_scientific_state",
]
ACCOUNTING_FIELDS = [
    "family", "status", "count", "family_total", "multiple_testing_denominator",
    "multiple_testing_rule", "threshold_name", "threshold", "observed_fraction",
    "gate_state", "notes", "terminal_scientific_state",
]
CANDIDATE_FIELDS = [
    "locus_index", "locus", "chromosome", "start", "stop", "n_snps",
    "eligible_models", "eligible_pairs",
]


class ExportError(RuntimeError):
    """A fail-closed diagnostic-export violation."""


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ExportError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


CONTRACT = _load_module("track_b_lava_export_contract", CONTRACT_PATH)
SUPERVISOR = _load_module("track_b_lava_export_supervisor", SUPERVISOR_PATH)
RESULT_VALIDATOR = _load_module("track_b_lava_export_result_validator", RESULT_VALIDATOR_PATH)


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError as error:
        raise ExportError(f"path escapes repository: {path}") from error


def _identity(path: Path) -> dict[str, object]:
    try:
        size, digest = CONTRACT.file_identity(path)
    except CONTRACT.ContinuationError as error:
        raise ExportError(str(error)) from error
    return {"path": _relative(path), "bytes": size, "sha256": digest}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ExportError(f"{label} is unreadable: {error}") from error
    if not isinstance(value, dict):
        raise ExportError(f"{label} is not a JSON object")
    return value


def read_tsv(path: Path, expected_fields: list[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
        raise ExportError(f"missing real nonempty TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if expected_fields is not None and fields != expected_fields:
        raise ExportError(
            f"TSV schema drifted for {path}: observed={fields} expected={expected_fields}"
        )
    return fields, rows


def tsv_bytes(fields: list[str], rows: Iterable[dict[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _number(value: str, label: str) -> float:
    try:
        observed = float(value)
    except (TypeError, ValueError) as error:
        raise ExportError(f"invalid numeric {label}: {value}") from error
    if not math.isfinite(observed):
        raise ExportError(f"non-finite numeric {label}: {value}")
    return observed


def _probability(value: str, label: str) -> float:
    observed = _number(value, label)
    if not 0 <= observed <= 1:
        raise ExportError(f"{label} outside [0,1]: {value}")
    return observed


def _boolean(value: str, label: str) -> bool:
    lowered = value.lower()
    if lowered not in {"true", "false"}:
        raise ExportError(f"invalid {label}: {value}")
    return lowered == "true"


def _fmt(value: float | int | None) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return format(float(value), ".12g")


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [0.0] * len(values)
    running = 1.0
    for rank0 in range(len(values) - 1, -1, -1):
        index = order[rank0]
        running = min(running, min(1.0, values[index] * len(values) / (rank0 + 1)))
        adjusted[index] = running
    return adjusted


def _assert_close(observed: str, expected: float, label: str) -> None:
    if not math.isclose(_number(observed, label), expected, rel_tol=1e-8, abs_tol=1e-12):
        raise ExportError(f"incorrect {label}: observed={observed} expected={expected}")


@dataclass(frozen=True)
class SealedTerminal:
    source_fingerprint: str
    continuation_fingerprint: str
    terminal_target: Path
    terminal_receipt: dict[str, Any]
    terminal_attestation: dict[str, Any]
    aggregate_target: Path
    aggregate_receipt: dict[str, Any]
    conditional_checkpoints: tuple[dict[str, object], ...]
    conditional_family_sha256: str
    candidate_path: Path


@dataclass(frozen=True)
class ValidatedFamily:
    policy: dict[str, Any]
    pairs: tuple[dict[str, str], ...]
    status: tuple[dict[str, str], ...]
    univ: tuple[dict[str, str], ...]
    bivar: tuple[dict[str, str], ...]
    conditional: tuple[dict[str, str], ...]
    candidates: tuple[dict[str, str], ...]
    conditional_checkpoint_qc: tuple[dict[str, object], ...]
    terminal_attestation: dict[str, Any]
    profile: dict[str, object]


def _require_terminal_attestation(
    payload: dict[str, Any], source_fingerprint: str, continuation_fingerprint: str,
) -> None:
    required = {
        "state": "TERMINAL_FAILED_QC",
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint,
        "full_family_complete": True,
        "scientific_validation_passed": False,
        "canonical_publication_allowed": False,
        "terminal": True,
    }
    if payload.get("schema_version") != RESULT_VALIDATOR.TERMINAL_QC_SCHEMA:
        raise ExportError("terminal QC attestation schema drifted")
    for key, value in required.items():
        if payload.get(key) != value:
            raise ExportError(f"terminal QC attestation is not a complete failed-QC family ({key})")
    qc = payload.get("qc")
    if not isinstance(qc, dict) or not isinstance(qc.get("reasons"), list) or not qc["reasons"]:
        raise ExportError("terminal QC attestation lacks frozen-gate failure reasons")
    for reason in qc["reasons"]:
        try:
            numerator = int(reason["numerator"])
            denominator = int(reason["denominator"])
            observed = float(reason["observed_fraction"])
            threshold = float(reason["maximum_allowed_fraction"])
        except (KeyError, TypeError, ValueError) as error:
            raise ExportError("terminal QC reason is malformed") from error
        if denominator <= 0 or numerator < 0 or not math.isclose(
            observed, numerator / denominator, rel_tol=1e-12, abs_tol=1e-15,
        ) or reason.get("comparison") != ">" or not observed > threshold:
            raise ExportError("terminal QC reason does not prove a frozen-gate failure")


def _compact_checkpoint(
    receipt: dict[str, Any], target: Path, locus_index: int,
) -> dict[str, object]:
    result = next(
        (item for item in receipt["artifacts"] if item["path"] == "result.rds"), None,
    )
    if result is None:
        raise ExportError(f"conditional checkpoint lacks result.rds: {locus_index}")
    return {
        "locus_index": locus_index,
        "locus": str(receipt["marker"]["locus"]),
        "chromosome": str(receipt["marker"]["chromosome"]),
        "pair": str(receipt["marker"]["pair"]),
        "qc": str(receipt["marker"]["qc"]),
        "receipt_sha256": _identity(target / "receipt.json")["sha256"],
        "result_sha256": str(result["sha256"]),
    }


def deep_validate_terminal(
    source_fingerprint: str,
    continuation_fingerprint: str,
    *,
    progress: Callable[[int, int], None] | None = None,
) -> SealedTerminal:
    """Deeply validate the terminal bundle and every upstream V2 checkpoint."""

    if not FINGERPRINT_RE.fullmatch(source_fingerprint) or not FINGERPRINT_RE.fullmatch(
        continuation_fingerprint
    ):
        raise ExportError("source and continuation fingerprints must be lowercase SHA-256 values")
    if FORBIDDEN_CANONICAL_PROVENANCE.exists() or FORBIDDEN_CANONICAL_PROVENANCE.is_symlink():
        raise ExportError(
            "success-only canonical LAVA provenance exists; failed-QC diagnostic export is forbidden"
        )
    try:
        terminal_bundle = SUPERVISOR.bundle_path(
            "terminal-qc", None, continuation_fingerprint,
        )
    except (SUPERVISOR.ExecutionError, CONTRACT.ContinuationError) as error:
        raise ExportError(str(error)) from error
    if not terminal_bundle.exists() and not terminal_bundle.is_symlink():
        raise ExportError(
            "V2 TERMINAL_FAILED_QC READY bundle does not exist; diagnostic preflight fails closed"
        )
    try:
        CONTRACT.validate_lineage(source_fingerprint, continuation_fingerprint, deep=True)
        terminal_receipt = SUPERVISOR.validate_bundle(
            "terminal-qc", None, source_fingerprint, continuation_fingerprint,
            revalidate_semantics=True, reconcile_active=True,
        )
    except (CONTRACT.ContinuationError, SUPERVISOR.ExecutionError) as error:
        raise ExportError(str(error)) from error
    if terminal_receipt is None:
        raise ExportError("V2 terminal READY receipt is absent")
    terminal_target = SUPERVISOR._safe_attempt_target(
        terminal_bundle, continuation_fingerprint,
    )
    attestation_path = terminal_target / "terminal_qc.json"
    attestation = _read_json(attestation_path, "terminal QC attestation")
    _require_terminal_attestation(attestation, source_fingerprint, continuation_fingerprint)

    artifact_names = {str(item["path"]) for item in terminal_receipt["artifacts"]}
    required_terminal_artifacts = {
        "result.rds", "semantic_validation.txt", "terminal_qc.json",
        "lava_locus_status.tsv", "lava_univariate.tsv", "lava_bivariate.tsv",
        "lava_conditional.tsv",
    }
    if not required_terminal_artifacts.issubset(artifact_names):
        raise ExportError("terminal READY receipt lacks the complete staged diagnostic family")
    if terminal_receipt["marker"].get("qc") != "TERMINAL_FAILED_QC":
        raise ExportError("terminal READY receipt is not marked TERMINAL_FAILED_QC")

    try:
        aggregate_receipt = SUPERVISOR.validate_bundle(
            "aggregate-discovery", None, source_fingerprint, continuation_fingerprint,
            revalidate_semantics=True, reconcile_active=True,
        )
    except SUPERVISOR.ExecutionError as error:
        raise ExportError(str(error)) from error
    if aggregate_receipt is None:
        raise ExportError("terminal family lacks its aggregate-discovery READY checkpoint")
    aggregate_bundle = SUPERVISOR.bundle_path(
        "aggregate-discovery", None, continuation_fingerprint,
    )
    aggregate_target = SUPERVISOR._safe_attempt_target(
        aggregate_bundle, continuation_fingerprint,
    )
    candidate_path = aggregate_target / "conditional_candidates.tsv"
    if "conditional_candidates.tsv" not in {
        str(item["path"]) for item in aggregate_receipt["artifacts"]
    }:
        raise ExportError("aggregate receipt lacks its conditional-candidate artifact")
    read_tsv(candidate_path, CANDIDATE_FIELDS)

    try:
        policy = CONTRACT.legacy_contract().load_policy()
    except Exception as error:
        raise ExportError(f"could not load frozen LAVA policy: {error}") from error
    expected_loci = int(policy["expected_loci"])
    locus_attestation = attestation["qc"]["counts"].get("locus")
    if not isinstance(locus_attestation, dict) or int(locus_attestation.get("family", -1)) != expected_loci:
        raise ExportError("terminal attestation locus-family size differs from frozen policy")

    compact: list[dict[str, object]] = []
    for locus_index in range(1, expected_loci + 1):
        try:
            receipt = SUPERVISOR.validate_bundle(
                "conditional", locus_index, source_fingerprint, continuation_fingerprint,
                revalidate_semantics=True, reconcile_active=True,
            )
        except SUPERVISOR.ExecutionError as error:
            raise ExportError(f"conditional checkpoint {locus_index} failed deep validation: {error}") from error
        if receipt is None:
            raise ExportError(f"missing conditional checkpoint {locus_index}/{expected_loci}")
        bundle = SUPERVISOR.bundle_path("conditional", locus_index, continuation_fingerprint)
        target = SUPERVISOR._safe_attempt_target(bundle, continuation_fingerprint)
        compact.append(_compact_checkpoint(receipt, target, locus_index))
        if progress is not None:
            progress(locus_index, expected_loci)
    if [int(item["locus_index"]) for item in compact] != list(range(1, expected_loci + 1)):
        raise ExportError("conditional checkpoint family is incomplete, duplicated, or unordered")

    upstream = attestation.get("upstream_continuation_family")
    family_hash = _digest([
        {
            "locus_index": item["locus_index"],
            "receipt_sha256": item["receipt_sha256"],
            "result_sha256": item["result_sha256"],
            "qc": item["qc"],
        }
        for item in compact
    ])
    if (
        not isinstance(upstream, dict)
        or int(upstream.get("conditional_checkpoint_count", -1)) != expected_loci
        or upstream.get("conditional_checkpoint_family_sha256") != family_hash
        or upstream.get("complete") is not True
    ):
        raise ExportError("deep conditional checkpoint family differs from terminal attestation")

    return SealedTerminal(
        source_fingerprint=source_fingerprint,
        continuation_fingerprint=continuation_fingerprint,
        terminal_target=terminal_target,
        terminal_receipt=terminal_receipt,
        terminal_attestation=attestation,
        aggregate_target=aggregate_target,
        aggregate_receipt=aggregate_receipt,
        conditional_checkpoints=tuple(compact),
        conditional_family_sha256=family_hash,
        candidate_path=candidate_path,
    )


def _validate_bh(
    rows: list[dict[str, str]], *, tested_status: str, failed_statuses: set[str],
    status_field: str, p_field: str, adjusted_field: str, family_field: str,
    label: str,
) -> int:
    denominator = len(rows)
    if rows and {int(row[family_field]) for row in rows} != {denominator}:
        raise ExportError(f"{label} multiple-testing denominator drifted")
    family_p: list[float] = []
    for row in rows:
        status = row[status_field]
        if status == tested_status:
            family_p.append(_probability(row[p_field], f"{label} P"))
        elif status in failed_statuses:
            family_p.append(1.0)
        else:
            raise ExportError(f"unknown {label} family status: {status}")
    for row, adjusted in zip(rows, bh(family_p), strict=True):
        if row[status_field] == tested_status:
            _assert_close(row[adjusted_field], adjusted, f"{label} BH FDR")
        elif row[adjusted_field] not in NA_VALUES:
            raise ExportError(f"failed {label} row has a reported adjusted P")
    return denominator


def validate_family(sealed: SealedTerminal) -> ValidatedFamily:
    """Read and cross-check all diagnostic tables against the terminal attestation."""

    legacy = RESULT_VALIDATOR.LEGACY
    terminal = sealed.terminal_target
    status_fields, status = read_tsv(terminal / "lava_locus_status.tsv", legacy.STATUS_FIELDS)
    univ_fields, univ = read_tsv(terminal / "lava_univariate.tsv", legacy.UNIV_FIELDS)
    bivar_fields, bivar = read_tsv(terminal / "lava_bivariate.tsv", legacy.BIVAR_FIELDS)
    conditional_fields, conditional = read_tsv(
        terminal / "lava_conditional.tsv", legacy.RAW_CONDITIONAL_FIELDS,
    )
    del status_fields, univ_fields, bivar_fields, conditional_fields
    _, candidates = read_tsv(sealed.candidate_path, CANDIDATE_FIELDS)

    policy = CONTRACT.legacy_contract().load_policy()
    _, pairs = read_tsv(ROOT / "results/track_b/lava_pair_manifest.tsv")
    pair_order = [str(value) for value in policy["pair_order"]]
    if [row["pair_id"] for row in pairs] != pair_order:
        raise ExportError("frozen LAVA pair order drifted")
    expected_loci = int(policy["expected_loci"])
    expected_traits = [str(value) for value in policy["trait_order"]]
    if len(status) != expected_loci or len(sealed.conditional_checkpoints) != expected_loci:
        raise ExportError("terminal locus/checkpoint family is not complete")
    if len(univ) != expected_loci * len(expected_traits):
        raise ExportError("terminal univariate family is not the complete locus x trait grid")
    if len(pairs) * expected_loci != int(policy["planned_bivariate_pair_locus_family_max"]):
        raise ExportError("frozen pair-locus grid size drifted")

    qc_counts = sealed.terminal_attestation["qc"]["counts"]
    failed_loci = sum(row["status"] != "PROCESSED" for row in status)
    untested_univ = sum(row["analysis_status"] != "TESTED" for row in univ)
    failed_bivar = sum(row["analysis_status"] == "BIVARIATE_FAILED" for row in bivar)
    failed_conditional = sum(
        row["analysis_status"] in {"CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2"}
        for row in conditional
    )
    conditional_ineligible = sum(
        row["analysis_status"] == "CONDITIONER_LOCAL_H2_INELIGIBLE" for row in conditional
    )
    conditional_eligible = len(conditional) - conditional_ineligible
    assertions = {
        "locus": (failed_loci, len(status)),
        "univariate": (untested_univ, len(univ)),
        "bivariate": (failed_bivar, len(bivar)),
        "conditional": (failed_conditional, conditional_eligible),
    }
    for family, (failed, denominator) in assertions.items():
        terminal_count = qc_counts.get(family)
        if not isinstance(terminal_count, dict):
            raise ExportError(f"terminal attestation lacks {family} accounting")
        failed_key = "failed_or_untested" if family == "univariate" else "failed"
        family_key = "family" if family in {"locus", "univariate"} else "eligible_family"
        if int(terminal_count.get(failed_key, -1)) != failed or int(
            terminal_count.get(family_key, -1)
        ) != denominator:
            raise ExportError(f"derived {family} counts differ from terminal attestation")
    conditional_count = qc_counts["conditional"]
    if (
        int(conditional_count.get("ineligible", -1)) != conditional_ineligible
        or int(conditional_count.get("reported_family", -1)) != len(conditional)
    ):
        raise ExportError("derived conditional ineligible/reported counts differ from attestation")

    univ_denominator = _validate_bh(
        univ, tested_status="TESTED",
        failed_statuses={"PHENOTYPE_DROPPED", "LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED"},
        status_field="analysis_status", p_field="p", adjusted_field="p_fdr",
        family_field="univariate_test_family_n", label="univariate",
    )
    bivar_denominator = _validate_bh(
        bivar, tested_status="TESTED", failed_statuses={"BIVARIATE_FAILED"},
        status_field="analysis_status", p_field="p", adjusted_field="p_fdr",
        family_field="bivariate_test_family_n", label="bivariate",
    )
    conditional_denominator = _validate_bh(
        [row for row in conditional if row["analysis_status"] != "CONDITIONER_LOCAL_H2_INELIGIBLE"],
        tested_status="TESTED",
        failed_statuses={"CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2"},
        status_field="analysis_status", p_field="p", adjusted_field="p_fdr",
        family_field="conditional_test_family_n", label="conditional",
    )
    if any(
        int(row["conditional_test_family_n"]) != conditional_denominator for row in conditional
    ):
        raise ExportError("conditioner-ineligible rows do not retain the exact eligible denominator")

    process_failures = sum(row["status"] == "PROCESS_FAILED" for row in status)
    phenotype_drops = sum(row["analysis_status"] == "PHENOTYPE_DROPPED" for row in univ)
    profile: dict[str, object] = {
        "official_loci": expected_loci,
        "validated_conditional_checkpoints": len(sealed.conditional_checkpoints),
        "processed_loci": sum(row["status"] == "PROCESSED" for row in status),
        "failed_loci": failed_loci,
        "process_failed_loci": process_failures,
        "univariate_family_n": len(univ),
        "univariate_tested": sum(row["analysis_status"] == "TESTED" for row in univ),
        "univariate_failed_or_untested": untested_univ,
        "phenotype_dropped": phenotype_drops,
        "bivariate_pair_locus_grid_n": len(pairs) * expected_loci,
        "bivariate_bh_family_n": bivar_denominator,
        "bivariate_tested": sum(row["analysis_status"] == "TESTED" for row in bivar),
        "bivariate_failed": failed_bivar,
        "conditional_candidate_loci": len(candidates),
        "conditional_reported_rows": len(conditional),
        "conditional_eligible_bh_family_n": conditional_denominator,
        "conditional_ineligible_rows": conditional_ineligible,
        "terminal_qc_reason_codes": [
            str(reason["code"]) for reason in sealed.terminal_attestation["qc"]["reasons"]
        ],
    }
    if univ_denominator != len(univ):
        raise ExportError("univariate denominator does not span the complete frozen family")

    return ValidatedFamily(
        policy=policy,
        pairs=tuple(pairs),
        status=tuple(status),
        univ=tuple(univ),
        bivar=tuple(bivar),
        conditional=tuple(conditional),
        candidates=tuple(candidates),
        conditional_checkpoint_qc=sealed.conditional_checkpoints,
        terminal_attestation=sealed.terminal_attestation,
        profile=profile,
    )


def build_univariate_rows(family: ValidatedFamily) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in family.univ:
        tested = source["analysis_status"] == "TESTED"
        rows.append(dict(source) | {
            "correction_only_family_p": source["p"] if tested else "1",
            "correction_only_role": (
                "OBSERVED_TEST_P" if tested
                else "IMPUTED_ONE_FOR_COMPLETE_FAMILY_BH_ONLY_NOT_AN_OBSERVED_NULL"
            ),
            "terminal_scientific_state": SCIENTIFIC_STATE,
        })
    return rows


def build_bivariate_grid(family: ValidatedFamily) -> list[dict[str, str]]:
    status_index = {row["LOC"]: row for row in family.status}
    univ_index = {(row["LOC"], row["phen"]): row for row in family.univ}
    bivar_index = {(row["LOC"], row["pair_id"]): row for row in family.bivar}
    loci = [row["LOC"] for row in family.status]
    threshold = float(family.policy["univariate_p_threshold"])
    eligible_denominator = len(family.bivar)
    grid_denominator = len(loci) * len(family.pairs)
    rows: list[dict[str, str]] = []
    for pair in family.pairs:
        for locus in loci:
            status = status_index[locus]
            first = univ_index[(locus, pair["trait1"])]
            second = univ_index[(locus, pair["trait2"])]
            eligible = all(
                row["analysis_status"] == "TESTED"
                and row["p"] not in NA_VALUES
                and _probability(row["p"], "local-h2 eligibility P") <= threshold
                for row in (first, second)
            )
            raw = bivar_index.get((locus, pair["pair_id"]))
            if (raw is not None) != eligible:
                raise ExportError(
                    f"bivariate row/local-h2 eligibility mismatch: {locus}/{pair['pair_id']}"
                )
            if raw is not None:
                analysis_status = raw["analysis_status"]
                reason = "BOTH_TRAITS_PASSED_FROZEN_LOCAL_H2_GATE"
                correction_p = raw["p"] if analysis_status == "TESTED" else "1"
                correction_role = (
                    "OBSERVED_TEST_P" if analysis_status == "TESTED"
                    else "IMPUTED_ONE_FOR_ELIGIBLE_FAILED_TEST_BH_ONLY_NOT_AN_OBSERVED_NULL"
                )
                stability = "PASS" if analysis_status == "TESTED" else "BIVARIATE_FAILURE_RETAINED"
            elif status["status"] == "PROCESS_FAILED":
                analysis_status = "LOCUS_PROCESS_FAILED"
                reason = "WHOLE_LOCUS_PROCESS_FAILED_BEFORE_LOCAL_H2_ELIGIBILITY"
                correction_p = correction_role = "NA"
                stability = "NOT_EVALUATED"
            elif status["status"] == "UNIVARIATE_FAILED":
                analysis_status = "LOCUS_UNIVARIATE_FAILED"
                reason = "WHOLE_LOCUS_UNIVARIATE_STAGE_FAILED"
                correction_p = correction_role = "NA"
                stability = "NOT_EVALUATED"
            else:
                analysis_status = "NOT_ELIGIBLE_LOCAL_H2"
                reason = (
                    f"TARGET_LOCAL_H2_STATES:{first['analysis_status']};{second['analysis_status']}"
                )
                correction_p = correction_role = "NA"
                stability = "NOT_EVALUATED"
            raw_value = lambda field: raw[field] if raw is not None else "NA"
            error_parts = [] if raw is not None else [
                value for value in (first.get("error", ""), second.get("error", ""))
                if value not in NA_VALUES
            ]
            error = raw_value("error") if raw is not None else (
                " | ".join(dict.fromkeys(error_parts)) or reason
            )
            rows.append({
                "LOC": locus, "CHR": status["CHR"], "START": status["START"],
                "STOP": status["STOP"], "pair_id": pair["pair_id"],
                "trait1": pair["trait1"], "trait2": pair["trait2"],
                "locus_status": status["status"],
                "trait1_univariate_status": first["analysis_status"],
                "trait2_univariate_status": second["analysis_status"],
                "trait1_local_h2_obs": first["h2.obs"],
                "trait2_local_h2_obs": second["h2.obs"],
                "both_traits_local_h2_eligible": str(eligible).upper(),
                "analysis_status": analysis_status, "eligibility_reason": reason,
                "local_covariance": raw_value("local_covariance"), "rho": raw_value("rho"),
                "rho.lower": raw_value("rho.lower"), "rho.upper": raw_value("rho.upper"),
                "r2": raw_value("r2"), "r2.lower": raw_value("r2.lower"),
                "r2.upper": raw_value("r2.upper"), "p": raw_value("p"),
                "correction_only_family_p": correction_p,
                "correction_only_role": correction_role,
                "bivariate_test_family_n": str(eligible_denominator),
                "predeclared_pair_locus_family_n": str(grid_denominator),
                "p_fdr": raw_value("p_fdr"),
                "fdr_significant": raw_value("fdr_significant"),
                "n_snps": status["n_snps"], "K": status["K"], "error": error,
                "numerical_stability": stability,
                "terminal_scientific_state": SCIENTIFIC_STATE,
            })
    if len(rows) != grid_denominator:
        raise ExportError("full bivariate pair-locus grid was not preserved")
    return rows


def build_bivariate_fdr_rows(family: ValidatedFamily) -> list[dict[str, str]]:
    rule = str(family.policy["bivariate_multiple_testing"])
    rows: list[dict[str, str]] = []
    for source in family.bivar:
        tested = source["analysis_status"] == "TESTED"
        rows.append(dict(source) | {
            "correction_only_family_p": source["p"] if tested else "1",
            "correction_only_role": (
                "OBSERVED_TEST_P" if tested
                else "IMPUTED_ONE_FOR_ELIGIBLE_FAILED_TEST_BH_ONLY_NOT_AN_OBSERVED_NULL"
            ),
            "multiple_testing_rule": rule,
            "terminal_scientific_state": SCIENTIFIC_STATE,
        })
    return rows


def build_failure_rows(
    family: ValidatedFamily, grid: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    def append(
        record_class: str, stage: str, row: dict[str, str], *, pair_id: str = "NA",
        phenotype: str = "NA", model: str = "NA", error: str = "NA",
        correction: str = "NOT_IN_MULTIPLE_TESTING_FAMILY",
    ) -> None:
        rows.append({
            "record_class": record_class, "stage": stage, "pair_id": pair_id,
            "LOC": row.get("LOC", "NA"), "CHR": row.get("CHR", "NA"),
            "START": row.get("START", "NA"), "STOP": row.get("STOP", "NA"),
            "phenotype": phenotype, "conditional_model_id": model,
            "analysis_status": row["analysis_status"], "error": error,
            "correction_treatment": correction,
            "terminal_scientific_state": SCIENTIFIC_STATE,
        })

    for row in family.status:
        if row["status"] != "PROCESSED":
            append(
                "FAILURE_RETAINED", "LOCUS", dict(row, analysis_status=row["status"]),
                error="Locus terminal execution status retained; see phenotype rows for diagnostics",
            )
    for row in family.univ:
        if row["analysis_status"] != "TESTED":
            append(
                "FAILURE_RETAINED", "UNIVARIATE", row, phenotype=row["phen"],
                error=row["error"],
                correction="P_EQUALS_1_FOR_COMPLETE_UNIVARIATE_BH_ONLY_NOT_AN_OBSERVED_NULL",
            )
    for row in grid:
        if row["analysis_status"] != "TESTED":
            record_class = (
                "FAILURE_RETAINED" if row["analysis_status"] == "BIVARIATE_FAILED"
                else "INELIGIBLE_NOT_NULL"
            )
            correction = (
                "P_EQUALS_1_FOR_ELIGIBLE_BIVARIATE_BH_ONLY_NOT_AN_OBSERVED_NULL"
                if row["analysis_status"] == "BIVARIATE_FAILED"
                else "EXCLUDED_FROM_ELIGIBLE_BH_FAMILY_NOT_AN_OBSERVED_NULL"
            )
            append(
                record_class, "BIVARIATE_PAIR_LOCUS_GRID", row,
                pair_id=row["pair_id"], error=row["error"], correction=correction,
            )
    for row in family.conditional:
        if row["analysis_status"] != "TESTED":
            eligible_failure = row["analysis_status"] in {
                "CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2",
            }
            append(
                "FAILURE_RETAINED" if eligible_failure else "INELIGIBLE_NOT_NULL",
                "CONDITIONAL_MODEL", row, pair_id=row["pair_id"],
                model=row["conditional_model_id"], error=row["error"],
                correction=(
                    "P_EQUALS_1_FOR_ELIGIBLE_CONDITIONAL_BH_ONLY_NOT_AN_OBSERVED_NULL"
                    if eligible_failure
                    else "EXCLUDED_FROM_ELIGIBLE_CONDITIONAL_BH_NOT_AN_OBSERVED_NULL"
                ),
            )
    status_by_locus = {row["LOC"]: row for row in family.status}
    for checkpoint in family.conditional_checkpoint_qc:
        if checkpoint["qc"] != "CONDITIONAL_LOCUS_PROCESSED":
            locus_status = status_by_locus[str(checkpoint["locus"])]
            append(
                "INAPPLICABLE_NOT_NULL", "CONDITIONAL_CHECKPOINT",
                dict(locus_status, analysis_status=str(checkpoint["qc"])),
                pair_id=str(checkpoint["pair"]),
                error="Immutable per-locus conditional checkpoint state",
            )
    return rows


def build_accounting_rows(
    family: ValidatedFamily, grid: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    reasons = {
        str(reason["metric"]): reason for reason in family.terminal_attestation["qc"]["reasons"]
    }

    def add_counts(
        name: str, values: Iterable[str], total: int, denominator: int | None,
        rule: str, metric: str | None = None, notes: str = "",
    ) -> None:
        counts = Counter(values)
        if not counts:
            counts["NONE"] = 0
        reason = reasons.get(metric or "")
        for status, count in sorted(counts.items()):
            rows.append({
                "family": name, "status": status, "count": str(count),
                "family_total": str(total),
                "multiple_testing_denominator": "NA" if denominator is None else str(denominator),
                "multiple_testing_rule": rule, "threshold_name": metric or "NA",
                "threshold": (
                    "NA" if reason is None else str(reason["maximum_allowed_fraction_decimal"])
                ),
                "observed_fraction": (
                    "NA" if reason is None else str(reason["observed_fraction_decimal"])
                ),
                "gate_state": "FAILED" if reason is not None else "NOT_FAILED_OR_NOT_APPLICABLE",
                "notes": notes, "terminal_scientific_state": SCIENTIFIC_STATE,
            })

    add_counts(
        "LOCUS_EXECUTION", (row["status"] for row in family.status), len(family.status), None,
        "NO_MULTIPLE_TESTING", "locus_failure_fraction", "All official loci retained",
    )
    add_counts(
        "LOCAL_UNIVARIATE", (row["analysis_status"] for row in family.univ), len(family.univ),
        len(family.univ), "BH_OVER_COMPLETE_LOCUS_X_TRAIT_FAMILY_FAILED_OR_UNTESTED_AS_ONE",
        "univariate_untested_fraction", "Correction-only P=1 is never an observed null",
    )
    add_counts(
        "LOCAL_BIVARIATE_PAIR_LOCUS_GRID", (row["analysis_status"] for row in grid), len(grid),
        len(family.bivar), str(family.policy["bivariate_multiple_testing"]), None,
        "Grid includes every frozen pair x official locus; only locally h2-eligible rows enter BH",
    )
    add_counts(
        "LOCAL_BIVARIATE_ELIGIBLE_BH", (row["analysis_status"] for row in family.bivar),
        len(family.bivar), len(family.bivar), str(family.policy["bivariate_multiple_testing"]),
        "bivariate_failure_fraction", "Failed eligible tests contribute correction-only P=1",
    )
    add_counts(
        "CONDITIONAL_CANDIDATE_LOCI", ("ELIGIBLE_CANDIDATE" for _ in family.candidates),
        len(family.candidates), None, "NO_CORRECTION_AT_CANDIDATE_MANIFEST", None,
        "Candidate requires bivariate FDR support and same-locus conditioner h2 eligibility",
    )
    add_counts(
        "CONDITIONAL_REPORTED_MODELS", (row["analysis_status"] for row in family.conditional),
        len(family.conditional), int(family.profile["conditional_eligible_bh_family_n"]),
        str(family.policy["conditional_multiple_testing"]), "conditional_failure_fraction",
        "Conditioner-ineligible rows are retained but excluded from eligible conditional BH",
    )
    add_counts(
        "CONDITIONAL_CHECKPOINTS", (str(row["qc"]) for row in family.conditional_checkpoint_qc),
        len(family.conditional_checkpoint_qc), None, "NO_MULTIPLE_TESTING", None,
        "One immutable checkpoint per official locus",
    )
    for reason in family.terminal_attestation["qc"]["reasons"]:
        rows.append({
            "family": "TERMINAL_QC_GATE", "status": str(reason["code"]),
            "count": str(reason["numerator"]), "family_total": str(reason["denominator"]),
            "multiple_testing_denominator": "NA", "multiple_testing_rule": "NO_MULTIPLE_TESTING",
            "threshold_name": str(reason["metric"]),
            "threshold": str(reason["maximum_allowed_fraction_decimal"]),
            "observed_fraction": str(reason["observed_fraction_decimal"]),
            "gate_state": "FAILED", "notes": "Frozen QC threshold exceeded",
            "terminal_scientific_state": SCIENTIFIC_STATE,
        })
    return rows


def _pair_architecture(
    pair_id: str, tested: list[dict[str, str]], significant: list[dict[str, str]],
) -> str:
    del pair_id
    if not tested:
        return "INSUFFICIENT_VALID_LOCAL_RG_TO_CLASSIFY"
    if not significant:
        return "NO_FDR_SUPPORTED_LOCAL_CONCENTRATION_DETECTED"
    signs = {1 if float(row["rho"]) > 0 else -1 if float(row["rho"]) < 0 else 0 for row in significant}
    if 1 in signs and -1 in signs:
        return "OPPOSING_FDR_SUPPORTED_LOCAL_EFFECTS"
    if len(significant) == 1:
        return "ONE_DETECTABLE_FDR_SUPPORTED_LOCUS"
    return "SEVERAL_SAME_DIRECTION_FDR_SUPPORTED_LOCI"


def build_summary_markdown(
    family: ValidatedFamily, sealed: SealedTerminal, grid: list[dict[str, str]],
) -> str:
    profile = family.profile
    lines = [
        "# Track B LAVA discovery diagnostic summary",
        "",
        "> **Terminal state: `TERMINAL_FAILED_QC`.** This is a complete-family diagnostic export, "
        "not a canonical scientific PASS. The success-only `lava_results.provenance.json` was not "
        "created, and these tables must not be cited as a QC-passing LAVA result set.",
        "",
        "## Sealed identity",
        "",
        f"- Source discovery fingerprint: `{sealed.source_fingerprint}`",
        f"- Continuation execution fingerprint: `{sealed.continuation_fingerprint}`",
        f"- Official loci accounted for: {profile['validated_conditional_checkpoints']}/{profile['official_loci']}",
        f"- Full bivariate pair-locus grid: {profile['bivariate_pair_locus_grid_n']}",
        "",
        "## Completion and QC accounting",
        "",
        "| Metric | Derived value |",
        "|---|---:|",
        f"| Processed loci | {profile['processed_loci']} |",
        f"| Failed loci (all locus failure states) | {profile['failed_loci']} |",
        f"| `PROCESS_FAILED` loci | {profile['process_failed_loci']} |",
        f"| Univariate family | {profile['univariate_family_n']} |",
        f"| Tested univariate rows | {profile['univariate_tested']} |",
        f"| Failed/untested univariate rows | {profile['univariate_failed_or_untested']} |",
        f"| `PHENOTYPE_DROPPED` rows | {profile['phenotype_dropped']} |",
        f"| Locally h2-eligible bivariate BH denominator | {profile['bivariate_bh_family_n']} |",
        f"| Valid bivariate tests | {profile['bivariate_tested']} |",
        f"| Failed eligible bivariate tests | {profile['bivariate_failed']} |",
        f"| Conditional candidate loci | {profile['conditional_candidate_loci']} |",
        f"| Reported conditional model rows | {profile['conditional_reported_rows']} |",
        f"| Eligible conditional BH denominator | {profile['conditional_eligible_bh_family_n']} |",
        f"| Conditioner-local-h2-ineligible rows | {profile['conditional_ineligible_rows']} |",
        "",
        "Frozen multiple-testing rules were retained exactly. Untested univariate rows and failed "
        "eligible tests contribute a value of 1 only inside the predeclared correction family; that "
        "correction-only value is not an observed null P-value. Pair-locus rows that never passed "
        "the local-h2 gate are excluded from the eligible bivariate BH family and remain explicitly "
        "labeled `NOT_ELIGIBLE_LOCAL_H2`.",
        "",
        "### Failed frozen QC gates",
        "",
        "| Code | Numerator | Denominator | Observed | Maximum allowed |",
        "|---|---:|---:|---:|---:|",
    ]
    for reason in family.terminal_attestation["qc"]["reasons"]:
        lines.append(
            f"| `{reason['code']}` | {reason['numerator']} | {reason['denominator']} | "
            f"{float(reason['observed_fraction']):.6%} | "
            f"{float(reason['maximum_allowed_fraction']):.6%} |"
        )

    lines.extend([
        "",
        "## Local architecture by frozen pair",
        "",
        "| Pair | Genome-wide rg | Eligible local rg rows | Valid | FDR-supported | Architecture label |",
        "|---|---:|---:|---:|---:|---|",
    ])
    raw_by_pair = {pair["pair_id"]: [] for pair in family.pairs}
    for row in family.bivar:
        raw_by_pair[row["pair_id"]].append(row)
    for pair in family.pairs:
        eligible = raw_by_pair[pair["pair_id"]]
        tested = [row for row in eligible if row["analysis_status"] == "TESTED"]
        significant = [
            row for row in tested if _boolean(row["fdr_significant"], "bivariate significance")
        ]
        lines.append(
            f"| {pair['pair_id']} ({pair['trait1']} ↔ {pair['trait2']}) | "
            f"{pair['discovery_rg']} | {len(eligible)} | {len(tested)} | {len(significant)} | "
            f"`{_pair_architecture(pair['pair_id'], tested, significant)}` |"
        )

    lines.extend([
        "",
        "Because the frozen QC gates failed, an architecture label describes only the detectable "
        "sealed rows; it is not evidence that untested loci are null or that sharing is absent.",
        "",
        "## Strongest diagnostic local signals for Pair A and Pair B",
        "",
    ])
    grid_index = {(row["LOC"], row["pair_id"]): row for row in grid}
    for pair_id in ("A", "B"):
        pair = next(row for row in family.pairs if row["pair_id"] == pair_id)
        tested = [row for row in raw_by_pair[pair_id] if row["analysis_status"] == "TESTED"]
        tested.sort(key=lambda row: (
            not _boolean(row["fdr_significant"], "bivariate significance"),
            _probability(row["p_fdr"], "bivariate FDR"),
            -abs(_number(row["rho"], "local rg")),
        ))
        lines.extend([f"### Pair {pair_id}: {pair['trait1']} ↔ {pair['trait2']}", ""])
        if not tested:
            lines.extend([
                "No valid local-rg row was available. This is an underpowered/ineligible state, not a null result.",
                "",
            ])
            continue
        lines.extend([
            "| Rank | Locus | Chr | Local h2 trait 1 | Local h2 trait 2 | Local rg | P | BH FDR | Stability | Opposes genome-wide direction? |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
        ])
        global_rg = _number(pair["discovery_rg"], "genome-wide rg")
        for rank, row in enumerate(tested[:10], start=1):
            full = grid_index[(row["LOC"], pair_id)]
            opposite = _number(row["rho"], "local rg") * global_rg < 0
            lines.append(
                f"| {rank} | {row['LOC']} | {row['CHR']} | {full['trait1_local_h2_obs']} | "
                f"{full['trait2_local_h2_obs']} | {row['rho']} | {row['p']} | {row['p_fdr']} | "
                f"{full['numerical_stability']} | {'YES' if opposite else 'NO'} |"
            )
        lines.append("")

    lines.extend([
        "A large or FDR-supported local genetic correlation does **not** establish causality. "
        "Opposition to the genome-wide direction is flagged explicitly above.",
        "",
        "## Conditional pass",
        "",
        f"The sealed aggregate produced {profile['conditional_candidate_loci']} fully eligible "
        "conditional candidate loci. All frozen per-locus conditional checkpoints are nevertheless "
        "present, and conditioner-ineligible/failed states remain in the failure ledger. No skipped "
        "locus is called null.",
        "",
        "## Permitted interpretation",
        "",
        "These artifacts support transparent accounting, method debugging, and planning of independent "
        "analyses. They do not override the frozen QC failure, validate a LAVA discovery claim, or "
        "authorize canonical LAVA-dependent triangulation.",
        "",
    ])
    return "\n".join(lines)


def render_outputs(
    family: ValidatedFamily, sealed: SealedTerminal,
) -> tuple[dict[Path, bytes], dict[str, object]]:
    univ = build_univariate_rows(family)
    grid = build_bivariate_grid(family)
    fdr = build_bivariate_fdr_rows(family)
    failures = build_failure_rows(family, grid)
    accounting = build_accounting_rows(family, grid)
    outputs = {
        OUTPUT_DIR / "local_univariate_results.tsv": tsv_bytes(
            RESULT_VALIDATOR.LEGACY.UNIV_FIELDS + UNIV_EXTRA_FIELDS, univ,
        ),
        OUTPUT_DIR / "local_bivariate_results.tsv": tsv_bytes(BIVARIATE_GRID_FIELDS, grid),
        OUTPUT_DIR / "local_bivariate_fdr.tsv": tsv_bytes(
            RESULT_VALIDATOR.LEGACY.BIVAR_FIELDS + BIVARIATE_FDR_EXTRA_FIELDS, fdr,
        ),
        OUTPUT_DIR / "local_failures.tsv": tsv_bytes(FAILURE_FIELDS, failures),
        OUTPUT_DIR / "local_family_accounting.tsv": tsv_bytes(ACCOUNTING_FIELDS, accounting),
        SUMMARY_PATH: build_summary_markdown(family, sealed, grid).encode("utf-8"),
    }
    derived = dict(family.profile) | {
        "bivariate_grid_status_counts": dict(Counter(row["analysis_status"] for row in grid)),
        "failure_ledger_rows": len(failures),
        "accounting_rows": len(accounting),
    }
    return outputs, derived


def versioned_lock_path(continuation_fingerprint: str) -> Path:
    if not FINGERPRINT_RE.fullmatch(continuation_fingerprint):
        raise ExportError("invalid continuation fingerprint for versioned export lock")
    return OUTPUT_DIR / f"diagnostic_export.{continuation_fingerprint}.provenance.lock.json"


def _output_record(path: Path, content: bytes) -> dict[str, object]:
    return {
        "path": _relative(path), "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def build_lock_payload(
    sealed: SealedTerminal, outputs: dict[Path, bytes], derived: dict[str, object],
) -> dict[str, object]:
    terminal_receipt_path = sealed.terminal_target / "receipt.json"
    aggregate_receipt_path = sealed.aggregate_target / "receipt.json"
    return {
        "schema_version": EXPORT_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "export_kind": "COMPLETE_FAMILY_TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY",
        "scientific_state": "TERMINAL_FAILED_QC",
        "scientific_validation_passed": False,
        "canonical_publication_allowed": False,
        "source_discovery_fingerprint": sealed.source_fingerprint,
        "continuation_execution_fingerprint": sealed.continuation_fingerprint,
        "generator": _identity(Path(__file__).resolve()),
        "terminal_bundle": {
            "path": _relative(SUPERVISOR.bundle_path(
                "terminal-qc", None, sealed.continuation_fingerprint,
            )),
            "receipt": _identity(terminal_receipt_path),
            "attestation": _identity(sealed.terminal_target / "terminal_qc.json"),
            "result": _identity(sealed.terminal_target / "result.rds"),
        },
        "aggregate_bundle": {
            "path": _relative(SUPERVISOR.bundle_path(
                "aggregate-discovery", None, sealed.continuation_fingerprint,
            )),
            "receipt": _identity(aggregate_receipt_path),
            "result": _identity(sealed.aggregate_target / "result.rds"),
            "conditional_candidates": _identity(sealed.candidate_path),
        },
        "conditional_checkpoint_family": {
            "count": len(sealed.conditional_checkpoints),
            "family_sha256": sealed.conditional_family_sha256,
            "deep_semantic_validation": True,
            "active_input_reconciliation": True,
        },
        "staged_low_level_inputs": [
            _identity(sealed.terminal_target / name) for name in (
                "lava_locus_status.tsv", "lava_univariate.tsv", "lava_bivariate.tsv",
                "lava_conditional.tsv",
            )
        ],
        "derived_accounting": derived,
        "multiple_testing": {
            "univariate": "COMPLETE_FROZEN_LOCUS_X_TRAIT_FAMILY_FAILED_OR_UNTESTED_AS_CORRECTION_ONLY_ONE",
            "bivariate": "ALL_LOCALLY_H2_ELIGIBLE_ROWS_ACROSS_THREE_FROZEN_PAIRS_FAILED_ELIGIBLE_AS_CORRECTION_ONLY_ONE",
            "conditional": "ALL_ELIGIBLE_MODEL_LOCUS_ROWS_CONDITIONER_INELIGIBLE_EXCLUDED_FAILED_ELIGIBLE_AS_CORRECTION_ONLY_ONE",
            "unobserved_rows_called_null": False,
        },
        "outputs": [_output_record(path, content) for path, content in outputs.items()],
        "forbidden_canonical_provenance": {
            "path": _relative(FORBIDDEN_CANONICAL_PROVENANCE),
            "required_absent": True,
            "created_by_exporter": False,
        },
    }


def _publish_or_verify(path: Path, content: bytes) -> None:
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
            raise ExportError(f"refusing to replace differing diagnostic artifact: {path}")
        return
    try:
        CONTRACT.publish_no_replace(path, content)
    except CONTRACT.ContinuationError as error:
        if path.is_file() and not path.is_symlink() and path.read_bytes() == content:
            return
        raise ExportError(str(error)) from error


def publish_outputs(
    sealed: SealedTerminal, outputs: dict[Path, bytes], derived: dict[str, object],
) -> Path:
    if FORBIDDEN_CANONICAL_PROVENANCE.exists() or FORBIDDEN_CANONICAL_PROVENANCE.is_symlink():
        raise ExportError("canonical success provenance appeared before diagnostic publication")
    for path, content in outputs.items():
        if path.exists() or path.is_symlink():
            if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
                raise ExportError(f"pre-existing diagnostic output differs; overwrite forbidden: {path}")
    for path, content in outputs.items():
        _publish_or_verify(path, content)
    payload = build_lock_payload(sealed, outputs, derived)
    lock = versioned_lock_path(sealed.continuation_fingerprint)
    encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    _publish_or_verify(lock, encoded)
    if FORBIDDEN_CANONICAL_PROVENANCE.exists() or FORBIDDEN_CANONICAL_PROVENANCE.is_symlink():
        raise ExportError("canonical success provenance appeared during diagnostic publication")
    return lock


def verify_outputs(
    sealed: SealedTerminal, outputs: dict[Path, bytes], derived: dict[str, object],
) -> Path:
    for path, expected in outputs.items():
        if not path.is_file() or path.is_symlink() or path.read_bytes() != expected:
            raise ExportError(f"diagnostic output is absent or differs from sealed family: {path}")
    payload = build_lock_payload(sealed, outputs, derived)
    lock = versioned_lock_path(sealed.continuation_fingerprint)
    expected_lock = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    if not lock.is_file() or lock.is_symlink() or lock.read_bytes() != expected_lock:
        raise ExportError("versioned diagnostic provenance lock is absent or differs")
    return lock


def _progress(index: int, total: int) -> None:
    if index == total or index % 250 == 0:
        print(
            f"TRACK_B_LAVA_DIAGNOSTIC_DEEP_VALIDATION conditional={index}/{total}",
            file=sys.stderr, flush=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--export", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument(
        "--source-fingerprint", default=CONTRACT.PINNED_SOURCE_FINGERPRINT,
    )
    parser.add_argument("--continuation-fingerprint", required=True)
    parser.add_argument("--quiet-progress", action="store_true")
    args = parser.parse_args()

    sealed = deep_validate_terminal(
        args.source_fingerprint, args.continuation_fingerprint,
        progress=None if args.quiet_progress else _progress,
    )
    family = validate_family(sealed)
    outputs, derived = render_outputs(family, sealed)
    if args.preflight:
        print(
            "TRACK_B_LAVA_DIAGNOSTIC_PREFLIGHT_PASS "
            f"source={sealed.source_fingerprint} continuation={sealed.continuation_fingerprint} "
            f"state=TERMINAL_FAILED_QC loci={derived['official_loci']} "
            f"conditional_checkpoints={derived['validated_conditional_checkpoints']} "
            f"bivariate_bh_family={derived['bivariate_bh_family_n']} "
            f"conditional_candidates={derived['conditional_candidate_loci']}"
        )
    elif args.export:
        lock = publish_outputs(sealed, outputs, derived)
        print(
            "TRACK_B_LAVA_DIAGNOSTIC_EXPORT_PUBLISHED "
            f"state=TERMINAL_FAILED_QC canonical_pass=FALSE lock={_relative(lock)}"
        )
    else:
        lock = verify_outputs(sealed, outputs, derived)
        print(
            "TRACK_B_LAVA_DIAGNOSTIC_EXPORT_VERIFIED "
            f"state=TERMINAL_FAILED_QC canonical_pass=FALSE lock={_relative(lock)}"
        )


if __name__ == "__main__":
    try:
        main()
    except (ExportError, CONTRACT.ContinuationError, SUPERVISOR.ExecutionError) as error:
        raise SystemExit(f"ERROR: {error}") from error
