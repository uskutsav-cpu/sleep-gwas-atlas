#!/usr/bin/env python3
"""Validate and publish the additive Track B LAVA V2 continuation.

The validator deliberately distinguishes a complete scientific family that
fails frozen QC gates from an incomplete or malformed family.  The former is
sealed as an immutable ``TERMINAL_FAILED_QC`` attestation and is never eligible
for canonical publication.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "scripts/134_track_b_lava_continuation_contract.py"
SUPERVISOR_PATH = ROOT / "scripts/136_run_track_b_lava_continuation.py"
LEGACY_VALIDATOR_PATH = ROOT / "scripts/121_validate_track_b_lava.py"
TERMINAL_QC_SCHEMA = "track-b-lava-terminal-qc.1"
RESULT_LOCK_SCHEMA = "track-b-lava-continuation-results.1"
TERMINAL_EXIT_STATUS = 78


class ValidationError(RuntimeError):
    """A malformed, incomplete, stale, or otherwise invalid result family."""


class TerminalQCOutcome(RuntimeError):
    """A complete result family that fails one or more frozen QC gates."""


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ValidationError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


CONTRACT = _load_module("track_b_lava_continuation_contract_for_results", CONTRACT_PATH)
LEGACY = _load_module("track_b_lava_legacy_result_validator", LEGACY_VALIDATOR_PATH)


def supervisor():
    return _load_module("track_b_lava_continuation_supervisor_for_results", SUPERVISOR_PATH)


def _legacy_call(function: Callable[..., Any], *args: object, **kwargs: object) -> Any:
    try:
        return function(*args, **kwargs)
    except SystemExit as error:
        raise ValidationError(str(error)) from error


def _identity(path: Path) -> dict[str, object]:
    try:
        size, digest = CONTRACT.file_identity(path)
        relative = str(path.relative_to(ROOT))
    except (CONTRACT.ContinuationError, ValueError) as error:
        raise ValidationError(f"could not identify continuation artifact {path}: {error}") from error
    return {"path": relative, "bytes": size, "sha256": digest}


def _digest(value: object) -> str:
    return CONTRACT._digest(value)


def _quick_contract(source_fingerprint: str, continuation_fingerprint: str) -> None:
    try:
        CONTRACT._require_source_fingerprint(source_fingerprint)
        CONTRACT.validate_lineage(source_fingerprint, continuation_fingerprint, deep=False)
    except CONTRACT.ContinuationError as error:
        raise ValidationError(str(error)) from error


def _attempt_directory(directory: Path, continuation_fingerprint: str) -> Path:
    if directory.is_symlink():
        raise ValidationError("staging directory may not be a symlink")
    try:
        resolved = directory.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ValidationError(f"staging directory is unreadable: {error}") from error
    expected_parent = (CONTRACT.continuation_run_root(continuation_fingerprint) / ".attempts").resolve()
    if (
        resolved.parent != expected_parent
        or not resolved.name.startswith("finalize_unit_all.")
    ):
        raise ValidationError("staging directory is not a supervised V2 finalize attempt")
    return resolved


def _staged_paths(directory: Path) -> list[Path]:
    return [directory / path.name for path in LEGACY.contract.RESULTS]


def _fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def _reason(
    code: str, metric: str, numerator: int, denominator: int, threshold: float,
) -> dict[str, object]:
    observed = _fraction(numerator, denominator)
    return {
        "code": code,
        "metric": metric,
        "numerator": numerator,
        "denominator": denominator,
        "observed_fraction": observed,
        "observed_fraction_decimal": format(observed, ".17g"),
        "comparison": ">",
        "maximum_allowed_fraction": threshold,
        "maximum_allowed_fraction_decimal": format(threshold, ".17g"),
    }


def _diagnostic_validate(paths: list[Path]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Run every legacy scientific check while collecting all QC gate failures.

    The four failure-fraction gates are deferred in memory only so the legacy
    validator can reach its schema, topology, full-family BH, and numerical
    checks.  Acceptance is never relaxed: a passing family is subsequently run
    through the unmodified validator with the frozen thresholds.
    """

    before = {path: _identity(path) for path in paths[:4]}
    original_load_policy = LEGACY.contract.load_policy

    def diagnostic_policy() -> dict[str, Any]:
        policy = copy.deepcopy(original_load_policy())
        for key in (
            "maximum_locus_failure_fraction",
            "maximum_univariate_untested_fraction",
            "maximum_bivariate_failure_fraction",
            "maximum_conditional_failure_fraction",
        ):
            policy[key] = 1.0
        return policy

    LEGACY.contract.load_policy = diagnostic_policy
    try:
        context = _legacy_call(LEGACY.validate_low_level, paths[:4])
    finally:
        LEGACY.contract.load_policy = original_load_policy
    after = {path: _identity(path) for path in paths[:4]}
    if before != after:
        raise ValidationError("staged low-level family changed during diagnostic validation")

    policy = original_load_policy()
    status = context["status"]
    univ = context["univ"]
    bivar = context["bivar"]
    conditional = context["conditional"]
    counts = {
        "locus": {
            "failed": sum(row["status"] != "PROCESSED" for row in status),
            "family": len(status),
        },
        "univariate": {
            "failed_or_untested": sum(row["analysis_status"] != "TESTED" for row in univ),
            "family": len(univ),
        },
        "bivariate": {
            "failed": sum(row["analysis_status"] == "BIVARIATE_FAILED" for row in bivar),
            "eligible_family": len(bivar),
        },
        "conditional": {
            "failed": sum(
                row["analysis_status"] in {"CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2"}
                for row in conditional
            ),
            "eligible_family": sum(
                row["analysis_status"] != "CONDITIONER_LOCAL_H2_INELIGIBLE"
                for row in conditional
            ),
            "ineligible": sum(
                row["analysis_status"] == "CONDITIONER_LOCAL_H2_INELIGIBLE"
                for row in conditional
            ),
            "reported_family": len(conditional),
        },
    }
    gates = [
        (
            "LOCUS_FAILURE_FRACTION_EXCEEDED", "locus_failure_fraction",
            counts["locus"]["failed"], counts["locus"]["family"],
            float(policy["maximum_locus_failure_fraction"]),
        ),
        (
            "UNIVARIATE_UNTESTED_FRACTION_EXCEEDED", "univariate_untested_fraction",
            counts["univariate"]["failed_or_untested"], counts["univariate"]["family"],
            float(policy["maximum_univariate_untested_fraction"]),
        ),
        (
            "BIVARIATE_FAILURE_FRACTION_EXCEEDED", "bivariate_failure_fraction",
            counts["bivariate"]["failed"], counts["bivariate"]["eligible_family"],
            float(policy["maximum_bivariate_failure_fraction"]),
        ),
        (
            "CONDITIONAL_FAILURE_FRACTION_EXCEEDED", "conditional_failure_fraction",
            counts["conditional"]["failed"], counts["conditional"]["eligible_family"],
            float(policy["maximum_conditional_failure_fraction"]),
        ),
    ]
    reasons = [
        _reason(code, metric, numerator, denominator, threshold)
        for code, metric, numerator, denominator, threshold in gates
        if denominator and _fraction(numerator, denominator) > threshold
    ]
    qc = {
        "thresholds": {
            "maximum_locus_failure_fraction": float(policy["maximum_locus_failure_fraction"]),
            "maximum_univariate_untested_fraction": float(policy["maximum_univariate_untested_fraction"]),
            "maximum_bivariate_failure_fraction": float(policy["maximum_bivariate_failure_fraction"]),
            "maximum_conditional_failure_fraction": float(policy["maximum_conditional_failure_fraction"]),
        },
        "counts": counts,
        "reasons": reasons,
    }
    return context, qc


def _bundle_summary(
    phase: str, index: int | None, source_fingerprint: str,
    continuation_fingerprint: str, engine: Any,
) -> dict[str, object]:
    receipt = engine.validate_bundle(
        phase, index, source_fingerprint, continuation_fingerprint,
        revalidate_semantics=False, reconcile_active=False,
    )
    if receipt is None:
        raise ValidationError(f"missing V2 prerequisite checkpoint: {phase}/{index}")
    bundle = engine.bundle_path(phase, index, continuation_fingerprint)
    target = engine._safe_attempt_target(bundle, continuation_fingerprint)
    result = next(
        (item for item in receipt["artifacts"] if item["path"] == "result.rds"), None,
    )
    if result is None:
        raise ValidationError(f"V2 prerequisite lacks result.rds: {phase}/{index}")
    return {
        "phase": phase,
        "locus_index": index,
        "receipt": _identity(target / "receipt.json"),
        "result": {
            "path": str((target / "result.rds").relative_to(ROOT)),
            "bytes": int(result["bytes"]),
            "sha256": str(result["sha256"]),
        },
        "qc": str(receipt["marker"]["qc"]),
    }


def _upstream_family(
    source_fingerprint: str, continuation_fingerprint: str,
) -> dict[str, object]:
    engine = supervisor()
    aggregate = _bundle_summary(
        "aggregate-discovery", None, source_fingerprint, continuation_fingerprint, engine,
    )
    conditional = [
        _bundle_summary(
            "conditional", index, source_fingerprint, continuation_fingerprint, engine,
        )
        for index in range(1, 2496)
    ]
    compact = [
        {
            "locus_index": item["locus_index"],
            "receipt_sha256": item["receipt"]["sha256"],
            "result_sha256": item["result"]["sha256"],
            "qc": item["qc"],
        }
        for item in conditional
    ]
    return {
        "aggregate": aggregate,
        "conditional_checkpoint_count": len(conditional),
        "conditional_checkpoint_family_sha256": _digest(compact),
        "complete": len(conditional) == 2495,
    }


def _terminal_payload(
    directory: Path, source_fingerprint: str, continuation_fingerprint: str,
) -> dict[str, object]:
    _quick_contract(source_fingerprint, continuation_fingerprint)
    staged = _staged_paths(directory)
    _, qc = _diagnostic_validate(staged)
    if not qc["reasons"]:
        raise ValidationError("terminal QC attestation requested for a family that passes frozen gates")
    upstream = _upstream_family(source_fingerprint, continuation_fingerprint)
    if not upstream["complete"]:
        raise ValidationError("terminal QC attestation requires all 2,495 conditional checkpoints")
    source_lock = CONTRACT.validate_source_lock(
        source_fingerprint, validate_bundles=False, require_live_fingerprint=False,
    )
    return {
        "schema_version": TERMINAL_QC_SCHEMA,
        "state": "TERMINAL_FAILED_QC",
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint,
        "full_family_complete": True,
        "scientific_validation_passed": False,
        "canonical_publication_allowed": False,
        "terminal": True,
        "validation_basis": (
            "UNMODIFIED_LEGACY_SCHEMA_TOPOLOGY_NUMERICAL_AND_FULL_FAMILY_BH_CHECKS_"
            "WITH_QC_THRESHOLDS_DIAGNOSTICALLY_DEFERRED_THEN_EXACT_FROZEN_GATE_EVALUATION"
        ),
        "qc": qc,
        "staged_low_level_results": [_identity(path) for path in staged[:4]],
        "upstream_continuation_family": upstream,
        "source_family_lock": _identity(CONTRACT.source_lock_path(source_fingerprint)),
        "source_checkpoint_family_sha256": source_lock["checkpoint_family_sha256"],
        "lineage": _identity(CONTRACT.lineage_path(continuation_fingerprint)),
        "validators": {
            "legacy_scientific_validator": _identity(LEGACY_VALIDATOR_PATH),
            "continuation_result_validator": _identity(Path(__file__).resolve()),
        },
    }


def validate_staged(
    directory: Path, source_fingerprint: str, continuation_fingerprint: str,
) -> None:
    _quick_contract(source_fingerprint, continuation_fingerprint)
    directory = _attempt_directory(directory, continuation_fingerprint)
    staged = _staged_paths(directory)
    _, qc = _diagnostic_validate(staged)
    if qc["reasons"]:
        payload = _terminal_payload(directory, source_fingerprint, continuation_fingerprint)
        path = directory / "terminal_qc.json"
        encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
        if path.exists():
            if path.read_bytes() != encoded:
                raise ValidationError("existing terminal QC attestation differs; overwrite is forbidden")
        else:
            try:
                CONTRACT.publish_no_replace(path, encoded)
            except CONTRACT.ContinuationError as error:
                raise ValidationError(str(error)) from error
        codes = ",".join(item["code"] for item in payload["qc"]["reasons"])
        print(
            "TRACK_B_LAVA_V2_TERMINAL_FAILED_QC "
            f"source={source_fingerprint} continuation={continuation_fingerprint} "
            f"full_family_complete=TRUE reasons={codes} attestation={path}"
        )
        raise TerminalQCOutcome("complete Track B LAVA family failed frozen QC gates")

    # Passing acceptance always uses the original, byte-frozen validator and
    # its original thresholds.  The diagnostic deferral above is not an
    # acceptance path.
    _legacy_call(LEGACY.validate, staged, build_summary_files=True)
    print(
        "TRACK_B_LAVA_V2_STAGED_RESULTS_SEMANTICALLY_VALIDATED "
        f"source={source_fingerprint} continuation={continuation_fingerprint}"
    )


def verify_terminal_qc(
    directory: Path, source_fingerprint: str, continuation_fingerprint: str,
) -> None:
    directory = _attempt_directory(directory, continuation_fingerprint)
    path = directory / "terminal_qc.json"
    try:
        observed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValidationError(f"terminal QC attestation is unreadable: {error}") from error
    expected = _terminal_payload(directory, source_fingerprint, continuation_fingerprint)
    if observed != expected:
        raise ValidationError("terminal QC attestation differs from current immutable inputs")
    print(
        "TRACK_B_LAVA_V2_TERMINAL_QC_VERIFIED "
        f"state=TERMINAL_FAILED_QC reasons={len(observed['qc']['reasons'])}"
    )


def _finalized_staged(
    source_fingerprint: str, continuation_fingerprint: str,
) -> tuple[Any, Path, list[Path]]:
    _quick_contract(source_fingerprint, continuation_fingerprint)
    engine = supervisor()
    receipt = engine.validate_bundle(
        "finalize", None, source_fingerprint, continuation_fingerprint,
        revalidate_semantics=True, reconcile_active=True,
    )
    if receipt is None:
        raise ValidationError("canonical publication requires a validated finalize READY bundle")
    terminal = engine.validate_bundle(
        "terminal-qc", None, source_fingerprint, continuation_fingerprint,
        revalidate_semantics=False, reconcile_active=False,
    )
    if terminal is not None:
        raise ValidationError("TERMINAL_FAILED_QC family may not be canonically published")
    bundle = engine.bundle_path("finalize", None, continuation_fingerprint)
    target = engine._safe_attempt_target(bundle, continuation_fingerprint)
    return receipt, target, _staged_paths(target)


def _copy_results_no_replace(
    staged: list[Path], identities: dict[Path, tuple[int, str]],
) -> None:
    temporaries: list[Path] = []
    try:
        for index, (source, destination) in enumerate(
            zip(staged, LEGACY.contract.RESULTS, strict=True), start=1,
        ):
            if CONTRACT.file_identity(source) != identities[source]:
                raise ValidationError("staged result changed after semantic validation")
            if destination.exists():
                if CONTRACT.file_identity(destination) != identities[source]:
                    raise ValidationError(f"canonical output differs; overwrite forbidden: {destination}")
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = destination.with_name(
                f".{destination.name}.{os.getpid()}.{index}.v2-publishing"
            )
            with source.open("rb") as source_handle, temporary.open("xb") as target_handle:
                shutil.copyfileobj(source_handle, target_handle, length=1024 * 1024)
                target_handle.flush()
                os.fsync(target_handle.fileno())
            temporaries.append(temporary)
            if CONTRACT.file_identity(temporary) != identities[source]:
                raise ValidationError("canonical publication copy changed in transit")
            try:
                os.link(temporary, destination)
            except FileExistsError as error:
                raise ValidationError(f"canonical output appeared during publication: {destination}") from error
            CONTRACT.fsync_directory(destination.parent)
            temporary.unlink()
            CONTRACT.fsync_directory(destination.parent)
        for source, destination in zip(staged, LEGACY.contract.RESULTS, strict=True):
            if (
                CONTRACT.file_identity(source) != identities[source]
                or CONTRACT.file_identity(destination) != identities[source]
            ):
                raise ValidationError("result changed between validation and sealing")
    finally:
        for path in temporaries:
            path.unlink(missing_ok=True)


def _result_payload(
    source_fingerprint: str, continuation_fingerprint: str, finalize_target: Path,
) -> dict[str, object]:
    engine = supervisor()
    final_receipt = json.loads((finalize_target / "receipt.json").read_text(encoding="utf-8"))
    return {
        "schema_version": RESULT_LOCK_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint,
        "canonical_publication_allowed": True,
        "terminal_qc": None,
        "source_family_lock": _identity(CONTRACT.source_lock_path(source_fingerprint)),
        "lineage": _identity(CONTRACT.lineage_path(continuation_fingerprint)),
        "finalize_bundle": {
            "path": str(engine.bundle_path("finalize", None, continuation_fingerprint).relative_to(ROOT)),
            "receipt": _identity(finalize_target / "receipt.json"),
            "receipt_phase": final_receipt["phase"],
        },
        "results": [_identity(path) for path in LEGACY.contract.RESULTS],
        "result_validator": _identity(Path(__file__).resolve()),
    }


def publish_results(source_fingerprint: str, continuation_fingerprint: str) -> None:
    _, target, staged = _finalized_staged(source_fingerprint, continuation_fingerprint)
    identities = _legacy_call(LEGACY.validate, staged, build_summary_files=False)
    # Revalidate the expensive scientific inputs at the final publication
    # boundary, matching the original fail-closed publication discipline.
    _legacy_call(LEGACY.contract.fully_verify_chromosome_inputs)
    _legacy_call(LEGACY.contract.validate_reference, rehash=True)
    _copy_results_no_replace(staged, identities)
    payload = _result_payload(source_fingerprint, continuation_fingerprint, target)
    encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    lock = LEGACY.contract.RESULT_LOCK
    if lock.exists():
        if lock.read_bytes() != encoded:
            raise ValidationError("canonical result provenance differs; overwrite is forbidden")
    else:
        try:
            CONTRACT.publish_no_replace(lock, encoded)
        except CONTRACT.ContinuationError as error:
            raise ValidationError(str(error)) from error
    print(
        "TRACK_B_LAVA_V2_RESULTS_VALIDATED_AND_PUBLISHED "
        f"source={source_fingerprint} continuation={continuation_fingerprint}"
    )


def verify_results(source_fingerprint: str, continuation_fingerprint: str) -> None:
    _, target, _ = _finalized_staged(source_fingerprint, continuation_fingerprint)
    _legacy_call(LEGACY.validate, LEGACY.contract.RESULTS, build_summary_files=False)
    expected = _result_payload(source_fingerprint, continuation_fingerprint, target)
    try:
        observed = json.loads(LEGACY.contract.RESULT_LOCK.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValidationError(f"canonical V2 result provenance is unreadable: {error}") from error
    if observed != expected:
        raise ValidationError("canonical V2 results differ from immutable dual-fingerprint provenance")
    print("TRACK_B_LAVA_V2_RESULTS_VERIFIED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--validate-staged-results", action="store_true")
    mode.add_argument("--verify-staged-results", action="store_true")
    mode.add_argument("--verify-terminal-qc", action="store_true")
    mode.add_argument("--publish-staged-results", action="store_true")
    mode.add_argument("--verify-results", action="store_true")
    parser.add_argument("--staging-dir", type=Path)
    parser.add_argument("--terminal-qc-output", type=Path)
    parser.add_argument("--source-fingerprint", required=True)
    parser.add_argument("--continuation-fingerprint", required=True)
    args = parser.parse_args()

    if args.validate_staged_results or args.verify_staged_results or args.verify_terminal_qc:
        if args.staging_dir is None:
            parser.error("selected staging mode requires --staging-dir")
        directory = _attempt_directory(args.staging_dir, args.continuation_fingerprint)
        if args.terminal_qc_output is not None:
            expected_terminal = directory / "terminal_qc.json"
            if args.terminal_qc_output != expected_terminal:
                try:
                    supplied_terminal = args.terminal_qc_output.resolve(strict=False)
                except (OSError, RuntimeError) as error:
                    raise ValidationError(f"invalid terminal-QC output path: {error}") from error
                if supplied_terminal != expected_terminal:
                    raise ValidationError(
                        "--terminal-qc-output must name terminal_qc.json in the supervised attempt"
                    )
        if args.validate_staged_results:
            validate_staged(directory, args.source_fingerprint, args.continuation_fingerprint)
        elif args.verify_terminal_qc:
            verify_terminal_qc(directory, args.source_fingerprint, args.continuation_fingerprint)
        else:
            _quick_contract(args.source_fingerprint, args.continuation_fingerprint)
            staged = _staged_paths(directory)
            _legacy_call(LEGACY.validate, staged, build_summary_files=False)
            print("TRACK_B_LAVA_V2_STAGED_RESULTS_READ_ONLY_VERIFIED")
    elif args.publish_staged_results:
        if args.staging_dir is not None or args.terminal_qc_output is not None:
            parser.error("publication uses only the immutable finalize READY bundle")
        publish_results(args.source_fingerprint, args.continuation_fingerprint)
    else:
        if args.staging_dir is not None or args.terminal_qc_output is not None:
            parser.error("canonical verification does not accept staging paths")
        verify_results(args.source_fingerprint, args.continuation_fingerprint)


if __name__ == "__main__":
    try:
        main()
    except TerminalQCOutcome as error:
        print(f"TRACK_B_LAVA_V2_TERMINAL_FAILED_QC message={error}", file=sys.stderr)
        raise SystemExit(TERMINAL_EXIT_STATUS) from error
    except (ValidationError, CONTRACT.ContinuationError) as error:
        raise SystemExit(f"ERROR: {error}") from error
