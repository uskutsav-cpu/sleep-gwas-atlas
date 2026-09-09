#!/usr/bin/env python3
"""Seal and verify the additive post-LAVA Track B PLACO+ execution gate.

This continuation never changes the frozen V1 pleiotropy policy or its input
gate.  It recognizes only a deeply validated LAVA V2 TERMINAL_FAILED_QC bundle
as satisfying the V1 ordering requirement.  The failed LAVA family is ordering
evidence only and is never an input to PLACO+.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import fcntl
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import stat
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
POLICY = Path("config/track_b_pleiotropy_policy.json")
V1_CONTRACT_LOCK = Path("results/track_b/pleiotropy/contract.lock.json")
V1_INPUT_GATE = Path("results/track_b/pleiotropy/input_gate.tsv")
V1_READINESS_GATE = Path("results/track_b/pleiotropy/readiness_gate.tsv")
V1_INPUT_LOCK = Path("results/track_b/pleiotropy/input_gate.lock.json")
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
PAIR_MANIFEST_LOCK = Path("results/track_b/pair_manifest.lock.json")
DENSE_QC = Path("results/track_b/03_dense_input_qc.tsv")
DENSE_QC_LOCK = Path("results/track_b/03_dense_input_qc.lock.json")
GATE_ROOT = Path("results/track_b/pleiotropy/continuations/post_lava_terminal_v2")
READINESS_GATE = GATE_ROOT / "readiness_gate.tsv"
TERMINAL_GATE_LOCK = GATE_ROOT / "terminal_gate.lock.json"
GATE_PUBLICATION_FLOCK = GATE_ROOT / "publication.flock"
V1_MATERIALIZER = Path("scripts/125_materialize_track_b_placo_pair.py")
V1_RUNNER = Path("scripts/126_run_track_b_placo_pair.R")
V2_RUNNER = Path("scripts/141_run_track_b_placo_pair_v2.R")
V1_CONTRACT_SCRIPT = Path("scripts/123_track_b_pleiotropy_contract.py")
V1_INPUT_GATE_SCRIPT = Path("scripts/124_build_track_b_pleiotropy_input_gate.py")
V1_LD_MATERIALIZER = Path("scripts/127_prepare_track_b_pleiotropy_ld.py")
R_RUNTIME = Path(".r-env/bin/Rscript")
GATE_SCRIPT = Path("scripts/139_build_track_b_placo_terminal_gate_v2.py")
BRIDGE_SCRIPT = Path("scripts/140_materialize_track_b_placo_pair_v2.py")
ARCHIVE_STATE_GUARD_SCRIPT = Path("scripts/146_manage_lava_reference_archives.py")
LAVA_RESULT_LOCK = Path("results/track_b/local/lava_results.provenance.json")
LAVA_CANONICAL_RESULTS = [
    Path("results/track_b/local/lava_locus_status.tsv"),
    Path("results/track_b/local/lava_univariate.tsv"),
    Path("results/track_b/local/lava_bivariate.tsv"),
    Path("results/track_b/local/lava_conditional.tsv"),
    Path("results/track_b/04_lava_local_results.tsv"),
    Path("results/track_b/05_local_conditional_results.tsv"),
]
PLACO_WORK_ROOT = Path("work/track_b_pleiotropy/placo")
LAVA_V2_CODE = [
    Path("scripts/134_track_b_lava_continuation_contract.py"),
    Path("scripts/135_run_track_b_lava_postdiscovery_v2.R"),
    Path("scripts/136_run_track_b_lava_continuation.py"),
    Path("scripts/137_validate_track_b_lava_checkpoint_v2.R"),
    Path("scripts/138_validate_track_b_lava_results_v2.py"),
]

GATE_SCHEMA = "sleep-atlas-track-b-pleiotropy-terminal-gate.2"
READINESS_SCHEMA = "sleep-atlas-track-b-pleiotropy-terminal-readiness.2"
RAM_SAFETY_RESERVE_BYTES = 1 * 1024**3

READINESS_FIELDS = [
    "schema_version", "component_id", "pair_scope", "dense_input_gate",
    "replication_gate", "local_analysis_gate", "lava_scientific_evidence",
    "implementation_gate", "software_gate", "full_p_scan_reference_gate",
    "compute_gate", "scan_gate", "locus_publication_gate", "overall_status",
    "blockers", "claim_status",
]
V1_INPUT_GATE_FIELDS = [
    "pair_id", "family_role", "trait1", "trait2", "ancestry", "analysis_build",
    "trait1_file", "trait1_bytes", "trait1_sha256", "trait1_rows", "trait1_schema",
    "trait1_INFO_status", "trait2_file", "trait2_bytes", "trait2_sha256", "trait2_rows",
    "trait2_schema", "trait2_INFO_status", "minimum_rows_per_trait",
    "full_genome_wide_dense_required", "hapmap3_only_forbidden", "gzip_integrity",
    "minimum_aligned_eligible_variants", "minimum_fraction_of_smaller_dense_input",
    "required_autosomes", "per_autosome_provenance_required", "pair_alignment_status",
    "input_gate_status", "blocker",
]
V1_INPUT_LOCK_FIELDS = {
    "schema_version", "analysis_id", "selection_timing",
    "pleiotropy_results_accessed_before_input_gate_freeze", "policy_sha256",
    "contract_lock_sha256", "contract_script_sha256", "input_gate_script_sha256",
    "ld_materializer_script_sha256", "pair_manifest_sha256", "pair_manifest_lock_sha256",
    "dense_qc_sha256", "dense_qc_lock_sha256", "input_gate_sha256",
    "readiness_gate_sha256", "pair_count", "primary_pair_count", "control_pair_count",
    "unique_dense_trait_count", "all_pair_dense_input_gates_pass",
    "full_gzip_integrity_scans_completed", "live_dense_inputs", "scientific_result_count",
    "scientific_result_substitution_policy",
}


class GateError(RuntimeError):
    """A fail-closed V2 terminal-gate violation."""


def load_module(name: str, relative: str):
    path = ROOT / relative
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise GateError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


PLEIOTROPY = load_module("track_b_pleiotropy_contract_for_terminal_gate", "scripts/123_track_b_pleiotropy_contract.py")
LAVA_CONTRACT = load_module("track_b_lava_continuation_for_placo_gate", "scripts/134_track_b_lava_continuation_contract.py")
LAVA_SUPERVISOR = load_module("track_b_lava_supervisor_for_placo_gate", "scripts/136_run_track_b_lava_continuation.py")
ARCHIVE_STATE_GUARD = load_module(
    "track_b_lava_reference_archive_guard_for_placo_gate",
    "scripts/146_manage_lava_reference_archives.py",
)

ARCHIVES_PRESENT = "ARCHIVES_PRESENT_FULLY_VERIFIED"
ARCHIVES_EVICTED = (
    "ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_"
    "EXTRACTED_PAYLOADS_REHASHED"
)
ACCEPTED_ARCHIVE_STATES = {ARCHIVES_PRESENT, ARCHIVES_EVICTED}


def sha256(path: Path) -> str:
    return file_identity(path)[1]


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def require_real_repository_path(path: Path, label: str) -> Path:
    """Reject path escape and every symlink in an existing in-repository prefix."""
    try:
        relative = path.relative_to(ROOT)
        root_resolved = ROOT.resolve(strict=True)
    except (OSError, RuntimeError, ValueError) as error:
        raise GateError(f"{label} escapes repository: {path}") from error
    current = ROOT
    for offset, part in enumerate(relative.parts):
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise GateError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise GateError(f"{label} contains a symbolic link: {current}")
        if offset < len(relative.parts) - 1 and not stat.S_ISDIR(observed.st_mode):
            raise GateError(f"{label} has a non-directory ancestor: {current}")
    try:
        path.resolve(strict=False).relative_to(root_resolved)
    except (OSError, RuntimeError, ValueError) as error:
        raise GateError(f"{label} escapes repository: {path}") from error
    return path


def file_identity(path: Path) -> tuple[int, str]:
    require_real_repository_path(path, "artifact")
    digest = hashlib.sha256()
    try:
        path_stat = os.lstat(path)
        if not stat.S_ISREG(path_stat.st_mode) or path_stat.st_size <= 0:
            raise GateError(f"artifact is not a real non-empty file: {path}")
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise GateError(f"artifact is not a regular file: {path}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = os.lstat(path)
    except OSError as error:
        raise GateError(f"could not hash artifact {path}: {error}") from error
    if (
        not stat.S_ISREG(current.st_mode)
        or _stat_identity(path_stat) != _stat_identity(before)
        or _stat_identity(before) != _stat_identity(after)
        or _stat_identity(after) != _stat_identity(current)
    ):
        raise GateError(f"artifact changed while hashing or became a symlink: {path}")
    return int(after.st_size), digest.hexdigest()


def identity(path: Path) -> dict[str, object]:
    size, digest = file_identity(path)
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": size,
        "sha256": digest,
    }


def read_json(path: Path) -> dict[str, Any]:
    require_real_repository_path(path, "JSON artifact")
    if path.is_symlink():
        raise GateError(f"JSON artifact may not be a symbolic link: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise GateError(f"unreadable JSON artifact {path}: {error}") from error
    if not isinstance(value, dict):
        raise GateError(f"JSON artifact must contain one object: {path}")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require_real_repository_path(path, "TSV artifact")
    if path.is_symlink():
        raise GateError(f"TSV artifact may not be a symbolic link: {path}")
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise GateError(f"unreadable TSV artifact {path}: {error}") from error
    if not fields or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise GateError(f"malformed TSV artifact: {path}")
    return fields, rows


def tsv_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def canonical_lava_result_evidence() -> list[Path]:
    evidence: list[Path] = []
    for relative in [*LAVA_CANONICAL_RESULTS, LAVA_RESULT_LOCK]:
        path = ROOT / relative
        require_real_repository_path(path, "canonical LAVA result path")
        if path.exists() or path.is_symlink():
            evidence.append(path)
    return evidence


def assert_no_canonical_lava_results() -> None:
    evidence = canonical_lava_result_evidence()
    if evidence:
        names = ", ".join(str(path.relative_to(ROOT)) for path in evidence[:3])
        raise GateError(
            "TERMINAL_FAILED_QC cannot coexist with a partial or sealed canonical LAVA result family: "
            + names
        )


def placo_result_evidence(policy: dict[str, Any]) -> list[Path]:
    """Find canonical and resumable scientific PLACO artifacts before gate freeze."""
    try:
        evidence = list(PLEIOTROPY.result_evidence(ROOT, policy))
    except SystemExit as error:
        raise GateError(str(error)) from error
    work = ROOT / PLACO_WORK_ROOT
    if work.is_symlink():
        return evidence + [work]
    if not work.exists():
        return evidence
    require_real_repository_path(work, "PLACO work namespace")
    if not work.is_dir():
        return evidence + [work]
    scientific_names = {
        "benchmark.raw.tsv", "benchmark.tsv", "nuisance.rds", "nuisance.sha256.tsv",
        "staged.full.tsv.gz", "run.summary.tsv", "staged.provenance.json",
    }

    def scientific_name(name: str) -> bool:
        return any(
            name == frozen
            or name.startswith(frozen + ".")
            or name.startswith("." + frozen + ".")
            for frozen in scientific_names
        )

    for path in sorted(work.rglob("*")):
        if path.is_symlink():
            evidence.append(path)
            continue
        relative = path.relative_to(work)
        if path.is_file() and (
            scientific_name(path.name)
            or "checkpoints" in relative.parts
            or path.name.startswith("shard_")
        ):
            evidence.append(path)
    return evidence


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


@contextmanager
def gate_publication_lock():
    path = ROOT / GATE_PUBLICATION_FLOCK
    require_real_repository_path(path.parent, "PLACO V2 gate publication-lock parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        before = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(before.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (before.st_dev, before.st_ino) != (current.st_dev, current.st_ino)
            or current.st_nlink != 1
        ):
            raise GateError("PLACO V2 gate publication flock is not a private regular file")
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        locked = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(locked.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (locked.st_dev, locked.st_ino) != (current.st_dev, current.st_ino)
            or current.st_nlink != 1
        ):
            raise GateError("PLACO V2 gate publication flock changed during acquisition")
        yield
    finally:
        os.close(descriptor)


def reconcile_gate_temporaries(encoded: list[bytes]) -> None:
    gate_root = ROOT / GATE_ROOT
    if not gate_root.exists():
        return
    require_real_repository_path(gate_root, "PLACO V2 gate namespace")
    patterns = {
        READINESS_GATE.name: (ROOT / READINESS_GATE, encoded[0]),
        TERMINAL_GATE_LOCK.name: (ROOT / TERMINAL_GATE_LOCK, encoded[1]),
    }
    stale: list[dict[str, object]] = []
    for path in sorted(gate_root.iterdir()):
        matched: tuple[Path, bytes] | None = None
        for basename, target in patterns.items():
            prefix = f".{basename}."
            if not path.name.startswith(prefix):
                continue
            if re.fullmatch(re.escape(prefix) + r"[0-9]+\.tmp", path.name) is None:
                raise GateError(f"malformed PLACO V2 gate temporary name: {path.name}")
            matched = target
            break
        if matched is None:
            continue
        target_path, wanted = matched
        observed = os.lstat(path)
        if not stat.S_ISREG(observed.st_mode):
            raise GateError(f"PLACO V2 gate temporary is not a regular artifact: {path}")
        content = path.read_bytes()
        current = os.lstat(path)
        if (current.st_dev, current.st_ino, current.st_size) != (
            observed.st_dev, observed.st_ino, observed.st_size,
        ):
            raise GateError("PLACO V2 gate temporary changed while hashing")
        size, digest = len(content), hashlib.sha256(content).hexdigest()
        if target_path.exists() or target_path.is_symlink():
            target_size, target_hash = file_identity(target_path)
            if target_size != len(wanted) or target_hash != hashlib.sha256(wanted).hexdigest():
                raise GateError("gate temporary coexists with a noncanonical final member")
        stale.append({
            "path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest,
            "device": int(observed.st_dev), "inode": int(observed.st_ino),
            "target": str(target_path.relative_to(ROOT)),
        })
    if not stale:
        return
    receipt = {
        "schema_version": "sleep-atlas-track-b-placo-gate-temp-reconciliation.1",
        "artifacts": stale,
        "rule": "ONLY_EXACT_EXPECTED_PID_TEMPORARIES_UNLINKED_UNDER_GATE_PUBLICATION_FLOCK",
        "scientific_result": False,
    }
    content = json.dumps(receipt, indent=2, sort_keys=True).encode() + b"\n"
    digest = hashlib.sha256(content).hexdigest()
    audit = gate_root / "publication_reconciliation" / f"gate_temporaries.{digest}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    require_real_repository_path(audit.parent, "PLACO V2 gate reconciliation audit parent")
    if audit.exists() or audit.is_symlink():
        if file_identity(audit) != (len(content), digest):
            raise GateError("PLACO V2 gate reconciliation audit drifted")
    else:
        prefix = f".{audit.name}."
        for stale_audit in sorted(audit.parent.iterdir()):
            if not stale_audit.name.startswith(prefix) or not stale_audit.name.endswith(".tmp"):
                continue
            if re.fullmatch(re.escape(prefix) + r"[0-9a-f]{64}\.[0-9]+\.tmp", stale_audit.name) is None:
                raise GateError("malformed gate-reconciliation audit temporary")
            observed = os.lstat(stale_audit)
            if not stat.S_ISREG(observed.st_mode):
                raise GateError("gate-reconciliation audit temporary is not regular")
            stale_audit.unlink()
        temporary = audit.with_name(f".{audit.name}.{digest}.{os.getpid()}.tmp")
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        if hashlib.sha256(temporary.read_bytes()).hexdigest() != digest:
            raise GateError("gate-reconciliation audit temporary changed before commit")
        os.rename(temporary, audit)
        fsync_directory(audit.parent)
    for record in stale:
        path = ROOT / str(record["path"])
        current = os.lstat(path)
        if (current.st_dev, current.st_ino, current.st_size) != (
            record["device"], record["inode"], record["bytes"],
        ):
            raise GateError("PLACO V2 gate temporary changed after reconciliation audit")
        path.unlink()
        fsync_directory(path.parent)


def unlink_if_identity(path: Path, expected: tuple[int, int]) -> None:
    """Rollback only the exact inode published by this process."""
    try:
        observed = os.lstat(path)
    except FileNotFoundError:
        return
    except OSError as error:
        raise GateError(f"could not inspect rollback artifact {path}: {error}") from error
    if (observed.st_dev, observed.st_ino) == expected:
        path.unlink()
        fsync_directory(path.parent)


def validate_terminal_payload(payload: dict[str, Any]) -> None:
    """Require a complete failed-QC family, never a successful LAVA result."""
    if (
        payload.get("schema_version") != "track-b-lava-terminal-qc.1"
        or payload.get("state") != "TERMINAL_FAILED_QC"
        or payload.get("analysis_id") != "track-b-v1.0-local"
        or payload.get("full_family_complete") is not True
        or payload.get("scientific_validation_passed") is not False
        or payload.get("canonical_publication_allowed") is not False
        or payload.get("terminal") is not True
    ):
        raise GateError("LAVA terminal attestation does not encode a complete TERMINAL_FAILED_QC family")
    qc = payload.get("qc")
    if not isinstance(qc, dict) or not isinstance(qc.get("reasons"), list) or not qc["reasons"]:
        raise GateError("LAVA terminal attestation lacks one or more frozen-QC failure reasons")
    counts = qc.get("counts")
    thresholds = qc.get("thresholds")
    if not isinstance(counts, dict) or not isinstance(thresholds, dict):
        raise GateError("LAVA terminal attestation lacks exact QC counts or thresholds")
    locus = counts.get("locus")
    univariate = counts.get("univariate")
    if not isinstance(locus, dict) or not isinstance(univariate, dict):
        raise GateError("LAVA terminal attestation has malformed frozen-family counts")
    try:
        locus_family = int(locus.get("family", -1))
        univariate_family = int(univariate.get("family", -1))
        maximum_locus = float(thresholds.get("maximum_locus_failure_fraction", -1))
        maximum_univariate = float(thresholds.get("maximum_univariate_untested_fraction", -1))
    except (TypeError, ValueError) as error:
        raise GateError("LAVA terminal attestation has nonnumeric frozen-family metadata") from error
    if (
        locus_family != 2495
        or univariate_family != 19960
        or not math.isclose(maximum_locus, 0.01)
        or not math.isclose(maximum_univariate, 0.05)
    ):
        raise GateError("LAVA terminal attestation differs from the frozen local-analysis family")


def validate_v1_pleiotropy() -> dict[str, Any]:
    try:
        policy = PLEIOTROPY.validate_policy(ROOT)
        contract_lock = PLEIOTROPY.verify_contract(ROOT)
    except SystemExit as error:
        raise GateError(str(error)) from error
    lock = read_json(ROOT / V1_INPUT_LOCK)
    if (
        set(lock) != V1_INPUT_LOCK_FIELDS
        or lock.get("schema_version") != "sleep-atlas-track-b-pleiotropy-input-gate.1"
        or lock.get("analysis_id") != policy["analysis_id"]
        or lock.get("selection_timing") != policy["selection_timing"]
        or lock.get("policy_sha256") != sha256(ROOT / POLICY)
        or lock.get("contract_lock_sha256") != sha256(ROOT / V1_CONTRACT_LOCK)
        or lock.get("contract_script_sha256")
        != contract_lock["script_sha256"][str(PLEIOTROPY.CONTRACT_SCRIPT)]
        or lock.get("input_gate_script_sha256")
        != contract_lock["script_sha256"][str(PLEIOTROPY.INPUT_GATE_SCRIPT)]
        or lock.get("ld_materializer_script_sha256")
        != contract_lock["script_sha256"][str(PLEIOTROPY.LD_MATERIALIZER_SCRIPT)]
        or lock.get("input_gate_sha256") != sha256(ROOT / V1_INPUT_GATE)
        or lock.get("readiness_gate_sha256") != sha256(ROOT / V1_READINESS_GATE)
        or lock.get("pair_count") != 3
        or lock.get("primary_pair_count") != 2
        or lock.get("control_pair_count") != 1
        or lock.get("unique_dense_trait_count") != 5
        or lock.get("all_pair_dense_input_gates_pass") is not True
        or lock.get("full_gzip_integrity_scans_completed") != 5
        or lock.get("scientific_result_count") != 0
        or lock.get("pleiotropy_results_accessed_before_input_gate_freeze") is not False
        or lock.get("scientific_result_substitution_policy")
        != "FORBIDDEN; readiness gates are not pleiotropy results"
    ):
        raise GateError("frozen V1 pleiotropy input gate is incomplete or drifted")
    upstream_identities = {
        "pair_manifest": identity(ROOT / PAIR_MANIFEST),
        "pair_manifest_lock": identity(ROOT / PAIR_MANIFEST_LOCK),
        "dense_qc": identity(ROOT / DENSE_QC),
        "dense_qc_lock": identity(ROOT / DENSE_QC_LOCK),
    }
    for name, record in upstream_identities.items():
        if record["sha256"] != lock.get(f"{name}_sha256"):
            raise GateError(f"frozen V1 upstream identity differs from input lock: {name}")
    _, readiness = read_tsv(ROOT / V1_READINESS_GATE)
    placo_rows = [row for row in readiness if row.get("component_id") == "PLACO_PLUS"]
    if len(placo_rows) != 1:
        raise GateError("frozen V1 readiness lacks exactly one PLACO_PLUS row")
    row = placo_rows[0]
    if (
        row.get("dense_input_gate") != "PASS_ALL_THREE_PAIRS"
        or not row.get("replication_gate", "").startswith("PASS_TERMINAL")
        or row.get("local_analysis_gate") != "BLOCKED_LOCAL_ANALYSIS_NOT_TERMINAL"
        or row.get("software_gate") != "READY_PINNED_PLACO_SOURCE"
        or row.get("claim_status") != "NO_SCIENTIFIC_RESULT"
    ):
        raise GateError("frozen V1 readiness is not the expected pre-LAVA blocked snapshot")
    dense = lock.get("live_dense_inputs")
    if not isinstance(dense, dict) or set(dense) != set(PLEIOTROPY.EXPECTED_TRAITS):
        raise GateError("V1 input lock lacks the exact five-trait dense family")
    dense_identities: dict[str, dict[str, object]] = {}
    for trait in PLEIOTROPY.EXPECTED_TRAITS:
        record = dense[trait]
        if (
            not isinstance(record, dict)
            or record.get("schema") != PLEIOTROPY.EXPECTED_DENSE_SCHEMA
            or record.get("gzip_integrity") != "FULL_DECOMPRESSION_CRC_PASS"
        ):
            raise GateError(f"V1 input lock has malformed dense-source metadata: {trait}")
        try:
            path = PLEIOTROPY.root_path(ROOT, str(record.get("path", "")))
            rows = int(record["rows"])
        except (SystemExit, KeyError, TypeError, ValueError) as error:
            raise GateError(f"V1 input lock has unsafe or invalid dense-source metadata: {trait}") from error
        observed = identity(path)
        if observed["bytes"] != record.get("bytes") or observed["sha256"] != record.get("sha256"):
            raise GateError(f"live dense source differs from V1 frozen identity: {trait}")
        dense_identities[trait] = observed | {
            "rows": rows, "schema": str(record["schema"]),
        }
    input_fields, pair_rows = read_tsv(ROOT / V1_INPUT_GATE)
    if (
        input_fields != V1_INPUT_GATE_FIELDS
        or len(pair_rows) != 3
        or [item.get("pair_id") for item in pair_rows] != [item[0] for item in PLEIOTROPY.EXPECTED_PAIRS]
    ):
        raise GateError("frozen V1 pair input gate schema, count, or order drifted")
    floors = policy["dense_input_contract"]["pair_materialization_plausibility"]
    for pair_row, (pair_id, trait1, trait2, role) in zip(
        pair_rows, PLEIOTROPY.EXPECTED_PAIRS, strict=True,
    ):
        source1, source2 = dense[trait1], dense[trait2]
        expected = {
            "pair_id": pair_id, "family_role": role, "trait1": trait1, "trait2": trait2,
            "ancestry": "EUR", "analysis_build": "hg19",
            "trait1_file": source1["path"], "trait1_bytes": str(source1["bytes"]),
            "trait1_sha256": source1["sha256"], "trait1_rows": str(source1["rows"]),
            "trait1_schema": source1["schema"], "trait2_file": source2["path"],
            "trait2_bytes": str(source2["bytes"]), "trait2_sha256": source2["sha256"],
            "trait2_rows": str(source2["rows"]), "trait2_schema": source2["schema"],
            "minimum_rows_per_trait": str(policy["dense_input_contract"]["minimum_variant_rows_per_trait"]),
            "full_genome_wide_dense_required": "TRUE", "hapmap3_only_forbidden": "TRUE",
            "gzip_integrity": "FULL_DECOMPRESSION_CRC_PASS_BOTH",
            "minimum_aligned_eligible_variants": str(floors["minimum_aligned_eligible_variants"]),
            "minimum_fraction_of_smaller_dense_input": str(floors["minimum_fraction_of_smaller_dense_input"]),
            "required_autosomes": ",".join(str(value) for value in range(1, 23)),
            "per_autosome_provenance_required": "TRUE",
            "pair_alignment_status": (
                "NOT_RUN_PRE_RESULT;EXACT_GLOBAL_AND_22_AUTOSOME_COUNTS_REQUIRED_AT_MATERIALIZATION"
            ),
            "input_gate_status": "READY_DENSE_INPUTS", "blocker": "NONE",
        }
        if any(pair_row.get(key) != value for key, value in expected.items()):
            raise GateError(f"frozen V1 pair input gate differs from exact identity: {pair_id}")
    placo_source = ROOT / policy["placo_plus"]["source_path"]
    if sha256(placo_source) != policy["placo_plus"]["source_sha256"]:
        raise GateError("pinned PLACO+ source differs from frozen policy")
    return {
        "policy": policy,
        "v1_input_lock": lock,
        "v1_readiness_row": row,
        "upstream_identities": upstream_identities,
        "dense_identities": dense_identities,
        "placo_source": identity(placo_source),
    }


def verified_archive_state() -> dict[str, Any]:
    """Rehash the scientific reference and accept only the two audited states."""

    try:
        observed = ARCHIVE_STATE_GUARD.verify_reference_state(ROOT)
    except (ARCHIVE_STATE_GUARD.ArchiveEvictionError, OSError) as error:
        raise GateError(f"LAVA reference archive state failed verification: {error}") from error
    if not isinstance(observed, dict) or observed.get("state") not in ACCEPTED_ARCHIVE_STATES:
        raise GateError("LAVA reference archive guard returned a non-permitted state")
    if (
        observed.get("extracted_payload_count") != 44
        or not re.fullmatch(r"[0-9a-f]{64}", str(observed.get("archive_family_sha256", "")))
        or not re.fullmatch(r"[0-9a-f]{64}", str(observed.get("extracted_family_sha256", "")))
        or not isinstance(observed.get("reference_provenance_identity"), dict)
    ):
        raise GateError("LAVA reference archive guard returned incomplete scientific identity")
    return observed


def validated_lava_continuation(reference_state: dict[str, Any]) -> tuple[str, str]:
    """Validate the frozen lineage without asking an evicted ZIP to be live."""

    source = LAVA_CONTRACT.PINNED_SOURCE_FINGERPRINT
    state = reference_state.get("state")
    if state == ARCHIVES_PRESENT:
        continuation = LAVA_CONTRACT.continuation_fingerprint(
            source, validate_source_bundles=True,
        )
        LAVA_CONTRACT.validate_lineage(source, continuation, deep=True)
    elif state == ARCHIVES_EVICTED:
        # The final eviction receipt and all 44 extracted payloads were deeply
        # checked above.  Source checkpoints remain fully and semantically
        # revalidated, but the legacy live fingerprint is intentionally not
        # invoked because it includes the redundant acquisition ZIP family.
        LAVA_CONTRACT.validate_source_lock(
            source, validate_bundles=True, revalidate_semantics=True,
            require_live_fingerprint=False,
        )
        continuation = LAVA_CONTRACT.continuation_fingerprint(
            source, validate_source_bundles=False,
        )
        LAVA_CONTRACT.validate_lineage_record(source, continuation)
    else:
        raise GateError("LAVA continuation requested for a non-permitted archive state")
    return source, continuation


def archive_state_invariant(reference_state: dict[str, Any]) -> dict[str, Any]:
    """Return only identities that are invariant across present/evicted states."""

    return {
        "accepted_states": [ARCHIVES_PRESENT, ARCHIVES_EVICTED],
        "guard_implementation": identity(ROOT / ARCHIVE_STATE_GUARD_SCRIPT),
        "reference_provenance": reference_state["reference_provenance_identity"],
        "archive_family_sha256": reference_state["archive_family_sha256"],
        "extracted_payload_count": reference_state["extracted_payload_count"],
        "extracted_family_sha256": reference_state["extracted_family_sha256"],
        "verification": (
            "CURRENT_STATE_REVALIDATED_EVERY_CALL;MUTABLE_PRESENT_OR_EVICTED_TOKEN_"
            "INTENTIONALLY_EXCLUDED_FROM_THE_IMMUTABLE_TERMINAL_GATE"
        ),
    }


def deep_lava_terminal() -> dict[str, Any]:
    reference_state = verified_archive_state()
    try:
        source, continuation = validated_lava_continuation(reference_state)
        receipt = LAVA_SUPERVISOR.validate_bundle(
            "terminal-qc", None, source, continuation,
            revalidate_semantics=True, reconcile_active=True,
        )
    except (LAVA_CONTRACT.ContinuationError, LAVA_SUPERVISOR.ExecutionError, OSError) as error:
        raise GateError(f"LAVA V2 terminal bundle failed deep validation: {error}") from error
    if receipt is None:
        raise GateError("LAVA V2 TERMINAL_FAILED_QC READY bundle does not yet exist")
    if (
        receipt.get("phase") != "terminal-qc"
        or receipt.get("exit_status") != 78
        or receipt.get("marker", {}).get("qc") != "TERMINAL_FAILED_QC"
    ):
        raise GateError("LAVA V2 terminal receipt has the wrong phase, exit status, or state")
    bundle = LAVA_SUPERVISOR.bundle_path("terminal-qc", None, continuation)
    try:
        target = LAVA_SUPERVISOR._safe_attempt_target(bundle, continuation)
    except LAVA_SUPERVISOR.ExecutionError as error:
        raise GateError(str(error)) from error
    terminal_path = target / "terminal_qc.json"
    terminal = read_json(terminal_path)
    validate_terminal_payload(terminal)
    if (
        terminal.get("source_discovery_fingerprint") != source
        or terminal.get("continuation_execution_fingerprint") != continuation
    ):
        raise GateError("LAVA terminal attestation fingerprint lineage drifted")
    assert_no_canonical_lava_results()
    artifacts = [identity(path) for path in sorted(target.iterdir()) if path.is_file()]
    return {
        "source_discovery_fingerprint": source,
        "continuation_execution_fingerprint": continuation,
        "ready_link": str(bundle.relative_to(ROOT)),
        "ready_target": str(target.relative_to(ROOT)),
        "ready_bundle_artifacts": artifacts,
        "terminal_attestation": identity(terminal_path),
        "terminal_payload": terminal,
        "source_family_lock": identity(LAVA_CONTRACT.source_lock_path(source)),
        "lineage": identity(LAVA_CONTRACT.lineage_path(continuation)),
        "reference_archive_state_guard": archive_state_invariant(reference_state),
        # The deeply revalidated receipt can contain thousands of active-input
        # records.  Its immutable byte identity binds that entire family
        # without copying it into a second oversized lock.
        "receipt": identity(target / "receipt.json"),
        "receipt_active_input_count": len(receipt["active_inputs"]),
        "continuation_code": {
            str(path): identity(ROOT / path) for path in LAVA_V2_CODE
        },
    }


def readiness_text() -> str:
    row: dict[str, object] = {
        "schema_version": READINESS_SCHEMA,
        "component_id": "PLACO_PLUS_FULL_P_V2",
        "pair_scope": "A;B;CONTROL",
        "dense_input_gate": "PASS_INHERITED_FROZEN_FIVE_TRAIT_GATE",
        "replication_gate": "PASS_TERMINAL_PRIMARY_REPLICATION_FAMILY",
        "local_analysis_gate": "PASS_TERMINAL_FAILED_QC_ATTESTED",
        "lava_scientific_evidence": "FAILED_QC_NOT_CONSUMED",
        "implementation_gate": "READY_PINNED_MATERIALIZER_ENGINE_RUNNER_AND_V2_BRIDGE",
        "software_gate": "READY_PINNED_PLACO_PLUS_0.2.0",
        "full_p_scan_reference_gate": "NOT_REQUIRED_FOR_FULL_P_SCAN",
        "compute_gate": "PER_PAIR_MATERIALIZATION_PREFLIGHT_AND_MATCHED_WORKER_BENCHMARK_REQUIRED",
        "scan_gate": "READY_CONDITIONAL_FULL_P_ONLY",
        "locus_publication_gate": "BLOCKED_FULL_EUR_HG19_LD_AND_FAMILY_COLLATOR",
        "overall_status": "READY_FOR_MATERIALIZATION_BENCHMARK_AND_CONDITIONAL_FULL_P_SCAN",
        "blockers": "LD_REFERENCE_AND_COLLATOR_FOR_LOCUS_PUBLICATION",
        "claim_status": "FULL_P_METHOD_RESULT_ONLY;NO_LAVA_PASS;NO_LD_LOCUS_CLAIM",
    }
    return tsv_text(READINESS_FIELDS, [row])


def lock_payload(v1: dict[str, Any], lava: dict[str, Any], readiness: str) -> dict[str, Any]:
    policy = v1["policy"]
    return {
        "schema_version": GATE_SCHEMA,
        "analysis_id": policy["analysis_id"],
        "selection_timing": "AFTER_LAVA_TERMINAL_ATTESTATION_BEFORE_ANY_PLACO_RESULT_ACCESS",
        "pleiotropy_results_accessed_before_terminal_gate_freeze": False,
        "pair_scope": ["A", "B", "CONTROL"],
        "permitted_scope": "FULL_P_MATERIALIZATION_BENCHMARK_AND_SCAN_ONLY",
        "lava_scientific_evidence": {
            "status": "FAILED_QC_NOT_CONSUMED",
            "scientific_validation_passed": False,
            "consumption": "ORDERING_TERMINALITY_ATTESTATION_ONLY",
            "passed_lava_results_used_by_placo": False,
            "canonical_lava_result_lock_present": False,
        },
        "lava_terminal": lava,
        "frozen_v1": {
            "policy": identity(ROOT / POLICY),
            "contract_lock": identity(ROOT / V1_CONTRACT_LOCK),
            "input_gate": identity(ROOT / V1_INPUT_GATE),
            "readiness_gate": identity(ROOT / V1_READINESS_GATE),
            "input_gate_lock": identity(ROOT / V1_INPUT_LOCK),
            "upstream": v1["upstream_identities"],
            "dense_sources": v1["dense_identities"],
            "placo_source": v1["placo_source"],
        },
        "execution_code": {
            "v1_contract": identity(ROOT / V1_CONTRACT_SCRIPT),
            "v1_input_gate": identity(ROOT / V1_INPUT_GATE_SCRIPT),
            "v1_materialization_validation_engine": identity(ROOT / V1_MATERIALIZER),
            "v1_full_p_runner": identity(ROOT / V1_RUNNER),
            "v2_worker_bound_full_p_runner": identity(ROOT / V2_RUNNER),
            "r_runtime_launcher": identity(ROOT / R_RUNTIME),
            "v1_ld_materializer_not_executed": identity(ROOT / V1_LD_MATERIALIZER),
            "v2_terminal_gate": identity(ROOT / GATE_SCRIPT),
            "v2_materializer_bridge": identity(ROOT / BRIDGE_SCRIPT),
            "lava_reference_archive_state_guard": identity(ROOT / ARCHIVE_STATE_GUARD_SCRIPT),
        },
        "resource_envelope": {
            "minimum_free_bytes_at_runner_start": 8 * 1024**3,
            "materialization_estimate_rule": "8_GIB_PLUS_600_BYTES_TIMES_SUM_OF_PAIR_SOURCE_ROWS",
            "maximum_workers": 4,
            "pair_concurrency_limit": 1,
            "shard_size": 20_000,
            "ram_safety_reserve_bytes": RAM_SAFETY_RESERVE_BYTES,
            "matched_worker_benchmark_required_before_full_scan": True,
            "multiworker_benchmark_requires_aggregate_process_tree_rss": True,
            "live_rss_guard_boundary": "FRESH_SESSION_PGID_INCLUDING_REPARENTED_DESCENDANTS",
            "worker_identity_bound_in_nuisance_shards_stage_and_publication": True,
            "inherited_worker_lease_authorization": "PAIR_RUN_FINGERPRINT_DEVICE_INODE_AND_COORDINATOR_CODE",
            "repo_scoped_single_pair_flock_required": True,
            "global_nuisance_rule": "OFFICIAL_VAR_PLACO_AND_COR_PEARSON_ON_ALL_VALID_ALIGNED_VARIANTS_BEFORE_SHARDING",
        },
        "reference_and_publication": {
            "full_p_scan_ld_reference_required": False,
            "ld_locus_publication_allowed": False,
            "required_to_unblock": "SEALED_FULL_EUR_HG19_LD_AND_VERSIONED_FAMILY_COLLATOR",
        },
        "readiness_gate": str(READINESS_GATE),
        "readiness_gate_sha256": hashlib.sha256(readiness.encode()).hexdigest(),
    }


def build_expected() -> tuple[str, dict[str, Any]]:
    v1 = validate_v1_pleiotropy()
    lava = deep_lava_terminal()
    readiness = readiness_text()
    return readiness, lock_payload(v1, lava, readiness)


def verify_gate() -> dict[str, Any]:
    readiness, expected = build_expected()
    readiness_path = ROOT / READINESS_GATE
    expected_readiness = readiness.encode()
    try:
        observed_size, observed_hash = file_identity(readiness_path)
    except GateError as error:
        raise GateError("additive PLACO V2 readiness gate is missing or drifted") from error
    if (
        observed_size != len(expected_readiness)
        or observed_hash != hashlib.sha256(expected_readiness).hexdigest()
    ):
        raise GateError("additive PLACO V2 readiness gate is missing or drifted")
    observed = read_json(ROOT / TERMINAL_GATE_LOCK)
    if observed != expected:
        raise GateError("additive PLACO V2 terminal-gate lock is missing or drifted")
    return observed


def _publish_family_locked(readiness: str, payload: dict[str, Any]) -> None:
    paths = [ROOT / READINESS_GATE, ROOT / TERMINAL_GATE_LOCK]
    encoded = [readiness.encode(), json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"]
    reconcile_gate_temporaries(encoded)
    for path in paths:
        require_real_repository_path(path.parent, "PLACO V2 gate publication parent")
    existing = [path.exists() or path.is_symlink() for path in paths]
    if existing == [False, True]:
        raise GateError("partial additive PLACO V2 gate is not an exact publication prefix")
    for present, path, content in zip(existing, paths, encoded, strict=True):
        if not present:
            continue
        observed_size, observed_hash = file_identity(path)
        if observed_size != len(content) or observed_hash != hashlib.sha256(content).hexdigest():
            raise GateError("partial additive PLACO V2 gate contains noncanonical bytes")
    if all(existing):
        verify_gate()
        return
    try:
        policy = PLEIOTROPY.validate_policy(ROOT)
        evidence = placo_result_evidence(policy)
    except SystemExit as error:
        raise GateError(str(error)) from error
    if evidence:
        raise GateError(
            "refusing to freeze PLACO V2 terminal gate after result access: "
            + ", ".join(str(path.relative_to(ROOT)) for path in evidence[:3])
        )
    assert_no_canonical_lava_results()
    # A SIGKILL can occur between the two no-replace links.  Exact bytes from
    # the already-linked prefix are retained; only missing members are created.
    # Any noncanonical or out-of-order state remains a hard failure.
    temporaries: list[tuple[Path, int, int, Path]] = []
    published: list[tuple[Path, int, int]] = []
    try:
        for present, path, content in zip(existing, paths, encoded, strict=True):
            if present:
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            require_real_repository_path(path.parent, "PLACO V2 gate publication parent")
            temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
            with temporary.open("xb") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            temporary_stat = os.lstat(temporary)
            if not stat.S_ISREG(temporary_stat.st_mode):
                raise GateError(f"PLACO V2 gate temporary is not a regular file: {temporary}")
            temporaries.append((temporary, temporary_stat.st_dev, temporary_stat.st_ino, path))
        # Recheck at publication time: building the deeply validated terminal
        # payload can be slow, and neither this continuation nor V1 owns a
        # shared publication mutex.
        fresh_readiness, fresh_payload = build_expected()
        if fresh_readiness != readiness or fresh_payload != payload:
            raise GateError("LAVA terminal or frozen V1 identity drifted while freezing the PLACO V2 gate")
        evidence = placo_result_evidence(policy)
        assert_no_canonical_lava_results()
        if evidence:
            raise GateError("scientific result evidence appeared while freezing the PLACO V2 gate")
        for temporary, _, _, path in temporaries:
            require_real_repository_path(path.parent, "PLACO V2 gate publication parent")
            source_stat = os.lstat(temporary)
            os.link(temporary, path)
            destination_stat = os.lstat(path)
            published.append((path, source_stat.st_dev, source_stat.st_ino))
            if (
                not stat.S_ISREG(destination_stat.st_mode)
                or (destination_stat.st_dev, destination_stat.st_ino)
                != (source_stat.st_dev, source_stat.st_ino)
            ):
                raise GateError(f"PLACO V2 no-replace publication identity mismatch: {path}")
            fsync_directory(path.parent)
        fresh_readiness, fresh_payload = build_expected()
        evidence = placo_result_evidence(policy)
        assert_no_canonical_lava_results()
        if fresh_readiness != readiness or fresh_payload != payload or evidence:
            raise GateError("scientific result evidence appeared during PLACO V2 gate publication")
        verify_gate()
    except BaseException:
        for path, device, inode in published:
            unlink_if_identity(path, (device, inode))
        raise
    finally:
        for path, device, inode, _ in temporaries:
            unlink_if_identity(path, (device, inode))


def publish_family(readiness: str, payload: dict[str, Any]) -> None:
    with gate_publication_lock():
        _publish_family_locked(readiness, payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--seal", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        payload = verify_gate()
        print(
            "TRACK_B_PLACO_TERMINAL_GATE_V2_VERIFIED "
            f"continuation={payload['lava_terminal']['continuation_execution_fingerprint']} "
            "lava_evidence=FAILED_QC_NOT_CONSUMED pairs=A,B,CONTROL"
        )
    else:
        readiness, payload = build_expected()
        publish_family(readiness, payload)
        print(
            "TRACK_B_PLACO_TERMINAL_GATE_V2_SEALED "
            f"continuation={payload['lava_terminal']['continuation_execution_fingerprint']} "
            "lava_evidence=FAILED_QC_NOT_CONSUMED pairs=A,B,CONTROL"
        )


if __name__ == "__main__":
    try:
        main()
    except GateError as error:
        raise SystemExit(f"ERROR: {error}") from error
