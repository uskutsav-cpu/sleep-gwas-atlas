#!/usr/bin/env python3
"""Freeze/verify revision 2 of the result-blind Track B mechanism contract.

Revision 1 remains byte-exact for audit.  Revision 2 applies one narrowly
specified correction: cell-sQTL payloads follow the already-frozen molecular
policy's leafcutter quantification family, not exon quantification.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import io
import json
import os
import stat
import sys
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

sys.dont_write_bytecode = True

ROOT_AT_IMPORT = Path(__file__).resolve().parents[1]
BASE_REL = "scripts/160_build_track_b_mechanism_readiness.py"
BASE_SHA256 = "7283de7830ed6a9602f32b3ebc9e38768efea97ccbb4b48d543dee4aa9629b27"


def _bootstrap_hash(path: Path) -> str:
    info = os.lstat(path)
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise RuntimeError(f"base implementation is not a real regular file: {path}")
    digest = hashlib.sha256()
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(path, flags)
    try:
        before = os.fstat(fd)
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise RuntimeError("base implementation drifted during bootstrap")
    return digest.hexdigest()


if _bootstrap_hash(ROOT_AT_IMPORT / BASE_REL) != BASE_SHA256:
    raise RuntimeError("revision-1 base implementation differs from the revision-2 pin")
_BASE_SPEC = importlib.util.spec_from_file_location("track_b_mechanism_readiness_v1_base", ROOT_AT_IMPORT / BASE_REL)
if _BASE_SPEC is None or _BASE_SPEC.loader is None:
    raise RuntimeError("cannot load revision-1 base implementation")
B = importlib.util.module_from_spec(_BASE_SPEC)
sys.modules[_BASE_SPEC.name] = B
_BASE_SPEC.loader.exec_module(B)

SCHEMA = "track-b-mechanism-pre-result-readiness.2"
POLICY_REL = "config/track_b_mechanism_followup_policy_v2.json"
SCRIPT_REL = "scripts/161_build_track_b_mechanism_readiness_v2.py"
OUTPUT_DIR_REL = "results/track_b/mechanism_followup/v2"
OUTPUT_RELS = (
    f"{OUTPUT_DIR_REL}/PRE_RESULT_RESOURCE_EVIDENCE.tsv",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.tsv",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.json",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS_REPORT.md",
)
SUPERSESSION_REL = "results/track_b/mechanism_followup/PRE_RESULT_READINESS_SUPERSESSION.json"
LOCK_REL = f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.lock.json"
OLD_OUTPUT_RELS = B.OUTPUT_RELS
OLD_LOCK_REL = B.LOCK_REL
EXPECTED_OVERLAY_KEYS = {
    "schema_version", "contract_revision", "base_policy_path", "base_policy_sha256",
    "base_implementation_path", "base_implementation_sha256", "supersedes",
    "exact_correction", "forbidden_change", "result_blind", "future_results_accessed",
}


def load_effective_policy(root: Path) -> tuple[dict[str, object], dict[str, object], str]:
    overlay = B.read_json(root, POLICY_REL)
    if set(overlay) != EXPECTED_OVERLAY_KEYS or overlay.get("schema_version") != SCHEMA or overlay.get("contract_revision") != 2:
        raise B.ContractError("revision-2 overlay schema is invalid")
    if overlay.get("result_blind") is not True or overlay.get("future_results_accessed") is not False:
        raise B.ContractError("revision-2 overlay is not result-blind")
    if overlay.get("base_policy_path") != B.POLICY_REL or overlay.get("base_policy_sha256") != "503d34e1cdebac5c1e427c41a71b9c320bde948828ddf18df9767d7ac4e06c7e":
        raise B.ContractError("revision-2 base-policy lineage differs from its exact pin")
    if overlay.get("base_implementation_path") != BASE_REL or overlay.get("base_implementation_sha256") != BASE_SHA256:
        raise B.ContractError("revision-2 base-implementation lineage differs from its exact pin")
    if B.stable_file(root, B.POLICY_REL).sha256 != overlay["base_policy_sha256"] or B.stable_file(root, BASE_REL).sha256 != BASE_SHA256:
        raise B.ContractError("revision-1 policy/implementation drifted")
    base = B.read_json(root, B.POLICY_REL)
    B.validate_policy_shape(base)
    correction = overlay.get("exact_correction")
    if not isinstance(correction, dict) or set(correction) != {"required_qtl_datasets", "analysis_id", "scope", "required_inputs"}:
        raise B.ContractError("revision-2 exact correction schema is invalid")
    expected_qtl = {"eQTL": ["QTD000559", "QTD000569"], "sQTL": ["QTD000563", "QTD000573"]}
    if correction.get("required_qtl_datasets") != expected_qtl or correction.get("analysis_id") != "P18_CELL_SQTL":
        raise B.ContractError("revision-2 QTL correction is not the exact leafcutter-only change")
    if base.get("required_qtl_datasets") != {
        "eQTL": ["QTD000559", "QTD000569"],
        "sQTL": ["QTD000560", "QTD000563", "QTD000570", "QTD000573"],
    }:
        raise B.ContractError("revision-1 QTL family is not the known superseded family")
    effective = copy.deepcopy(base)
    effective["required_qtl_datasets"] = copy.deepcopy(expected_qtl)
    p18 = [row for row in effective["analyses"] if row["analysis_id"] == "P18_CELL_SQTL"]
    if len(p18) != 1:
        raise B.ContractError("revision-1 P18 row is not unique")
    p18[0]["scope"] = correction["scope"]
    p18[0]["required_inputs"] = correction["required_inputs"]
    B.validate_policy_shape(effective)
    return overlay, effective, B.bytes_sha256(B.canonical_json(effective))


def validate_superseded_family(root: Path, overlay: Mapping[str, object]) -> dict[str, object]:
    supersedes = overlay.get("supersedes")
    if not isinstance(supersedes, dict) or set(supersedes) != {"lock_path", "lock_sha256", "retention", "reason"}:
        raise B.ContractError("supersession lineage record is malformed")
    if supersedes.get("lock_path") != OLD_LOCK_REL or supersedes.get("lock_sha256") != "4ab3c004d109fa856936a5411f4495975b176333916bbe3f70b7974aac9de40f":
        raise B.ContractError("superseded lock identity differs from its exact pin")
    if supersedes.get("retention") != "PRESERVE_BYTE_EXACT_AS_AUDIT_TRAIL_DO_NOT_USE_FOR_EXECUTION":
        raise B.ContractError("superseded freeze is not explicitly retained")
    lock_payload = B.read_stable_bytes(root, OLD_LOCK_REL)
    if hashlib.sha256(lock_payload).hexdigest() != supersedes["lock_sha256"]:
        raise B.ContractError("superseded lock bytes differ from revision-2 lineage")
    try:
        lock = json.loads(lock_payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise B.ContractError(f"superseded lock is invalid JSON: {exc}") from exc
    if B.canonical_json(lock) != lock_payload or lock.get("script_sha256") != BASE_SHA256 or lock.get("policy_sha256") != overlay["base_policy_sha256"]:
        raise B.ContractError("superseded lock canonical policy/script lineage is invalid")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != set(OLD_OUTPUT_RELS):
        raise B.ContractError("superseded output family is incomplete")
    for relative in OLD_OUTPUT_RELS:
        item = B.stable_file(root, relative)
        expected = outputs[relative]
        if item.observed_bytes != expected.get("bytes") or item.sha256 != expected.get("sha256"):
            raise B.ContractError(f"superseded output is not byte-exact: {relative}")
    return lock


def _tsv_bytes(fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> bytes:
    return B.tsv_bytes(fields, rows)


def supersession_receipt(overlay: Mapping[str, object], effective_hash: str) -> dict[str, object]:
    supersedes = overlay["supersedes"]
    assert isinstance(supersedes, dict)
    return {
        "schema_version": SCHEMA,
        "record_kind": "IMMUTABLE_SUPERSESSION_RECEIPT",
        "superseded_contract": {
            "lock_path": supersedes["lock_path"],
            "lock_sha256": supersedes["lock_sha256"],
            "status": "SUPERSEDED_DO_NOT_USE_FOR_EXECUTION",
            "retention": supersedes["retention"],
        },
        "replacement_contract": {
            "policy_path": POLICY_REL,
            "policy_sha256": B.stable_file(ROOT_AT_IMPORT, POLICY_REL).sha256,
            "effective_policy_sha256": effective_hash,
            "script_path": SCRIPT_REL,
            "script_sha256": B.stable_file(ROOT_AT_IMPORT, SCRIPT_REL).sha256,
            "lock_path": LOCK_REL,
        },
        "scientific_correction": supersedes["reason"],
        "result_blind": True,
        "future_results_accessed": False,
    }


def render_report_v2(
    policy_hash: str,
    evidence_rows: Sequence[Mapping[str, object]],
    readiness: Sequence[Mapping[str, object]],
    supersedes: Mapping[str, object],
) -> bytes:
    base = B.render_report(policy_hash, evidence_rows, readiness).decode("utf-8")
    notice = (
        "# Track B mechanism follow-up: pre-result local readiness — revision 2\n\n"
        f"This contract supersedes `{supersedes['lock_path']}` (SHA-256 `{supersedes['lock_sha256']}`). "
        "Revision 1 is retained byte-exact as an audit trail and must not be used for execution. "
        "Revision 2 corrects the cell-sQTL family to the predeclared leafcutter datasets QTD000563 and QTD000573 only.\n\n"
    )
    return (notice + base.replace("# Track B mechanism follow-up: pre-result local readiness\n\n", "", 1)).encode("utf-8")


def construct_artifacts(root: Path) -> tuple[dict[str, bytes], dict[str, object]]:
    overlay, effective, effective_hash = load_effective_policy(root)
    old_lock = validate_superseded_family(root, overlay)
    policy_hash = B.stable_file(root, POLICY_REL).sha256
    script_hash = B.stable_file(root, SCRIPT_REL).sha256
    registry, manifests = B.collect_evidence(root, effective)
    registry.add(POLICY_REL, expected_sha256=policy_hash, detail="revision-2 correction policy")
    registry.add(SCRIPT_REL, expected_sha256=script_hash, detail="revision-2 contract implementation")
    registry.add(OLD_LOCK_REL, expected_sha256=str(overlay["supersedes"]["lock_sha256"]), detail="retained superseded lock")
    for relative, expected in sorted(old_lock["outputs"].items()):
        registry.add(relative, expected_bytes=int(expected["bytes"]), expected_sha256=str(expected["sha256"]), detail="retained byte-exact superseded output")
    facts = B.build_facts(root, effective, registry, manifests)
    registry.recheck()
    evidence_rows = registry.rows()
    readiness = B.evaluate_analyses(effective, facts, registry)
    p18 = next(row for row in readiness if row["analysis_id"] == "P18_CELL_SQTL")
    if "QTD000560" in p18["evidence_paths"] or "QTD000570" in p18["evidence_paths"] or "QTD000563" not in p18["evidence_paths"] or "QTD000573" not in p18["evidence_paths"]:
        raise B.ContractError("revision-2 P18 evidence is not exactly leafcutter-only")
    evidence_digest = B.bytes_sha256(B.canonical_json(evidence_rows))
    resource_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    for row in readiness:
        resource_counts[str(row["resource_status"])] = resource_counts.get(str(row["resource_status"]), 0) + 1
        status_counts[str(row["status"])] = status_counts.get(str(row["status"]), 0) + 1
    document = {
        "schema_version": SCHEMA,
        "contract_revision": 2,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS",
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": POLICY_REL,
        "policy_sha256": policy_hash,
        "effective_policy_sha256": effective_hash,
        "script_path": SCRIPT_REL,
        "script_sha256": script_hash,
        "supersedes": overlay["supersedes"],
        "analysis_family_count": len(readiness),
        "evidence_record_count": len(evidence_rows),
        "evidence_bundle_sha256": evidence_digest,
        "resource_status_counts": dict(sorted(resource_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "facts": [facts[name].record() for name in sorted(facts)],
        "analyses": readiness,
    }
    receipt = {
        "schema_version": SCHEMA,
        "record_kind": "IMMUTABLE_SUPERSESSION_RECEIPT",
        "superseded_contract": {
            "lock_path": overlay["supersedes"]["lock_path"],
            "lock_sha256": overlay["supersedes"]["lock_sha256"],
            "status": "SUPERSEDED_DO_NOT_USE_FOR_EXECUTION",
            "retention": overlay["supersedes"]["retention"],
        },
        "replacement_contract": {
            "policy_path": POLICY_REL,
            "policy_sha256": policy_hash,
            "effective_policy_sha256": effective_hash,
            "script_path": SCRIPT_REL,
            "script_sha256": script_hash,
            "lock_path": LOCK_REL,
        },
        "scientific_correction": overlay["supersedes"]["reason"],
        "result_blind": True,
        "future_results_accessed": False,
    }
    artifacts = {
        OUTPUT_RELS[0]: _tsv_bytes(B.EVIDENCE_FIELDS, evidence_rows),
        OUTPUT_RELS[1]: _tsv_bytes(B.READINESS_FIELDS, readiness),
        OUTPUT_RELS[2]: B.canonical_json(document),
        OUTPUT_RELS[3]: render_report_v2(policy_hash, evidence_rows, readiness, overlay["supersedes"]),
        SUPERSESSION_REL: B.canonical_json(receipt),
    }
    lock = {
        "schema_version": SCHEMA,
        "contract_revision": 2,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK",
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": POLICY_REL,
        "policy_sha256": policy_hash,
        "effective_policy_sha256": effective_hash,
        "script_sha256": script_hash,
        "superseded_lock_path": OLD_LOCK_REL,
        "superseded_lock_sha256": overlay["supersedes"]["lock_sha256"],
        "analysis_family_count": len(readiness),
        "evidence_record_count": len(evidence_rows),
        "evidence_bundle_sha256": evidence_digest,
        "outputs": {
            relative: {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
            for relative, payload in sorted(artifacts.items())
        },
    }
    artifacts[LOCK_REL] = B.canonical_json(lock)
    return artifacts, lock


def freeze_no_replace(root: Path, artifacts: Mapping[str, bytes]) -> None:
    output_dir = B._secure_mkdirs(root, OUTPUT_DIR_REL)
    ordered = list(OUTPUT_RELS) + [SUPERSESSION_REL, LOCK_REL]
    if set(artifacts) != set(ordered):
        raise B.ContractError("revision-2 artifact family is incomplete")
    for relative in ordered:
        try:
            os.lstat(root / relative)
        except FileNotFoundError:
            continue
        raise B.ContractError(f"no-replace revision-2 freeze refused existing target: {relative}")
    stage = Path(tempfile.mkdtemp(prefix=".readiness-v2-freeze-", dir=output_dir))
    staged: dict[str, Path] = {}
    try:
        for index, relative in enumerate(ordered):
            path = stage / f"{index:02d}.stage"
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
            try:
                view = memoryview(artifacts[relative])
                while view:
                    written = os.write(fd, view)
                    if written <= 0:
                        raise B.ContractError(f"short write while staging {relative}")
                    view = view[written:]
                os.fsync(fd)
            finally:
                os.close(fd)
            staged[relative] = path
        for relative in ordered:
            parent = (root / relative).parent.relative_to(root).as_posix()
            B._secure_mkdirs(root, parent)
            os.link(staged[relative], root / relative, follow_symlinks=False)
            item = B.stable_file(root, relative)
            if item.observed_bytes != len(artifacts[relative]) or item.sha256 != hashlib.sha256(artifacts[relative]).hexdigest():
                raise B.ContractError(f"published revision-2 artifact mismatch: {relative}")
        dir_fd = os.open(output_dir, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    finally:
        for path in staged.values():
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
        try:
            os.rmdir(stage)
        except FileNotFoundError:
            pass


def verify_frozen(root: Path) -> None:
    lock_payload = B.read_stable_bytes(root, LOCK_REL)
    try:
        lock = json.loads(lock_payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise B.ContractError(f"invalid revision-2 lock: {exc}") from exc
    if not isinstance(lock, dict) or B.canonical_json(lock) != lock_payload or lock.get("schema_version") != SCHEMA:
        raise B.ContractError("revision-2 lock is noncanonical or has the wrong schema")
    expected_paths = set(OUTPUT_RELS) | {SUPERSESSION_REL}
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != expected_paths:
        raise B.ContractError("revision-2 lock output family is incomplete")
    for relative in sorted(expected_paths):
        item = B.stable_file(root, relative)
        expected = outputs[relative]
        if item.observed_bytes != expected.get("bytes") or item.sha256 != expected.get("sha256"):
            raise B.ContractError(f"revision-2 output differs from lock: {relative}")
    recomputed, expected_lock = construct_artifacts(root)
    for relative in sorted(expected_paths):
        if B.read_stable_bytes(root, relative) != recomputed[relative]:
            raise B.ContractError(f"current resources no longer reproduce revision-2 output: {relative}")
    if lock != expected_lock or lock_payload != recomputed[LOCK_REL]:
        raise B.ContractError("current resources no longer reproduce revision-2 lock")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument("--root", type=Path, default=ROOT_AT_IMPORT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if root != ROOT_AT_IMPORT and _bootstrap_hash(root / BASE_REL) != BASE_SHA256:
        raise B.ContractError("selected repository root lacks the exact revision-1 base")
    if args.freeze:
        artifacts, _lock = construct_artifacts(root)
        freeze_no_replace(root, artifacts)
        print(f"FROZEN_REVISION_2 {len(B.EXPECTED_ANALYSIS_IDS)} analyses at {OUTPUT_DIR_REL}")
    else:
        verify_frozen(root)
        print(f"VERIFIED_REVISION_2 {len(B.EXPECTED_ANALYSIS_IDS)} analyses at {OUTPUT_DIR_REL}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except B.ContractError as exc:
        raise SystemExit(f"ERROR: {exc}")
