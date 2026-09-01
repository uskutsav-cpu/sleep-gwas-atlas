#!/usr/bin/env python3
"""Verify and seal the dedicated Track B LAVA execution contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import lava_contract


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "config/track_b_local_analysis_policy.json"
INPUT_LOCK = ROOT / "results/track_b/local_analysis_input.lock.json"
REFERENCE = ROOT / "ref/lava/ukb_v1.1/reference.provenance.json"
CHECKPOINT_DIR = ROOT / "results/track_b/checkpoints/lava"
RESULT_DIR = ROOT / "results/track_b/local"
RESULTS = [
    RESULT_DIR / "lava_locus_status.tsv", RESULT_DIR / "lava_univariate.tsv",
    RESULT_DIR / "lava_bivariate.tsv", RESULT_DIR / "lava_conditional.tsv",
]
RESULT_LOCK = RESULT_DIR / "lava_results.provenance.json"


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path.relative_to(ROOT)}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_policy() -> dict[str, object]:
    try:
        policy = json.loads(POLICY.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B local policy is unreadable: {error}")
    if (
        policy.get("analysis_id") != "track-b-v1.0-local"
        or policy.get("lava_version") != "0.1.5"
        or policy.get("expected_loci") != 2495
        or policy.get("expected_analysis_traits") != 8
        or policy.get("expected_discovery_pairs") != 3
        or policy.get("planned_univariate_tests") != 2495 * 8
        or policy.get("planned_bivariate_pair_locus_family_max") != 2495 * 3
    ):
        fail("Track B local policy scope drifted")
    return policy


def validate_inputs(policy: dict[str, object]) -> dict[str, object]:
    try:
        lock = json.loads(INPUT_LOCK.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B local input lock is unreadable: {error}")
    if (
        lock.get("analysis_id") != policy["analysis_id"]
        or lock.get("policy_sha256") != sha256(POLICY)
        or lock.get("local_results_accessed_before_input_freeze") is not False
    ):
        fail("Track B local input lock differs from policy or pre-result state")
    hashes = lock.get("input_artifact_sha256", {})
    if not isinstance(hashes, dict) or len(hashes) != 7:
        fail("Track B local input artifact family is incomplete")
    for relative, expected in hashes.items():
        if sha256(ROOT / relative) != expected:
            fail(f"Track B local input differs from lock: {relative}")
    return lock


def validate_reference() -> dict[str, object]:
    generic_path, generic_policy = lava_contract.load_policy(ROOT)
    del generic_path
    return lava_contract.validate_reference(ROOT, generic_policy, rehash=False)


def run_fingerprint() -> str:
    policy = load_policy()
    validate_inputs(policy)
    validate_reference()
    payload = {
        "analysis_id": policy["analysis_id"], "policy_sha256": sha256(POLICY),
        "input_lock_sha256": sha256(INPUT_LOCK), "reference_sha256": sha256(REFERENCE),
        "runner_sha256": sha256(ROOT / "scripts/120_run_track_b_lava.R"),
        "validator_sha256": sha256(ROOT / "scripts/121_validate_track_b_lava.py"),
        "contract_sha256": sha256(Path(__file__)),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def result_payload() -> dict[str, object]:
    policy = load_policy()
    fingerprint = run_fingerprint()
    checkpoints = []
    for index in range(1, int(policy["expected_loci"]) + 1):
        path = CHECKPOINT_DIR / f"locus_{index:04d}.rds"
        checkpoints.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    results = [
        {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in RESULTS
    ]
    return {
        "schema_version": 1, "analysis_id": policy["analysis_id"],
        "lava_version": policy["lava_version"], "run_fingerprint": fingerprint,
        "input_lock_sha256": sha256(INPUT_LOCK), "reference_sha256": sha256(REFERENCE),
        "checkpoint_count": len(checkpoints), "checkpoints": checkpoints,
        "results": results,
    }


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def seal_results() -> None:
    expected = result_payload()
    if RESULT_LOCK.is_file():
        validate_results()
        return
    atomic_json(RESULT_LOCK, expected)
    print(f"TRACK_B_LAVA_RESULTS_SEALED checkpoints={expected['checkpoint_count']}")


def validate_results() -> None:
    try:
        observed = json.loads(RESULT_LOCK.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B LAVA result provenance is unreadable: {error}")
    if observed != result_payload():
        fail("Track B LAVA results differ from immutable provenance")
    print("TRACK_B_LAVA_RESULTS_VALIDATED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--preflight", action="store_true")
    group.add_argument("--run-fingerprint", action="store_true")
    group.add_argument("--seal-results", action="store_true")
    group.add_argument("--verify-results", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        fingerprint = run_fingerprint()
        print(f"TRACK_B_LAVA_PREFLIGHT_PASS fingerprint={fingerprint}")
    elif args.run_fingerprint:
        print(run_fingerprint())
    elif args.seal_results:
        seal_results()
    else:
        validate_results()


if __name__ == "__main__":
    main()
