#!/usr/bin/env python3
"""Read-only, full-family integrity and frozen failure-gate audit for FI×sleep LAVA."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SERIAL_PATH = ROOT / "frailty_paper/scripts/45_run_lava_sensitivity.py"
INTERIM_AUDIT_PATH = ROOT / "frailty_paper/scripts/47_audit_lava_interim_receipts.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import project module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SERIAL = _load(SERIAL_PATH, "frailty_lava_serial_for_full_audit")
INTERIM = _load(INTERIM_AUDIT_PATH, "frailty_lava_interim_audit_for_full_audit")
TRAITS = tuple(SERIAL.TRAITS)
LOCUS_LIMIT = 2495
EXPECTED_RECEIPTS = len(TRAITS) * LOCUS_LIMIT
FAILURE_LIMIT = 0.01
RECEIPT_RE = re.compile(r"^locus_(\d{4})\.json$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _state(run_root: Path) -> dict[str, Any] | None:
    path = run_root / "runner_state.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _claim_snapshot(run_root: Path) -> tuple[list[dict[str, Any]], list[str]]:
    claims: list[dict[str, Any]] = []
    issues: list[str] = []
    claim_root = run_root / "claims"
    if not claim_root.exists():
        return claims, issues
    for path in sorted(claim_root.glob("*/*.claim")):
        if path.name.startswith("._"):
            continue
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
            identity = (value.get("trait"), value.get("locus_index"), value.get("pair_order"))
            if not all(isinstance(x, (str, int)) for x in identity):
                raise ValueError("claim identity fields are missing")
            claims.append({"path": path.relative_to(run_root).as_posix(), **value})
        except (OSError, ValueError, json.JSONDecodeError) as error:
            issues.append(f"{path.relative_to(run_root).as_posix()}: {error}")
    return claims, issues


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def audit_run(run_root: Path) -> dict[str, Any]:
    """Audit a fixed file-list snapshot while the append-only runner may continue."""
    run_root = run_root.resolve()
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    state_start = _state(run_root)
    receipt_list = sorted(
        ((trait_index, trait, path)
         for trait_index, trait in enumerate(TRAITS, start=1)
         for path in (run_root / "pairs" / trait / "loci").glob("locus_*.json")),
        key=lambda row: (row[0], row[2].name),
    )

    issues: list[str] = []
    duplicate_identities: list[dict[str, Any]] = []
    seen: dict[tuple[str, int], str] = {}
    by_trait: dict[str, Counter[str]] = {trait: Counter() for trait in TRAITS}
    failure_categories: Counter[str] = Counter()
    failure_counts: Counter[str] = Counter()
    digest = hashlib.sha256()
    for trait_order, trait, receipt_path in receipt_list:
        match = RECEIPT_RE.fullmatch(receipt_path.name)
        if match is None:
            issues.append(f"invalid receipt filename: {receipt_path}")
            continue
        locus = int(match.group(1))
        if not 1 <= locus <= LOCUS_LIMIT:
            issues.append(f"locus index outside frozen family: {receipt_path}")
            continue
        log_path = run_root / "pairs" / trait / "logs" / f"locus_{locus:04d}.log"
        identity = (trait, locus)
        rel = receipt_path.relative_to(run_root).as_posix()
        if identity in seen:
            duplicate_identities.append({"trait": trait, "locus_index": locus,
                                         "paths": [seen[identity], rel]})
            issues.append(f"duplicate receipt identity: {trait}/{locus}")
            continue
        seen[identity] = rel
        try:
            value = SERIAL.verify_receipt(receipt_path, trait, trait_order, locus)
        except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
            issues.append(f"invalid receipt {rel}: {error}")
            continue
        if not log_path.is_file() or log_path.stat().st_size == 0:
            issues.append(f"missing or empty log: {log_path.relative_to(run_root).as_posix()}")
            continue
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
        flagged, counts, unexpected = INTERIM.classify_receipt(value, log_text, trait)
        by_trait[trait]["receipt_files_validated"] += 1
        by_trait[trait][value["process_status"].lower()] += 1
        if value.get("fi_status") != "TESTED":
            by_trait[trait]["fi_status_not_tested"] += 1
        if value.get("sleep_status") != "TESTED":
            by_trait[trait]["sleep_status_not_tested"] += 1
        for key, count in counts.items():
            # process_status was already counted once from the receipt field.
            if key != value["process_status"].lower():
                by_trait[trait][key] += count
        if flagged:
            failure_counts[trait] += 1
        if value["process_status"] == "PROCESS_FAILED":
            if "none of specified SNP IDs are present in reference data" in log_text:
                failure_categories["no_specified_snps_in_reference"] += 1
            elif "Negative variance estimate for all phenotypes" in log_text:
                failure_categories["all_phenotypes_negative_variance"] += 1
            else:
                failure_categories["other_process_failure"] += 1
        # Same ordered, path-bound receipt hash stream used by the earlier full audit.
        digest.update(rel.encode("utf-8") + b"\0")
        digest.update(receipt_path.read_bytes())
        digest.update(b"\0")

    state_end = _state(run_root)
    claims, claim_issues = _claim_snapshot(run_root)
    issues.extend(claim_issues)
    duplicated_claims: list[dict[str, Any]] = []
    claim_seen: dict[tuple[str, int], str] = {}
    for claim in claims:
        identity = (claim["trait"], int(claim["locus_index"]))
        if identity in claim_seen:
            duplicated_claims.append({"trait": identity[0], "locus_index": identity[1],
                                      "paths": [claim_seen[identity], claim["path"]]})
        else:
            claim_seen[identity] = claim["path"]
    active_jobs = (state_end or {}).get("active_jobs", [])
    active_workers = sorted({int(job["worker_pid"]) for job in active_jobs
                             if isinstance(job.get("worker_pid"), int)})
    runner_receipts = (state_end or {}).get("verified_receipts")
    actual_after = sum(1 for trait in TRAITS
                       for _ in (run_root / "pairs" / trait / "loci").glob("locus_*.json"))
    trait_summary = {}
    for trait in TRAITS:
        failed = failure_counts[trait]
        trait_summary[trait] = {
            **dict(sorted(by_trait[trait].items())),
            "expected_loci": LOCUS_LIMIT,
            "process_or_univariate_failure_lower_bound": failed,
            "failure_fraction": failed / LOCUS_LIMIT,
            "passes_1pct_gate": failed / LOCUS_LIMIT <= FAILURE_LIMIT,
        }
    inputs_manifest = run_root / "input_manifest.json"
    return {
        "schema_version": 1,
        "audit_time_start_utc": started,
        "audit_time_end_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "analysis_id": (state_end or state_start or {}).get("analysis_id"),
        "expected_receipts": EXPECTED_RECEIPTS,
        "receipt_files_at_scan_start": len(receipt_list),
        "receipt_files_after_scan": actual_after,
        "receipt_files_validated": sum(v["receipt_files_validated"] for v in by_trait.values()),
        "runner_verified_receipts_at_start": (state_start or {}).get("verified_receipts"),
        "runner_verified_receipts_at_end": runner_receipts,
        "receipt_file_list_changed_during_scan": actual_after != len(receipt_list),
        "receipt_hash_stream_sha256": digest.hexdigest(),
        "receipt_issues": issues,
        "duplicate_trait_locus_identities": duplicate_identities,
        "duplicate_claim_identities": duplicated_claims,
        "claims": claims,
        "claim_issues": claim_issues,
        "active_job_pids": active_workers,
        "active_job_pid_liveness": {str(pid): _pid_alive(pid) for pid in active_workers},
        "effective_parallel_workers": (state_end or {}).get("effective_parallel_workers"),
        "requested_parallel_workers": (state_end or {}).get("requested_parallel_workers"),
        "launch_failures_this_invocation": (state_end or {}).get("launch_failures_this_invocation"),
        "stale_claims_recovered": (state_end or {}).get("stale_claims_recovered"),
        "runner_state_updated_utc": (state_end or {}).get("updated_utc"),
        "runner_state_lock_sha256": (state_end or {}).get("lock_sha256"),
        "resource_fallback_reason": (state_end or {}).get("resource_fallback_reason"),
        "process_failure_categories": dict(sorted(failure_categories.items())),
        "process_failures_by_trait": trait_summary,
        "locked_max_locus_failure_fraction": FAILURE_LIMIT,
        "all_traits_pass_1pct_gate": all(row["passes_1pct_gate"] for row in trait_summary.values()),
        "lock_sha256": sha256(SERIAL.LOCK),
        "input_manifest_sha256": sha256(inputs_manifest) if inputs_manifest.is_file() else None,
        "serial_runner_sha256": sha256(SERIAL_PATH),
        "interim_classifier_sha256": sha256(INTERIM_AUDIT_PATH),
        "auditor_sha256": sha256(Path(__file__)),
        "external_run_root": str(run_root),
        "mutates_run_or_receipts": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=SERIAL.OUTPUT_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {output}")
    result = audit_run(args.run_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "audit_time_start_utc", "audit_time_end_utc", "receipt_files_at_scan_start",
        "receipt_files_after_scan", "receipt_files_validated", "receipt_issues",
        "effective_parallel_workers", "requested_parallel_workers",
        "process_failure_categories", "all_traits_pass_1pct_gate")}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
