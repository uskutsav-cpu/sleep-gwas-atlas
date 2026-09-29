#!/usr/bin/env python3
"""Read-only integrity and QC lower-bound audit for a partial canonical LAVA v3 run.

This tool never changes the production run. It verifies every observed receipt,
counts statuses over the exact frozen 7 × 2,495 family, and reports whether the
frozen untested-cell ceiling is already impossible to meet.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
import run_lava_canonical_v3 as canonical

DEFAULT_FAMILY = ROOT / "brain6/config/lava_family_canonical_v3.json"
DEFAULT_EXECUTION = ROOT / "brain6/config/lava_execution_canonical_v3.json"


def summarize(valid_loci: int, planned_loci: int, traits: tuple[str, ...],
              counts_by_trait: dict[str, collections.Counter[str]],
              maximum_untested_fraction: float) -> dict[str, Any]:
    """Summarize valid cells using the locked full-family denominator."""
    if planned_loci < 1 or not 0 <= valid_loci <= planned_loci or not traits:
        raise ValueError("Invalid canonical family dimensions or valid-locus count")
    if set(counts_by_trait) != set(traits):
        raise ValueError("Status summary must contain every canonical trait exactly once")
    intended = planned_loci * len(traits)
    expected_observed_per_trait = valid_loci
    for trait in traits:
        observed = sum(counts_by_trait[trait].values())
        if observed != expected_observed_per_trait:
            raise ValueError(f"Status counts for {trait} do not match verified receipts")

    totals: collections.Counter[str] = collections.Counter()
    by_trait: dict[str, dict[str, int]] = {}
    for trait in traits:
        counts = counts_by_trait[trait]
        if set(counts) - canonical.VALID_STATUSES:
            raise ValueError(f"Unknown canonical status for {trait}")
        totals.update(counts)
        by_trait[trait] = {status: counts.get(status, 0)
                           for status in sorted(canonical.VALID_STATUSES)}

    missing_cells = intended - valid_loci * len(traits)
    observed_untested = totals.get("NOT_RUN", 0) + totals.get("FAILED", 0)
    current_uncompleted = missing_cells + observed_untested
    max_fraction = maximum_untested_fraction
    if not 0 <= max_fraction <= 1:
        raise ValueError("Maximum untested fraction must lie in [0, 1]")
    max_allowed = int(intended * max_fraction)
    complete = valid_loci == planned_loci
    return {
        "state": "COMPLETE" if complete else "PARTIAL",
        "planned_loci": planned_loci,
        "verified_loci": valid_loci,
        "planned_cells": intended,
        "verified_cells": valid_loci * len(traits),
        "missing_or_invalid_cells": missing_cells,
        "currently_uncompleted_or_untested_cells": current_uncompleted,
        "statuses": {status: totals.get(status, 0)
                     for status in sorted(canonical.VALID_STATUSES)},
        "status_by_trait": by_trait,
        "untested_lower_bound_cells": observed_untested,
        "untested_lower_bound_fraction": observed_untested / intended,
        "maximum_untested_fraction": max_fraction,
        "maximum_allowed_untested_cells": max_allowed,
        "family_qc_guaranteed_fail": observed_untested > max_allowed,
        "family_qc_pass": complete and current_uncompleted <= max_allowed,
    }


def audit(run_dir: Path, input_root: Path, reference_provenance: Path,
          family_path: Path, execution_path: Path) -> dict[str, Any]:
    family = canonical.read_json(family_path)
    execution = canonical.read_json(execution_path)
    family_sha, execution_sha = canonical.sha256(family_path), canonical.sha256(execution_path)
    if family.get("analysis_id") != "brain6-lava-canonical-v3":
        raise ValueError("Unexpected canonical family identity")
    if execution.get("family_lock_sha256") != family_sha:
        raise ValueError("Execution lock does not bind the selected family lock")
    if (execution.get("worker_count") != 4 or
            execution.get("execution_policy", {}).get("worker_count") != 4):
        raise ValueError("Canonical execution lock does not preserve the four-worker policy")
    loci_path = ROOT / family["locus_definition"]["path"]
    loci = canonical.read_loci(loci_path, family["locus_definition"]["sha256"], 2495)
    traits = tuple(family["trait_ids"])
    if traits != canonical.TRAITS:
        raise ValueError("Canonical trait order differs from the frozen seven-trait family")

    identity_path = run_dir / "run_identity.json"
    identity = canonical.read_json(identity_path)
    run_id = run_dir.name
    if identity.get("run_id") != run_id:
        raise ValueError("Run identity document does not match the run directory")
    if identity.get("analysis_id") != family["analysis_id"]:
        raise ValueError("Run identity has a different analysis id")
    if identity.get("family_lock_sha256") != family_sha or identity.get("execution_lock_sha256") != execution_sha:
        raise ValueError("Run identity is not bound to the selected frozen locks")
    script_sha = canonical.sha256(canonical.SCRIPT)
    if identity.get("script_sha256") != script_sha:
        raise ValueError("Run identity is not bound to the current canonical worker script")
    if identity.get("run_identity_sha256") != run_id:
        raise ValueError("Run identity self-hash differs from the run directory")
    if canonical.read_json(run_dir / "family_lock.json") != family:
        raise ValueError("Run-local family lock differs from the frozen family lock")
    if canonical.read_json(run_dir / "execution_lock.json") != execution:
        raise ValueError("Run-local execution lock differs from the frozen execution lock")
    manifest_path = run_dir / "canonical_family_manifest.tsv"
    if canonical.sha256(manifest_path) != identity.get("canonical_manifest_sha256"):
        raise ValueError("Run-local canonical family manifest hash mismatch")
    input_provenance = input_root / "provenance.json"
    input_inventory = canonical.read_json(input_provenance)
    if (input_inventory.get("n_loci") != 2495 or input_inventory.get("n_traits") != len(traits) or
            len(input_inventory.get("records", [])) != 2495 * (len(traits) + 1)):
        raise ValueError("Materialized input provenance is not the exact frozen locus/trait inventory")
    if not reference_provenance.is_file():
        raise FileNotFoundError(reference_provenance)

    locus_by_id = {row["LOC"]: row for row in loci}
    observed_paths = sorted((run_dir / "receipts").glob("locus_*.json"),
                            key=lambda p: p.stem.removeprefix("locus_"))
    unknown_receipts = []
    valid_loci: list[dict[str, str]] = []
    invalid_receipts = []
    for path in observed_paths:
        locus_id = path.stem.removeprefix("locus_")
        locus = locus_by_id.get(locus_id)
        if locus is None:
            unknown_receipts.append(locus_id)
            continue
        if canonical.verify_locus(run_dir, locus, traits, family["analysis_id"], family_sha,
                                  execution_sha, input_root, reference_provenance, script_sha):
            valid_loci.append(locus)
        else:
            invalid_receipts.append(locus_id)

    counts_by_trait = {trait: collections.Counter() for trait in traits}
    reasons = collections.Counter()
    for locus in valid_loci:
        rows = canonical.load_cells(canonical.result_path(run_dir, locus["LOC"]))
        for row in rows:
            counts_by_trait[row["phen"]][row["status"]] += 1
            if row["status"] != "TESTED":
                reasons[f"{row['status']}:{row.get('reason', '')}"] += 1
    summary = summarize(len(valid_loci), len(loci), traits, counts_by_trait,
                        float(family["execution"]["maximum_univariate_untested_fraction"]))
    summary.update({
        "analysis_id": family["analysis_id"],
        "run_id": run_id,
        "worker_count": execution["worker_count"],
        "observed_receipt_files": len(observed_paths),
        "verified_locus_ids": [locus["LOC"] for locus in valid_loci],
        "invalid_receipts": sorted(invalid_receipts),
        "unknown_receipt_ids": sorted(unknown_receipts),
        "invalid_receipt_count": len(invalid_receipts) + len(unknown_receipts),
        "non_tested_reasons": dict(sorted(reasons.items())),
        "family_lock_sha256": family_sha,
        "execution_lock_sha256": execution_sha,
        "canonical_script_sha256": script_sha,
        "input_provenance_sha256": canonical.sha256(input_provenance),
        "reference_provenance_sha256": canonical.sha256(reference_provenance),
        "read_only": True,
    })
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference-provenance", type=Path, required=True)
    parser.add_argument("--family-lock", type=Path, default=DEFAULT_FAMILY)
    parser.add_argument("--execution-lock", type=Path, default=DEFAULT_EXECUTION)
    args = parser.parse_args()
    result = audit(args.run_dir, args.input_root, args.reference_provenance,
                   args.family_lock, args.execution_lock)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 2 if result["invalid_receipt_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
