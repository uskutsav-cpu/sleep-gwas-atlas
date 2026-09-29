#!/usr/bin/env python3
"""Read-only production throughput profile for canonical LAVA v3 receipts."""
from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
import audit_lava_canonical_v3_partial as partial
import run_lava_canonical_v3 as canonical

METRIC = re.compile(
    r"BRAIN6_CANONICAL_LOCUS_METRICS\s+(\S+)\s+"
    r"compute_wall_seconds=([0-9.]+)\s+write_seconds=([0-9.]+)"
)
READY = re.compile(r"BRAIN6_WORKER_READY\s+.*?reference_load_seconds=([0-9.]+)")


def percentile90(values: list[float]) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    return statistics.quantiles(values, n=10)[8]


def window_rates(receipt_times: list[float], remaining_loci: int,
                 now_timestamp: float, windows_minutes: tuple[int, ...] = (15, 60, 120)
                 ) -> list[dict[str, Any]]:
    result = []
    for minutes in windows_minutes:
        if minutes < 1:
            raise ValueError("Window lengths must be positive")
        count = sum(ts >= now_timestamp - minutes * 60 for ts in receipt_times)
        per_hour = count * 60 / minutes
        result.append({
            "window_minutes": minutes,
            "verified_receipts_in_window": count,
            "receipts_per_hour": per_hour,
            "remaining_loci_projection_hours": remaining_loci / per_hour if per_hour else None,
        })
    return result


def profile(run_dir: Path, audit_result: dict[str, Any], now_timestamp: float | None = None
            ) -> dict[str, Any]:
    """Summarize metrics only for receipts accepted by the full verifier."""
    if now_timestamp is None:
        now_timestamp = time.time()
    metric_rows: dict[str, list[tuple[float, float]]] = defaultdict(list)
    ready_seconds: list[float] = []
    bound_logs: set[Path] = set()
    receipt_times: list[float] = []
    missing_metric_loci: list[str] = []
    duplicate_metric_loci: list[str] = []
    valid_receipt_loci: list[str] = []

    # Use the verifier's accepted-ID snapshot. Rescanning the live directory
    # here can observe a new atomic receipt that was not included in the audit's
    # cell/status totals, producing internally inconsistent progress counts.
    for locus_id in audit_result["verified_locus_ids"]:
        receipt_path = canonical.receipt_path(run_dir, locus_id)
        receipt = canonical.read_json(receipt_path)
        valid_receipt_loci.append(locus_id)
        receipt_times.append(receipt_path.stat().st_mtime)
        config_path = Path(receipt["worker_config_path"])
        log_path = config_path.with_suffix(".log")
        if not log_path.is_file():
            missing_metric_loci.append(locus_id)
            continue
        bound_logs.add(log_path)
        for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
            metric = METRIC.search(line)
            if metric and metric.group(1) == locus_id:
                metric_rows[locus_id].append((float(metric.group(2)), float(metric.group(3))))

    for log_path in bound_logs:
        for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
            match = READY.search(line)
            if match:
                ready_seconds.append(float(match.group(1)))

    for locus_id in valid_receipt_loci:
        count = len(metric_rows.get(locus_id, []))
        if count == 0:
            missing_metric_loci.append(locus_id)
        elif count > 1:
            duplicate_metric_loci.append(locus_id)

    unique_rows = [rows[0] for rows in metric_rows.values() if len(rows) == 1]
    compute = [row[0] for row in unique_rows]
    write = [row[1] for row in unique_rows]
    remaining = max(0, int(audit_result["planned_loci"]) - len(valid_receipt_loci))
    first = min(receipt_times) if receipt_times else None
    last = max(receipt_times) if receipt_times else None
    return {
        "analysis_id": audit_result["analysis_id"],
        "run_id": audit_result["run_id"],
        "state": audit_result["state"],
        "workers": audit_result["worker_count"],
        "verified_receipts": len(valid_receipt_loci),
        "verified_cells": audit_result["verified_cells"],
        "planned_cells": audit_result["planned_cells"],
        "status_counts": audit_result["statuses"],
        "status_by_trait": audit_result["status_by_trait"],
        "untested_lower_bound_cells": audit_result["untested_lower_bound_cells"],
        "maximum_allowed_untested_cells": audit_result["maximum_allowed_untested_cells"],
        "family_qc_guaranteed_fail": audit_result["family_qc_guaranteed_fail"],
        "invalid_receipts": audit_result["invalid_receipt_count"],
        "remaining_loci": remaining,
        "receipt_time_range_utc": {
            "first": datetime.fromtimestamp(first, timezone.utc).isoformat() if first else None,
            "last": datetime.fromtimestamp(last, timezone.utc).isoformat() if last else None,
        },
        "throughput_windows": window_rates(receipt_times, remaining, now_timestamp),
        "receipt_bound_locus_metrics": len(unique_rows),
        "missing_metric_loci": sorted(set(missing_metric_loci)),
        "duplicate_metric_loci": sorted(set(duplicate_metric_loci)),
        "compute_wall_seconds": {
            "mean": statistics.mean(compute) if compute else None,
            "median": statistics.median(compute) if compute else None,
            "p90": percentile90(compute),
            "sum_worker_hours": sum(compute) / 3600,
        },
        "result_write_seconds_median": statistics.median(write) if write else None,
        "reference_load_seconds_median": statistics.median(ready_seconds) if ready_seconds else None,
        "bound_worker_logs": len(bound_logs),
        "receipt_metrics_join": "receipt.worker_config_path -> its sibling worker log; only verified receipt locus ids counted",
        "read_only": True,
        "projection_note": "Window-based arithmetic projections are descriptive; runtime varies by locus and workload.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference-provenance", type=Path, required=True)
    parser.add_argument("--family-lock", type=Path, default=partial.DEFAULT_FAMILY)
    parser.add_argument("--execution-lock", type=Path, default=partial.DEFAULT_EXECUTION)
    args = parser.parse_args()
    audit_result = partial.audit(args.run_dir, args.input_root, args.reference_provenance,
                                 args.family_lock, args.execution_lock)
    report = profile(args.run_dir, audit_result)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 2 if report["invalid_receipts"] or report["missing_metric_loci"] or report["duplicate_metric_loci"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
