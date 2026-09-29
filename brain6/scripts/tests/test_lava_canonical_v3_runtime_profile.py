from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/profile_lava_canonical_v3_runtime.py"
spec = importlib.util.spec_from_file_location("lava_v3_runtime_profile", SCRIPT)
profile = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(profile)


def test_window_rates_use_receipt_timestamps_and_remaining_loci():
    now = 10_000.0
    result = profile.window_rates([now - 30, now - 600, now - 8_000], 900, now,
                                  windows_minutes=(15, 60))
    assert result[0]["verified_receipts_in_window"] == 2
    assert result[0]["receipts_per_hour"] == 8
    assert result[0]["remaining_loci_projection_hours"] == 112.5
    assert result[1]["verified_receipts_in_window"] == 2
    assert result[1]["receipts_per_hour"] == 2


def test_percentile_handles_small_and_empty_samples():
    assert profile.percentile90([]) is None
    assert profile.percentile90([4.5]) == 4.5
    assert profile.percentile90([1.0, 2.0, 3.0, 4.0]) > 3.0


def test_nonpositive_throughput_has_no_projection():
    result = profile.window_rates([], 100, 10_000, windows_minutes=(5,))[0]
    assert result["receipts_per_hour"] == 0
    assert result["remaining_loci_projection_hours"] is None


def test_profile_uses_the_audited_receipt_snapshot(tmp_path):
    run_dir = tmp_path / "run"
    receipts = run_dir / "receipts"
    receipts.mkdir(parents=True)
    config = tmp_path / "worker.json"
    config.write_text("{}\n")
    (config.with_suffix(".log")).write_text(
        "BRAIN6_CANONICAL_LOCUS_METRICS 1 compute_wall_seconds=2.0 write_seconds=0.1\n"
    )
    (receipts / "locus_1.json").write_text(json.dumps({"worker_config_path": str(config)}))
    # Simulate a receipt arriving after the audit snapshot. The profiler must
    # ignore it until a subsequent audit verifies its seven cell results.
    (receipts / "locus_2.json").write_text("{}\n")
    audit_result = {
        "invalid_receipts": [],
        "unknown_receipt_ids": [],
        "verified_locus_ids": ["1"],
        "planned_loci": 2,
        "analysis_id": "test",
        "run_id": "test-run",
        "state": "PARTIAL",
        "worker_count": 4,
        "verified_cells": 7,
        "planned_cells": 14,
        "statuses": {"TESTED": 7, "NOT_RUN": 0, "FAILED": 0},
        "status_by_trait": {},
        "untested_lower_bound_cells": 0,
        "maximum_allowed_untested_cells": 0,
        "family_qc_guaranteed_fail": False,
        "invalid_receipt_count": 0,
    }

    result = profile.profile(run_dir, audit_result, now_timestamp=1000.0)

    assert result["verified_receipts"] == 1
    assert result["verified_cells"] == 7
    assert result["remaining_loci"] == 1
