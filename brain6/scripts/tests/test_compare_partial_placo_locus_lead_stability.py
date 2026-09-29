from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "compare_partial_placo_locus_lead_stability.py"
SPEC = importlib.util.spec_from_file_location("placo_lead_stability", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_lead_stability_comparison_is_receipt_bound_and_diagnostic_only() -> None:
    provenance = MODULE.build()
    assert provenance["status"] == "PASS_DIAGNOSTIC_ONLY_NO_PROMOTION"
    assert provenance["scope"]["baseline_leads"] == 21
    assert provenance["scope"]["factor_normalized_sensitivity_leads"] == 21
    assert provenance["scope"]["range_edges_excluded_leads"] == 25
    assert provenance["scope"]["stable_all_three"] == 21
    assert provenance["scope"]["range_filter_only"] == 4
    assert provenance["scope"]["full_five_track_family"] is False
    assert provenance["scope"]["genotype_level_ld_validation"] is False
    with MODULE.OUT_PATH.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 25
    range_only = {row["lead_variant"] for row in rows if row["stability_class"] == "RANGE_FILTER_ONLY"}
    assert range_only == {"rs6430538", "rs199535", "rs199534", "rs4938023"}
    assert all(row["interpretation"].startswith("Diagnostic PLACO candidate lead only") for row in rows)
