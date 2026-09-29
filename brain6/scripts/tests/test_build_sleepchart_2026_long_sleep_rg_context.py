from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "build_sleepchart_2026_long_sleep_rg_context.py"
SPEC = importlib.util.spec_from_file_location("sleepchart_context", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_sleepchart_context_is_rounded_prior_context_not_replication_or_power_upgrade():
    rows, provenance = MODULE.build()
    assert {row["brain_disorder"] for row in rows} == {"adhd", "mdd", "scz", "bipolar"}
    assert len(rows) == 4
    assert all(row["direction_concordant"] == "True" for row in rows)
    assert all(row["independent_replication"] == "False" for row in rows)
    assert all(row["power_replacement_eligible"] == "False" for row in rows)
    assert provenance["scope"]["novelty_or_fdr_changed"] is False
    assert provenance["source"]["sleep_gwas"]["cases"] == 25049
    assert provenance["source"]["sleep_gwas"]["summary_statistics_downloaded"] is False

