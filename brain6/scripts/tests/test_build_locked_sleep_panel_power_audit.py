import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "build_locked_sleep_panel_power_audit.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("locked_sleep_panel_power_audit", SCRIPT)
audit = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(audit)


def test_locked_sleep_panel_is_read_from_family_config():
    traits = audit.locked_sleep_traits()
    assert len(traits) == 12
    assert traits[:4] == ["insomnia", "sleepdur", "shortsleep", "longsleep"]
    assert len(set(traits)) == 12


def test_effective_sample_size_is_only_computed_for_case_control_traits():
    assert audit.effective_n("34184", "305742", "339926") == "122985.4"
    assert audit.effective_n("NA", "NA", "446118") == "446118"
    assert audit.effective_n("0", "305742", "305742") == "305742"


def test_lava_summary_separates_low_h2_from_reference_min_k():
    rows = [{"status": "TESTED", "h2.obs": "0.001", "reason": ""}]
    rows += [{"status": "NOT_RUN", "h2.obs": "NA", "reason": "LOW_LOCAL_H2_UNDERPOWERED"}]
    rows += [{"status": "NOT_RUN", "h2.obs": "NA", "reason": "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS"}]
    rows += [{"status": "TESTED", "h2.obs": "0.003", "reason": ""}] * (2495 - len(rows))

    result = audit.summarize_cells(rows)

    assert result["tested"] == 2493
    assert result["not_run"] == 2
    assert result["low_h2"] == 1
    assert result["shared_reference_min_k"] == 1
    assert result["failed"] == 0
    assert result["h2_observed_n"] == 2493
    assert result["h2_median"] == 0.003


def test_lava_summary_rejects_incomplete_locus_set():
    with pytest.raises(ValueError, match="2,495 locus rows"):
        audit.summarize_cells([])
