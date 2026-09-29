from __future__ import annotations

import collections
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/audit_lava_canonical_v3_partial.py"
spec = importlib.util.spec_from_file_location("lava_v3_partial_audit", SCRIPT)
audit = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(audit)


def counts(**values: int) -> collections.Counter[str]:
    return collections.Counter(values)


def test_partial_summary_uses_full_family_denominator_and_trait_counts():
    traits = ("a", "b")
    result = audit.summarize(2, 4, traits, {
        "a": counts(NOT_RUN=2),
        "b": counts(TESTED=1, NOT_RUN=1),
    }, 0.25)
    assert result["state"] == "PARTIAL"
    assert result["planned_cells"] == 8
    assert result["verified_cells"] == 4
    assert result["missing_or_invalid_cells"] == 4
    assert result["untested_lower_bound_cells"] == 3
    assert result["currently_uncompleted_or_untested_cells"] == 7
    assert result["maximum_allowed_untested_cells"] == 2
    assert result["family_qc_guaranteed_fail"] is True
    assert result["family_qc_pass"] is False
    assert result["status_by_trait"]["a"]["NOT_RUN"] == 2


def test_partial_under_cap_is_not_reported_as_family_pass():
    result = audit.summarize(1, 2, ("a",), {"a": counts(TESTED=1)}, 0.5)
    assert result["untested_lower_bound_cells"] == 0
    assert result["currently_uncompleted_or_untested_cells"] == 1
    assert result["family_qc_guaranteed_fail"] is False
    assert result["family_qc_pass"] is False


def test_complete_family_can_pass_frozen_untested_gate():
    result = audit.summarize(2, 2, ("a", "b"), {
        "a": counts(TESTED=2),
        "b": counts(TESTED=1, NOT_RUN=1),
    }, 0.25)
    assert result["state"] == "COMPLETE"
    assert result["untested_lower_bound_cells"] == 1
    assert result["maximum_allowed_untested_cells"] == 1
    assert result["family_qc_pass"] is True


def test_partial_status_summary_rejects_incomplete_trait_coverage():
    try:
        audit.summarize(1, 2, ("a", "b"), {"a": counts(TESTED=1)}, 0.05)
    except ValueError as error:
        assert "every canonical trait" in str(error)
    else:
        raise AssertionError("Incomplete trait status summary must be rejected")
