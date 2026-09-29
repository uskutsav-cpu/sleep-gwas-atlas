import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "compare_lava_baseline_roundoff.py"
spec = importlib.util.spec_from_file_location("lava_roundoff_comparison", SCRIPT)
comparison = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(comparison)


def test_number_treats_missing_values_as_none_and_parses_finite_values():
    assert comparison.number("NA") is None
    assert comparison.number("") is None
    assert comparison.number("0.125") == 0.125


def test_number_rejects_non_finite_values():
    with pytest.raises(ValueError, match="non-finite"):
        comparison.number("inf")


def test_stable_file_inventory_hash_is_deterministic_and_order_sensitive():
    records = [{"path": "a", "sha256": "1"}, {"path": "b", "sha256": "2"}]
    assert comparison.stable_file_inventory_hash(records) == comparison.stable_file_inventory_hash(records)
    assert comparison.stable_file_inventory_hash(records) != comparison.stable_file_inventory_hash(list(reversed(records)))


def test_validate_audit_requires_complete_valid_receipt_families(tmp_path):
    audit = {
        "audit_status": "PASS_RECEIPT_INTEGRITY",
        "lava": {"run_id": "base", "receipt_errors": 0, "complete": True, "verified_loci": 2, "planned_loci": 2},
        "lava_roundoff": {"run_id": "fix", "receipt_errors": 0, "complete": True, "verified_loci": 2, "planned_loci": 2},
    }
    path = tmp_path / "audit.json"
    path.write_text(json.dumps(audit))
    assert comparison.validate_audit(path, "base", "fix")["audit_status"] == "PASS_RECEIPT_INTEGRITY"
    audit["lava_roundoff"]["complete"] = False
    path.write_text(json.dumps(audit))
    with pytest.raises(ValueError, match="not a complete receipt-verified run"):
        comparison.validate_audit(path, "base", "fix")
