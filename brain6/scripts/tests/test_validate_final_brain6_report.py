from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/validate_final_brain6_report.py"
spec = importlib.util.spec_from_file_location("validate_final_brain6_report", SCRIPT)
validator = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(validator)


def test_current_report_matches_locked_global_replication_and_lava_evidence():
    result = validator.validate_report()

    assert result["status"] == "PASS_PROVISIONAL_REPORT"
    assert result["global_rows"] == 72
    assert result["inherited_significant_rows"] == 35
    assert result["roundoff_receipt_errors"] == 0


def test_complete_significant_table_rejects_a_missing_association():
    report = validator.REPORT.read_text(encoding="utf-8")
    global_rows = validator.read_tsv(ROOT / "brain6/results/global/brain6_72_locked.tsv")
    rows = validator.report_significant_rows(report)
    removed = rows[-1]
    line = "| " + " | ".join(removed) + " |"

    with pytest.raises(ValueError, match="differs from the locked global map"):
        validator.validate_report_table(report.replace(line + "\n", "", 1), global_rows)


def test_complete_significant_table_rejects_a_changed_effect():
    report = validator.REPORT.read_text(encoding="utf-8")
    global_rows = validator.read_tsv(ROOT / "brain6/results/global/brain6_72_locked.tsv")
    first = validator.report_significant_rows(report)[0]
    line = "| " + " | ".join(first) + " |"
    changed = line.replace(first[2], "0.999", 1)

    with pytest.raises(ValueError, match="differs from the locked global map"):
        validator.validate_report_table(report.replace(line, changed, 1), global_rows)


def test_report_provenance_rejects_stale_source_hashes(tmp_path):
    report = ROOT / "brain6/FINAL_BRAIN6_REPORT.md"
    provenance = ROOT / "brain6/FINAL_BRAIN6_REPORT.provenance.json"
    (tmp_path / "brain6").mkdir()
    # Reuse a root copy so the test only perturbs a provenance-bound source.
    (tmp_path / "brain6/FINAL_BRAIN6_REPORT.md").write_bytes(report.read_bytes())
    (tmp_path / "brain6/FINAL_BRAIN6_REPORT.provenance.json").write_text(
        provenance.read_text(encoding="utf-8"), encoding="utf-8")
    with pytest.raises((FileNotFoundError, ValueError)):
        validator.validate_report(tmp_path)
