from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/validate_current_outputs.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("validate_current_outputs", SCRIPT)
validator = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(validator)


def test_verified_mdd2025_sidecar_reports_documented_pgc29_overlap():
    result = validator.validate_pgc_mdd2025_cohort_overlap()

    assert result["status"] == "PASS_COHORT_OVERLAP_NOT_INDEPENDENT_REPLICATION"
    assert result["sidecar_rows"] == 74
    assert result["matched_pgc29_cohorts"] == 9


def test_sidecar_tampering_is_rejected(tmp_path):
    source = ROOT / "brain6/qc/replication_sources"
    target = tmp_path / "brain6/qc/replication_sources"
    target.mkdir(parents=True)
    for name in ("pgc-mdd2025_no23andMe-noUKBB_eur_v3.49.24.11.txt",
                 "pgc_mdd2025_noUKBB_overlap_audit.json"):
        shutil.copy2(source / name, target / name)
    sidecar = target / "pgc-mdd2025_no23andMe-noUKBB_eur_v3.49.24.11.txt"
    sidecar.write_bytes(sidecar.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="source checksum is invalid"):
        validator.validate_pgc_mdd2025_cohort_overlap(tmp_path)


def test_audit_cannot_drop_a_documented_overlap(tmp_path):
    source = ROOT / "brain6/qc/replication_sources"
    target = tmp_path / "brain6/qc/replication_sources"
    target.mkdir(parents=True)
    for name in ("pgc-mdd2025_no23andMe-noUKBB_eur_v3.49.24.11.txt",
                 "pgc_mdd2025_noUKBB_overlap_audit.json"):
        shutil.copy2(source / name, target / name)
    audit_path = target / "pgc_mdd2025_noUKBB_overlap_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["matched_legacy_pgc29_cohorts"] = audit["matched_legacy_pgc29_cohorts"][:-1]
    audit_path.write_text(json.dumps(audit), encoding="utf-8")

    with pytest.raises(ValueError, match="match_ids|match_count"):
        validator.validate_pgc_mdd2025_cohort_overlap(tmp_path)
