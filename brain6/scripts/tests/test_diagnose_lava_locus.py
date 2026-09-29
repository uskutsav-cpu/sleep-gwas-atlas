import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6" / "scripts"))

import diagnose_lava_locus as diagnostic
from diagnose_lava_locus import diagnose


def test_diagnostic_requires_a_complete_receipt(tmp_path):
    (tmp_path / "receipt.json").write_text(json.dumps({"locus_id": "1997", "run_id": "x"}))
    with pytest.raises(ValueError, match="not complete|identity mismatch"):
        diagnose(tmp_path)


def test_diagnostic_selects_requested_pair_without_mutating_saved_config(tmp_path, monkeypatch):
    from run_lava_family import EXPECTED_PAIRS

    (tmp_path / "receipt.json").write_text(json.dumps({"locus_id": "1997", "run_id": "x"}))
    (tmp_path / "worker_config.json").write_text(json.dumps({
        "locus_id": "1997", "pairs": [{"pair_id": p} for p in EXPECTED_PAIRS],
    }))
    monkeypatch.setattr(diagnostic, "verify_unit", lambda *args: {"status": "COMPLETE"})

    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        observed["config"] = json.loads(Path(command[-1]).read_text())
        return type("Result", (), {"returncode": 0})()

    monkeypatch.setattr(diagnostic.subprocess, "run", fake_run)
    assert diagnose(tmp_path, "longsleep__scz") == 0
    assert [p["pair_id"] for p in observed["config"]["pairs"]] == ["longsleep__scz"]
    assert json.loads((tmp_path / "worker_config.json").read_text())["pairs"] == [
        {"pair_id": p} for p in EXPECTED_PAIRS
    ]
