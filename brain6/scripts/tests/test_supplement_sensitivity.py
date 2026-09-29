import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/build_supplement_sensitivity.py"
spec = importlib.util.spec_from_file_location("supplement_sensitivity", SCRIPT)
supplement = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(supplement)


def test_table_s16_is_a_provenance_bound_interim_copy_of_claim_ledger(monkeypatch, tmp_path):
    source = tmp_path / "brain6/results/sensitivity/headline_claim_sensitivity.tsv"
    source.parent.mkdir(parents=True)
    source.write_bytes(supplement.SOURCE.read_bytes())
    script = tmp_path / "brain6/scripts/build_supplement_sensitivity.py"
    script.parent.mkdir(parents=True)
    script.write_bytes(supplement.SCRIPT.read_bytes())
    monkeypatch.setattr(supplement, "ROOT", tmp_path)
    monkeypatch.setattr(supplement, "SOURCE", source)
    monkeypatch.setattr(supplement, "SCRIPT", script)
    monkeypatch.setattr(supplement, "OUTPUT", tmp_path / "table_S16.tsv")
    monkeypatch.setattr(supplement, "PROVENANCE", tmp_path / "table_S16.provenance.json")
    provenance = supplement.main()
    source = supplement.SOURCE.read_bytes()
    assert supplement.OUTPUT.read_bytes() == source
    assert provenance["status"] == "INTERIM_PARTIAL_SENSITIVITY_COVERAGE"
    assert provenance["claim_rows"] == 12
    assert provenance["source_sha256"] == supplement.sha256(source)
    assert provenance["script_sha256"] == supplement.sha256(SCRIPT.read_bytes())
    stored = json.loads(supplement.PROVENANCE.read_text())
    assert stored == provenance


def configure_temporary_ledger(monkeypatch, tmp_path, mutate):
    source_rows = list(csv.DictReader(supplement.SOURCE.open(encoding="utf-8"), delimiter="\t"))
    mutate(source_rows)
    source = tmp_path / "sensitivity.tsv"
    with source.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=source_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(source_rows)
    script = tmp_path / "builder.py"
    script.write_text("test builder identity\n", encoding="utf-8")
    monkeypatch.setattr(supplement, "ROOT", tmp_path)
    monkeypatch.setattr(supplement, "SOURCE", source)
    monkeypatch.setattr(supplement, "SCRIPT", script)
    monkeypatch.setattr(supplement, "OUTPUT", tmp_path / "table_S16.tsv")
    monkeypatch.setattr(supplement, "PROVENANCE", tmp_path / "table_S16.provenance.json")
    return source


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda rows: [row.pop("LAVA_UKB_v1_1") for row in rows], "unexpected schema"),
        (lambda rows: rows.__setitem__(1, {**rows[1], "claim_id": rows[0]["claim_id"]}),
         "duplicate claim identities"),
        (lambda rows: rows[0].__setitem__("overall_sensitivity_status", "ROBUST"),
         "must not promote"),
    ],
    ids=("missing-column", "duplicate-claim", "robust-status"),
)
def test_table_s16_rejects_invalid_or_promoted_claim_ledgers(
    monkeypatch, tmp_path, mutate, message
):
    configure_temporary_ledger(monkeypatch, tmp_path, mutate)
    with pytest.raises(ValueError, match=message):
        supplement.main()
    assert not supplement.OUTPUT.exists()
    assert not supplement.PROVENANCE.exists()
