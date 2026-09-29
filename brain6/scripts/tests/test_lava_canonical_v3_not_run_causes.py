from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/build_lava_canonical_v3_not_run_causes.py"
spec = importlib.util.spec_from_file_location("lava_v3_not_run_causes", SCRIPT)
builder = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(builder)


def test_extracts_one_reason_bound_row_per_not_run_cell():
    rows = [
        {"phen": "adhd", "locus_id": "12", "chromosome": "1", "status": "NOT_RUN",
         "reason": "LOW_LOCAL_H2_UNDERPOWERED"},
        {"phen": "bipolar", "locus_id": "12", "chromosome": "1", "status": "NOT_RUN",
         "reason": "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS"},
        {"phen": "insomnia", "locus_id": "2", "chromosome": "1", "status": "TESTED", "reason": ""},
    ]
    causes = builder.extract_cause_rows(rows)
    assert len(causes) == 2
    assert {row["cause_category"] for row in causes} == {
        "LOCAL_H2_SUPPORT_GATE", "INSUFFICIENT_SHARED_REFERENCE_VARIANTS"
    }
    assert all("not evidence for a null association" in row["interpretation"] for row in causes)
    assert builder.render_tsv(causes) == builder.render_tsv(causes)


def test_rejects_unknown_not_run_reason():
    with pytest.raises(ValueError, match="Unclassified canonical NOT_RUN reason"):
        builder.extract_cause_rows([
            {"phen": "adhd", "locus_id": "1", "status": "NOT_RUN", "reason": "NEW_REASON"}
        ])


def test_rejects_duplicate_trait_locus_identity():
    row = {"phen": "adhd", "locus_id": "1", "status": "TESTED", "reason": ""}
    with pytest.raises(ValueError, match="duplicate cell identity"):
        builder.extract_cause_rows([row, dict(row)])


def test_immutable_writer_refuses_to_replace_different_bytes(tmp_path):
    path = tmp_path / "cause.tsv"
    builder.write_immutable(path, b"first\n")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        builder.write_immutable(path, b"changed\n")
    assert path.read_bytes() == b"first\n"
