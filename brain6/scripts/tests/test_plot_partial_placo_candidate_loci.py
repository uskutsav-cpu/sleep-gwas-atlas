from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/plot_partial_placo_candidate_loci.py"
spec = importlib.util.spec_from_file_location("partial_placo_candidate_figure", SCRIPT)
builder = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(builder)


def test_source_bound_partial_loci_have_all_four_available_pairs():
    rows = builder.load_rows()
    assert len(rows) == 19
    assert {row["pair_id"] for row in rows} == set(builder.EXPECTED_PAIRS)
    assert {row["locus_status"] for row in rows} == {"PLACO_ONLY_CANDIDATE_LOCUS"}


def test_refuses_full_family_claim_when_source_provenance_changes(tmp_path):
    source = json.loads(builder.PROVENANCE_PATH.read_text(encoding="utf-8"))
    source["status"] = "PASS"
    path = tmp_path / "provenance.json"
    path.write_text(json.dumps(source), encoding="utf-8")
    with pytest.raises(ValueError, match="partial family"):
        builder.load_rows(provenance_path=path)


def test_render_is_source_driven_and_never_overwrites(tmp_path):
    png, pdf = builder.render(builder.load_rows(), tmp_path)
    assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert pdf.read_bytes().startswith(b"%PDF-")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        builder.render(builder.load_rows(), tmp_path)
