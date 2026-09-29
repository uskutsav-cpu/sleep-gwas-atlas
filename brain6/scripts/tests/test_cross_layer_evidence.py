from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/build_cross_layer_evidence.py"
spec = importlib.util.spec_from_file_location("cross_layer_evidence", SCRIPT)
builder = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(builder)


@pytest.mark.parametrize("values,expected", [
    ({"global_significant": True, "local_promoted": True, "replication_admissible": True,
      "variant_or_molecular_layer": True, "any_evidence": True}, "MULTI_LAYER_SUPPORTED"),
    ({"global_significant": True, "local_promoted": False, "replication_admissible": True,
      "variant_or_molecular_layer": False, "any_evidence": True}, "GENETIC_ONLY"),
    ({"global_significant": False, "local_promoted": True, "replication_admissible": False,
      "variant_or_molecular_layer": False, "any_evidence": True}, "LOCAL_ONLY"),
    ({"global_significant": False, "local_promoted": False, "replication_admissible": True,
      "variant_or_molecular_layer": False, "any_evidence": True}, "REPLICATION_ONLY"),
    ({"global_significant": False, "local_promoted": False, "replication_admissible": False,
      "variant_or_molecular_layer": False, "any_evidence": True}, "SUGGESTIVE"),
    ({"global_significant": False, "local_promoted": False, "replication_admissible": False,
      "variant_or_molecular_layer": False, "any_evidence": False}, "UNRESOLVED"),
])
def test_frozen_evidence_category_precedence(values, expected):
    assert builder.assign_category(**values) == expected


def test_primary_pair_table_preserves_partial_and_blocked_layers():
    rows, inputs = builder.build_rows()
    assert len(rows) == 5
    assert len({row["pair_id"] for row in rows}) == 5
    assert {row["overall_evidence_category"] for row in rows} == {"GENETIC_ONLY"}
    adhd = next(row for row in rows if row["pair_id"] == "insomnia__adhd")
    assert adhd["replication_class"] == "DIRECTIONAL_REPLICATION"
    assert adhd["placo_pair_status"] == "PROTECTED_TRACK_B_ABSENT_NOT_RERUN"
    assert all(row["canonical_v3_local_rg_family_status"] == "FAILED_QC_NOT_PROMOTED" for row in rows)
    assert all(row["fine_mapping_status"].startswith("NOT_RUN_") for row in rows)
    assert "brain6/config/cross_layer_evidence_rules_v1.json" in inputs


def test_table_render_is_deterministic():
    rows, _ = builder.build_rows()
    assert builder.render(rows) == builder.render(rows)
