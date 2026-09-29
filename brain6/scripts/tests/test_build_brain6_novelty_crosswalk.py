from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_brain6_novelty_crosswalk as crosswalk


def test_crosswalk_covers_exactly_the_inherited_significant_pairs():
    rows, summary = crosswalk.build_rows()
    assert len(rows) == 35
    assert len({(row["sleep_trait"], row["brain_disorder"]) for row in rows}) == 35
    assert all(row["significance_under_original_396_family"].lower() == "true" for row in rows)
    assert summary["global_pair_count"] == 72
    assert summary["inherited_significant_pair_count"] == 35


def test_crosswalk_keeps_novelty_classes_and_overlap_limit_explicit():
    rows, summary = crosswalk.build_rows()
    assert Counter(row["novelty_classification"] for row in rows) == {
        "DIRECT_RG_PREVIOUSLY_REPORTED": 25,
        "DIRECT_RG_REPLICATION_DIFFERENT_DATASET": 9,
        "MR_ONLY": 1,
    }
    assert summary["apparent_novelty_count"] == 0
    assert summary["parent_audit_search_as_of_date"] == "2026-08-28"
    assert "not a fresh current-date search" in summary["interpretation"]
    assert all("first-ever claim" in row["claim_limit"].lower() for row in rows)


def test_persisted_crosswalk_hashes_and_rows_validate():
    result = crosswalk.validate()
    assert result["status"] == "PASS_DATED_LITERATURE_CROSSWALK_NO_NOVELTY_CLAIM"
    assert result["direct_rg_count"] == 34
    assert result["output_sha256"] == crosswalk.sha256(crosswalk.OUTPUT)
