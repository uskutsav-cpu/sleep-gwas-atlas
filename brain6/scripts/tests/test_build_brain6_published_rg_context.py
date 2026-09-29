from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from build_brain6_published_rg_context import build_rows


def test_published_rows_map_to_unique_locked_brain6_pairs() -> None:
    rows = build_rows()
    keys = [(row["brain6_sleep_trait"], row["brain6_disorder"]) for row in rows]
    assert len(rows) == 21
    assert len(set(keys)) == 21
    assert sum(row["brain6_original_family_significant"].lower() == "true" for row in rows) == 18


def test_published_comparison_is_descriptive_and_not_independent_replication() -> None:
    rows = build_rows()
    assert all("NOT_ESTABLISHED" in row["comparison_independence"] for row in rows)
    assert all("not independent replication" in row["interpretation"] for row in rows)
    assert sum(row["intercept_overlap_flag"] == "true" for row in rows) == 6


def test_published_rg_direction_agrees_for_all_21_pairs() -> None:
    assert all(row["direction_concordant"] == "true" for row in build_rows())
