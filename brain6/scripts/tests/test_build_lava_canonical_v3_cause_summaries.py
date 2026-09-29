import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "build_lava_canonical_v3_cause_summaries.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("lava_v3_cause_summaries", SCRIPT)
summaries = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(summaries)


def full_family_rows():
    traits = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
    rows = []
    cause_remaining = dict(zip(summaries.CAUSES, (3564, 154, 2), strict=True))
    for trait in traits:
        for locus_num in range(1, 2496):
            status, reason = "TESTED", ""
            next_cause = next((cause for cause in summaries.CAUSES if cause_remaining[cause]), None)
            if next_cause is not None:
                status = "NOT_RUN"
                reason = next_cause
                cause_remaining[next_cause] -= 1
            rows.append({
                "phen": trait, "locus_id": str(locus_num),
                "chromosome": str((locus_num - 1) // 114 + 1),
                "status": status, "reason": reason,
            })
    return rows


def test_trait_and_locus_breakdowns_account_for_every_frozen_cell():
    traits, loci, summary = summaries.summarize(full_family_rows())
    assert len(traits) == 7
    assert len(loci) == 2495
    assert summary["cells"] == 17465
    assert summary["not_run_cells"] == 3720
    assert summary["tested_cells"] == 13745
    assert summary["failed_cells"] == 0
    assert summary["cause_counts"] == dict(zip(summaries.CAUSES, (3564, 154, 2), strict=True))
    assert sum(row["planned_cells"] for row in traits) == 17465
    assert sum(row["planned_cells"] for row in loci) == 17465
    assert all(row["planned_cells"] == 7 for row in loci)


def test_summaries_reject_missing_or_duplicate_cells():
    rows = full_family_rows()
    with pytest.raises(ValueError, match="dimensions"):
        summaries.summarize(rows[:-1])
    rows.append(rows[0].copy())
    with pytest.raises(ValueError, match="duplicate"):
        summaries.summarize(rows)


def test_immutable_write_refuses_a_different_existing_payload(tmp_path):
    output = tmp_path / "table.tsv"
    summaries.immutable_write(output, b"frozen\n")
    summaries.immutable_write(output, b"frozen\n")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        summaries.immutable_write(output, b"changed\n")
