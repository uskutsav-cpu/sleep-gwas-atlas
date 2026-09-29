import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import clump_placo_candidates as clump
from rebuild_partial_placo_clumps_valid_ld import exception_rows, split_ld_edges


def test_out_of_range_ld_is_excluded_without_clipping():
    valid, invalid = split_ld_edges({("a", "a"): 1.0, ("a", "b"): 1.001,
                                     ("b", "c"): -1.0001, ("c", "c"): -1.0})
    assert valid == {("a", "a"): 1.0, ("c", "c"): -1.0}
    assert invalid == {("a", "b"): 1.001, ("b", "c"): -1.0001}


def test_range_exception_is_reported_for_in_window_pair():
    rows = {"pair": [
        {"SNP": "a", "CHR": 1, "BP": 10_000},
        {"SNP": "b", "CHR": 1, "BP": 100_000},
    ]}
    reference = {(1, "a"): {"BP": 10_000}, (1, "b"): {"BP": 100_000}}
    found = exception_rows(rows, reference, {("a", "b"): 1.001}, 500_000)
    assert len(found) == 1
    assert float(found[0]["lava_r_raw"]) == 1.001
    assert found[0]["disposition"].startswith("EXCLUDED_FROM_R2_CLUMPING")


def test_candidate_clumps_do_not_use_out_of_range_edges():
    candidates = {"pair": [
        {"pair_id": "pair", "sleep_trait": "sleep", "brain_disorder": "brain",
         "SNP": "a", "CHR": 1, "BP": 10_000, "P_PLACO": 1e-10},
        {"pair_id": "pair", "sleep_trait": "sleep", "brain_disorder": "brain",
         "SNP": "b", "CHR": 1, "BP": 100_000, "P_PLACO": 1e-9},
    ]}
    reference = {(1, "a"): {"BP": 10_000, "A1": "A", "A2": "C"},
                 (1, "b"): {"BP": 100_000, "A1": "G", "A2": "T"}}
    rule = json.loads((ROOT / "brain6/config/shared_locus_rule_v1.json").read_text())
    # The raw invalid edge is intentionally absent from the valid-edge map.
    variants, loci = clump.build_outputs(candidates, reference, {("a", "a"): 1.0, ("b", "b"): 1.0}, rule)
    assert len(loci) == 1
    assert loci[0]["n_lead_signals"] == 2
    assert set(loci[0]["lead_variants"].split(";")) == {"a", "b"}
    assert {row["candidate_status"] for row in variants} == {"LEAD"}
    assert all(row["lead_SNP"] == row["SNP"] for row in variants)
