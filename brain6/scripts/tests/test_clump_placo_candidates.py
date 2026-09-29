import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from clump_placo_candidates import greedy_clumps, merge_overlapping_regions


def candidate(snp, bp, p):
    return {"SNP": snp, "CHR": 1, "BP": bp, "P_PLACO": p}


def test_greedy_clumping_uses_ld_window_and_p_order_deterministically():
    rows = [candidate("rs2", 200_000, 1e-9), candidate("rs1", 100_000, 1e-10),
            candidate("rs3", 300_000, 1e-8), candidate("rs4", 900_001, 1e-7)]
    r2 = {("rs1", "rs2"): 0.8, ("rs1", "rs3"): 0.4, ("rs2", "rs3"): 0.4}
    clumps = greedy_clumps(rows, r2, window_bp=500_000, r2_threshold=0.1)
    assert [(x["lead"]["SNP"], [m[0]["SNP"] for m in x["members"]]) for x in clumps] == [
        ("rs1", ["rs1", "rs2", "rs3"]), ("rs4", ["rs4"]) ]


def test_r2_below_threshold_does_not_clump_and_outside_window_does_not_clump():
    rows = [candidate("rs1", 100_000, 1e-10), candidate("rs2", 200_000, 1e-9),
            candidate("rs3", 700_001, 1e-8)]
    clumps = greedy_clumps(rows, {("rs1", "rs2"): 0.2, ("rs1", "rs3"): 0.99},
                           window_bp=500_000, r2_threshold=0.1)
    assert [x["lead"]["SNP"] for x in clumps] == ["rs1", "rs2", "rs3"]


def test_overlapping_regions_merge_but_keep_all_independent_leads():
    clumps = [{"chr": 1, "lead": candidate("rs1", 100_000, 1e-10), "members": []},
              {"chr": 1, "lead": candidate("rs2", 900_000, 1e-9), "members": []},
              {"chr": 2, "lead": candidate("rs3", 500_000, 1e-8), "members": []}]
    regions = merge_overlapping_regions(clumps, 500_000)
    assert len(regions) == 2
    chr1 = next(x for x in regions if x["chr"] == 1)
    assert chr1["start"] == 1 and chr1["stop"] == 1_400_000
    assert [c["lead"]["SNP"] for c in chr1["clumps"]] == ["rs1", "rs2"]
