from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "build_partial_placo_signed_ld.py"
SPEC = importlib.util.spec_from_file_location("partial_placo_signed_ld", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(module)


def row(snp: str, lead: str, r2: str) -> dict[str, str]:
    return {
        "pair_id": "pair", "SNP": snp, "CHR": "1", "BP": "100",
        "candidate_status": "LEAD" if snp == lead else "LD_CLUMPED",
        "lead_SNP": lead, "reference_status": "EXACT_MATCH",
        "r2_to_lead": r2, "locus_id": "locus1",
    }


def test_signed_ld_preserves_sign_and_matches_frozen_r2():
    refs = {
        (1, "rs1"): {"A1": "A", "A2": "G"},
        (1, "rs2"): {"A1": "C", "A2": "T"},
    }
    rows = module.make_signed_rows(
        [row("rs1", "rs1", "1"), row("rs2", "rs1", "0.64")],
        refs,
        {("rs1", "rs1"): 1.0, ("rs1", "rs2"): -0.8},
    )
    assert rows[1]["signed_r_ref_A1"] == "-0.8"
    assert rows[1]["r2_ref"] == "0.64"
    assert rows[1]["ld_qc_status"] == "PASS_VALID_SIGNED_R"
    assert rows[1]["A1_REF"] == "C"
    assert rows[1]["LEAD_A1_REF"] == "A"


def test_exact_match_with_missing_ld_is_rejected():
    refs = {(1, "rs1"): {"A1": "A", "A2": "G"},
            (1, "rs2"): {"A1": "C", "A2": "T"}}
    try:
        module.make_signed_rows([row("rs2", "rs1", "0.5")], refs, {})
    except ValueError as error:
        assert "Signed LD is absent" in str(error)
    else:
        raise AssertionError("Missing signed LD must not be silently treated as zero")


def test_out_of_range_lava_value_is_not_reported_as_a_correlation():
    refs = {(1, "rs1"): {"A1": "A", "A2": "G"},
            (1, "rs2"): {"A1": "C", "A2": "T"}}
    out = module.make_signed_rows([row("rs2", "rs1", "1.0006")], refs,
                                  {("rs1", "rs2"): 1.0003})
    assert out[0]["lava_r_raw"] == "1.0003"
    assert out[0]["signed_r_ref_A1"] == "NA"
    assert out[0]["r2_ref"] == "NA"
    assert out[0]["ld_qc_status"] == "OUT_OF_RANGE_RAW_LAVA_R"
    assert out[0]["r2_matches_existing_clump"] == "RAW_R_SQUARED_MATCHES_CLUMP_ONLY"


def test_unmatched_candidate_is_retained_without_invented_ld():
    candidate = row("rs3", "NA", "NA")
    candidate["reference_status"] = "NO_REFERENCE_VARIANT"
    candidate["candidate_status"] = "NOT_CLUMPED_NO_EXACT_REFERENCE_MATCH"
    out = module.make_signed_rows([candidate], {}, {})
    assert out[0]["signed_r_ref_A1"] == "NA"
    assert out[0]["reference_status"] == "NO_REFERENCE_VARIANT"
