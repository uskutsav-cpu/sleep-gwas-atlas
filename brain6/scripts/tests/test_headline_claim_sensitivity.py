import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
import build_headline_claim_sensitivity as builder


def test_interim_sensitivity_matrix_traces_evidence_and_keeps_gaps_explicit(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "OUT", tmp_path / "headline_claim_sensitivity.tsv")
    builder.main()
    with builder.OUT.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert set(rows[0]) == set(builder.FIELDS)
    by_id = {row["claim_id"]: row for row in rows}
    with (ROOT / "brain6/results/placo/placo_master.tsv").open(newline="", encoding="utf-8") as stream:
        placo_pairs = list(csv.DictReader(stream, delimiter="\t"))
    expected_pair_claims = {f"placo_{row['pair_id']}_headline_variants" for row in placo_pairs
                            if row["stage_status"].startswith("COMPLETE_QC_PASS")}
    assert {claim_id for claim_id in by_id if claim_id.startswith("placo_")} == \
        expected_pair_claims | {"placo_partial_candidate_regions"}
    assert "35 retain significance" in by_id["global_map_35_of_72"]["manuscript_claim"]
    assert "0/12 significant" in by_id["alzheimer_no_inherited_signals"]["manuscript_claim"]
    assert by_id["profile_adhd_mdd"]["profile_leave_one_sleep_trait_out"] == \
        "12 omissions; Pearson r range 0.961–0.983; descriptive only"
    assert by_id["profile_scz_bipolar"]["profile_leave_one_sleep_trait_out"] == \
        "12 omissions; Pearson r range 0.830–0.917; descriptive only"
    assert all("LOO_RANGES_REPORTED" in by_id[claim_id]["overall_sensitivity_status"]
               for claim_id in ("profile_adhd_mdd", "profile_scz_bipolar"))
    assert by_id["insomnia_adhd_directional_replication"]["independent_replication"] == \
        "PARTIAL_LEVEL_B_DIRECTIONAL_ONLY"
    context = by_id["longsleep_scz_published_context"]
    assert context["evidence_status"] == "SUPPORTED_PUBLISHED_CONTEXT_ONLY_LEVEL_D"
    assert context["independent_replication"].startswith("LEVEL_D_CONTEXT_ONLY;")
    assert "NOT_REPLICATION" in context["overall_sensitivity_status"]
    assert by_id["placo_partial_candidate_regions"]["overall_sensitivity_status"] == \
        "PRELIMINARY_CANDIDATE_GROUPING_ONLY; MHC_NO_CANDIDATES_IN_INTERVAL_FOR_OBSERVED_PARTIAL_SET; NOT_FINAL_LOCI"
    assert by_id["placo_partial_candidate_regions"]["MHC_exclusion"].startswith(
        "NO_CANDIDATES_IN_INTERVAL; 0/2273 position-validated candidates")
    assert "across 4 published QC-passed PLACO pairs" in by_id["placo_partial_candidate_regions"]["manuscript_claim"]
    assert "protected Track B is absent" in by_id["placo_partial_candidate_regions"]["claim_scope"]
    assert by_id["placo_longsleep__scz_headline_variants"]["sample_overlap_sensitivity"].startswith("UNKNOWN;")
    bipolar_overlap = by_id["placo_longsleep__bipolar_headline_variants"]["sample_overlap_sensitivity"]
    assert bipolar_overlap.startswith("DOCUMENTED;")
    assert "exact participant-level overlap" in bipolar_overlap
    pair_sensitivity = [row for claim_id, row in by_id.items()
                        if claim_id.startswith("placo_") and claim_id.endswith("_headline_variants")]
    assert len(pair_sensitivity) == 4
    assert all(row["MHC_exclusion"].startswith("NO_CANDIDATES_IN_INTERVAL;") for row in pair_sensitivity)
    assert all("MHC_NO_CANDIDATES_IN_INTERVAL_FOR_PARTIAL_CANDIDATE_SET" in row["overall_sensitivity_status"]
               for row in pair_sensitivity)
    assert all(row["overall_sensitivity_status"] != "PASS" for row in rows)
    for row in rows:
        for record in row["evidence_sources"].split(";"):
            relative, expected = record.split("@sha256:", 1)
            path = ROOT / relative
            assert path.is_file() and builder.digest(path) == expected
