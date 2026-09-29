import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from validate_manuscript_claims import (
    read_tsv,
    validate_citations,
    validate_claims,
    validate_local_feasibility_claims,
)


def inputs():
    manuscript = (ROOT / "brain6/paper/global_map_manuscript.md").read_text(encoding="utf-8")
    tables = [
        read_tsv(ROOT / "brain6/results/global/brain6_72_locked.tsv"),
        read_tsv(ROOT / "brain6/results/global/disorder_profile_similarity.tsv"),
        read_tsv(ROOT / "brain6/results/global/disorder_profile_leave_one_sleep_trait_out.tsv"),
        read_tsv(ROOT / "brain6/results/replication/replication_master.tsv"),
        read_tsv(ROOT / "brain6/results/supplement/table_S18_published_rg_context.tsv"),
    ]
    return manuscript, tables


def test_manuscript_key_numbers_match_frozen_result_tables():
    manuscript, tables = inputs()
    result = validate_claims(manuscript, *tables)
    assert result == {
        "global_pairs": 72,
        "original_family_significant": 35,
        "profile_leave_one_out_pairs": 2,
        "replication_rows": 1,
        "published_context_rows": 1,
    }


def test_local_feasibility_claims_match_receipt_derived_tables():
    manuscript, _ = inputs()
    base = ROOT / "brain6/results/power_optimized_sensitivity_v1"
    result = validate_local_feasibility_claims(
        manuscript,
        read_tsv(base / "canonical_family_trait_power_audit_v1.tsv"),
        read_tsv(base / "completed_trait_only_screen_comparison_v1.tsv"),
        json.loads((base / "sensitivity_family_feasibility_v1.json").read_text(encoding="utf-8")),
    )
    assert result == {
        "canonical_cells": 17465,
        "canonical_tested": 13745,
        "canonical_not_run": 3720,
        "candidate_family_not_run": 3305,
    }


def test_local_feasibility_validator_rejects_stale_family_counts():
    manuscript, _ = inputs()
    manuscript = manuscript.replace("3,720 as `NOT_RUN`", "3,719 as `NOT_RUN`")
    base = ROOT / "brain6/results/power_optimized_sensitivity_v1"
    with pytest.raises(ValueError, match="canonical local-analysis family status counts"):
        validate_local_feasibility_claims(
            manuscript,
            read_tsv(base / "canonical_family_trait_power_audit_v1.tsv"),
            read_tsv(base / "completed_trait_only_screen_comparison_v1.tsv"),
            json.loads((base / "sensitivity_family_feasibility_v1.json").read_text(encoding="utf-8")),
        )


def test_all_manuscript_citations_resolve_to_doi_bearing_bib_entries():
    manuscript, _ = inputs()
    bibliography = (ROOT / "brain6/paper/references.bib").read_text(encoding="utf-8")
    assert validate_citations(manuscript, bibliography) == {
        "cited_references": 15,
        "doi_bearing_entries": 15,
    }


def test_citation_validator_rejects_missing_and_unused_keys():
    with pytest.raises(ValueError, match="missing from references"):
        validate_citations("A claim cites [@missing2026].", "")
    with pytest.raises(ValueError, match="uncited entries"):
        validate_citations("No citations.", "@article{unused,\n doi = {10.1234/example}\n}\n")


def test_citation_validator_rejects_duplicate_keys_and_invalid_doi():
    duplicate = "@article{same,\n doi = {10.1234/a}\n}\n@article{same,\n doi = {10.1234/b}\n}\n"
    with pytest.raises(ValueError, match="duplicate citation keys"):
        validate_citations("[@same]", duplicate)
    with pytest.raises(ValueError, match="valid DOI"):
        validate_citations("[@bad_doi]", "@article{bad_doi,\n doi = {not-a-doi}\n}\n")


@pytest.mark.parametrize(
    ("old", "new", "claim"),
    [
        (
            "The corrected canonical LAVA v3 family completed its full receipt audit but failed the frozen untested-cell limit, so it is not promoted",
            "The corrected canonical LAVA v3 family is still running and remains under evaluation",
            "canonical LAVA v3 QC disposition",
        ),
        (
            "This global-map manuscript does not report local-rg, PLACO, or other downstream results as final analyses.",
            "This global-map manuscript reports local-rg and PLACO results as final analyses.",
            "downstream-analysis evidence boundary",
        ),
    ],
)
def test_manuscript_validator_rejects_stale_lava_status_or_overclaim(old, new, claim):
    manuscript, tables = inputs()
    assert old in manuscript
    manuscript = manuscript.replace(old, new)
    with pytest.raises(ValueError, match=claim):
        validate_claims(manuscript, *tables)


@pytest.mark.parametrize(
    ("old", "new", "claim"),
    [
        ("original *q* = 8.988 × 10<sup>−67</sup>", "original *q* = 8.987 × 10<sup>−67</sup>", "largest absolute global estimate"),
        ("The map contained 72 unique pairs", "The map contained 73 unique pairs", "global map totals"),
        ("FDR procedure: 9 for attention-deficit", "FDR procedure: 8 for attention-deficit", "per-disorder inherited significance counts"),
        ("Pearson correlation of 0.967", "Pearson correlation of 0.966", "abstract disorder-profile similarities"),
        ("ranged from 0.961 to 0.983", "ranged from 0.962 to 0.983", "leave-one-sleep-trait-out profile sensitivity"),
        ("*r*<sub>g</sub> = 0.3817", "*r*<sub>g</sub> = 0.3818", "insomnia–ADHD replication estimate"),
        ("*P* = 1.58 × 10<sup>−10</sup>", "*P* = 1.59 × 10<sup>−10</sup>", "published long-sleep–SCZ context"),
    ],
)
def test_manuscript_validator_rejects_stale_claim(old, new, claim):
    manuscript, tables = inputs()
    assert old in manuscript
    manuscript = manuscript.replace(old, new)
    with pytest.raises(ValueError, match=claim):
        validate_claims(manuscript, *tables)
