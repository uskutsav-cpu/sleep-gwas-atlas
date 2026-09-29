#!/usr/bin/env python3
"""Bind the post-admission Brain6 PLACO, locus, and LAVA gate states."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/brain6_unblock_parallel_v2"
ADMISSION = ROOT / "brain6/results/placo/insomnia__adhd/admission_receipt.json"
FIVE = ROOT / "brain6/results/track_b_admission_v1/five_track_family_audit.json"
LOCI = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/provenance.json"
LD = ROOT / "brain6/results/loci/admitted_B_1000g_ld_v1/provenance.json"
GATE = ROOT / "brain6/results/loci/five_track_candidate_region_gate_v1/provenance.json"
DECISIONS = ROOT / "brain6/results/loci/five_track_candidate_region_gate_v1/decisions.tsv"
LAVA = ROOT / "brain6/results/lava_rescue_v1/lava_rescue_v1_qc.json"
FEASIBILITY = ROOT / "brain6/results/power_optimized_sensitivity_v1/sensitivity_family_feasibility_v1.json"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    paths = (ADMISSION,FIVE,LOCI,LD,GATE,DECISIONS,LAVA,FEASIBILITY)
    admission,five,loci,ld,gate,lava,feasibility = (json.loads(p.read_text()) for p in (ADMISSION,FIVE,LOCI,LD,GATE,LAVA,FEASIBILITY))
    with DECISIONS.open(newline="") as stream:
        decisions = list(csv.DictReader(stream,delimiter="\t"))
    if admission["status"] != "PROMOTED_AFTER_ALL_CHECKS_PASS" or five["status"] != "PASS_FIVE_PAIR_QC_LOCUS_DEDUP_PENDING":
        raise ValueError("Track B/five-track QC has not passed")
    if loci["status"] != "PASS_COMPLETE_PAIR_FAMILY_CANDIDATE_LOCI_ONLY" or loci["n_candidate_loci"] != 25:
        raise ValueError("five-track PLACO loci incomplete")
    if ld["status"] != "PASS_DIAGNOSTIC_1000G_LD_AUDIT" or not ld["lead_sets_identical"]:
        raise ValueError("independent B lead evidence changed")
    if gate["status"] != "25_CANDIDATE_LOCI_BLOCKED_BY_LAVA" or gate["decisions_sha256"] != sha(DECISIONS):
        raise ValueError("region gate changed")
    if len(decisions) != 25 or {r["Decision"] for r in decisions} != {"BLOCKED_EXTERNAL"}:
        raise ValueError("25-row decision classification changed")
    if lava["classification"] != "FAIL_QC" or lava["family_complete"] or lava["promotion_permitted"]:
        raise ValueError("LAVA classification changed")
    if feasibility["status"] != "PASS_SLEEP_ONLY_REPLACEMENTS_CANNOT_PASS_FROZEN_FAMILY_QC":
        raise ValueError("sleep-only feasibility changed")
    hashes = {str(p.relative_to(ROOT)):sha(p) for p in paths}
    receipt = {
        "analysis_id":"brain6-parallel-unblock-v2-post-B-admission",
        "track_B":"PROMOTED_AFTER_ALL_CHECKS_PASS",
        "five_track_PLACO":"PASS",
        "five_track_candidate_loci":25,
        "original_candidate_regions":19,
        "new_B_candidate_loci":6,
        "independent_B_1000G_genotypes":ld["n_exact_allele_genotypes"],
        "independent_B_1000G_lead_concordance":ld["lead_sets_identical"],
        "LAVA_family":"FAIL_QC",
        "LAVA_not_run":lava["canonical_observed_status_counts"]["NOT_RUN"],
        "LAVA_maximum_not_run":lava["frozen_maximum_untested"],
        "five_retained_disorder_not_run_even_if_both_sleep_traits_perfect":feasibility["best_case_lower_bounds"]["both_sleep_traits_zero_not_run"]["retained_five_disorder_traits_not_run"],
        "region_tiering":"BLOCKED",
        "blocked_candidate_regions":len(decisions),
        "source_sha256":hashes,
    }
    report = f"""# Brain6 combined unblock report v2

## Track B and five-track PLACO

The exact historical insomnia–ADHD Track B package passed the authorized legacy validator, full 5,514,399-row adapter/BH audit, and scientific compatibility check before promotion. Its original 17-field ledger and provenance remain beside the checked ten-field view in the protected slot. The five-track pair QC covers **20,370,452 tested rows** and passes the frozen five-track correction. No earlier four-track, v2, roundoff, or canonical LAVA output was overwritten.

The frozen UKB/LAVA LD rule gives **25 complete-family PLACO candidate loci**: the prior four pairs retain their exact **19** lead/boundary identities, and admitted B adds **six**. All **389** B headline candidates match the locked UKB reference by coordinate and compatible allele. An independent 503-sample 1000 Genomes EUR diagnostic measures exact allele genotypes for **388/389** B candidates and reproduces the same six lead identities; the single unmatched nonlead variant is retained in the audit. These are shared-association candidate loci, not final local-sharing or causal claims.

## LAVA rescue and final gates

Canonical LAVA v3 remains `FAILED_QC_NOT_PROMOTED`. **{receipt['LAVA_not_run']}/17,465** univariate cells are `NOT_RUN`, above the frozen maximum **{receipt['LAVA_maximum_not_run']}**. The rescue-v1 pilot found no broad technical correction; the full 12,475-slot pair-locus family was not launched or promoted. Low local-h² accounts for 3,564 of the 3,720 untested canonical cells. Even hypothetically eliminating every NOT_RUN cell in both sleep traits leaves **{receipt['five_retained_disorder_not_run_even_if_both_sleep_traits_perfect']}** untested disorder cells, still above the ceiling. The observed continuous-duration long-sleep sensitivity screen therefore cannot rescue the frozen seven-trait family by itself.

| Gate | State |
|---|---|
| Protected Track B | PASS, promoted under the approved narrow amendment |
| Five-track PLACO pair family and candidate loci | PASS; 25 candidate loci |
| Complete LAVA family | FAIL_QC |
| Final region tiering | BLOCKED; all 25 candidate loci remain `BLOCKED_EXTERNAL` |

The 19 original candidate regions are unchanged as PLACO loci; [the versioned 25-row decision table](../loci/five_track_candidate_region_gate_v1/decisions.tsv) includes the six new B regions. No final tier, fine-mapping, colocalization, gene, tissue, cell-type, regulatory, or therapeutic claim is promoted from the failed LAVA family.

## Exact next requirement

Identify and source-verify a scientifically equivalent, better-powered **multi-trait** GWAS strategy that can plausibly bring the frozen seven-trait family to at most 873 untested cells, then predeclare and validate a separate pilot before any full rescue run. A sleep-only substitution cannot meet that ceiling. If the intended scientific question instead requires a different family/QC framework, obtain explicit authorization for a prospective amendment; do not relabel or overwrite canonical v3.
"""
    OUT.mkdir(parents=True,exist_ok=False)
    (OUT / "combined_unblock_report.md").write_text(report)
    receipt["report_sha256"] = sha(OUT / "combined_unblock_report.md")
    (OUT / "terminal_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:receipt[k] for k in ("five_track_PLACO","five_track_candidate_loci","LAVA_family","region_tiering","blocked_candidate_regions")}))


if __name__ == "__main__":
    main()
