#!/usr/bin/env python3
"""Versioned 25-locus gate after five-track PLACO and independent B LD audit."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "brain6/results/loci/candidate_region_decisions_v2/decisions.tsv"
LOCI = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1"
INDEPENDENT = ROOT / "brain6/results/loci/admitted_B_1000g_ld_v1"
LAVA = ROOT / "brain6/results/lava_rescue_v1/lava_rescue_v1_qc.json"
OUT = ROOT / "brain6/results/loci/five_track_candidate_region_gate_v1"
OLD_REASON = "Terminal-valid protected insomnia–ADHD Track B is unavailable, so the frozen complete five-track PLACO prerequisite cannot be satisfied; canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED and no eligible independent exact-pair replication is established."
NEW_REASON = "The complete five-track PLACO candidate locus has the same lead and boundaries as this prior four-pair locus; canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED and no eligible independent exact-pair replication is established."


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def table(path: Path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream,delimiter="\t")
        return reader.fieldnames,list(reader)


def main() -> None:
    locus_receipt = json.loads((LOCI / "provenance.json").read_text())
    independent = json.loads((INDEPENDENT / "provenance.json").read_text())
    lava = json.loads(LAVA.read_text())
    if locus_receipt["status"] != "PASS_COMPLETE_PAIR_FAMILY_CANDIDATE_LOCI_ONLY" or locus_receipt["n_candidate_loci"] != 25:
        raise ValueError("five-track PLACO locus family incomplete")
    if independent["status"] != "PASS_DIAGNOSTIC_1000G_LD_AUDIT" or independent["five_track_loci_receipt_sha256"] != sha(LOCI / "provenance.json"):
        raise ValueError("B independent 1000G audit incomplete")
    if lava["classification"] != "FAIL_QC" or lava["family_complete"] or lava["promotion_permitted"]:
        raise ValueError("LAVA gate changed; full tiering must be reevaluated")
    old_fields,old_rows = table(OLD)
    locus_fields,loci = table(LOCI / "candidate_loci.tsv")
    _,b_matches = table(INDEPENDENT / "variant_matches.tsv")
    if old_fields is None or locus_fields is None or len(old_rows) != 19 or len(loci) != 25:
        raise ValueError("locus/decision row count changed")
    old_loci = {x["locus_id"] for x in old_rows}
    four_loci = {x["locus_id"] for x in loci if x["pair_id"] != "insomnia__adhd"}
    if old_loci != four_loci or len(old_loci) != 19:
        raise ValueError("prior 19 lead/boundary identities are not stable")
    match_status = {x["SNP"]:x["status"] for x in b_matches}
    if len(match_status) != 389 or independent["n_candidates"] != 389:
        raise ValueError("B independent genotype accounting changed")
    rows = []
    for row in old_rows:
        if row["Decision"] != "BLOCKED_EXTERNAL" or not row["Reason"].startswith(OLD_REASON):
            raise ValueError("prior decision content changed")
        copy = row.copy()
        copy["Reason"] = NEW_REASON + copy["Reason"][len(OLD_REASON):]
        copy["promotion_gate"] = "CANONICAL_LAVA_V3_FAILED_QC;FINAL_TIER_NOT_ASSIGNED"
        rows.append(copy)
    for locus in loci:
        if locus["pair_id"] != "insomnia__adhd":
            continue
        snps = locus["candidate_variants"].split(";")
        measured = sum(match_status.get(s) == "EXACT_ALLELE_GENOTYPE" for s in snps)
        if len(snps) != int(locus["n_candidate_variants"]):
            raise ValueError("B locus candidate count differs")
        concordance = "YES" if independent["lead_sets_identical"] else "NO"
        row = {x:"NA" for x in old_fields}
        row.update({
            "Region":f"chr{locus['CHR']}:{locus['START']}-{locus['STOP']}",
            "Pair":"insomnia__adhd","Decision":"BLOCKED_EXTERNAL",
            "Reason":"Five-track PLACO candidate locus is complete under the frozen UKB LD rule; canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED, so no final local-sharing tier or downstream mechanistic claim can be assigned.",
            "Fine-map":"NOT_RUN_GATE","Coloc":"NOT_RUN_GATE","Gene":"NOT_PRIORITIZED",
            "Tissue":"NOT_PRIORITIZED","Cell type":"NOT_PRIORITIZED","Therapeutic target":"NONE_DEFENSIBLE",
            "locus_id":locus["locus_id"],"lead_variants":locus["lead_variants"],
            "PLACO_lead_P":locus["lead_P_PLACO"],
            "LD_check":f"frozen UKB/LAVA exact coordinate and compatible alleles; independent 1000G EUR genotypes {measured}/{len(snps)}; six-lead set concordance {concordance}",
            "stable_partial_leads":"NOT_APPLICABLE_NEW_B",
            "signed_effect_relation":"NOT_ASSESSED","catalog_trait_context_rows":"NOT_ASSESSED",
            "external_replication_adjudication":"NOT_ESTABLISHED",
            "candidate_variants":str(len(snps)),"frozen_reference_eligible_variants":str(len(snps)),
            "measured_1000g_genotypes":str(measured),"excluded_no_locked_reference":"0",
            "phase3_rescued_variants":"NA","diagnostic_leads":locus["n_lead_signals"],
            "promotion_gate":"CANONICAL_LAVA_V3_FAILED_QC;FINAL_TIER_NOT_ASSIGNED",
        })
        rows.append(row)
    if len(rows) != 25 or len({x["locus_id"] for x in rows}) != 25 or {x["Decision"] for x in rows} != {"BLOCKED_EXTERNAL"}:
        raise ValueError("25-locus decision accounting failed")
    rows.sort(key=lambda x:(x["Pair"],int(x["Region"].split(":")[0][3:]),int(x["Region"].split(":")[1].split("-")[0])))
    if OUT.exists():
        raise FileExistsError("five-track region gate directory already exists")
    OUT.mkdir(parents=True,exist_ok=False)
    path = OUT / "decisions.tsv"
    with path.open("x",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=old_fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "analysis_id":"brain6-five-track-candidate-region-gate-v1",
        "status":"25_CANDIDATE_LOCI_BLOCKED_BY_LAVA",
        "prior_four_pair_region_count":19,"admitted_B_new_region_count":6,
        "five_track_candidate_locus_count":25,"blocked_count":25,
        "five_track_loci_receipt_sha256":sha(LOCI / "provenance.json"),
        "independent_B_1000G_receipt_sha256":sha(INDEPENDENT / "provenance.json"),
        "lava_rescue_qc_sha256":sha(LAVA),"prior_decisions_sha256":sha(OLD),
        "decisions_sha256":sha(path),
        "independent_B_lead_sets_identical":independent["lead_sets_identical"],
        "B_exact_1000G_genotypes":independent["n_exact_allele_genotypes"],
        "lava_not_run_cells":lava["canonical_observed_status_counts"]["NOT_RUN"],
        "lava_maximum_untested":lava["frozen_maximum_untested"],
        "interpretation":"All prior 19 PLACO loci are unchanged under complete five-track input; six protected B loci are added. No final evidence tier is assigned because the canonical LAVA family fails QC. The 1000G check is independent diagnostic evidence and does not replace the frozen UKB clumping rule."
    }
    (OUT / "provenance.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:report[k] for k in ("status","five_track_candidate_locus_count","blocked_count","independent_B_lead_sets_identical")}))


if __name__ == "__main__":
    main()
