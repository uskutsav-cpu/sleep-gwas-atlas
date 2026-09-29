#!/usr/bin/env python3
"""Compose the verified four-track QC and newly admitted protected B into five tracks."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/track_b_admission_v1"
PLACOROOT = ROOT / "brain6/results/placo"
MASTER4 = PLACOROOT / "placo_master.tsv"
QC4 = PLACOROOT / "placo_v3_pair_qc_validation.json"
SLOT = PLACOROOT / "insomnia__adhd"
FAMILY = ROOT / "brain6/config/placo_family_v3/family_lock.json"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def require(value: bool, reason: str) -> None:
    if not value:
        raise ValueError(reason)


def write_once(path: Path, value: str) -> None:
    if path.exists() and path.read_text() != value:
        raise FileExistsError(f"refusing to overwrite different five-track output: {path}")
    if not path.exists():
        path.write_text(value)


def main() -> None:
    criteria = read(OUT / "criteria.lock.json")
    adapter = read(OUT / "adapter_validation.json")
    receipt_path = SLOT / "admission_receipt.json"
    receipt = read(receipt_path)
    four = read(QC4)
    family = read(FAMILY)
    require(receipt["status"] == "PROMOTED_AFTER_ALL_CHECKS_PASS", "protected B is not admitted")
    require(receipt["criteria_sha256"] == sha(OUT / "criteria.lock.json") and
            receipt["adapter_validation_sha256"] == sha(OUT / "adapter_validation.json"), "admission identity drifted")
    require(sha(FAMILY) == receipt["family_lock_sha256"] == four["family_lock_sha256"] == criteria["expected"]["family_lock_sha256"], "frozen family lock drifted")
    require(family["n_selected_tracks"] == 5 and family["protected_legacy_pair"] == "insomnia__adhd" and
            family["headline_threshold"] == criteria["frozen_statistics"]["headline_p_strictly_less_than"], "five-track policy drifted")
    for name, info in receipt["files"].items():
        path = SLOT / name
        require(path.stat().st_size == info["bytes"] and sha(path) == info["sha256"], f"admitted B file drifted: {name}")
    require(receipt["files"]["variants.tsv.gz"]["sha256"] == adapter["staged_view_sha256"] and
            receipt["files"]["B.full.tsv.gz"]["sha256"] == criteria["expected"]["ledger_sha256"], "protected B view or legacy source drifted")
    require(adapter["status"] == "PASS" and adapter["rows"] == adapter["independent_bh_checked_rows"] == criteria["expected"]["rows"] and
            adapter["independent_bh_q_mismatches"] == adapter["five_track_adjustment_mismatches"] == 0, "B multiplicity audit failed")
    require(four["status"] == "PASS_FOUR_PAIR_QC; FIVE_TRACK_FAMILY_INCOMPLETE" and
            four["n_required_tracks"] == 5 and four["full_family_complete"] is False and
            len(four["published_v3_pairs"]) == len(four["published_output_sha256"]) == len(four["pairs"]) == 4,
            "prior four-track QC cannot be composed")
    with MASTER4.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    require(fields is not None and len(rows) == 4 and {r["pair_id"] for r in rows} == set(four["published_v3_pairs"]), "prior four-track master changed")
    pair_qc = {x["pair_id"]:x for x in four["pairs"]}
    for row in rows:
        pair = row["pair_id"]
        require(pair in family["pairs"] and row["output_sha256"] == four["published_output_sha256"][pair]
                and sha(ROOT / row["output_path"]) == row["output_sha256"], f"frozen pair output drifted: {pair}")
        require(row["family_lock_sha256"] == receipt["family_lock_sha256"] and
                row["stage_status"].startswith("COMPLETE_QC_PASS") and
                int(row["total_rows"]) == pair_qc[pair]["total_rows"] and
                int(row["tested_rows"]) == pair_qc[pair]["tested"] and
                int(row["headline_variants"]) == pair_qc[pair]["headline_variants"] and
                int(row["excluded_extreme_z"]) == pair_qc[pair]["excluded_extreme_z"] and
                int(row["numerical_failures"]) == pair_qc[pair]["numerical_failures"] and
                pair_qc[pair]["adjustment_tracks"] == 5 and
                float(row["across_track_p_threshold"]) == family["headline_threshold"], f"frozen pair QC drifted: {pair}")
    require(four["q_values_independently_verified"] == sum(x["tested"] for x in four["pairs"]), "prior BH audit denominator changed")
    b = {
        "pair_id":"insomnia__adhd","sleep_trait":"insomnia","brain_disorder":"adhd",
        "stage_status":"COMPLETE_QC_PASS_LEGACY_ADMISSION_V1_LOCUS_DEDUP_PENDING",
        "total_rows":str(adapter["rows"]),"tested_rows":str(adapter["rows"]),
        "excluded_extreme_z":"0","numerical_failures":"0","failure_rate":"0.0",
        "across_track_p_threshold":str(family["headline_threshold"]),
        "headline_variants":str(adapter["headline_variants"]),
        "independent_loci_status":"NOT_YET_ESTIMATED",
        "output_path":"brain6/results/placo/insomnia__adhd/variants.tsv.gz",
        "output_sha256":receipt["files"]["variants.tsv.gz"]["sha256"],
        "source_artifact_fingerprint":receipt["historical_run_fingerprint"],
        "source_receipt_sha256":sha(receipt_path),
        "family_lock_sha256":receipt["family_lock_sha256"],
    }
    require(set(b) == set(fields), "five-track master schema differs")
    all_rows = sorted(rows + [b], key=lambda x:x["pair_id"])
    require(len(all_rows) == 5 and len({r["pair_id"] for r in all_rows}) == 5, "five-track set is incomplete")
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(all_rows)
    master_path = OUT / "placo_master_five_track_v1.tsv"
    write_once(master_path, output.getvalue())
    report = {
        "analysis_id":"brain6-placo-five-track-after-legacy-B-admission-v1", "status":"PASS_FIVE_PAIR_QC_LOCUS_DEDUP_PENDING",
        "family_lock_sha256":sha(FAMILY), "protected_B_admission_receipt_sha256":sha(receipt_path),
        "prior_four_pair_qc_receipt_sha256":sha(QC4),"prior_four_pair_master_sha256":sha(MASTER4),
        "five_pair_master_sha256":sha(master_path),
        "pairs":[r["pair_id"] for r in all_rows],
        "output_sha256":{r["pair_id"]:r["output_sha256"] for r in all_rows},
        "rows_total":sum(int(r["total_rows"]) for r in all_rows),
        "tested_rows_total":four["q_values_independently_verified"] + adapter["independent_bh_checked_rows"],
        "q_values_independently_verified_in_original_four_pair_audit":four["q_values_independently_verified"],
        "q_values_independently_verified_in_B_adapter":adapter["independent_bh_checked_rows"],
        "five_track_adjustment_mismatches_B":adapter["five_track_adjustment_mismatches"],
        "headline_variants_by_pair":{r["pair_id"]:int(r["headline_variants"]) for r in all_rows},
        "within_pair_q_below_0_05_B":adapter["within_pair_q_below_0_05"],
        "independent_loci_status":"NOT_YET_ESTIMATED",
        "candidate_region_decision_status":"BLOCKED_BY_LAVA_AND_LD_QC",
        "audit_method":"Composition of the hash-pinned prior independent four-pair BH and five-track Bonferroni audit with the new full-row B independent BH and five-track Bonferroni audit. Each pair correction is independent and uses the same frozen multiplier 5.",
    }
    report_path = OUT / "five_track_family_audit.json"
    write_once(report_path,json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":report["status"],"pairs":report["pairs"],"tested_rows_total":report["tested_rows_total"],"headline_variants_by_pair":report["headline_variants_by_pair"]}))


if __name__ == "__main__":
    main()
