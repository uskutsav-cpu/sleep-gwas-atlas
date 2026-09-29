#!/usr/bin/env python3
"""Refresh the 19 carried-forward candidate gate after protected B admission."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "brain6/results/loci/candidate_region_decisions_v2/decisions.tsv"
DEST = ROOT / "brain6/results/loci/candidate_region_decisions_v3_post_b_admission"
AUDIT = ROOT / "brain6/results/track_b_admission_v1/five_track_family_audit.json"
LAVA = ROOT / "brain6/results/lava_rescue_v1/lava_rescue_v1_qc.json"
OLD = "Terminal-valid protected insomnia–ADHD Track B is unavailable, so the frozen complete five-track PLACO prerequisite cannot be satisfied; canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED and no eligible independent exact-pair replication is established."
NEW = "Protected insomnia–ADHD Track B and the five-track PLACO pair family now pass QC, while full-family locus deduplication remains pending; canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED and no eligible independent exact-pair replication is established."


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    a, lava = json.loads(AUDIT.read_text()), json.loads(LAVA.read_text())
    if a["status"] != "PASS_FIVE_PAIR_QC_LOCUS_DEDUP_PENDING" or len(a["pairs"]) != 5:
        raise ValueError("five-track PLACO QC is not complete")
    if lava["classification"] != "FAIL_QC" or lava["family_complete"] or lava["promotion_permitted"]:
        raise ValueError("LAVA gate state changed; reevaluate rather than copy decisions")
    if not lava["canonical_observed_status_counts"]["NOT_RUN"] > lava["frozen_maximum_untested"]:
        raise ValueError("LAVA failure count changed")
    with SOURCE.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if fields is None or len(rows) != 19 or {r["Decision"] for r in rows} != {"BLOCKED_EXTERNAL"}:
        raise ValueError("frozen 19-row candidate state changed")
    if {r["promotion_gate"] for r in rows} != {"FIVE_TRACK_PLACO_INCOMPLETE;CANONICAL_LAVA_V3_FAILED_QC"}:
        raise ValueError("candidate gate text changed")
    for row in rows:
        if not row["Reason"].startswith(OLD):
            raise ValueError(f"candidate rationale changed: {row['Region']}")
        row["Reason"] = NEW + row["Reason"][len(OLD):]
        row["promotion_gate"] = "FIVE_TRACK_PLACO_PAIR_QC_PASS_LOCUS_DEDUP_PENDING;CANONICAL_LAVA_V3_FAILED_QC"
    DEST.mkdir(parents=True, exist_ok=False)
    result_path = DEST / "decisions.tsv"
    with result_path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "analysis_id":"brain6-19-region-gate-refresh-after-protected-B-admission-v1",
        "status":"19_CARRIED_FORWARD_CANDIDATES_BLOCKED",
        "source_decisions_sha256":sha(SOURCE),"five_track_audit_sha256":sha(AUDIT),
        "lava_rescue_qc_sha256":sha(LAVA),"result_sha256":sha(result_path),
        "rows":len(rows),"blocked_rows":sum(r["Decision"] == "BLOCKED_EXTERNAL" for r in rows),
        "placo_pair_family_qc":"PASS_FIVE_PAIR_QC_LOCUS_DEDUP_PENDING",
        "placo_full_family_locus_status":"NOT_YET_ESTIMATED",
        "lava_family":"FAIL_QC","lava_not_run_cells":lava["canonical_observed_status_counts"]["NOT_RUN"],
        "lava_allowed_not_run_cells":lava["frozen_maximum_untested"],
        "explanation":"This versioned gate refresh carries forward the 19 prior partial-family candidate regions without recomputing or claiming new loci. The B admission clears the pair-family PLACO QC deficit, but exact-reference locus deduplication and a complete admissible LAVA family remain unavailable. No final evidence tier is assigned."
    }
    (DEST / "provenance.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":receipt["status"],"blocked_rows":receipt["blocked_rows"],"lava_not_run_cells":receipt["lava_not_run_cells"]}))


if __name__ == "__main__":
    main()
