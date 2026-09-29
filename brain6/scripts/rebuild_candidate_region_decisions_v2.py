#!/usr/bin/env python3
"""Bind 19 partial candidate decisions to complete eligible 1000G LD checks."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "brain6/results/loci/candidate_region_decisions_v1"
LD = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v2"
PLACO_QC = ROOT / "brain6/results/placo/placo_v3_pair_qc_validation.json"
LAVA = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json"
TRACK_B = ROOT / "results/track_b/pleiotropy/results/placo/B.full.tsv.gz"
OUT = ROOT / "brain6/results/loci/candidate_region_decisions_v2"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    old_receipt = json.loads((OLD / "provenance.json").read_text())
    ld_receipt = json.loads((LD / "provenance.json").read_text())
    for name, expected in old_receipt["outputs_sha256"].items():
        assert sha(OLD / name) == expected
    for name, expected in ld_receipt["outputs_sha256"].items():
        assert sha(LD / name) == expected
    qc = json.loads(PLACO_QC.read_text())
    lava = json.loads(LAVA.read_text())
    assert qc["full_family_complete"] is False and qc["protected_track_b_published"] is False
    assert lava["overall_status"] == "FAILED_QC_NOT_PROMOTED" and lava["promotion_permitted"] is False
    assert not TRACK_B.exists()
    prior = tsv(OLD / "decisions.tsv")
    regions = tsv(LD / "diagnostic_regions.tsv")
    assignments = tsv(LD / "variant_assignments.tsv")
    assert len(prior) == len(regions) == 19 and len(assignments) == 2297
    region_map = {(r["pair_id"], r["CHR"], int(r["START"]), int(r["STOP"])): r for r in regions}
    assert len(region_map) == 19
    groups: dict[tuple[str, str, int, int], list[dict[str, str]]] = defaultdict(list)
    for variant in assignments:
        matches = [key for key in region_map if key[0] == variant["pair_id"] and key[1] == variant["CHR"]
                   and key[2] <= int(variant["BP"]) <= key[3]]
        assert len(matches) == 1, f"candidate has {len(matches)} regions: {variant['SNP']}"
        groups[matches[0]].append(variant)
    assert sum(map(len, groups.values())) == 2297
    output = []
    for old in prior:
        chrom, bounds = old["Region"].split(":", 1)
        start, stop = map(int, bounds.split("-", 1))
        key = (old["Pair"], chrom.removeprefix("chr"), start, stop)
        region = region_map[key]
        variants = groups[key]
        status = Counter(x["status"] for x in variants)
        exact = status["LEAD"] + status["LD_CLUMPED"]
        excluded = status["NOT_CLUMPED_NO_EXACT_LOCKED_REFERENCE"]
        rescued = sum(x["reference"] == "PHASE3_GRCH37_CANDIDATE_TO_30X_GRCH38_LEAD" for x in variants)
        assert exact + excluded == len(variants)
        assert old["Decision"] == "BLOCKED_EXTERNAL"
        assert old["Pair"] == region["pair_id"]
        assert set(old["lead_variants"].split(";")) == set(region["lead_variants"].split(";"))
        row = dict(old)
        row["LD_check"] = (f"partial-family 1000G EUR exact genotypes {exact}/{exact} frozen-reference-eligible; "
                           f"Phase3-to-30x rescues {rescued}; excluded without exact frozen UKB reference {excluded}; "
                           f"diagnostic leads {region['n_leads']}; region unchanged")
        row.update(candidate_variants=len(variants), frozen_reference_eligible_variants=exact,
                   measured_1000g_genotypes=exact, excluded_no_locked_reference=excluded,
                   phase3_rescued_variants=rescued, diagnostic_leads=region["n_leads"],
                   promotion_gate="FIVE_TRACK_PLACO_INCOMPLETE;CANONICAL_LAVA_V3_FAILED_QC")
        output.append(row)
    assert len(output) == 19
    assert sum(int(x["measured_1000g_genotypes"]) for x in output) == 2273
    assert sum(int(x["excluded_no_locked_reference"]) for x in output) == 24
    assert sum(int(x["phase3_rescued_variants"]) for x in output) == 8
    OUT.mkdir(parents=True, exist_ok=False)
    table = OUT / "decisions.tsv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    inputs = [OLD / "decisions.tsv", OLD / "provenance.json", LD / "variant_assignments.tsv",
              LD / "diagnostic_regions.tsv", LD / "provenance.json", PLACO_QC, LAVA, Path(__file__)]
    receipt = {"analysis_id": "brain6_candidate_region_decisions_v2",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "scope": "Four-pair partial candidates; complete eligible genotype diagnostics but no full-family PLACO or admissible local LAVA decision",
               "status_counts": dict(Counter(x["Decision"] for x in output)),
               "genotype_counts": {"candidate_variants": 2297, "measured_eligible": 2273,
                                   "excluded_no_locked_reference": 24, "phase3_rescued": 8},
               "frozen_result_mutation": False,
               "inputs_sha256": {str(p): sha(p) for p in inputs},
               "outputs_sha256": {table.name: sha(table)}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"regions": len(output), "status_counts": receipt["status_counts"],
                      "genotype_counts": receipt["genotype_counts"]}))


if __name__ == "__main__":
    main()
