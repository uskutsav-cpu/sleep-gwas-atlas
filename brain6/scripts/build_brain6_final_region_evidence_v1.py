#!/usr/bin/env python3
"""Join frozen candidate geography, canonical LAVA gate, and exploratory lookups."""

import csv
import json
import hashlib
from pathlib import Path


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    repo = Path(__file__).resolve().parents[2]
    core = repo / "brain6/results/loci"
    source = repo / "brain6/results/lava_longsleep_source_rescue_v1"
    paths = {
        "geography": core / "five_track_cross_pair_region_reconciliation_v1/geographic_region_groups.tsv",
        "candidates": core / "five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv",
        "lava": core / "five_track_candidate_lava_gate_v1/candidate_lava_gate.tsv",
        "replication": source / "exploratory_long_sleep_replication_v1/pair_candidate_replication.tsv",
    }
    regions = {r["region_group"]: r for r in rows(paths["geography"])}
    candidates = rows(paths["candidates"])
    lava = {r["candidate_locus_id"]: r for r in rows(paths["lava"])}
    replication = {r["candidate_locus_id"]: r for r in rows(paths["replication"])}
    assert len(regions) == 20 and len(candidates) == 25
    output = []
    for c in candidates:
        l = lava[c["candidate_locus_id"]]
        e = replication[c["candidate_locus_id"]]
        g = regions[c["region_group"]]
        output.append({
            "geographic_region": c["region_group"], "chromosome": g["chromosome"],
            "start_grch37": g["start"], "stop_grch37": g["stop"],
            "pair_id": c["pair_id"], "candidate_locus_id": c["candidate_locus_id"],
            "lead_variants": c["lead_variants"], "placo_lead_p": c["lead_p_placo"],
            "canonical_lava_locus_id": l["lava_locus_id"],
            "canonical_trait1_status": l["trait1_status"],
            "canonical_trait2_status": l["trait2_status"],
            "canonical_both_strict_univariate_gates_pass": l["both_strict_gates_pass"],
            "confirmatory_family_status": "FAILED_QC_NOT_PROMOTED",
            "confirmatory_promotion": "NO",
            "exploratory_source": "Austin-Zimmerman_2023_UKB_MVP_EUR_ge10_vs_7to8",
            "exploratory_classification": e["classification"],
            "exploratory_best_exact_lead": e["best_exact_lead"],
            "exploratory_best_exact_p": e["best_exact_p"],
            "exploratory_best_exact_direction": e["best_exact_direction"],
            "exploratory_min_region_p": e["minimum_region_p"],
            "independent_replication": "NO_UKB_OVERLAP",
        })
    assert len(output) == 25 and len({r["geographic_region"] for r in output}) == 20
    target = source / "BRAIN6_FINAL_REGION_EVIDENCE.tsv"
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    receipt = {"analysis_id": "brain6_final_region_evidence_v1", "input_sha256": {k: sha(p) for k, p in paths.items()},
               "output_sha256": sha(target), "geographic_regions": 20, "pair_specific_candidates": 25,
               "confirmatory_status": "FAILED_QC_NOT_PROMOTED"}
    (source / "BRAIN6_FINAL_REGION_EVIDENCE.provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
