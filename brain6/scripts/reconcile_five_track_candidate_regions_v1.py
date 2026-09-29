#!/usr/bin/env python3
"""Group overlapping five-track PLACO candidate intervals across pairs.

The frozen per-pair clumps remain unchanged. This cross-pair region grouping is
descriptive and cannot promote a candidate through the failed LAVA family.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1"
RULE = ROOT / "brain6/config/shared_locus_rule_v1.json"
OUT = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1"


def sha(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path: Path, records: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite reconciliation: {OUT}")
    source_receipt = json.loads((SOURCE / "provenance.json").read_text())
    if source_receipt["status"] != "PASS_COMPLETE_PAIR_FAMILY_CANDIDATE_LOCI_ONLY":
        raise ValueError("Five-track candidate source is not admitted")
    loci_path, variants_path = SOURCE / "candidate_loci.tsv", SOURCE / "candidate_variants.tsv"
    if (sha(loci_path) != source_receipt["candidate_loci_sha256"]
            or sha(variants_path) != source_receipt["candidate_variants_sha256"]):
        raise ValueError("Five-track candidate files differ from their source receipt")
    rule = json.loads(RULE.read_text())
    if (rule["clumping"]["r2_threshold"] != 0.1
            or rule["region_definition"]["lead_flank_kb"] != 500):
        raise ValueError("Frozen source locus rule changed")
    loci = rows(loci_path)
    if len(loci) != source_receipt["n_candidate_loci"] or len(loci) != 25:
        raise ValueError("Expected exactly 25 five-track pair-specific candidates")
    by_id = {}
    for row in loci:
        key = row["locus_id"]
        if key in by_id:
            raise ValueError(f"Duplicate candidate locus: {key}")
        if int(row["START"]) > int(row["STOP"]):
            raise ValueError(f"Invalid candidate interval: {key}")
        by_id[key] = row
    variant_rows = rows(variants_path)
    if len(variant_rows) != source_receipt["n_candidate_variants"]:
        raise ValueError("Candidate-variant count differs from source receipt")
    variants = defaultdict(list)
    unassigned_no_reference = 0
    for row in variant_rows:
        if row["locus_id"] not in by_id:
            if (row["locus_id"] == "NA"
                    and row["candidate_status"] == "NOT_CLUMPED_NO_EXACT_REFERENCE_MATCH"):
                unassigned_no_reference += 1
                continue
            raise ValueError("Unexpected candidate variant without a pair-specific locus")
        variants[row["locus_id"]].append(row)
    if len(variant_rows) - unassigned_no_reference != source_receipt["n_exact_reference_matched_variants"]:
        raise ValueError("Exact-reference variant accounting differs from source receipt")

    def lead_ld(a: dict, b: dict) -> list[str]:
        observations = set()
        for lead_a in a["lead_variants"].split(";"):
            for lead_b in b["lead_variants"].split(";"):
                if lead_a == lead_b:
                    observations.add(f"{lead_a}|{lead_b}|1")
                for source, snp, assigned in ((a, lead_b, lead_a), (b, lead_a, lead_b)):
                    for row in variants[source["locus_id"]]:
                        if (row["SNP"] == snp and row["lead_SNP"] == assigned
                                and row["candidate_status"] in {"LEAD", "LD_CLUMPED"}):
                            r2 = float(row["r2_to_lead"])
                            if not 0 <= r2 <= 1.000001:
                                raise ValueError("Invalid frozen UKB LD value")
                            observations.add(f"{lead_a}|{lead_b}|{r2:.8g}")
        return sorted(observations)

    parents = {key: key for key in by_id}

    def find(key: str) -> str:
        while parents[key] != key:
            parents[key] = parents[parents[key]]
            key = parents[key]
        return key

    edges = []
    for a, b in combinations(loci, 2):
        if a["CHR"] != b["CHR"]:
            continue
        overlap = min(int(a["STOP"]), int(b["STOP"])) - max(int(a["START"]), int(b["START"])) + 1
        if overlap <= 0:
            continue
        if a["pair_id"] == b["pair_id"]:
            raise ValueError("Within-pair overlapping regions survived the frozen pair-specific merge")
        root_a, root_b = find(a["locus_id"]), find(b["locus_id"])
        if root_a != root_b:
            parents[root_b] = root_a
        edges.append({"locus_a": a["locus_id"], "locus_b": b["locus_id"],
                      "pair_a": a["pair_id"], "pair_b": b["pair_id"],
                      "chromosome": a["CHR"], "interval_overlap_bp": overlap,
                      "observed_frozen_ukb_lead_r2": ";".join(lead_ld(a, b)),
                      "grouping_basis": "OVERLAPPING_FROZEN_CANDIDATE_INTERVALS"})
    groups = defaultdict(list)
    for locus in loci:
        groups[find(locus["locus_id"])].append(locus)
    group_rows = []
    member_rows = []
    for members in groups.values():
        members.sort(key=lambda r: (r["pair_id"], int(r["START"]), r["locus_id"]))
        chrom = members[0]["CHR"]
        start = min(int(r["START"]) for r in members)
        stop = max(int(r["STOP"]) for r in members)
        group_id = f"chr{chrom}:{start}-{stop}"
        pairs = sorted({r["pair_id"] for r in members})
        group_rows.append({"region_group": group_id, "chromosome": chrom,
                           "start": start, "stop": stop,
                           "n_pair_specific_candidates": len(members),
                           "n_distinct_pairs": len(pairs), "pairs": ";".join(pairs),
                           "candidate_locus_ids": ";".join(r["locus_id"] for r in members),
                           "lead_variants": ";".join(sorted({lead for r in members
                                                              for lead in r["lead_variants"].split(";")})),
                           "status": "CROSS_PAIR_GEOGRAPHIC_GROUP_ONLY"})
        for locus in members:
            member_rows.append({"region_group": group_id, "pair_id": locus["pair_id"],
                                "candidate_locus_id": locus["locus_id"],
                                "lead_variants": locus["lead_variants"],
                                "lead_p_placo": locus["lead_P_PLACO"],
                                "status": "PAIR_SPECIFIC_CANDIDATE_UNCHANGED"})
    group_rows.sort(key=lambda r: (int(r["chromosome"]), int(r["start"]), int(r["stop"])))
    member_rows.sort(key=lambda r: (r["region_group"], r["pair_id"], r["candidate_locus_id"]))
    edges.sort(key=lambda r: (int(r["chromosome"]), r["locus_a"], r["locus_b"]))
    if len(group_rows) + len(edges) < len(loci):
        raise ValueError("Unexpected overlap graph accounting")

    OUT.mkdir(parents=True)
    write(OUT / "geographic_region_groups.tsv", group_rows)
    write(OUT / "pair_candidate_members.tsv", member_rows)
    write(OUT / "cross_pair_overlap_edges.tsv", edges)
    result = {
        "analysis_id": "brain6-five-track-cross-pair-region-reconciliation-v1",
        "status": "PASS_GEOGRAPHIC_RECONCILIATION_DIAGNOSTIC_ONLY",
        "input_candidate_loci": len(loci),
        "geographic_region_groups": len(group_rows),
        "multi_pair_groups": sum(int(row["n_distinct_pairs"]) > 1 for row in group_rows),
        "cross_pair_overlap_edges": len(edges),
        "unassigned_no_exact_reference_variants": unassigned_no_reference,
        "pair_specific_counts": dict(Counter(r["pair_id"] for r in loci)),
        "source_sha256": {str(path.relative_to(ROOT)): sha(path)
                          for path in (loci_path, variants_path, SOURCE / "provenance.json", RULE)},
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {path.name: sha(path) for path in OUT.glob("*.tsv")},
        "interpretation": "Cross-pair coordinate grouping of frozen candidate intervals only. The pair-specific UKB clumps and all region promotion gates remain unchanged; grouping does not prove independent signals or valid local rg.",
    }
    (OUT / "provenance.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("status", "input_candidate_loci", "geographic_region_groups",
                       "multi_pair_groups", "cross_pair_overlap_edges")}, sort_keys=True))


if __name__ == "__main__":
    main()
