#!/usr/bin/env python3
"""Freeze 88 long-sleep pilot loci from canonical status strata without associations."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCUS = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
CANDIDATES = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv"
CANONICAL = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
DECISION = CANONICAL.parents[1] / "canonical_family_decision.json"
OUT = ROOT / "brain6/results/lava_confirmatory_pilots_v1"
EXPECTED = {
    LOCUS: "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882",
    CANDIDATES: "fffc0e23ab8c111385d78d00168d8dbdbf20b935e7a8fa4b7b6a77e754f3180c",
    CANONICAL: "ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352",
    DECISION: "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22",
}
PREFIX = "brain6-longsleep-linear-pilot-v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main() -> None:
    for path, expected in EXPECTED.items():
        if sha(path) != expected:
            raise ValueError(f"Frozen source changed: {path}")
    loci = read_tsv(LOCUS, " ")
    candidates = read_tsv(CANDIDATES)
    if len(loci) != 2495 or len(candidates) != 25:
        raise ValueError("Frozen locus or candidate family changed")
    by_id = {row["LOC"]: row for row in loci}
    if len(by_id) != len(loci):
        raise ValueError("Duplicate locus ID")
    excluded = {row["LOC"] for row in loci for candidate in candidates
                if row["CHR"] == candidate["CHR"]
                and int(row["START"]) <= int(candidate["STOP"])
                and int(candidate["START"]) <= int(row["STOP"])}
    if len(excluded) != 38:
        raise ValueError("Candidate-overlap exclusion changed")
    canonical = [row for row in read_tsv(CANONICAL) if row["phen"] == "longsleep"]
    if len(canonical) != 2495 or Counter(row["status"] for row in canonical) != {"TESTED": 1204, "NOT_RUN": 1291}:
        raise ValueError("Canonical long-sleep stratum changed")
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in canonical:
        locus_id = row["locus_id"]
        if locus_id not in by_id or row["chromosome"] != by_id[locus_id]["CHR"]:
            raise ValueError("Canonical result and locus definition disagree")
        if locus_id not in excluded:
            groups[(row["chromosome"], row["status"])].append(locus_id)
    selected = []
    for chromosome in map(str, range(1, 23)):
        for status in ("TESTED", "NOT_RUN"):
            values = groups[(chromosome, status)]
            if len(values) < 2:
                raise ValueError("Insufficient loci in stratum")
            ranked = sorted(values, key=lambda locus_id:
                (hashlib.sha256(f"{PREFIX}|{chromosome}|{status}|{locus_id}".encode()).hexdigest(), int(locus_id)))
            for rank, locus_id in enumerate(ranked[:2], 1):
                locus = by_id[locus_id]
                selected.append({"LOC": locus_id, "CHR": chromosome,
                                 "START": locus["START"], "STOP": locus["STOP"],
                                 "canonical_status": status, "stratum_rank": rank})
    if len(selected) != 88 or len({row["LOC"] for row in selected}) != 88:
        raise ValueError("Selected locus family incomplete or duplicated")
    OUT.mkdir(parents=True, exist_ok=True)
    table = OUT / "longsleep_linear_88_loci.tsv"
    receipt = OUT / "longsleep_linear_88_loci.selection.json"
    if table.exists() or receipt.exists():
        raise FileExistsError("Prospective selection already exists")
    with table.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=("LOC", "CHR", "START", "STOP", "canonical_status", "stratum_rank"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    metadata = {"analysis_id": PREFIX, "n_loci": len(selected), "n_candidate_intervals": 25,
                "n_candidate_overlapping_loci_excluded": len(excluded),
                "selection": "two smallest SHA-256 ranks per chromosome and canonical long-sleep TESTED/NOT_RUN status",
                "association_effect_or_p_fields_used": False,
                "input_sha256": {str(path): value for path, value in EXPECTED.items()},
                "table_sha256": sha(table),
                "status_counts": dict(Counter(row["canonical_status"] for row in selected))}
    receipt.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"selected": len(selected), "excluded": len(excluded),
                      "table_sha256": sha(table), "receipt_sha256": sha(receipt)}, sort_keys=True))


if __name__ == "__main__":
    main()
