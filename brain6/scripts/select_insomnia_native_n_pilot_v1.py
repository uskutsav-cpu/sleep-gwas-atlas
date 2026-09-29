#!/usr/bin/env python3
"""Select 88 insomnia diagnostic loci without association values or candidates."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCI = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
CANDIDATES = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv"
CANONICAL = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
DECISION = CANONICAL.parents[1] / "canonical_family_decision.json"
OUT = ROOT / "brain6/results/confirmatory_source_rescue_20260927/other_traits"
EXPECTED = {
    LOCI: "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882",
    CANDIDATES: "fffc0e23ab8c111385d78d00168d8dbdbf20b935e7a8fa4b7b6a77e754f3180c",
    CANONICAL: "ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352",
    DECISION: "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22",
}
PREFIX = "brain6-insomnia-native-n-pilot-v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path, sep: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter=sep))


def main() -> None:
    for path, expected in EXPECTED.items():
        if sha(path) != expected:
            raise ValueError(f"Frozen input changed: {path}")
    loci = read(LOCI, " ")
    candidates = read(CANDIDATES)
    if len(loci) != 2495 or len(candidates) != 25:
        raise ValueError("Locus or candidate count changed")
    locus_by_id = {row["LOC"]: row for row in loci}
    if len(locus_by_id) != 2495:
        raise ValueError("Duplicate locus ID")
    excluded = {row["LOC"] for row in loci for candidate in candidates
                if row["CHR"] == candidate["CHR"]
                and int(row["START"]) <= int(candidate["STOP"])
                and int(candidate["START"]) <= int(row["STOP"])}
    if len(excluded) != 38:
        raise ValueError("Protected candidate overlap changed")
    original = [r for r in read(CANONICAL) if r["phen"] == "insomnia"]
    if len(original) != 2495 or Counter(r["status"] for r in original) != {"TESTED": 1824, "NOT_RUN": 671}:
        raise ValueError("Canonical insomnia stratum changed")
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in original:
        loc = locus_by_id[row["locus_id"]]
        if loc["CHR"] != row["chromosome"]:
            raise ValueError("Canonical chromosome mismatch")
        if row["locus_id"] not in excluded:
            groups[(loc["CHR"], row["status"])].append(row["locus_id"])
    chosen = []
    for chrom in map(str, range(1, 23)):
        for status in ("TESTED", "NOT_RUN"):
            eligible = groups[(chrom, status)]
            if len(eligible) < 2:
                raise ValueError(f"Insufficient loci in chr{chrom} {status}")
            ranked = sorted(eligible, key=lambda loc: (
                hashlib.sha256(f"{PREFIX}|{chrom}|{status}|{loc}".encode()).hexdigest(), int(loc)))
            for rank, loc in enumerate(ranked[:2], 1):
                source = locus_by_id[loc]
                chosen.append({"LOC": loc, "CHR": chrom, "START": source["START"],
                               "STOP": source["STOP"], "canonical_status": status,
                               "stratum_rank": rank})
    if len(chosen) != 88 or len({r["LOC"] for r in chosen}) != 88:
        raise ValueError("Selected family invalid")
    table = OUT / "insomnia_native_n_88_loci.tsv"
    receipt = OUT / "insomnia_native_n_88_loci.selection.json"
    with table.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(chosen[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(chosen)
    metadata = {"analysis_id": PREFIX, "n_loci": 88, "n_candidates": 25,
                "n_candidate_overlapping_loci_excluded": 38,
                "selection": "two lowest SHA-256 ranks per chromosome and canonical insomnia TESTED/NOT_RUN stratum",
                "association_effect_or_p_fields_used": False,
                "input_sha256": {str(p): digest for p, digest in EXPECTED.items()},
                "table_sha256": sha(table),
                "status_counts": dict(Counter(r["canonical_status"] for r in chosen))}
    with receipt.open("x", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"selected": 88, "excluded": 38, "selection_sha256": sha(table),
                      "receipt_sha256": sha(receipt)}, sort_keys=True))


if __name__ == "__main__":
    main()
