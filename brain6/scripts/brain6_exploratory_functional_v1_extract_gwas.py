#!/usr/bin/env python3
"""Extract raw, variant-N-bearing insomnia/ADHD GWAS at six frozen windows."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path

from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[2]
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
CHAIN = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz")
OUT = Path("/Volumes/Extreme SSD/brain6-work/brain6-exploratory-functional-v1/gwas_regional")
SOURCES = {
    "insomnia": (Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/Insomnia_sumstats_Jansenetal.txt.gz"),
                 "32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b"),
    "adhd": (Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/ADHD2022_iPSYCH_deCODE_PGC.meta.gz"),
             "c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5"),
}
COLUMNS = ("trait", "region_grch37", "snp", "chr37", "bp37", "a1", "a2", "chr38", "bp38",
           "beta_log_or", "se_log_or", "p", "n_total", "n_case", "n_control", "n_semantics", "mapping_status")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with MEMBERS.open(newline="") as f:
        selected = [r for r in csv.DictReader(f, delimiter="\t") if r["pair_id"] == "insomnia__adhd"]
    if len(selected) != 6:
        raise ValueError("Frozen insomnia–ADHD candidate count changed")
    bounds = {}
    by_chr = defaultdict(list)
    for item in selected:
        region = item["region_group"]
        chrom, span = region.removeprefix("chr").split(":")
        start, stop = map(int, span.split("-"))
        bounds[region] = (chrom, start, stop)
        by_chr[chrom].append(region)
    lift = LiftOver(str(CHAIN))
    for trait, (source, expected_hash) in SOURCES.items():
        if digest(source) != expected_hash:
            raise ValueError(f"Raw GWAS hash mismatch: {source}")
        destination = {region: OUT / f"{trait}_{region.replace(':', '_').replace('-', '_')}.tsv.gz" for region in bounds}
        receipt_path = {region: OUT / f"{trait}_{region.replace(':', '_').replace('-', '_')}.receipt.json" for region in bounds}
        if any(p.exists() for p in destination.values()) or any(p.exists() for p in receipt_path.values()):
            if all(p.exists() for p in destination.values()) and all(p.exists() for p in receipt_path.values()):
                for region in bounds:
                    receipt = json.loads(receipt_path[region].read_text())
                    if receipt["sha256"] != digest(destination[region]):
                        raise ValueError(f"Existing GWAS slice changed: {destination[region]}")
                continue
            raise FileExistsError("Incomplete GWAS slice family; preserve it for audit")
        streams = {}
        writers = {}
        counts = defaultdict(int)
        best = defaultdict(lambda: 1.0)
        try:
            for region, path in destination.items():
                stream = gzip.open(str(path) + ".partial", "wt", newline="")
                streams[region] = stream
                writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writers[region] = writer
            with gzip.open(source, "rt") as input_file:
                header = input_file.readline().strip().split()
                for line in input_file:
                    values = line.strip().split()
                    if len(values) != len(header):
                        raise ValueError(f"Malformed raw GWAS line for {trait}")
                    r = dict(zip(header, values))
                    chrom = r["CHR"]
                    if chrom not in by_chr:
                        continue
                    pos = int(r["BP"])
                    targets = [region for region in by_chr[chrom] if bounds[region][1] <= pos <= bounds[region][2]]
                    if not targets:
                        continue
                    try:
                        odds, se, p = float(r["OR"]), float(r["SE"]), float(r["P"])
                        ncase = float(r["Nca"]) if trait == "adhd" else float("nan")
                        ncontrol = float(r["Nco"]) if trait == "adhd" else float("nan")
                        n = ncase + ncontrol if trait == "adhd" else float(r["N"])
                    except (KeyError, ValueError):
                        continue
                    if not (odds > 0 and se > 0 and 0 <= p <= 1 and n > 0 and
                            all(math.isfinite(v) for v in (odds, se, p, n))):
                        continue
                    a1, a2 = r["A1"].upper(), r["A2"].upper()
                    hits = lift.convert_coordinate("chr" + chrom, pos - 1)
                    valid = [h for h in hits if h[0] == "chr" + chrom and h[2] == "+"]
                    mapping = "UNIQUE_PLUS_SNP" if len(valid) == 1 and len(a1) == len(a2) == 1 and {a1, a2} <= set("ACGT") else "UNRESOLVED"
                    b38 = valid[0][1] + 1 if mapping == "UNIQUE_PLUS_SNP" else "NA"
                    row = {"trait": trait, "snp": r["SNP"], "chr37": chrom, "bp37": pos,
                           "a1": a1, "a2": a2, "chr38": chrom if mapping == "UNIQUE_PLUS_SNP" else "NA",
                           "bp38": b38, "beta_log_or": math.log(odds), "se_log_or": se,
                           "p": p, "n_total": n, "n_case": ncase if trait == "adhd" else "NA",
                           "n_control": ncontrol if trait == "adhd" else "NA",
                           "n_semantics": "RAW_NCASE_PLUS_NCONTROL" if trait == "adhd" else "RAW_LITERAL_PER_SNP_N",
                           "mapping_status": mapping}
                    for region in targets:
                        row["region_grch37"] = region
                        writers[region].writerow(row)
                        counts[region] += 1
                        best[region] = min(best[region], p)
        finally:
            for stream in streams.values():
                stream.close()
        for region, path in destination.items():
            os.replace(str(path) + ".partial", path)
            receipt = {"trait": trait, "region_grch37": region, "source_path": str(source),
                       "source_sha256": expected_hash, "chain_sha256": digest(CHAIN),
                       "member_sha256": digest(MEMBERS), "n_rows": counts[region],
                       "min_p": best[region], "sha256": digest(path),
                       "analysis_label": "EXPLORATORY_SOURCE_INPUT_ONLY"}
            with receipt_path[region].open("x") as f:
                json.dump(receipt, f, indent=2, sort_keys=True)
                f.write("\n")
            print(trait, region, counts[region], best[region], flush=True)


if __name__ == "__main__":
    main()
