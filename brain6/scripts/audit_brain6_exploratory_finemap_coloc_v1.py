#!/usr/bin/env python3
"""Pre-outcome dense-GWAS and LD coverage audit for the 25 protected PLACO loci.

Reads association columns only to test that they are finite; never writes effect,
standard-error, P value, or Z values. This is an input-admission audit, not a
fine-mapping or colocalization result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
import statistics
from collections import Counter, defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CANDIDATES = REPO / "brain6/results/brain6_bounded_manuscript_v1/candidate_evidence_25.tsv"
SSD = Path("/Volumes/Extreme SSD/brain6-work")
GWAS = SSD / "preparation-v1"
LD = SSD / "lava-ukb-v1.1"
OUT = REPO / "brain6/results/brain6_exploratory_finemap_coloc_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def orient(a1: str, a2: str, ref1: str, ref2: str) -> str:
    alleles = {"A", "C", "G", "T"}
    a1, a2, ref1, ref2 = (v.upper() for v in (a1, a2, ref1, ref2))
    if not all(x in alleles for x in (a1, a2, ref1, ref2)) or a1 == a2:
        return "INVALID_ALLELE"
    if {a1, a2} in ({"A", "T"}, {"C", "G"}):
        return "PALINDROMIC_EXCLUDED"
    if (a1, a2) == (ref1, ref2):
        return "EXACT"
    if (a1, a2) == (ref2, ref1):
        return "SWAP"
    comp = str.maketrans("ACGT", "TGCA")
    if (a1.translate(comp), a2.translate(comp)) == (ref1, ref2):
        return "STRAND_COMPLEMENT"
    if (a1.translate(comp), a2.translate(comp)) == (ref2, ref1):
        return "STRAND_COMPLEMENT_SWAP"
    return "ALLELE_MISMATCH"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    candidates = list(csv.DictReader(CANDIDATES.open(), delimiter="\t"))
    if len(candidates) != 25 or len({c["candidate_locus_id"] for c in candidates}) != 25:
        raise RuntimeError("Expected 25 unique protected pair-specific candidates")
    loci = []
    by_chr = defaultdict(list)
    for c in candidates:
        match = re.search(r"_chr(\d+)_(\d+)_(\d+)$", c["candidate_locus_id"])
        if not match:
            raise RuntimeError(f"Cannot parse frozen candidate bounds: {c['candidate_locus_id']}")
        chrom, start, stop = map(int, match.groups())
        if chrom not in range(1, 23) or start >= stop:
            raise RuntimeError("Invalid candidate interval")
        obj = dict(c, chromosome=chrom, start=start, stop=stop)
        loci.append(obj)
        by_chr[chrom].append(obj)
    ref = defaultdict(dict)
    ref_files = {}
    for chrom, chr_loci in by_chr.items():
        path = LD / f"lava-ukb-v1.1_chr{chrom}.info"
        ref_files[str(path)] = sha256(path)
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"SNP", "CHR", "POS", "A1", "A2"}
            if not required.issubset(reader.fieldnames or []):
                raise RuntimeError(f"Invalid reference INFO header: {path}")
            for row in reader:
                pos = int(row["POS"])
                if any(loc["start"] <= pos <= loc["stop"] for loc in chr_loci):
                    snp = row["SNP"].lower()
                    if snp in ref[chrom]:
                        raise RuntimeError(f"Duplicate LD reference SNP {snp}")
                    ref[chrom][snp] = (pos, row["A1"], row["A2"])

    connections = {}
    trait_paths = {}
    for trait in sorted({t for c in loci for t in c["pair_id"].split("__")}):
        db = GWAS / f"normalize_{trait}/variants.sqlite"
        if not db.exists():
            raise RuntimeError(f"Missing normalized dense GWAS: {db}")
        connections[trait] = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        trait_paths[trait] = db

    rows = []
    oriented_sets = {}
    for locus in loci:
        chrom, start, stop = (locus[k] for k in ("chromosome", "start", "stop"))
        window_ref = {s: x for s, x in ref[chrom].items() if start <= x[0] <= stop}
        for trait in locus["pair_id"].split("__"):
            conn = connections[trait]
            source = conn.execute(
                "SELECT snp,chr,bp,a1,a2,beta,se,p,n FROM variants "
                "WHERE chr=? AND bp>=? AND bp<=?", (chrom, start, stop)
            )
            counts = Counter()
            admitted = set()
            admitted_n = []
            for snp, chr_val, pos, a1, a2, beta, se, p, n in source:
                counts["gwas_rows"] += 1
                if not all(math.isfinite(x) for x in (beta, se, p, n)) or se <= 0 or n <= 0 or not (0 < p <= 1):
                    counts["invalid_statistics"] += 1
                    continue
                counts["valid_statistics"] += 1
                entry = window_ref.get(snp.lower())
                if entry is None:
                    counts["reference_absent"] += 1
                    continue
                if entry[0] != pos:
                    counts["position_mismatch"] += 1
                    continue
                status = orient(a1, a2, entry[1], entry[2])
                counts[status] += 1
                if status in {"EXACT", "SWAP", "STRAND_COMPLEMENT", "STRAND_COMPLEMENT_SWAP"}:
                    admitted.add(snp.lower())
                    admitted_n.append(n)
            key = (locus["candidate_locus_id"], trait)
            oriented_sets[key] = admitted
            n_median = statistics.median(admitted_n) if admitted_n else None
            rows.append({
                "candidate_locus_id": locus["candidate_locus_id"],
                "pair_id": locus["pair_id"],
                "trait_id": trait,
                "chromosome": chrom,
                "start": start,
                "stop": stop,
                "reference_variants": len(window_ref),
                "gwas_rows": counts["gwas_rows"],
                "valid_statistics": counts["valid_statistics"],
                "invalid_statistics": counts["invalid_statistics"],
                "reference_absent": counts["reference_absent"],
                "position_mismatch": counts["position_mismatch"],
                "palindromic_excluded": counts["PALINDROMIC_EXCLUDED"],
                "allele_mismatch": counts["ALLELE_MISMATCH"] + counts["INVALID_ALLELE"],
                "exact": counts["EXACT"],
                "swap": counts["SWAP"],
                "strand_complement": counts["STRAND_COMPLEMENT"],
                "strand_complement_swap": counts["STRAND_COMPLEMENT_SWAP"],
                "oriented_ld_variants": len(admitted),
                "admitted_n_min": min(admitted_n) if admitted_n else "",
                "admitted_n_median": n_median if n_median is not None else "",
                "admitted_n_max": max(admitted_n) if admitted_n else "",
                "admitted_n_within_20pct_fraction": round(
                    sum(0.8 * n_median <= x <= 1.2 * n_median for x in admitted_n) / len(admitted_n), 6
                ) if admitted_n else 0,
                "reference_coverage_fraction": round(len(admitted) / len(window_ref), 6) if window_ref else 0,
                "gwas_coverage_fraction": round(len(admitted) / counts["gwas_rows"], 6) if counts["gwas_rows"] else 0,
            })
    pair_rows = []
    for locus in loci:
        t1, t2 = locus["pair_id"].split("__")
        shared = oriented_sets[(locus["candidate_locus_id"], t1)] & oriented_sets[(locus["candidate_locus_id"], t2)]
        pair_rows.append({
            "candidate_locus_id": locus["candidate_locus_id"],
            "pair_id": locus["pair_id"],
            "chromosome": locus["chromosome"],
            "start": locus["start"],
            "stop": locus["stop"],
            "trait1_ld_variants": len(oriented_sets[(locus["candidate_locus_id"], t1)]),
            "trait2_ld_variants": len(oriented_sets[(locus["candidate_locus_id"], t2)]),
            "shared_oriented_ld_variants": len(shared),
            "lead_in_shared": "YES" if any(v.lower() in shared for v in locus["lead_variants"].split(";")) else "NO",
        })
    write_tsv(args.out / "trait_locus_input_audit.tsv", rows, list(rows[0]))
    write_tsv(args.out / "pair_locus_input_audit.tsv", pair_rows, list(pair_rows[0]))
    provenance = {
        "status": "PRE_OUTCOME_INPUT_AUDIT",
        "candidate_sha256": sha256(CANDIDATES),
        "script_sha256": sha256(Path(__file__)),
        "reference_version_sha256": sha256(LD / "VERSION_INFO"),
        "reference_info_sha256": ref_files,
        "normalized_sqlite_sha256": {trait: sha256(path) for trait, path in trait_paths.items()},
        "source_receipt_sha256": {trait: sha256(path.parent / "receipt.json") for trait, path in trait_paths.items()},
        "output_sha256": {
            name: sha256(args.out / name)
            for name in ("trait_locus_input_audit.tsv", "pair_locus_input_audit.tsv")
        },
        "limitations": [
            "No LD matrix read or numerical conditioning assessed by this script",
            "No association outcome, credible set, or colocalization posterior inspected",
            "Palindromic alleles excluded irrespective of allele frequency",
        ],
    }
    (args.out / "input_audit_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    for conn in connections.values():
        conn.close()
    print(f"audited {len(loci)} candidates; {len(rows)} trait-locus rows")


if __name__ == "__main__":
    main()
