#!/usr/bin/env python3
"""Source-only MDD2025 coverage screen against the frozen Brain6 LAVA panel.

Reads no association statistics. The immutable config fixes source identity and
filters before coverage is observed. Writes a new diagnostic receipt only.
"""
from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "brain6/config/lava_multitrait_feasibility_v1/mdd2025_source_screen.json"
OUTPUT = ROOT / "brain6/results/lava_multitrait_feasibility_v1/mdd2025_source_coverage_v1.json"
BASES = set("ACGT")
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
SOURCE_FIELDS = {"#CHROM", "POS", "ID", "EA", "NEA", "FCAS", "FCON", "IMPINFO", "NEFF", "NCAS", "NCON"}
REFERENCE_FIELDS = ("SNP", "CHR", "POS", "A1", "A2", "NOBS", "MISS", "FREQ", "NCORRS")


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_reference(path: Path, chromosome: int) -> dict[str, tuple[int, frozenset[str]]]:
    reference = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REFERENCE_FIELDS:
            raise ValueError(f"Unexpected reference schema: {path}")
        for row in reader:
            if int(row["CHR"]) != chromosome:
                raise ValueError("Reference chromosome mismatch")
            key = row["SNP"].lower()
            if key in reference:
                raise ValueError(f"Duplicate reference SNP: {key}")
            reference[key] = (int(row["POS"]), frozenset((row["A1"].upper(), row["A2"].upper())))
    return reference


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    source = Path(config["source_path"])
    if digest(source, "md5") != config["source_md5"] or digest(source) != config["source_sha256"]:
        raise ValueError("Source archive checksum mismatch")
    sidecar = source.parent / "pgc-mdd2025_no23andMe_eur_v3.49.24.11.txt"
    if digest(sidecar, "md5") != config["cohort_sidecar_md5"]:
        raise ValueError("Cohort sidecar checksum mismatch")
    locus_path = ROOT / config["locus_file"]
    loci_by_chr: dict[int, list[tuple[int, int, str]]] = {chrom: [] for chrom in range(1, 23)}
    with locus_path.open(encoding="utf-8") as stream:
        if stream.readline().split() != ["LOC", "CHR", "START", "STOP"]:
            raise ValueError("Frozen locus schema mismatch")
        for line in stream:
            if not line.strip():
                continue
            locus_id, chrom, start, stop = line.split()
            loci_by_chr[int(chrom)].append((int(start), int(stop), locus_id))
    if sum(map(len, loci_by_chr.values())) != config["locus_count"]:
        raise ValueError("Frozen locus count mismatch")
    for loci in loci_by_chr.values():
        loci.sort()
    counters = Counter()
    chromosome_counts: dict[str, dict[str, int]] = {}
    locus_counts: dict[str, int] = {loc[2]: 0 for loci in loci_by_chr.values() for loc in loci}
    current_chr = 0
    reference: dict[str, tuple[int, frozenset[str]]] = {}
    seen: set[str] = set()
    with gzip.open(source, "rt", encoding="utf-8", newline="") as stream:
        for line in stream:
            if not line.startswith("##"):
                header = line.rstrip("\r\n").split("\t")
                break
        else:
            raise ValueError("Missing source table header")
        if not SOURCE_FIELDS.issubset(header):
            raise ValueError(f"Missing source fields: {sorted(SOURCE_FIELDS - set(header))}")
        idx = {name: header.index(name) for name in SOURCE_FIELDS}
        for line in stream:
            counters["source_rows"] += 1
            parts = line.rstrip("\r\n").split("\t")
            try:
                chrom = int(parts[idx["#CHROM"]])
                pos = int(parts[idx["POS"]])
            except (ValueError, IndexError):
                counters["invalid_coordinates"] += 1
                continue
            if not 1 <= chrom <= 22:
                counters["nonautosomal_rows"] += 1
                continue
            if chrom != current_chr:
                if chrom <= current_chr:
                    raise ValueError("Source chromosome order is not monotonic")
                current_chr = chrom
                ref_path = ROOT / f"{config['reference_prefix']}{chrom}.info"
                reference = load_reference(ref_path, chrom)
                chromosome_counts[str(chrom)] = {"reference_ids": len(reference), "source_rows": 0,
                                                   "exact_id_and_allele": 0, "strict_qc_exact": 0}
                seen = set()
            chr_counts = chromosome_counts[str(chrom)]
            chr_counts["source_rows"] += 1
            snp = parts[idx["ID"]].lower()
            if snp in seen:
                counters["duplicate_source_ids"] += 1
                continue
            seen.add(snp)
            hit = reference.get(snp)
            if hit is None:
                counters["missing_reference_id"] += 1
                continue
            a1, a2 = parts[idx["EA"]].upper(), parts[idx["NEA"]].upper()
            alleles = frozenset((a1, a2))
            if pos != hit[0] or alleles != hit[1] or len(alleles) != 2 or not alleles <= BASES:
                counters["position_or_allele_mismatch"] += 1
                continue
            chr_counts["exact_id_and_allele"] += 1
            counters["exact_id_and_allele"] += 1
            if config["exclude_palindromic"] and alleles in PALINDROMIC:
                counters["palindromic"] += 1
                continue
            try:
                fcase = float(parts[idx["FCAS"]])
                fcontrol = float(parts[idx["FCON"]])
                ncase = float(parts[idx["NCAS"]])
                ncontrol = float(parts[idx["NCON"]])
                neff = float(parts[idx["NEFF"]])
                info = float(parts[idx["IMPINFO"]])
                freq = (fcase * ncase + fcontrol * ncontrol) / (ncase + ncontrol)
            except (ValueError, ZeroDivisionError):
                counters["invalid_qc_fields"] += 1
                continue
            if not all(map(math.isfinite, (freq, neff, info))) or not 0 < freq < 1:
                counters["invalid_qc_fields"] += 1
                continue
            if info < config["minimum_info"] or min(freq, 1 - freq) < config["minimum_maf"] or neff < config["minimum_neff"]:
                counters["below_qc_threshold"] += 1
                continue
            chr_counts["strict_qc_exact"] += 1
            counters["strict_qc_exact"] += 1
            loci = loci_by_chr[chrom]
            starts = [x[0] for x in loci] if "starts" not in chr_counts else chr_counts["starts"]
            # Cache starts only within the current chromosome, without placing
            # that internal list in the final JSON receipt.
            chr_counts["starts"] = starts
            index = bisect.bisect_right(starts, pos) - 1
            if index >= 0 and pos <= loci[index][1]:
                locus_counts[loci[index][2]] += 1
    for counts in chromosome_counts.values():
        counts.pop("starts", None)
    result = {
        "analysis_id": config["analysis_id"],
        "status": "PASS_SOURCE_COVERAGE_DIAGNOSTIC_ONLY",
        "source_sha256": config["source_sha256"],
        "source_md5": config["source_md5"],
        "cohort_sidecar_sha256": digest(sidecar),
        "config_sha256": digest(CONFIG),
        "locus_file_sha256": digest(locus_path),
        "counts": dict(counters),
        "chromosomes": chromosome_counts,
        "locus_coverage": {"loci": len(locus_counts), "with_at_least_one_strict_qc_variant": sum(v > 0 for v in locus_counts.values()),
                           "below_30_strict_qc_variants": sum(v < 30 for v in locus_counts.values()),
                           "minimum": min(locus_counts.values()), "maximum": max(locus_counts.values())},
        "interpretation": config["interpretation_rule"],
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "counts": result["counts"],
                      "locus_coverage": result["locus_coverage"]}, sort_keys=True))


if __name__ == "__main__":
    main()
