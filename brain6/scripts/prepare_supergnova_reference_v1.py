#!/usr/bin/env python3
"""Make an auditable GRCh37 EUR PLINK reference for SUPERGNOVA.

This only reads reference genotypes and a published genetic map. It never
opens GWAS association statistics or candidate outcomes.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import pathlib
import shutil
import zipfile

import numpy as np


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def genetic_maps(archive: pathlib.Path) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    result = {}
    with zipfile.ZipFile(archive) as zf:
        for chrom in range(1, 23):
            name = f"plink.chr{chrom}.GRCh37.map"
            bp, cm = [], []
            with zf.open(name) as stream:
                for line in stream:
                    fields = line.split()
                    if len(fields) != 4:
                        raise ValueError(f"Malformed genetic-map row: {name}")
                    bp.append(int(fields[3]))
                    cm.append(float(fields[2]))
            positions = np.asarray(bp, dtype=np.int64)
            distances = np.asarray(cm, dtype=np.float64)
            if not np.all(np.diff(positions) > 0) or not np.all(np.diff(distances) >= 0):
                raise ValueError(f"Nonmonotone genetic map: {name}")
            result[chrom] = positions, distances
    return result


def convert_partitions(source: pathlib.Path, target: pathlib.Path) -> int:
    count = 0
    last = {}
    with source.open() as inp, target.open("w") as out:
        out.write("chr\tstart\tend\n")
        header = next(inp).split()
        if header != ["chr", "start", "stop"]:
            raise ValueError(f"Unexpected LDetect header: {header}")
        for line in inp:
            chrom_s, start_s, end_s = line.split()
            chrom = int(chrom_s.removeprefix("chr"))
            start, end = int(start_s) + 1, int(end_s)
            if not 1 <= chrom <= 22 or start > end:
                raise ValueError(f"Invalid LDetect interval: {line.strip()}")
            if chrom in last and start <= last[chrom]:
                raise ValueError(f"Overlapping LDetect intervals: {line.strip()}")
            last[chrom] = end
            if chrom == 6 and start <= 34_000_000 and end >= 25_000_000:
                continue  # published SUPERGNOVA application excluded MHC
            out.write(f"{chrom}\t{start}\t{end}\n")
            count += 1
    return count


def prepare(source: pathlib.Path, output: pathlib.Path) -> dict:
    output.mkdir(parents=True, exist_ok=False)
    prefix = source / "1000G_EUR"
    bed, bim, fam, afreq = (prefix.with_suffix(s) for s in (".bed", ".bim", ".fam", ".afreq"))
    mapzip = source / "plink.GRCh37.map.zip"
    blocks = source / "LDetect_EUR_GRCh37.bed"
    inputs = [bed, bim, fam, afreq, mapzip, blocks]
    if any(not p.is_file() for p in inputs):
        raise FileNotFoundError([str(p) for p in inputs if not p.is_file()])
    map_by_chr = genetic_maps(mapzip)
    n_samples = sum(1 for _ in fam.open())
    stride = math.ceil(n_samples / 4)
    if n_samples != 503 or bed.stat().st_size != 3 + stride * sum(1 for _ in bim.open()):
        raise ValueError("Reference BED/BIM/FAM dimensions disagree")

    # The archive's AFREQ and BIM files share order and SNP ID. Count IDs first
    # so duplicated IDs are removed from both the reference and GWAS merge.
    ids = collections.Counter()
    with bim.open() as handle:
        for line in handle:
            ids[line.split()[1]] += 1
    kept = 0
    exclusion = collections.Counter()
    dest = output / "eur_maf05_nonpal"
    with bed.open("rb") as bed_in, bim.open() as bim_in, afreq.open() as freq_in, \
            dest.with_suffix(".bed").open("wb") as bed_out, \
            dest.with_suffix(".bim").open("w") as bim_out:
        magic = bed_in.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError("Input BED is not variant-major PLINK 1")
        bed_out.write(magic)
        if next(freq_in).split() != ["#CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"]:
            raise ValueError("Unexpected allele-frequency header")
        for i, (bline, fline) in enumerate(zip(bim_in, freq_in), start=1):
            genotype = bed_in.read(stride)
            if len(genotype) != stride:
                raise ValueError("Truncated BED")
            b = bline.split()
            f = fline.split()
            if len(b) != 6 or len(f) != 6 or b[1] != f[1] or b[0] != f[0]:
                raise ValueError(f"BIM/AFREQ mismatch on row {i}")
            chrom, snp, _, bp, a1, a2 = b
            if ids[snp] != 1 or not snp.startswith("rs"):
                exclusion["duplicate_or_non_rs"] += 1
                continue
            if len(a1) != 1 or len(a2) != 1 or a1 not in "ACGT" or a2 not in "ACGT":
                exclusion["non_snp"] += 1
                continue
            if {a1, a2} in ({"A", "T"}, {"C", "G"}):
                exclusion["strand_ambiguous"] += 1
                continue
            frequency = float(f[4])
            if not 0.05 <= frequency <= 0.95:
                exclusion["maf_below_05"] += 1
                continue
            chrom_i = int(chrom)
            positions, cm_values = map_by_chr[chrom_i]
            cm = float(np.interp(int(bp), positions, cm_values))
            if not math.isfinite(cm):
                raise ValueError(f"Missing cM: {snp}")
            bed_out.write(genotype)
            bim_out.write(f"{chrom}\t{snp}\t{cm:.9f}\t{bp}\t{a1}\t{a2}\n")
            kept += 1
        if next(bim_in, None) is not None or next(freq_in, None) is not None or bed_in.read(1):
            raise ValueError("Reference input row count disagrees")
    shutil.copyfile(fam, dest.with_suffix(".fam"))
    n_blocks = convert_partitions(blocks, output / "ldetect_eur_nonmhc.tsv")
    if dest.with_suffix(".bed").stat().st_size != 3 + stride * kept:
        raise ValueError("Filtered BED length mismatch")
    manifest = {
        "n_samples": n_samples,
        "n_input_variants": sum(ids.values()),
        "n_kept_variants": kept,
        "exclusion_counts": dict(exclusion),
        "n_ld_blocks_nonmhc": n_blocks,
        "source_sha256": {p.name: sha256(p) for p in inputs},
        "output_sha256": {p.name: sha256(p) for p in output.iterdir()
                          if p.is_file() and not p.name.startswith("._")},
        "genetic_map": "Beagle GRCh37 HapMap II map, linear interpolation of cM at reference BP",
        "partition_conversion": "LDetect EUR GRCh37 BED 0-based half-open to 1-based inclusive; exclude blocks overlapping chr6:25-34Mb",
    }
    (output / "reference_prep_provenance.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.output), indent=2, sort_keys=True))
