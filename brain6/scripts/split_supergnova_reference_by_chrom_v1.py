#!/usr/bin/env python3
"""Split one variant-major PLINK panel into chromosomes without changing SNPs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import shutil


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def split(prefix: pathlib.Path, output: pathlib.Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    bim = prefix.with_suffix(".bim")
    fam = prefix.with_suffix(".fam")
    bed = prefix.with_suffix(".bed")
    samples = sum(1 for _ in fam.open())
    stride = math.ceil(samples / 4)
    counts = {str(c): 0 for c in range(1, 23)}
    paths = {}
    with bim.open() as bim_in, bed.open("rb") as bed_in:
        if bed_in.read(3) != b"\x6c\x1b\x01":
            raise ValueError("Not variant-major PLINK BED")
        current = 0
        bim_out = bed_out = None
        for line in bim_in:
            chrom = int(line.split()[0])
            if chrom < current or not 1 <= chrom <= 22:
                raise ValueError("BIM not sorted by autosomal chromosome")
            if chrom != current:
                if bim_out:
                    bim_out.close()
                    bed_out.close()
                current = chrom
                dst = output / f"eur_chr{chrom}_maf05_nonpal"
                bim_out = dst.with_suffix(".bim").open("w")
                bed_out = dst.with_suffix(".bed").open("wb")
                bed_out.write(b"\x6c\x1b\x01")
                shutil.copyfile(fam, dst.with_suffix(".fam"))
                paths[str(chrom)] = str(dst)
            chunk = bed_in.read(stride)
            if len(chunk) != stride:
                raise ValueError("Truncated source BED")
            bim_out.write(line)
            bed_out.write(chunk)
            counts[str(chrom)] += 1
        if bim_out:
            bim_out.close()
            bed_out.close()
        if bed_in.read(1):
            raise ValueError("Source BED has extra variant data")
    if any(v <= 0 for v in counts.values()):
        raise ValueError("At least one chromosome absent")
    output_hashes = {}
    for chrom, dst_s in paths.items():
        dst = pathlib.Path(dst_s)
        for ext in (".bed", ".bim", ".fam"):
            p = dst.with_suffix(ext)
            output_hashes[p.name] = sha256(p)
            if ext == ".bed" and p.stat().st_size != 3 + stride * counts[chrom]:
                raise ValueError(f"Wrong BED length for chr{chrom}")
    receipt = {
        "source_prefix": str(prefix),
        "source_sha256": {ext: sha256(prefix.with_suffix(ext)) for ext in (".bed", ".bim", ".fam")},
        "n_samples": samples,
        "n_input_variants": sum(counts.values()),
        "variants_per_chr": counts,
        "output_prefix_pattern": str(output / "eur_chr@_maf05_nonpal"),
        "output_sha256": output_hashes,
        "interpretation": "Exact byte-preserving chromosome split of variant-major genotype records and BIM rows; original per-variant genotypes, alleles, positions and cM unchanged",
    }
    (output / "split_reference_provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    print(json.dumps(split(args.prefix, args.output), indent=2, sort_keys=True))
