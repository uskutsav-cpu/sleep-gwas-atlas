#!/usr/bin/env python3
"""Filter frozen LDSC summary files to a fixed EUR reference, without N imputation."""

from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import io
import json
import math
import pathlib


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def alleles_match(ref1: str, ref2: str, a1: str, a2: str) -> bool:
    pair = (a1, a2)
    return pair in {
        (ref1, ref2), (ref2, ref1),
        (ref1.translate(COMPLEMENT), ref2.translate(COMPLEMENT)),
        (ref2.translate(COMPLEMENT), ref1.translate(COMPLEMENT)),
    }


def prepare(reference_bim: pathlib.Path, source: pathlib.Path,
            output: pathlib.Path, trait: str, n_total: int) -> dict:
    if output.exists():
        raise FileExistsError(output)
    ref = {}
    with reference_bim.open() as handle:
        for line in handle:
            chrom, snp, cm, bp, a1, a2 = line.split()
            if snp in ref:
                raise ValueError(f"Reference duplicate: {snp}")
            ref[snp] = (a1, a2)
    counts = collections.Counter()
    kept = {}
    seen = set()
    duplicates = set()
    with gzip.open(source, "rt") as handle:
        header = next(handle).split()
        if header != ["SNP", "A1", "A2", "Z", "N"]:
            raise ValueError(f"Unexpected frozen munged header: {header}")
        for line in handle:
            counts["input_rows"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                counts["malformed"] += 1
                continue
            snp, a1, a2, z_s, _historical_n = fields
            if snp not in ref:
                counts["not_in_reference"] += 1
                continue
            if snp in seen:
                duplicates.add(snp)
                kept.pop(snp, None)
                counts["duplicate_rows"] += 1
                continue
            seen.add(snp)
            a1, a2 = a1.upper(), a2.upper()
            if len(a1) != 1 or len(a2) != 1 or a1 not in "ACGT" or a2 not in "ACGT":
                counts["invalid_alleles"] += 1
                continue
            if not alleles_match(*ref[snp], a1, a2):
                counts["allele_mismatch"] += 1
                continue
            try:
                z = float(z_s)
            except ValueError:
                counts["missing_z"] += 1
                continue
            if not math.isfinite(z):
                counts["nonfinite_z"] += 1
                continue
            kept[snp] = (a1, a2, z_s)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6,
                                                   mtime=0, filename="") as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as out:
        out.write("SNP\tA1\tA2\tZ\n")
        with reference_bim.open() as handle:
            for line in handle:
                snp = line.split()[1]
                if snp in kept:
                    a1, a2, z = kept[snp]
                    out.write(f"{snp}\t{a1}\t{a2}\t{z}\n")
    counts["output_rows"] = len(kept)
    receipt = {
        "trait": trait,
        "study_total_n_used_by_supergnova": n_total,
        "historical_munged_n_column": "ignored; no N column written to adapter output",
        "reference_bim": str(reference_bim),
        "reference_bim_sha256": sha256(reference_bim),
        "input": str(source),
        "input_sha256": sha256(source),
        "output": str(output),
        "output_sha256": sha256(output),
        "counts": dict(counts),
    }
    receipt_path = output.parent / f"{trait}.adapter_receipt.json"
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-bim", type=pathlib.Path, required=True)
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--trait", required=True)
    parser.add_argument("--n-total", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.reference_bim, args.source, args.output,
                             args.trait, args.n_total), indent=2, sort_keys=True))
