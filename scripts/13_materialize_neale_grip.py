#!/usr/bin/env python3
"""Materialize the public Neale Lab round-2 left-grip-strength release.

The association file identifies variants as ``chr:pos:ref:alt``.  The companion
variant annotation supplied in the same public release is therefore required to
obtain documented rsIDs and the alternate-allele frequency/INFO fields.  This
program streams the two bgzip files in lockstep, refuses a row-order mismatch,
and never performs a liftover or inferred coordinate-to-rsID lookup.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import os
from pathlib import Path
import re
import sys
import tempfile


RESULTS_BYTES = 539_248_977
RESULTS_MD5 = "e93db6d09ffa591e5a27a77492dbd21d"
RESULTS_SHA256 = "da2ceb5e59c2d0edf2c405d0e6e918c7a07876d3dfb3f42b1b3c43133ca65b55"
VARIANTS_BYTES = 916_143_947
VARIANTS_MD5 = "4fc9936ba4b4fd446cb4325dd63b6e72"
VARIANTS_SHA256 = "e7f035e9264536b416e7b7a514dd863bd508dc380f40e48f12ffbc7c06891ea5"

RESULTS_COLUMNS = [
    "variant", "minor_allele", "minor_AF", "low_confidence_variant",
    "n_complete_samples", "AC", "ytx", "beta", "se", "tstat", "pval",
]
VARIANT_COLUMNS = [
    "variant", "chr", "pos", "ref", "alt", "rsid", "varid",
    "consequence", "consequence_category", "info", "call_rate", "AC",
    "AF", "minor_allele", "minor_AF", "p_hwe", "n_called",
    "n_not_called", "n_hom_ref", "n_het", "n_hom_var", "n_non_ref",
    "r_heterozygosity", "r_het_hom_var", "r_expected_het_frequency",
]
RSID = re.compile(r"rs[0-9]+$", re.IGNORECASE)
ALLELES = {"A", "C", "G", "T"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def digest(path: Path) -> tuple[int, str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            md5.update(chunk)
            sha256.update(chunk)
    return size, md5.hexdigest(), sha256.hexdigest()


def verify_file(path: Path, expected_bytes: int, expected_md5: str, expected_sha256: str) -> None:
    if not path.is_file():
        fail(f"missing required source file: {path}")
    actual_bytes, actual_md5, actual_sha256 = digest(path)
    if actual_bytes != expected_bytes:
        fail(f"unexpected byte count for {path.name}: expected {expected_bytes}, got {actual_bytes}")
    if actual_md5 != expected_md5:
        fail(f"MD5 mismatch for {path.name}")
    if actual_sha256 != expected_sha256:
        fail(f"SHA-256 mismatch for {path.name}")


def checked_reader(path: Path, expected_columns: list[str], label: str) -> csv.DictReader:
    handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames != expected_columns:
        handle.close()
        fail(
            f"unexpected {label} header in {path.name}: expected {expected_columns}, "
            f"got {reader.fieldnames}"
        )
    # Keep a reference so the file remains open for the caller's full stream.
    reader._source_handle = handle  # type: ignore[attr-defined]
    return reader


def materialize(results: Path, variants: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    results_reader = checked_reader(results, RESULTS_COLUMNS, "results")
    variants_reader = checked_reader(variants, VARIANT_COLUMNS, "variant annotation")
    result_handle = results_reader._source_handle  # type: ignore[attr-defined]
    variant_handle = variants_reader._source_handle  # type: ignore[attr-defined]

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".partial", dir=output.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    total = kept = non_rsid = non_snp = low_confidence = 0
    try:
        with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N", "INFO"])
            for line_number, pair in enumerate(
                itertools.zip_longest(results_reader, variants_reader), start=2
            ):
                result, variant = pair
                if result is None or variant is None:
                    fail(f"source files have different row counts near line {line_number}")
                total += 1
                if result["variant"] != variant["variant"]:
                    fail(
                        "source files are not in the documented shared variant order at "
                        f"line {line_number}: {result['variant']} != {variant['variant']}"
                    )
                rsid = variant["rsid"].strip().lower()
                if not RSID.fullmatch(rsid):
                    non_rsid += 1
                    continue
                ref, alt = variant["ref"].strip().upper(), variant["alt"].strip().upper()
                if ref not in ALLELES or alt not in ALLELES or ref == alt:
                    non_snp += 1
                    continue
                if result["low_confidence_variant"].strip().lower() == "true":
                    low_confidence += 1
                    continue
                writer.writerow([
                    rsid, variant["chr"], variant["pos"], alt, ref, variant["AF"],
                    result["beta"], result["se"], result["pval"],
                    result["n_complete_samples"], variant["info"],
                ])
                kept += 1
        # Read the complete compressed stream so gzip validates its trailer before
        # a result becomes visible to the Phase 0 harmonizer.
        with gzip.open(temporary, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        result_handle.close()
        variant_handle.close()

    with output.open("rb") as handle:
        output_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    print(f"Materialized {output}")
    print(f"  source rows: {total:,}")
    print(f"  retained rsID SNPs: {kept:,}")
    print(f"  skipped non-rsID: {non_rsid:,}")
    print(f"  skipped non-SNP: {non_snp:,}")
    print(f"  skipped source low-confidence: {low_confidence:,}")
    print(f"  output sha256: {output_sha256}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if args.verify_only == (args.out is not None):
        fail("supply exactly one of --verify-only or --out")

    verify_file(args.results, RESULTS_BYTES, RESULTS_MD5, RESULTS_SHA256)
    verify_file(args.variants, VARIANTS_BYTES, VARIANTS_MD5, VARIANTS_SHA256)
    print("Verified Neale Lab left-grip-strength source files")
    if args.verify_only:
        return
    materialize(args.results, args.variants, args.out)


if __name__ == "__main__":
    main()
