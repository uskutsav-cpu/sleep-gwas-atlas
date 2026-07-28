#!/usr/bin/env python3
"""Materialize selected Jones 2019 accelerometer sleep GWAS traits.

The GWAS Catalog release for GCST007803 is a gzip-wrapped ZIP archive whose
single member is a tab-delimited table containing several phenotypes. This
utility reads that member once and writes three standardized, gzipped inputs
for this atlas without loading the release into memory. It deliberately
selects the README-recommended rank-normalized (``RAW_SIN``) effects and
retains the source alleles/coordinates unchanged.

It is intentionally source-specific: adding another multi-phenotype release
requires its own explicit, reviewed mapping rather than guessing columns from
phenotype labels.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import sys
import zipfile
from pathlib import Path
from typing import TextIO


OUTPUTS = {
    "sleep_efficiency.txt.gz": ("ACC_SLEEP_EFF_RAW_SIN", 84_810),
    "accel_sleep_duration.txt.gz": ("ACC_SLEEP_DUR_RAW_SIN", 85_449),
    "sleep_timing.txt.gz": ("ACC_SLEEP_MIDP_RAW_SIN", 84_810),
}
BASE_COLUMNS = ("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "INFO")
OUTPUT_HEADER = ("SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "INFO")
EXPECTED_MEMBER = "accel_GWAS_all_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt"


def require_columns(header: list[str]) -> dict[str, int]:
    """Map explicit source columns, refusing a release-layout surprise."""
    positions = {name: index for index, name in enumerate(header)}
    required = set(BASE_COLUMNS)
    for prefix, _ in OUTPUTS.values():
        required.update(f"{prefix}_{field}" for field in ("A1FREQ", "BETA", "SE", "P"))
    missing = sorted(required.difference(positions))
    if missing:
        raise ValueError(
            "Unexpected GCST007803 header; required columns absent: " + ", ".join(missing)
        )
    return positions


def source_row(row: list[str], positions: dict[str, int], prefix: str) -> list[str]:
    return [
        row[positions["SNP"]],
        row[positions["CHR"]],
        row[positions["BP"]],
        row[positions["ALLELE1"]],
        row[positions["ALLELE0"]],
        row[positions[f"{prefix}_A1FREQ"]],
        row[positions[f"{prefix}_BETA"]],
        row[positions[f"{prefix}_SE"]],
        row[positions[f"{prefix}_P"]],
        row[positions["INFO"]],
    ]


def source_member(release: zipfile.ZipFile) -> str:
    members = [info.filename for info in release.infolist() if not info.is_dir()]
    if members != [EXPECTED_MEMBER]:
        raise ValueError(
            "Unexpected GCST007803 ZIP member list: " + ", ".join(members or ["<empty>"])
        )
    return EXPECTED_MEMBER


def verify_archive(source: Path) -> None:
    """Read every nested compressed byte, checking both container CRCs."""
    with gzip.open(source, "rb") as gzip_stream:
        with zipfile.ZipFile(gzip_stream) as release:
            source_member(release)
            failed_member = release.testzip()
            if failed_member:
                raise ValueError(f"ZIP CRC check failed for member {failed_member}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        required=True,
        type=Path,
        help="GCST007803 direct gzip stream downloaded from the registered URL.",
    )
    parser.add_argument("--out-dir", default="data/raw", type=Path)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Validate the gzip wrapper and ZIP member/CRC without writing raw inputs.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    source = args.source
    out_dir = args.out_dir
    if not source.is_file() or source.stat().st_size == 0:
        raise SystemExit(f"ERROR: source is missing or empty: {source}")
    if args.verify_only:
        try:
            verify_archive(source)
        except (OSError, ValueError, zipfile.BadZipFile) as exc:
            raise SystemExit(f"ERROR: archive verification failed: {exc}") from exc
        print(f"Verified gzip-wrapped ZIP source: {source}")
        return 0
    out_dir.mkdir(parents=True, exist_ok=True)

    final_paths = {name: out_dir / name for name in OUTPUTS}
    present = [str(path) for path in final_paths.values() if path.exists()]
    if present:
        raise SystemExit(
            "ERROR: refusing to overwrite existing raw input(s): " + ", ".join(present)
        )
    temporary_paths = {
        name: path.with_name(f".{path.name}.partial.{os.getpid()}")
        for name, path in final_paths.items()
    }

    writers: dict[str, csv.writer] = {}
    handles: list[TextIO] = []
    row_count = 0
    try:
        for name, temporary in temporary_paths.items():
            handle = gzip.open(temporary, "wt", encoding="utf-8", newline="")
            handles.append(handle)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(OUTPUT_HEADER)
            writers[name] = writer

        # Fully reading the sole ZIP member checks its CRC; the enclosing
        # GzipFile reaches EOF after ZIP central-directory handling, checking
        # the outer CRC too.
        with gzip.open(source, "rb") as gzip_stream:
            with zipfile.ZipFile(gzip_stream) as release:
                member = source_member(release)
                with release.open(member, "r") as member_stream:
                    with io.TextIOWrapper(member_stream, encoding="utf-8", newline="") as stream:
                        reader = csv.reader(stream, delimiter="\t")
                        try:
                            header = next(reader)
                        except StopIteration as exc:
                            raise ValueError("Source ZIP member is empty") from exc
                        positions = require_columns(header)
                        maximum_column = max(positions.values())
                        for row_number, row in enumerate(reader, start=2):
                            if len(row) <= maximum_column:
                                raise ValueError(
                                    f"Malformed source row {row_number}: expected at least "
                                    f"{maximum_column + 1} columns, found {len(row)}"
                                )
                            for name, (prefix, _) in OUTPUTS.items():
                                writers[name].writerow(source_row(row, positions, prefix))
                            row_count += 1
    except (OSError, UnicodeError, ValueError, csv.Error, zipfile.BadZipFile) as exc:
        for handle in handles:
            handle.close()
        for temporary in temporary_paths.values():
            temporary.unlink(missing_ok=True)
        raise SystemExit(f"ERROR: materialization failed: {exc}") from exc
    else:
        for handle in handles:
            handle.close()

    try:
        # Re-read each output fully: this validates output gzip CRCs before
        # atomically publishing any of the selected trait inputs.
        for temporary in temporary_paths.values():
            with gzip.open(temporary, "rb") as stream:
                while stream.read(1024 * 1024):
                    pass
        for name, final in final_paths.items():
            os.replace(temporary_paths[name], final)
    except OSError as exc:
        for temporary in temporary_paths.values():
            temporary.unlink(missing_ok=True)
        raise SystemExit(f"ERROR: output validation/publish failed: {exc}") from exc

    print(f"Materialized {row_count:,} source variants from {source}")
    for name, (prefix, n_total) in OUTPUTS.items():
        print(f"  {name}\t{prefix}\tN={n_total:,}\t{final_paths[name]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
