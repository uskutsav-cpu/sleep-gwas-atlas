#!/usr/bin/env python3
"""Materialize the PRACTICAL 2018 prostate-cancer release fail-closed.

The checksum-pinned public TSV contains a small, reproducible set of damaged
physical records.  They have the wrong number of tab-separated fields and
cannot be reconstructed without inventing association data.  This helper
retains only exact 16-field records and refuses to produce an output unless the
source-wide row counts match the reviewed release exactly.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

from panel_guard import require_locked_traits


ARCHIVE_BYTES = 2_534_497_125
ARCHIVE_SHA256 = "649921348dcc20f22a86e17bf897a2f8ae93139613a393f492d97e006a923390"
EXPECTED_FIELD_COUNT = 16
EXPECTED_DATA_ROWS = 20_065_528
EXPECTED_MALFORMED_ROWS = 367
EXPECTED_NON_RSID_ROWS = 5_887_673
EXPECTED_NON_SNP_ROWS = 1_104_805
EXPECTED_NONAUTOSOMAL_ROWS = 391_229
EXPECTED_RETAINED_ROWS = 12_681_454
HEADER = (
    b"MarkerName\trs_id\tSNP\tChr\tposition\tAllele1\tAllele2\tFreq1\tFreqSE\t"
    b"MinFreq\tMaxFreq\tEffect\tStdErr\tPvalue\tDirection\tOncoArray_imputation_r2"
)
OUTPUT_COLUMN_INDEXES = (2, 3, 4, 5, 6, 7, 11, 12, 13, 15)
OUTPUT_HEADER = (
    b"SNP\tChr\tposition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tPvalue\t"
    b"OncoArray_imputation_r2"
)
RSID = re.compile(rb"rs[0-9]+$", re.IGNORECASE)
VALID_ALLELES = {b"A", b"C", b"G", b"T"}
AUTOSOMES = {str(chromosome).encode("ascii") for chromosome in range(1, 23)}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_archive(source: Path) -> None:
    if not source.is_file():
        fail(f"source file not found: {source}")
    if source.stat().st_size != ARCHIVE_BYTES:
        fail(f"source byte count is {source.stat().st_size}, expected {ARCHIVE_BYTES}")
    if sha256(source) != ARCHIVE_SHA256:
        fail("source SHA-256 does not match the registered PRACTICAL release")


def scan_or_materialize(
    source: Path,
    output: Path | None,
    *,
    expected_data_rows: int = EXPECTED_DATA_ROWS,
    expected_malformed_rows: int = EXPECTED_MALFORMED_ROWS,
    expected_field_count: int = EXPECTED_FIELD_COUNT,
    expected_non_rsid_rows: int = EXPECTED_NON_RSID_ROWS,
    expected_non_snp_rows: int = EXPECTED_NON_SNP_ROWS,
    expected_nonautosomal_rows: int = EXPECTED_NONAUTOSOMAL_ROWS,
    expected_retained_rows: int = EXPECTED_RETAINED_ROWS,
) -> dict[str, int]:
    """Audit the source and optionally project harmonization-eligible rows."""
    temporary: Path | None = None
    raw_output = None
    gzip_output = None
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary_name = tempfile.mkstemp(
            prefix=f".{output.name}.", suffix=".partial", dir=output.parent
        )
        os.close(fd)
        temporary = Path(temporary_name)
        raw_output = temporary.open("wb")
        gzip_output = gzip.GzipFile(
            filename="", mode="wb", compresslevel=6, mtime=0, fileobj=raw_output
        )

    source_rows = retained_rows = malformed_rows = 0
    non_rsid_rows = non_snp_rows = nonautosomal_rows = 0
    malformed_examples: list[dict[str, int]] = []
    try:
        with source.open("rb") as source_handle:
            header = source_handle.readline().rstrip(b"\r\n")
            if header != HEADER:
                fail(f"unexpected PRACTICAL header: {header[:200]!r}")
            if gzip_output is not None:
                gzip_output.write(OUTPUT_HEADER + b"\n")
            for line_number, line in enumerate(source_handle, start=2):
                source_rows += 1
                stripped = line.rstrip(b"\r\n")
                field_count = stripped.count(b"\t") + 1
                if field_count != expected_field_count:
                    malformed_rows += 1
                    if len(malformed_examples) < 10:
                        malformed_examples.append(
                            {"line_number": line_number, "field_count": field_count}
                        )
                    continue
                fields = stripped.split(b"\t")
                if RSID.fullmatch(fields[2]) is None:
                    non_rsid_rows += 1
                    continue
                if fields[5].upper() not in VALID_ALLELES or fields[6].upper() not in VALID_ALLELES:
                    non_snp_rows += 1
                    continue
                if fields[3] not in AUTOSOMES:
                    nonautosomal_rows += 1
                    continue
                retained_rows += 1
                if gzip_output is not None:
                    gzip_output.write(
                        b"\t".join(fields[index] for index in OUTPUT_COLUMN_INDEXES) + b"\n"
                    )
        if source_rows != expected_data_rows:
            fail(f"source has {source_rows:,} data rows, expected {expected_data_rows:,}")
        if malformed_rows != expected_malformed_rows:
            fail(
                f"source has {malformed_rows:,} malformed rows, "
                f"expected {expected_malformed_rows:,}"
            )
        expected_categories = {
            "non-rsID": (non_rsid_rows, expected_non_rsid_rows),
            "non-SNP": (non_snp_rows, expected_non_snp_rows),
            "non-autosomal": (nonautosomal_rows, expected_nonautosomal_rows),
            "retained": (retained_rows, expected_retained_rows),
        }
        for label, (observed, expected) in expected_categories.items():
            if observed != expected:
                fail(f"source has {observed:,} {label} rows, expected {expected:,}")
    except BaseException:
        if gzip_output is not None:
            gzip_output.close()
        if raw_output is not None:
            raw_output.close()
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
    else:
        if gzip_output is not None:
            gzip_output.close()
        if raw_output is not None:
            raw_output.close()

    counts = {
        "source_rows": source_rows,
        "retained_rows": retained_rows,
        "malformed_rows": malformed_rows,
        "non_rsid_rows": non_rsid_rows,
        "non_snp_rows": non_snp_rows,
        "nonautosomal_rows": nonautosomal_rows,
    }
    if output is not None and temporary is not None:
        os.replace(temporary, output)
        provenance = {
            "schema_version": "sleep-gwas-atlas.raw-materialization.v1",
            "source_path": str(source),
            "source_bytes": source.stat().st_size,
            "source_sha256": sha256(source),
            "output_path": str(output),
            "output_bytes": output.stat().st_size,
            "output_sha256": sha256(output),
            "expected_field_count": expected_field_count,
            "materialization_policy": (
                "drop_unrecoverable_physical_rows_then_project_explicit_rsid_"
                "single_base_autosomal_records"
            ),
            "counts": counts,
            "malformed_examples": malformed_examples,
        }
        provenance_path = Path(f"{output}.provenance.json")
        provenance_temp = provenance_path.with_name(f".{provenance_path.name}.partial")
        provenance_temp.write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(provenance_temp, provenance_path)
    return counts


def main() -> None:
    require_locked_traits({"prostate_cancer": "schumacher_2018_prostate_cancer"})
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if args.verify_only == (args.out is not None):
        fail("supply exactly one of --verify-only or --out")
    verify_archive(args.source)
    counts = scan_or_materialize(args.source, None if args.verify_only else args.out)
    action = "Verified" if args.verify_only else "Materialized"
    print(f"{action} PRACTICAL prostate-cancer source")
    print(f"  source rows: {counts['source_rows']:,}")
    print(f"  skipped unrecoverable malformed rows: {counts['malformed_rows']:,}")
    print(f"  skipped non-rsID rows: {counts['non_rsid_rows']:,}")
    print(f"  skipped non-SNP rows: {counts['non_snp_rows']:,}")
    print(f"  skipped non-autosomal rows: {counts['nonautosomal_rows']:,}")
    print(f"  retained harmonization-eligible rows: {counts['retained_rows']:,}")
    if args.out is not None:
        print(f"  output: {args.out}")
        print(f"  output SHA-256: {sha256(args.out)}")
        print(f"  provenance: {args.out}.provenance.json")


if __name__ == "__main__":
    main()
