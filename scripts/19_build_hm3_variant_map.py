#!/usr/bin/env python3
"""Build an audited HapMap3 rsID/GRCh37/allele map from the pinned LDSC bundle.

The bundled ``w_hm3.snplist`` supplies rsIDs and alleles; the 22
``*.l2.ldscore.gz`` files supply GRCh37 chromosome and position. Only their
intersection is emitted. The output is deterministic and accompanied by a
JSON provenance record that pins every input and the output SHA-256.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tarfile
import tempfile


SCHEMA_VERSION = "atlas.hm3-grch37-variant-map.v1"
VALID_ALLELES = {"A", "C", "G", "T"}
RSID = re.compile(r"rs[0-9]+$", re.IGNORECASE)
DEFAULT_CHROMOSOMES = tuple(range(1, 23))
DEFAULT_EXPECTED_HM3 = 1_217_311
DEFAULT_MIN_COVERAGE = 0.95


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def file_hash(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_chromosomes(value: str) -> tuple[int, ...]:
    chromosomes = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            chromosomes.extend(range(int(start), int(end) + 1))
        else:
            chromosomes.append(int(token))
    result = tuple(dict.fromkeys(chromosomes))
    if not result or any(chromosome < 1 or chromosome > 22 for chromosome in result):
        fail(f"invalid autosomal chromosome selection: {value!r}")
    return result


class ReferenceReader:
    def __init__(self, archive: Path | None, directory: Path | None):
        self.archive = archive
        self.directory = directory
        self.tar = tarfile.open(archive, "r:gz") if archive else None

    def close(self) -> None:
        if self.tar is not None:
            self.tar.close()

    def binary(self, relative: str):
        if self.tar is not None:
            member = f"eur_w_ld_chr/{relative}"
            try:
                handle = self.tar.extractfile(member)
            except KeyError:
                handle = None
            if handle is None:
                fail(f"reference archive is missing {member}")
            return handle
        path = self.directory / relative
        if not path.is_file():
            fail(f"reference directory is missing {path}")
        return path.open("rb")


def read_hm3(reader: ReferenceReader) -> dict[str, tuple[str, str]]:
    with reader.binary("w_hm3.snplist") as binary:
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as text_handle:
            rows = csv.DictReader(text_handle, delimiter="\t")
            if rows.fieldnames != ["SNP", "A1", "A2"]:
                fail(f"unexpected w_hm3 header: {rows.fieldnames}")
            hm3: dict[str, tuple[str, str]] = {}
            for line_number, row in enumerate(rows, start=2):
                snp = row["SNP"].strip().lower()
                a1, a2 = row["A1"].strip().upper(), row["A2"].strip().upper()
                if not RSID.fullmatch(snp):
                    fail(f"invalid rsID in w_hm3 at line {line_number}: {snp!r}")
                if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES or a1 == a2:
                    fail(f"invalid SNP alleles in w_hm3 at line {line_number}: {a1}/{a2}")
                if snp in hm3:
                    fail(f"duplicate rsID in w_hm3: {snp}")
                hm3[snp] = (a1, a2)
    if not hm3:
        fail("w_hm3 contains no variants")
    return hm3


def read_coordinates(
    reader: ReferenceReader,
    chromosomes: tuple[int, ...],
) -> dict[str, tuple[int, int]]:
    coordinates: dict[str, tuple[int, int]] = {}
    required = {"CHR", "SNP", "BP"}
    for expected_chromosome in chromosomes:
        with reader.binary(f"{expected_chromosome}.l2.ldscore.gz") as compressed:
            with gzip.GzipFile(fileobj=compressed, mode="rb") as binary:
                with io.TextIOWrapper(binary, encoding="utf-8", newline="") as text_handle:
                    rows = csv.DictReader(text_handle, delimiter="\t")
                    if not required.issubset(rows.fieldnames or []):
                        fail(
                            f"chromosome {expected_chromosome} LD-score header lacks "
                            f"{sorted(required)}: {rows.fieldnames}"
                        )
                    for line_number, row in enumerate(rows, start=2):
                        snp = row["SNP"].strip().lower()
                        try:
                            chromosome = int(row["CHR"])
                            position = int(row["BP"])
                        except ValueError:
                            fail(
                                f"invalid coordinate in chromosome {expected_chromosome} "
                                f"LD-score file at line {line_number}"
                            )
                        if chromosome != expected_chromosome or position <= 0:
                            fail(
                                f"unexpected coordinate for {snp}: {chromosome}:{position}; "
                                f"expected chromosome {expected_chromosome}"
                            )
                        if not RSID.fullmatch(snp):
                            fail(f"invalid LD-score rsID at line {line_number}: {snp!r}")
                        previous = coordinates.get(snp)
                        if previous is not None and previous != (chromosome, position):
                            fail(f"rsID maps to multiple GRCh37 coordinates: {snp}")
                        if previous is not None:
                            fail(f"duplicate LD-score rsID: {snp}")
                        coordinates[snp] = (chromosome, position)
    if not coordinates:
        fail("LD-score files contain no coordinates")
    return coordinates


def input_provenance(
    archive: Path | None,
    directory: Path | None,
    chromosomes: tuple[int, ...],
) -> dict:
    if archive is not None:
        return {
            "mode": "archive",
            "path": str(archive.resolve()),
            "bytes": archive.stat().st_size,
            "md5": file_hash(archive, "md5"),
            "sha256": file_hash(archive),
        }
    paths = [directory / "w_hm3.snplist"] + [
        directory / f"{chromosome}.l2.ldscore.gz" for chromosome in chromosomes
    ]
    files = []
    combined = hashlib.sha256()
    for path in paths:
        digest = file_hash(path)
        relative = str(path.relative_to(directory))
        combined.update(f"{relative}\t{digest}\n".encode())
        files.append({"path": relative, "bytes": path.stat().st_size, "sha256": digest})
    return {
        "mode": "directory",
        "path": str(directory.resolve()),
        "combined_sha256": combined.hexdigest(),
        "files": files,
    }


def write_map(path: Path, rows: list[tuple[str, int, int, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent
    )
    os.close(file_descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as raw:
            # Do not let GzipFile embed the randomized temporary filename in
            # the gzip header. Together with mtime=0 this makes the map bytes
            # reproducible across output paths and repeated setup runs.
            with gzip.GzipFile(
                filename="", fileobj=raw, mode="wb", mtime=0
            ) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                    writer = csv.writer(text_handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["SNP", "CHR", "BP", "A1", "A2"])
                    writer.writerows(rows)
        os.replace(temporary, path)
        path.chmod(0o644)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive", type=Path)
    source.add_argument("--reference-dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--chromosomes", default="1-22")
    parser.add_argument("--expected-hm3-count", type=int, default=DEFAULT_EXPECTED_HM3)
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE)
    args = parser.parse_args()

    if args.archive is not None and not args.archive.is_file():
        fail(f"reference archive not found: {args.archive}")
    if args.reference_dir is not None and not args.reference_dir.is_dir():
        fail(f"reference directory not found: {args.reference_dir}")
    if not 0 < args.min_coverage <= 1:
        fail("--min-coverage must be in (0, 1]")
    chromosomes = parse_chromosomes(args.chromosomes)

    reader = ReferenceReader(args.archive, args.reference_dir)
    try:
        hm3 = read_hm3(reader)
        coordinates = read_coordinates(reader, chromosomes)
    finally:
        reader.close()

    if args.expected_hm3_count and len(hm3) != args.expected_hm3_count:
        fail(
            f"expected {args.expected_hm3_count:,} w_hm3 variants, found {len(hm3):,}"
        )
    mapped_ids = hm3.keys() & coordinates.keys()
    coverage = len(mapped_ids) / len(hm3)
    if coverage < args.min_coverage:
        fail(
            f"only {len(mapped_ids):,}/{len(hm3):,} w_hm3 variants have coordinates "
            f"({coverage:.2%}); required {args.min_coverage:.2%}"
        )
    rows = [
        (snp, coordinates[snp][0], coordinates[snp][1], hm3[snp][0], hm3[snp][1])
        for snp in mapped_ids
    ]
    rows.sort(key=lambda row: (row[1], row[2], int(row[0][2:])))
    write_map(args.out, rows)

    provenance_path = Path(f"{args.out}.provenance.json")
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "resource": "eur_w_ld_chr bundled w_hm3 + LD-score coordinates",
        "genome_build": "GRCh37/hg19",
        "chromosomes": list(chromosomes),
        "source": input_provenance(args.archive, args.reference_dir, chromosomes),
        "hm3_rows": len(hm3),
        "ldscore_rows": len(coordinates),
        "mapped_rows": len(rows),
        "hm3_rows_without_eur_ldscore_coordinate": len(hm3) - len(rows),
        "coverage": coverage,
        "map_path": str(args.out),
        "map_bytes": args.out.stat().st_size,
        "map_sha256": file_hash(args.out),
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    provenance_path.chmod(0o644)
    print(
        f"Wrote {len(rows):,} GRCh37 HapMap3 mappings to {args.out} "
        f"({coverage:.2%} of w_hm3)"
    )
    print(f"Wrote provenance: {provenance_path}")


if __name__ == "__main__":
    main()
