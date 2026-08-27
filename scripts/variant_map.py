"""Validated loading and application of the pinned GRCh37 HapMap3 map."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import re


SCHEMA_VERSION = "atlas.hm3-grch37-variant-map.v1"
MAP_HEADER = ["SNP", "CHR", "BP", "A1", "A2"]
VALID_ALLELES = {"A", "C", "G", "T"}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
RSID = re.compile(r"rs[0-9]+$", re.IGNORECASE)
STRATEGIES = {"BY_COORD_ALLELES", "BY_RSID_ALLELES"}
CONFLICT = object()


class VariantMapError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def allele_key(a1: str, a2: str) -> tuple[str, str] | None:
    if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES or a1 == a2:
        return None
    return tuple(sorted((a1, a2)))


def complement_key(a1: str, a2: str) -> tuple[str, str] | None:
    if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES:
        return None
    return allele_key(COMPLEMENT[a1], COMPLEMENT[a2])


def validate_provenance(
    map_path: str | Path,
    expected_sha256: str | None = None,
    expected_bytes: int | None = None,
) -> dict:
    path = Path(map_path)
    provenance_path = Path(f"{path}.provenance.json")
    if not path.is_file():
        raise VariantMapError(f"variant map not found: {path}")
    if not provenance_path.is_file():
        raise VariantMapError(f"variant-map provenance not found: {provenance_path}")
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise VariantMapError(f"invalid variant-map provenance: {error}") from error
    if provenance.get("schema_version") != SCHEMA_VERSION:
        raise VariantMapError(
            f"variant-map schema is {provenance.get('schema_version')!r}, expected {SCHEMA_VERSION}"
        )
    if provenance.get("genome_build") != "GRCh37/hg19":
        raise VariantMapError("variant map does not declare GRCh37/hg19")
    actual_hash = sha256(path)
    if expected_sha256 is not None and actual_hash != expected_sha256:
        raise VariantMapError(
            f"variant-map SHA-256 differs from the registered map: "
            f"expected {expected_sha256}, got {actual_hash}"
        )
    if expected_bytes is not None and path.stat().st_size != expected_bytes:
        raise VariantMapError(
            f"variant-map bytes differ from the registered map: "
            f"expected {expected_bytes}, got {path.stat().st_size}"
        )
    if provenance.get("map_sha256") != actual_hash:
        raise VariantMapError(
            f"variant-map SHA-256 mismatch: expected {provenance.get('map_sha256')}, got {actual_hash}"
        )
    if provenance.get("map_bytes") != path.stat().st_size:
        raise VariantMapError("variant-map byte count does not match provenance")
    return provenance


def load_variant_map(
    map_path: str | Path,
    strategy: str,
    expected_sha256: str | None = None,
    expected_bytes: int | None = None,
):
    if strategy not in STRATEGIES:
        raise VariantMapError(f"invalid variant-map strategy: {strategy}")
    provenance = validate_provenance(map_path, expected_sha256, expected_bytes)
    by_rsid = {} if strategy == "BY_RSID_ALLELES" else None
    by_locus = {} if strategy == "BY_COORD_ALLELES" else None
    row_count = 0
    with gzip.open(map_path, "rt", encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames != MAP_HEADER:
            raise VariantMapError(f"unexpected variant-map header: {rows.fieldnames}")
        for line_number, row in enumerate(rows, start=2):
            snp = row["SNP"].strip().lower()
            a1, a2 = row["A1"].strip().upper(), row["A2"].strip().upper()
            try:
                chromosome, position = int(row["CHR"]), int(row["BP"])
            except ValueError as error:
                raise VariantMapError(f"invalid map coordinate at line {line_number}") from error
            key = allele_key(a1, a2)
            if (
                not RSID.fullmatch(snp)
                or chromosome < 1
                or chromosome > 22
                or position <= 0
                or key is None
            ):
                raise VariantMapError(f"invalid variant-map row at line {line_number}")
            if by_rsid is not None:
                if snp in by_rsid:
                    raise VariantMapError(f"duplicate rsID in variant map: {snp}")
                by_rsid[snp] = (chromosome, position, key)
            else:
                locus = (chromosome, position, *key)
                previous = by_locus.get(locus)
                if previous is None and locus not in by_locus:
                    by_locus[locus] = snp
                elif previous != snp:
                    by_locus[locus] = CONFLICT
            row_count += 1
    if row_count != provenance.get("mapped_rows"):
        raise VariantMapError(
            f"variant-map row count {row_count:,} does not match provenance "
            f"{provenance.get('mapped_rows')!r}"
        )
    return (by_rsid if by_rsid is not None else by_locus), provenance


def coordinate_match(index, chromosome, position, a1, a2):
    if chromosome is None or position is None:
        return "unmatched_reference", None
    direct = allele_key(a1, a2)
    reverse = complement_key(a1, a2)
    if direct is None:
        return "allele_conflict", None
    try:
        chromosome, position = int(chromosome), int(position)
    except (TypeError, ValueError, OverflowError):
        return "unmatched_reference", None
    candidates = []
    conflict = False
    for key in dict.fromkeys((direct, reverse)):
        hit = index.get((chromosome, position, *key))
        if hit is CONFLICT:
            conflict = True
        elif hit is not None:
            candidates.append(hit)
    unique = set(candidates)
    if conflict or len(unique) > 1:
        return "ambiguous_reference", None
    if not unique:
        return "unmatched_reference", None
    return "mapped", unique.pop()


def rsid_match(index, snp, chromosome, position, a1, a2):
    snp = str(snp).strip().lower()
    if not RSID.fullmatch(snp):
        return "unmatched_reference", None
    hit = index.get(snp)
    if hit is None:
        return "unmatched_reference", None
    ref_chromosome, ref_position, ref_alleles = hit
    direct = allele_key(a1, a2)
    reverse = complement_key(a1, a2)
    if direct is None or ref_alleles not in {direct, reverse}:
        return "allele_conflict", None
    if chromosome is not None and position is not None:
        try:
            existing = (int(chromosome), int(position))
        except (TypeError, ValueError, OverflowError):
            existing = None
        if existing is not None and existing != (ref_chromosome, ref_position):
            return "coordinate_conflict", None
    return "mapped", (snp, ref_chromosome, ref_position)
