"""Validated loading and application of checksum-sealed GRCh37 identity maps."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import re


HM3_SCHEMA_VERSION = "atlas.hm3-grch37-variant-map.v1"
DENSE_SCHEMA_VERSION = "atlas.grch37-variant-map.v2"
SCHEMA_VERSION = HM3_SCHEMA_VERSION
SCHEMA_SCOPES = {
    HM3_SCHEMA_VERSION: "HAPMAP3_ONLY",
    DENSE_SCHEMA_VERSION: "GENOME_WIDE_IMPUTED_VARIANT_IDENTITY",
}
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
    schema_version = provenance.get("schema_version")
    if schema_version not in SCHEMA_SCOPES:
        raise VariantMapError(
            f"unsupported variant-map schema: {schema_version!r}"
        )
    expected_scope = SCHEMA_SCOPES[schema_version]
    observed_scope = provenance.get("map_scope", expected_scope)
    if observed_scope != expected_scope:
        raise VariantMapError(
            f"variant-map scope is {observed_scope!r}, expected {expected_scope!r}"
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
    if schema_version == DENSE_SCHEMA_VERSION:
        validate_dense_policy_seal(path, provenance)
    return {**provenance, "map_scope": observed_scope}


def safe_relative(value, label: str) -> Path:
    relative = Path(str(value))
    if not str(value) or relative.is_absolute() or ".." in relative.parts:
        raise VariantMapError(f"unsafe {label} in dense variant-map provenance")
    return relative


def root_from_relative_artifact(path: Path, relative: Path) -> Path:
    absolute = path.resolve()
    root = absolute
    for _ in relative.parts:
        root = root.parent
    if (root / relative).resolve() != absolute:
        raise VariantMapError("dense variant map is not at its provenance-declared path")
    return root


def validate_dense_policy_seal(path: Path, provenance: dict) -> None:
    """Bind a genome-wide map to its live policy and builder, not a QC label."""
    map_relative = safe_relative(provenance.get("map_path", ""), "map path")
    policy_relative = safe_relative(provenance.get("policy_path", ""), "policy path")
    builder_relative = safe_relative(provenance.get("builder_path", ""), "builder path")
    root = root_from_relative_artifact(path, map_relative)
    policy_path, builder_path = root / policy_relative, root / builder_relative
    if not policy_path.is_file() or not builder_path.is_file():
        raise VariantMapError("dense variant-map policy or builder is absent")
    try:
        policy = json.loads(policy_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise VariantMapError(f"invalid dense variant-map policy: {error}") from error
    source = provenance.get("source", {})
    if not isinstance(source, dict):
        raise VariantMapError("dense variant-map source provenance is invalid")
    required_matches = (
        provenance.get("analysis_panel") == "atlas-v1.0",
        provenance.get("policy_sha256") == sha256(policy_path),
        provenance.get("builder_sha256") == sha256(builder_path),
        provenance.get("builder_path") == policy.get("builder_path"),
        provenance.get("builder_sha256") == policy.get("builder_sha256"),
        provenance.get("schema_version") == policy.get("map_schema_version"),
        provenance.get("map_scope") == policy.get("map_scope"),
        provenance.get("map_path") == policy.get("map_path"),
        str(Path(f"{provenance.get('map_path', '')}.provenance.json"))
        == policy.get("map_provenance_path"),
        provenance.get("resource") == policy.get("resource"),
        provenance.get("identity_rule") == policy.get("identity_rule"),
        provenance.get("role_limit") == policy.get("role_limit"),
        source.get("path") == policy.get("source_path"),
        source.get("url") == policy.get("source_url"),
        source.get("bytes") == policy.get("source_bytes"),
        source.get("md5") == policy.get("source_md5"),
        source.get("s3_version_id") == policy.get("source_s3_version_id"),
        isinstance(source.get("sha256"), str)
        and bool(re.fullmatch(r"[0-9a-f]{64}", source["sha256"])),
        isinstance(provenance.get("mapped_rows"), int)
        and provenance["mapped_rows"] >= int(policy.get("minimum_mapped_rows", 0)),
    )
    if not all(required_matches):
        raise VariantMapError("dense variant map differs from its live policy/builder seal")


def load_variant_map(
    map_path: str | Path,
    strategy: str,
    expected_sha256: str | None = None,
    expected_bytes: int | None = None,
    required_keys: set | None = None,
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
                wanted = required_keys is None or snp in required_keys
                if snp in by_rsid:
                    raise VariantMapError(f"duplicate rsID in variant map: {snp}")
                if wanted:
                    by_rsid[snp] = (chromosome, position, key)
            else:
                locus = (chromosome, position, *key)
                if required_keys is not None and locus not in required_keys:
                    row_count += 1
                    continue
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


def required_coordinate_keys(chromosome, position, a1, a2):
    """Return conservative direct/complement lookup keys for one source row."""
    try:
        chromosome, position = int(chromosome), int(position)
    except (TypeError, ValueError, OverflowError):
        return ()
    keys = []
    for alleles in dict.fromkeys((allele_key(a1, a2), complement_key(a1, a2))):
        if alleles is not None:
            keys.append((chromosome, position, *alleles))
    return tuple(keys)


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
