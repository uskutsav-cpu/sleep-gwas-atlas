#!/usr/bin/env python3
"""Materialize and verify the checksum-sealed Track B EUR PLINK LD reference."""

from __future__ import annotations

import argparse
import ctypes
import ctypes.util
import csv
import errno
import hashlib
import importlib.util
import io
import json
import math
import os
import resource
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from typing import Any


PRODUCTION_POLICY = Path("config/track_b_pleiotropy_policy.json")
CONTRACT_SCRIPT = Path("scripts/123_track_b_pleiotropy_contract.py")
SCRIPT = Path("scripts/127_prepare_track_b_pleiotropy_ld.py")
COMPONENT_ORDER = ("bed", "bim", "fam")
MANIFEST_FIELDS = [
    "component_id", "path", "bytes", "sha256", "records", "validation_status",
]
TEST_SCHEMA = "sleep-atlas-track-b-pleiotropy-ld-fixture.1"
PRODUCTION_SCHEMA = "sleep-atlas-track-b-pleiotropy-policy.1"


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def safe_path(root: Path, value: str | Path) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        fail(f"unsafe repository-relative path: {value}")
    return root / relative


def relative(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError as error:
        fail(f"artifact is outside repository root: {path}")
        raise AssertionError from error


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        fail(f"unreadable JSON artifact {path}: {error}")
    if not isinstance(value, dict):
        fail(f"JSON artifact must contain one object: {path}")
    return value


def atomic_text(path: Path, payload: str) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def load_contract(root: Path):
    path = safe_path(root, CONTRACT_SCRIPT)
    spec = importlib.util.spec_from_file_location("track_b_pleiotropy_contract_for_ld", path)
    if spec is None or spec.loader is None:
        fail("could not load the Track B pleiotropy contract")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_reference_policy(ld: dict[str, Any], *, test_fixture: bool) -> None:
    required_strings = [
        "source_archive_path", "source_archive_sha256", "source_release",
        "materialized_directory", "materialized_prefix", "materialized_manifest",
        "materialized_provenance", "staging_root", "materializer_script",
        "plink_path", "plink_sha256", "plink_version", "ancestry", "build",
        "raw_reference_fidelity_rule", "analysis_chromosome_rule",
        "materialization_transform_rule",
    ]
    if any(not isinstance(ld.get(field), str) or not ld[field] for field in required_strings):
        fail("LD reference policy lacks a required non-empty string")
    for field in ("source_archive_sha256", "plink_sha256"):
        if len(ld[field]) != 64 or any(character not in "0123456789abcdef" for character in ld[field]):
            fail(f"invalid SHA-256 in LD reference policy: {field}")
    for field in (
        "source_archive_expected_bytes", "sample_count", "variant_count",
        "expected_autosomal_variant_count", "staging_safety_bytes", "minimum_free_storage_bytes",
    ):
        if type(ld.get(field)) is not int or ld[field] <= 0:
            fail(f"invalid positive integer in LD reference policy: {field}")
    members = ld.get("source_members")
    if not isinstance(members, dict) or list(members) != list(COMPONENT_ORDER):
        fail("LD reference policy must pin BED/BIM/FAM in exact order")
    for component in COMPONENT_ORDER:
        member = members[component]
        if (
            not isinstance(member, dict)
            or not isinstance(member.get("name"), str)
            or member["name"] != f"g1000_eur.{component}"
            or type(member.get("expected_bytes")) is not int
            or member["expected_bytes"] <= 0
            or type(member.get("records")) is not int
            or member["records"] <= 0
            or not isinstance(member.get("sha256"), str)
            or len(member["sha256"]) != 64
            or any(character not in "0123456789abcdef" for character in member["sha256"])
        ):
            fail(f"incomplete LD member pin: {component}")
    if (
        members["bed"]["records"] != ld["variant_count"]
        or members["bim"]["records"] != ld["variant_count"]
        or members["fam"]["records"] != ld["sample_count"]
    ):
        fail("PLINK member record counts differ from reference totals")
    bytes_per_variant = math.ceil(ld["sample_count"] / 4)
    expected_bed_bytes = 3 + bytes_per_variant * ld["variant_count"]
    if members["bed"]["expected_bytes"] != expected_bed_bytes:
        fail("BED byte count is inconsistent with SNP-major PLINK dimensions")
    allowed = ld.get("archive_allowed_members")
    member_names = [members[component]["name"] for component in COMPONENT_ORDER]
    if (
        not isinstance(allowed, list)
        or any(
            not isinstance(value, str) or not value
            or Path(value).is_absolute() or ".." in Path(value).parts
            for value in allowed
        )
        or len(allowed) != len(set(allowed))
        or not set(member_names).issubset(allowed)
    ):
        fail("archive member allowlist is invalid or incomplete")
    reference_chromosomes = ld.get("reference_chromosomes")
    analysis_autosomes = ld.get("analysis_autosomes")
    nonanalysis = ld.get("nonanalysis_chromosomes_retained_for_source_fidelity")
    counts_raw = ld.get("expected_reference_chromosome_variant_counts")
    lengths_raw = ld.get("grch37_chromosome_lengths")
    if (
        not isinstance(reference_chromosomes, list)
        or not isinstance(analysis_autosomes, list)
        or not isinstance(nonanalysis, list)
        or any(type(value) is not int or value <= 0 for value in reference_chromosomes)
        or any(type(value) is not int or value <= 0 for value in analysis_autosomes)
        or any(type(value) is not int or value <= 0 for value in nonanalysis)
        or not isinstance(counts_raw, dict)
        or not isinstance(lengths_raw, dict)
        or len(reference_chromosomes) != len(set(reference_chromosomes))
        or not set(analysis_autosomes).isdisjoint(nonanalysis)
        or set(reference_chromosomes) != set(analysis_autosomes) | set(nonanalysis)
        or set(counts_raw) != {str(value) for value in reference_chromosomes}
        or set(lengths_raw) != {str(value) for value in reference_chromosomes}
    ):
        fail("reference/analysis chromosome policy is inconsistent")
    try:
        counts = {int(key): int(value) for key, value in counts_raw.items()}
        lengths = {int(key): int(value) for key, value in lengths_raw.items()}
    except (TypeError, ValueError) as error:
        fail(f"invalid chromosome count/length policy: {error}")
    if (
        any(value <= 0 for value in counts.values())
        or any(value <= 0 for value in lengths.values())
        or sum(counts.values()) != ld["variant_count"]
        or sum(counts[chromosome] for chromosome in analysis_autosomes)
        != ld["expected_autosomal_variant_count"]
    ):
        fail("chromosome counts do not reconcile to reference/autosomal totals")
    extracted_bytes = sum(members[component]["expected_bytes"] for component in COMPONENT_ORDER)
    if ld["minimum_free_storage_bytes"] < extracted_bytes + ld["staging_safety_bytes"]:
        fail("LD disk gate is smaller than extracted payload plus staging safety")
    if ld.get("materialization_resource_metrics_required") != [
        "staging_extraction_validation_wall_seconds", "peak_rss_bytes",
        "streaming_chunk_bytes", "extracted_uncompressed_bytes",
    ]:
        fail("LD materialization resource-metric schema drifted")
    markers = ld.get("conflicting_activity_markers")
    if (
        not isinstance(markers, list)
        or any(
            not isinstance(marker, str) or not marker
            or Path(marker).is_absolute() or ".." in Path(marker).parts
            for marker in markers
        )
        or len(markers) != len(set(markers))
    ):
        fail("LD conflicting-activity marker policy is invalid")
    destination = Path(ld["materialized_directory"])
    prefix = Path(ld["materialized_prefix"])
    manifest = Path(ld["materialized_manifest"])
    provenance = Path(ld["materialized_provenance"])
    if (
        prefix.parent != destination or manifest.parent != destination
        or provenance.parent != destination or prefix.name != "g1000_eur"
        or manifest.name != "g1000_eur.manifest.tsv"
        or provenance.name != "g1000_eur.provenance.json"
        or ld["materializer_script"] != str(SCRIPT)
    ):
        fail("permanent LD directory/prefix/manifest/provenance paths are inconsistent")
    if not test_fixture and (
        ld.get("reference_chromosomes") != list(range(1, 24))
        or ld.get("analysis_autosomes") != list(range(1, 23))
        or ld.get("nonanalysis_chromosomes_retained_for_source_fidelity") != [23]
        or ld.get("variant_count") != 22_665_064
        or ld.get("expected_autosomal_variant_count") != 22_132_657
        or ld.get("sample_count") != 503
        or ld.get("build") != "GRCh37/hg19"
        or ld.get("ancestry") != "EUR"
        or "chr23/X" not in ld.get("analysis_chromosome_rule", "")
        or "no LD-block splitting" not in ld.get("materialization_transform_rule", "")
        or "SNP reduction" not in ld.get("materialization_transform_rule", "")
        or ld.get("materialization_resource_metrics_required") != [
            "staging_extraction_validation_wall_seconds", "peak_rss_bytes",
            "streaming_chunk_bytes", "extracted_uncompressed_bytes",
        ]
    ):
        fail("production raw-reference/autosome split drifted")


def materializer_path(root: Path, *, test_fixture: bool) -> Path:
    observed = Path(__file__).resolve()
    expected = safe_path(root, SCRIPT).resolve()
    if not test_fixture and observed != expected:
        fail("production LD materializer is not running from its policy-pinned repository path")
    return observed


def load_policy_context(
    root: Path, policy_relative: Path, *, test_fixture: bool,
) -> tuple[dict[str, Any], dict[str, Any], Path, str]:
    policy_path = safe_path(root, policy_relative)
    policy = read_json(policy_path)
    if test_fixture:
        if policy.get("schema_version") != TEST_SCHEMA or policy.get("analysis_id") != "track-b-test-pleiotropy-ld":
            fail("--test-fixture requires the exact test-only LD policy schema and analysis ID")
        contract_lock_sha256 = "NOT_APPLICABLE_TEST_FIXTURE"
    else:
        if policy_relative != PRODUCTION_POLICY or policy.get("schema_version") != PRODUCTION_SCHEMA:
            fail("production LD operations require config/track_b_pleiotropy_policy.json")
        contract = load_contract(root)
        live_policy = contract.validate_policy(root)
        if live_policy != policy:
            fail("live Track B pleiotropy policy differs from its contract loader")
        contract.verify_contract(root)
        contract_lock_sha256 = sha256(root / contract.CONTRACT_LOCK)
    ld = policy.get("ld_reference_policy")
    if not isinstance(ld, dict):
        fail("policy lacks an LD reference object")
    validate_reference_policy(ld, test_fixture=test_fixture)
    return policy, ld, policy_path, contract_lock_sha256


def plink_identity(root: Path, ld: dict[str, Any]) -> dict[str, Any]:
    path = safe_path(root, ld["plink_path"])
    if not path.is_file() or not os.access(path, os.X_OK):
        fail(f"BLOCKED_BY_SOFTWARE pinned PLINK is absent or not executable: {path}")
    observed_sha256 = sha256(path)
    if observed_sha256 != ld["plink_sha256"]:
        fail("BLOCKED_BY_SOFTWARE pinned PLINK SHA-256 drifted")
    result = subprocess.run([str(path), "--version"], check=False, capture_output=True, text=True)
    version = (result.stdout + result.stderr).strip()
    if result.returncode or not version.startswith(ld["plink_version"]):
        fail(f"BLOCKED_BY_SOFTWARE pinned PLINK version drifted: {version!r}")
    return {
        "path": ld["plink_path"], "bytes": path.stat().st_size,
        "sha256": observed_sha256, "version": version,
    }


def archive_identity(root: Path, ld: dict[str, Any]) -> tuple[Path, dict[str, zipfile.ZipInfo]]:
    path = safe_path(root, ld["source_archive_path"])
    if (
        not path.is_file() or path.stat().st_size != ld["source_archive_expected_bytes"]
        or sha256(path) != ld["source_archive_sha256"]
    ):
        fail("BLOCKED_BY_DATA Track B EUR LD source archive size or SHA-256 drifted")
    try:
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
    except (OSError, zipfile.BadZipFile) as error:
        fail(f"BLOCKED_BY_DATA unreadable Track B EUR LD ZIP: {error}")
    names = [info.filename for info in infos]
    if len(names) != len(set(names)) or set(names) != set(ld["archive_allowed_members"]):
        fail("BLOCKED_BY_DATA Track B EUR LD ZIP member family drifted or contains duplicates")
    by_name = {info.filename: info for info in infos}
    for component in COMPONENT_ORDER:
        member = ld["source_members"][component]
        info = by_name[member["name"]]
        if info.is_dir() or info.flag_bits & 0x1 or info.file_size != member["expected_bytes"]:
            fail(f"BLOCKED_BY_DATA invalid/encrypted/wrong-size ZIP member: {member['name']}")
    return path, by_name


def conflict_markers(root: Path, ld: dict[str, Any]) -> list[Path]:
    return [
        safe_path(root, marker) for marker in ld.get("conflicting_activity_markers", [])
        if safe_path(root, marker).exists()
    ]


def reference_paths(root: Path, ld: dict[str, Any]) -> dict[str, Path]:
    directory = safe_path(root, ld["materialized_directory"])
    paths = {
        component: directory / ld["source_members"][component]["name"]
        for component in COMPONENT_ORDER
    }
    paths["manifest"] = safe_path(root, ld["materialized_manifest"])
    paths["provenance"] = safe_path(root, ld["materialized_provenance"])
    return paths


def expected_destination_names(ld: dict[str, Any]) -> set[str]:
    return {
        *(ld["source_members"][component]["name"] for component in COMPONENT_ORDER),
        Path(ld["materialized_manifest"]).name,
        Path(ld["materialized_provenance"]).name,
    }


def destination_state(root: Path, ld: dict[str, Any]) -> str:
    directory = safe_path(root, ld["materialized_directory"])
    if not os.path.lexists(directory):
        return "ABSENT"
    if directory.is_symlink() or not directory.is_dir():
        return "PARTIAL_OR_TAMPERED"
    children = list(directory.iterdir())
    observed = {path.name for path in children}
    return (
        "COMPLETE_SHAPE"
        if observed == expected_destination_names(ld)
        and all(path.is_file() and not path.is_symlink() for path in children)
        else "PARTIAL_OR_TAMPERED"
    )


def fingerprint(
    root: Path, policy: dict[str, Any], ld: dict[str, Any], policy_path: Path,
    contract_lock_sha256: str, *, test_fixture: bool,
) -> str:
    payload = {
        "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path),
        "contract_lock_sha256": contract_lock_sha256,
        "materializer_sha256": sha256(materializer_path(root, test_fixture=test_fixture)),
        "archive": {
            "path": ld["source_archive_path"], "bytes": ld["source_archive_expected_bytes"],
            "sha256": ld["source_archive_sha256"], "members": ld["source_members"],
        },
        "destination": ld["materialized_directory"],
        "plink_sha256": ld["plink_sha256"],
        "sample_count": ld["sample_count"], "variant_count": ld["variant_count"],
        "analysis_autosomes": ld["analysis_autosomes"],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def staging_directory(root: Path, ld: dict[str, Any], run_fingerprint: str) -> Path:
    if len(run_fingerprint) != 64 or any(value not in "0123456789abcdef" for value in run_fingerprint):
        fail("invalid LD materialization fingerprint")
    return safe_path(root, ld["staging_root"]) / f"g1000_eur_{run_fingerprint}"


def preflight(
    root: Path, policy: dict[str, Any], ld: dict[str, Any], policy_path: Path,
    contract_lock_sha256: str, *, test_fixture: bool,
) -> dict[str, Any]:
    state = destination_state(root, ld)
    if state == "PARTIAL_OR_TAMPERED":
        fail("BLOCKED_BY_DATA partial/preexisting Track B LD destination; overwrite is forbidden")
    run_fingerprint = fingerprint(
        root, policy, ld, policy_path, contract_lock_sha256, test_fixture=test_fixture,
    )
    stage = staging_directory(root, ld, run_fingerprint)
    if state == "ABSENT" and stage.exists():
        fail("BLOCKED_BY_DATA fingerprinted Track B LD staging directory already exists")
    if state == "ABSENT":
        conflicts = conflict_markers(root, ld)
        if conflicts:
            fail(
                "BLOCKED_ACTIVE_LAVA_DOWNLOAD conflicting marker(s): "
                + ", ".join(relative(root, path) for path in conflicts)
            )
        free_bytes = shutil.disk_usage(root).free
        if free_bytes < ld["minimum_free_storage_bytes"]:
            fail(
                "BLOCKED_BY_COMPUTE insufficient disk for Track B LD staging: "
                f"observed={free_bytes} required={ld['minimum_free_storage_bytes']}"
            )
    else:
        free_bytes = shutil.disk_usage(root).free
    archive_path, archive_infos = archive_identity(root, ld)
    plink = plink_identity(root, ld)
    return {
        "state": state, "fingerprint": run_fingerprint,
        "staging_directory": relative(root, stage), "free_bytes": free_bytes,
        "required_free_bytes": ld["minimum_free_storage_bytes"],
        "archive_path": relative(root, archive_path), "archive_infos": archive_infos,
        "plink": plink,
    }


def copy_zip_member(
    archive: zipfile.ZipFile, info: zipfile.ZipInfo, destination: Path,
    expected_size: int, expected_sha256: str,
) -> dict[str, Any]:
    digest = hashlib.sha256()
    written = 0
    try:
        with archive.open(info, "r") as source, destination.open("xb") as target:
            while True:
                block = source.read(8 * 1024 * 1024)
                if not block:
                    break
                target.write(block)
                digest.update(block)
                written += len(block)
            target.flush()
            os.fsync(target.fileno())
    except (OSError, EOFError, zipfile.BadZipFile) as error:
        fail(f"LD member extraction failed for {info.filename}: {error}")
    observed_sha256 = digest.hexdigest()
    if written != expected_size or observed_sha256 != expected_sha256:
        fail(
            f"extracted LD member identity drifted: {info.filename} "
            f"bytes={written} sha256={observed_sha256}"
        )
    return {"bytes": written, "sha256": observed_sha256}


def validate_bed(path: Path, sample_count: int, variant_count: int) -> dict[str, Any]:
    expected_bytes = 3 + math.ceil(sample_count / 4) * variant_count
    if path.stat().st_size != expected_bytes:
        fail("PLINK BED size is inconsistent with the exact sample/variant dimensions")
    with path.open("rb") as handle:
        magic = handle.read(3)
    if magic != b"\x6c\x1b\x01":
        fail("PLINK BED lacks the SNP-major magic header")
    return {
        "magic_hex": magic.hex(), "variant_major": True,
        "bytes_per_variant": math.ceil(sample_count / 4),
    }


def validate_fam(path: Path, sample_count: int) -> dict[str, Any]:
    observed = 0
    identities: set[tuple[bytes, bytes]] = set()
    try:
        with path.open("rb") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.split()
                if len(fields) != 6:
                    fail(f"FAM schema drift at line {line_number}")
                identity = (fields[0], fields[1])
                if identity in identities:
                    fail(f"duplicate FAM FID/IID at line {line_number}")
                identities.add(identity)
                if fields[4] not in {b"0", b"1", b"2"}:
                    fail(f"invalid FAM sex code at line {line_number}")
                try:
                    phenotype = float(fields[5])
                except ValueError:
                    fail(f"invalid FAM phenotype at line {line_number}")
                if not math.isfinite(phenotype):
                    fail(f"non-finite FAM phenotype at line {line_number}")
                observed += 1
    except OSError as error:
        fail(f"could not validate FAM: {error}")
    if observed != sample_count:
        fail(f"FAM sample count drifted: observed={observed} expected={sample_count}")
    return {"sample_count": observed, "unique_fid_iid_count": len(identities)}


def validate_bim(path: Path, ld: dict[str, Any]) -> dict[str, Any]:
    expected_counts = {
        int(chromosome): int(count)
        for chromosome, count in ld["expected_reference_chromosome_variant_counts"].items()
    }
    chromosome_lengths = {
        int(chromosome): int(length)
        for chromosome, length in ld["grch37_chromosome_lengths"].items()
    }
    counts = {chromosome: 0 for chromosome in ld["reference_chromosomes"]}
    minima: dict[int, int] = {}
    maxima: dict[int, int] = {}
    observed = 0
    previous_chromosome = -1
    previous_position = -1
    allowed_alleles = {b"A", b"C", b"G", b"T"}
    try:
        with path.open("rb") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.split()
                if len(fields) != 6:
                    fail(f"BIM schema drift at line {line_number}")
                chromosome_raw, variant_id, cm_raw, position_raw, allele1, allele2 = fields
                try:
                    chromosome = int(chromosome_raw)
                    genetic_distance = float(cm_raw)
                    position = int(position_raw)
                except ValueError:
                    fail(f"BIM numeric field drift at line {line_number}")
                if (
                    chromosome not in counts or position <= 0
                    or position > chromosome_lengths[chromosome]
                    or not math.isfinite(genetic_distance) or genetic_distance < 0
                    or not variant_id or variant_id == b"."
                    or allele1 not in allowed_alleles or allele2 not in allowed_alleles
                    or allele1 == allele2
                ):
                    fail(f"BIM coordinate/ID/allele drift at line {line_number}")
                if (
                    chromosome < previous_chromosome
                    or (chromosome == previous_chromosome and position < previous_position)
                ):
                    fail(f"BIM coordinate order drift at line {line_number}")
                if chromosome != previous_chromosome:
                    previous_position = -1
                previous_chromosome, previous_position = chromosome, position
                counts[chromosome] += 1
                minima.setdefault(chromosome, position)
                maxima[chromosome] = position
                observed += 1
    except OSError as error:
        fail(f"could not validate BIM: {error}")
    if observed != ld["variant_count"] or counts != expected_counts:
        fail(
            "BIM chromosome/variant counts drifted: "
            f"observed_total={observed} expected_total={ld['variant_count']}"
        )
    autosomal = sum(counts[chromosome] for chromosome in ld["analysis_autosomes"])
    if autosomal != ld["expected_autosomal_variant_count"]:
        fail("BIM eligible autosomal count drifted")
    return {
        "variant_count": observed, "autosomal_variant_count": autosomal,
        "chromosome_variant_counts": {str(key): value for key, value in counts.items()},
        "chromosome_min_position": {str(key): value for key, value in minima.items()},
        "chromosome_max_position": {str(key): value for key, value in maxima.items()},
        "analysis_autosomes": ld["analysis_autosomes"],
        "excluded_source_fidelity_chromosomes": ld["nonanalysis_chromosomes_retained_for_source_fidelity"],
        "allele_schema": "BIALLELIC_ACGT_DISTINCT",
        "coordinate_order": "NONDECREASING_WITHIN_ORDERED_CHROMOSOMES",
        "build_validation": "GRCH37_RELEASE_PIN_PLUS_CHROMOSOME_LENGTH_BOUNDS",
    }


def validate_plink_family(paths: dict[str, Path], ld: dict[str, Any]) -> dict[str, Any]:
    return {
        "bed": validate_bed(paths["bed"], ld["sample_count"], ld["variant_count"]),
        "bim": validate_bim(paths["bim"], ld),
        "fam": validate_fam(paths["fam"], ld["sample_count"]),
    }


def artifact_rows(
    root: Path, paths: dict[str, Path], ld: dict[str, Any],
    *, canonical_paths: dict[str, Path] | None = None,
) -> list[dict[str, Any]]:
    statuses = {
        "bed": "PLINK_BED_SNP_MAJOR_DIMENSIONS_VALID",
        "bim": "BIM_GRCH37_SCHEMA_ALLELES_CHROMOSOME_COUNTS_VALID",
        "fam": f"FAM_SIX_FIELD_{ld['sample_count']}_SAMPLE_FAMILY_VALID",
    }
    return [
        {
            "component_id": component.upper(),
            "path": relative(root, (canonical_paths or paths)[component]),
            "bytes": paths[component].stat().st_size, "sha256": sha256(paths[component]),
            "records": ld["source_members"][component]["records"],
            "validation_status": statuses[component],
        }
        for component in COMPONENT_ORDER
    ]


def manifest_text(rows: list[dict[str, Any]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def provenance_payload(
    root: Path, policy: dict[str, Any], ld: dict[str, Any], policy_path: Path,
    contract_lock_sha256: str, run_fingerprint: str, plink: dict[str, Any],
    paths: dict[str, Path], manifest: str, validation: dict[str, Any],
    resource_metrics: dict[str, Any], *, test_fixture: bool,
) -> dict[str, Any]:
    rows = artifact_rows(root, paths, ld, canonical_paths=reference_paths(root, ld))
    return {
        "schema_version": "sleep-atlas-track-b-pleiotropy-ld-reference.1",
        "analysis_id": policy["analysis_id"], "test_fixture": test_fixture,
        "run_fingerprint": run_fingerprint, "policy": relative(root, policy_path),
        "policy_sha256": sha256(policy_path), "contract_lock_sha256": contract_lock_sha256,
        "materializer": str(SCRIPT),
        "materializer_sha256": sha256(materializer_path(root, test_fixture=test_fixture)),
        "source": {
            "release": ld["source_release"], "archive": ld["source_archive_path"],
            "archive_bytes": ld["source_archive_expected_bytes"],
            "archive_sha256": ld["source_archive_sha256"], "member_pins": ld["source_members"],
        },
        "reference": {
            "directory": ld["materialized_directory"], "prefix": ld["materialized_prefix"],
            "ancestry": ld["ancestry"], "build": ld["build"],
            "sample_count": ld["sample_count"], "variant_count": ld["variant_count"],
            "analysis_autosomes": ld["analysis_autosomes"],
            "eligible_autosomal_variant_count": ld["expected_autosomal_variant_count"],
            "excluded_source_fidelity_chromosomes": ld["nonanalysis_chromosomes_retained_for_source_fidelity"],
            "analysis_chromosome_rule": ld["analysis_chromosome_rule"],
        },
        "plink": plink, "validation": validation, "outputs": rows,
        "materialization_resource_metrics": resource_metrics,
        "materialization_transform_rule": ld["materialization_transform_rule"],
        "manifest": {
            "path": ld["materialized_manifest"], "bytes": len(manifest.encode("utf-8")),
            "sha256": hashlib.sha256(manifest.encode("utf-8")).hexdigest(),
            "rows": len(rows), "schema": MANIFEST_FIELDS,
        },
        "publication": {
            "mode": "ATOMIC_DIRECTORY_RENAME_WITH_KERNEL_NO_REPLACE",
            "commit_marker": ld["materialized_provenance"],
            "partial_or_preexisting_rule": "FAIL_CLOSED_NO_OVERWRITE_NO_AUTOMATIC_RECOVERY",
        },
        "scientific_result": False,
        "claim_limit": "An ancestry/build-matched LD reference is infrastructure, not pleiotropic or mechanistic evidence",
    }


def live_reference_payload(
    root: Path, policy: dict[str, Any], ld: dict[str, Any], policy_path: Path,
    contract_lock_sha256: str, run_fingerprint: str, plink: dict[str, Any],
    resource_metrics: dict[str, Any], *, test_fixture: bool,
) -> tuple[str, dict[str, Any]]:
    paths = reference_paths(root, ld)
    validation = validate_plink_family(paths, ld)
    rows = artifact_rows(root, paths, ld)
    manifest = manifest_text(rows)
    provenance = provenance_payload(
        root, policy, ld, policy_path, contract_lock_sha256, run_fingerprint,
        plink, paths, manifest, validation, resource_metrics, test_fixture=test_fixture,
    )
    return manifest, provenance


def verify_reference(
    root: Path, policy: dict[str, Any], ld: dict[str, Any], policy_path: Path,
    contract_lock_sha256: str, *, test_fixture: bool, quiet: bool = False,
) -> dict[str, Any]:
    if destination_state(root, ld) == "ABSENT":
        fail("BLOCKED_BY_DATA materialized Track B EUR LD reference is absent")
    if destination_state(root, ld) != "COMPLETE_SHAPE":
        fail("BLOCKED_BY_DATA materialized Track B EUR LD reference is partial or tampered")
    archive_identity(root, ld)
    plink = plink_identity(root, ld)
    run_fingerprint = fingerprint(
        root, policy, ld, policy_path, contract_lock_sha256, test_fixture=test_fixture,
    )
    paths = reference_paths(root, ld)
    observed_provenance = read_json(paths["provenance"])
    resource_metrics = observed_provenance.get("materialization_resource_metrics")
    required_metrics = ld["materialization_resource_metrics_required"]
    if (
        not isinstance(resource_metrics, dict) or set(resource_metrics) != set(required_metrics)
        or type(resource_metrics.get("staging_extraction_validation_wall_seconds")) not in {int, float}
        or resource_metrics["staging_extraction_validation_wall_seconds"] < 0
        or type(resource_metrics.get("peak_rss_bytes")) is not int
        or resource_metrics["peak_rss_bytes"] <= 0
        or resource_metrics.get("streaming_chunk_bytes") != 8 * 1024 * 1024
        or resource_metrics.get("extracted_uncompressed_bytes")
        != sum(ld["source_members"][component]["expected_bytes"] for component in COMPONENT_ORDER)
    ):
        fail("materialized Track B LD resource metrics are absent or invalid")
    expected_manifest, expected_provenance = live_reference_payload(
        root, policy, ld, policy_path, contract_lock_sha256, run_fingerprint,
        plink, resource_metrics, test_fixture=test_fixture,
    )
    if paths["manifest"].read_text(encoding="utf-8") != expected_manifest:
        fail("materialized Track B LD manifest differs from live full hashes/counts/schema")
    if observed_provenance != expected_provenance:
        fail("materialized Track B LD provenance differs from the sealed live family")
    if not quiet:
        print(
            "TRACK_B_PLEIOTROPY_LD_VERIFIED "
            f"samples={ld['sample_count']} variants={ld['variant_count']} "
            f"autosomal={ld['expected_autosomal_variant_count']}"
        )
    return observed_provenance


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def publish_directory_no_replace(stage: Path, destination: Path) -> None:
    """Atomically rename one complete same-filesystem stage without replacement."""
    staging_parent = stage.parent
    fsync_directory(stage)
    try:
        library_name = ctypes.util.find_library("c")
        libc = ctypes.CDLL(library_name or None, use_errno=True)
    except OSError as error:
        fail(f"BLOCKED_BY_SOFTWARE could not load libc for atomic publication: {error}")
    source_raw, destination_raw = os.fsencode(stage), os.fsencode(destination)
    if sys.platform == "darwin":
        try:
            rename = libc.renamex_np
        except AttributeError:
            fail("BLOCKED_BY_SOFTWARE renamex_np is unavailable for atomic no-replace publication")
        rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        result = rename(source_raw, destination_raw, 0x00000004)  # RENAME_EXCL
    elif sys.platform.startswith("linux") and hasattr(libc, "renameat2"):
        rename = libc.renameat2
        rename.argtypes = [
            ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint,
        ]
        rename.restype = ctypes.c_int
        result = rename(-100, source_raw, -100, destination_raw, 0x00000001)  # RENAME_NOREPLACE
    else:
        fail("BLOCKED_BY_SOFTWARE kernel atomic no-replace directory publication is unavailable")
    if result != 0:
        error_number = ctypes.get_errno()
        if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
            fail("Track B LD destination appeared during staging; atomic no-replace publication refused")
        fail(
            "atomic no-replace Track B LD publication failed; fingerprinted stage retained: "
            + os.strerror(error_number)
        )
    fsync_directory(staging_parent)
    fsync_directory(destination.parent)


def materialize_reference(
    root: Path, policy: dict[str, Any], ld: dict[str, Any], policy_path: Path,
    contract_lock_sha256: str, *, test_fixture: bool,
) -> dict[str, Any]:
    state = destination_state(root, ld)
    if state == "COMPLETE_SHAPE":
        return verify_reference(
            root, policy, ld, policy_path, contract_lock_sha256,
            test_fixture=test_fixture,
        )
    checks = preflight(
        root, policy, ld, policy_path, contract_lock_sha256, test_fixture=test_fixture,
    )
    resource_started = time.monotonic()
    stage = safe_path(root, checks["staging_directory"])
    stage.parent.mkdir(parents=True, exist_ok=True)
    try:
        stage.mkdir()
    except FileExistsError:
        fail("fingerprinted Track B LD staging directory already exists; automatic reuse is forbidden")
    archive_path = safe_path(root, checks["archive_path"])
    stage_paths = {
        component: stage / ld["source_members"][component]["name"]
        for component in COMPONENT_ORDER
    }
    stage_paths["manifest"] = stage / Path(ld["materialized_manifest"]).name
    stage_paths["provenance"] = stage / Path(ld["materialized_provenance"]).name
    try:
        with zipfile.ZipFile(archive_path) as archive:
            for component in COMPONENT_ORDER:
                member = ld["source_members"][component]
                info = checks["archive_infos"][member["name"]]
                copy_zip_member(
                    archive, info, stage_paths[component],
                    member["expected_bytes"], member["sha256"],
                )
        validation = validate_plink_family(stage_paths, ld)
        canonical_paths = reference_paths(root, ld)
        rows = artifact_rows(root, stage_paths, ld, canonical_paths=canonical_paths)
        manifest = manifest_text(rows)
        atomic_text(stage_paths["manifest"], manifest)
        peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != "darwin":
            peak_rss *= 1024
        resource_metrics = {
            "staging_extraction_validation_wall_seconds": round(time.monotonic() - resource_started, 6),
            "peak_rss_bytes": int(peak_rss),
            "streaming_chunk_bytes": 8 * 1024 * 1024,
            "extracted_uncompressed_bytes": sum(
                ld["source_members"][component]["expected_bytes"] for component in COMPONENT_ORDER
            ),
        }
        provenance = provenance_payload(
            root, policy, ld, policy_path, contract_lock_sha256, checks["fingerprint"],
            checks["plink"], stage_paths, manifest, validation, resource_metrics,
            test_fixture=test_fixture,
        )
        atomic_text(stage_paths["provenance"], json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    except BaseException:
        # Preserve an incomplete fingerprinted stage for audit; never auto-promote or reuse it.
        raise

    destination = safe_path(root, ld["materialized_directory"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    publish_directory_no_replace(stage, destination)
    verified = verify_reference(
        root, policy, ld, policy_path, contract_lock_sha256,
        test_fixture=test_fixture, quiet=True,
    )
    print(
        "TRACK_B_PLEIOTROPY_LD_MATERIALIZED "
        f"fingerprint={checks['fingerprint']} samples={ld['sample_count']} "
        f"variants={ld['variant_count']}"
    )
    return verified


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default=str(PRODUCTION_POLICY))
    parser.add_argument(
        "--test-fixture", action="store_true",
        help="accept only the explicit test-only fixture schema; never valid for production",
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--materialize", action="store_true")
    action.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_relative = Path(args.policy)
    policy, ld, policy_path, contract_lock_sha256 = load_policy_context(
        root, policy_relative, test_fixture=args.test_fixture,
    )
    if args.preflight:
        checks = preflight(
            root, policy, ld, policy_path, contract_lock_sha256,
            test_fixture=args.test_fixture,
        )
        if checks["state"] == "COMPLETE_SHAPE":
            verify_reference(
                root, policy, ld, policy_path, contract_lock_sha256,
                test_fixture=args.test_fixture,
            )
        else:
            print(
                "TRACK_B_PLEIOTROPY_LD_PREFLIGHT_READY "
                f"fingerprint={checks['fingerprint']} free_bytes={checks['free_bytes']} "
                f"required_free_bytes={checks['required_free_bytes']}"
            )
    elif args.materialize:
        materialize_reference(
            root, policy, ld, policy_path, contract_lock_sha256,
            test_fixture=args.test_fixture,
        )
    else:
        verify_reference(
            root, policy, ld, policy_path, contract_lock_sha256,
            test_fixture=args.test_fixture,
        )


if __name__ == "__main__":
    main()
