#!/usr/bin/env python3
"""Freeze the Track B fine-mapping family before any locus materialization.

The sole scientific input to this gate is the stable, side-effect-free
``verify_results(root=...)`` API in script 147. This module never opens a raw
PLACO ledger and never treats the terminal failed-QC LAVA family as support or
ranking evidence. Verified independent PLACO leads are grouped by pair and
official LAVA block; an official block is only a work partition, not evidence
that its co-members are one signal or share a causal variant.

The public :func:`verify_results` API is side-effect free and intentionally
stable for scripts 148--151. Publication is explicit, no-replace, resumable,
and commits ``DEEP_LOCUS_MANIFEST.lock.json`` last.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import stat
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, NamedTuple, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_REL = Path("scripts/145_build_track_b_finemapping_continuation_gate.py")
ARCHIVE_MANAGER_REL = Path("scripts/146_manage_lava_reference_archives.py")
PLEIOTROPY_VERIFIER_REL = Path("scripts/147_build_track_b_pleiotropy_results_v2.py")
PLEIOTROPY_RESULT_PROVENANCE_REL = Path(
    "results/track_b/pleiotropy/results/results.provenance.json"
)
POLICY_REL = Path("config/track_b_finemapping_policy.json")
OFFICIAL_BLOCKS_REL = Path("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
BASE_CONTRACT_SCRIPT_REL = Path("scripts/130_track_b_finemapping_contract.py")
BASE_CONTRACT_LOCK_REL = Path("results/track_b/finemapping/contract.lock.json")
BASE_READINESS_REL = Path("results/track_b/finemapping/readiness.tsv")
PLEIOTROPY_CONTRACT_LOCK_REL = Path("results/track_b/pleiotropy/contract.lock.json")
PLEIOTROPY_INPUT_GATE_LOCK_REL = Path("results/track_b/pleiotropy/input_gate.lock.json")
LOCAL_POLICY_REL = Path("config/track_b_local_analysis_policy.json")
LOCAL_INPUT_LOCK_REL = Path("results/track_b/local_analysis_input.lock.json")

FAMILY_REL = Path("results/track_b/DEEP_LOCUS_MANIFEST.tsv")
FAMILY_LOCK_REL = Path("results/track_b/DEEP_LOCUS_MANIFEST.lock.json")
UNAVAILABLE_REL = Path("results/track_b/DEEP_LOCUS_UNAVAILABLE.tsv")
ZERO_PROVENANCE_REL = Path("results/track_b/DEEP_LOCUS_ZERO_FAMILY.provenance.json")

ANALYSIS_ID = "track-b-v1.0-finemapping-trait-coloc"
UPSTREAM_ANALYSIS_ID = "track-b-v1.0-pleiotropy"
SCHEMA = "sleep-atlas-track-b-finemapping-deep-locus-family.1"
ARTIFACT_ROLE = "PRE_MATERIALIZATION_PAIR_BY_OFFICIAL_BLOCK_FAMILY"
ZERO_SCHEMA = "sleep-atlas-track-b-finemapping-deep-locus-zero-family.1"
UPSTREAM_VERIFIER_SCHEMA = "sleep-atlas-track-b-post-placo-verified-family.2"
UPSTREAM_VERIFIER_SHA256 = "cc0cbb0b9b3f8c6045b02367a4d4d4ed2a8b21804f376318efce429ea26518a5"
OFFICIAL_BLOCK_SHA256 = "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882"

UPSTREAM_ZERO = "ZERO_ELIGIBLE_SIGNAL_FAMILY"
UPSTREAM_BLOCKED = "BLOCKED_BY_LD_REFERENCE_COVERAGE"
UPSTREAM_COMPLETE = "COMPLETE_WITH_INDEPENDENT_LOCI"
UPSTREAM_STATES = {UPSTREAM_ZERO, UPSTREAM_BLOCKED, UPSTREAM_COMPLETE}
READY_STATE = "READY_WHOLE_LOCUS_MATERIALIZATION"
ZERO_STATE = "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY"
CONJFDR_BLOCKED = "NOT_TESTED_METHOD_BLOCKED"

# Literal constants and the executable verify_reference_state call are required
# by script 146's downstream-consumer audit.
ARCHIVES_PRESENT = "ARCHIVES_PRESENT_FULLY_VERIFIED"
ARCHIVES_EVICTED = (
    "ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_"
    "EXTRACTED_PAYLOADS_REHASHED"
)
ACCEPTED_ARCHIVE_STATES = {ARCHIVES_PRESENT, ARCHIVES_EVICTED}

PAIR_ORDER = ("A", "B", "CONTROL")
PAIR_RANK = {pair: index for index, pair in enumerate(PAIR_ORDER)}
PAIR_IDENTITIES = {
    "A": ("snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    "B": ("insomnia", "adhd", "PRIMARY_DISCOVERY"),
    "CONTROL": ("insomnia", "frailty", "POSITIVE_CONTROL_NON_NOVELTY"),
}
UPSTREAM_PAIR_ROLES = {
    "A": "PRIMARY_DISCOVERY",
    "B": "PRIMARY_DISCOVERY",
    "CONTROL": "POSITIVE_CONTROL",
}

FAMILY_FIELDS = [
    "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
    "CHR", "START", "STOP", "ld_block_id", "inclusion_sources", "priority_tier",
    "local_evidence_ids", "local_lead_variants", "pleiotropy_lead_evidence_ids",
    "pleiotropy_lead_variants", "pleiotropy_clump_evidence_ids",
    "method_availability", "partition_role", "selection_rule", "claim_limit",
]
UNAVAILABLE_FIELDS = [
    "analysis_id", "unavailable_entry_id", "pair_id", "family_role", "trait1",
    "trait2", "CHR", "BP", "lead_variant", "ld_block_id", "evidence_ids",
    "inclusion_sources", "terminal_status", "blocker_reason", "claim_limit",
]
UPSTREAM_LEAD_FIELDS = [
    "pair_id", "lead_variant", "chr", "position", "PLACO_P", "FDR",
    "trait1_P", "trait2_P", "locus_start", "locus_end",
    "independent_signal", "annotations", "analysis_id", "family_role",
    "trait1", "trait2", "lead_A1", "lead_A2", "lead_Z1", "lead_Z2",
    "evidence_labels", "primary_headline", "eligibility_basis", "clump_id",
    "clump_evidence_ids", "lava_block_id", "lava_block_start", "lava_block_end",
    "pair_block_evidence_ids", "ld_reference_status", "plink_clump_status",
    "pairwise_r2_status", "max_pairwise_r2_within_1mb", "lava_terminal_state",
    "lava_failed_qc_accounting", "conditional_status", "method_availability",
    "claim_limit",
]
UPSTREAM_TOP_LEVEL_KEYS = {
    "schema_version", "state", "fingerprint", "zero_family", "primary_leads",
    "control_leads", "complete_family", "archive_state", "archive_reference",
    "ld_reference_coverage", "method_availability",
}
LOCK_KEYS = {
    "schema_version", "analysis_id", "artifact_role", "state",
    "results_accessed_before_lock", "pair_order", "manifest", "locus_count",
    "locus_entry_ids_in_order", "evidence_count", "evidence_ids_sha256", "policy",
    "official_blocks", "upstream_family_provenance", "generator_scripts",
    "selection_contract", "zero_family_provenance", "unavailable_manifest",
    "unavailable_count", "unavailable_entry_ids_in_order",
    "unavailable_evidence_count", "unavailable_evidence_ids_sha256",
    "upstream_evidence_count", "upstream_evidence_ids_sha256",
}
SELECTION_CONTRACT = {
    "pair_by_official_block_grouping": True,
    "all_eligible_evidence_retained": True,
    "cross_pair_collapse": False,
    "consensus_only": False,
    "result_ranked_filtering": False,
    "maximum_locus_variants": None,
    "locus_splitting": False,
    "variant_thinning": False,
    "lead_centered_truncation": False,
}
PARTITION_ROLE = "OFFICIAL_LAVA_BLOCK_COORDINATES_AS_WORK_PARTITION_ONLY"
SELECTION_RULE = (
    "ALL_MAPPED_PAIR_BY_OFFICIAL_BLOCK_ENTRIES_CONTINUE;"
    "PRIORITY_IS_PRE_FINEMAP_NON_FILTERING_LABEL_ONLY"
)
METHOD_AVAILABILITY = f"PLACO:COMPLETE_WITH_HITS;CONJFDR:{CONJFDR_BLOCKED}"
PRIMARY_CLAIM = (
    "STATISTICAL_PLEIOTROPY_ONLY_NO_SHARED_CAUSAL_VARIANT_NO_MECHANISM_OR_MEDIATION"
)
CONTROL_CLAIM = "CONTROL_RECOVERY_ONLY_NO_NOVELTY_NO_MECHANISM"

DOWNSTREAM_OR_COMPETING_PATHS = {
    Path("results/track_b/finemapping/locus_manifest.tsv"),
    Path("results/track_b/finemapping/locus_manifest.lock.json"),
    Path("results/track_b/finemapping/zero_family.provenance.json"),
    Path("results/track_b/finemapping/results.provenance.json"),
    Path("results/track_b/finemapping/diagnostics.tsv"),
    Path("results/track_b/finemapping/resource_metrics.tsv"),
    Path("results/track_b/finemapping/run_index.tsv"),
    Path("results/track_b/finemapping/RAM_BY_LOCUS.namespace.json"),
    Path("results/track_b/finemapping/materialization_metrics"),
    Path("results/track_b/finemapping/materialization_attempts"),
    Path("results/track_b/finemapping/runs"),
    Path("results/track_b/finemapping/publications"),
    Path("results/track_b/10_finemap_trait1.tsv"),
    Path("results/track_b/11_finemap_trait2.tsv"),
    Path("results/track_b/12_trait_trait_coloc.tsv"),
    Path("work/track_b_finemapping/inputs"),
    Path("work/track_b_finemapping/runs"),
    Path("work/track_b_finemapping/quarantine"),
    Path("work/track_b_finemapping/leases"),
    Path("work/track_b_finemapping/publications"),
}
OWN_PATHS = {FAMILY_REL, FAMILY_LOCK_REL, UNAVAILABLE_REL, ZERO_PROVENANCE_REL}

HEX64 = re.compile(r"^[0-9a-f]{64}$")
RSID = re.compile(r"^rs[0-9]+$", re.IGNORECASE)
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
BASES = frozenset("ACGT")
MAX_STRUCTURED_BYTES = 128 * 1024**2


class GateError(RuntimeError):
    """A fail-closed input, scientific-contract, or publication violation."""


class UpstreamBlocked(GateError):
    """The verified upstream family cannot produce a runnable frozen family."""


class Block(NamedTuple):
    locus: int
    chromosome: int
    start: int
    stop: int


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _stat_tuple(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def safe_path(
    root: Path, relative: Path | str, label: str, *, must_exist: bool = False,
) -> Path:
    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise GateError(f"unsafe {label} path: {relative}")
    try:
        root_real = root.resolve(strict=True)
    except OSError as error:
        raise GateError(f"repository root is unavailable: {error}") from error
    current = root
    for offset, part in enumerate(relative.parts):
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise GateError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise GateError(f"{label} contains a symbolic link: {current}")
        if offset < len(relative.parts) - 1 and not stat.S_ISDIR(observed.st_mode):
            raise GateError(f"{label} has a non-directory ancestor: {current}")
    try:
        path = root / relative
        path.resolve(strict=False).relative_to(root_real)
    except (OSError, RuntimeError, ValueError) as error:
        raise GateError(f"{label} escapes repository: {relative}") from error
    if must_exist and not path.exists():
        raise GateError(f"required {label} is absent: {relative}")
    return path


def stable_bytes(
    root: Path, relative: Path | str, *, allow_empty: bool = False,
) -> tuple[bytes, dict[str, object]]:
    path = safe_path(root, relative, "artifact", must_exist=True)
    relative_string = str(Path(relative))
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        named_before = os.lstat(path)
        if (
            not stat.S_ISREG(named_before.st_mode)
            or (named_before.st_size == 0 and not allow_empty)
            or named_before.st_size > MAX_STRUCTURED_BYTES
        ):
            raise GateError(f"artifact is not a permitted regular file: {relative_string}")
        descriptor = os.open(path, flags)
        try:
            opened = os.fstat(descriptor)
            chunks: list[bytes] = []
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                chunks.append(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        named_after = os.lstat(path)
    except OSError as error:
        raise GateError(f"could not read stable artifact {relative_string}: {error}") from error
    if len({_stat_tuple(item) for item in (named_before, opened, after, named_after)}) != 1:
        raise GateError(f"artifact changed while reading: {relative_string}")
    content = b"".join(chunks)
    return content, {
        "path": relative_string,
        "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def stable_identity(
    root: Path, relative: Path | str, *, allow_empty: bool = False,
) -> dict[str, object]:
    """Hash a file without retaining it in memory, detecting named-inode change."""

    path = safe_path(root, relative, "artifact", must_exist=True)
    relative_string = str(Path(relative))
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    digest = hashlib.sha256()
    try:
        named_before = os.lstat(path)
        if not stat.S_ISREG(named_before.st_mode) or (
            named_before.st_size == 0 and not allow_empty
        ):
            raise GateError(f"artifact is not a permitted regular file: {relative_string}")
        descriptor = os.open(path, flags)
        try:
            opened = os.fstat(descriptor)
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        named_after = os.lstat(path)
    except OSError as error:
        raise GateError(f"could not hash stable artifact {relative_string}: {error}") from error
    if len({_stat_tuple(item) for item in (named_before, opened, after, named_after)}) != 1:
        raise GateError(f"artifact changed while hashing: {relative_string}")
    return {
        "path": relative_string,
        "bytes": int(named_after.st_size),
        "sha256": digest.hexdigest(),
    }


def identity_for_content(relative: Path | str, content: bytes) -> dict[str, object]:
    return {
        "path": str(Path(relative)),
        "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def identity_matches(observed: Mapping[str, Any], expected: Mapping[str, Any]) -> bool:
    return all(observed.get(key) == expected.get(key) for key in ("path", "bytes", "sha256"))


def _identity_record(value: Any, label: str) -> dict[str, object]:
    if (
        not isinstance(value, Mapping)
        or set(value) != {"path", "bytes", "sha256"}
        or not isinstance(value.get("path"), str)
        or Path(str(value["path"])).is_absolute()
        or ".." in Path(str(value["path"])).parts
        or not isinstance(value.get("bytes"), int)
        or isinstance(value.get("bytes"), bool)
        or int(value["bytes"]) <= 0
        or not isinstance(value.get("sha256"), str)
        or not HEX64.fullmatch(str(value["sha256"]))
    ):
        raise GateError(f"malformed immutable identity for {label}")
    return {key: value[key] for key in ("path", "bytes", "sha256")}


def read_json(root: Path, relative: Path | str) -> tuple[dict[str, Any], dict[str, object]]:
    content, identity = stable_bytes(root, relative)
    try:
        value = json.loads(content.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise GateError(f"invalid JSON artifact {relative}: {error}") from error
    if not isinstance(value, dict):
        raise GateError(f"JSON artifact must contain one object: {relative}")
    return value, identity


def read_tsv(
    root: Path, relative: Path | str, fields: Sequence[str], *, allow_header_only: bool,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    content, identity = stable_bytes(root, relative)
    try:
        reader = csv.DictReader(io.StringIO(content.decode("utf-8")), delimiter="\t")
        observed_fields = list(reader.fieldnames or [])
        rows = list(reader)
    except (UnicodeError, csv.Error) as error:
        raise GateError(f"invalid TSV artifact {relative}: {error}") from error
    if (
        observed_fields != list(fields)
        or any(None in row or any(value is None for value in row.values()) for row in rows)
        or (not rows and not allow_header_only)
    ):
        raise GateError(f"TSV schema or row family drifted: {relative}")
    return rows, identity


def tsv_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return buffer.getvalue().encode("utf-8")


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise GateError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        specification.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


def _integer(value: Any, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise GateError(f"invalid integer {label}: {value!r}")
    try:
        parsed = int(str(value))
    except (TypeError, ValueError, OverflowError) as error:
        raise GateError(f"invalid integer {label}: {value!r}") from error
    if str(parsed) != str(value) or parsed < minimum:
        raise GateError(f"invalid integer {label}: {value!r}")
    return parsed


def _probability(value: Any, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise GateError(f"invalid probability {label}: {value!r}") from error
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        raise GateError(f"probability outside [0,1] for {label}: {value!r}")
    return parsed


def _finite(value: Any, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise GateError(f"invalid number {label}: {value!r}") from error
    if not math.isfinite(parsed):
        raise GateError(f"non-finite number {label}: {value!r}")
    return parsed


def _tokens(value: Any, label: str, *, allow_na: bool = False) -> list[str]:
    if value == "NA" and allow_na:
        return []
    if not isinstance(value, str) or value == "NA":
        raise GateError(f"{label} is absent")
    tokens = value.split(";")
    if any(not token or token.strip() != token for token in tokens):
        raise GateError(f"{label} contains an empty or whitespace-drifted token")
    if len(tokens) != len(set(tokens)):
        raise GateError(f"{label} contains duplicate identities")
    return tokens


def _parse_evidence_id(value: str, pair_id: str) -> tuple[str, int, int]:
    fields = value.split(":")
    if len(fields) != 5 or fields[0] != "PLACO" or fields[1] != pair_id:
        raise GateError("PLACO evidence ID does not encode its exact method and pair")
    snp = fields[2]
    if not RSID.fullmatch(snp):
        raise GateError("PLACO evidence ID does not contain an rsID")
    chromosome = _integer(fields[3], "evidence chromosome", minimum=1)
    position = _integer(fields[4], "evidence position", minimum=1)
    if chromosome > 22:
        raise GateError("PLACO evidence chromosome is not autosomal")
    return snp, chromosome, position


def parse_blocks_text(
    text: str, *, expected_count: int | None = None,
) -> tuple[list[Block], dict[int, tuple[list[int], list[Block]]]]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or lines[0].split() != ["LOC", "CHR", "START", "STOP"]:
        raise GateError("official LAVA block schema drifted")
    blocks: list[Block] = []
    for row_number, line in enumerate(lines[1:], start=2):
        values = line.split()
        if len(values) != 4:
            raise GateError(f"malformed official LAVA block row {row_number}")
        try:
            locus, chromosome, start, stop = map(int, values)
        except ValueError as error:
            raise GateError(f"nonnumeric official LAVA block row {row_number}") from error
        if locus < 1 or chromosome not in range(1, 23) or start < 1 or stop < start:
            raise GateError(f"invalid official LAVA block row {row_number}")
        blocks.append(Block(locus, chromosome, start, stop))
    if expected_count is not None and len(blocks) != expected_count:
        raise GateError(f"official LAVA block count is {len(blocks)}, expected {expected_count}")
    if [block.locus for block in blocks] != list(range(1, len(blocks) + 1)):
        raise GateError("official LAVA block IDs are not the exact ordered 1..N family")
    if blocks != sorted(blocks, key=lambda item: (item.chromosome, item.start, item.stop)):
        raise GateError("official LAVA block family is not in genomic order")
    by_chromosome: dict[int, list[Block]] = defaultdict(list)
    for block in blocks:
        by_chromosome[block.chromosome].append(block)
    if expected_count == 2495 and set(by_chromosome) != set(range(1, 23)):
        raise GateError("official LAVA block family does not cover all autosomes")
    index: dict[int, tuple[list[int], list[Block]]] = {}
    for chromosome, family in by_chromosome.items():
        previous_stop = 0
        for block in family:
            if block.start <= previous_stop:
                raise GateError(f"official LAVA blocks overlap on chromosome {chromosome}")
            previous_stop = block.stop
        index[chromosome] = ([block.start for block in family], family)
    return blocks, index


def map_to_block(
    chromosome: int, position: int,
    index: Mapping[int, tuple[list[int], list[Block]]],
) -> Block | None:
    family = index.get(chromosome)
    if family is None:
        return None
    starts, blocks = family
    offset = bisect.bisect_right(starts, position) - 1
    if offset < 0 or position > blocks[offset].stop:
        return None
    return blocks[offset]


def load_official_blocks(
    root: Path = ROOT,
) -> tuple[list[Block], dict[int, tuple[list[int], list[Block]]], dict[str, object]]:
    content, identity = stable_bytes(root, OFFICIAL_BLOCKS_REL)
    if identity["sha256"] != OFFICIAL_BLOCK_SHA256:
        raise GateError("official LAVA block identity differs from the frozen reference")
    try:
        text = content.decode("utf-8")
    except UnicodeError as error:
        raise GateError("official LAVA block file is not UTF-8") from error
    blocks, index = parse_blocks_text(text, expected_count=2495)
    return blocks, index, identity


def validate_policy(root: Path = ROOT) -> tuple[dict[str, Any], dict[str, object]]:
    policy, identity = read_json(root, POLICY_REL)
    pairs = policy.get("pair_scope", {})
    pair_records = pairs.get("pairs") if isinstance(pairs, Mapping) else None
    expected_pairs = [
        {
            "pair_id": pair,
            "trait1": PAIR_IDENTITIES[pair][0],
            "trait2": PAIR_IDENTITIES[pair][1],
            "family_role": PAIR_IDENTITIES[pair][2],
        }
        for pair in PAIR_ORDER
    ]
    locus = policy.get("locus_entry_contract", {})
    dense = policy.get("dense_input_contract", {})
    signed = policy.get("signed_ld_contract", {})
    execution = policy.get("ram_aware_execution_contract", {})
    outputs = policy.get("required_science_outputs", {})
    if (
        policy.get("schema_version") != "sleep-atlas-track-b-finemapping-policy.1"
        or policy.get("analysis_id") != ANALYSIS_ID
        or policy.get("analysis_build") != "hg19"
        or policy.get("ancestry") != "EUR"
        or policy.get("results_accessed_before_contract_freeze") is not False
        or not isinstance(pairs, Mapping)
        or pairs.get("pair_order") != list(PAIR_ORDER)
        or pair_records != expected_pairs
        or not isinstance(locus, Mapping)
        or locus.get("consensus_only_selection_forbidden") is not True
        or locus.get("result_ranked_maximum_forbidden") is not True
        or locus.get("intersection_required") is not False
        or not isinstance(dense, Mapping)
        or dense.get("maximum_locus_variants") is not None
        or "no count cap" not in str(dense.get("locus_size_rule", ""))
        or not isinstance(signed, Mapping)
        or signed.get("reference_id") != "LAVA_UKB_v1.1_EUR_GRCh37"
        or signed.get("ancestry") != "EUR"
        or signed.get("build") != "GRCh37/hg19"
        or not isinstance(execution, Mapping)
        or execution.get("execution_unit") != "ONE_VALIDATED_PAIR_LOCUS_PER_FRESH_OS_PROCESS"
        or execution.get("full_locus_variant_universe_required") is not True
        or execution.get("full_ancestry_matched_signed_ld_required") is not True
        or execution.get("locus_splitting_forbidden") is not True
        or execution.get("variant_thinning_for_compute_forbidden") is not True
        or not isinstance(outputs, Mapping)
        or [outputs.get(key, {}).get("path") for key in ("output_10", "output_11", "output_12")]
        != [
            "results/track_b/10_finemap_trait1.tsv",
            "results/track_b/11_finemap_trait2.tsv",
            "results/track_b/12_trait_trait_coloc.tsv",
        ]
    ):
        raise GateError("fine-mapping policy lost its pair, whole-locus, or signed-LD contract")
    return policy, identity


def validate_base_contract_snapshot(root: Path = ROOT) -> dict[str, Any]:
    """Read-only compatibility audit for the immutable historical contract.

    The historical lock is not rewritten when additive policy inputs change.
    The corrected freezer binds the current policy and scripts in its new lock.
    """

    policy, policy_identity = validate_policy(root)
    lock, lock_identity = read_json(root, BASE_CONTRACT_LOCK_REL)
    contract_identity = stable_identity(root, BASE_CONTRACT_SCRIPT_REL)
    readiness_identity = stable_identity(root, BASE_READINESS_REL)
    if (
        lock.get("schema_version") != "sleep-atlas-track-b-finemapping-contract-lock.1"
        or lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("artifact_role")
        != "PRE_RESULT_FINE_MAPPING_AND_TRAIT_COLOC_CONTRACT_NOT_SCIENTIFIC_RESULT"
        or lock.get("results_accessed_before_contract_freeze") is not False
        or lock.get("science_outputs_created_by_contract") is not False
        or lock.get("readiness_sha256") != readiness_identity["sha256"]
        or lock.get("zero_family_terminal_status") != ZERO_STATE
        or lock.get("zero_family_creates_10_11_12") is not False
        or lock.get("consensus_only_selection") != "FORBIDDEN"
        or lock.get("result_ranked_maximum") != "FORBIDDEN"
        or lock.get("locus_splitting") != "FORBIDDEN"
    ):
        raise GateError("historical fine-mapping contract identity or invariants drifted")
    historical_inputs = lock.get("immutable_input_sha256")
    if not isinstance(historical_inputs, Mapping):
        raise GateError("historical fine-mapping lock lacks immutable input identities")
    permitted_additive_lineage = {
        str(POLICY_REL), str(PLEIOTROPY_CONTRACT_LOCK_REL),
        str(PLEIOTROPY_INPUT_GATE_LOCK_REL),
    }
    dependency_lineage: dict[str, dict[str, Any]] = {}
    for relative, expected_hash in historical_inputs.items():
        if not isinstance(relative, str) or not isinstance(expected_hash, str):
            raise GateError("historical fine-mapping dependency identity is malformed")
        current_hash = stable_identity(root, relative)["sha256"]
        matches = current_hash == expected_hash
        if not matches and relative not in permitted_additive_lineage:
            raise GateError(f"historically locked fine-mapping dependency drifted: {relative}")
        dependency_lineage[relative] = {
            "historical_sha256": expected_hash,
            "current_sha256": current_hash,
            "matches_historical_lock": matches,
            "handling": (
                "UNCHANGED" if matches
                else "PERMITTED_ADDITIVE_PRE_RESULT_LINEAGE_BOUND_BY_THIS_CONTINUATION"
            ),
        }
    current_contract = _load_module(
        "track_b_finemapping_historical_contract_for_145",
        safe_path(root, BASE_CONTRACT_SCRIPT_REL, "base contract", must_exist=True),
    )
    try:
        validated_policy = current_contract.load_policy()
    except SystemExit as error:
        raise GateError(f"current fine-mapping policy failed semantic validation: {error}") from error
    if validated_policy != policy:
        raise GateError("base contract validator returned a different current policy")
    return {
        "policy": policy,
        "contract_lock": lock,
        "lineage": {
            "historical_lock_preserved": True,
            "historical_policy_sha256": lock.get("policy_sha256"),
            "current_policy_sha256": policy_identity["sha256"],
            "policy_identity_matches_historical_lock": (
                lock.get("policy_sha256") == policy_identity["sha256"]
            ),
            "historical_contract_script_sha256": lock.get("script_sha256"),
            "current_contract_script_sha256": contract_identity["sha256"],
            "contract_script_identity_matches_historical_lock": (
                lock.get("script_sha256") == contract_identity["sha256"]
            ),
            "historical_dependency_lineage": dependency_lineage,
            "mismatch_handling": (
                "DISCLOSED_AND_BOUND_IN_ADDITIVE_PRE_RESULT_CONTINUATION;"
                "HISTORICAL_LOCK_NOT_REWRITTEN"
            ),
        },
        "identities": {
            "policy": policy_identity,
            "contract_script": contract_identity,
            "contract_lock": lock_identity,
            "readiness": readiness_identity,
            "pleiotropy_contract_lock": stable_identity(root, PLEIOTROPY_CONTRACT_LOCK_REL),
            "pleiotropy_input_gate_lock": stable_identity(root, PLEIOTROPY_INPUT_GATE_LOCK_REL),
            "local_analysis_policy": stable_identity(root, LOCAL_POLICY_REL),
            "local_analysis_input_lock": stable_identity(root, LOCAL_INPUT_LOCK_REL),
        },
    }


def _expected_labels(pair_id: str, p_value: float, q_value: float) -> list[str]:
    labels: list[str] = []
    if pair_id == "CONTROL":
        if p_value <= 5e-8:
            labels.append("CONTROL_PLACO_GWS_RECOVERED")
        if q_value <= 0.05:
            labels.append("CONTROL_PLACO_WITHIN_PAIR_FDR")
    else:
        if p_value <= 2.5e-8:
            labels.append("PLACO_PRIMARY_HEADLINE")
        if p_value <= 5e-8:
            labels.append("PLACO_PAIRWISE_GWS")
        if q_value <= 0.05:
            labels.append("PLACO_WITHIN_PAIR_FDR")
    return labels


def _expected_eligibility(p_value: float, q_value: float) -> list[str]:
    values: list[str] = []
    if p_value <= 5e-8:
        values.append("PLACO_P_LE_5E_8")
    if q_value <= 0.05:
        values.append("PLACO_COMPLETE_WITHIN_PAIR_BH_Q_LE_0_05")
    return values


def _validate_complete_family(value: Any) -> dict[str, str]:
    if not isinstance(value, Mapping) or set(value) != {"order", "rows", "canonical_07", "pairs"}:
        raise GateError("147 complete-family payload has an unexpected shape")
    total = _integer(value.get("rows"), "complete-family row count", minimum=0)
    if value.get("order") != list(PAIR_ORDER):
        raise GateError("147 complete-family order is not exact A/B/CONTROL")
    canonical = _identity_record(value.get("canonical_07"), "147 canonical 07")
    if canonical["path"] != "results/track_b/07_placo_plus_variants.tsv":
        raise GateError("147 complete-family canonical 07 path drifted")
    pairs = value.get("pairs")
    if not isinstance(pairs, Mapping) or set(pairs) != set(PAIR_ORDER):
        raise GateError("147 complete-family pair map is not exact A/B/CONTROL")
    observed_total = 0
    pair_states: dict[str, str] = {}
    count_fields = {
        "rows", "failures", "primary", "pairwise", "bh", "minimum_p", "terminal_status",
    }
    for pair_id in PAIR_ORDER:
        pair = pairs[pair_id]
        if not isinstance(pair, Mapping) or set(pair) != {
            "terminal_state", "rows", "ledger", "provenance", "counts",
        }:
            raise GateError(f"147 complete-family record shape drifted for {pair_id}")
        rows = _integer(pair.get("rows"), f"{pair_id} complete rows", minimum=1)
        observed_total += rows
        _identity_record(pair.get("ledger"), f"{pair_id} ledger")
        _identity_record(pair.get("provenance"), f"{pair_id} provenance")
        counts = pair.get("counts")
        if not isinstance(counts, Mapping) or set(counts) != count_fields:
            raise GateError(f"147 count family shape drifted for {pair_id}")
        if _integer(counts.get("rows"), f"{pair_id} count rows", minimum=1) != rows:
            raise GateError(f"147 pair row denominator drifted for {pair_id}")
        parsed_counts: dict[str, int] = {}
        for key in ("failures", "primary", "pairwise", "bh"):
            count = _integer(counts.get(key), f"{pair_id} {key}", minimum=0)
            if count > rows:
                raise GateError(f"147 {pair_id} {key} exceeds its complete row family")
            parsed_counts[key] = count
        _probability(counts.get("minimum_p"), f"{pair_id} minimum P")
        terminal = counts.get("terminal_status")
        if terminal not in {"COMPLETE_WITH_HITS", "TESTED_NO_HIT"}:
            raise GateError(f"147 pair terminal state is invalid for {pair_id}")
        if pair.get("terminal_state") != terminal:
            raise GateError(f"147 pair terminal state/count accounting differs for {pair_id}")
        expected_terminal = (
            "COMPLETE_WITH_HITS"
            if any(parsed_counts[key] for key in ("primary", "pairwise", "bh"))
            else "TESTED_NO_HIT"
        )
        if terminal != expected_terminal:
            raise GateError(f"147 pair hit counts do not imply its terminal state for {pair_id}")
        pair_states[pair_id] = str(terminal)
    if observed_total != total:
        raise GateError("147 complete-family total does not equal its pair denominators")
    return pair_states


def _validate_archive_payload(payload: Mapping[str, Any]) -> None:
    state = payload.get("archive_state")
    reference = payload.get("archive_reference")
    if (
        state not in ACCEPTED_ARCHIVE_STATES
        or not isinstance(reference, Mapping)
        or reference.get("state") != state
        or reference.get("extracted_payload_count") != 44
        or not HEX64.fullmatch(str(reference.get("archive_family_sha256", "")))
        or not HEX64.fullmatch(str(reference.get("extracted_family_sha256", "")))
    ):
        raise GateError("147 did not return one fully verified script-146 archive state")


def _validate_lead(
    raw: Any, expected_collection: str,
    block_index: Mapping[int, tuple[list[int], list[Block]]],
) -> dict[str, Any]:
    if not isinstance(raw, Mapping) or set(raw) != set(UPSTREAM_LEAD_FIELDS):
        raise GateError("147 lead row does not have the exact verified schema")
    row = {field: str(raw[field]) for field in UPSTREAM_LEAD_FIELDS}
    pair_id = row["pair_id"]
    if pair_id not in PAIR_IDENTITIES:
        raise GateError("147 lead belongs to an unknown pair")
    if (expected_collection == "primary") != (pair_id in {"A", "B"}):
        raise GateError("147 A/B and CONTROL lead families are not separated")
    trait1, trait2, _ = PAIR_IDENTITIES[pair_id]
    chromosome = _integer(row["chr"], "lead chromosome", minimum=1)
    position = _integer(row["position"], "lead position", minimum=1)
    locus_start = _integer(row["locus_start"], "lead locus start", minimum=1)
    locus_end = _integer(row["locus_end"], "lead locus end", minimum=locus_start)
    if chromosome > 22 or not locus_start <= position <= locus_end:
        raise GateError("147 lead has invalid autosomal locus coordinates")
    snp = row["lead_variant"]
    if not RSID.fullmatch(snp):
        raise GateError("147 lead variant is not an rsID")
    p_value = _probability(row["PLACO_P"], "lead PLACO P")
    q_value = _probability(row["FDR"], "lead PLACO FDR")
    _probability(row["trait1_P"], "lead trait1 P")
    _probability(row["trait2_P"], "lead trait2 P")
    _finite(row["lead_Z1"], "lead Z1")
    _finite(row["lead_Z2"], "lead Z2")
    if (
        row["analysis_id"] != UPSTREAM_ANALYSIS_ID
        or row["family_role"] != UPSTREAM_PAIR_ROLES[pair_id]
        or row["trait1"] != trait1
        or row["trait2"] != trait2
        or row["lead_A1"] not in BASES
        or row["lead_A2"] not in BASES
        or row["lead_A1"] == row["lead_A2"]
        or row["independent_signal"] != "TRUE"
        or row["annotations"]
        != "STATISTICAL_CROSS_TRAIT_PLEIOTROPY_ONLY;NO_SHARED_CAUSAL_VARIANT_NO_MEDIATION_NO_CAUSAL_DIRECTION"
        or row["evidence_labels"] != ";".join(_expected_labels(pair_id, p_value, q_value))
        or row["eligibility_basis"] != ";".join(_expected_eligibility(p_value, q_value))
        or row["clump_id"] != f"PLACO:{pair_id}:{snp}"
        or row["ld_reference_status"]
        != "EXACT_RSID_CHR_BP_UNORDERED_NONPALINDROMIC_ALLELES_MATCHED"
        or row["plink_clump_status"] != "INDEPENDENT_LEAD_REAL_PLINK_1_9"
        or row["lava_terminal_state"] != "TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY"
        or row["method_availability"] != METHOD_AVAILABILITY
        or row["claim_limit"] != (CONTROL_CLAIM if pair_id == "CONTROL" else PRIMARY_CLAIM)
    ):
        raise GateError("147 lead identity, eligibility, method, or claim fields drifted")
    expected_headline = (
        "NOT_APPLICABLE_CONTROL" if pair_id == "CONTROL"
        else str(p_value <= 2.5e-8).upper()
    )
    if row["primary_headline"] != expected_headline:
        raise GateError("147 lead primary-headline flag differs from its PLACO P")
    clump_evidence = _tokens(row["clump_evidence_ids"], "147 clump evidence")
    pair_block_leads = _tokens(row["pair_block_evidence_ids"], "147 pair/block leads")
    for evidence_id in [*clump_evidence, *pair_block_leads]:
        _parse_evidence_id(evidence_id, pair_id)
    lead_evidence_id = f"PLACO:{pair_id}:{snp}:{chromosome}:{position}"
    if lead_evidence_id not in clump_evidence:
        raise GateError("147 independent lead is absent from its ordered clump-evidence family")
    block = map_to_block(chromosome, position, block_index)
    expected_block_id = "UNMAPPED" if block is None else f"LOC{block.locus}"
    if row["lava_block_id"] != expected_block_id:
        raise GateError("147 lead official LAVA block mapping drifted")
    if block is None:
        if (
            row["lava_block_start"] != "NA"
            or row["lava_block_end"] != "NA"
            or row["lava_failed_qc_accounting"]
            != "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;OFFICIAL_BLOCK=UNMAPPED"
            or row["conditional_status"] != "BLOCKED_BY_OFFICIAL_LAVA_PARTITION_COVERAGE"
        ):
            raise GateError("147 unmapped lead lost its exact no-block accounting")
    else:
        if (
            row["lava_block_start"] != str(block.start)
            or row["lava_block_end"] != str(block.stop)
            or not row["lava_failed_qc_accounting"].startswith(
                "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;LOCUS_STATUS="
            )
            or "SCIENTIFIC_SUPPORT=FORBIDDEN" not in row["lava_failed_qc_accounting"]
            or not row["conditional_status"].startswith(
                "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;FAILURE_STATUSES="
            )
        ):
            raise GateError("147 mapped lead lost failed-QC LAVA accounting semantics")
    status = row["pairwise_r2_status"]
    maximum = row["max_pairwise_r2_within_1mb"]
    if status == "PASS_NO_WITHIN_WINDOW_LEAD_PAIR":
        if maximum != "NA":
            raise GateError("147 no-pair LD audit unexpectedly reports an r2 maximum")
    elif status == "PASS_STRICT_R2_LT_0.1_WITHIN_1_MB":
        observed = _finite(maximum, "147 maximum pairwise r2")
        if not 0 <= observed < 0.1:
            raise GateError("147 independent leads do not satisfy strict r2 < 0.1")
    else:
        raise GateError("147 lead lacks a passing independent-lead r2 audit")
    return {
        "row": row, "pair_id": pair_id, "chromosome": chromosome,
        "position": position, "snp": snp, "p": p_value, "q": q_value,
        "primary_headline": expected_headline == "TRUE", "block": block,
        "block_id": expected_block_id, "lead_evidence_id": lead_evidence_id,
        "clump_evidence": clump_evidence, "pair_block_leads": pair_block_leads,
    }


def validate_upstream_payload(
    payload: Any,
    block_index: Mapping[int, tuple[list[int], list[Block]]],
) -> dict[str, Any]:
    """Deeply validate the exact public value returned by script 147."""

    if not isinstance(payload, Mapping) or set(payload) != UPSTREAM_TOP_LEVEL_KEYS:
        raise GateError("147 verifier returned an unexpected top-level schema")
    state = payload.get("state")
    fingerprint = payload.get("fingerprint")
    if (
        payload.get("schema_version") != UPSTREAM_VERIFIER_SCHEMA
        or state not in UPSTREAM_STATES
        or not isinstance(fingerprint, str)
        or not HEX64.fullmatch(fingerprint)
    ):
        raise GateError("147 verifier identity or state drifted")
    zero = payload.get("zero_family")
    if zero != {"is_zero": state == UPSTREAM_ZERO, "state": state}:
        raise GateError("147 zero-family receipt is inconsistent with its state")
    if payload.get("method_availability") != {
        "PLACO": "COMPLETE", "CONJFDR": CONJFDR_BLOCKED,
    }:
        raise GateError("147 method availability does not prove complete PLACO and blocked conjFDR")
    _validate_archive_payload(payload)
    pair_states = _validate_complete_family(payload.get("complete_family"))
    primary_raw = payload.get("primary_leads")
    control_raw = payload.get("control_leads")
    if not isinstance(primary_raw, list) or not isinstance(control_raw, list):
        raise GateError("147 lead collections are not lists")
    leads = [
        *[_validate_lead(row, "primary", block_index) for row in primary_raw],
        *[_validate_lead(row, "control", block_index) for row in control_raw],
    ]
    expected_order = sorted(
        leads,
        key=lambda item: (
            PAIR_RANK[item["pair_id"]], item["chromosome"], item["position"], item["snp"],
        ),
    )
    if leads != expected_order:
        raise GateError("147 independent-lead family is not in canonical pair/coordinate order")
    clump_ids: set[str] = set()
    all_clump_evidence: set[str] = set()
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for lead in leads:
        if pair_states[lead["pair_id"]] != "COMPLETE_WITH_HITS":
            raise GateError("147 emitted a lead from a TESTED_NO_HIT pair")
        clump_id = lead["row"]["clump_id"]
        if clump_id in clump_ids:
            raise GateError("147 independent lead identity is duplicated")
        clump_ids.add(clump_id)
        overlap = all_clump_evidence.intersection(lead["clump_evidence"])
        if overlap:
            raise GateError("147 clump-evidence families overlap across independent leads")
        all_clump_evidence.update(lead["clump_evidence"])
        grouped[(lead["pair_id"], lead["block_id"])].append(lead)
    for key, family in grouped.items():
        ordered_leads = [lead["lead_evidence_id"] for lead in family]
        recorded_order = family[0]["pair_block_leads"]
        if len(recorded_order) != len(ordered_leads) or set(recorded_order) != set(ordered_leads):
            raise GateError(f"147 pair/block lead-evidence membership drifted for {key}")
        for lead in family:
            if lead["pair_block_leads"] != recorded_order:
                raise GateError(f"147 pair/block lead-evidence accounting drifted for {key}")

    coverage = payload.get("ld_reference_coverage")
    if not isinstance(coverage, Mapping) or set(coverage) != {
        "eligible_rows", "status_counts", "coverage_failure_rows",
    }:
        raise GateError("147 LD-reference coverage payload has an unexpected shape")
    eligible = _integer(coverage.get("eligible_rows"), "eligible evidence rows", minimum=0)
    failures = _integer(coverage.get("coverage_failure_rows"), "coverage failures", minimum=0)
    counts = coverage.get("status_counts")
    if not isinstance(counts, Mapping) or any(not isinstance(key, str) or not key for key in counts):
        raise GateError("147 LD-reference status counts are malformed")
    parsed_counts = {
        key: _integer(value, f"LD-reference status {key}", minimum=0)
        for key, value in counts.items()
    }
    if any(value == 0 for value in parsed_counts.values()):
        raise GateError("147 LD-reference status counts contain an impossible zero group")
    matched = parsed_counts.get("MATCHED", 0)
    if sum(parsed_counts.values()) != eligible or failures != eligible - matched:
        raise GateError("147 LD-reference denominator does not reconcile")
    lead_pairs = {lead["pair_id"] for lead in leads}
    hit_pairs = {pair for pair, pair_state in pair_states.items() if pair_state == "COMPLETE_WITH_HITS"}
    if state == UPSTREAM_ZERO:
        if eligible != 0 or failures != 0 or leads or hit_pairs:
            raise GateError("147 zero family contains eligible evidence or independent leads")
    elif state == UPSTREAM_BLOCKED:
        if eligible < 1 or matched != 0 or failures != eligible or leads or not hit_pairs:
            raise GateError("147 LD-reference-blocked state is not an exact all-failure family")
    else:
        if eligible < 1 or not leads:
            raise GateError("147 complete state lacks eligible evidence or independent leads")
        if failures != 0 or matched != eligible:
            raise UpstreamBlocked(
                "147 complete loci coexist with LD-reference coverage failures; the complete "
                "fine-mapping evidence denominator cannot be frozen"
            )
        if len(all_clump_evidence) != matched:
            raise GateError("147 clump-evidence union does not equal its MATCHED denominator")
        if lead_pairs != hit_pairs:
            raise GateError("147 independent-lead pairs do not equal COMPLETE_WITH_HITS pairs")
    return {
        "state": state, "fingerprint": fingerprint, "leads": leads,
        "matched_evidence_count": matched, "eligible_evidence_count": eligible,
        "coverage_failure_rows": failures,
    }


def _verify_archive_state_live(root: Path, upstream: Mapping[str, Any]) -> dict[str, object]:
    manager_identity = stable_identity(root, ARCHIVE_MANAGER_REL)
    manager = _load_module(
        "track_b_lava_archive_state_for_finemapping_145",
        safe_path(root, ARCHIVE_MANAGER_REL, "archive manager", must_exist=True),
    )
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        observed = manager.verify_reference_state(root)
    except BaseException as error:
        raise GateError(f"script 146 reference-state guard failed: {error}") from error
    finally:
        sys.dont_write_bytecode = previous
    if (
        not isinstance(observed, dict)
        or observed.get("state") not in ACCEPTED_ARCHIVE_STATES
        or observed != upstream.get("archive_reference")
        or observed.get("state") != upstream.get("archive_state")
    ):
        raise GateError("script 146 live reference state differs from the 147 verified snapshot")
    return manager_identity


def load_verified_upstream(
    root: Path = ROOT,
    block_index: Mapping[int, tuple[list[int], list[Block]]] | None = None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, object], dict[str, object]]:
    """Call only script 147's stable verifier; never open its raw PLACO ledgers here."""

    verifier_identity = stable_identity(root, PLEIOTROPY_VERIFIER_REL)
    if verifier_identity["sha256"] != UPSTREAM_VERIFIER_SHA256:
        raise GateError("script 147 differs from its pinned stable API implementation")
    verifier = _load_module(
        "track_b_pleiotropy_verified_family_for_finemapping_145",
        safe_path(root, PLEIOTROPY_VERIFIER_REL, "script 147", must_exist=True),
    )
    if not callable(getattr(verifier, "verify_results", None)):
        raise GateError("script 147 lacks verify_results(root=...) API")
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        payload = verifier.verify_results(root=root)
    except BaseException as error:
        raise GateError(f"script 147 verified-family API failed: {error}") from error
    finally:
        sys.dont_write_bytecode = previous
    if stable_identity(root, PLEIOTROPY_VERIFIER_REL) != verifier_identity:
        raise GateError("script 147 changed while its verified-family API was running")
    if block_index is None:
        _, block_index, _ = load_official_blocks(root)
    validated = validate_upstream_payload(payload, block_index)
    archive_manager_identity = _verify_archive_state_live(root, payload)
    return dict(payload), validated, verifier_identity, archive_manager_identity


def _priority_key(group: Sequence[Mapping[str, Any]]) -> tuple[Any, ...]:
    return (
        0 if any(bool(lead["primary_headline"]) for lead in group) else 1,
        min(float(lead["p"]) for lead in group),
        min(float(lead["q"]) for lead in group),
        int(group[0]["block"].chromosome),
        int(group[0]["block"].start),
        int(group[0]["block"].stop),
        tuple(str(lead["lead_evidence_id"]) for lead in group),
    )


def build_family_rows(validated: Mapping[str, Any]) -> dict[str, Any]:
    """Build the full pair-by-block family from already validated 147 rows."""

    state = validated.get("state")
    leads = validated.get("leads")
    if state not in UPSTREAM_STATES or not isinstance(leads, list):
        raise GateError("validated upstream family is malformed")
    if state == UPSTREAM_BLOCKED:
        raise UpstreamBlocked(
            "147 state BLOCKED_BY_LD_REFERENCE_COVERAGE is upstream-blocked, never zero"
        )
    if state == UPSTREAM_ZERO:
        return {
            "state": ZERO_STATE, "mapped_rows": [], "unavailable_rows": [],
            "mapped_evidence_ids": [], "unavailable_evidence_ids": [],
            "upstream_evidence_ids": [],
        }

    mapped_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    unavailable_leads: list[dict[str, Any]] = []
    for lead in leads:
        if lead["block"] is None:
            unavailable_leads.append(lead)
        else:
            mapped_groups[(lead["pair_id"], lead["block_id"])].append(lead)

    tier_by_group: dict[tuple[str, str], str] = {}
    for pair_id in ("A", "B"):
        candidates = [
            (key, family) for key, family in mapped_groups.items() if key[0] == pair_id
        ]
        ranked = sorted(candidates, key=lambda item: _priority_key(item[1]))
        primary_keys = {key for key, _ in ranked[:min(3, len(ranked))]}
        for key, _ in candidates:
            tier_by_group[key] = "PRIMARY_TIER" if key in primary_keys else "SECONDARY_TIER"
    for key in mapped_groups:
        if key[0] == "CONTROL":
            tier_by_group[key] = "CONTROL_TIER"

    mapped_rows: list[dict[str, str]] = []
    mapped_evidence_ids: list[str] = []
    group_order = sorted(
        mapped_groups,
        key=lambda key: (
            PAIR_RANK[key[0]], mapped_groups[key][0]["block"].chromosome,
            mapped_groups[key][0]["block"].start,
            mapped_groups[key][0]["block"].stop, key[1],
        ),
    )
    for key in group_order:
        pair_id, block_id = key
        family = mapped_groups[key]
        block = family[0]["block"]
        if block is None or any(lead["block"] != block for lead in family):
            raise GateError("mapped pair/block family has inconsistent official coordinates")
        trait1, trait2, role = PAIR_IDENTITIES[pair_id]
        lead_ids = [str(lead["lead_evidence_id"]) for lead in family]
        lead_variants = [str(lead["snp"]) for lead in family]
        clump_union = [
            evidence_id for lead in family for evidence_id in lead["clump_evidence"]
        ]
        if len(clump_union) != len(set(clump_union)) or not set(lead_ids).issubset(clump_union):
            raise GateError("mapped lead identities are not a subset of a disjoint clump union")
        if len(lead_ids) != len(lead_variants):
            raise GateError("mapped lead evidence IDs and lead variants are not parallel")
        if any(lead["row"]["method_availability"] != METHOD_AVAILABILITY for lead in family):
            raise GateError("mapped pair/block method availability is inconsistent")
        claim = CONTROL_CLAIM if pair_id == "CONTROL" else PRIMARY_CLAIM
        if any(lead["row"]["claim_limit"] != claim for lead in family):
            raise GateError("mapped pair/block claim limits are inconsistent")
        mapped_rows.append({
            "analysis_id": ANALYSIS_ID,
            "locus_entry_id": f"TBFM_{pair_id}_{block_id}",
            "pair_id": pair_id,
            "family_role": role,
            "trait1": trait1,
            "trait2": trait2,
            "CHR": str(block.chromosome),
            "START": str(block.start),
            "STOP": str(block.stop),
            "ld_block_id": block_id,
            "inclusion_sources": "PLACO",
            "priority_tier": tier_by_group[key],
            "local_evidence_ids": "NA",
            "local_lead_variants": "NA",
            "pleiotropy_lead_evidence_ids": ";".join(lead_ids),
            "pleiotropy_lead_variants": ";".join(lead_variants),
            "pleiotropy_clump_evidence_ids": ";".join(clump_union),
            "method_availability": METHOD_AVAILABILITY,
            "partition_role": PARTITION_ROLE,
            "selection_rule": SELECTION_RULE,
            "claim_limit": claim,
        })
        mapped_evidence_ids.extend(clump_union)

    unavailable_rows: list[dict[str, str]] = []
    unavailable_evidence_ids: list[str] = []
    for lead in sorted(
        unavailable_leads,
        key=lambda item: (
            PAIR_RANK[item["pair_id"]], item["chromosome"], item["position"],
            item["lead_evidence_id"],
        ),
    ):
        pair_id = str(lead["pair_id"])
        trait1, trait2, role = PAIR_IDENTITIES[pair_id]
        evidence = [str(value) for value in lead["clump_evidence"]]
        if str(lead["lead_evidence_id"]) not in evidence:
            raise GateError("unmapped independent lead is absent from its clump evidence")
        unavailable_id = (
            f"TBFM_UNAVAILABLE_{pair_id}_{lead['snp']}_{lead['chromosome']}_{lead['position']}"
        )
        if not SAFE_ID.fullmatch(unavailable_id):
            raise GateError("derived unavailable identity is unsafe")
        unavailable_rows.append({
            "analysis_id": ANALYSIS_ID,
            "unavailable_entry_id": unavailable_id,
            "pair_id": pair_id,
            "family_role": role,
            "trait1": trait1,
            "trait2": trait2,
            "CHR": str(lead["chromosome"]),
            "BP": str(lead["position"]),
            "lead_variant": str(lead["snp"]),
            "ld_block_id": "UNMAPPED",
            "evidence_ids": ";".join(evidence),
            "inclusion_sources": "PLACO",
            "terminal_status": "BLOCKED_BY_DATA",
            "blocker_reason": "NO_PREDECLARED_SIGNED_LD_BLOCK",
            "claim_limit": CONTROL_CLAIM if pair_id == "CONTROL" else PRIMARY_CLAIM,
        })
        unavailable_evidence_ids.extend(evidence)

    upstream_evidence_ids = [*mapped_evidence_ids, *unavailable_evidence_ids]
    if len(upstream_evidence_ids) != len(set(upstream_evidence_ids)):
        raise GateError("mapped and unavailable clump-evidence families are not disjoint")
    denominator = _integer(
        validated.get("matched_evidence_count"), "validated MATCHED denominator", minimum=0,
    )
    if len(upstream_evidence_ids) != denominator:
        raise GateError("mapped plus unavailable evidence does not equal the 147 MATCHED denominator")
    if not mapped_rows and not unavailable_rows:
        raise GateError("147 complete state produced no mapped or unavailable independent lead")
    return {
        "state": READY_STATE, "mapped_rows": mapped_rows,
        "unavailable_rows": unavailable_rows, "mapped_evidence_ids": mapped_evidence_ids,
        "unavailable_evidence_ids": unavailable_evidence_ids,
        "upstream_evidence_ids": upstream_evidence_ids,
    }


def _upstream_provenance_identity(root: Path) -> dict[str, object]:
    return stable_identity(root, PLEIOTROPY_RESULT_PROVENANCE_REL)


def assemble_expected(root: Path = ROOT) -> dict[str, Any]:
    """Assemble exact bytes from stable APIs without publishing anything."""

    policy, policy_identity = validate_policy(root)
    del policy
    _, block_index, blocks_identity = load_official_blocks(root)
    upstream, validated, verifier_identity, archive_manager_identity = load_verified_upstream(
        root, block_index,
    )
    del upstream
    if validated["state"] == UPSTREAM_BLOCKED:
        return {
            "publication_allowed": False, "upstream_state": UPSTREAM_BLOCKED,
            "upstream_fingerprint": validated["fingerprint"],
            "reason": "UPSTREAM_LD_REFERENCE_COVERAGE_BLOCKED_NOT_ZERO",
        }
    built = build_family_rows(validated)
    mapped_rows = built["mapped_rows"]
    unavailable_rows = built["unavailable_rows"]
    manifest_content = tsv_bytes(FAMILY_FIELDS, mapped_rows)
    unavailable_content = (
        tsv_bytes(UNAVAILABLE_FIELDS, unavailable_rows) if unavailable_rows else None
    )
    zero_content: bytes | None = None
    if built["state"] == ZERO_STATE:
        zero_payload = {
            "schema_version": ZERO_SCHEMA,
            "analysis_id": ANALYSIS_ID,
            "artifact_role": "SEALED_PRE_MATERIALIZATION_ZERO_FAMILY_NOT_SCIENTIFIC_RESULT",
            "state": ZERO_STATE,
            "upstream_state": UPSTREAM_ZERO,
            "upstream_fingerprint": validated["fingerprint"],
            "upstream_matched_evidence_count": 0,
            "mapped_locus_count": 0,
            "unavailable_count": 0,
            "local_lava_evidence_used": False,
            "science_outputs_10_11_12_created": False,
            "header_only_science_outputs_created": False,
        }
        zero_content = json.dumps(zero_payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"

    manifest_identity = identity_for_content(FAMILY_REL, manifest_content)
    unavailable_identity = (
        identity_for_content(UNAVAILABLE_REL, unavailable_content)
        if unavailable_content is not None else None
    )
    zero_identity = (
        identity_for_content(ZERO_PROVENANCE_REL, zero_content)
        if zero_content is not None else None
    )
    mapped_evidence = list(built["mapped_evidence_ids"])
    unavailable_evidence = list(built["unavailable_evidence_ids"])
    upstream_evidence = list(built["upstream_evidence_ids"])
    lock = {
        "schema_version": SCHEMA,
        "analysis_id": ANALYSIS_ID,
        "artifact_role": ARTIFACT_ROLE,
        "state": built["state"],
        "results_accessed_before_lock": False,
        "pair_order": list(PAIR_ORDER),
        "manifest": manifest_identity,
        "locus_count": len(mapped_rows),
        "locus_entry_ids_in_order": [row["locus_entry_id"] for row in mapped_rows],
        "evidence_count": len(mapped_evidence),
        "evidence_ids_sha256": digest_json(mapped_evidence),
        "policy": policy_identity,
        "official_blocks": blocks_identity,
        "upstream_family_provenance": {
            "post_placo_results": _upstream_provenance_identity(root),
        },
        "generator_scripts": {
            "family_freezer": stable_identity(root, SCRIPT_REL),
            "archive_state_guard": archive_manager_identity,
            "pleiotropy_family_verifier": verifier_identity,
        },
        "selection_contract": dict(SELECTION_CONTRACT),
        "zero_family_provenance": zero_identity,
        "unavailable_manifest": unavailable_identity,
        "unavailable_count": len(unavailable_rows),
        "unavailable_entry_ids_in_order": [
            row["unavailable_entry_id"] for row in unavailable_rows
        ],
        "unavailable_evidence_count": len(unavailable_evidence),
        "unavailable_evidence_ids_sha256": digest_json(unavailable_evidence),
        "upstream_evidence_count": len(upstream_evidence),
        "upstream_evidence_ids_sha256": digest_json(upstream_evidence),
    }
    if set(lock) != LOCK_KEYS:
        raise GateError("internal corrected lock schema drifted")
    lock_content = json.dumps(lock, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    return {
        "publication_allowed": True,
        "upstream_state": validated["state"],
        "upstream_fingerprint": validated["fingerprint"],
        "state": built["state"],
        "mapped_rows": mapped_rows,
        "unavailable_rows": unavailable_rows,
        "upstream_evidence_ids": upstream_evidence,
        "manifest_content": manifest_content,
        "unavailable_content": unavailable_content,
        "zero_content": zero_content,
        "lock": lock,
        "lock_content": lock_content,
    }


def _artifact_present(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def assert_publication_boundary_clear(root: Path = ROOT) -> None:
    offenders = [
        str(relative) for relative in sorted(DOWNSTREAM_OR_COMPETING_PATHS, key=str)
        if _artifact_present(safe_path(root, relative, "downstream publication"))
    ]
    if offenders:
        raise GateError(
            "fine-mapping materialization/result artifact exists before family freeze: "
            + ",".join(offenders)
        )


def assert_preflight_targets_absent(root: Path = ROOT) -> None:
    offenders = [
        str(relative) for relative in sorted(OWN_PATHS, key=str)
        if _artifact_present(safe_path(root, relative, "pre-materialization output"))
    ]
    if offenders:
        raise GateError("pre-materialization output already exists: " + ",".join(offenders))


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _ensure_parent(root: Path, relative: Path) -> Path:
    path = safe_path(root, relative, "publication")
    current = root
    for part in relative.parent.parts:
        current /= part
        if current.exists():
            observed = os.lstat(current)
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISDIR(observed.st_mode):
                raise GateError(f"publication parent is not a real directory: {current}")
        else:
            os.mkdir(current, 0o755)
            _fsync_directory(current.parent)
    return path


def publish_no_replace(root: Path, relative: Path, content: bytes) -> dict[str, object]:
    """Publish exact bytes by exclusive hard-link; accept only an identical partial."""

    expected = identity_for_content(relative, content)
    path = _ensure_parent(root, relative)
    if _artifact_present(path):
        observed = stable_identity(root, relative)
        if not identity_matches(observed, expected):
            raise GateError(f"existing no-replace artifact differs: {relative}")
        return observed
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{os.urandom(8).hex()}.tmp")
    created_temporary = False
    try:
        with temporary.open("xb") as handle:
            created_temporary = True
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError:
            observed = stable_identity(root, relative)
            if not identity_matches(observed, expected):
                raise GateError(f"concurrent no-replace publication differs: {relative}")
        _fsync_directory(path.parent)
    finally:
        if created_temporary:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
    observed = stable_identity(root, relative)
    if not identity_matches(observed, expected):
        raise GateError(f"published artifact differs immediately after commit: {relative}")
    return observed


def _same_expected(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    keys = (
        "publication_allowed", "upstream_state", "upstream_fingerprint", "state",
        "manifest_content", "unavailable_content", "zero_content", "lock_content",
    )
    return all(left.get(key) == right.get(key) for key in keys)


def preflight(root: Path = ROOT) -> dict[str, Any]:
    """Validate the publication boundary and stable upstream API without writes."""

    assert_publication_boundary_clear(root)
    assert_preflight_targets_absent(root)
    expected = assemble_expected(root)
    assert_publication_boundary_clear(root)
    assert_preflight_targets_absent(root)
    if not expected.get("publication_allowed"):
        return {
            "state": UPSTREAM_BLOCKED,
            "publication_allowed": False,
            "upstream_fingerprint": expected["upstream_fingerprint"],
            "mapped_count": 0,
            "unavailable_count": 0,
            "science_outputs_10_11_12_created": False,
        }
    return {
        "state": expected["state"],
        "publication_allowed": True,
        "upstream_state": expected["upstream_state"],
        "upstream_fingerprint": expected["upstream_fingerprint"],
        "mapped_count": len(expected["mapped_rows"]),
        "unavailable_count": len(expected["unavailable_rows"]),
        "upstream_evidence_denominator": len(expected["upstream_evidence_ids"]),
        "science_outputs_10_11_12_created": False,
    }


def freeze_family(root: Path = ROOT) -> dict[str, Any]:
    """Exclusively publish the pre-materialization family, with provenance last."""

    if _artifact_present(safe_path(root, FAMILY_LOCK_REL, "family lock")):
        return verify_results(root)
    assert_publication_boundary_clear(root)
    expected = assemble_expected(root)
    if not expected.get("publication_allowed"):
        raise UpstreamBlocked(
            "147 is BLOCKED_BY_LD_REFERENCE_COVERAGE; no zero or runnable family was published"
        )
    if expected["unavailable_content"] is None and _artifact_present(
        safe_path(root, UNAVAILABLE_REL, "unavailable ledger")
    ):
        raise GateError("unexpected unavailable ledger exists for an empty unavailable family")
    if expected["zero_content"] is None and _artifact_present(
        safe_path(root, ZERO_PROVENANCE_REL, "zero receipt")
    ):
        raise GateError("zero-family receipt exists for a nonzero family")
    assert_publication_boundary_clear(root)
    publish_no_replace(root, FAMILY_REL, expected["manifest_content"])
    if expected["unavailable_content"] is not None:
        assert_publication_boundary_clear(root)
        publish_no_replace(root, UNAVAILABLE_REL, expected["unavailable_content"])
    if expected["zero_content"] is not None:
        assert_publication_boundary_clear(root)
        publish_no_replace(root, ZERO_PROVENANCE_REL, expected["zero_content"])

    # Re-run the stable upstream verifier and all identities at the final commit
    # boundary. The lock is the only usability marker and is linked last.
    assert_publication_boundary_clear(root)
    final_expected = assemble_expected(root)
    if not final_expected.get("publication_allowed") or not _same_expected(expected, final_expected):
        raise GateError("pre-materialization inputs changed before provenance-last publication")
    assert_publication_boundary_clear(root)
    publish_no_replace(root, FAMILY_LOCK_REL, expected["lock_content"])
    return verify_results(root)


def _verify_absence(root: Path, relative: Path, label: str) -> None:
    if _artifact_present(safe_path(root, relative, label)):
        raise GateError(f"unexpected {label} exists: {relative}")


def verify_results(root: Path = ROOT) -> dict[str, Any]:
    """Side-effect-free stable API for scripts 148--151.

    Deliberately no competitor/result-absence check is performed here: once
    this immutable family is frozen, authorized downstream materialization and
    canonical results must not invalidate its verifier.
    """

    expected = assemble_expected(root)
    if not expected.get("publication_allowed"):
        raise UpstreamBlocked(
            "147 is BLOCKED_BY_LD_REFERENCE_COVERAGE; no frozen family is verifiable"
        )
    lock, lock_identity = read_json(root, FAMILY_LOCK_REL)
    if set(lock) != LOCK_KEYS or lock != expected["lock"]:
        raise GateError("DEEP_LOCUS_MANIFEST lock differs from current exact upstream family")
    rows, family_identity = read_tsv(
        root, FAMILY_REL, FAMILY_FIELDS, allow_header_only=(expected["state"] == ZERO_STATE),
    )
    if rows != expected["mapped_rows"] or not identity_matches(family_identity, lock["manifest"]):
        raise GateError("DEEP_LOCUS_MANIFEST differs from its provenance-last lock")
    unavailable_rows: list[dict[str, str]] = []
    unavailable_identity: dict[str, object] | None = None
    if expected["unavailable_content"] is not None:
        unavailable_rows, unavailable_identity = read_tsv(
            root, UNAVAILABLE_REL, UNAVAILABLE_FIELDS, allow_header_only=False,
        )
        if (
            unavailable_rows != expected["unavailable_rows"]
            or not identity_matches(unavailable_identity, lock["unavailable_manifest"])
        ):
            raise GateError("DEEP_LOCUS_UNAVAILABLE differs from the main family lock")
    else:
        _verify_absence(root, UNAVAILABLE_REL, "unavailable ledger")
    zero_identity: dict[str, object] | None = None
    if expected["zero_content"] is not None:
        content, zero_identity = stable_bytes(root, ZERO_PROVENANCE_REL)
        if content != expected["zero_content"] or not identity_matches(
            zero_identity, lock["zero_family_provenance"],
        ):
            raise GateError("pre-materialization zero-family receipt drifted")
    else:
        _verify_absence(root, ZERO_PROVENANCE_REL, "pre-materialization zero receipt")
    return {
        "state": expected["state"],
        "rows": rows,
        "unavailable_rows": unavailable_rows,
        "family_path": str(FAMILY_REL),
        "family_lock_path": str(FAMILY_LOCK_REL),
        "family_identity": family_identity,
        "family_lock_identity": lock_identity,
        "unavailable_path": str(UNAVAILABLE_REL) if unavailable_identity is not None else None,
        "unavailable_lock_path": None,
        "unavailable_identity": unavailable_identity,
        "unavailable_lock_identity": None,
        "zero_provenance_path": (
            str(ZERO_PROVENANCE_REL) if zero_identity is not None else None
        ),
        "zero_provenance_identity": zero_identity,
        "mapped_count": len(rows),
        "unavailable_count": len(unavailable_rows),
        "upstream_evidence_denominator": lock["upstream_evidence_count"],
        "family_sha256": family_identity["sha256"],
        "unavailable_sha256": (
            unavailable_identity["sha256"] if unavailable_identity is not None else None
        ),
    }


verify_frozen_family = verify_results
deep_verify_published_family = verify_results


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--freeze", action="store_true")
    action.add_argument("--seal", action="store_true", help="alias for --freeze")
    action.add_argument("--run", action="store_true", help="alias for --freeze")
    action.add_argument("--verify", action="store_true")
    parser.add_argument("--fingerprint", help="require this exact script-147 family fingerprint")
    parser.add_argument(
        "--execute", action="store_true",
        help="required for --freeze, --seal, or --run; never required for read-only checks",
    )
    arguments = parser.parse_args(argv)
    try:
        mutating = arguments.freeze or arguments.seal or arguments.run
        if mutating and not arguments.execute:
            raise GateError("explicit --execute is required for every family mutation")
        if arguments.fingerprint:
            fingerprint_check = assemble_expected(ROOT)
            if arguments.fingerprint != fingerprint_check.get("upstream_fingerprint"):
                raise GateError("requested fingerprint differs from the verified script-147 family")
        if arguments.preflight:
            payload = preflight(ROOT)
        elif arguments.verify:
            payload = verify_results(ROOT)
        else:
            payload = freeze_family(ROOT)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    except (GateError, OSError, ValueError, KeyError) as error:
        print(f"TRACK_B_FINEMAPPING_FAMILY_ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
