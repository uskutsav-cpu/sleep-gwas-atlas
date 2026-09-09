#!/usr/bin/env python3
"""Run the frozen P14/P16 human brain broad-class MAGMA analysis.

The command is deliberately inert until an exact upstream G4/G5 brain-PASS
eligible locus/gene family validates.  It never selects loci, genes, classes,
or thresholds from cell-enrichment results.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import stat
import subprocess
import sys
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence

sys.dont_write_bytecode = True

ROOT_AT_IMPORT = Path(__file__).resolve().parents[1]
POLICY_REL = "config/track_b_mechanism_followup_policy_v4.json"
CONTRACT_LOCK_REL = "results/track_b/mechanism_followup/v4/P14_P16_EXECUTION_AMENDMENT.lock.json"
MAPPING_REL = "results/track_b/mechanism_followup/v4/ALLEN_MTG_75_TO_7_CLASS_MAP.tsv"
ALIGNMENT_REL = "results/track_b/mechanism_followup/v4/ALLEN_GSE67835_CLASS_ALIGNMENT.tsv"
AMENDMENT_REL = "results/track_b/mechanism_followup/v4/P14_P16_EXECUTION_AMENDMENT.json"
READINESS_REL = "results/track_b/mechanism_followup/v4/P14_P16_READINESS.tsv"
REPORT_REL = "results/track_b/mechanism_followup/v4/P14_P16_EXECUTION_AMENDMENT_REPORT.md"
SUPERSESSION_REL = "results/track_b/mechanism_followup/PRE_RESULT_READINESS_SUPERSESSION_V4.json"
CONTRACT_OUTPUT_RELS = (
    MAPPING_REL, ALIGNMENT_REL, AMENDMENT_REL, READINESS_REL, REPORT_REL,
    SUPERSESSION_REL,
)
SELF_REL = "scripts/164_run_track_b_cell_class_enrichment.py"
CLASSIFIER_REL = "scripts/165_classify_track_b_cell_replication.py"
BUILDER_REL = "scripts/163_build_track_b_cell_class_contract.py"
PARENT_LOCK_REL = "results/track_b/mechanism_followup/v3/PRE_RESULT_READINESS.lock.json"
PARENT_LOCK_SHA256 = "28ba9e3b4be939657391f9ad9e90e89fdacd6e2be834e7602ced3e67c72eae0e"
POLICY_SCHEMA = "track-b-p14-p16-execution-policy.1"
LOCK_SCHEMA = "track-b-p14-p16-execution-amendment.1"
GATE_SCHEMA = "track-b-p14-p16-upstream-gate.1"
GATE_STATUS = "ELIGIBLE_FOR_BRAIN_CELL_CLASS"
UPSTREAM_SOURCE_SCHEMA = "track-b-p14-p16-upstream-eligibility-source.1"
UPSTREAM_SOURCE_KIND = "G4_G5_BRAIN_TISSUE_ELIGIBILITY_SOURCE"
UPSTREAM_SOURCE_STATUS = "COMPLETE_WITH_FAILURES_PRESERVED"
ELIGIBLE_FIELDS = [
    "pair_id", "target_trait_id", "locus_id", "gene_id",
    "genetic_evidence_tier", "brain_tissue_gate", "eligibility_status",
]
CLASS_ORDER = (
    "ASTROCYTES", "ENDOTHELIAL", "EXCITATORY_NEURONS",
    "INHIBITORY_NEURONS", "MICROGLIA", "OLIGODENDROCYTES", "OPC",
)
PREFIX_TO_CLASS = (
    ("Astro_", "ASTROCYTES"), ("Endo_", "ENDOTHELIAL"),
    ("Exc_", "EXCITATORY_NEURONS"), ("Inh_", "INHIBITORY_NEURONS"),
    ("Micro_", "MICROGLIA"), ("Oligo_", "OLIGODENDROCYTES"),
    ("OPC_", "OPC"),
)
EXPECTED_CLASS_COUNTS = {
    "ASTROCYTES": 2, "ENDOTHELIAL": 1, "EXCITATORY_NEURONS": 24,
    "INHIBITORY_NEURONS": 45, "MICROGLIA": 1,
    "OLIGODENDROCYTES": 1, "OPC": 1,
}
GSE_ORDER = (
    "astrocytes", "endothelial", "hybrid", "microglia", "neurons",
    "oligodendrocytes", "OPC",
)
MAPPING_FIELDS = [
    "mapping_version", "dataset_id", "source_position", "source_column",
    "broad_class", "class_member_index", "class_member_count",
    "aggregation_rule",
]
ALIGNMENT_FIELDS = [
    "alignment_version", "gse_source_position", "gse_class",
    "allen_broad_classes", "alignment_status", "replication_resolution",
    "gse_bh_denominator_inclusion", "independent_replication_credit_rule",
]
RESULT_FIELDS = [
    "dataset_id", "target_trait_id", "class_id", "source_column",
    "aggregation_rule", "member_count", "analysis_status", "n_genes",
    "beta", "se", "z", "p_value", "p_value_for_bh", "bh_fdr",
    "direction_pass",
]
GENE_RE = re.compile(r"ENSG[0-9]{11}")
SAFE_TOKEN_RE = re.compile(r"[A-Za-z0-9_.:+-]+")


class ContractError(RuntimeError):
    pass


@dataclass(frozen=True)
class FileRecord:
    path: str
    bytes: int
    sha256: str
    identity: tuple[int, int, int, int, int, int]


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def safe_relative(value: str) -> str:
    pure = PurePosixPath(value)
    if (
        not value or pure.is_absolute() or "\\" in value
        or any(part in {"", ".", ".."} for part in pure.parts)
    ):
        raise ContractError(f"unsafe relative path: {value!r}")
    return pure.as_posix()


def _identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        value.st_dev, value.st_ino, value.st_mode, value.st_size,
        value.st_mtime_ns, value.st_ctime_ns,
    )


def _check_parents(root: Path, relative: str) -> None:
    current = root
    for part in PurePosixPath(safe_relative(relative)).parts[:-1]:
        current = current / part
        info = os.lstat(current)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise ContractError(f"non-directory or symlink in path: {relative}")


def stable_file(root: Path, relative: str) -> FileRecord:
    relative = safe_relative(relative)
    _check_parents(root, relative)
    path = root / relative
    before_path = os.lstat(path)
    if stat.S_ISLNK(before_path.st_mode) or not stat.S_ISREG(before_path.st_mode):
        raise ContractError(f"not a real regular file: {relative}")
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        digest = hashlib.sha256()
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    if _identity(before_path) != _identity(before) or _identity(before) != _identity(after) or _identity(after) != _identity(final):
        raise ContractError(f"file drifted while hashing: {relative}")
    return FileRecord(relative, after.st_size, digest.hexdigest(), _identity(after))


def stable_absolute_file(path: Path) -> FileRecord:
    """Hash an owned temporary file without buffering it or following symlinks."""
    before_path = os.lstat(path)
    if stat.S_ISLNK(before_path.st_mode) or not stat.S_ISREG(before_path.st_mode):
        raise ContractError(f"not a real regular file: {path}")
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        digest = hashlib.sha256()
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    if (
        _identity(before_path) != _identity(before)
        or _identity(before) != _identity(after)
        or _identity(after) != _identity(final)
    ):
        raise ContractError(f"temporary file drifted while hashing: {path}")
    return FileRecord(str(path), after.st_size, digest.hexdigest(), _identity(after))


def read_absolute_bytes(path: Path, *, limit: int = 64 * 1024 * 1024) -> bytes:
    record = stable_absolute_file(path)
    if record.bytes > limit:
        raise ContractError(f"refusing to buffer oversized runtime output: {path}")
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        chunks = []
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            chunks.append(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    payload = b"".join(chunks)
    if (
        _identity(before) != record.identity or _identity(after) != record.identity
        or _identity(final) != record.identity or len(payload) != record.bytes
        or hashlib.sha256(payload).hexdigest() != record.sha256
    ):
        raise ContractError(f"runtime output drifted while reading: {path}")
    return payload


def read_stable_bytes(root: Path, relative: str, *, limit: int = 64 * 1024 * 1024) -> bytes:
    record = stable_file(root, relative)
    if record.bytes > limit:
        raise ContractError(f"refusing to buffer large file: {relative}")
    path = root / relative
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        chunks: list[bytes] = []
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            chunks.append(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    payload = b"".join(chunks)
    if (
        _identity(before) != record.identity or _identity(after) != record.identity
        or _identity(final) != record.identity or len(payload) != record.bytes
        or hashlib.sha256(payload).hexdigest() != record.sha256
    ):
        raise ContractError(f"file drifted while reading: {relative}")
    return payload


def read_canonical_json(root: Path, relative: str) -> dict[str, object]:
    payload = read_stable_bytes(root, relative)
    try:
        value = json.loads(payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"invalid JSON: {relative}: {exc}") from exc
    if not isinstance(value, dict) or canonical_json(value) != payload:
        raise ContractError(f"JSON is not a canonical object: {relative}")
    return value


def tsv_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: str(row.get(field, "")).replace("\t", " ").replace("\n", " ").replace("\r", " ") for field in fields})
    return output.getvalue().encode()


def parse_tsv_payload(payload: bytes, fields: Sequence[str], label: str) -> list[dict[str, str]]:
    try:
        text = payload.decode("utf-8")
    except UnicodeError as exc:
        raise ContractError(f"{label} is not UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames != list(fields):
        raise ContractError(f"{label} has an unexpected header")
    rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ContractError(f"{label} has a malformed row")
    return rows


def class_mapping_rows(columns: Sequence[str]) -> list[dict[str, object]]:
    if len(columns) != 75 or len(set(columns)) != 75 or "Average" in columns:
        raise ContractError("Allen MTG source must contain exactly 75 unique subtype columns excluding Average")
    grouped: dict[str, list[str]] = {name: [] for name in CLASS_ORDER}
    mapped: list[tuple[str, str]] = []
    for column in columns:
        matches = [broad for prefix, broad in PREFIX_TO_CLASS if column.startswith(prefix)]
        if len(matches) != 1:
            raise ContractError(f"Allen MTG subtype does not map exactly once: {column}")
        grouped[matches[0]].append(column)
        mapped.append((column, matches[0]))
    counts = {name: len(values) for name, values in grouped.items()}
    if counts != EXPECTED_CLASS_COUNTS:
        raise ContractError(f"Allen MTG class counts differ from the frozen family: {counts}")
    indexes = {name: 0 for name in CLASS_ORDER}
    rows: list[dict[str, object]] = []
    for position, (column, broad) in enumerate(mapped, 1):
        indexes[broad] += 1
        rows.append({
            "mapping_version": "ALLEN_MTG_75_TO_7_V1",
            "dataset_id": "Allen_Human_MTG_level2",
            "source_position": position,
            "source_column": column,
            "broad_class": broad,
            "class_member_index": indexes[broad],
            "class_member_count": counts[broad],
            "aggregation_rule": "ARITHMETIC_MEAN_ALL_MAPPED_SUBTYPES_MATH_FSUM_SOURCE_ORDER",
        })
    return rows


def alignment_rows() -> list[dict[str, object]]:
    alignments = (
        ("astrocytes", "ASTROCYTES", "ONE_TO_ONE", "BROAD_CLASS", "ONE_IF_REPLICATED"),
        ("endothelial", "ENDOTHELIAL", "ONE_TO_ONE", "BROAD_CLASS", "ONE_IF_REPLICATED"),
        ("hybrid", "NA", "UNALIGNED_HYBRID_PRESERVED", "UNAVAILABLE", "ZERO"),
        ("microglia", "MICROGLIA", "ONE_TO_ONE", "BROAD_CLASS", "ONE_IF_REPLICATED"),
        ("neurons", "EXCITATORY_NEURONS;INHIBITORY_NEURONS", "MANY_TO_ONE", "COARSE_NEURON_ONLY", "ONE_SHARED_TARGET_MAXIMUM"),
        ("oligodendrocytes", "OLIGODENDROCYTES", "ONE_TO_ONE", "BROAD_CLASS", "ONE_IF_REPLICATED"),
        ("OPC", "OPC", "ONE_TO_ONE", "BROAD_CLASS", "ONE_IF_REPLICATED"),
    )
    return [{
        "alignment_version": "ALLEN_GSE67835_7_CLASS_V1",
        "gse_source_position": index,
        "gse_class": gse,
        "allen_broad_classes": allen,
        "alignment_status": status,
        "replication_resolution": resolution,
        "gse_bh_denominator_inclusion": "YES",
        "independent_replication_credit_rule": credit,
    } for index, (gse, allen, status, resolution, credit) in enumerate(alignments, 1)]


def bh_adjust(labels: Sequence[str], p_values: Sequence[float]) -> dict[str, float]:
    if len(labels) != 7 or len(p_values) != 7 or len(set(labels)) != 7:
        raise ContractError("BH requires the exact unique seven-hypothesis class family")
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in p_values):
        raise ContractError("BH received an invalid probability")
    ordered = sorted(zip(p_values, labels), key=lambda value: (value[0], value[1]))
    adjusted: dict[str, float] = {}
    running = 1.0
    for rank in range(7, 0, -1):
        p_value, label = ordered[rank - 1]
        running = min(running, p_value * 7 / rank, 1.0)
        adjusted[label] = running
    return adjusted


def load_policy(root: Path) -> tuple[dict[str, object], FileRecord]:
    policy = read_canonical_json(root, POLICY_REL)
    required = {
        "schema_version", "contract_revision", "parent_v3", "amended_analysis_ids",
        "preserved_analysis_ids", "static_inputs", "source_datasets", "executors",
        "mapping", "statistics", "upstream_gate", "production_outputs",
        "result_blind", "future_results_accessed",
    }
    if (
        set(policy) != required or policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("contract_revision") != 4
        or policy.get("result_blind") is not True
        or policy.get("future_results_accessed") is not False
        or policy.get("amended_analysis_ids") != [
            "P14_BRAIN_CELL_CLASS_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION",
        ]
        or policy.get("preserved_analysis_ids") != [
            "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY", "P15_BRAIN_CELL_SUBTYPE_DISCOVERY",
            "P16_HUMAN_CELL_SUBTYPE_REPLICATION",
        ]
    ):
        raise ContractError("V4 P14/P16 policy schema or scope is invalid")
    parent = policy.get("parent_v3")
    if not isinstance(parent, dict) or parent != {
        "lock_path": PARENT_LOCK_REL, "lock_sha256": PARENT_LOCK_SHA256,
        "retention": "PRESERVE_V1_V2_V3_BYTE_EXACT",
    }:
        raise ContractError("V4 parent lineage differs from the frozen V3 lock")
    statistics = policy.get("statistics")
    if not isinstance(statistics, dict) or statistics != {
        "tested_statistic": "MAGMA_V1.10_GENE_PROPERTY_BETA_OVER_SE",
        "model": "condition-hide=Average direction=greater",
        "average_role": "COVARIATE_ONLY_EXCLUDED_FROM_SEVEN_HYPOTHESES",
        "allen_bh_family": list(CLASS_ORDER),
        "gse_bh_family": list(GSE_ORDER),
        "bh_denominator_per_atlas_per_trait": 7,
        "fdr_alpha": 0.05,
        "direction_rule": "BETA_GT_ZERO_IN_BOTH_ATLASES",
        "failed_p_for_bh": 1.0,
        "zero_p_rule": "VALID_AND_PRESERVED",
        "na_rule": "PRESERVE_ROW_MARK_FAILED_AND_USE_P_EQ_1_FOR_BH_ONLY",
        "replication_rule": "BOTH_ATLAS_FDR_LE_0.05_AND_DIRECTION_CONCORDANT",
        "many_to_one_rule": "EXC_AND_INH_SHARE_GSE_NEURONS_MAX_PARTIAL_COARSE_NEURON_ONE_SHARED_TARGET_NO_DOUBLE_CREDIT",
        "hybrid_rule": "TEST_AND_INCLUDE_IN_GSE_BH_DENOMINATOR_BUT_NEVER_ALIGN",
    }:
        raise ContractError("V4 statistical contract differs from the hard-coded model")
    mapping = policy.get("mapping")
    if not isinstance(mapping, dict) or mapping != {
        "allen_mapping_version": "ALLEN_MTG_75_TO_7_V1",
        "allen_source_subtype_count": 75,
        "allen_broad_class_order": list(CLASS_ORDER),
        "allen_expected_class_counts": EXPECTED_CLASS_COUNTS,
        "allen_aggregation_rule": "ARITHMETIC_MEAN_ALL_MAPPED_SUBTYPES_MATH_FSUM_SOURCE_ORDER",
        "allen_average_rule": "PRESERVE_SOURCE_AVERAGE_AS_COVARIATE_ONLY",
        "alignment_version": "ALLEN_GSE67835_7_CLASS_V1",
        "gse_class_order": list(GSE_ORDER),
        "hybrid_alignment": "UNALIGNED_HYBRID_PRESERVED",
        "neurons_alignment": "MANY_TO_ONE_EXCITATORY_AND_INHIBITORY_WITHOUT_DOUBLE_CREDIT",
    }:
        raise ContractError("V4 mapping contract differs from the hard-coded complete families")
    gate = policy.get("upstream_gate")
    if not isinstance(gate, dict) or gate != {
        "gate_schema": GATE_SCHEMA,
        "gate_terminal_status": GATE_STATUS,
        "source_contract_schema": UPSTREAM_SOURCE_SCHEMA,
        "source_contract_kind": UPSTREAM_SOURCE_KIND,
        "source_contract_terminal_status": UPSTREAM_SOURCE_STATUS,
        "pair_manifest_path": "results/track_b/pair_manifest.tsv",
        "analysis_panel_path": "config/analysis_panel.tsv",
        "complete_candidate_fields": ELIGIBLE_FIELDS,
        "eligible_family_fields": ELIGIBLE_FIELDS,
        "eligible_genetic_tiers": ["G4", "G5"],
        "selection_rule": "PREDECLARED_G4_G5_AND_BRAIN_TISSUE_PASS_ONLY",
        "maximum_family_bytes": 67108864,
        "gene_results_path_template": "work/interpretation_automatic/fuma_scrna/gene_results/{trait_id}/{trait_id}.genes.raw",
        "gene_results_provenance_template": "work/interpretation_automatic/fuma_scrna/gene_results/{trait_id}/provenance.json",
        "minimum_gene_rows": 10000,
        "magma_timeout_sec": 86400,
    }:
        raise ContractError("V4 upstream gate differs from the hard-coded complete-family rule")
    production = policy.get("production_outputs")
    if not isinstance(production, dict) or production != {
        "namespace_prefix": "results/track_b/mechanism_followup/science/cell_class/",
        "enrichment_dir_template": "results/track_b/mechanism_followup/science/cell_class/{pair_id}/{trait_id}/enrichment",
        "replication_dir_template": "results/track_b/mechanism_followup/science/cell_class/{pair_id}/{trait_id}/replication",
        "enrichment_required_files": [
            "P14_ALLEN_BROAD_CLASS_RESULTS.tsv",
            "P16_GSE67835_CLASS_RESULTS.tsv",
            "Allen_Human_MTG_level2.gsa.out",
            "Allen_Human_MTG_level2.log.out",
            "Allen_Human_MTG_level2.stdout.txt",
            "Allen_Human_MTG_level2.stderr.txt",
            "GSE67835_Human_Cortex_woFetal.gsa.out",
            "GSE67835_Human_Cortex_woFetal.log.out",
            "GSE67835_Human_Cortex_woFetal.stdout.txt",
            "GSE67835_Human_Cortex_woFetal.stderr.txt",
            "P14_P16_ENRICHMENT.provenance.json",
        ],
        "replication_required_files": [
            "P16_CLASS_REPLICATION.tsv", "P16_CLASS_REPLICATION_SUMMARY.json",
            "P16_CLASS_REPLICATION.provenance.json",
        ],
        "publication_rule": "NO_REPLACE_DIRECTORY_WITH_PROVENANCE_LAST",
    }:
        raise ContractError("V4 production output family differs from the hard-coded no-replace family")
    return policy, stable_file(root, POLICY_REL)


def _static_by_id(policy: Mapping[str, object]) -> dict[str, dict[str, object]]:
    rows = policy.get("static_inputs")
    if not isinstance(rows, list) or not rows:
        raise ContractError("V4 static input family is absent")
    by_id: dict[str, dict[str, object]] = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"input_id", "path", "bytes", "sha256"}:
            raise ContractError("V4 static input record is malformed")
        identity = str(row["input_id"])
        if identity in by_id or not SAFE_TOKEN_RE.fullmatch(identity):
            raise ContractError("V4 static input identity is duplicated or unsafe")
        safe_relative(str(row["path"]))
        if not isinstance(row["bytes"], int) or row["bytes"] <= 0 or not re.fullmatch(r"[0-9a-f]{64}", str(row["sha256"])):
            raise ContractError(f"V4 static input pin is invalid: {identity}")
        by_id[identity] = row
    expected = {
        "V3_LOCK", "FUMA_MANIFEST", "FUMA_ARCHIVE", "INTERPRETATION_POLICY",
        "MAGMA_BINARY", "MAGMA_GENE_PRODUCER", "PAIR_MANIFEST", "ANALYSIS_PANEL",
    }
    if set(by_id) != expected:
        raise ContractError("V4 static input family is incomplete")
    return by_id


def validate_static_inputs(root: Path, policy: Mapping[str, object]) -> dict[str, FileRecord]:
    records: dict[str, FileRecord] = {}
    for identity, pin in _static_by_id(policy).items():
        record = stable_file(root, str(pin["path"]))
        if record.bytes != pin["bytes"] or record.sha256 != pin["sha256"]:
            raise ContractError(f"V4 static input differs from its pin: {identity}")
        records[identity] = record
    if records["V3_LOCK"].path != PARENT_LOCK_REL or records["V3_LOCK"].sha256 != PARENT_LOCK_SHA256:
        raise ContractError("V4 static family does not retain the exact V3 lock")
    return records


def validate_parent_v3_family(root: Path) -> None:
    try:
        lock = json.loads(read_stable_bytes(root, PARENT_LOCK_REL).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError("frozen V3 parent lock is invalid JSON") from exc
    if not isinstance(lock, dict):
        raise ContractError("frozen V3 parent lock is not an object")
    if (
        lock.get("schema_version") != "track-b-mechanism-pre-result-readiness.3"
        or lock.get("contract_revision") != 3
        or lock.get("contract_kind") != "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK"
        or lock.get("authoritative_for_execution") is not True
        or lock.get("result_blind") is not True
        or lock.get("future_results_accessed") is not False
        or lock.get("policy_sha256") != "30f6cd67b1227421664f55b3cb20f5c2139f9f0293365d687c42fa03309de0a4"
        or lock.get("script_sha256") != "28569632ecbeb842f07d0214a6de7912f18d089ce7f813f6aec16f867c7c1f1b"
    ):
        raise ContractError("frozen V3 parent lock semantics differ")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or len(outputs) != 5:
        raise ContractError("frozen V3 parent output family is incomplete")
    for relative, expected in outputs.items():
        if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
            raise ContractError("frozen V3 parent output identity is malformed")
        record = stable_file(root, str(relative))
        if record.bytes != expected["bytes"] or record.sha256 != expected["sha256"]:
            raise ContractError(f"frozen V3 parent output differs: {relative}")


def _dataset_specs(policy: Mapping[str, object]) -> dict[str, dict[str, object]]:
    datasets = policy.get("source_datasets")
    if not isinstance(datasets, list) or len(datasets) != 2:
        raise ContractError("V4 source dataset family must contain exactly two human matrices")
    by_id: dict[str, dict[str, object]] = {}
    for row in datasets:
        if not isinstance(row, dict):
            raise ContractError("V4 source dataset record is malformed")
        identity = str(row.get("dataset_id"))
        if identity in by_id:
            raise ContractError("V4 source dataset identity is duplicated")
        by_id[identity] = row
    if set(by_id) != {"Allen_Human_MTG_level2", "GSE67835_Human_Cortex_woFetal"}:
        raise ContractError("V4 source family substituted an unapproved atlas")
    for identity, expected in {
        "Allen_Human_MTG_level2": (
            "human", "middle_temporal_gyrus", 75, 29155, 30505558,
            "eae19362b53e8d6eec3d2a50582fd48e8738b2c7fc49ee3b18988c2940d72124",
            13781051, "f20bdb4e9a5d759063a1dd28e7fe8658f2893042e346077f89d97ca1151a565b",
            "FUMA_scRNA_data-dd526163ea80af1a80a6cdc80db167144500694b/processed_data/Allen_Human_MTG_level2.txt.gz",
        ),
        "GSE67835_Human_Cortex_woFetal": (
            "human", "cerebral_cortex", 7, 19749, 2438789,
            "35676120adf207a8153b5a40947cfe8ca0b6fd23bb77d0865847ec94be9c2c57",
            1048690, "aaf3a2ad3e0e58a4b275c082189b6c486e5594a36c3b140084a806c16ac62780",
            "FUMA_scRNA_data-dd526163ea80af1a80a6cdc80db167144500694b/processed_data/GSE67835_Human_Cortex_woFetal.txt.gz",
        ),
    }.items():
        row = by_id[identity]
        if (
            row.get("species") != expected[0] or row.get("tissue") != expected[1]
            or row.get("cell_type_count") != expected[2] or row.get("gene_rows") != expected[3]
            or row.get("text_bytes") != expected[4] or row.get("text_sha256") != expected[5]
            or row.get("compressed_bytes") != expected[6] or row.get("compressed_sha256") != expected[7]
            or row.get("archive_member") != expected[8]
            or set(row) != {
                "dataset_id", "archive_member", "species", "tissue", "gene_annotation",
                "gene_rows", "cell_type_count", "compressed_bytes", "compressed_sha256",
                "text_bytes", "text_sha256",
            }
            or row.get("gene_annotation") != "Ensembl_v92_GRCh37_20260_coding_genes"
        ):
            raise ContractError(f"V4 human source semantics differ from their pin: {identity}")
    return by_id


class _HashingReader(io.RawIOBase):
    def __init__(self, source: object) -> None:
        self.source = source
        self.digest = hashlib.sha256()
        self.total = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: bytearray) -> int:
        data = self.source.read(len(buffer))  # type: ignore[attr-defined]
        if not data:
            return 0
        self.digest.update(data)
        self.total += len(data)
        buffer[:len(data)] = data
        return len(data)


def archive_headers(root: Path, policy: Mapping[str, object]) -> dict[str, list[str]]:
    """Deep-stream both exact source payloads and return their validated headers."""
    pins = _static_by_id(policy)
    archive_record = stable_file(root, str(pins["FUMA_ARCHIVE"]["path"]))
    if archive_record.bytes != pins["FUMA_ARCHIVE"]["bytes"] or archive_record.sha256 != pins["FUMA_ARCHIVE"]["sha256"]:
        raise ContractError("FUMA archive differs from the V4 pin")
    specs = _dataset_specs(policy)
    headers: dict[str, list[str]] = {}
    try:
        with tarfile.open(root / archive_record.path, "r:gz") as archive:
            members = {member.name: member for member in archive.getmembers()}
            for identity, spec in specs.items():
                name = str(spec["archive_member"])
                member = members.get(name)
                if member is None or not member.isfile() or member.size != spec["compressed_bytes"]:
                    raise ContractError(f"FUMA archive member identity differs: {identity}")
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise ContractError(f"cannot read FUMA archive member: {identity}")
                raw = _HashingReader(extracted)
                decoded_digest = hashlib.sha256()
                decoded_bytes = 0
                row_count = 0
                genes: set[str] = set()
                with extracted, io.BufferedReader(raw) as buffered, gzip.GzipFile(fileobj=buffered, mode="rb") as decoded:
                    raw_header = decoded.readline()
                    decoded_digest.update(raw_header)
                    decoded_bytes += len(raw_header)
                    try:
                        header = raw_header.decode("utf-8").split()
                    except UnicodeError as exc:
                        raise ContractError(f"non-UTF-8 FUMA source header: {identity}") from exc
                    expected_width = int(spec["cell_type_count"]) + 2
                    if len(header) != expected_width:
                        raise ContractError(f"wrong FUMA source header width: {identity}")
                    for line_number, payload in enumerate(decoded, 2):
                        decoded_digest.update(payload)
                        decoded_bytes += len(payload)
                        try:
                            fields = payload.decode("utf-8").split()
                        except UnicodeError as exc:
                            raise ContractError(f"non-UTF-8 FUMA source row: {identity}/{line_number}") from exc
                        if (
                            len(fields) != expected_width or not GENE_RE.fullmatch(fields[0])
                            or fields[0] in genes
                        ):
                            raise ContractError(f"malformed/duplicate FUMA source row: {identity}/{line_number}")
                        try:
                            numeric = [float(value) for value in fields[1:]]
                        except ValueError as exc:
                            raise ContractError(f"nonnumeric FUMA source row: {identity}/{line_number}") from exc
                        if any(not math.isfinite(value) or value < 0 for value in numeric):
                            raise ContractError(f"invalid expression in FUMA source row: {identity}/{line_number}")
                        genes.add(fields[0])
                        row_count += 1
                if (
                    raw.total != spec["compressed_bytes"]
                    or raw.digest.hexdigest() != spec["compressed_sha256"]
                    or decoded_bytes != spec["text_bytes"]
                    or decoded_digest.hexdigest() != spec["text_sha256"]
                    or row_count != spec["gene_rows"]
                ):
                    raise ContractError(f"deep FUMA source identity/count differs: {identity}")
                headers[identity] = header
    except (tarfile.TarError, OSError, UnicodeError) as exc:
        raise ContractError(f"cannot validate FUMA archive headers: {exc}") from exc
    allen = headers["Allen_Human_MTG_level2"]
    gse = headers["GSE67835_Human_Cortex_woFetal"]
    if allen[:1] != ["GENE"] or allen[-1:] != ["Average"] or len(allen) != 77:
        raise ContractError("Allen MTG header is not GENE + 75 subtypes + Average")
    if gse != ["GENE", *GSE_ORDER, "Average"]:
        raise ContractError("GSE67835 header is not the exact seven-class-plus-Average family")
    class_mapping_rows(allen[1:-1])
    return headers


def validate_runtime_contract(root: Path) -> dict[str, object]:
    policy, policy_record = load_policy(root)
    records = validate_static_inputs(root, policy)
    validate_parent_v3_family(root)
    executors = policy.get("executors")
    if not isinstance(executors, dict) or set(executors) != {"enrichment", "classifier"}:
        raise ContractError("V4 executor family is malformed")
    for name, relative in (("enrichment", SELF_REL), ("classifier", CLASSIFIER_REL)):
        spec = executors.get(name)
        if not isinstance(spec, dict) or set(spec) != {"path", "sha256"} or spec.get("path") != relative:
            raise ContractError(f"V4 executor pin is malformed: {name}")
        record = stable_file(root, relative)
        if record.sha256 != spec.get("sha256"):
            raise ContractError(f"V4 executor differs from its checksum pin: {name}")
    lock = read_canonical_json(root, CONTRACT_LOCK_REL)
    if (
        set(lock) != {
            "schema_version", "contract_revision", "contract_kind",
            "authoritative_for_execution", "result_blind", "future_results_accessed",
            "policy_path", "policy_sha256", "builder_path", "builder_sha256",
            "parent_v3_lock_path", "parent_v3_lock_sha256", "analysis_ids",
            "mapping_row_count", "alignment_row_count", "outputs",
        }
        or lock.get("schema_version") != LOCK_SCHEMA
        or lock.get("contract_revision") != 4
        or lock.get("contract_kind") != "PRE_RESULT_P14_P16_EXECUTION_AMENDMENT_LOCK"
        or lock.get("authoritative_for_execution") is not True
        or lock.get("policy_sha256") != policy_record.sha256
        or lock.get("policy_path") != POLICY_REL
        or lock.get("parent_v3_lock_sha256") != PARENT_LOCK_SHA256
        or lock.get("parent_v3_lock_path") != PARENT_LOCK_REL
        or lock.get("result_blind") is not True
        or lock.get("future_results_accessed") is not False
        or lock.get("analysis_ids") != [
            "P14_BRAIN_CELL_CLASS_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION",
        ]
        or lock.get("mapping_row_count") != 75
        or lock.get("alignment_row_count") != 7
    ):
        raise ContractError("V4 execution-amendment lock lineage is invalid")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != set(CONTRACT_OUTPUT_RELS):
        raise ContractError("V4 execution-amendment output family is incomplete")
    for relative, expected in outputs.items():
        if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
            raise ContractError("V4 locked output identity is malformed")
        record = stable_file(root, str(relative))
        if record.bytes != expected["bytes"] or record.sha256 != expected["sha256"]:
            raise ContractError(f"V4 locked output differs: {relative}")
    if lock.get("builder_path") != BUILDER_REL:
        raise ContractError("V4 lock names an unexpected builder")
    builder_record = stable_file(root, BUILDER_REL)
    if builder_record.sha256 != lock.get("builder_sha256"):
        raise ContractError("V4 builder differs from the lock pin")
    headers = archive_headers(root, policy)
    expected_map = tsv_bytes(MAPPING_FIELDS, class_mapping_rows(headers["Allen_Human_MTG_level2"][1:-1]))
    expected_alignment = tsv_bytes(ALIGNMENT_FIELDS, alignment_rows())
    if read_stable_bytes(root, MAPPING_REL) != expected_map or read_stable_bytes(root, ALIGNMENT_REL) != expected_alignment:
        raise ContractError("V4 mapping/alignment artifacts are not reproduced by the frozen human headers")
    binary = root / records["MAGMA_BINARY"].path
    version = subprocess.run([binary, "--version"], cwd=root, capture_output=True, text=True, check=False, timeout=30)
    if version.returncode or (version.stdout + version.stderr).strip() != "MAGMA version: v1.10 (custom)":
        raise ContractError("V4 MAGMA binary does not report the frozen v1.10 build")
    return policy


def _read_tsv_stable(
    root: Path, relative: str, *, maximum_bytes: int = 64 * 1024 * 1024,
) -> tuple[list[str], list[dict[str, str]], FileRecord]:
    record = stable_file(root, relative)
    if record.bytes > maximum_bytes:
        raise ContractError(f"TSV exceeds its RAM-bounded contract: {relative}")
    path = root / relative
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        with os.fdopen(fd, "r", encoding="utf-8", newline="", closefd=False) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            rows = list(reader)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    if _identity(before) != record.identity or _identity(after) != record.identity or _identity(final) != record.identity:
        raise ContractError(f"TSV drifted while parsing: {relative}")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ContractError(f"TSV contains a malformed row: {relative}")
    return fields, rows, record


def validate_upstream_gate(root: Path, gate_relative: str, policy: Mapping[str, object]) -> dict[str, object]:
    gate_relative = safe_relative(gate_relative)
    gate = read_canonical_json(root, gate_relative)
    required = {
        "schema_version", "gate_kind", "terminal_status", "pair_id",
        "target_trait_id", "target_trait_role", "brain_tissue_gate",
        "eligible_genetic_tiers", "selection_rule", "complete_candidate_family",
        "eligible_family", "upstream_contract", "failed_units_preserved",
    }
    if (
        set(gate) != required or gate.get("schema_version") != GATE_SCHEMA
        or gate.get("gate_kind") != "REAL_UPSTREAM_ELIGIBLE_LOCUS_GENE_FAMILY"
        or gate.get("terminal_status") != GATE_STATUS
        or gate.get("brain_tissue_gate") != "PASS"
        or gate.get("eligible_genetic_tiers") != ["G4", "G5"]
        or gate.get("selection_rule") != "PREDECLARED_G4_G5_AND_BRAIN_TISSUE_PASS_ONLY"
        or gate.get("failed_units_preserved") is not True
        or gate.get("target_trait_role") not in {"sleep", "external"}
    ):
        raise ContractError("upstream P14/P16 gate is absent, nonterminal, or semantically invalid")
    pair_id, trait_id = str(gate.get("pair_id")), str(gate.get("target_trait_id"))
    if pair_id not in {"A", "B", "CONTROL"} or not SAFE_TOKEN_RE.fullmatch(trait_id):
        raise ContractError("upstream P14/P16 pair/trait identity is unsafe")
    gate_policy = policy.get("upstream_gate")
    if not isinstance(gate_policy, dict):
        raise ContractError("V4 upstream gate policy is absent")
    maximum_bytes = int(gate_policy["maximum_family_bytes"])
    pair_fields, pairs, _pair_record = _read_tsv_stable(
        root, str(gate_policy["pair_manifest_path"]), maximum_bytes=maximum_bytes,
    )
    if not {"pair_id", "sleep_trait", "external_trait"}.issubset(pair_fields):
        raise ContractError("frozen pair manifest lacks required identities")
    pair_matches = [row for row in pairs if row.get("pair_id") == pair_id]
    if len(pair_matches) != 1 or trait_id != pair_matches[0].get(str(gate["target_trait_role"]) + "_trait"):
        raise ContractError("upstream gate target is not the exact frozen pair trait/role")
    panel_fields, panel, _panel_record = _read_tsv_stable(
        root, str(gate_policy["analysis_panel_path"]), maximum_bytes=maximum_bytes,
    )
    if not {"trait_id", "source_id", "ancestry", "source_status"}.issubset(panel_fields):
        raise ContractError("analysis panel lacks trait validation fields")
    panel_matches = [row for row in panel if row.get("trait_id") == trait_id]
    if len(panel_matches) != 1 or panel_matches[0].get("ancestry") != "EUR" or panel_matches[0].get("source_status") != "SOURCE_VERIFIED":
        raise ContractError("upstream gate trait is not a frozen source-verified EUR panel member")
    upstream = gate.get("upstream_contract")
    candidates = gate.get("complete_candidate_family")
    family = gate.get("eligible_family")
    if not isinstance(upstream, dict) or set(upstream) != {
        "path", "bytes", "sha256", "schema_version", "contract_kind", "terminal_status",
    }:
        raise ContractError("upstream gate contract identity is malformed")
    if (
        upstream.get("schema_version") != UPSTREAM_SOURCE_SCHEMA
        or upstream.get("contract_kind") != UPSTREAM_SOURCE_KIND
        or upstream.get("terminal_status") != UPSTREAM_SOURCE_STATUS
    ):
        raise ContractError("upstream source-contract semantic pin is invalid")
    for label, block in (("complete candidate family", candidates), ("eligible family", family)):
        if not isinstance(block, dict) or set(block) != {"path", "bytes", "sha256", "rows", "fields"}:
            raise ContractError(f"{label} identity is malformed")
        if block.get("fields") != ELIGIBLE_FIELDS or not isinstance(block.get("rows"), int) or int(block["rows"]) < 1:
            raise ContractError(f"{label} schema/count is invalid")
    for label, block in (
        ("upstream contract", upstream),
        ("complete candidate family", candidates),
        ("eligible family", family),
    ):
        assert isinstance(block, dict)
        relative = safe_relative(str(block["path"]))
        if not relative.startswith("results/track_b/") or relative.startswith("results/track_b/mechanism_followup/science/cell_class/"):
            raise ContractError(f"{label} path is outside the allowed upstream namespace")
        record = stable_file(root, relative)
        if record.bytes > maximum_bytes or record.bytes != block["bytes"] or record.sha256 != block["sha256"]:
            raise ContractError(f"{label} differs from the gate identity")
    source = read_canonical_json(root, str(upstream["path"]))
    if (
        set(source) != {
            "schema_version", "contract_kind", "terminal_status", "pair_id",
            "target_trait_id", "target_trait_role", "selection_rule",
            "complete_candidate_family", "eligible_family", "failed_units_preserved",
        }
        or source.get("schema_version") != UPSTREAM_SOURCE_SCHEMA
        or source.get("contract_kind") != UPSTREAM_SOURCE_KIND
        or source.get("terminal_status") != UPSTREAM_SOURCE_STATUS
        or source.get("pair_id") != pair_id
        or source.get("target_trait_id") != trait_id
        or source.get("target_trait_role") != gate.get("target_trait_role")
        or source.get("selection_rule") != gate.get("selection_rule")
        or source.get("complete_candidate_family") != candidates
        or source.get("eligible_family") != family
        or source.get("failed_units_preserved") is not True
    ):
        raise ContractError("upstream source contract does not bind the complete and eligible families")
    candidate_fields, candidate_rows, _candidate_record = _read_tsv_stable(
        root, str(candidates["path"]), maximum_bytes=maximum_bytes,
    )
    if candidate_fields != ELIGIBLE_FIELDS or len(candidate_rows) != candidates["rows"]:
        raise ContractError("upstream complete candidate family differs from its schema/count")
    identities: set[tuple[str, str]] = set()
    selected_rows: list[dict[str, str]] = []
    for row in candidate_rows:
        key = (row.get("locus_id", ""), row.get("gene_id", ""))
        if (
            row.get("pair_id") != pair_id or row.get("target_trait_id") != trait_id
            or row.get("genetic_evidence_tier") not in {"G1", "G2", "G3", "G4", "G5", "NA", "FAILED"}
            or row.get("brain_tissue_gate") not in {"PASS", "FAIL", "NA", "FAILED"}
            or row.get("eligibility_status") not in {"ELIGIBLE", "INELIGIBLE", "FAILED_UPSTREAM"}
            or not SAFE_TOKEN_RE.fullmatch(key[0])
            or not (GENE_RE.fullmatch(key[1]) or key[1] == "NA")
            or key in identities
        ):
            raise ContractError("upstream complete candidate family contains an invalid or duplicate unit")
        meets_rule = row["genetic_evidence_tier"] in {"G4", "G5"} and row["brain_tissue_gate"] == "PASS"
        if meets_rule != (row["eligibility_status"] == "ELIGIBLE"):
            raise ContractError("upstream complete candidate family does not apply the frozen selection rule exactly")
        if row["eligibility_status"] == "FAILED_UPSTREAM" and meets_rule:
            raise ContractError("an upstream failure was relabelled eligible")
        identities.add(key)
        if meets_rule:
            if not GENE_RE.fullmatch(key[1]):
                raise ContractError("eligible upstream unit lacks an Ensembl gene identity")
            selected_rows.append(row)
    if not selected_rows:
        raise ContractError("upstream complete candidate family has no G4/G5 brain-PASS unit")
    candidate_payload = tsv_bytes(ELIGIBLE_FIELDS, candidate_rows)
    if read_stable_bytes(root, str(candidates["path"]), limit=maximum_bytes) != candidate_payload:
        raise ContractError("upstream complete candidate family is not canonical TSV")
    eligible_payload = tsv_bytes(ELIGIBLE_FIELDS, selected_rows)
    eligible_fields, eligible_rows, _eligible_record = _read_tsv_stable(
        root, str(family["path"]), maximum_bytes=maximum_bytes,
    )
    if (
        eligible_fields != ELIGIBLE_FIELDS or len(eligible_rows) != family["rows"]
        or read_stable_bytes(root, str(family["path"]), limit=maximum_bytes) != eligible_payload
    ):
        raise ContractError("eligible family is not the exact deterministic subset of the complete candidate family")
    return {
        **gate, "gate_path": gate_relative,
        "gate_sha256": stable_file(root, gate_relative).sha256,
        "target_source_id": panel_matches[0]["source_id"],
    }


def _materialize_dataset(root: Path, policy: Mapping[str, object], dataset_id: str, output: Path) -> None:
    specs = _dataset_specs(policy)
    spec = specs[dataset_id]
    archive_pin = _static_by_id(policy)["FUMA_ARCHIVE"]
    archive_path = root / str(archive_pin["path"])
    compressed_path = output.with_suffix(".member.gz")
    try:
        with tarfile.open(archive_path, "r:gz") as archive:
            matches = [member for member in archive.getmembers() if member.name == spec["archive_member"]]
            if len(matches) != 1 or not matches[0].isfile() or matches[0].size != spec["compressed_bytes"]:
                raise ContractError(f"source member is absent/duplicated/wrong-size: {dataset_id}")
            source = archive.extractfile(matches[0])
            if source is None:
                raise ContractError(f"cannot extract source member: {dataset_id}")
            digest = hashlib.sha256()
            total = 0
            with source, compressed_path.open("xb") as target:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(block)
                    total += len(block)
                    target.write(block)
            if total != spec["compressed_bytes"] or digest.hexdigest() != spec["compressed_sha256"]:
                raise ContractError(f"compressed source member differs from its pin: {dataset_id}")
        digest = hashlib.sha256()
        total = 0
        with gzip.open(compressed_path, "rb") as source, output.open("xb") as target:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
                total += len(block)
                target.write(block)
        if total != spec["text_bytes"] or digest.hexdigest() != spec["text_sha256"]:
            raise ContractError(f"decoded source matrix differs from its pin: {dataset_id}")
    except (tarfile.TarError, OSError) as exc:
        raise ContractError(f"cannot materialize frozen source matrix {dataset_id}: {exc}") from exc
    finally:
        try:
            compressed_path.unlink()
        except FileNotFoundError:
            pass


def aggregate_allen_matrix(source: Path, mapping: Sequence[Mapping[str, str]], output: Path) -> dict[str, object]:
    ordered_columns = [row["source_column"] for row in mapping]
    if len(mapping) != 75 or [int(row["source_position"]) for row in mapping] != list(range(1, 76)):
        raise ContractError("frozen Allen mapping is incomplete or reordered")
    expected_rows = parse_tsv_payload(
        tsv_bytes(MAPPING_FIELDS, class_mapping_rows(ordered_columns)),
        MAPPING_FIELDS,
        "recomputed Allen class mapping",
    )
    if list(mapping) != expected_rows:
        raise ContractError("Allen mapping semantics differ from the exhaustive prefix rule")
    by_class = {
        broad: [index for index, row in enumerate(mapping) if row["broad_class"] == broad]
        for broad in CLASS_ORDER
    }
    digest = hashlib.sha256()
    gene_ids: set[str] = set()
    rows = 0
    with source.open(encoding="utf-8", newline="") as source_handle, output.open("x", encoding="utf-8", newline="") as target:
        header = source_handle.readline().split()
        if header != ["GENE", *ordered_columns, "Average"]:
            raise ContractError("Allen source matrix header differs from the frozen mapping")
        header_payload = "\t".join(["GENE", *CLASS_ORDER, "Average"]) + "\n"
        target.write(header_payload)
        digest.update(header_payload.encode())
        for line_number, line in enumerate(source_handle, 2):
            fields = line.split()
            if len(fields) != 77 or not GENE_RE.fullmatch(fields[0]) or fields[0] in gene_ids:
                raise ContractError(f"Allen matrix row is malformed or duplicated at line {line_number}")
            try:
                values = [float(value) for value in fields[1:76]]
                average = float(fields[76])
            except ValueError as exc:
                raise ContractError(f"Allen matrix row is nonnumeric at line {line_number}") from exc
            if any(not math.isfinite(value) or value < 0 for value in [*values, average]):
                raise ContractError(f"Allen matrix row contains invalid expression at line {line_number}")
            aggregates = [math.fsum(values[index] for index in by_class[broad]) / len(by_class[broad]) for broad in CLASS_ORDER]
            payload = "\t".join([fields[0], *(format(value, ".17g") for value in aggregates), fields[76]]) + "\n"
            target.write(payload)
            digest.update(payload.encode())
            gene_ids.add(fields[0])
            rows += 1
    if rows != 29155:
        raise ContractError(f"Allen matrix gene family differs from its 29,155-row pin: {rows}")
    output_record = stable_absolute_file(output)
    if output_record.sha256 != digest.hexdigest():
        raise ContractError("Allen aggregate matrix drifted while it was constructed")
    return {"gene_rows": rows, "sha256": output_record.sha256, "bytes": output_record.bytes}


def _value_or_na(value: str) -> float | None:
    if value.strip().upper() in {"NA", "NAN", ".", ""}:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def normalize_gsa(
    path: Path,
    expected: Sequence[str],
    dataset_id: str,
    trait_id: str,
    aggregation: str,
    member_counts: Mapping[str, int],
    *,
    minimum_gene_rows: int = 10000,
    forced_failure: str | None = None,
) -> list[dict[str, object]]:
    if len(expected) != 7 or len(set(expected)) != 7 or set(member_counts) != set(expected):
        raise ContractError("MAGMA normalization requires one exact seven-class family")
    header: list[str] | None = None
    observed: dict[str, dict[str, str]] = {}
    schema_failed = False
    if path.is_file():
        try:
            with path.open(encoding="utf-8") as handle:
                for line in handle:
                    values = line.split()
                    if not values:
                        continue
                    if values[0] == "VARIABLE":
                        if header is not None:
                            schema_failed = True
                            break
                        header = values
                        continue
                    if header is None or values[0].startswith("#"):
                        continue
                    if len(values) != len(header):
                        schema_failed = True
                        break
                    row = dict(zip(header, values))
                    variable = row.get("VARIABLE", "")
                    if row.get("TYPE") != "COVAR" or variable not in expected or variable in observed:
                        schema_failed = True
                        break
                    observed[variable] = row
        except (OSError, UnicodeError):
            schema_failed = True
    if header is None or not {"VARIABLE", "TYPE", "NGENES", "BETA", "SE", "P"}.issubset(header):
        schema_failed = True
    if schema_failed:
        observed = {}
    preliminary: list[dict[str, object]] = []
    p_for_bh: list[float] = []
    for variable in expected:
        raw = observed.get(variable)
        status = forced_failure or "VALID"
        n_genes: int | None = None
        beta = se = p_value = None
        if schema_failed:
            status = forced_failure or "FAILED_OUTPUT_SCHEMA"
        elif raw is None:
            status = forced_failure or "FAILED_MISSING_RESULT_ROW"
        else:
            beta, se, p_value = (_value_or_na(raw.get(field, "NA")) for field in ("BETA", "SE", "P"))
            try:
                n_genes = int(raw.get("NGENES", ""))
            except ValueError:
                n_genes = None
            if forced_failure is not None:
                status = forced_failure
            elif n_genes is None or n_genes < minimum_gene_rows:
                status = "FAILED_GENE_COUNT"
            elif beta is None or se is None or p_value is None:
                status = "FAILED_NA_OR_NONFINITE"
            elif not 0 <= p_value <= 1:
                status = "FAILED_INVALID_P"
            elif se <= 0:
                status = "FAILED_NONPOSITIVE_SE"
        bh_p = p_value if status == "VALID" and p_value is not None and 0 <= p_value <= 1 else 1.0
        p_for_bh.append(bh_p)
        z = beta / se if status == "VALID" and beta is not None and se is not None else None
        preliminary.append({
            "dataset_id": dataset_id,
            "target_trait_id": trait_id,
            "class_id": variable,
            "source_column": variable,
            "aggregation_rule": aggregation,
            "member_count": member_counts[variable],
            "analysis_status": status,
            "n_genes": "NA" if n_genes is None else n_genes,
            "beta": "NA" if beta is None else format(beta, ".17g"),
            "se": "NA" if se is None else format(se, ".17g"),
            "z": "NA" if z is None else format(z, ".17g"),
            "p_value": "NA" if p_value is None else format(p_value, ".17g"),
            "p_value_for_bh": format(bh_p, ".17g"),
            "bh_fdr": "PENDING",
            "direction_pass": "YES" if status == "VALID" and beta is not None and beta > 0 else "NO",
        })
    adjusted = bh_adjust(list(expected), p_for_bh)
    for row in preliminary:
        row["bh_fdr"] = format(adjusted[str(row["class_id"])], ".17g")
    return preliminary


def _secure_mkdirs(root: Path, relative: str) -> Path:
    current = root
    for part in PurePosixPath(safe_relative(relative)).parts:
        current = current / part
        try:
            info = os.lstat(current)
        except FileNotFoundError:
            os.mkdir(current, 0o755)
            info = os.lstat(current)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise ContractError(f"unsafe output directory component: {current}")
    return current


def publish_directory_no_replace(
    root: Path,
    out_relative: str,
    payloads: Mapping[str, bytes],
    *,
    completion_name: str,
) -> None:
    out_relative = safe_relative(out_relative)
    if not out_relative.startswith("results/track_b/mechanism_followup/science/cell_class/"):
        raise ContractError("cell-class science output escapes its dedicated namespace")
    target = root / out_relative
    try:
        os.lstat(target)
    except FileNotFoundError:
        pass
    else:
        raise ContractError("cell-class science output directory already exists; refusing replacement")
    if completion_name not in payloads:
        raise ContractError("cell-class completion artifact is absent")
    parent_rel = PurePosixPath(out_relative).parent.as_posix()
    parent = _secure_mkdirs(root, parent_rel)
    stage = Path(tempfile.mkdtemp(prefix=".cell-class-stage-", dir=parent))
    try:
        for name, payload in payloads.items():
            if PurePosixPath(name).name != name or name in {"", ".", ".."}:
                raise ContractError(f"unsafe cell-class output name: {name}")
            fd = os.open(stage / name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
            try:
                view = memoryview(payload)
                while view:
                    written = os.write(fd, view)
                    if written <= 0:
                        raise ContractError(f"short write while staging cell-class output: {name}")
                    view = view[written:]
                os.fsync(fd)
            finally:
                os.close(fd)
            os.chmod(stage / name, 0o444)
        try:
            os.mkdir(target, 0o755)
        except FileExistsError as exc:
            raise ContractError("cell-class publication race found an existing target") from exc
        ordered = [name for name in sorted(payloads) if name != completion_name] + [completion_name]
        for name in ordered:
            try:
                os.link(stage / name, target / name, follow_symlinks=False)
            except FileExistsError as exc:
                raise ContractError(f"cell-class publication target appeared concurrently: {name}") from exc
            published = stable_absolute_file(target / name)
            if published.bytes != len(payloads[name]) or published.sha256 != hashlib.sha256(payloads[name]).hexdigest():
                raise ContractError(f"cell-class publication identity mismatch: {name}")
        target_fd = os.open(target, os.O_RDONLY)
        try:
            os.fsync(target_fd)
        finally:
            os.close(target_fd)
        fd = os.open(parent, os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    finally:
        if stage.exists():
            for child in stage.iterdir():
                child.unlink()
            stage.rmdir()


def magma_gene_raw_count(root: Path, relative: str) -> int:
    record = stable_file(root, relative)
    path = root / relative
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        header_seen = False
        rows = 0
        genes: set[str] = set()
        with os.fdopen(fd, "r", encoding="utf-8", newline="", closefd=False) as handle:
            for line in handle:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "GENE":
                    if header_seen:
                        raise ContractError("MAGMA gene-results family has a duplicate header")
                    header_seen = True
                    continue
                if not header_seen or not GENE_RE.fullmatch(fields[0]) or fields[0] in genes:
                    raise ContractError("MAGMA gene-results family has a malformed row")
                genes.add(fields[0])
                rows += 1
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    if (
        _identity(before) != record.identity or _identity(after) != record.identity
        or _identity(final) != record.identity
    ):
        raise ContractError("MAGMA gene-results family drifted during validation")
    if not header_seen:
        raise ContractError("MAGMA gene-results family lacks its GENE header")
    return rows


def validate_gene_results(
    root: Path,
    trait_id: str,
    target_source_id: str,
    policy: Mapping[str, object],
) -> tuple[str, FileRecord, str, FileRecord]:
    gene_spec = policy.get("upstream_gate")
    if not isinstance(gene_spec, dict):
        raise ContractError("V4 upstream-gate policy is malformed")
    genes_rel = safe_relative(str(gene_spec["gene_results_path_template"]).format(trait_id=trait_id))
    provenance_rel = safe_relative(str(gene_spec["gene_results_provenance_template"]).format(trait_id=trait_id))
    genes_record = stable_file(root, genes_rel)
    provenance_record = stable_file(root, provenance_rel)
    try:
        provenance = json.loads(read_stable_bytes(root, provenance_rel).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError("MAGMA gene-results provenance is invalid") from exc
    try:
        interpretation = json.loads(
            read_stable_bytes(root, str(_static_by_id(policy)["INTERPRETATION_POLICY"]["path"])).decode("utf-8")
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError("frozen interpretation policy is invalid JSON") from exc
    if not isinstance(interpretation, dict):
        raise ContractError("frozen interpretation policy is not an object")
    minimum_rows = int(gene_spec["minimum_gene_rows"])
    observed_rows = magma_gene_raw_count(root, genes_rel)
    source_path = Path(str(provenance.get("genes_raw_path", "")))
    if not source_path.is_absolute():
        source_path = root / source_path
    try:
        source_path = source_path.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise ContractError("MAGMA gene-results provenance names an invalid source path") from exc
    if (
        not isinstance(provenance, dict)
        or provenance.get("schema_version") != interpretation.get("schema_version")
        or provenance.get("analysis_id") != interpretation.get("analysis_id")
        or provenance.get("trait_id") != trait_id
        or provenance.get("trait_source_id") != target_source_id
        or provenance.get("policy_sha256") != _static_by_id(policy)["INTERPRETATION_POLICY"]["sha256"]
        or provenance.get("component_manifest_sha256") != _static_by_id(policy)["FUMA_MANIFEST"]["sha256"]
        or provenance.get("magma_version") != "1.10"
        or provenance.get("magma_binary_sha256") != _static_by_id(policy)["MAGMA_BINARY"]["sha256"]
        or provenance.get("genes_raw_sha256") != genes_record.sha256
        or provenance.get("gene_rows") != observed_rows
        or observed_rows < minimum_rows
        or source_path != (root / genes_rel).resolve(strict=True)
    ):
        raise ContractError("MAGMA gene results are not the exact policy-bound full genome-wide gene family")
    return genes_rel, genes_record, provenance_rel, provenance_record


def _magma_failure_code(returncode: int, exception_name: str | None) -> str | None:
    if exception_name == "TIMEOUT":
        return "FAILED_MAGMA_TIMEOUT"
    if exception_name is not None:
        return "FAILED_MAGMA_EXECUTION"
    if returncode != 0:
        return "FAILED_MAGMA_EXIT"
    return None


def run_production(root: Path, gate_relative: str, out_relative: str) -> dict[str, object]:
    policy = validate_runtime_contract(root)
    # This gate validation is intentionally before tempfile, mkdir, or output access.
    gate = validate_upstream_gate(root, gate_relative, policy)
    trait_id = str(gate["target_trait_id"])
    pair_id = str(gate["pair_id"])
    production = policy["production_outputs"]
    assert isinstance(production, dict)
    expected_out = str(production["enrichment_dir_template"]).format(
        pair_id=pair_id, trait_id=trait_id,
    )
    if safe_relative(out_relative) != expected_out:
        raise ContractError(f"enrichment output must use the frozen pair/trait namespace: {expected_out}")
    try:
        os.lstat(root / expected_out)
    except FileNotFoundError:
        pass
    else:
        raise ContractError("enrichment output directory already exists; refusing to rerun science")
    genes_rel, genes_record, provenance_rel, provenance_record = validate_gene_results(
        root, trait_id, str(gate["target_source_id"]), policy,
    )
    gene_spec = policy["upstream_gate"]
    assert isinstance(gene_spec, dict)
    mapping_payload = read_stable_bytes(root, MAPPING_REL)
    mapping = parse_tsv_payload(mapping_payload, MAPPING_FIELDS, "Allen class mapping")
    binary = root / str(_static_by_id(policy)["MAGMA_BINARY"]["path"])
    source_specs = _dataset_specs(policy)
    with tempfile.TemporaryDirectory(prefix="track_b_p14_p16_") as temporary:
        temp = Path(temporary)
        allen_source, gse_source = temp / "allen.txt", temp / "gse.txt"
        _materialize_dataset(root, policy, "Allen_Human_MTG_level2", allen_source)
        _materialize_dataset(root, policy, "GSE67835_Human_Cortex_woFetal", gse_source)
        allen_matrix = temp / "allen_7_class.txt"
        aggregate_meta = aggregate_allen_matrix(allen_source, mapping, allen_matrix)
        run_specs = (
            ("Allen_Human_MTG_level2", allen_matrix, CLASS_ORDER, "ARITHMETIC_MEAN_ALL_MAPPED_SUBTYPES_MATH_FSUM_SOURCE_ORDER", EXPECTED_CLASS_COUNTS),
            ("GSE67835_Human_Cortex_woFetal", gse_source, GSE_ORDER, "SOURCE_CLASS_COLUMN_NO_AGGREGATION", {name: 1 for name in GSE_ORDER}),
        )
        families: dict[str, list[dict[str, object]]] = {}
        run_records: list[dict[str, object]] = []
        raw_payloads: dict[str, bytes] = {}
        for dataset_id, matrix, variables, aggregation, counts in run_specs:
            prefix = temp / dataset_id
            command = [
                str(binary), "--gene-results", str(root / genes_rel),
                "--gene-covar", str(matrix), "--model", "condition-hide=Average",
                "direction=greater", "--settings", "abbreviate=0", "--out", str(prefix),
            ]
            stdout = ""
            stderr = ""
            returncode = -1
            exception_name: str | None = None
            try:
                result = subprocess.run(
                    command, cwd=root, capture_output=True, text=True, check=False,
                    timeout=int(gene_spec["magma_timeout_sec"]),
                )
                stdout, stderr, returncode = result.stdout, result.stderr, result.returncode
            except subprocess.TimeoutExpired as exc:
                exception_name = "TIMEOUT"
                stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
                stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
                returncode = -9
            except OSError as exc:
                exception_name = type(exc).__name__
                stderr = str(exc)
                returncode = -1
            gsa, log = Path(str(prefix) + ".gsa.out"), Path(str(prefix) + ".log")
            forced_failure = _magma_failure_code(returncode, exception_name)
            if forced_failure is None and (not gsa.is_file() or not log.is_file()):
                forced_failure = "FAILED_MISSING_MAGMA_ARTIFACT"
            rows = normalize_gsa(
                gsa, variables, dataset_id, trait_id, aggregation, counts,
                minimum_gene_rows=int(gene_spec["minimum_gene_rows"]),
                forced_failure=forced_failure,
            )
            families[dataset_id] = rows
            for path, label in ((gsa, "gsa"), (log, "log")):
                payload = read_absolute_bytes(path) if path.is_file() else b""
                raw_payloads[f"{dataset_id}.{label}.out"] = payload
            raw_payloads[f"{dataset_id}.stdout.txt"] = stdout.encode("utf-8", errors="replace")
            raw_payloads[f"{dataset_id}.stderr.txt"] = stderr.encode("utf-8", errors="replace")
            matrix_record = stable_absolute_file(matrix)
            run_records.append({
                "dataset_id": dataset_id, "command": command,
                "source_text_sha256": source_specs[dataset_id]["text_sha256"],
                "matrix_bytes": matrix_record.bytes, "matrix_sha256": matrix_record.sha256,
                "result_rows": len(rows), "valid_rows": sum(row["analysis_status"] == "VALID" for row in rows),
                "returncode": returncode, "execution_exception": exception_name or "NONE",
                "failure_code": forced_failure or "NONE",
            })
        allen_payload = tsv_bytes(RESULT_FIELDS, families["Allen_Human_MTG_level2"])
        gse_payload = tsv_bytes(RESULT_FIELDS, families["GSE67835_Human_Cortex_woFetal"])
    outputs: dict[str, bytes] = {
        "P14_ALLEN_BROAD_CLASS_RESULTS.tsv": allen_payload,
        "P16_GSE67835_CLASS_RESULTS.tsv": gse_payload,
        **raw_payloads,
    }
    valid_rows = sum(
        row["analysis_status"] == "VALID"
        for rows in families.values() for row in rows
    )
    provenance_out = {
        "schema_version": "track-b-p14-p16-enrichment-run.1",
        "terminal_status": (
            "COMPLETED_ALL_14_ROWS_VALID" if valid_rows == 14
            else "COMPLETED_WITH_FAILURES_PRESERVED"
        ),
        "pair_id": pair_id, "target_trait_id": trait_id,
        "upstream_gate_path": gate["gate_path"], "upstream_gate_sha256": gate["gate_sha256"],
        "eligible_family": gate["eligible_family"], "contract_lock_path": CONTRACT_LOCK_REL,
        "contract_lock_sha256": stable_file(root, CONTRACT_LOCK_REL).sha256,
        "policy_path": POLICY_REL, "policy_sha256": stable_file(root, POLICY_REL).sha256,
        "genes_raw_path": genes_rel, "genes_raw_sha256": genes_record.sha256,
        "gene_results_provenance_path": provenance_rel,
        "gene_results_provenance_sha256": provenance_record.sha256,
        "allen_aggregation": aggregate_meta, "runs": run_records,
        "multiple_testing": {
            "allen_denominator": 7, "gse_denominator": 7,
            "average_is_hypothesis": False, "failed_or_na_p_for_bh": 1.0,
        },
        "valid_result_rows": valid_rows, "total_result_rows": 14,
        "output_sha256": {name: hashlib.sha256(payload).hexdigest() for name, payload in sorted(outputs.items())},
        "claim_limit": "Broad-class gene-property enrichment is not locus mediation, a causal-cell claim, subtype replication, or pericyte evidence.",
    }
    outputs["P14_P16_ENRICHMENT.provenance.json"] = canonical_json(provenance_out)
    if set(outputs) != set(production["enrichment_required_files"]):
        raise ContractError("enrichment output family differs from the frozen production family")
    publish_directory_no_replace(
        root, out_relative, outputs,
        completion_name="P14_P16_ENRICHMENT.provenance.json",
    )
    return provenance_out


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT_AT_IMPORT)
    parser.add_argument("--upstream-gate", required=True, help="repository-relative canonical G4/G5 brain-PASS gate")
    parser.add_argument("--out-dir", required=True, help="repository-relative no-replace science output directory")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        raise ContractError("repository root is not a directory")
    result = run_production(root, args.upstream_gate, args.out_dir)
    print(f"P14_P16_ENRICHMENT_OK pair={result['pair_id']} trait={result['target_trait_id']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ContractError as exc:
        raise SystemExit(f"ERROR: {exc}")
