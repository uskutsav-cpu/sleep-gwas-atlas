#!/usr/bin/env python3
"""Freeze or verify a result-blind Track B mechanism resource contract.

This program deliberately inspects only local resources and immutable upstream
metadata.  Fine-mapping, PLACO/pleiotropy, and future mechanism results are
forbidden inputs.  READY means resource/software readiness; the current
pre-result execution state remains NOT_APPLICABLE_UNTIL_UPSTREAM until the
predeclared upstream gates exist.
"""
from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Callable, Iterable, Mapping, Sequence

sys.dont_write_bytecode = True

SCHEMA = "track-b-mechanism-pre-result-readiness.1"
POLICY_REL = "config/track_b_mechanism_followup_policy.json"
OUTPUT_DIR_REL = "results/track_b/mechanism_followup"
OUTPUT_RELS = (
    f"{OUTPUT_DIR_REL}/PRE_RESULT_RESOURCE_EVIDENCE.tsv",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.tsv",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.json",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS_REPORT.md",
)
LOCK_REL = f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.lock.json"
ALLOWED_STATUSES = {
    "READY", "BLOCKED_BY_DATA", "BLOCKED_BY_SOFTWARE",
    "NOT_APPLICABLE_UNTIL_UPSTREAM",
}
EVIDENCE_FIELDS = (
    "evidence_id", "path", "kind", "observed_bytes", "sha256",
    "file_count", "stable_before_after", "validation_status", "detail",
)
READINESS_FIELDS = (
    "analysis_id", "phase", "analysis", "scope", "resource_status", "status",
    "blocker_types", "status_reason", "exact_unblock_condition", "required_inputs",
    "required_software", "reference", "build", "ancestry", "evidence_paths",
    "evidence_sha256", "ram_decomposition_safety", "planned_output", "claim_limit",
)
EXPECTED_ANALYSIS_IDS = (
    "P07_CONJFDR", "P09_HDL_L", "P09_RHO_HESS", "P13_MAGMA_BROAD_TISSUE",
    "P13_SLDSC_GTEX", "P13_GTEX_TISSUE_QTL", "P14_BRAIN_CELL_CLASS_DISCOVERY",
    "P15_BRAIN_CELL_SUBTYPE_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION",
    "P16_HUMAN_CELL_SUBTYPE_REPLICATION", "P17_CELL_EQTL", "P18_CELL_SQTL",
    "P19_THREE_WAY_COLOC", "P20_CATLAS_SCATAC", "P20_SCREEN_CHROMATIN",
    "P20_HOCOMOCO_MOTIF", "P21_ABC_ENHANCER_GENE", "P21_PCHIC_ENHANCER_GENE",
    "P21_TWAS", "P22_REGULATORY_CHAIN", "P23_SPATIAL",
    "P24_LOCKED_ATLAS_CROSS_SLEEP", "P25_PAIR_A_VS_B", "P26_POSITIVE_CONTROL",
    "P27_PATHWAYS", "P28_GRN", "P29_BIDIRECTIONAL_MR",
)
ANALYSIS_REQUIRED_KEYS = {
    "analysis_id", "phase", "analysis", "scope", "requirements", "required_inputs",
    "required_software", "reference", "build", "ancestry", "ram_decomposition",
    "output", "unblock", "claim_limit",
}
FORBIDDEN_RESULT_PREFIXES = (
    "results/track_b/finemap", "results/track_b/fine_mapping",
    "results/track_b/placo", "results/track_b/pleiotropy",
    "results/track_b/mechanism_followup/science",
)
SHA256_RE = re.compile(r"[0-9a-f]{64}")


class ContractError(RuntimeError):
    """A fail-closed readiness-contract violation."""


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")


def bytes_sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def safe_relative(value: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value:
        raise ContractError(f"unsafe repository-relative path: {value!r}")
    pure = PurePosixPath(value)
    if pure.is_absolute() or any(part in {"", ".", ".."} for part in pure.parts):
        raise ContractError(f"unsafe repository-relative path: {value!r}")
    normalized = pure.as_posix()
    if normalized != value.rstrip("/"):
        raise ContractError(f"non-canonical repository-relative path: {value!r}")
    for forbidden in FORBIDDEN_RESULT_PREFIXES:
        if normalized == forbidden or normalized.startswith(forbidden + "/"):
            raise ContractError(f"future science/result path is forbidden: {normalized}")
    return normalized


def _check_components(root: Path, relative: str, *, allow_missing_leaf: bool = False) -> Path:
    relative = safe_relative(relative)
    current = root
    parts = PurePosixPath(relative).parts
    for index, part in enumerate(parts):
        current = current / part
        try:
            info = os.lstat(current)
        except FileNotFoundError:
            if allow_missing_leaf and index == len(parts) - 1:
                return current
            raise
        if stat.S_ISLNK(info.st_mode):
            raise ContractError(f"symlink is forbidden in evidence path: {relative}")
        if index < len(parts) - 1 and not stat.S_ISDIR(info.st_mode):
            raise ContractError(f"non-directory path component in evidence path: {relative}")
    return current


def _identity(info: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        info.st_dev, info.st_ino, info.st_mode, info.st_size,
        info.st_mtime_ns, info.st_ctime_ns,
    )


@dataclass(frozen=True)
class StableObject:
    kind: str
    observed_bytes: int
    sha256: str
    file_count: int
    token: tuple[object, ...]


def stable_file(
    root: Path,
    relative: str,
    *,
    chunk_hook: Callable[[int], None] | None = None,
) -> StableObject:
    """Hash a real regular file and prove its identity did not change while read."""
    path = _check_components(root, relative)
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode):
            raise ContractError(f"evidence is not a regular file: {relative}")
        digest = hashlib.sha256()
        block_count = 0
        while True:
            block = os.read(fd, 8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
            block_count += 1
            if chunk_hook is not None:
                chunk_hook(block_count)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    if not stat.S_ISREG(final.st_mode) or _identity(before) != _identity(after) or _identity(after) != _identity(final):
        raise ContractError(f"evidence drifted while hashing: {relative}")
    return StableObject("FILE", before.st_size, digest.hexdigest(), 1, _identity(before))


def stable_tree(root: Path, relative: str) -> StableObject:
    """Hash a directory tree without following links, streaming every file."""
    directory = _check_components(root, relative)
    info_before = os.lstat(directory)
    if not stat.S_ISDIR(info_before.st_mode):
        raise ContractError(f"tree probe is not a directory: {relative}")
    members: list[tuple[str, StableObject]] = []
    for base, dirs, files in os.walk(directory, topdown=True, followlinks=False):
        dirs.sort()
        files.sort()
        for name in dirs:
            child = Path(base) / name
            if stat.S_ISLNK(os.lstat(child).st_mode):
                raise ContractError(f"symlink is forbidden in evidence tree: {child}")
        for name in files:
            child = Path(base) / name
            rel = child.relative_to(root).as_posix()
            members.append((rel, stable_file(root, rel)))
    info_after = os.lstat(directory)
    if _identity(info_before) != _identity(info_after):
        raise ContractError(f"evidence tree drifted while hashing: {relative}")
    digest = hashlib.sha256()
    total = 0
    token_members: list[tuple[object, ...]] = []
    for rel, item in members:
        digest.update(f"{rel}\t{item.observed_bytes}\t{item.sha256}\n".encode())
        total += item.observed_bytes
        token_members.append((rel, *item.token, item.sha256))
    token = (_identity(info_before), tuple(token_members))
    return StableObject("TREE", total, digest.hexdigest(), len(members), token)


def stable_object(root: Path, relative: str) -> StableObject | None:
    try:
        path = _check_components(root, relative, allow_missing_leaf=True)
    except FileNotFoundError:
        return None
    if not path.exists():
        return None
    mode = os.lstat(path).st_mode
    if stat.S_ISREG(mode):
        return stable_file(root, relative)
    if stat.S_ISDIR(mode):
        return stable_tree(root, relative)
    if stat.S_ISLNK(mode):
        raise ContractError(f"symlink is forbidden in evidence path: {relative}")
    raise ContractError(f"unsupported filesystem object in evidence path: {relative}")


def read_stable_bytes(root: Path, relative: str) -> bytes:
    item = stable_file(root, relative)
    path = _check_components(root, relative)
    flags = os.O_RDONLY | (getattr(os, "O_NOFOLLOW", 0))
    fd = os.open(path, flags)
    try:
        before = os.fstat(fd)
        payload = bytearray()
        while True:
            block = os.read(fd, 1024 * 1024)
            if not block:
                break
            payload.extend(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    if _identity(before) != _identity(after) or len(payload) != item.observed_bytes or bytes_sha256(payload) != item.sha256:
        raise ContractError(f"file drifted while reading: {relative}")
    return bytes(payload)


def read_json(root: Path, relative: str) -> dict[str, object]:
    try:
        value = json.loads(read_stable_bytes(root, relative).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"invalid JSON in {relative}: {exc}") from exc
    if not isinstance(value, dict):
        raise ContractError(f"JSON root must be an object: {relative}")
    return value


def read_tsv(root: Path, relative: str) -> tuple[list[str], list[dict[str, str]]]:
    try:
        payload = read_stable_bytes(root, relative).decode("utf-8")
        reader = csv.DictReader(io.StringIO(payload), delimiter="\t")
        rows = list(reader)
    except UnicodeError as exc:
        raise ContractError(f"invalid UTF-8 TSV in {relative}: {exc}") from exc
    return list(reader.fieldnames or []), rows


@dataclass
class Evidence:
    evidence_id: str
    path: str
    kind: str
    observed_bytes: int | str
    sha256: str
    file_count: int | str
    stable_before_after: str
    validation_status: str
    detail: str
    token: tuple[object, ...] | None

    def row(self) -> dict[str, object]:
        return {field: getattr(self, field) for field in EVIDENCE_FIELDS}


class EvidenceRegistry:
    def __init__(self, root: Path):
        self.root = root
        self._by_path: dict[str, Evidence] = {}

    def add(
        self,
        relative: str,
        *,
        expected_bytes: int | None = None,
        expected_sha256: str | None = None,
        expected_kind: str | None = None,
        expected_file_count: int | None = None,
        required: bool = True,
        detail: str = "checksum-and-identity validated",
    ) -> str:
        relative = safe_relative(relative)
        evidence_id = "EVIDENCE::" + relative
        if relative in self._by_path:
            existing = self._by_path[relative]
            if expected_bytes is not None and existing.observed_bytes != expected_bytes:
                raise ContractError(f"conflicting size pin for {relative}")
            if expected_sha256 is not None and existing.sha256 != expected_sha256:
                raise ContractError(f"conflicting SHA-256 pin for {relative}")
            if expected_kind is not None and existing.kind != expected_kind:
                raise ContractError(f"conflicting object-kind pin for {relative}")
            if expected_file_count is not None and existing.file_count != expected_file_count:
                raise ContractError(f"conflicting file-count pin for {relative}")
            return existing.evidence_id
        item = stable_object(self.root, relative)
        if item is None:
            if required:
                raise ContractError(f"required evidence is missing: {relative}")
            evidence = Evidence(evidence_id, relative, "MISSING", "NA", "NA", "NA", "YES", "MISSING", detail, None)
        else:
            if expected_bytes is not None and item.observed_bytes != expected_bytes:
                raise ContractError(f"byte-size mismatch for {relative}: {item.observed_bytes} != {expected_bytes}")
            if expected_kind is not None and item.kind != expected_kind:
                raise ContractError(f"object-kind mismatch for {relative}: {item.kind} != {expected_kind}")
            if expected_file_count is not None and item.file_count != expected_file_count:
                raise ContractError(f"file-count mismatch for {relative}: {item.file_count} != {expected_file_count}")
            if expected_sha256 is not None:
                if not isinstance(expected_sha256, str) or not SHA256_RE.fullmatch(expected_sha256):
                    raise ContractError(f"invalid expected SHA-256 for {relative}")
                if item.sha256 != expected_sha256:
                    raise ContractError(f"SHA-256 mismatch for {relative}: {item.sha256} != {expected_sha256}")
            evidence = Evidence(
                evidence_id, relative, item.kind, item.observed_bytes, item.sha256,
                item.file_count, "YES", "VALIDATED" if required else "PRESENT_UNVALIDATED",
                detail, item.token,
            )
        self._by_path[relative] = evidence
        return evidence_id

    def get(self, evidence_id: str) -> Evidence:
        for value in self._by_path.values():
            if value.evidence_id == evidence_id:
                return value
        raise ContractError(f"unknown evidence id: {evidence_id}")

    def id_for(self, relative: str) -> str:
        return self._by_path[safe_relative(relative)].evidence_id

    def rows(self) -> list[dict[str, object]]:
        return [self._by_path[path].row() for path in sorted(self._by_path)]

    def recheck(self) -> None:
        for relative in sorted(self._by_path):
            prior = self._by_path[relative]
            current = stable_object(self.root, relative)
            if prior.kind == "MISSING":
                if current is not None:
                    raise ContractError(f"previously absent evidence appeared during validation: {relative}")
                continue
            if current is None:
                raise ContractError(f"evidence disappeared during validation: {relative}")
            if (
                current.kind != prior.kind or current.observed_bytes != prior.observed_bytes
                or current.sha256 != prior.sha256 or current.file_count != prior.file_count
                or current.token != prior.token
            ):
                raise ContractError(f"evidence drifted between pre/post snapshots: {relative}")


@dataclass(frozen=True)
class Fact:
    name: str
    ready: bool
    blocker_type: str
    reason: str
    evidence_ids: tuple[str, ...]

    def record(self) -> dict[str, object]:
        return {
            "fact_id": self.name,
            "ready": self.ready,
            "blocker_type": self.blocker_type,
            "reason": self.reason,
            "evidence_ids": list(self.evidence_ids),
        }


FACT_TYPES = {
    "abc_bundle": "DATA",
    "analysis_panel": "DATA",
    "brain_discovery_human": "DATA",
    "brain_replication_human_class": "DATA",
    "broad_tissue_magma": "DATA",
    "catlas_bundle": "DATA",
    "causal_bundle": "SOFTWARE",
    "conjfdr_runtime": "SOFTWARE",
    "cross_sleep_executor": "SOFTWARE",
    "dense_sleep_family": "DATA",
    "fuma_bundle": "DATA",
    "grn_executor": "SOFTWARE",
    "hdl_l_reference": "DATA",
    "hocomoco_bundle": "DATA",
    "hocomoco_sequence": "DATA",
    "human_subtype_replication": "DATA",
    "ldsc_source_bundle": "DATA",
    "ldsc_runtime": "SOFTWARE",
    "ldsc_static_reference": "DATA",
    "matched_multimodal_grn": "DATA",
    "pair_manifest": "DATA",
    "pathway_bundle": "DATA",
    "pathway_executor": "SOFTWARE",
    "pchic_bundle": "DATA",
    "qtl_metadata": "DATA",
    "qtl_payloads": "DATA",
    "regulatory_chain_executor": "SOFTWARE",
    "rho_hess_reference": "DATA",
    "screen_bundle": "DATA",
    "spatial": "DATA",
    "tabix": "SOFTWARE",
    "three_way_coloc_executor": "SOFTWARE",
    "track_b_mr_executor": "SOFTWARE",
    "twas_models": "DATA",
    "twas_runtime": "SOFTWARE",
}


def validate_policy_shape(policy: Mapping[str, object]) -> None:
    if policy.get("schema_version") != SCHEMA:
        raise ContractError("unexpected mechanism-readiness policy schema")
    if policy.get("analysis_build") != "GRCh37" or policy.get("analysis_ancestry") != "European":
        raise ContractError("policy analysis build/ancestry must remain GRCh37/European")
    forbidden = policy.get("forbidden_inputs")
    if forbidden != list(FORBIDDEN_RESULT_PREFIXES):
        raise ContractError("policy forbidden-input family differs from the hard-coded result-blind boundary")
    expected = policy.get("expected_analysis_ids")
    if expected != list(EXPECTED_ANALYSIS_IDS):
        raise ContractError("policy expected analysis family differs from the hard-coded complete family")
    analyses = policy.get("analyses")
    if not isinstance(analyses, list):
        raise ContractError("policy analyses must be a list")
    identifiers = [row.get("analysis_id") for row in analyses if isinstance(row, dict)]
    if identifiers != list(EXPECTED_ANALYSIS_IDS) or len(set(identifiers)) != len(identifiers):
        raise ContractError("policy analysis rows are incomplete, duplicated, or out of order")
    for row in analyses:
        if not isinstance(row, dict) or set(row) != ANALYSIS_REQUIRED_KEYS:
            raise ContractError(f"analysis row has an unexpected schema: {getattr(row, 'get', lambda *_: 'UNKNOWN')('analysis_id')}")
        if "status" in row or "resource_status" in row:
            raise ContractError("policy may not assert its own readiness status")
        requirements = row["requirements"]
        if not isinstance(requirements, list) or not requirements or any(value not in FACT_TYPES for value in requirements):
            raise ContractError(f"unknown or empty requirement family: {row['analysis_id']}")
        if len(requirements) != len(set(requirements)):
            raise ContractError(f"duplicate analysis requirement: {row['analysis_id']}")
        output = row["output"]
        if not isinstance(output, str) or not output.startswith("results/track_b/mechanism_followup/science/"):
            raise ContractError(f"planned science output escapes its dedicated future namespace: {row['analysis_id']}")
        for field in ANALYSIS_REQUIRED_KEYS - {"requirements"}:
            if not isinstance(row[field], str) or not row[field].strip():
                raise ContractError(f"empty/non-string analysis field {field}: {row['analysis_id']}")
    pinned = policy.get("pinned_inputs")
    if not isinstance(pinned, dict) or not pinned:
        raise ContractError("pinned input family is absent")
    for relative, digest in pinned.items():
        safe_relative(relative)
        if not isinstance(digest, str) or not SHA256_RE.fullmatch(digest):
            raise ContractError(f"invalid policy SHA-256 pin: {relative}")
    probes = policy.get("absence_probes")
    if not isinstance(probes, dict) or set(probes) != {
        "broad_tissue_magma", "human_subtype_replication", "qtl_payload_root",
        "psychencode_payload_root", "tabix", "hocomoco_sequence",
        "ldsc_static_reference", "ldsc_static_reference_root", "spatial",
        "matched_multimodal_grn", "track_b_mr_executor", "conjfdr_runtime",
        "hdl_l_reference", "rho_hess_reference", "three_way_coloc_executor",
        "regulatory_chain_executor", "cross_sleep_executor", "grn_executor",
    }:
        raise ContractError("absence-probe family is incomplete or has an unknown member")
    for relative in probes.values():
        safe_relative(str(relative))
    runtime_evidence = policy.get("runtime_evidence")
    if not isinstance(runtime_evidence, list) or len(runtime_evidence) != 4:
        raise ContractError("runtime evidence family must contain the exact four pinned records")
    expected_runtime_ids = {
        "LDSC_PYTHON_SITE_PACKAGES", "CAUSAL_R_LIBRARY", "BASE_R_LIBRARY", "RSCRIPT_ENTRYPOINT",
    }
    if {row.get("runtime_id") for row in runtime_evidence if isinstance(row, dict)} != expected_runtime_ids:
        raise ContractError("runtime evidence identities are incomplete")
    for row in runtime_evidence:
        if not isinstance(row, dict) or set(row) != {"runtime_id", "path", "kind", "bytes", "file_count", "sha256"}:
            raise ContractError("runtime evidence row schema is invalid")
        safe_relative(str(row["path"]))
        if row["kind"] not in {"FILE", "TREE"} or not isinstance(row["bytes"], int) or row["bytes"] <= 0 or not isinstance(row["file_count"], int) or row["file_count"] <= 0 or not isinstance(row["sha256"], str) or not SHA256_RE.fullmatch(row["sha256"]):
            raise ContractError(f"runtime evidence pin is invalid: {row['runtime_id']}")


def _add_manifest_components(registry: EvidenceRegistry, root: Path, manifest_rel: str) -> list[str]:
    manifest = read_json(root, manifest_rel)
    ids = [registry.id_for(manifest_rel)]
    components = manifest.get("components", [])
    if not isinstance(components, list):
        raise ContractError(f"manifest components are not a list: {manifest_rel}")
    for component in components:
        if not isinstance(component, dict):
            raise ContractError(f"malformed component in {manifest_rel}")
        ids.append(registry.add(
            str(component.get("path", "")),
            expected_bytes=int(component["bytes"]),
            expected_sha256=str(component["sha256"]),
            detail=f"component {component.get('component_id', 'UNKNOWN')} pinned by {manifest_rel}",
        ))
    return ids


def collect_evidence(
    root: Path, policy: Mapping[str, object]
) -> tuple[EvidenceRegistry, dict[str, dict[str, object]]]:
    registry = EvidenceRegistry(root)
    registry.add(POLICY_REL, detail="mechanism readiness policy")
    registry.add("scripts/160_build_track_b_mechanism_readiness.py", detail="mechanism readiness contract implementation")
    pinned = policy["pinned_inputs"]
    assert isinstance(pinned, dict)
    for relative in sorted(pinned):
        registry.add(relative, expected_sha256=str(pinned[relative]), detail="immutable contract/upstream pin")
    runtime_evidence = policy.get("runtime_evidence")
    if not isinstance(runtime_evidence, list) or not runtime_evidence:
        raise ContractError("runtime evidence family is absent")
    for row in runtime_evidence:
        if not isinstance(row, dict) or set(row) != {"runtime_id", "path", "kind", "bytes", "file_count", "sha256"}:
            raise ContractError("runtime evidence record has an unexpected schema")
        registry.add(
            str(row["path"]), expected_kind=str(row["kind"]),
            expected_bytes=int(row["bytes"]), expected_file_count=int(row["file_count"]),
            expected_sha256=str(row["sha256"]), detail=f"complete runtime identity: {row['runtime_id']}",
        )

    interpretation = read_json(root, "config/interpretation_analysis_policy.json")
    molecular = read_json(root, "config/molecular_analysis_policy.json")
    manifests: dict[str, dict[str, object]] = {
        "interpretation": interpretation,
        "molecular": molecular,
    }

    component_manifest_rels = (
        "config/interpretation_screen_registry_v4.json",
        "config/interpretation_fuma_scrna.json",
        "config/interpretation_catlas_adult_v4.json",
        "config/interpretation_abc_2021.json",
        "config/interpretation_pchic_2016.json",
        "config/interpretation_hocomoco_v14.json",
        "config/interpretation_causal_runtime.json",
    )
    for relative in component_manifest_rels:
        manifests[Path(relative).stem] = read_json(root, relative)
        _add_manifest_components(registry, root, relative)

    pathway_rel = "config/interpretation_pathway_sources.json"
    pathway = read_json(root, pathway_rel)
    manifests["pathway"] = pathway
    mapping = pathway.get("identifier_mapping")
    if not isinstance(mapping, dict):
        raise ContractError("pathway identifier mapping is malformed")
    registry.add(str(mapping["path"]), expected_bytes=int(mapping["bytes"]), expected_sha256=str(mapping["sha256"]), detail="pathway identifier reference")
    resources = pathway.get("resources")
    if not isinstance(resources, dict):
        raise ContractError("pathway resource family is malformed")
    for source_id in sorted(resources):
        row = resources[source_id]
        if not isinstance(row, dict):
            raise ContractError(f"malformed pathway resource: {source_id}")
        candidates = [row]
        candidates.extend(value for value in row.values() if isinstance(value, dict))
        for candidate in candidates:
            if all(key in candidate for key in ("path", "bytes", "sha256")):
                registry.add(str(candidate["path"]), expected_bytes=int(candidate["bytes"]), expected_sha256=str(candidate["sha256"]), detail=f"pathway resource {source_id}")

    ldsc_source_rel = str(interpretation["ldsc_seg_gtex"]["source_manifest"])
    registry.add(
        ldsc_source_rel,
        expected_bytes=int(interpretation["ldsc_seg_gtex"]["source_manifest_bytes"]),
        expected_sha256=str(interpretation["ldsc_seg_gtex"]["source_manifest_sha256"]),
        detail="complete S-LDSC GTEx source manifest",
    )
    ldsc_source = read_json(root, ldsc_source_rel)
    manifests["ldsc_source"] = ldsc_source
    files = ldsc_source.get("files")
    if not isinstance(files, list) or len(files) != 375:
        raise ContractError("S-LDSC source manifest is not the exact 375-file family")
    for row in files:
        if not isinstance(row, dict):
            raise ContractError("malformed S-LDSC source file record")
        registry.add(str(row["local_path"]), expected_bytes=int(row["bytes"]), expected_sha256=str(row["sha256"]), detail="S-LDSC source annotation family")

    causal = manifests["interpretation_causal_runtime"]
    runtime = causal.get("runtime")
    if not isinstance(runtime, dict):
        raise ContractError("causal runtime block is malformed")
    for key, hash_key in (
        ("genomicsem_environment_manifest", "genomicsem_environment_manifest_sha256"),
        ("dependency_source_manifest", "dependency_source_manifest_sha256"),
        ("installed_package_manifest", "installed_package_manifest_sha256"),
    ):
        registry.add(str(runtime[key]), expected_sha256=str(runtime[hash_key]), detail=f"causal runtime {key}")
    _, dependency_rows = read_tsv(root, str(runtime["dependency_source_manifest"]))
    for row in dependency_rows:
        registry.add(row["archive"], expected_bytes=int(row["bytes"]), expected_sha256=row["sha256"], detail=f"causal R source archive {row['package']}")
    ld_reference = causal.get("ld_reference")
    if not isinstance(ld_reference, dict):
        raise ContractError("causal LD reference block is malformed")
    registry.add(str(ld_reference["source_manifest"]), expected_sha256=str(ld_reference["source_manifest_sha256"]), detail="causal LD source manifest")
    registry.add(str(ld_reference["archive_path"]), expected_bytes=int(ld_reference["archive_bytes"]), expected_sha256=str(ld_reference["archive_sha256"]), detail="causal EUR LD archive")
    registry.add(str(ld_reference["instrument_universe_path"]), expected_bytes=int(ld_reference["instrument_universe_bytes"]), expected_sha256=str(ld_reference["instrument_universe_sha256"]), detail="causal HapMap3 instrument universe")

    for asset in molecular.get("metadata_assets", []):
        if not isinstance(asset, dict):
            raise ContractError("malformed molecular metadata asset")
        registry.add(str(asset["path"]), expected_bytes=int(asset["bytes"]), expected_sha256=str(asset["sha256"]), detail=f"molecular metadata {asset['id']}")
    for block_name in ("reference_build", "source_to_analysis_chain"):
        block = molecular.get(block_name)
        if not isinstance(block, dict):
            raise ContractError(f"malformed molecular {block_name}")
        registry.add(str(block["path"] if "path" in block else block["chain_path"]), expected_bytes=int(block["bytes"] if "bytes" in block else block["chain_bytes"]), expected_sha256=str(block["sha256"] if "sha256" in block else block["chain_sha256"]), detail=f"molecular {block_name}")
    generic = molecular.get("generic_engine")
    if not isinstance(generic, dict):
        raise ContractError("molecular generic engine is malformed")
    registry.add(str(generic["path"]), expected_sha256=str(generic["sha256"]), detail="molecular SuSiE/coloc engine")
    twas = molecular.get("twas")
    if not isinstance(twas, dict):
        raise ContractError("TWAS policy block is malformed")
    registry.add(str(twas["source_archive_path"]), expected_bytes=int(twas["source_archive_bytes"]), expected_sha256=str(twas["source_archive_sha256"]), detail="MetaXcan source archive")
    registry.add(str(twas["entrypoint_path"]), expected_sha256=str(twas["entrypoint_sha256"]), detail="MetaXcan entrypoint")

    _, inventory = read_tsv(root, "results/tables/twas_phi_model_inventory.tsv")
    model_root = str(twas["phi_model_source"]["install_dir"])
    registry.add(model_root, required=False, detail="complete TWAS model install-root probe")
    for row in inventory:
        registry.add(f"{model_root}/{row['filename']}", required=False, detail=f"locked TWAS model inventory file {row['file_id']}")

    probes = policy["absence_probes"]
    assert isinstance(probes, dict)
    for probe_id in sorted(probes):
        registry.add(str(probes[probe_id]), required=False, detail=f"fail-closed readiness probe: {probe_id}")
    qtl_root = str(probes["qtl_payload_root"])
    required_qtl = policy.get("required_qtl_datasets")
    if not isinstance(required_qtl, dict):
        raise ContractError("required QTL dataset family is malformed")
    for modality in sorted(required_qtl):
        datasets = required_qtl[modality]
        if not isinstance(datasets, list):
            raise ContractError("required QTL dataset list is malformed")
        for dataset_id in datasets:
            registry.add(f"{qtl_root}/{dataset_id}.tsv.gz", required=False, detail=f"required {modality} payload")
            registry.add(f"{qtl_root}/{dataset_id}.tsv.gz.tbi", required=False, detail=f"required {modality} tabix index")

    return registry, manifests


def validate_fuma_species(manifest: Mapping[str, object]) -> dict[str, Mapping[str, object]]:
    datasets = manifest.get("datasets")
    if not isinstance(datasets, list):
        raise ContractError("FUMA dataset family is malformed")
    by_id: dict[str, Mapping[str, object]] = {}
    for row in datasets:
        if not isinstance(row, dict) or not isinstance(row.get("dataset_id"), str):
            raise ContractError("FUMA dataset record is malformed")
        identity = str(row["dataset_id"])
        if identity in by_id:
            raise ContractError(f"duplicate FUMA dataset: {identity}")
        by_id[identity] = row
    expected = {
        "Allen_Human_MTG_level2": ("brain", "human", "middle_temporal_gyrus", 75),
        "GSE67835_Human_Cortex_woFetal": ("brain", "human", "cerebral_cortex", 7),
        "GSE89232_Human_Blood": ("immune", "human", "blood_dendritic_cell_compartment", 4),
        "PBMC_10x_68k": ("immune", "human", "peripheral_blood_mononuclear_cells", 11),
        "GSE81547_Human_Pancreas": ("metabolic", "human", "pancreas", 7),
        "GSE84133_Human_Pancreas": ("metabolic", "human", "pancreas", 14),
        "GSE98816_Mouse_Brain_Vascular": ("vascular", "mouse_mapped_to_human_Ensembl_gene_ID", "brain_vasculature", 15),
        "GSE99235_Mouse_Lung_Vascular": ("vascular", "mouse_mapped_to_human_Ensembl_gene_ID", "lung_vasculature", 17),
    }
    if set(by_id) != set(expected):
        raise ContractError("FUMA dataset family is not the exact eight-dataset pre-result family")
    for identity, values in expected.items():
        observed = by_id[identity]
        fields = (observed.get("domain"), observed.get("species"), observed.get("tissue"), observed.get("cell_type_count"))
        if fields != values:
            raise ContractError(f"FUMA species/context/count mismatch: {identity}: {fields!r} != {values!r}")
    if by_id["GSE67835_Human_Cortex_woFetal"]["species"] != "human":
        raise ContractError("independent human class replication is not human")
    if any(by_id[key]["species"] == "human" for key in (
        "GSE98816_Mouse_Brain_Vascular", "GSE99235_Mouse_Lung_Vascular"
    )):
        raise ContractError("mouse vascular matrices were mislabeled as human")
    return by_id


def validate_build_semantics(manifests: Mapping[str, Mapping[str, object]]) -> None:
    molecular = manifests["molecular"]
    reference_build = molecular.get("reference_build")
    if not isinstance(reference_build, dict) or reference_build.get("analysis_build") != "GRCh37" or reference_build.get("molecular_source_build") != "GRCh38":
        raise ContractError("molecular source/analysis build contract is not GRCh38-to-GRCh37")
    expected = {
        "interpretation_screen_registry_v4": "GRCh38",
        "interpretation_catlas_adult_v4": "GRCh38",
        "interpretation_abc_2021": "GRCh37",
        "interpretation_pchic_2016": "GRCh37",
    }
    for key, build in expected.items():
        if manifests[key].get("genome_build") != build:
            raise ContractError(f"wrong genome build for {key}: expected {build}")
    catlas = manifests["interpretation_catlas_adult_v4"]
    if catlas.get("life_stage") != "Adult" or catlas.get("source_counts", {}).get("adult_nuclei") != 615998:
        raise ContractError("CATlas is not the frozen human adult 615,998-nucleus atlas")
    causal = manifests["interpretation_causal_runtime"]
    ld_reference = causal.get("ld_reference")
    if not isinstance(ld_reference, dict) or ld_reference.get("build") != "GRCh37" or ld_reference.get("ancestry") != "European" or ld_reference.get("sample_count") != 503:
        raise ContractError("causal LD reference is not the frozen 503-sample European GRCh37 panel")


def validate_pair_and_panel(
    root: Path, policy: Mapping[str, object]
) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    panel_fields, panel = read_tsv(root, "config/analysis_panel.tsv")
    required_panel_fields = {
        "panel_version", "trait_id", "domain", "raw_file", "build",
        "source_status", "ancestry", "n_total",
    }
    if not required_panel_fields.issubset(panel_fields) or len(panel) != 45:
        raise ContractError("analysis panel schema/count differs from the 45-trait lock")
    if len({row["trait_id"] for row in panel}) != 45 or any(row["panel_version"] != "atlas-v1.0" for row in panel):
        raise ContractError("analysis panel IDs/version are not unique and frozen")
    sleep = [row for row in panel if row["domain"] == "sleep"]
    expected_sleep = policy.get("expected_sleep_trait_ids")
    if [row["trait_id"] for row in sleep] != expected_sleep or len(sleep) != 12:
        raise ContractError("sleep-trait family/order differs from the exact 12-trait lock")
    if any(row["ancestry"] != "EUR" or row["source_status"] != "SOURCE_VERIFIED" for row in sleep):
        raise ContractError("sleep-trait ancestry/source status differs from the lock")
    builds = {row["trait_id"]: row["build"] for row in sleep}
    if builds.get("sleep_apnea") != "hg38" or any(value != "hg19" for key, value in builds.items() if key != "sleep_apnea"):
        raise ContractError("sleep-trait mixed-build family differs from the frozen 11 hg19 + sleep_apnea hg38 design")

    pair_fields, pair_rows = read_tsv(root, "results/track_b/pair_manifest.tsv")
    required_pair_fields = {
        "pair_id", "sleep_trait", "external_trait", "ancestry", "build",
        "dense_data_readiness", "primary_scientific_role",
    }
    if not required_pair_fields.issubset(pair_fields) or len(pair_rows) != 3:
        raise ContractError("pair manifest schema/count differs from the exact three-pair lock")
    by_pair = {row["pair_id"]: row for row in pair_rows}
    expected_pairs = policy.get("expected_pair_identities")
    if not isinstance(expected_pairs, dict) or set(by_pair) != set(expected_pairs):
        raise ContractError("pair manifest identities are incomplete")
    for pair_id, traits in expected_pairs.items():
        if [by_pair[pair_id]["sleep_trait"], by_pair[pair_id]["external_trait"]] != traits:
            raise ContractError(f"frozen pair identity mismatch: {pair_id}")
        if by_pair[pair_id]["ancestry"] != "sleep=EUR;external=EUR" or by_pair[pair_id]["build"] != "sleep=hg19;external=hg19":
            raise ContractError(f"pair ancestry/build mismatch: {pair_id}")
        if by_pair[pair_id]["dense_data_readiness"] != "READY_FULL_SUMSTATS_BOTH":
            raise ContractError(f"pair dense-data readiness is not locked READY: {pair_id}")
    return by_pair, sleep


def _qc_rows_out(root: Path, trait_id: str) -> tuple[int, str]:
    relative = f"data/harmonized/{trait_id}.qc.txt"
    payload = read_stable_bytes(root, relative).decode("utf-8")
    rows_out: int | None = None
    output_build: str | None = None
    for line in payload.splitlines():
        parts = line.split("\t")
        if len(parts) == 2 and parts[0] == "rows_out":
            rows_out = int(parts[1])
        if len(parts) == 2 and parts[0] == "output_build":
            output_build = parts[1]
    if rows_out is None or rows_out <= 0 or output_build not in {"hg19", "hg38"}:
        raise ContractError(f"harmonized QC lacks valid rows_out/output_build: {trait_id}")
    return rows_out, output_build


def validate_dense_gzip(root: Path, trait_id: str, registry: EvidenceRegistry) -> dict[str, object]:
    relative = f"data/harmonized/{trait_id}.harmonized.tsv.gz"
    qc_rel = f"data/harmonized/{trait_id}.qc.txt"
    registry.add(relative, detail=f"dense harmonized summary statistics: {trait_id}")
    registry.add(qc_rel, detail=f"harmonized QC ledger: {trait_id}")
    expected_rows, output_build = _qc_rows_out(root, trait_id)
    path = _check_components(root, relative)
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(path, flags)
    try:
        before = os.fstat(fd)
        raw = os.fdopen(fd, "rb", closefd=False)
        try:
            with gzip.GzipFile(fileobj=raw, mode="rb") as handle:
                header = handle.readline()
                if header != b"SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\n":
                    raise ContractError(f"dense harmonized header mismatch: {trait_id}")
                rows = 0
                ended_newline = True
                while True:
                    block = handle.read(16 * 1024 * 1024)
                    if not block:
                        break
                    rows += block.count(b"\n")
                    ended_newline = block.endswith(b"\n")
                if not ended_newline or rows != expected_rows:
                    raise ContractError(f"dense harmonized row count/final newline mismatch: {trait_id}: {rows} != {expected_rows}")
        except (OSError, EOFError, gzip.BadGzipFile) as exc:
            raise ContractError(f"invalid/truncated dense gzip for {trait_id}: {exc}") from exc
        finally:
            raw.close()
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    if _identity(before) != _identity(after) or _identity(after) != _identity(final):
        raise ContractError(f"dense gzip drifted while validating: {trait_id}")
    return {"trait_id": trait_id, "rows": rows, "output_build": output_build}


def validate_qtl_metadata(root: Path, policy: Mapping[str, object]) -> dict[str, object]:
    fields, dataset_rows = read_tsv(root, "ref/molecular/eqtl_catalogue_r7/dataset_metadata_r7.tsv")
    expected_fields = [
        "study_id", "dataset_id", "study_label", "sample_group", "tissue_id",
        "tissue_label", "condition_label", "sample_size", "quant_method", "pmid", "study_type",
    ]
    if fields != expected_fields or len(dataset_rows) != 758 or len({row["dataset_id"] for row in dataset_rows}) != 758:
        raise ContractError("eQTL Catalogue r7 dataset metadata schema/count is invalid")
    path_fields, path_rows = read_tsv(root, "ref/molecular/eqtl_catalogue_r7/tabix_ftp_paths.tsv")
    expected_path_fields = expected_fields[:-2] + ["ftp_path", "ftp_cs_path", "ftp_lbf_path"]
    if path_fields != expected_path_fields or len(path_rows) != 758:
        raise ContractError("eQTL Catalogue r7 path metadata schema/count is invalid")
    dataset_by_id = {row["dataset_id"]: row for row in dataset_rows}
    paths_by_id = {row["dataset_id"]: row for row in path_rows}
    expected = {
        "QTD000559": ("Young_2019", "microglia", "ge", "104"),
        "QTD000560": ("Young_2019", "microglia", "exon", "104"),
        "QTD000563": ("Young_2019", "microglia", "leafcutter", "104"),
        "QTD000569": ("Aygun_2021", "neuron", "ge", "73"),
        "QTD000570": ("Aygun_2021", "neuron", "exon", "73"),
        "QTD000573": ("Aygun_2021", "neuron", "leafcutter", "73"),
    }
    for identity, values in expected.items():
        if identity not in dataset_by_id or identity not in paths_by_id:
            raise ContractError(f"required cell QTL metadata is absent: {identity}")
        row = dataset_by_id[identity]
        observed = (row["study_label"], row["tissue_label"], row["quant_method"], row["sample_size"])
        if observed != values:
            raise ContractError(f"required cell QTL context/method mismatch: {identity}")
        if not paths_by_id[identity]["ftp_path"].endswith(".tsv.gz"):
            raise ContractError(f"required cell QTL payload URL is malformed: {identity}")
    imported_fields, imported = read_tsv(root, "ref/molecular/eqtl_catalogue_r7/tabix_ftp_paths_imported.tsv")
    if imported_fields != ["study", "qtl_group", "tissue_ontology_id", "tissue_ontology_term", "tissue_label", "condition_label", "quant_method", "ftp_path"] or len(imported) != 49:
        raise ContractError("GTEx v8 imported QTL metadata is not the exact 49-context family")
    return {"release7_datasets": len(dataset_rows), "release7_path_rows": len(path_rows), "gtex_v8_contexts": len(imported), "required_cell_datasets": sorted(expected)}


def validate_twas_inventory(root: Path, molecular: Mapping[str, object]) -> dict[str, object]:
    fields, rows = read_tsv(root, "results/tables/twas_phi_model_inventory.tsv")
    expected_fields = [
        "file_id", "filename", "context", "file_role", "bytes", "download_url",
        "results_accessed_before_lock",
    ]
    if fields != expected_fields or len(rows) != 98:
        raise ContractError("TWAS inventory is not the exact 98-file family")
    if len({row["file_id"] for row in rows}) != 98 or len({row["filename"] for row in rows}) != 98:
        raise ContractError("TWAS inventory contains duplicate file IDs/names")
    contexts = {row["context"] for row in rows}
    if len(contexts) != 49 or any(row["file_role"] not in {"COVARIANCE", "MODEL_DB"} for row in rows):
        raise ContractError("TWAS inventory context/role family is invalid")
    if any(row["results_accessed_before_lock"] != "NO" for row in rows):
        raise ContractError("TWAS inventory was not frozen result-blind")
    if any(sum(item["context"] == context for item in rows) != 2 for context in contexts):
        raise ContractError("TWAS inventory does not contain one covariance and one model DB per context")
    total = sum(int(row["bytes"]) for row in rows)
    twas = molecular.get("twas")
    if not isinstance(twas, dict) or not isinstance(twas.get("phi_model_source"), dict):
        raise ContractError("TWAS policy/model source block is malformed")
    source = twas["phi_model_source"]
    if len(rows) != source.get("expected_file_count") or len(contexts) != source.get("expected_context_count") or total != source.get("expected_total_bytes"):
        raise ContractError("TWAS inventory counts/bytes differ from molecular policy")
    lock = read_json(root, "results/tables/twas_phi_model_inventory.lock.json")
    if lock.get("inventory_sha256") != stable_file(root, "results/tables/twas_phi_model_inventory.tsv").sha256 or lock.get("file_count") != 98 or lock.get("context_count") != 49 or lock.get("total_bytes") != total or lock.get("results_accessed_before_lock") is not False:
        raise ContractError("TWAS inventory lock is invalid")
    return {"files": 98, "contexts": 49, "expected_bytes": total}


def run_deep_validators(root: Path, interpretation: Mapping[str, object]) -> dict[str, tuple[bool, str]]:
    script = root / "scripts/74_interpretation_preflight.py"
    scripts_dir = str(root / "scripts")
    added = scripts_dir not in sys.path
    if added:
        sys.path.insert(0, scripts_dir)
    try:
        spec = importlib.util.spec_from_file_location("track_b_readiness_interpretation_preflight", script)
        if spec is None or spec.loader is None:
            raise ContractError("cannot load pinned interpretation preflight module")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if added and sys.path and sys.path[0] == scripts_dir:
            sys.path.pop(0)
    validators = {
        "screen_bundle": module.validate_screen_bundle,
        "hocomoco_bundle": module.validate_hocomoco_bundle,
        "abc_bundle": module.validate_abc_bundle,
        "pchic_bundle": module.validate_pchic_bundle,
        "fuma_bundle": module.validate_fuma_bundle,
        "catlas_bundle": module.validate_catlas_bundle,
        "ldsc_source_bundle": module.validate_ldsc_seg_gtex_bundle,
        "pathway_bundle": module.validate_public_pathway_bundle,
        "causal_bundle": module.validate_causal_runtime_bundle,
    }
    results: dict[str, tuple[bool, str]] = {}
    for name, validator in validators.items():
        try:
            ok, reason, _observed = validator(root, interpretation)
            results[name] = (bool(ok), str(reason) if reason else "complete deep semantic validation passed")
        except (Exception, SystemExit) as exc:  # retain a blocker instead of manufacturing readiness
            results[name] = (False, f"deep validator failed: {type(exc).__name__}: {exc}")
        finally:
            gc.collect()
    return results


def validate_ldsc_runtime(root: Path) -> tuple[bool, str]:
    python = _check_components(root, ".ldsc-env/bin/python3.9")
    script = _check_components(root, "ldsc/ldsc.py")
    check = subprocess.run(
        [str(python), "-c", "import numpy,pandas,scipy; print('\\t'.join([numpy.__version__,pandas.__version__,scipy.__version__]))"],
        cwd=root, capture_output=True, text=True, check=False, timeout=30,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    if check.returncode != 0 or check.stdout.strip() != "1.21.5\t1.3.3\t1.7.3":
        return False, "LDSC Python dependency versions failed exact validation"
    help_check = subprocess.run(
        [str(python), str(script), "--help"], cwd=root, capture_output=True,
        text=True, check=False, timeout=30, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    if help_check.returncode != 0 or "--h2-cts" not in (help_check.stdout + help_check.stderr):
        return False, "pinned LDSC entrypoint does not expose --h2-cts"
    return True, "Python 3.9.23-compatible runtime with numpy 1.21.5, pandas 1.3.3, scipy 1.7.3 and LDSC --h2-cts entrypoint validated"


def _evidence_ids_matching(registry: EvidenceRegistry, *needles: str) -> tuple[str, ...]:
    ids = [row["evidence_id"] for row in registry.rows() if any(needle in str(row["path"]) for needle in needles)]
    return tuple(sorted(set(str(value) for value in ids)))


def _probe_fact(
    name: str,
    registry: EvidenceRegistry,
    relative: str,
    blocker_type: str,
    absent_reason: str,
    present_reason: str | None = None,
) -> Fact:
    evidence_id = registry.id_for(relative)
    evidence = registry.get(evidence_id)
    reason = absent_reason if evidence.kind == "MISSING" else (present_reason or f"{relative} exists but has no pre-result checksum/semantic validation and cannot confer readiness")
    return Fact(name, False, blocker_type, reason, (evidence_id,))


def build_facts(
    root: Path,
    policy: Mapping[str, object],
    registry: EvidenceRegistry,
    manifests: Mapping[str, Mapping[str, object]],
) -> dict[str, Fact]:
    validate_build_semantics(manifests)
    fuma_by_id = validate_fuma_species(manifests["interpretation_fuma_scrna"])
    pair_rows, sleep_rows = validate_pair_and_panel(root, policy)
    qtl_observed = validate_qtl_metadata(root, policy)
    twas_observed = validate_twas_inventory(root, manifests["molecular"])

    dense_observed = []
    trait_ids = [row["trait_id"] for row in sleep_rows]
    trait_ids.extend(row["external_trait"] for row in pair_rows.values())
    for trait_id in dict.fromkeys(trait_ids):
        dense_observed.append(validate_dense_gzip(root, trait_id, registry))

    deep = run_deep_validators(root, manifests["interpretation"])
    ldsc_runtime_ok, ldsc_runtime_reason = validate_ldsc_runtime(root)
    probes = policy["absence_probes"]
    assert isinstance(probes, dict)

    facts: dict[str, Fact] = {}
    pinned_panel_ids = _evidence_ids_matching(registry, "config/analysis_panel")
    pinned_pair_ids = _evidence_ids_matching(registry, "results/track_b/pair_manifest")
    dense_ids = _evidence_ids_matching(registry, "data/harmonized/")
    facts["analysis_panel"] = Fact("analysis_panel", True, "DATA", "45-trait panel and exact ordered 12-sleep-trait family validated; mixed build is explicitly 11 hg19 plus sleep_apnea hg38", pinned_panel_ids)
    facts["pair_manifest"] = Fact("pair_manifest", True, "DATA", "exact frozen A, B and CONTROL identities are EUR/hg19 with full-summary readiness", pinned_pair_ids)
    facts["dense_sleep_family"] = Fact("dense_sleep_family", True, "DATA", f"gzip CRC/header/count validation passed for {len(dense_observed)} unique sleep/pair traits; 12-sleep family complete", dense_ids)

    bundle_needles = {
        "screen_bundle": ("interpretation_screen_registry_v4", "regulatory/screen_"),
        "hocomoco_bundle": ("interpretation_hocomoco_v14", "regulatory/hocomoco_"),
        "abc_bundle": ("interpretation_abc_2021", "regulatory/abc_2021"),
        "pchic_bundle": ("interpretation_pchic_2016", "regulatory/javierre_pchic"),
        "fuma_bundle": ("interpretation_fuma_scrna", "cell_type/fuma_scrna", "cell_type/magma_v1.10"),
        "catlas_bundle": ("interpretation_catlas_adult_v4", "cell_type/catlas_adult_v4"),
        "ldsc_source_bundle": ("interpretation_ldsc_seg_", "cell_type/ldsc_seg_gtex_"),
        "pathway_bundle": ("interpretation_pathway_sources", "interpretation/pathway/", "gencode.v26.GRCh38.genes.gtf"),
        "causal_bundle": ("interpretation_causal_", "interpretation/causal/", "ref/eur_w_ld_chr/w_hm3.snplist", "magma_v1.10/g1000_eur.zip", ".mr-env/library", ".r-env/lib/R/library", ".r-env/bin/Rscript"),
    }
    for name, needles in bundle_needles.items():
        ok, reason = deep[name]
        facts[name] = Fact(name, ok, FACT_TYPES[name], reason, _evidence_ids_matching(registry, *needles))
    facts["ldsc_runtime"] = Fact(
        "ldsc_runtime", ldsc_runtime_ok, "SOFTWARE", ldsc_runtime_reason,
        _evidence_ids_matching(registry, ".ldsc-env/bin/python3.9", ".ldsc-env/lib/python3.9/site-packages", "ldsc/ldsc.py", "scripts/98_run_ldsc_seg_task.py"),
    )
    facts["brain_discovery_human"] = Fact(
        "brain_discovery_human", fuma_by_id["Allen_Human_MTG_level2"]["species"] == "human", "DATA",
        "Allen Human MTG level 2 is an exact human middle-temporal-gyrus discovery matrix with 75 cell columns",
        _evidence_ids_matching(registry, "interpretation_fuma_scrna", "FUMA_scRNA_data-"),
    )
    facts["brain_replication_human_class"] = Fact(
        "brain_replication_human_class", fuma_by_id["GSE67835_Human_Cortex_woFetal"]["species"] == "human", "DATA",
        "GSE67835 Human Cortex without fetal samples is an independent human seven-class matrix; it is class-level only",
        _evidence_ids_matching(registry, "interpretation_fuma_scrna", "FUMA_scRNA_data-"),
    )
    facts["qtl_metadata"] = Fact(
        "qtl_metadata", True, "DATA",
        f"eQTL Catalogue r7 metadata validated: {qtl_observed['release7_datasets']} datasets, {qtl_observed['gtex_v8_contexts']} imported GTEx contexts and exact selected microglia/neuron IDs",
        _evidence_ids_matching(registry, "eqtl_catalogue_r7/dataset_metadata", "eqtl_catalogue_r7/tabix_ftp_paths"),
    )
    facts["pathway_executor"] = Fact(
        "pathway_executor", True, "SOFTWARE", "checksum-pinned automatic pathway executor is present",
        _evidence_ids_matching(registry, "scripts/92_run_pathway_task.py"),
    )

    # A file merely appearing at a probe path never upgrades readiness.  A new
    # policy pin and semantic validator are required, preventing pseudo-READY.
    probe_specs = {
        "broad_tissue_magma": ("DATA", "no complete human broad-tissue MAGMA expression matrix is local; FUMA scRNA is not a broad-tissue substitute"),
        "human_subtype_replication": ("DATA", "no independent human subtype-resolved replication atlas is local; GSE67835 has only seven broad classes and mouse is ineligible"),
        "hocomoco_sequence": ("DATA", "no local checksum-pinned hg38 FASTA/sequence family exists; network sequence calls are forbidden for this local contract"),
        "ldsc_static_reference": ("DATA", "the baselineLD v2.2/weights/HapMap3/static EUR reference provenance is absent"),
        "spatial": ("DATA", "no authorized local human spatial payload/metadata manifest exists"),
        "matched_multimodal_grn": ("DATA", "no matched human expression-accessibility dataset/metadata manifest exists for GRN inference"),
        "track_b_mr_executor": ("SOFTWARE", "MR packages and EUR reference are locally valid, but no validated automatic Track B bidirectional MR science executor exists"),
        "conjfdr_runtime": ("SOFTWARE", "no checksum-pinned MATLAB-compatible conjFDR runtime/reference/executor manifest exists"),
        "hdl_l_reference": ("DATA", "no checksum-pinned ancestry/build-matched HDL-L reference manifest exists"),
        "rho_hess_reference": ("DATA", "no separate sealed Track B rho-HESS reference/contract exists"),
        "three_way_coloc_executor": ("SOFTWARE", "no checksum-pinned multi-signal three-way coloc executor exists; pairwise coloc cannot substitute"),
        "regulatory_chain_executor": ("SOFTWARE", "no checksum-pinned same-locus/same-context regulatory-chain integrator exists"),
        "cross_sleep_executor": ("SOFTWARE", "no checksum-pinned complete 12-trait cross-sleep integration executor exists"),
        "grn_executor": ("SOFTWARE", "no validated executor matched to a frozen multimodal GRN design exists"),
        "tabix": ("SOFTWARE", "no repository-local checksum-pinned tabix entrypoint exists"),
    }
    for name, (blocker_type, reason) in probe_specs.items():
        facts[name] = _probe_fact(name, registry, str(probes[name]), blocker_type, reason)
    static_fact = facts["ldsc_static_reference"]
    facts["ldsc_static_reference"] = Fact(
        static_fact.name, static_fact.ready, static_fact.blocker_type, static_fact.reason,
        tuple(sorted(set(static_fact.evidence_ids + (registry.id_for(str(probes["ldsc_static_reference_root"])),)))),
    )

    qtl_ids = _evidence_ids_matching(registry, str(probes["qtl_payload_root"]))
    qtl_present = [registry.get(value) for value in qtl_ids if registry.get(value).kind != "MISSING"]
    facts["qtl_payloads"] = Fact(
        "qtl_payloads", False, "DATA",
        "complete selected QTL payload plus .tbi family is absent" if not qtl_present else "one or more QTL probe files exist but no complete checksum-pinned and semantically validated payload/index family exists",
        qtl_ids + (registry.id_for(str(probes["psychencode_payload_root"])),),
    )
    model_ids = _evidence_ids_matching(registry, str(manifests["molecular"]["twas"]["phi_model_source"]["install_dir"]))
    model_present = [registry.get(value) for value in model_ids if registry.get(value).kind != "MISSING"]
    facts["twas_models"] = Fact(
        "twas_models", False, "DATA",
        f"locked result-blind inventory is valid ({twas_observed['files']} files/{twas_observed['contexts']} contexts/{twas_observed['expected_bytes']} bytes) but model payloads are absent" if not model_present else "TWAS files appeared but lack a complete post-download SHA-256 lock and cannot confer readiness",
        _evidence_ids_matching(registry, "twas_phi_model_inventory", str(manifests["molecular"]["twas"]["phi_model_source"]["install_dir"])),
    )
    runtime_rel = str(manifests["molecular"]["twas"]["runtime_environment"])
    runtime_probe = f"{runtime_rel}" if runtime_rel.endswith("python") else f"{runtime_rel}/python"
    # The molecular runtime path is not an absence-probe key, but adding it now
    # freezes its missing/present-unvalidated state in the evidence family.
    runtime_eid = registry.add(runtime_probe, required=False, detail="required exact MetaXcan Python runtime")
    facts["twas_runtime"] = Fact(
        "twas_runtime", False, "SOFTWARE",
        "MetaXcan source/entrypoint are pinned, but the required exact .molecular-env/python runtime is absent or unvalidated",
        (runtime_eid,) + _evidence_ids_matching(registry, ".molecular-env/MetaXcan", ".molecular-env/source_archives/MetaXcan"),
    )

    missing = set(FACT_TYPES) - set(facts)
    extra = set(facts) - set(FACT_TYPES)
    if missing or extra:
        raise ContractError(f"fact family mismatch; missing={sorted(missing)} extra={sorted(extra)}")
    return facts


def evaluate_analyses(
    policy: Mapping[str, object], facts: Mapping[str, Fact], registry: EvidenceRegistry
) -> list[dict[str, object]]:
    analyses = policy["analyses"]
    assert isinstance(analyses, list)
    result: list[dict[str, object]] = []
    for declared in analyses:
        assert isinstance(declared, dict)
        requirements = [str(value) for value in declared["requirements"]]
        failed = [facts[name] for name in requirements if not facts[name].ready]
        blocker_types = sorted({fact.blocker_type for fact in failed}, key=lambda value: (value != "DATA", value))
        if failed:
            resource_status = "BLOCKED_BY_DATA" if "DATA" in blocker_types else "BLOCKED_BY_SOFTWARE"
            status = resource_status
        else:
            resource_status = "READY"
            status = "NOT_APPLICABLE_UNTIL_UPSTREAM"
        if resource_status not in ALLOWED_STATUSES or status not in ALLOWED_STATUSES:
            raise ContractError(f"derived invalid status for {declared['analysis_id']}")
        evidence_ids: list[str] = []
        for name in requirements:
            evidence_ids.extend(facts[name].evidence_ids)
        evidence_ids = sorted(set(evidence_ids))
        if not evidence_ids:
            raise ContractError(f"analysis has no exact evidence: {declared['analysis_id']}")
        evidence = [registry.get(value) for value in evidence_ids]
        evidence_paths = ";".join(item.path for item in evidence)
        evidence_sha = ";".join(f"{item.path}={item.sha256}" for item in evidence)
        reasons = [fact.reason for fact in failed]
        row = {
            "analysis_id": declared["analysis_id"],
            "phase": declared["phase"],
            "analysis": declared["analysis"],
            "scope": declared["scope"],
            "resource_status": resource_status,
            "status": status,
            "blocker_types": ";".join(blocker_types) if blocker_types else "NONE",
            "status_reason": " | ".join(reasons) if reasons else "local resource/software contract is READY; scientific execution awaits its predeclared upstream gate",
            "exact_unblock_condition": declared["unblock"],
            "required_inputs": declared["required_inputs"],
            "required_software": declared["required_software"],
            "reference": declared["reference"],
            "build": declared["build"],
            "ancestry": declared["ancestry"],
            "requirements": requirements,
            "evidence_ids": evidence_ids,
            "evidence_paths": evidence_paths,
            "evidence_sha256": evidence_sha,
            "ram_decomposition_safety": declared["ram_decomposition"],
            "planned_output": declared["output"],
            "claim_limit": declared["claim_limit"],
        }
        result.append(row)
    if [row["analysis_id"] for row in result] != list(EXPECTED_ANALYSIS_IDS):
        raise ContractError("derived readiness matrix lost or reordered an analysis")
    return result


def _clean_tsv(value: object) -> str:
    if isinstance(value, list):
        value = ";".join(str(item) for item in value)
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def tsv_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: _clean_tsv(row.get(field, "")) for field in fields})
    return output.getvalue().encode("utf-8")


def render_report(
    policy_hash: str,
    evidence_rows: Sequence[Mapping[str, object]],
    readiness: Sequence[Mapping[str, object]],
) -> bytes:
    resource_counts: dict[str, int] = {}
    current_counts: dict[str, int] = {}
    for row in readiness:
        resource_counts[str(row["resource_status"])] = resource_counts.get(str(row["resource_status"]), 0) + 1
        current_counts[str(row["status"])] = current_counts.get(str(row["status"]), 0) + 1
    lines = [
        "# Track B mechanism follow-up: pre-result local readiness",
        "",
        "This is a result-blind resource/software contract. It did not read fine-mapping, PLACO/pleiotropy, or future mechanism results and it makes no post-result scientific claim. Missing, blocked, skipped, and unavailable resources are never interpreted as null evidence.",
        "",
        f"- Policy SHA-256: `{policy_hash}`",
        f"- Complete analysis family: {len(readiness)} records",
        f"- Exact evidence family: {len(evidence_rows)} filesystem records",
        f"- Resource readiness: {', '.join(f'{key}={resource_counts[key]}' for key in sorted(resource_counts))}",
        f"- Current pre-result state: {', '.join(f'{key}={current_counts[key]}' for key in sorted(current_counts))}",
        "",
        "## Readiness matrix",
        "",
        "| Phase | Analysis ID | Resource status | Current status | Exact boundary |",
        "|---:|---|---|---|---|",
    ]
    for row in readiness:
        boundary = str(row["status_reason"]).replace("|", "/")
        lines.append(f"| {row['phase']} | {row['analysis_id']} | {row['resource_status']} | {row['status']} | {boundary} |")
    lines.extend([
        "",
        "## Scientifically important boundaries",
        "",
        "- Human brain discovery is available in Allen Human MTG level 2; independent human replication is available only at seven broad classes in GSE67835. No independent human subtype atlas is local. The two vascular FUMA matrices are mouse mapped to human gene IDs and cannot count as human replication.",
        "- The frozen S-LDSC annotation source family contains 16 GTEx tissues in four domains, but its static genotype/baselineLD/weights/HapMap3 reference family is absent. Lung/airway and skeletal-muscle coverage are not present in that 16-tissue selection.",
        "- eQTL Catalogue and GTEx metadata are present, but the molecular-QTL payloads and matching tabix indexes are not. PsychENCODE remains controlled supplementary access, not a local resource.",
        "- SCREEN, CATlas, ABC and immune-only promoter-capture Hi-C sources are local. HOCOMOCO motifs are local but no immutable local hg38 sequence family is available. No airway-specific or vascular human expression reference was invented.",
        "- The 98-file/49-context phi-enabled TWAS inventory is locked before results, but all model payloads and the exact MetaXcan Python runtime are absent.",
        "- No authorized spatial dataset or matched multimodal human GRN dataset/runtime exists locally. Those absences are blockers, never negative findings.",
        "- MR packages and the 503-sample European GRCh37 LD reference validate, but a production executor enforcing the full bidirectional estimator/overlap/Steiger/colocalization/FDR contract is absent.",
        "- HDL-L, rho-HESS and conjFDR remain blocked by their method-specific runtime/reference contracts.",
        "",
        "## RAM-equivalent decomposition",
        "",
        "Locus-dependent QTL, colocalization, regulatory linking and local-correlation work may run one complete locus per fresh process; an LD-dependent locus must never be split. Trait/cell/pathway analyses may run one frozen trait/dataset/resource unit at a time while preserving every predeclared family and correction denominator. Global nuisance estimation (for example conjFDR or any method that estimates genome-wide parameters) must remain global; chunked tests are permitted only after those parameters are frozen and mathematical identity is demonstrated.",
        "",
        "## Claim boundary",
        "",
        "A READY resource status means only that the local inputs and entrypoints passed this pre-result contract. It is not authorization to bypass an upstream gate, not evidence that an association exists, and not a biological or causal conclusion.",
        "",
    ])
    return "\n".join(lines).encode("utf-8")


def construct_artifacts(root: Path, policy_rel: str = POLICY_REL) -> tuple[dict[str, bytes], dict[str, object]]:
    if policy_rel != POLICY_REL:
        raise ContractError("production CLI accepts only the canonical mechanism policy path")
    policy = read_json(root, policy_rel)
    validate_policy_shape(policy)
    policy_hash = stable_file(root, policy_rel).sha256
    script_hash = stable_file(root, "scripts/160_build_track_b_mechanism_readiness.py").sha256
    registry, manifests = collect_evidence(root, policy)
    facts = build_facts(root, policy, registry, manifests)
    registry.recheck()
    evidence_rows = registry.rows()
    readiness = evaluate_analyses(policy, facts, registry)
    evidence_digest = bytes_sha256(canonical_json(evidence_rows))
    resource_summary: dict[str, int] = {}
    status_summary: dict[str, int] = {}
    for row in readiness:
        resource_summary[str(row["resource_status"])] = resource_summary.get(str(row["resource_status"]), 0) + 1
        status_summary[str(row["status"])] = status_summary.get(str(row["status"]), 0) + 1
    document = {
        "schema_version": SCHEMA,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS",
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": policy_rel,
        "policy_sha256": policy_hash,
        "script_path": "scripts/160_build_track_b_mechanism_readiness.py",
        "script_sha256": script_hash,
        "analysis_family_count": len(readiness),
        "evidence_record_count": len(evidence_rows),
        "evidence_bundle_sha256": evidence_digest,
        "resource_status_counts": dict(sorted(resource_summary.items())),
        "status_counts": dict(sorted(status_summary.items())),
        "facts": [facts[name].record() for name in sorted(facts)],
        "analyses": readiness,
    }
    readiness_tsv_rows = [
        {**row, "exact_unblock_condition": row["exact_unblock_condition"]}
        for row in readiness
    ]
    artifacts = {
        OUTPUT_RELS[0]: tsv_bytes(EVIDENCE_FIELDS, evidence_rows),
        OUTPUT_RELS[1]: tsv_bytes(READINESS_FIELDS, readiness_tsv_rows),
        OUTPUT_RELS[2]: canonical_json(document),
        OUTPUT_RELS[3]: render_report(policy_hash, evidence_rows, readiness),
    }
    lock = {
        "schema_version": SCHEMA,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK",
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": policy_rel,
        "policy_sha256": policy_hash,
        "script_sha256": script_hash,
        "analysis_family_count": len(readiness),
        "evidence_record_count": len(evidence_rows),
        "evidence_bundle_sha256": evidence_digest,
        "outputs": {
            relative: {"bytes": len(payload), "sha256": bytes_sha256(payload)}
            for relative, payload in sorted(artifacts.items())
        },
    }
    artifacts[LOCK_REL] = canonical_json(lock)
    return artifacts, lock


def _secure_mkdirs(root: Path, relative: str) -> Path:
    relative = safe_relative(relative)
    current = root
    for part in PurePosixPath(relative).parts:
        current = current / part
        try:
            info = os.lstat(current)
        except FileNotFoundError:
            os.mkdir(current, 0o755)
            info = os.lstat(current)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise ContractError(f"unsafe output directory component: {current}")
    return current


def freeze_no_replace(root: Path, artifacts: Mapping[str, bytes]) -> None:
    output_dir = _secure_mkdirs(root, OUTPUT_DIR_REL)
    ordered = list(OUTPUT_RELS) + [LOCK_REL]
    if set(artifacts) != set(ordered):
        raise ContractError("artifact family differs from the exact canonical output family")
    for relative in ordered:
        target = root / relative
        try:
            os.lstat(target)
        except FileNotFoundError:
            continue
        raise ContractError(f"no-replace freeze refused existing target: {relative}")
    stage = Path(tempfile.mkdtemp(prefix=".readiness-freeze-", dir=output_dir))
    staged: dict[str, Path] = {}
    try:
        for index, relative in enumerate(ordered):
            path = stage / f"{index:02d}.stage"
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
            try:
                payload = artifacts[relative]
                view = memoryview(payload)
                while view:
                    written = os.write(fd, view)
                    if written <= 0:
                        raise ContractError(f"short write while staging {relative}")
                    view = view[written:]
                os.fsync(fd)
            finally:
                os.close(fd)
            staged[relative] = path
        for relative in ordered:
            target = root / relative
            os.link(staged[relative], target, follow_symlinks=False)
            observed = stable_file(root, relative)
            if observed.observed_bytes != len(artifacts[relative]) or observed.sha256 != bytes_sha256(artifacts[relative]):
                raise ContractError(f"published artifact identity mismatch: {relative}")
        dir_fd = os.open(output_dir, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    finally:
        for path in staged.values():
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
        try:
            os.rmdir(stage)
        except FileNotFoundError:
            pass


def verify_frozen(root: Path) -> None:
    lock_payload = read_stable_bytes(root, LOCK_REL)
    try:
        lock = json.loads(lock_payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"invalid readiness lock JSON: {exc}") from exc
    if not isinstance(lock, dict) or canonical_json(lock) != lock_payload:
        raise ContractError("readiness lock is not canonical JSON")
    if lock.get("schema_version") != SCHEMA or lock.get("result_blind") is not True or lock.get("future_results_accessed") is not False:
        raise ContractError("readiness lock schema/result-blind invariants failed")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or list(sorted(outputs)) != list(sorted(OUTPUT_RELS)):
        raise ContractError("readiness lock output family is incomplete")
    for relative in OUTPUT_RELS:
        item = stable_file(root, relative)
        expected = outputs[relative]
        if not isinstance(expected, dict) or item.observed_bytes != expected.get("bytes") or item.sha256 != expected.get("sha256"):
            raise ContractError(f"frozen readiness output differs from lock: {relative}")
    recomputed, expected_lock = construct_artifacts(root)
    for relative in OUTPUT_RELS:
        frozen = read_stable_bytes(root, relative)
        if frozen != recomputed[relative]:
            raise ContractError(f"current local resources no longer reproduce frozen artifact: {relative}")
    if lock != expected_lock or lock_payload != recomputed[LOCK_REL]:
        raise ContractError("current local resources no longer reproduce the frozen lock")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true", help="create the canonical result-blind artifacts with no replacement")
    mode.add_argument("--verify", action="store_true", help="verify frozen artifacts and current resources without writes")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        raise ContractError(f"repository root is not a directory: {root}")
    if args.freeze:
        artifacts, _lock = construct_artifacts(root)
        freeze_no_replace(root, artifacts)
        print(f"FROZEN {len(EXPECTED_ANALYSIS_IDS)} analyses at {OUTPUT_DIR_REL}")
    else:
        verify_frozen(root)
        print(f"VERIFIED {len(EXPECTED_ANALYSIS_IDS)} analyses at {OUTPUT_DIR_REL}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ContractError as exc:
        raise SystemExit(f"ERROR: {exc}")
