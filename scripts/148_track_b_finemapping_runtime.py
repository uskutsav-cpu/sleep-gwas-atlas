#!/usr/bin/env python3
"""Strict Track B fine-mapping materialization and validation primitives.

This module is the additive production boundary between the frozen
pre-materialization pair x official-LAVA-block family emitted by script 145
and the fresh-process coordinator in script 151.  It deliberately does not
discover or rank loci and never reads a SuSiE/coloc result while constructing
the materialized manifest.

The public functions are intentionally importable by downstream verifiers:

* :func:`validate_pre_materialization_family`
* :func:`materialize_family`
* :func:`validate_materialized_manifest`
* :func:`deep_validate_locus_task`
* :func:`validate_engine_bundle`
* :func:`validate_canonical_family`

All publication operations are exclusive/no-replace and all verification
operations are side-effect free.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import stat
import subprocess
import sys
import tempfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
POLICY_REL = Path("config/track_b_finemapping_policy.json")
DENSE_QC_REL = Path("results/track_b/03_dense_input_qc.tsv")
DENSE_LOCK_REL = Path("results/track_b/03_dense_input_qc.lock.json")
PANEL_REL = Path("config/analysis_panel.tsv")
ARCHIVE_MANAGER_REL = Path("scripts/146_manage_lava_reference_archives.py")
PRE_GATE_REL = Path("scripts/145_build_track_b_finemapping_continuation_gate.py")
PRE_MANIFEST_REL = Path("results/track_b/DEEP_LOCUS_MANIFEST.tsv")
PRE_LOCK_REL = Path("results/track_b/DEEP_LOCUS_MANIFEST.lock.json")
PRE_UNAVAILABLE_REL = Path("results/track_b/DEEP_LOCUS_UNAVAILABLE.tsv")
MATERIALIZED_MANIFEST_REL = Path("results/track_b/finemapping/locus_manifest.tsv")
MATERIALIZED_LOCK_REL = Path("results/track_b/finemapping/locus_manifest.lock.json")
ZERO_PROVENANCE_REL = Path("results/track_b/finemapping/zero_family.provenance.json")
INPUT_ROOT_REL = Path("work/track_b_finemapping/inputs")

ANALYSIS_ID = "track-b-v1.0-finemapping-trait-coloc"
PRE_FAMILY_SCHEMA = "sleep-atlas-track-b-finemapping-deep-locus-family.1"
PRE_FAMILY_ROLE = "PRE_MATERIALIZATION_PAIR_BY_OFFICIAL_BLOCK_FAMILY"
MATERIALIZED_LOCK_SCHEMA = "sleep-atlas-track-b-finemapping-materialized-family.1"
TASK_SCHEMA = "sleep-atlas-track-b-finemapping-locus-task.2"
LD_VALIDATION_SCHEMA = "sleep-atlas-track-b-signed-ld-validation-attestation.2"
LD_RUNTIME_ATTESTATION_SCHEMA = "sleep-atlas-track-b-signed-ld-runtime-attestation.1"
ZERO_SCHEMA = "sleep-atlas-track-b-finemapping-zero-family-production.1"
READY_STATE = "READY_WHOLE_LOCUS_MATERIALIZATION"
ZERO_STATE = "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY"
REFERENCE_STATES = {
    "ARCHIVES_PRESENT_FULLY_VERIFIED",
    (
        "ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_"
        "EXTRACTED_PAYLOADS_REHASHED"
    ),
}
PAIR_ORDER = ("A", "B", "CONTROL")
PAIR_IDENTITIES = {
    "A": ("snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    "B": ("insomnia", "adhd", "PRIMARY_DISCOVERY"),
    "CONTROL": ("insomnia", "frailty", "POSITIVE_CONTROL_NON_NOVELTY"),
}
TRAIT_TYPES = {
    "snoring": "cc", "insomnia": "cc", "adhd": "cc",
    "parental_lifespan": "quant", "frailty": "quant",
}
BASES = frozenset("ACGT")
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
RSID = re.compile(r"^rs[0-9]+$", re.IGNORECASE)
HEX64 = re.compile(r"^[0-9a-f]{64}$")
MAX_STRUCTURED_BYTES = 128 * 1024**2
BASE_PHYSICAL_ENVELOPE_BYTES = 8 * 1024**3

PRE_FAMILY_FIELDS = [
    "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
    "CHR", "START", "STOP", "ld_block_id", "inclusion_sources", "priority_tier",
    "local_evidence_ids", "local_lead_variants", "pleiotropy_lead_evidence_ids",
    "pleiotropy_lead_variants", "pleiotropy_clump_evidence_ids",
    "method_availability", "partition_role",
    "selection_rule", "claim_limit",
]

PRE_UNAVAILABLE_FIELDS = [
    "analysis_id", "unavailable_entry_id", "pair_id", "family_role", "trait1",
    "trait2", "CHR", "BP", "lead_variant", "ld_block_id", "evidence_ids",
    "inclusion_sources", "terminal_status", "blocker_reason", "claim_limit",
]

SUMMARY_FIELDS = [
    "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "INFO", "N",
]
ORDER_FIELDS = ["SNP", "CHR", "BP", "A1", "A2"]
LD_QC_FIELDS = [
    "reference_id", "reference_sample_size", "variant_count", "symmetry_max_abs",
    "diagonal_max_abs", "minimum_eigenvalue", "maximum_abs_correlation",
    "matrix_kind", "variant_order_sha256", "status",
]
TASK_FIELDS = [
    "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
    "CHR", "START", "STOP", "ld_block_id", "variant_count", "summary1_path",
    "summary1_sha256", "summary2_path", "summary2_sha256", "variant_order_path",
    "variant_order_sha256", "signed_ld_path", "signed_ld_sha256", "ld_qc_path",
    "ld_qc_sha256", "ld_validation_path", "ld_validation_sha256",
    "ld_runtime_package_lock_path", "ld_runtime_package_lock_sha256",
    "ld_runtime_attestation_path", "ld_runtime_attestation_sha256",
    "runtime_package_lock_path", "runtime_package_lock_sha256",
    "execution_amendment_sha256", "trait1_type", "trait1_ncase", "trait1_ncontrol",
    "trait1_case_fraction", "trait1_sdY", "trait1_scalar_N", "trait1_scalar_N_rule",
    "trait1_per_snp_N_sha256", "trait1_per_snp_N_count", "trait1_per_snp_N_min",
    "trait1_per_snp_N_q1", "trait1_per_snp_N_median", "trait1_per_snp_N_q3",
    "trait1_per_snp_N_max", "trait1_per_snp_N_max_to_min_ratio",
    "trait1_per_snp_N_fraction_below_90pct_max",
    "trait1_per_snp_N_fraction_below_50pct_max", "trait1_N_dispersion_status",
    "trait1_INFO_status", "trait2_type", "trait2_ncase",
    "trait2_ncontrol", "trait2_case_fraction", "trait2_sdY", "trait2_scalar_N",
    "trait2_scalar_N_rule", "trait2_per_snp_N_sha256", "trait2_per_snp_N_count",
    "trait2_per_snp_N_min", "trait2_per_snp_N_q1", "trait2_per_snp_N_median",
    "trait2_per_snp_N_q3", "trait2_per_snp_N_max",
    "trait2_per_snp_N_max_to_min_ratio",
    "trait2_per_snp_N_fraction_below_90pct_max",
    "trait2_per_snp_N_fraction_below_50pct_max", "trait2_N_dispersion_status",
    "trait2_INFO_status",
    "reference_id", "reference_state", "reference_sample_size", "prior_method",
    "maximum_causal_signals", "credible_set_coverage", "minimum_absolute_correlation",
    "maximum_iterations", "estimate_residual_variance", "p1", "p2", "p12_grid",
    "single_signal_fallback", "maximum_locus_variants", "locus_splitting",
    "variant_thinning", "lead_centered_truncation", "claim_limit",
]

RUNTIME_PACKAGE_LOCK_FIELDS = [
    "record_type", "name", "version", "path", "file_count", "bytes", "sha256",
]
LD_RUNTIME_ATTESTATION_FIELDS = [
    "schema_version", "analysis_id", "locus_entry_id", "reference_id",
    "execution_amendment_sha256", "runtime_package_lock_sha256",
    "precheck_status", "postcheck_status", "loaded_namespace_count",
    "loaded_namespace_closure", "lava_version",
]

CREDIBLE_SET_FIELDS = [
    "analysis_id", "pair_id", "family_role", "locus_entry_id", "trait_id",
    "trait_role", "signal_id", "component_index", "lead_snp", "lead_pip",
    "credible_set_size", "credible_set_snps", "requested_coverage",
    "achieved_coverage", "min_abs_corr", "mean_abs_corr", "median_abs_corr",
    "cs_log10bf", "model_converged", "analysis_status", "claim_limit",
]

FINEMAP_QC_FIELDS = [
    "analysis_id", "pair_id", "family_role", "locus_entry_id", "trait_id",
    "trait_role", "variant_count", "model_converged", "niter", "credible_set_count",
    "max_pip", "rss_ld_s", "kriging_allele_switch_outlier_count",
    "ld_symmetry_max_abs", "ld_diagonal_max_abs", "ld_minimum_eigenvalue",
    "reference_sample_size", "analysis_status", "diagnostic_status", "error",
]

N_DIAGNOSTIC_FIELDS = [
    "analysis_id", "pair_id", "family_role", "locus_entry_id", "dataset_role",
    "trait_id", "per_snp_N_sha256", "per_snp_N_count", "per_snp_N_min",
    "per_snp_N_q1", "per_snp_N_median", "per_snp_N_q3", "per_snp_N_max",
    "per_snp_N_max_to_min_ratio", "per_snp_N_fraction_below_90pct_max",
    "per_snp_N_fraction_below_50pct_max", "scalar_N", "scalar_N_rule",
    "dispersion_status", "claim_cap",
]
RUNTIME_ATTESTATION_FIELDS = [
    "analysis_id", "locus_entry_id", "execution_amendment_sha256",
    "runtime_package_lock_sha256", "precheck_status", "postcheck_status",
    "loaded_namespace_count",
]

ENGINE_STATUS_FIELDS = [
    "analysis_id", "pair_id", "locus_entry_id", "engine_status", "trait1_status",
    "trait2_status", "coloc_status", "error",
]

SCALAR_N_RULE_CONSTANT = "EXACT_CONSTANT_ELIGIBLE_PER_SNP_N"
SCALAR_N_RULE_MEDIAN = "PRE_RESULT_MEDIAN_ELIGIBLE_PER_SNP_N_CLAIM_CAPPED"
N_DISPERSION_CONSTANT = "EXACT_CONSTANT_PER_SNP_N"
N_DISPERSION_VARIABLE = "VARIABLE_PER_SNP_N_MEDIAN_SCALAR_APPROXIMATION_CLAIM_CAPPED"
SCALAR_N_CONVENTION = (
    "Use the exact value when every eligible per-SNP N is numerically identical; otherwise "
    "use the ordinary sample median of all eligible per-SNP N values. The rule is fixed before "
    "result access and is never rounded or tuned using SuSiE or coloc output. Because pinned "
    "coloc::runsusie/susieR::susie_rss accepts only one scalar n, every nonconstant vector is "
    "retained and quantified but the resulting fine-mapping strength and coloc classification "
    "are claim-capped for sample-size uncertainty."
)
PINNED_RUNTIME_NAMESPACE_CLOSURE = (
    "base", "cli", "coloc", "compiler", "crayon", "data.table", "datasets",
    "dplyr", "farver", "generics", "ggplot2", "glue", "graphics", "grDevices",
    "grid", "gridExtra", "gtable", "irlba", "lattice", "lifecycle", "magrittr",
    "Matrix", "matrixStats", "methods", "mixsqp", "pillar", "pkgconfig", "plyr",
    "R6", "RColorBrewer", "Rcpp", "reshape", "rlang", "S7", "scales", "stats",
    "susieR", "tibble", "tidyselect", "tools", "utils", "vctrs", "viridis",
    "viridisLite",
)
PINNED_LD_RUNTIME_NAMESPACE_CLOSURE = (
    "base", "compiler", "data.table", "datasets", "graphics", "grDevices",
    "LAVA", "methods", "stats", "tools", "utils",
)
EXECUTION_AMENDMENT = {
    "schema_version": "sleep-atlas-track-b-finemapping-execution-amendment.1",
    "amends_policy_schema": "sleep-atlas-track-b-finemapping-policy.1",
    "timing": "BEFORE_ANY_TRACK_B_FINE_MAPPING_OR_COLOC_RESULT_ACCESS",
    "scientific_family_or_model_changed": False,
    "whole_official_block_and_full_eligible_snp_family_required": True,
    "signed_ld_spectral_validation": (
        "FULL_MATRIX_EXACT_EIGVALSH_ONLY_IN_MEASURED_FRESH_MATERIALIZATION_PROCESS;"
        "LATER_READ_ONLY_VERIFICATION_USES_IMMUTABLE_HASHED_ATTESTATION"
    ),
    "runtime_package_binding": (
        "RSCRIPT_IDENTITY_PLUS_EXACT_ACTUALLY_LOADED_LAVA_AND_SUSIER_COLOC_NAMESPACE_CLOSURES;"
        "EXACT_NAMESPACE_PATH_AND_DETERMINISTIC_EVERY_INSTALLED_FILE_TREE_DIGEST;"
        "EACH_PROCESS_CHECKS_BEFORE_SCIENTIFIC_WORK_AND_AFTER_OUTPUT_CONSTRUCTION;"
        "PYTHON_RECHECKS_EXACT_LOCKS_AFTER_PROCESS_EXIT"
    ),
    "nonconstant_per_snp_n": N_DISPERSION_VARIABLE,
    "sample_size_software_semantics_basis": (
        "PINNED_SUSIER_0.14.2_SUSIE_RSS_DOCUMENTS_ONE_SCALAR_N_AS_A_REASONABLE_ESTIMATE;"
        "PINNED_COLOC_5.2.3_RUNSUSIE_FORWARDS_ONE_DATASET_N;NO_PER_SNP_N_LIKELIHOOD_"
        "IS_CLAIMED_OR_INVENTED"
    ),
    "nonconstant_per_snp_n_claim_cap": (
        "RETAIN_FULL_VECTOR_AND_HASH;REPORT_COUNT_MIN_Q1_MEDIAN_Q3_MAX_MAX_TO_MIN_RATIO_"
        "AND_FRACTIONS_BELOW_90_AND_50_PERCENT_OF_LOCUS_MAX;USE_FIXED_MEDIAN_ONLY_BECAUSE_"
        "PINNED_RUNSUSIE_ACCEPTS_SCALAR_N;ADD_SAMPLE_SIZE_UNCERTAINTY_AND_FORCE_COLOC_"
        "INCONCLUSIVE"
    ),
    "complement_alleles": (
        "ONLY_UNAMBIGUOUS_NONPALINDROMIC_STRAND_EQUIVALENTS;COMPLEMENT_ONLY_PRESERVES_"
        "BETA_SIGN_AND_EAF;COMPLEMENT_PLUS_SWAP_NEGATES_BETA_AND_REPLACES_EAF_WITH_1_MINUS_EAF;"
        "TRUE_MISMATCH_FAILS_WHOLE_LOCUS"
    ),
    "resource_retry": (
        "BLOCKED_OR_OOM_PREDECESSORS_ONLY_FOR_MATERIALIZATION_OR_R_ENGINE;EXPLICIT_"
        "CONTINUATION;STRICTLY_LARGER_HOST_PHYSICAL_MEMORY_AND_EFFECTIVE_ENVELOPE;"
        "UNCHANGED_FINGERPRINT;IMMUTABLE_PREDECESSOR_SEAL_AND_PROVENANCE_IDENTITIES"
    ),
}
EXECUTION_AMENDMENT_SHA256 = hashlib.sha256(
    json.dumps(
        EXECUTION_AMENDMENT, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()
).hexdigest()
RSS_LD_S_MAXIMUM = 0.10
KRIGING_LOG_LR_THRESHOLD = 2.0
KRIGING_ABS_Z_THRESHOLD = 2.0

RAW_VARIANT_FIELDS = [
    "analysis_id", "pair_id", "family_role", "locus_entry_id", "locus_id", "CHR",
    "START", "STOP", "trait_id", "trait_role", "SNP", "BP", "A1", "A2", "BETA",
    "SE", "P", "EAF", "INFO", "N", "prior_role", "prior_weight", "PIP",
    "credible_set_ids", "credible_set_sizes", "credible_set_requested_coverage",
    "credible_set_achieved_coverage", "credible_set_min_abs_corr",
    "credible_set_mean_abs_corr", "max_alpha_component", "model_converged", "niter",
    "rss_ld_s", "kriging_allele_switch_outlier", "diagnostic_status", "error",
    "claim_limit",
]

RAW_COLOC_FIELDS = [
    "analysis_id", "pair_id", "family_role", "locus_entry_id", "locus_id", "CHR",
    "START", "STOP", "trait1", "trait2", "signal1", "signal2", "coloc_method", "p1",
    "p2", "p12", "prior_role", "nsnps", "PP_H0", "PP_H1", "PP_H2", "PP_H3",
    "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant",
    "top_shared_variant_PP_H4", "fine_mapping_qc", "ld_qc", "engine_status", "error",
    "claim_limit",
]


class FineMappingError(RuntimeError):
    """Fail-closed contract, validation, or publication error."""


class ComputeBlocked(FineMappingError):
    """The unchanged indivisible whole locus cannot be admitted on this host."""


def canonical_json_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest_json(value: object) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _relative(root: Path, path: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError as error:
        raise FineMappingError(f"path escapes repository: {path}") from error


def safe_path(root: Path, value: Path | str, label: str, *, must_exist: bool = False) -> Path:
    root = root.resolve(strict=True)
    raw = Path(value)
    if raw.is_absolute():
        try:
            relative = raw.relative_to(root)
        except ValueError as error:
            raise FineMappingError(f"unsafe {label} path outside repository: {raw}") from error
    else:
        relative = raw
    if not relative.parts or ".." in relative.parts:
        raise FineMappingError(f"unsafe {label} path: {value}")
    current = root
    for index, part in enumerate(relative.parts):
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            if must_exist:
                raise FineMappingError(f"missing {label}: {current}")
            break
        except OSError as error:
            raise FineMappingError(f"cannot inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise FineMappingError(f"{label} contains a symbolic link: {current}")
        if index < len(relative.parts) - 1 and not stat.S_ISDIR(observed.st_mode):
            raise FineMappingError(f"{label} has a non-directory ancestor: {current}")
    try:
        current = (root / relative).resolve(strict=must_exist)
        current.relative_to(root)
    except (OSError, RuntimeError, ValueError) as error:
        raise FineMappingError(f"unsafe {label} path: {value}") from error
    return root / relative


def stable_identity(
    root: Path, value: Path | str, *, allow_empty: bool = False,
    require_single_link: bool = True,
) -> dict[str, object]:
    path = safe_path(root, value, "artifact", must_exist=True)
    digest = hashlib.sha256()
    descriptor: int | None = None
    try:
        named_before = os.lstat(path)
        flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or (opened.st_size == 0 and not allow_empty)
            or (require_single_link and opened.st_nlink != 1)
        ):
            raise FineMappingError(f"artifact is not one private regular file: {path}")
        while True:
            block = os.read(descriptor, 4 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
        opened_after = os.fstat(descriptor)
        named_after = os.lstat(path)
    except FineMappingError:
        raise
    except OSError as error:
        raise FineMappingError(f"cannot hash artifact {path}: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
    keys = lambda item: (  # noqa: E731 - compact immutable-identity tuple
        item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns,
        item.st_nlink,
    )
    if len({keys(item) for item in (named_before, opened, opened_after, named_after)}) != 1:
        raise FineMappingError(f"artifact changed while hashing: {path}")
    return {
        "path": _relative(root, path), "bytes": int(opened.st_size),
        "sha256": digest.hexdigest(),
    }


def identity_matches(observed: object, expected: object) -> bool:
    return isinstance(observed, Mapping) and isinstance(expected, Mapping) and all(
        observed.get(key) == expected.get(key) for key in ("path", "bytes", "sha256")
    )


def stable_bytes(root: Path, value: Path | str, *, allow_empty: bool = False) -> tuple[bytes, dict[str, object]]:
    identity = stable_identity(root, value, allow_empty=allow_empty)
    if int(identity["bytes"]) > MAX_STRUCTURED_BYTES:
        raise FineMappingError(f"structured artifact is too large: {identity['path']}")
    path = safe_path(root, str(identity["path"]), "structured artifact", must_exist=True)
    try:
        content = path.read_bytes()
    except OSError as error:
        raise FineMappingError(f"cannot read structured artifact {path}: {error}") from error
    if hashlib.sha256(content).hexdigest() != identity["sha256"]:
        raise FineMappingError(f"structured artifact changed after hashing: {path}")
    return content, identity


def read_json(root: Path, value: Path | str) -> tuple[dict[str, Any], dict[str, object]]:
    content, identity = stable_bytes(root, value)
    try:
        payload = json.loads(content.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise FineMappingError(f"unreadable JSON artifact {value}: {error}") from error
    if not isinstance(payload, dict):
        raise FineMappingError(f"JSON artifact must contain one object: {value}")
    return payload, identity


def read_tsv(
    root: Path, value: Path | str, expected_fields: Sequence[str] | None = None,
    *, allow_header_only: bool = False,
) -> tuple[list[str], list[dict[str, str]], dict[str, object]]:
    content, identity = stable_bytes(root, value)
    try:
        reader = csv.DictReader(io.StringIO(content.decode("utf-8"), newline=""), delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    except (UnicodeError, csv.Error) as error:
        raise FineMappingError(f"unreadable TSV artifact {value}: {error}") from error
    if expected_fields is not None and fields != list(expected_fields):
        raise FineMappingError(f"TSV schema drifted for {value}: {fields}")
    if not fields or (not rows and not allow_header_only):
        raise FineMappingError(f"TSV is empty/header-only: {value}")
    if any(None in row or any(item is None for item in row.values()) for row in rows):
        raise FineMappingError(f"malformed TSV row in {value}")
    return fields, rows, identity


def tsv_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def deterministic_gzip(content: bytes) -> bytes:
    return gzip.compress(content, compresslevel=6, mtime=0)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_file(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        observed = os.fstat(descriptor)
        if not stat.S_ISREG(observed.st_mode):
            raise FineMappingError(f"cannot fsync non-regular file: {path}")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def ensure_directory(root: Path, relative: Path | str) -> Path:
    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise FineMappingError(f"unsafe output directory: {relative}")
    current = root.resolve(strict=True)
    for part in relative.parts:
        candidate = current / part
        try:
            os.mkdir(candidate)
        except FileExistsError:
            observed = os.lstat(candidate)
            if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
                raise FineMappingError(f"output path is not a real directory: {candidate}")
        else:
            fsync_directory(current)
        current = candidate
    fsync_directory(current)
    return current


def publish_bytes_no_replace(root: Path, relative: Path | str, content: bytes) -> dict[str, object]:
    relative = Path(relative)
    parent = ensure_directory(root, relative.parent)
    path = safe_path(root, relative, "publication")
    if path.exists() or path.is_symlink():
        raise FineMappingError(f"refusing to replace immutable artifact: {relative}")
    temporary = parent / f".{path.name}.{os.getpid()}.{os.urandom(6).hex()}.tmp"
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        descriptor = os.open(temporary, flags, 0o600)
        try:
            offset = 0
            while offset < len(content):
                offset += os.write(descriptor, content[offset:])
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        os.link(temporary, path, follow_symlinks=False)
        fsync_directory(parent)
    except FileExistsError as error:
        raise FineMappingError(f"refusing to replace immutable artifact: {relative}") from error
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return stable_identity(root, relative)


def _number(value: Any, label: str) -> float:
    try:
        observed = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise FineMappingError(f"invalid numeric {label}: {value}") from error
    if not math.isfinite(observed):
        raise FineMappingError(f"nonfinite numeric {label}: {value}")
    return observed


def _integer(value: Any, label: str, *, minimum: int | None = None) -> int:
    try:
        observed = int(str(value))
    except (TypeError, ValueError) as error:
        raise FineMappingError(f"invalid integer {label}: {value}") from error
    if str(observed) != str(value).strip() or (minimum is not None and observed < minimum):
        raise FineMappingError(f"invalid integer {label}: {value}")
    return observed


def _probability(value: Any, label: str) -> float:
    observed = _number(value, label)
    if observed < 0 or observed > 1:
        raise FineMappingError(f"probability outside [0,1] for {label}: {value}")
    return observed


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise FineMappingError(f"could not load module: {path}")
    module = importlib.util.module_from_spec(specification)
    # Contract verification must not mutate the repository merely because a
    # source module is imported (in particular, no new __pycache__ artifact).
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        specification.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


def verify_pre_gate_results(root: Path = ROOT) -> dict[str, Any]:
    """Invoke revised script 145's side-effect-free authoritative verifier.

    Directly hashing the two frozen files is necessary but not sufficient: the
    corrected 145 verifier also revalidates the complete, sealed 147 evidence
    family from which they were derived.  This adapter therefore fails closed
    until that revised API is present and returns the exact agreed interface.
    """

    gate = safe_path(root, PRE_GATE_REL, "pre-family gate", must_exist=True)
    module = _load_module("track_b_corrected_pre_family_gate_for_148", gate)
    verifier = getattr(module, "verify_results", None)
    if not callable(verifier):
        raise FineMappingError(
            "corrected script-145 verify_results(root=...) API is unavailable"
        )
    try:
        result = verifier(root=root)
    except BaseException as error:
        raise FineMappingError(f"corrected script-145 verification failed: {error}") from error
    required = {
        "state", "rows", "unavailable_rows", "family_path", "family_lock_path",
        "family_identity", "family_lock_identity", "unavailable_path",
        "unavailable_lock_path", "unavailable_identity", "unavailable_lock_identity",
        "zero_provenance_path", "zero_provenance_identity", "mapped_count",
        "unavailable_count", "upstream_evidence_denominator", "family_sha256",
        "unavailable_sha256",
    }
    if not isinstance(result, dict) or set(result) != required:
        raise FineMappingError("corrected script-145 verify_results return schema drifted")
    if result.get("state") not in {READY_STATE, ZERO_STATE}:
        raise FineMappingError("corrected script-145 verifier returned an unsupported state")
    return result


def load_policy(root: Path = ROOT) -> tuple[dict[str, Any], dict[str, object]]:
    policy, identity = read_json(root, POLICY_REL)
    fine = policy.get("fine_mapping_contract", {})
    coloc = policy.get("colocalization_contract", {})
    dense = policy.get("dense_input_contract", {})
    signed = policy.get("signed_ld_contract", {})
    execution = policy.get("ram_aware_execution_contract", {})
    if (
        policy.get("schema_version") != "sleep-atlas-track-b-finemapping-policy.1"
        or policy.get("analysis_id") != ANALYSIS_ID
        or policy.get("analysis_build") != "hg19"
        or policy.get("ancestry") != "EUR"
        or fine.get("primary_method") != "SuSiE-RSS"
        or fine.get("maximum_causal_signals") != 10
        or fine.get("credible_set_coverage") != 0.95
        or fine.get("minimum_absolute_correlation_for_credible_set") != 0.5
        or fine.get("maximum_iterations") != 1000
        or fine.get("primary_prior") != "FLAT"
        or fine.get("estimate_residual_variance") is not False
        or coloc.get("primary_method") != "coloc.susie"
        or coloc.get("p1") != 1e-4
        or coloc.get("p2") != 1e-4
        or coloc.get("p12_primary") != 1e-5
        or coloc.get("p12_sensitivity_grid") != [1e-6, 5e-6, 1e-5, 5e-5]
        or coloc.get("single_signal_fallback") != "FORBIDDEN for the primary Track B trait-trait family; absence of a valid SuSiE credible set is retained as INCONCLUSIVE or INVALID_INPUT according to diagnostics"
        or dense.get("minimum_maf") != 0.01
        or dense.get("minimum_locus_variants") != 50
        or dense.get("maximum_locus_variants") is not None
        or signed.get("reference_id") != "LAVA_UKB_v1.1_EUR_GRCh37"
        or signed.get("minimum_reference_sample_size") != 10000
        or signed.get("symmetry_tolerance") != 1e-8
        or signed.get("diagonal_tolerance") != 1e-6
        or signed.get("minimum_eigenvalue_tolerance") != -1e-6
        or execution.get("maximum_concurrent_loci_per_worker") != 1
        or execution.get("locus_splitting_forbidden") is not True
        or execution.get("variant_thinning_for_compute_forbidden") is not True
    ):
        raise FineMappingError("Track B fine-mapping policy/model identity drifted")
    software = policy.get("software_contract", {})
    for package in ("susieR", "coloc"):
        contract = software.get(package, {})
        archive = contract.get("source_archive")
        expected = contract.get("source_archive_sha256")
        if not isinstance(archive, str) or not isinstance(expected, str):
            raise FineMappingError(f"policy lacks pinned {package} source archive")
        if stable_identity(root, archive)["sha256"] != expected:
            raise FineMappingError(f"pinned {package} source archive hash drifted")
    return policy, identity


def verify_reference_state(root: Path = ROOT) -> tuple[dict[str, Any], dict[str, object]]:
    """Call script 146's only supported reference-state API and bind its hash."""

    manager = safe_path(root, ARCHIVE_MANAGER_REL, "reference-state manager", must_exist=True)
    manager_identity = stable_identity(root, ARCHIVE_MANAGER_REL)
    module = _load_module("track_b_reference_state_for_finemapping_148", manager)
    try:
        state = module.verify_reference_state(root)
    except BaseException as error:
        raise FineMappingError(f"LAVA reference state is unusable: {error}") from error
    if not isinstance(state, dict) or state.get("state") not in REFERENCE_STATES:
        raise FineMappingError("script 146 returned an unsupported LAVA reference state")
    if state.get("extracted_payload_count") != 44:
        raise FineMappingError("script 146 did not rehash the exact 44 extracted payload family")
    return state, manager_identity


def _identity_record(value: Any, label: str) -> Mapping[str, Any]:
    if (
        not isinstance(value, Mapping)
        or set(value) != {"path", "bytes", "sha256"}
        or not isinstance(value.get("path"), str)
        or not isinstance(value.get("bytes"), int)
        or isinstance(value.get("bytes"), bool)
        or int(value["bytes"]) <= 0
        or not isinstance(value.get("sha256"), str)
        or not HEX64.fullmatch(str(value["sha256"]))
    ):
        raise FineMappingError(f"malformed immutable identity for {label}")
    return value


def _split_tokens(value: str, label: str) -> list[str]:
    if value == "NA":
        return []
    tokens = value.split(";")
    if any(not token or token.strip() != token for token in tokens) or len(tokens) != len(set(tokens)):
        raise FineMappingError(f"{label} is empty, whitespace-drifted, or duplicated")
    return tokens


def _validate_selection_contract(value: Any) -> None:
    expected = {
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
    if value != expected:
        raise FineMappingError("pre-materialization selection contract drifted")


def validate_pre_materialization_family(
    root: Path = ROOT, manifest_rel: Path | str = PRE_MANIFEST_REL,
    lock_rel: Path | str = PRE_LOCK_REL,
) -> dict[str, Any]:
    """Validate only the revised script-145 frozen interface, never live results."""

    # This must precede acceptance of any caller-supplied path.  Script 145 is
    # the authority that revalidates frozen 147, while the checks below
    # independently rehash and semantically validate its canonical projection.
    gate_result = verify_pre_gate_results(root)
    policy, policy_identity = load_policy(root)
    lock, lock_identity = read_json(root, lock_rel)
    # Fail explicitly on the superseded/current 145 shape.  This prevents a
    # historical gate from being treated as the corrected production family.
    if lock.get("schema_version") != PRE_FAMILY_SCHEMA:
        raise FineMappingError(
            "unsupported script-145 interface; corrected frozen DEEP_LOCUS_MANIFEST is required"
        )
    required_keys = {
        "schema_version", "analysis_id", "artifact_role", "state",
        "results_accessed_before_lock", "pair_order", "manifest", "locus_count",
        "locus_entry_ids_in_order", "evidence_count", "evidence_ids_sha256", "policy",
        "official_blocks", "upstream_family_provenance", "generator_scripts",
        "selection_contract", "zero_family_provenance", "unavailable_manifest",
        "unavailable_count", "unavailable_entry_ids_in_order",
        "unavailable_evidence_count", "unavailable_evidence_ids_sha256",
        "upstream_evidence_count", "upstream_evidence_ids_sha256",
    }
    if set(lock) != required_keys:
        raise FineMappingError("corrected script-145 lock schema keys drifted")
    state = lock.get("state")
    if (
        lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("artifact_role") != PRE_FAMILY_ROLE
        or state not in {READY_STATE, ZERO_STATE}
        or lock.get("results_accessed_before_lock") is not False
        or lock.get("pair_order") != list(PAIR_ORDER)
        or not identity_matches(lock.get("policy"), policy_identity)
    ):
        raise FineMappingError("corrected script-145 family identity or timing drifted")
    _validate_selection_contract(lock.get("selection_contract"))
    for label in ("policy", "official_blocks"):
        record = _identity_record(lock[label], label)
        if not identity_matches(stable_identity(root, str(record["path"])), record):
            raise FineMappingError(f"corrected script-145 bound {label} drifted")
    for collection_name in ("upstream_family_provenance", "generator_scripts"):
        collection = lock.get(collection_name)
        if not isinstance(collection, Mapping) or not collection:
            raise FineMappingError(f"script-145 lock lacks {collection_name}")
        for label, raw_record in collection.items():
            record = _identity_record(raw_record, f"{collection_name}:{label}")
            if not identity_matches(stable_identity(root, str(record["path"])), record):
                raise FineMappingError(f"script-145 bound input drifted: {record['path']}")
    manifest_record = _identity_record(lock.get("manifest"), "pre-materialization manifest")
    if Path(str(manifest_record["path"])) != Path(manifest_rel):
        raise FineMappingError("script-145 lock points at a different DEEP_LOCUS_MANIFEST")
    fields, rows, manifest_identity = read_tsv(
        root, manifest_rel, PRE_FAMILY_FIELDS,
        allow_header_only=(state == ZERO_STATE or int(lock.get("locus_count", -1)) == 0),
    )
    del fields
    if not identity_matches(manifest_identity, manifest_record):
        raise FineMappingError("DEEP_LOCUS_MANIFEST differs from its corrected script-145 lock")
    if lock.get("locus_count") != len(rows):
        raise FineMappingError("script-145 locus count differs from DEEP_LOCUS_MANIFEST")
    ids = [row["locus_entry_id"] for row in rows]
    if lock.get("locus_entry_ids_in_order") != ids or len(ids) != len(set(ids)):
        raise FineMappingError("script-145 ordered locus family is duplicated or drifted")
    pair_rank = {pair: index for index, pair in enumerate(PAIR_ORDER)}
    block_keys: set[tuple[str, int, int, int]] = set()
    all_evidence: list[str] = []
    previous_sort: tuple[int, int, int, str] | None = None
    for row in rows:
        pair_id = row["pair_id"]
        if pair_id not in PAIR_IDENTITIES:
            raise FineMappingError("DEEP_LOCUS_MANIFEST contains an unknown pair")
        trait1, trait2, role = PAIR_IDENTITIES[pair_id]
        chromosome = _integer(row["CHR"], "locus chromosome", minimum=1)
        start = _integer(row["START"], "locus start", minimum=1)
        stop = _integer(row["STOP"], "locus stop", minimum=start)
        if (
            chromosome > 22 or row["analysis_id"] != ANALYSIS_ID
            or row["family_role"] != role or row["trait1"] != trait1 or row["trait2"] != trait2
            or not SAFE_ID.fullmatch(row["locus_entry_id"])
            or not re.fullmatch(r"LOC[0-9]+", row["ld_block_id"])
            or row["partition_role"] != "OFFICIAL_LAVA_BLOCK_COORDINATES_AS_WORK_PARTITION_ONLY"
            or (
                row["priority_tier"] not in {"PRIMARY_TIER", "SECONDARY_TIER"}
                if pair_id != "CONTROL" else row["priority_tier"] != "CONTROL_TIER"
            )
        ):
            raise FineMappingError("DEEP_LOCUS_MANIFEST row identity or partition drifted")
        block_key = (pair_id, chromosome, start, stop)
        if block_key in block_keys:
            raise FineMappingError("DEEP_LOCUS_MANIFEST contains duplicate pair x official block")
        block_keys.add(block_key)
        local_ids = _split_tokens(row["local_evidence_ids"], "local evidence IDs")
        local_leads = _split_tokens(row["local_lead_variants"], "local lead variants")
        pleio_lead_ids = _split_tokens(
            row["pleiotropy_lead_evidence_ids"], "pleiotropy lead evidence IDs",
        )
        pleio_leads = _split_tokens(row["pleiotropy_lead_variants"], "pleiotropy lead variants")
        pleio_ids = _split_tokens(
            row["pleiotropy_clump_evidence_ids"], "pleiotropy clump evidence IDs",
        )
        if local_ids or local_leads:
            raise FineMappingError("failed-QC LAVA/local evidence must remain literal NA")
        if len(local_ids) != len(local_leads) or len(pleio_lead_ids) != len(pleio_leads):
            raise FineMappingError("lead evidence IDs and retained lead lists are not parallel")
        if not set(pleio_lead_ids).issubset(pleio_ids):
            raise FineMappingError("lead evidence is absent from full clump-evidence union")
        if not local_ids and not pleio_ids:
            raise FineMappingError("nonzero locus row has no retained evidence")
        for evidence_id in [*local_ids, *pleio_ids]:
            if f":{pair_id}:" not in evidence_id:
                raise FineMappingError("evidence ID does not encode its frozen pair")
        for lead in [*local_leads, *pleio_leads]:
            if not RSID.fullmatch(lead):
                raise FineMappingError("retained evidence lead is not an rsID")
        all_evidence.extend([*local_ids, *pleio_ids])
        sort_key = (pair_rank[pair_id], chromosome, start, row["locus_entry_id"])
        if previous_sort is not None and sort_key <= previous_sort:
            raise FineMappingError("DEEP_LOCUS_MANIFEST is not in deterministic pair/block order")
        previous_sort = sort_key
    if len(all_evidence) != len(set(all_evidence)):
        raise FineMappingError("evidence was duplicated or cross-pair-collapsed")
    if lock.get("evidence_count") != len(all_evidence) or lock.get("evidence_ids_sha256") != digest_json(all_evidence):
        raise FineMappingError("complete evidence-family accounting drifted")
    unavailable_record = lock.get("unavailable_manifest")
    unavailable_count = lock.get("unavailable_count")
    if not isinstance(unavailable_count, int) or isinstance(unavailable_count, bool) or unavailable_count < 0:
        raise FineMappingError("script-145 unavailable count is invalid")
    unavailable_rows: list[dict[str, str]] = []
    unavailable_identity: dict[str, object] | None = None
    if unavailable_count:
        record = _identity_record(unavailable_record, "script-145 unavailable manifest")
        if Path(str(record["path"])) != PRE_UNAVAILABLE_REL:
            raise FineMappingError("script-145 unavailable manifest path drifted")
        _, unavailable_rows, unavailable_identity = read_tsv(
            root, PRE_UNAVAILABLE_REL, PRE_UNAVAILABLE_FIELDS,
        )
        if not identity_matches(unavailable_identity, record) or len(unavailable_rows) != unavailable_count:
            raise FineMappingError("script-145 unavailable manifest identity/count drifted")
    elif unavailable_record is not None:
        raise FineMappingError("zero unavailable family must use a null manifest identity")
    unavailable_ids = [row["unavailable_entry_id"] for row in unavailable_rows]
    if lock.get("unavailable_entry_ids_in_order") != unavailable_ids or len(unavailable_ids) != len(set(unavailable_ids)):
        raise FineMappingError("unavailable entry identities are duplicated or drifted")
    unavailable_evidence: list[str] = []
    previous_unavailable: tuple[int, int, int, str] | None = None
    for row in unavailable_rows:
        pair_id = row["pair_id"]
        if pair_id not in PAIR_IDENTITIES:
            raise FineMappingError("unavailable row has an unknown pair")
        trait1, trait2, role = PAIR_IDENTITIES[pair_id]
        chromosome = _integer(row["CHR"], "unavailable CHR", minimum=1)
        bp = _integer(row["BP"], "unavailable BP", minimum=1)
        evidence = _split_tokens(row["evidence_ids"], "unavailable evidence IDs")
        if (
            chromosome > 22 or row["analysis_id"] != ANALYSIS_ID
            or row["family_role"] != role or row["trait1"] != trait1 or row["trait2"] != trait2
            or not SAFE_ID.fullmatch(row["unavailable_entry_id"])
            or not RSID.fullmatch(row["lead_variant"])
            or row["ld_block_id"] != "UNMAPPED"
            or row["terminal_status"] != "BLOCKED_BY_DATA"
            or row["blocker_reason"] != "NO_PREDECLARED_SIGNED_LD_BLOCK"
            or not evidence
            or any(f":{pair_id}:" not in evidence_id for evidence_id in evidence)
        ):
            raise FineMappingError("unavailable no-block row identity/semantics drifted")
        sort_key = (pair_rank[pair_id], chromosome, bp, row["unavailable_entry_id"])
        if previous_unavailable is not None and sort_key <= previous_unavailable:
            raise FineMappingError("unavailable family is not in deterministic pair/position order")
        previous_unavailable = sort_key
        unavailable_evidence.extend(evidence)
    if len(unavailable_evidence) != len(set(unavailable_evidence)):
        raise FineMappingError("unavailable evidence IDs are duplicated")
    if set(all_evidence).intersection(unavailable_evidence):
        raise FineMappingError("one evidence ID appears in both mapped and unavailable families")
    if (
        lock.get("unavailable_evidence_count") != len(unavailable_evidence)
        or lock.get("unavailable_evidence_ids_sha256") != digest_json(unavailable_evidence)
        or lock.get("upstream_evidence_count") != len(all_evidence) + len(unavailable_evidence)
        or lock.get("upstream_evidence_ids_sha256")
        != digest_json([*all_evidence, *unavailable_evidence])
    ):
        raise FineMappingError("mapped/unavailable upstream evidence denominator drifted")
    zero_record = lock.get("zero_family_provenance")
    if state == ZERO_STATE:
        if rows or unavailable_rows or zero_record is None:
            raise FineMappingError("zero family lacks exclusive zero-family provenance")
        record = _identity_record(zero_record, "script-145 zero-family provenance")
        if not identity_matches(stable_identity(root, str(record["path"])), record):
            raise FineMappingError("script-145 zero-family provenance drifted")
    elif (not rows and not unavailable_rows) or zero_record is not None:
        raise FineMappingError("nonzero script-145 family is empty or has zero-family provenance")
    expected_gate = {
        "state": state,
        "rows": rows,
        "unavailable_rows": unavailable_rows,
        "family_path": str(Path(manifest_rel)),
        "family_lock_path": str(Path(lock_rel)),
        "family_identity": manifest_identity,
        "family_lock_identity": lock_identity,
        "unavailable_path": str(PRE_UNAVAILABLE_REL) if unavailable_count else None,
        "unavailable_lock_path": None,
        "unavailable_identity": unavailable_identity,
        "unavailable_lock_identity": None,
        "zero_provenance_path": (
            str(_identity_record(zero_record, "script-145 zero-family provenance")["path"])
            if state == ZERO_STATE else None
        ),
        "zero_provenance_identity": (
            dict(_identity_record(zero_record, "script-145 zero-family provenance"))
            if state == ZERO_STATE else None
        ),
        "mapped_count": len(rows),
        "unavailable_count": len(unavailable_rows),
        "upstream_evidence_denominator": len(all_evidence) + len(unavailable_evidence),
        "family_sha256": manifest_identity["sha256"],
        "unavailable_sha256": (
            unavailable_identity["sha256"] if unavailable_identity is not None else None
        ),
    }
    if gate_result != expected_gate:
        raise FineMappingError(
            "script-145 authoritative verification differs from its canonical adapter artifacts"
        )
    return {
        "state": state, "rows": rows, "lock": lock, "lock_identity": lock_identity,
        "manifest_identity": manifest_identity, "policy": policy,
        "policy_identity": policy_identity, "evidence_ids": all_evidence,
        "unavailable_rows": unavailable_rows, "unavailable_identity": unavailable_identity,
        "unavailable_evidence_ids": unavailable_evidence,
    }


def per_snp_n_diagnostics(values: Sequence[float]) -> dict[str, float | int | str]:
    """Quantify the complete N vector with a fixed, result-blind type-7 summary."""

    if not values:
        raise FineMappingError("cannot choose scalar N from an empty locus")
    parsed = [_number(value, "per-SNP N") for value in values]
    if any(value <= 0 for value in parsed):
        raise FineMappingError("per-SNP N must be positive")
    ordered = sorted(parsed)

    def type7(probability: float) -> float:
        position = (len(ordered) - 1) * probability
        lower = math.floor(position)
        upper = math.ceil(position)
        if lower == upper:
            return float(ordered[lower])
        weight = position - lower
        return float(ordered[lower] + weight * (ordered[upper] - ordered[lower]))

    minimum, maximum = float(ordered[0]), float(ordered[-1])
    median = type7(0.5)
    constant = all(value == parsed[0] for value in parsed[1:])
    return {
        "count": len(parsed), "min": minimum, "q1": type7(0.25),
        "median": median, "q3": type7(0.75), "max": maximum,
        "max_to_min_ratio": maximum / minimum,
        "fraction_below_90pct_max": sum(value < 0.90 * maximum for value in parsed) / len(parsed),
        "fraction_below_50pct_max": sum(value < 0.50 * maximum for value in parsed) / len(parsed),
        "dispersion_status": N_DISPERSION_CONSTANT if constant else N_DISPERSION_VARIABLE,
    }


def scalar_n(values: Sequence[float]) -> tuple[float, str]:
    """Apply the additive, result-blind scalar-N/claim-cap convention."""

    diagnostics = per_snp_n_diagnostics(values)
    if diagnostics["dispersion_status"] == N_DISPERSION_CONSTANT:
        return float(diagnostics["median"]), SCALAR_N_RULE_CONSTANT
    return float(diagnostics["median"]), SCALAR_N_RULE_MEDIAN


def per_snp_n_sha256(values: Sequence[float]) -> str:
    return hashlib.sha256(
        ("\n".join(format(float(value), ".17g") for value in values) + "\n").encode()
    ).hexdigest()


def oriented_alleles(
    source_a1: str, source_a2: str, source_eaf: float,
    reference_a1: str, reference_a2: str,
) -> tuple[int, float, str]:
    """Orient exact or unambiguous strand-equivalent alleles to reference A1."""

    source = (source_a1.upper(), source_a2.upper())
    reference = (reference_a1.upper(), reference_a2.upper())
    if (
        any(base not in BASES for base in (*source, *reference))
        or source[0] == source[1] or reference[0] == reference[1]
        or frozenset(source) in PALINDROMIC or frozenset(reference) in PALINDROMIC
    ):
        raise FineMappingError("alleles are mismatched, multibase, or palindromic")
    eaf = _probability(source_eaf, "source effect-allele frequency")
    if source == reference:
        return 1, eaf, "DIRECT"
    if source == (reference[1], reference[0]):
        return -1, 1.0 - eaf, "SWAP"
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    strand_equivalent = (complement[source[0]], complement[source[1]])
    if strand_equivalent == reference:
        return 1, eaf, "COMPLEMENT"
    if strand_equivalent == (reference[1], reference[0]):
        return -1, 1.0 - eaf, "COMPLEMENT_SWAP"
    raise FineMappingError("alleles are not an unambiguous strand-equivalent reference match")


def _open_text_table(path: Path) -> Iterator[io.TextIOBase]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def load_reference_block(
    path: Path, *, chromosome: int, start: int, stop: int, minimum_maf: float,
) -> tuple[list[str], dict[str, dict[str, str]], int]:
    required = {"SNP", "CHR", "POS", "A1", "A2", "FREQ", "NOBS"}
    order: list[str] = []
    output: dict[str, dict[str, str]] = {}
    sample_sizes: set[int] = set()
    seen_rsids: set[str] = set()
    try:
        with path.open("rt", encoding="utf-8", newline="") as handle:
            header = handle.readline().split()
            if not required.issubset(header):
                raise FineMappingError("LAVA .info schema lacks required identity fields")
            for line_number, line in enumerate(handle, start=2):
                values = line.split()
                if len(values) != len(header):
                    raise FineMappingError(f"malformed LAVA .info row {line_number}")
                row = dict(zip(header, values))
                row_chr = _integer(row["CHR"], "reference CHR", minimum=1)
                position = _integer(row["POS"], "reference POS", minimum=1)
                if row_chr != chromosome or not start <= position <= stop:
                    continue
                snp = row["SNP"].lower()
                if not RSID.fullmatch(snp):
                    raise FineMappingError("official block reference contains a non-rsID")
                if snp in seen_rsids:
                    raise FineMappingError("official block reference contains a duplicate rsID")
                seen_rsids.add(snp)
                a1, a2 = row["A1"].upper(), row["A2"].upper()
                frequency = _probability(row["FREQ"], f"reference frequency {snp}")
                nobs = _integer(row["NOBS"], f"reference NOBS {snp}", minimum=1)
                sample_sizes.add(nobs)
                # Palindromic and low-frequency reference SNPs are not eligible.
                if (
                    a1 not in BASES or a2 not in BASES or a1 == a2
                    or frozenset((a1, a2)) in PALINDROMIC
                    or min(frequency, 1 - frequency) < minimum_maf
                ):
                    continue
                output[snp] = {
                    "SNP": snp, "CHR": str(chromosome), "BP": str(position),
                    "A1": a1, "A2": a2, "REFERENCE_EAF": format(frequency, ".17g"),
                    "NOBS": str(nobs),
                }
                order.append(snp)
    except FineMappingError:
        raise
    except (OSError, UnicodeError) as error:
        raise FineMappingError(f"cannot read LAVA reference info: {error}") from error
    if not output:
        raise FineMappingError("official block has no eligible non-palindromic reference SNPs")
    if not sample_sizes:
        raise FineMappingError("LAVA reference lacks sample-size evidence")
    # Per-variant NOBS legitimately varies with UKB genotype missingness.  The
    # reference-wide sample size used for the >=10,000 gate is read from the
    # LAVA reference object and recorded by script 149, never inferred here.
    return order, output, min(sample_sizes)


def extract_dense_block(
    path: Path, reference: Mapping[str, Mapping[str, str]], *, chromosome: int,
    minimum_maf: float, info_status: str,
) -> tuple[dict[str, dict[str, str]], dict[str, int]]:
    expected = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
    output: dict[str, dict[str, str]] = {}
    counts: Counter[str] = Counter()
    try:
        with _open_text_table(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or []) != expected:
                raise FineMappingError(f"dense summary schema drifted: {path}")
            for row_index, row in enumerate(reader, start=1):
                if None in row or any(value is None for value in row.values()):
                    raise FineMappingError(f"malformed dense row {row_index}: {path}")
                snp = row["SNP"].lower()
                ref = reference.get(snp)
                if ref is None:
                    continue
                counts["reference_rsid_seen"] += 1
                if snp in output:
                    raise FineMappingError(f"duplicate dense locus rsID: {snp}")
                source_chr = _integer(row["CHR"], f"{snp} CHR", minimum=1)
                source_bp = _integer(row["BP"], f"{snp} BP", minimum=1)
                if source_chr != chromosome or source_bp != int(ref["BP"]):
                    raise FineMappingError(f"exact rsID/CHR/BP mismatch for {snp}")
                frequency = _probability(row["FRQ"], f"{snp} FRQ")
                if min(frequency, 1 - frequency) < minimum_maf:
                    raise FineMappingError(f"dense source violates MAF>=0.01 for {snp}")
                sign, eaf, orientation = oriented_alleles(
                    row["A1"], row["A2"], frequency, str(ref["A1"]), str(ref["A2"]),
                )
                beta = _number(row["BETA"], f"{snp} BETA") * sign
                se = _number(row["SE"], f"{snp} SE")
                p_value = _probability(row["P"], f"{snp} P")
                n_value = _number(row["N"], f"{snp} N")
                if se <= 0 or n_value <= 0:
                    raise FineMappingError(f"SE and N must be positive for {snp}")
                output[snp] = {
                    "SNP": snp, "CHR": str(chromosome), "BP": str(ref["BP"]),
                    "A1": str(ref["A1"]), "A2": str(ref["A2"]),
                    "BETA": format(beta, ".17g"), "SE": format(se, ".17g"),
                    "P": format(p_value, ".17g"), "EAF": format(eaf, ".17g"),
                    # The frozen dense files do not retain row-wise INFO.  NA is
                    # honest for both a source-level absence and an upstream
                    # row-wise filter whose individual values were not retained.
                    "INFO": "NA", "N": format(n_value, ".17g"),
                }
                counts["eligible"] += 1
                if orientation in {"SWAP", "COMPLEMENT_SWAP"}:
                    counts["allele_swaps"] += 1
                if orientation in {"COMPLEMENT", "COMPLEMENT_SWAP"}:
                    counts["allele_complements"] += 1
                counts[f"allele_orientation_{orientation.lower()}"] += 1
    except FineMappingError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        raise FineMappingError(f"cannot stream dense summary {path}: {error}") from error
    counts["info_source_level_absent"] = int("ABSENT_IN_RELEASE_SOURCE_LEVEL_ONLY" in info_status)
    counts["info_rowwise_values_not_retained"] = int("ROW_WISE_FILTERED" in info_status)
    return output, dict(counts)


def _panel_metadata(panel_row: Mapping[str, str]) -> dict[str, str]:
    kind = panel_row.get("type")
    if kind == "binary":
        ncase = _integer(panel_row.get("ncase"), "panel ncase", minimum=1)
        ncontrol = _integer(panel_row.get("ncontrol"), "panel ncontrol", minimum=1)
        fraction = ncase / (ncase + ncontrol)
        return {
            "type": "cc", "ncase": str(ncase), "ncontrol": str(ncontrol),
            "case_fraction": format(fraction, ".17g"), "sdY": "NA",
        }
    if kind == "continuous":
        _integer(panel_row.get("n_total"), "panel quantitative N", minimum=1)
        return {
            "type": "quant", "ncase": "NA", "ncontrol": "NA",
            "case_fraction": "NA", "sdY": "NA",
        }
    raise FineMappingError("trait type is neither frozen binary nor continuous")


def _load_dense_registry(root: Path, family: Mapping[str, Any]) -> dict[str, dict[str, str]]:
    _, rows, identity = read_tsv(root, DENSE_QC_REL)
    dense_lock, _ = read_json(root, DENSE_LOCK_REL)
    if dense_lock.get("dense_qc_sha256") != identity["sha256"]:
        raise FineMappingError("dense QC differs from its frozen lock")
    registry: dict[str, dict[str, str]] = {}
    for row in rows:
        trait = row.get("trait_id", "")
        if trait in registry:
            raise FineMappingError("dense QC contains duplicate trait rows")
        if (
            row.get("analysis_build") != "hg19" or row.get("ancestry") != "EUR"
            or row.get("observed_schema") != "SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N"
            or not row.get("fine_mapping_readiness", "").startswith("PASS_DENSE")
        ):
            raise FineMappingError(f"dense QC is not fine-mapping-ready for {trait}")
        registry[trait] = row
    required = {trait for pair in PAIR_IDENTITIES.values() for trait in pair[:2]}
    if set(registry) != required:
        raise FineMappingError("dense QC is not the exact five-trait family")
    del family
    return registry


def _load_panel(root: Path) -> dict[str, dict[str, str]]:
    _, rows, _ = read_tsv(root, PANEL_REL)
    output = {row["trait_id"]: row for row in rows}
    required = {trait for pair in PAIR_IDENTITIES.values() for trait in pair[:2]}
    if not required.issubset(output):
        raise FineMappingError("analysis panel lacks Track B traits")
    return output


def _find_reference_record(reference_state: Mapping[str, Any], path: str) -> Mapping[str, Any]:
    provenance = reference_state.get("reference_provenance")
    if not isinstance(provenance, Mapping):
        raise FineMappingError("reference-state API lacks immutable provenance")
    records = provenance.get("extracted_files")
    if not isinstance(records, list):
        raise FineMappingError("reference provenance lacks extracted payload family")
    matches = [record for record in records if isinstance(record, Mapping) and record.get("path") == path]
    if len(matches) != 1:
        raise FineMappingError(f"reference provenance lacks exact payload: {path}")
    return matches[0]


def estimate_peak_rss_bytes(variant_count: int, policy: Mapping[str, Any]) -> int:
    estimate = policy["ram_aware_execution_contract"]["admission_estimate"]
    count = _integer(variant_count, "variant count", minimum=1)
    calculated = (
        int(estimate["ld_scalar_bytes"]) * count * count
        * int(estimate["ld_live_copy_multiplier"])
        + int(estimate["fixed_process_overhead_bytes"])
    )
    return max(int(estimate["absolute_minimum_memory_bytes"]), calculated)


def _memory_available_bytes() -> int:
    # Linux MemAvailable is the most useful live admission metric.
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        for line in meminfo.read_text(encoding="ascii").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    # macOS: free + inactive + speculative pages are reclaimable.  This is a
    # live available-memory estimate, not physical RAM.
    try:
        page_size = int(subprocess.check_output(["sysctl", "-n", "hw.pagesize"], text=True).strip())
        vm_stat = subprocess.check_output(["vm_stat"], text=True)
        pages: dict[str, int] = {}
        for line in vm_stat.splitlines():
            match = re.match(r"([^:]+):\s+([0-9]+)\.", line)
            if match:
                pages[match.group(1)] = int(match.group(2))
        return page_size * sum(
            pages.get(key, 0) for key in (
                "Pages free", "Pages inactive", "Pages speculative", "Pages purgeable",
            )
        )
    except (OSError, subprocess.SubprocessError, ValueError):
        pass
    page_size = int(os.sysconf("SC_PAGE_SIZE"))
    available_pages = int(os.sysconf("SC_AVPHYS_PAGES"))
    return page_size * available_pages


def host_physical_memory_bytes() -> int:
    try:
        if Path("/proc/meminfo").is_file():
            for line in Path("/proc/meminfo").read_text(encoding="ascii").splitlines():
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) * 1024
        return int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True).strip())
    except (OSError, subprocess.SubprocessError, ValueError):
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))


def memory_admission(
    variant_count: int, policy: Mapping[str, Any], *, available_bytes: int | None = None,
    configured_envelope_bytes: int = BASE_PHYSICAL_ENVELOPE_BYTES,
    larger_host_continuation: bool = False,
) -> dict[str, int | str]:
    estimated = estimate_peak_rss_bytes(variant_count, policy)
    reserve = int(
        policy["ram_aware_execution_contract"]["admission_estimate"]
        ["minimum_unallocated_host_reserve_bytes"]
    )
    available = _memory_available_bytes() if available_bytes is None else int(available_bytes)
    physical = host_physical_memory_bytes()
    configured = int(configured_envelope_bytes)
    if configured < BASE_PHYSICAL_ENVELOPE_BYTES:
        raise FineMappingError("configured materialization envelope is below 8 GiB")
    if configured > BASE_PHYSICAL_ENVELOPE_BYTES and not larger_host_continuation:
        raise FineMappingError("larger materialization envelope requires explicit continuation")
    envelope = min(physical, configured)
    status = (
        "ADMITTED_WHOLE_LOCUS"
        if available >= estimated + reserve and estimated + reserve <= envelope
        else "BLOCKED_BY_COMPUTE"
    )
    return {
        "estimated_peak_rss_bytes": estimated,
        "minimum_unallocated_host_reserve_bytes": reserve,
        "available_memory_before_bytes": available,
        "host_physical_memory_bytes": physical,
        "configured_envelope_bytes": configured,
        "physical_envelope_bytes": envelope,
        "admission_status": status,
    }


def _write_new_file(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    descriptor = os.open(path, flags, 0o600)
    try:
        offset = 0
        while offset < len(content):
            offset += os.write(descriptor, content[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _task_row(task_path: Path) -> dict[str, str]:
    try:
        with task_path.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or []) != TASK_FIELDS:
                raise FineMappingError("locus task schema drifted")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise FineMappingError(f"cannot read locus task: {error}") from error
    if len(rows) != 1 or None in rows[0] or any(value is None for value in rows[0].values()):
        raise FineMappingError("locus task must contain exactly one complete row")
    return rows[0]


def _load_ld_qc(path: Path) -> dict[str, str]:
    try:
        with path.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or []) != LD_QC_FIELDS:
                raise FineMappingError("signed-LD QC schema drifted")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise FineMappingError(f"cannot read signed-LD QC: {error}") from error
    if len(rows) != 1:
        raise FineMappingError("signed-LD QC must contain exactly one row")
    return rows[0]


def require_measured_materialization_context() -> int:
    """Require coordinator authorization plus a fresh owned OS session."""

    if os.environ.get("TRACK_B_FINEMAP_MEASURED_PHASE") != "MATERIALIZATION":
        raise FineMappingError("materialization lacks the measured-phase authorization marker")
    if os.getsid(0) != os.getpid():
        raise FineMappingError("materialization is not the leader of a fresh measured OS session")
    lease_text = os.environ.get("TRACK_B_FINEMAP_WORKER_LEASE_FD")
    try:
        lease_fd = int(lease_text or "")
        lease_status = os.fstat(lease_fd)
    except (ValueError, OSError) as error:
        raise FineMappingError("materialization lacks an inherited live worker lease") from error
    if lease_fd <= 2 or not stat.S_ISREG(lease_status.st_mode):
        raise FineMappingError("materialization worker lease is unsafe")
    return lease_fd


def deep_validate_signed_ld(
    ld_path: Path, order_rows: Sequence[Mapping[str, str]], *, policy: Mapping[str, Any],
    scratch_directory: Path,
) -> dict[str, float | int | str]:
    """Full-stream LD validation, including exact order, symmetry and PSD.

    The matrix is staged in a disk-backed NumPy memmap while the gzip stream is
    consumed exactly once.  The eigenvalue calculation is intentionally exact
    for the complete indivisible matrix; no submatrix or approximate thinning
    is permitted.
    """

    require_measured_materialization_context()
    try:
        import numpy as np
    except ImportError as error:  # pragma: no cover - production dependency is bundled
        raise FineMappingError("NumPy is required for full signed-LD validation") from error
    names = [row["SNP"] for row in order_rows]
    if not names or len(names) != len(set(names)):
        raise FineMappingError("variant order is empty or duplicated")
    n = len(names)
    scratch = scratch_directory
    try:
        scratch_status = os.lstat(scratch)
    except OSError as error:
        raise FineMappingError("LD-validation scratch directory is absent") from error
    if not stat.S_ISDIR(scratch_status.st_mode) or stat.S_ISLNK(scratch_status.st_mode):
        raise FineMappingError("LD-validation scratch is not a real existing directory")
    descriptor, memmap_name = tempfile.mkstemp(prefix=".ld-validate-", suffix=".bin", dir=scratch)
    os.close(descriptor)
    matrix = np.memmap(memmap_name, dtype="float64", mode="w+", shape=(n, n))
    try:
        opener = gzip.open if ld_path.suffix == ".gz" else open
        try:
            with opener(ld_path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if header != ["SNP", *names]:
                    raise FineMappingError("signed-LD columns differ from exact variant order")
                observed_rows = 0
                for index, values in enumerate(reader):
                    if index >= n or len(values) != n + 1 or values[0] != names[index]:
                        raise FineMappingError("signed-LD row identity or dimension drifted")
                    try:
                        row = np.asarray(values[1:], dtype=np.float64)
                    except (TypeError, ValueError) as error:
                        raise FineMappingError("signed-LD contains a nonnumeric value") from error
                    if not np.isfinite(row).all():
                        raise FineMappingError("signed-LD contains a nonfinite value")
                    matrix[index, :] = row
                    observed_rows += 1
                if observed_rows != n:
                    raise FineMappingError("signed-LD matrix is truncated or has extra rows")
        except (OSError, EOFError, UnicodeError, csv.Error) as error:
            raise FineMappingError(f"cannot full-stream signed LD: {error}") from error
        matrix.flush()
        symmetry = float(np.max(np.abs(matrix - matrix.T)))
        diagonal = float(np.max(np.abs(np.diag(matrix) - 1.0)))
        maximum = float(np.max(np.abs(matrix)))
        signed = policy["signed_ld_contract"]
        if symmetry > float(signed["symmetry_tolerance"]):
            raise FineMappingError("signed-LD symmetry exceeds 1e-8")
        if diagonal > float(signed["diagonal_tolerance"]):
            raise FineMappingError("signed-LD diagonal differs from one by more than 1e-6")
        if maximum > 1.0 + float(signed["diagonal_tolerance"]):
            raise FineMappingError("signed-LD correlation lies outside [-1,1]")
        # eigvalsh operates on the complete matrix.  Converting the memmap to an
        # ndarray view does not alter order or select a subset.
        minimum_eigenvalue = float(np.linalg.eigvalsh(np.asarray(matrix))[0])
        if (
            not math.isfinite(minimum_eigenvalue)
            or minimum_eigenvalue < float(signed["minimum_eigenvalue_tolerance"])
        ):
            raise FineMappingError("signed-LD minimum eigenvalue is below -1e-6")
        return {
            "variant_count": n, "symmetry_max_abs": symmetry,
            "diagonal_max_abs": diagonal, "minimum_eigenvalue": minimum_eigenvalue,
            "maximum_abs_correlation": maximum, "matrix_kind": "SIGNED_PEARSON_CORRELATION",
        }
    finally:
        try:
            matrix._mmap.close()  # type: ignore[attr-defined]
        except (AttributeError, OSError):
            pass
        try:
            os.unlink(memmap_name)
        except FileNotFoundError:
            pass


def _validate_summary_file(
    path: Path, order_rows: Sequence[Mapping[str, str]], *, expected_n_hash: str,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    # Summary files are intentionally small enough to retain one locus, not one
    # chromosome/genome-wide dense family, in memory.
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or []) != SUMMARY_FIELDS:
                raise FineMappingError("materialized summary schema drifted")
            rows = list(reader)
    except (OSError, EOFError, UnicodeError, csv.Error) as error:
        raise FineMappingError(f"cannot read materialized summary: {error}") from error
    if len(rows) != len(order_rows):
        raise FineMappingError("summary and variant-order counts differ")
    n_values: list[float] = []
    for observed, expected in zip(rows, order_rows, strict=True):
        if any(observed[field] != expected[field] for field in ORDER_FIELDS):
            raise FineMappingError("summary SNP/order/coordinate/alleles drifted")
        if frozenset((observed["A1"], observed["A2"])) in PALINDROMIC:
            raise FineMappingError("materialized summary contains a palindromic SNP")
        beta = _number(observed["BETA"], "materialized BETA")
        se = _number(observed["SE"], "materialized SE")
        p_value = _probability(observed["P"], "materialized P")
        eaf = _probability(observed["EAF"], "materialized EAF")
        n_value = _number(observed["N"], "materialized N")
        del beta, p_value
        if se <= 0 or n_value <= 0 or min(eaf, 1 - eaf) < 0.01:
            raise FineMappingError("materialized SE/N/MAF constraint failed")
        if observed["INFO"] != "NA":
            info = _probability(observed["INFO"], "materialized INFO")
            if info <= 0:
                raise FineMappingError("materialized INFO must be positive when present")
        n_values.append(n_value)
    if per_snp_n_sha256(n_values) != expected_n_hash:
        raise FineMappingError("per-SNP N vector differs from its task hash")
    # The stable identity is calculated after semantic parsing; task validation
    # rechecks it against the hash frozen before any result access.
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    identity = {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}
    return rows, identity


def _validate_ld_runtime_attestation(
    root: Path, task: Mapping[str, str],
) -> dict[str, Any]:
    """Verify the extractor's exact pre/post loaded LAVA runtime attestation."""

    package_lock = verify_ld_runtime_package_lock(
        root, task["ld_runtime_package_lock_path"],
        expected_sha256=task["ld_runtime_package_lock_sha256"],
    )
    _, rows, identity = read_tsv(
        root, task["ld_runtime_attestation_path"], LD_RUNTIME_ATTESTATION_FIELDS,
    )
    expected = {
        "schema_version": LD_RUNTIME_ATTESTATION_SCHEMA,
        "analysis_id": ANALYSIS_ID,
        "locus_entry_id": task["locus_entry_id"],
        "reference_id": task["reference_id"],
        "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
        "runtime_package_lock_sha256": task["ld_runtime_package_lock_sha256"],
        "precheck_status": "PASS_EXACT_LOADED_CLOSURE_PATHS_AND_BYTES",
        "postcheck_status": "PASS_EXACT_LOADED_CLOSURE_PATHS_AND_BYTES",
        "loaded_namespace_count": str(len(PINNED_LD_RUNTIME_NAMESPACE_CLOSURE)),
        "loaded_namespace_closure": ",".join(PINNED_LD_RUNTIME_NAMESPACE_CLOSURE),
        "lava_version": "0.1.5",
    }
    if (
        len(rows) != 1 or rows[0] != expected
        or identity["sha256"] != task["ld_runtime_attestation_sha256"]
    ):
        raise FineMappingError("signed-LD extractor runtime attestation drifted")
    return {"row": rows[0], "identity": identity, "package_lock": package_lock}


def _validate_ld_attestation(
    root: Path, task: Mapping[str, str], ld_qc: Mapping[str, str],
    *, policy: Mapping[str, Any], variant_count: int,
) -> tuple[dict[str, Any], dict[str, object]]:
    attestation, identity = read_json(root, task["ld_validation_path"])
    required = {
        "schema_version", "analysis_id", "locus_entry_id", "validation_phase",
        "results_accessed", "algorithm", "signed_ld_sha256", "variant_order_sha256",
        "ld_qc_sha256", "reference_id", "reference_sample_size", "metrics",
        "tolerances", "validator", "execution_amendment_sha256",
        "ld_runtime_package_lock_sha256", "ld_runtime_attestation_sha256",
    }
    expected_tolerances = {
        "symmetry_tolerance": policy["signed_ld_contract"]["symmetry_tolerance"],
        "diagonal_tolerance": policy["signed_ld_contract"]["diagonal_tolerance"],
        "minimum_eigenvalue_tolerance": policy["signed_ld_contract"]["minimum_eigenvalue_tolerance"],
    }
    current_validator = stable_identity(root, SCRIPT)
    if (
        set(attestation) != required
        or identity["sha256"] != task["ld_validation_sha256"]
        or attestation.get("schema_version") != LD_VALIDATION_SCHEMA
        or attestation.get("analysis_id") != ANALYSIS_ID
        or attestation.get("locus_entry_id") != task["locus_entry_id"]
        or attestation.get("validation_phase") != "MEASURED_FRESH_MATERIALIZATION_PROCESS"
        or attestation.get("results_accessed") is not False
        or attestation.get("algorithm") != "FULL_MATRIX_EXACT_NUMPY_EIGVALSH_NO_SUBSETTING"
        or attestation.get("signed_ld_sha256") != task["signed_ld_sha256"]
        or attestation.get("variant_order_sha256") != task["variant_order_sha256"]
        or attestation.get("ld_qc_sha256") != task["ld_qc_sha256"]
        or attestation.get("reference_id") != task["reference_id"]
        or attestation.get("reference_sample_size") != int(task["reference_sample_size"])
        or attestation.get("tolerances") != expected_tolerances
        or attestation.get("validator") != current_validator
        or attestation.get("execution_amendment_sha256") != EXECUTION_AMENDMENT_SHA256
        or attestation.get("ld_runtime_package_lock_sha256")
        != task["ld_runtime_package_lock_sha256"]
        or attestation.get("ld_runtime_attestation_sha256")
        != task["ld_runtime_attestation_sha256"]
    ):
        raise FineMappingError("immutable measured signed-LD validation attestation drifted")
    metrics = attestation.get("metrics")
    if not isinstance(metrics, Mapping) or set(metrics) != {
        "variant_count", "symmetry_max_abs", "diagonal_max_abs",
        "minimum_eigenvalue", "maximum_abs_correlation", "matrix_kind",
    }:
        raise FineMappingError("signed-LD validation attestation metrics are malformed")
    if (
        metrics.get("matrix_kind") != "SIGNED_PEARSON_CORRELATION"
        or int(metrics.get("variant_count", -1)) != variant_count
        or _number(metrics.get("symmetry_max_abs"), "attested LD symmetry")
        > float(expected_tolerances["symmetry_tolerance"])
        or _number(metrics.get("diagonal_max_abs"), "attested LD diagonal")
        > float(expected_tolerances["diagonal_tolerance"])
        or _number(metrics.get("minimum_eigenvalue"), "attested LD minimum eigenvalue")
        < float(expected_tolerances["minimum_eigenvalue_tolerance"])
        or _number(metrics.get("maximum_abs_correlation"), "attested LD maximum correlation")
        > 1.0 + float(expected_tolerances["diagonal_tolerance"])
    ):
        raise FineMappingError("attested signed-LD metrics violate the locked QC thresholds")
    comparisons = (
        ("symmetry_max_abs", "symmetry_max_abs", 1e-12),
        ("diagonal_max_abs", "diagonal_max_abs", 1e-12),
        ("minimum_eigenvalue", "minimum_eigenvalue", 1e-8),
        ("maximum_abs_correlation", "maximum_abs_correlation", 1e-12),
    )
    if any(
        not math.isclose(
            _number(ld_qc[qc_field], f"LD-QC {qc_field}"),
            _number(metrics[metric_field], f"attested {metric_field}"),
            rel_tol=0, abs_tol=tolerance,
        )
        for metric_field, qc_field, tolerance in comparisons
    ):
        raise FineMappingError("R extraction QC differs from measured signed-LD attestation")
    return dict(attestation), identity


def deep_validate_locus_task(
    root: Path, task_rel: Path | str, *, expected_manifest_row: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Deeply revalidate an immutable materialized task without side effects."""

    policy, policy_identity = load_policy(root)
    task_identity = stable_identity(root, task_rel)
    task_path = safe_path(root, task_rel, "locus task", must_exist=True)
    task = _task_row(task_path)
    chromosome = _integer(task["CHR"], "task chromosome", minimum=1)
    start = _integer(task["START"], "task start", minimum=1)
    stop = _integer(task["STOP"], "task stop", minimum=start)
    if (
        task["analysis_id"] != ANALYSIS_ID
        or chromosome > 22
        or not SAFE_ID.fullmatch(task["locus_entry_id"])
        or not re.fullmatch(r"LOC[0-9]+", task["ld_block_id"])
        or task["pair_id"] not in PAIR_IDENTITIES
        or task["prior_method"] != "FLAT"
        or task["maximum_causal_signals"] != "10"
        or _number(task["credible_set_coverage"], "coverage") != 0.95
        or _number(task["minimum_absolute_correlation"], "min_abs_corr") != 0.5
        or task["maximum_iterations"] != "1000"
        or task["estimate_residual_variance"] != "FALSE"
        or _number(task["p1"], "p1") != 1e-4
        or _number(task["p2"], "p2") != 1e-4
        or [_number(value, "p12") for value in task["p12_grid"].split(",")]
        != [1e-6, 5e-6, 1e-5, 5e-5]
        or task["single_signal_fallback"] != "FORBIDDEN"
        or task["maximum_locus_variants"] != "NONE"
        or task["locus_splitting"] != "FORBIDDEN"
        or task["variant_thinning"] != "FORBIDDEN"
        or task["lead_centered_truncation"] != "FORBIDDEN"
        or task["reference_id"] != policy["signed_ld_contract"]["reference_id"]
        or task["reference_state"] not in REFERENCE_STATES
        or task["execution_amendment_sha256"] != EXECUTION_AMENDMENT_SHA256
    ):
        raise FineMappingError("locus task model, priors, or no-subsetting contract drifted")
    trait1, trait2, role = PAIR_IDENTITIES[task["pair_id"]]
    if (task["trait1"], task["trait2"], task["family_role"]) != (trait1, trait2, role):
        raise FineMappingError("locus task pair/trait/family identity drifted")
    variant_count = _integer(task["variant_count"], "task variant count", minimum=50)
    if expected_manifest_row is not None:
        comparison = {
            "analysis_id": task["analysis_id"], "locus_entry_id": task["locus_entry_id"],
            "pair_id": task["pair_id"], "family_role": task["family_role"],
            "trait1": task["trait1"], "trait2": task["trait2"], "CHR": task["CHR"],
            "START": task["START"], "STOP": task["STOP"], "ld_block_id": task["ld_block_id"],
            "variant_count": task["variant_count"], "variant_order_path": task["variant_order_path"],
            "signed_ld_path": task["signed_ld_path"], "trait1_summary_sha256": task["summary1_sha256"],
            "trait2_summary_sha256": task["summary2_sha256"],
            "variant_order_sha256": task["variant_order_sha256"],
            "signed_ld_sha256": task["signed_ld_sha256"], "task_path": str(Path(task_rel)),
        }
        if any(expected_manifest_row.get(key) != value for key, value in comparison.items()):
            raise FineMappingError("locus task differs from materialized manifest row")
    identities: dict[str, dict[str, object]] = {"task": task_identity, "policy": policy_identity}
    paths = {
        "summary1": (task["summary1_path"], task["summary1_sha256"]),
        "summary2": (task["summary2_path"], task["summary2_sha256"]),
        "variant_order": (task["variant_order_path"], task["variant_order_sha256"]),
        "signed_ld": (task["signed_ld_path"], task["signed_ld_sha256"]),
        "ld_qc": (task["ld_qc_path"], task["ld_qc_sha256"]),
        "ld_validation": (task["ld_validation_path"], task["ld_validation_sha256"]),
        "ld_runtime_package_lock": (
            task["ld_runtime_package_lock_path"], task["ld_runtime_package_lock_sha256"],
        ),
        "ld_runtime_attestation": (
            task["ld_runtime_attestation_path"], task["ld_runtime_attestation_sha256"],
        ),
        "runtime_package_lock": (
            task["runtime_package_lock_path"], task["runtime_package_lock_sha256"],
        ),
    }
    for label, (path_value, expected_hash) in paths.items():
        identity = stable_identity(root, path_value)
        if identity["sha256"] != expected_hash:
            raise FineMappingError(f"task-bound {label} hash drifted")
        identities[label] = identity
    _, order_rows, _ = read_tsv(root, task["variant_order_path"], ORDER_FIELDS)
    if len(order_rows) != variant_count:
        raise FineMappingError("task variant count differs from variant-order file")
    if len({row["SNP"] for row in order_rows}) != variant_count:
        raise FineMappingError("variant-order file contains duplicates")
    previous_position: int | None = None
    for row in order_rows:
        position = _integer(row["BP"], "variant-order BP", minimum=start)
        if (
            not RSID.fullmatch(row["SNP"]) or row["SNP"] != row["SNP"].lower()
            or _integer(row["CHR"], "variant-order CHR", minimum=1) != chromosome
            or position > stop or row["A1"] not in BASES or row["A2"] not in BASES
            or row["A1"] == row["A2"]
            or frozenset((row["A1"], row["A2"])) in PALINDROMIC
        ):
            raise FineMappingError("variant-order identity/coordinate/allele constraint failed")
        if previous_position is not None and position < previous_position:
            raise FineMappingError("variant-order positions are not in reference order")
        previous_position = position
    first_rows, first_identity = _validate_summary_file(
        safe_path(root, task["summary1_path"], "trait1 summary", must_exist=True), order_rows,
        expected_n_hash=task["trait1_per_snp_N_sha256"],
    )
    second_rows, second_identity = _validate_summary_file(
        safe_path(root, task["summary2_path"], "trait2 summary", must_exist=True), order_rows,
        expected_n_hash=task["trait2_per_snp_N_sha256"],
    )
    if first_identity["sha256"] != task["summary1_sha256"] or second_identity["sha256"] != task["summary2_sha256"]:
        raise FineMappingError("summary changed during deep validation")
    first_n = [_number(row["N"], "trait1 per-SNP N") for row in first_rows]
    second_n = [_number(row["N"], "trait2 per-SNP N") for row in second_rows]
    for prefix, values in (("trait1", first_n), ("trait2", second_n)):
        expected_scalar, expected_rule = scalar_n(values)
        expected_n_fields = _n_task_fields(prefix, per_snp_n_diagnostics(values))
        if (
            _number(task[f"{prefix}_scalar_N"], f"{prefix} scalar N") != expected_scalar
            or task[f"{prefix}_scalar_N_rule"] != expected_rule
            or any(task[field] != value for field, value in expected_n_fields.items())
        ):
            raise FineMappingError(f"{prefix} full-vector N diagnostics/scalar convention drifted")
        kind = task[f"{prefix}_type"]
        if kind == "quant":
            if any(task[f"{prefix}_{field}"] != "NA" for field in ("ncase", "ncontrol", "case_fraction", "sdY")):
                raise FineMappingError("quantitative task invents case fraction or sdY")
        elif kind == "cc":
            ncase = _integer(task[f"{prefix}_ncase"], f"{prefix} ncase", minimum=1)
            ncontrol = _integer(task[f"{prefix}_ncontrol"], f"{prefix} ncontrol", minimum=1)
            if (
                _number(task[f"{prefix}_case_fraction"], f"{prefix} case fraction")
                != ncase / (ncase + ncontrol) or task[f"{prefix}_sdY"] != "NA"
            ):
                raise FineMappingError("case-control metadata or no-sdY rule drifted")
        else:
            raise FineMappingError("task trait type is invalid")
        trait_id = task[prefix]
        if kind != TRAIT_TYPES.get(trait_id):
            raise FineMappingError(f"{trait_id} task type differs from frozen source metadata")
        if task[f"{prefix}_INFO_status"] not in {
            "SOURCE_LEVEL_INFO_ABSENT;ROW_VALUES_NA_NOT_IMPUTED",
            "ROW_WISE_INFO_FILTER_APPLIED_UPSTREAM;VALUES_NOT_RETAINED;ROW_VALUES_NA",
        }:
            raise FineMappingError(f"{prefix} INFO limitation provenance drifted")
    panel = _load_panel(root)
    for prefix in ("trait1", "trait2"):
        expected_metadata = _panel_metadata(panel[task[prefix]])
        observed_metadata = {
            "type": task[f"{prefix}_type"], "ncase": task[f"{prefix}_ncase"],
            "ncontrol": task[f"{prefix}_ncontrol"],
            "case_fraction": task[f"{prefix}_case_fraction"],
            "sdY": task[f"{prefix}_sdY"],
        }
        if observed_metadata != expected_metadata:
            raise FineMappingError(f"{prefix} metadata differs from frozen analysis panel")
    reference_state, manager_identity = verify_reference_state(root)
    if task["reference_state"] != reference_state["state"]:
        raise FineMappingError("task reference state differs from live fully verified state")
    identities["archive_manager"] = manager_identity
    ld_qc = _load_ld_qc(safe_path(root, task["ld_qc_path"], "LD QC", must_exist=True))
    if (
        ld_qc["status"] != "PASS" or ld_qc["matrix_kind"] != "SIGNED_PEARSON_CORRELATION"
        or ld_qc["reference_id"] != task["reference_id"]
        or _integer(ld_qc["variant_count"], "LD-QC variant count", minimum=1) != variant_count
        or _integer(ld_qc["reference_sample_size"], "reference sample size", minimum=1)
        != _integer(task["reference_sample_size"], "task reference sample size", minimum=1)
        or ld_qc["variant_order_sha256"] != task["variant_order_sha256"]
    ):
        raise FineMappingError("R extraction QC identity differs from the immutable locus task")
    if int(task["reference_sample_size"]) < 10000:
        raise FineMappingError("reference sample size is below 10,000")
    ld_attestation, attestation_identity = _validate_ld_attestation(
        root, task, ld_qc, policy=policy, variant_count=variant_count,
    )
    identities["ld_validation"] = attestation_identity
    ld_runtime_attestation = _validate_ld_runtime_attestation(root, task)
    identities["ld_runtime_package_lock"] = ld_runtime_attestation["package_lock"]["identity"]
    identities["ld_runtime_attestation"] = ld_runtime_attestation["identity"]
    package_lock = verify_runtime_package_lock(
        root, task["runtime_package_lock_path"],
        expected_sha256=task["runtime_package_lock_sha256"],
    )
    identities["runtime_package_lock"] = package_lock["identity"]
    ld_metrics = dict(ld_attestation["metrics"])
    return {
        "task": task, "task_identity": task_identity, "identities": identities,
        "order_rows": order_rows, "summary1_rows": first_rows, "summary2_rows": second_rows,
        "ld_metrics": ld_metrics, "reference_state": reference_state,
        "scalar_n_convention": SCALAR_N_CONVENTION, "policy": policy,
        "execution_amendment": EXECUTION_AMENDMENT,
        "runtime_package_lock": package_lock,
        "ld_runtime_package_lock": ld_runtime_attestation["package_lock"],
        "ld_runtime_attestation": ld_runtime_attestation["row"],
        "ld_validation_attestation": ld_attestation,
    }


def _description_version(path: Path, package: str) -> str:
    try:
        values = [
            line.split(":", 1)[1].strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith("Version:")
        ]
    except (OSError, UnicodeError) as error:
        raise FineMappingError(f"cannot read installed {package} DESCRIPTION") from error
    if len(values) != 1 or not values[0]:
        raise FineMappingError(f"installed {package} has no unique Version field")
    return values[0]


def _runtime_package_tree_identity(root: Path, package: str) -> dict[str, str]:
    package_rel = Path(".r-env/lib/R/library") / package
    package_path = safe_path(root, package_rel, f"installed {package}", must_exist=True)
    before = os.lstat(package_path)
    if not stat.S_ISDIR(before.st_mode) or stat.S_ISLNK(before.st_mode):
        raise FineMappingError(f"installed {package} path is not a real directory")
    relative_files: list[Path] = []
    for candidate in sorted(
        package_path.rglob("*"), key=lambda item: str(item.relative_to(package_path)),
    ):
        observed = os.lstat(candidate)
        if stat.S_ISLNK(observed.st_mode):
            raise FineMappingError(f"installed {package} contains a symbolic link")
        if stat.S_ISDIR(observed.st_mode):
            continue
        if not stat.S_ISREG(observed.st_mode):
            raise FineMappingError(f"installed {package} contains a non-regular file")
        relative_files.append(candidate.relative_to(package_path))
    if not relative_files:
        raise FineMappingError(f"installed {package} package tree is empty")
    tree_digest = hashlib.sha256()
    total_bytes = 0
    for relative_file in relative_files:
        # The managed R runtime deduplicates native libraries with hard links.
        # Hash their bytes through a stable descriptor and then re-enumerate the
        # complete tree; unlike scientific artifacts, link count is not an
        # exclusivity claim for installed runtime files.
        identity = stable_identity(
            root, package_rel / relative_file, allow_empty=True, require_single_link=False,
        )
        total_bytes += int(identity["bytes"])
        tree_digest.update(
            f"{relative_file}\t{identity['bytes']}\t{identity['sha256']}\n".encode("utf-8")
        )
    after = os.lstat(package_path)
    after_files: list[Path] = []
    for candidate in sorted(
        package_path.rglob("*"), key=lambda item: str(item.relative_to(package_path)),
    ):
        observed = os.lstat(candidate)
        if stat.S_ISLNK(observed.st_mode):
            raise FineMappingError(f"installed {package} acquired a symbolic link while hashing")
        if stat.S_ISREG(observed.st_mode):
            after_files.append(candidate.relative_to(package_path))
    if (
        (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino)
        or relative_files != after_files
    ):
        raise FineMappingError(f"installed {package} tree changed while hashing")
    return {
        "record_type": "PACKAGE_TREE", "name": package,
        "version": _description_version(package_path / "DESCRIPTION", package),
        "path": str(package_rel), "file_count": str(len(relative_files)),
        "bytes": str(total_bytes), "sha256": tree_digest.hexdigest(),
    }


def _runtime_package_rows_for_closure(
    root: Path, packages: Sequence[str],
) -> list[dict[str, str]]:
    """Bind Rscript plus one exact, ordered installed-package closure."""

    rscript = stable_identity(root, ".r-env/bin/Rscript")
    rows = [{
        "record_type": "R_EXECUTABLE", "name": "Rscript", "version": "R-4.3.3",
        "path": ".r-env/bin/Rscript", "file_count": "1",
        "bytes": str(rscript["bytes"]), "sha256": str(rscript["sha256"]),
    }]
    rows.extend(
        _runtime_package_tree_identity(root, package)
        for package in packages
    )
    return rows


def _runtime_package_rows(root: Path) -> list[dict[str, str]]:
    """Bind Rscript plus the exact loaded susieR/coloc namespace closure."""

    return _runtime_package_rows_for_closure(root, PINNED_RUNTIME_NAMESPACE_CLOSURE)


def runtime_package_lock_bytes(root: Path) -> bytes:
    return tsv_bytes(RUNTIME_PACKAGE_LOCK_FIELDS, _runtime_package_rows(root))


def ld_runtime_package_lock_bytes(root: Path) -> bytes:
    """Freeze the exact Rscript/LAVA dependency closure before LD extraction."""

    return tsv_bytes(
        RUNTIME_PACKAGE_LOCK_FIELDS,
        _runtime_package_rows_for_closure(root, PINNED_LD_RUNTIME_NAMESPACE_CLOSURE),
    )


def verify_runtime_package_lock(
    root: Path, lock_rel: Path | str, *, expected_sha256: str | None = None,
    packages: Sequence[str] = PINNED_RUNTIME_NAMESPACE_CLOSURE,
    label: str = "SuSiE/coloc runtime",
) -> dict[str, Any]:
    """Side-effect-free exact installed-package path and byte verification."""

    _, observed_rows, identity = read_tsv(root, lock_rel, RUNTIME_PACKAGE_LOCK_FIELDS)
    if expected_sha256 is not None and identity["sha256"] != expected_sha256:
        raise FineMappingError("task-bound runtime package lock hash drifted")
    expected_rows = _runtime_package_rows_for_closure(root, packages)
    if observed_rows != expected_rows:
        raise FineMappingError(
            f"Rscript or loaded {label} namespace-closure paths/bytes differ from runtime lock"
        )
    return {"identity": identity, "rows": observed_rows}


def verify_ld_runtime_package_lock(
    root: Path, lock_rel: Path | str, *, expected_sha256: str | None = None,
) -> dict[str, Any]:
    """Revalidate the exact Rscript/LAVA closure used for signed-LD extraction."""

    return verify_runtime_package_lock(
        root, lock_rel, expected_sha256=expected_sha256,
        packages=PINNED_LD_RUNTIME_NAMESPACE_CLOSURE, label="LAVA extractor runtime",
    )


def _runtime_script_identities(root: Path) -> dict[str, dict[str, object]]:
    relatives = (
        "config/track_b_finemapping_policy.json",
        "config/analysis_panel.tsv",
        "config/fine_mapping_sources.tsv",
        "discovery_extension/config/fine_mapping_method_references.tsv",
        ".r-env/bin/Rscript",
        ".r-env/share/finemapping/source_archives/susieR_0.14.2.tar.gz",
        ".r-env/share/finemapping/source_archives/coloc_v5.2.3.tar.gz",
        "scripts/130_track_b_finemapping_contract.py",
        "scripts/145_build_track_b_finemapping_continuation_gate.py",
        "scripts/146_manage_lava_reference_archives.py",
        "scripts/148_track_b_finemapping_runtime.py",
        "scripts/149_extract_track_b_signed_ld.R",
        "scripts/150_run_track_b_susie_coloc.R",
        "scripts/151_run_track_b_finemapping_sequential.py",
        "discovery_extension/scripts/35_run_susie_coloc.R",
        "scripts/fine_mapping_contract.py",
    )
    return {relative: stable_identity(root, relative) for relative in relatives}


def _publish_materialized_directory(
    root: Path, stage: Path, destination_rel: Path, deferred: Mapping[str, bytes],
) -> Path:
    """Claim a destination exclusively, then seal it last.

    The directory is not considered usable until ``bundle.lock.json`` exists.
    A process death leaves a visible, non-scientific partial that the
    coordinator must reconcile; it can never be mistaken for a complete task.
    """

    destination = safe_path(root, destination_rel, "materialized destination")
    ensure_directory(root, destination_rel.parent)
    try:
        os.mkdir(destination, 0o700)
    except FileExistsError as error:
        raise FineMappingError(
            f"materialized destination already exists (sealed or interrupted): {destination_rel}"
        ) from error
    fsync_directory(destination.parent)
    try:
        for source in sorted(stage.iterdir(), key=lambda item: item.name):
            if not source.is_file() or source.is_symlink():
                raise FineMappingError("materialization stage contains a non-regular entry")
            fsync_file(source)
            target = destination / source.name
            os.link(source, target, follow_symlinks=False)
            fsync_directory(destination)
            source.unlink()
        stage.rmdir()
        for name, content in deferred.items():
            if name == "bundle.lock.json":
                continue
            _write_new_file(destination / name, content)
            fsync_directory(destination)
        if "bundle.lock.json" not in deferred:
            raise FineMappingError("materialization publication lacks its final seal")
        _write_new_file(destination / "bundle.lock.json", deferred["bundle.lock.json"])
        fsync_directory(destination)
        fsync_directory(destination.parent)
        return destination
    except BaseException:
        # Never delete a claimed destination here: its partial state is durable
        # interruption evidence and must be reconciled explicitly by script 151.
        raise


def _info_semantics(status: str) -> str:
    if "ABSENT_IN_RELEASE_SOURCE_LEVEL_ONLY" in status:
        return "SOURCE_LEVEL_INFO_ABSENT;ROW_VALUES_NA_NOT_IMPUTED"
    if "ROW_WISE_FILTERED" in status:
        return "ROW_WISE_INFO_FILTER_APPLIED_UPSTREAM;VALUES_NOT_RETAINED;ROW_VALUES_NA"
    raise FineMappingError("dense INFO provenance is neither recognized source-level nor row-wise")


def _sha256_content(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _n_task_fields(prefix: str, diagnostics: Mapping[str, float | int | str]) -> dict[str, str]:
    return {
        f"{prefix}_per_snp_N_count": str(diagnostics["count"]),
        f"{prefix}_per_snp_N_min": format(float(diagnostics["min"]), ".17g"),
        f"{prefix}_per_snp_N_q1": format(float(diagnostics["q1"]), ".17g"),
        f"{prefix}_per_snp_N_median": format(float(diagnostics["median"]), ".17g"),
        f"{prefix}_per_snp_N_q3": format(float(diagnostics["q3"]), ".17g"),
        f"{prefix}_per_snp_N_max": format(float(diagnostics["max"]), ".17g"),
        f"{prefix}_per_snp_N_max_to_min_ratio": format(
            float(diagnostics["max_to_min_ratio"]), ".17g",
        ),
        f"{prefix}_per_snp_N_fraction_below_90pct_max": format(
            float(diagnostics["fraction_below_90pct_max"]), ".17g",
        ),
        f"{prefix}_per_snp_N_fraction_below_50pct_max": format(
            float(diagnostics["fraction_below_50pct_max"]), ".17g",
        ),
        f"{prefix}_N_dispersion_status": str(diagnostics["dispersion_status"]),
    }


def materialize_one_locus(
    root: Path, locus_entry_id: str, *, manifest_rel: Path | str = PRE_MANIFEST_REL,
    lock_rel: Path | str = PRE_LOCK_REL,
    destination_root_rel: Path | str = INPUT_ROOT_REL,
    available_memory_bytes: int | None = None,
    configured_envelope_bytes: int = BASE_PHYSICAL_ENVELOPE_BYTES,
    larger_host_continuation: bool = False,
) -> dict[str, Any]:
    """Materialize and seal one whole pair x official block, without results."""

    if not SAFE_ID.fullmatch(locus_entry_id):
        raise FineMappingError("unsafe locus_entry_id")
    family = validate_pre_materialization_family(root, manifest_rel, lock_rel)
    if family["state"] == ZERO_STATE:
        raise FineMappingError("cannot materialize a locus from a zero family")
    matches = [row for row in family["rows"] if row["locus_entry_id"] == locus_entry_id]
    if len(matches) != 1:
        raise FineMappingError("locus_entry_id does not select exactly one frozen pair x block")
    locus = matches[0]
    destination_root_rel = Path(destination_root_rel)
    destination_rel = destination_root_rel / locus_entry_id
    destination = safe_path(root, destination_rel, "materialized locus destination")
    if destination.exists() or destination.is_symlink():
        return verify_materialized_bundle(root, destination_rel, family=family)
    require_measured_materialization_context()

    policy = family["policy"]
    reference_state, archive_manager_identity = verify_reference_state(root)
    dense = _load_dense_registry(root, family)
    panel = _load_panel(root)
    chromosome = int(locus["CHR"])
    start, stop = int(locus["START"]), int(locus["STOP"])
    prefix = "ref/lava/ukb_v1.1/lava-ukb-v1.1"
    info_rel = f"{prefix}_chr{chromosome}.info"
    bcor_rel = f"{prefix}_chr{chromosome}.bcor"
    info_expected = _find_reference_record(reference_state, info_rel)
    bcor_expected = _find_reference_record(reference_state, bcor_rel)
    info_identity = stable_identity(root, info_rel)
    bcor_identity = stable_identity(root, bcor_rel)
    if not identity_matches(info_identity, info_expected) or not identity_matches(
        bcor_identity, bcor_expected
    ):
        raise FineMappingError("chromosome LAVA reference differs from verified script-146 state")
    reference_order, reference, _ = load_reference_block(
        safe_path(root, info_rel, "reference info", must_exist=True), chromosome=chromosome,
        start=start, stop=stop, minimum_maf=0.01,
    )

    source_identities: dict[str, dict[str, object]] = {}
    extracted: dict[str, dict[str, dict[str, str]]] = {}
    extraction_counts: dict[str, dict[str, int]] = {}
    for trait in (locus["trait1"], locus["trait2"]):
        registry_row = dense[trait]
        dense_rel = registry_row["dense_file"]
        identity = stable_identity(root, dense_rel)
        if (
            identity["sha256"] != registry_row["dense_file_sha256"]
            or identity["bytes"] != int(registry_row["compressed_size_bytes"])
        ):
            raise FineMappingError(f"dense source differs from frozen QC identity: {trait}")
        source_identities[trait] = identity
        extracted[trait], extraction_counts[trait] = extract_dense_block(
            safe_path(root, dense_rel, f"{trait} dense input", must_exist=True), reference,
            chromosome=chromosome, minimum_maf=0.01, info_status=registry_row["INFO_status"],
        )
        # Recheck the immutable source after full streaming to close the TOCTOU
        # window without accepting a same-path replacement.
        if stable_identity(root, dense_rel) != identity:
            raise FineMappingError(f"dense source changed during materialization: {trait}")

    joint_order = [
        snp for snp in reference_order
        if snp in extracted[locus["trait1"]] and snp in extracted[locus["trait2"]]
    ]
    if len(joint_order) < int(policy["dense_input_contract"]["minimum_locus_variants"]):
        raise FineMappingError(
            f"whole official block has only {len(joint_order)} jointly eligible SNPs (<50)"
        )
    if len(joint_order) != len(set(joint_order)):
        raise FineMappingError("joint whole-block variant family is duplicated")
    admission = memory_admission(
        len(joint_order), policy, available_bytes=available_memory_bytes,
        configured_envelope_bytes=configured_envelope_bytes,
        larger_host_continuation=larger_host_continuation,
    )
    if admission["admission_status"] != "ADMITTED_WHOLE_LOCUS":
        raise ComputeBlocked(
            "BLOCKED_BY_COMPUTE indivisible whole official block cannot fit the configured "
            "physical envelope with the locked 1-GiB reserve"
        )

    parent = ensure_directory(root, destination_root_rel)
    stage = parent / f".{locus_entry_id}.{os.getpid()}.{os.urandom(8).hex()}.building"
    try:
        os.mkdir(stage, 0o700)
    except FileExistsError as error:  # cryptographically improbable, still fail closed
        raise FineMappingError("materialization staging collision") from error
    try:
        trait1_rows = [extracted[locus["trait1"]][snp] for snp in joint_order]
        trait2_rows = [extracted[locus["trait2"]][snp] for snp in joint_order]
        order_rows = [{field: trait1_rows[index][field] for field in ORDER_FIELDS} for index in range(len(joint_order))]
        if any(
            any(trait1_rows[index][field] != trait2_rows[index][field] for field in ORDER_FIELDS)
            for index in range(len(joint_order))
        ):
            raise FineMappingError("trait summaries do not share identical SNP/order/alleles")

        summary1_content = deterministic_gzip(tsv_bytes(SUMMARY_FIELDS, trait1_rows))
        summary2_content = deterministic_gzip(tsv_bytes(SUMMARY_FIELDS, trait2_rows))
        order_content = tsv_bytes(ORDER_FIELDS, order_rows)
        summary1_stage = stage / "trait1.tsv.gz"
        summary2_stage = stage / "trait2.tsv.gz"
        order_stage = stage / "variant_order.tsv"
        _write_new_file(summary1_stage, summary1_content)
        _write_new_file(summary2_stage, summary2_content)
        _write_new_file(order_stage, order_content)

        ld_stage = stage / "signed_ld.tsv.gz"
        ld_qc_stage = stage / "signed_ld.qc.tsv"
        ld_runtime_package_stage = stage / "signed_ld.runtime_packages.tsv"
        ld_runtime_attestation_stage = stage / "signed_ld.runtime_attestation.tsv"
        ld_runtime_package_content = ld_runtime_package_lock_bytes(root)
        ld_runtime_package_sha256 = _sha256_content(ld_runtime_package_content)
        _write_new_file(ld_runtime_package_stage, ld_runtime_package_content)
        extractor = safe_path(
            root, "scripts/149_extract_track_b_signed_ld.R", "signed-LD extractor",
            must_exist=True,
        )
        rscript = safe_path(root, ".r-env/bin/Rscript", "pinned Rscript", must_exist=True)
        command = [
            str(rscript), "--vanilla", str(extractor), str(root / prefix), str(chromosome),
            str(order_stage), str(ld_stage), str(ld_qc_stage), "10000",
            "LAVA_UKB_v1.1_EUR_GRCh37", "SIGNED_PEARSON_FROM_LAVA_READ_LD",
            str(ld_runtime_package_stage), ld_runtime_package_sha256,
            str(ld_runtime_attestation_stage), ANALYSIS_ID, locus_entry_id,
            EXECUTION_AMENDMENT_SHA256,
        ]
        inherited_lease_fds: tuple[int, ...] = ()
        lease_fd_text = os.environ.get("TRACK_B_FINEMAP_WORKER_LEASE_FD")
        if lease_fd_text is not None:
            try:
                lease_fd = int(lease_fd_text)
                lease_status = os.fstat(lease_fd)
            except (ValueError, OSError) as error:
                raise FineMappingError("inherited worker-lease descriptor is invalid") from error
            if lease_fd <= 2 or not stat.S_ISREG(lease_status.st_mode):
                raise FineMappingError("inherited worker-lease descriptor is unsafe")
            inherited_lease_fds = (lease_fd,)
        completed = subprocess.run(
            command, cwd=root, capture_output=True, text=True, check=False,
            start_new_session=False, pass_fds=inherited_lease_fds,
        )
        if (
            completed.returncode != 0 or not ld_stage.is_file() or not ld_qc_stage.is_file()
            or not ld_runtime_attestation_stage.is_file()
        ):
            raise FineMappingError(
                "signed Pearson LD extraction failed: "
                + (completed.stdout + completed.stderr).strip()
            )
        ld_metrics = deep_validate_signed_ld(
            ld_stage, order_rows, policy=policy, scratch_directory=stage,
        )
        ld_qc = _load_ld_qc(ld_qc_stage)
        if (
            ld_qc["status"] != "PASS"
            or ld_qc["matrix_kind"] != "SIGNED_PEARSON_CORRELATION"
            or ld_qc["reference_id"] != policy["signed_ld_contract"]["reference_id"]
            or int(ld_qc["variant_count"]) != len(joint_order)
            or int(ld_qc["reference_sample_size"]) < 10000
            or ld_qc["variant_order_sha256"] != _sha256_content(order_content)
        ):
            raise FineMappingError("signed-LD extractor QC identity drifted")

        ld_sha256 = _file_sha256(ld_stage)
        ld_qc_sha256 = _file_sha256(ld_qc_stage)
        ld_runtime_attestation_sha256 = _file_sha256(ld_runtime_attestation_stage)
        ld_runtime_validation = _validate_ld_runtime_attestation(root, {
            "locus_entry_id": locus_entry_id,
            "reference_id": policy["signed_ld_contract"]["reference_id"],
            "ld_runtime_package_lock_path": str(ld_runtime_package_stage.relative_to(root)),
            "ld_runtime_package_lock_sha256": ld_runtime_package_sha256,
            "ld_runtime_attestation_path": str(ld_runtime_attestation_stage.relative_to(root)),
            "ld_runtime_attestation_sha256": ld_runtime_attestation_sha256,
        })
        ld_validation = {
            "schema_version": LD_VALIDATION_SCHEMA,
            "analysis_id": ANALYSIS_ID,
            "locus_entry_id": locus_entry_id,
            "validation_phase": "MEASURED_FRESH_MATERIALIZATION_PROCESS",
            "results_accessed": False,
            "algorithm": "FULL_MATRIX_EXACT_NUMPY_EIGVALSH_NO_SUBSETTING",
            "signed_ld_sha256": ld_sha256,
            "variant_order_sha256": _sha256_content(order_content),
            "ld_qc_sha256": ld_qc_sha256,
            "reference_id": policy["signed_ld_contract"]["reference_id"],
            "reference_sample_size": int(ld_qc["reference_sample_size"]),
            "metrics": ld_metrics,
            "tolerances": {
                "symmetry_tolerance": policy["signed_ld_contract"]["symmetry_tolerance"],
                "diagonal_tolerance": policy["signed_ld_contract"]["diagonal_tolerance"],
                "minimum_eigenvalue_tolerance": policy["signed_ld_contract"]["minimum_eigenvalue_tolerance"],
            },
            "validator": stable_identity(root, SCRIPT),
            "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
            "ld_runtime_package_lock_sha256": ld_runtime_package_sha256,
            "ld_runtime_attestation_sha256": ld_runtime_attestation_sha256,
        }
        ld_validation_content = json.dumps(
            ld_validation, indent=2, sort_keys=True,
        ).encode() + b"\n"
        ld_validation_stage = stage / "signed_ld.validation.json"
        _write_new_file(ld_validation_stage, ld_validation_content)

        runtime_package_content = runtime_package_lock_bytes(root)
        runtime_package_stage = stage / "runtime_packages.tsv"
        _write_new_file(runtime_package_stage, runtime_package_content)

        final_prefix = destination_rel
        summary1_rel = final_prefix / "trait1.tsv.gz"
        summary2_rel = final_prefix / "trait2.tsv.gz"
        order_rel = final_prefix / "variant_order.tsv"
        ld_rel = final_prefix / "signed_ld.tsv.gz"
        ld_qc_rel = final_prefix / "signed_ld.qc.tsv"
        ld_validation_rel = final_prefix / "signed_ld.validation.json"
        ld_runtime_package_rel = final_prefix / "signed_ld.runtime_packages.tsv"
        ld_runtime_attestation_rel = final_prefix / "signed_ld.runtime_attestation.tsv"
        runtime_package_rel = final_prefix / "runtime_packages.tsv"
        task_rel = final_prefix / "task.tsv"
        first_n = [float(row["N"]) for row in trait1_rows]
        second_n = [float(row["N"]) for row in trait2_rows]
        first_scalar, first_rule = scalar_n(first_n)
        second_scalar, second_rule = scalar_n(second_n)
        first_n_diagnostics = per_snp_n_diagnostics(first_n)
        second_n_diagnostics = per_snp_n_diagnostics(second_n)
        first_meta = _panel_metadata(panel[locus["trait1"]])
        second_meta = _panel_metadata(panel[locus["trait2"]])
        task = {
            "analysis_id": ANALYSIS_ID, "locus_entry_id": locus_entry_id,
            "pair_id": locus["pair_id"], "family_role": locus["family_role"],
            "trait1": locus["trait1"], "trait2": locus["trait2"],
            "CHR": locus["CHR"], "START": locus["START"], "STOP": locus["STOP"],
            "ld_block_id": locus["ld_block_id"], "variant_count": str(len(joint_order)),
            "summary1_path": str(summary1_rel), "summary1_sha256": _sha256_content(summary1_content),
            "summary2_path": str(summary2_rel), "summary2_sha256": _sha256_content(summary2_content),
            "variant_order_path": str(order_rel), "variant_order_sha256": _sha256_content(order_content),
            "signed_ld_path": str(ld_rel), "signed_ld_sha256": ld_sha256,
            "ld_qc_path": str(ld_qc_rel), "ld_qc_sha256": ld_qc_sha256,
            "ld_validation_path": str(ld_validation_rel),
            "ld_validation_sha256": _sha256_content(ld_validation_content),
            "ld_runtime_package_lock_path": str(ld_runtime_package_rel),
            "ld_runtime_package_lock_sha256": ld_runtime_package_sha256,
            "ld_runtime_attestation_path": str(ld_runtime_attestation_rel),
            "ld_runtime_attestation_sha256": ld_runtime_attestation_sha256,
            "runtime_package_lock_path": str(runtime_package_rel),
            "runtime_package_lock_sha256": _sha256_content(runtime_package_content),
            "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
            "trait1_type": first_meta["type"], "trait1_ncase": first_meta["ncase"],
            "trait1_ncontrol": first_meta["ncontrol"],
            "trait1_case_fraction": first_meta["case_fraction"], "trait1_sdY": first_meta["sdY"],
            "trait1_scalar_N": format(first_scalar, ".17g"), "trait1_scalar_N_rule": first_rule,
            "trait1_per_snp_N_sha256": per_snp_n_sha256(first_n),
            **_n_task_fields("trait1", first_n_diagnostics),
            "trait1_INFO_status": _info_semantics(dense[locus["trait1"]]["INFO_status"]),
            "trait2_type": second_meta["type"], "trait2_ncase": second_meta["ncase"],
            "trait2_ncontrol": second_meta["ncontrol"],
            "trait2_case_fraction": second_meta["case_fraction"], "trait2_sdY": second_meta["sdY"],
            "trait2_scalar_N": format(second_scalar, ".17g"), "trait2_scalar_N_rule": second_rule,
            "trait2_per_snp_N_sha256": per_snp_n_sha256(second_n),
            **_n_task_fields("trait2", second_n_diagnostics),
            "trait2_INFO_status": _info_semantics(dense[locus["trait2"]]["INFO_status"]),
            "reference_id": policy["signed_ld_contract"]["reference_id"],
            "reference_state": reference_state["state"],
            "reference_sample_size": ld_qc["reference_sample_size"], "prior_method": "FLAT",
            "maximum_causal_signals": "10", "credible_set_coverage": ".95",
            "minimum_absolute_correlation": ".5", "maximum_iterations": "1000",
            "estimate_residual_variance": "FALSE", "p1": ".0001", "p2": ".0001",
            "p12_grid": ".000001,.000005,.00001,.00005",
            "single_signal_fallback": "FORBIDDEN", "maximum_locus_variants": "NONE",
            "locus_splitting": "FORBIDDEN", "variant_thinning": "FORBIDDEN",
            "lead_centered_truncation": "FORBIDDEN", "claim_limit": locus["claim_limit"],
        }
        task_content = tsv_bytes(TASK_FIELDS, [task])
        manifest_row = {
            "analysis_id": ANALYSIS_ID, "locus_entry_id": locus_entry_id,
            "pair_id": locus["pair_id"], "family_role": locus["family_role"],
            "trait1": locus["trait1"], "trait2": locus["trait2"], "CHR": locus["CHR"],
            "START": locus["START"], "STOP": locus["STOP"], "ld_block_id": locus["ld_block_id"],
            "inclusion_sources": locus["inclusion_sources"],
            "local_evidence_ids": locus["local_evidence_ids"],
            "pleiotropy_evidence_ids": locus["pleiotropy_clump_evidence_ids"],
            "method_availability": locus["method_availability"],
            "trait1_dense_input": dense[locus["trait1"]]["dense_file"],
            "trait2_dense_input": dense[locus["trait2"]]["dense_file"],
            "variant_order_path": str(order_rel), "signed_ld_path": str(ld_rel),
            "variant_count": str(len(joint_order)),
            "trait1_summary_sha256": task["summary1_sha256"],
            "trait2_summary_sha256": task["summary2_sha256"],
            "variant_order_sha256": task["variant_order_sha256"],
            "signed_ld_sha256": task["signed_ld_sha256"],
            "estimated_peak_rss_bytes": str(admission["estimated_peak_rss_bytes"]),
            "memory_admission_status": str(admission["admission_status"]),
            "task_path": str(task_rel), "claim_limit": locus["claim_limit"],
        }
        provenance = {
            "schema_version": "sleep-atlas-track-b-finemapping-materialization.1",
            "analysis_id": ANALYSIS_ID, "artifact_role": "PRE_RESULT_WHOLE_LOCUS_INPUT",
            "results_accessed": False, "locus_entry_id": locus_entry_id,
            "pre_family_lock": family["lock_identity"],
            "pre_family_manifest": family["manifest_identity"],
            "policy": family["policy_identity"], "source_dense_inputs": source_identities,
            "reference_state": reference_state, "archive_manager": archive_manager_identity,
            "reference_info": info_identity, "reference_bcor": bcor_identity,
            "eligible_reference_variant_count": len(reference_order),
            "joint_eligible_variant_count": len(joint_order),
            "extraction_counts": extraction_counts,
            "retained_lead_evidence": {
                "evidence_ids": locus["pleiotropy_lead_evidence_ids"],
                "lead_variants": locus["pleiotropy_lead_variants"],
            },
            "priority_tier": locus["priority_tier"],
            "priority_tier_is_non_filtering": True,
            "info_semantics": {
                locus["trait1"]: task["trait1_INFO_status"],
                locus["trait2"]: task["trait2_INFO_status"],
            },
            "scalar_n_convention": SCALAR_N_CONVENTION,
            "trait_scalar_n": {
                locus["trait1"]: {
                    "value": task["trait1_scalar_N"], "rule": first_rule,
                    "full_vector_sha256": task["trait1_per_snp_N_sha256"],
                    "dispersion": first_n_diagnostics,
                },
                locus["trait2"]: {
                    "value": task["trait2_scalar_N"], "rule": second_rule,
                    "full_vector_sha256": task["trait2_per_snp_N_sha256"],
                    "dispersion": second_n_diagnostics,
                },
            },
            "execution_amendment": EXECUTION_AMENDMENT,
            "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
            "admission": admission, "ld_validation_attestation": ld_validation,
            "ld_runtime_attestation": ld_runtime_validation["row"],
            "ld_runtime_package_lock_sha256": task["ld_runtime_package_lock_sha256"],
            "runtime_package_lock_sha256": task["runtime_package_lock_sha256"],
            "ld_extractor_stdout": completed.stdout.strip(),
            "ld_extractor_stderr": completed.stderr.strip(),
            "script_identities": _runtime_script_identities(root),
            "no_subsetting": {
                "maximum_locus_variants": None, "locus_splitting": False,
                "variant_thinning": False, "lead_centered_truncation": False,
            },
            "manifest_row": manifest_row,
        }
        provenance_content = json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n"
        output_records = {
            "trait1.tsv.gz": {"bytes": len(summary1_content), "sha256": task["summary1_sha256"]},
            "trait2.tsv.gz": {"bytes": len(summary2_content), "sha256": task["summary2_sha256"]},
            "variant_order.tsv": {"bytes": len(order_content), "sha256": task["variant_order_sha256"]},
            "signed_ld.tsv.gz": {"bytes": ld_stage.stat().st_size, "sha256": task["signed_ld_sha256"]},
            "signed_ld.qc.tsv": {"bytes": ld_qc_stage.stat().st_size, "sha256": task["ld_qc_sha256"]},
            "signed_ld.validation.json": {
                "bytes": len(ld_validation_content), "sha256": task["ld_validation_sha256"],
            },
            "signed_ld.runtime_packages.tsv": {
                "bytes": len(ld_runtime_package_content),
                "sha256": task["ld_runtime_package_lock_sha256"],
            },
            "signed_ld.runtime_attestation.tsv": {
                "bytes": ld_runtime_attestation_stage.stat().st_size,
                "sha256": task["ld_runtime_attestation_sha256"],
            },
            "runtime_packages.tsv": {
                "bytes": len(runtime_package_content),
                "sha256": task["runtime_package_lock_sha256"],
            },
            "task.tsv": {"bytes": len(task_content), "sha256": _sha256_content(task_content)},
            "materialization.provenance.json": {
                "bytes": len(provenance_content), "sha256": _sha256_content(provenance_content),
            },
        }
        bundle_lock = {
            "schema_version": TASK_SCHEMA, "analysis_id": ANALYSIS_ID,
            "artifact_role": "IMMUTABLE_PRE_RESULT_MATERIALIZED_LOCUS",
            "locus_entry_id": locus_entry_id, "results_accessed_before_lock": False,
            "output_records": output_records, "task_sha256": _sha256_content(task_content),
            "manifest_row_sha256": digest_json(manifest_row),
            "pre_family_lock": family["lock_identity"], "policy": family["policy_identity"],
            "scalar_n_convention": SCALAR_N_CONVENTION,
            "execution_amendment": EXECUTION_AMENDMENT,
            "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
            "ld_runtime_package_lock_sha256": task["ld_runtime_package_lock_sha256"],
            "ld_runtime_attestation_sha256": task["ld_runtime_attestation_sha256"],
            "reference_state": reference_state["state"],
            "archive_manager": archive_manager_identity,
            "script_identities": provenance["script_identities"],
        }
        bundle_content = json.dumps(bundle_lock, indent=2, sort_keys=True).encode() + b"\n"
        _publish_materialized_directory(
            root, stage, destination_rel,
            {
                "task.tsv": task_content,
                "materialization.provenance.json": provenance_content,
                "bundle.lock.json": bundle_content,
            },
        )
        return verify_materialized_bundle(root, destination_rel, family=family)
    except BaseException:
        # Never erase a failed materialization.  The measured coordinator owns
        # the pair lease and atomically quarantines this complete attempt-local
        # tree (including LD/QC/runtime evidence) after child exit or on restart.
        if stage.exists():
            try:
                fsync_directory(stage)
                fsync_directory(stage.parent)
            except OSError:
                # Preserve the original scientific/materialization failure;
                # the coordinator's quarantine pass fails closed if it cannot
                # inventory or durably move the retained tree.
                pass
        raise


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_materialized_bundle(
    root: Path, bundle_rel: Path | str, *, family: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    bundle_rel = Path(bundle_rel)
    bundle = safe_path(root, bundle_rel, "materialized bundle", must_exist=True)
    observed_dir = os.lstat(bundle)
    if not stat.S_ISDIR(observed_dir.st_mode) or stat.S_ISLNK(observed_dir.st_mode):
        raise FineMappingError("materialized bundle is not a real directory")
    lock, lock_identity = read_json(root, bundle_rel / "bundle.lock.json")
    if (
        lock.get("schema_version") != TASK_SCHEMA or lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("artifact_role") != "IMMUTABLE_PRE_RESULT_MATERIALIZED_LOCUS"
        or lock.get("results_accessed_before_lock") is not False
        or lock.get("scalar_n_convention") != SCALAR_N_CONVENTION
        or lock.get("execution_amendment") != EXECUTION_AMENDMENT
        or lock.get("execution_amendment_sha256") != EXECUTION_AMENDMENT_SHA256
    ):
        raise FineMappingError("materialized bundle lock identity drifted")
    records = lock.get("output_records")
    expected_names = {
        "trait1.tsv.gz", "trait2.tsv.gz", "variant_order.tsv", "signed_ld.tsv.gz",
        "signed_ld.qc.tsv", "signed_ld.validation.json",
        "signed_ld.runtime_packages.tsv", "signed_ld.runtime_attestation.tsv",
        "runtime_packages.tsv",
        "task.tsv", "materialization.provenance.json",
    }
    if not isinstance(records, Mapping) or set(records) != expected_names:
        raise FineMappingError("materialized bundle file family is incomplete")
    observed_names = {path.name for path in bundle.iterdir()}
    if observed_names != expected_names | {"bundle.lock.json"}:
        raise FineMappingError("materialized bundle contains missing or unexpected files")
    identities: dict[str, dict[str, object]] = {}
    for name, expected in records.items():
        identity = stable_identity(root, bundle_rel / name)
        if identity["bytes"] != expected.get("bytes") or identity["sha256"] != expected.get("sha256"):
            raise FineMappingError(f"materialized bundle artifact drifted: {name}")
        identities[name] = identity
    provenance, _ = read_json(root, bundle_rel / "materialization.provenance.json")
    manifest_row = provenance.get("manifest_row")
    if (
        provenance.get("schema_version") != "sleep-atlas-track-b-finemapping-materialization.1"
        or provenance.get("results_accessed") is not False
        or not isinstance(manifest_row, Mapping)
        or lock.get("manifest_row_sha256") != digest_json(manifest_row)
        or lock.get("script_identities") != provenance.get("script_identities")
        or provenance.get("execution_amendment") != EXECUTION_AMENDMENT
        or provenance.get("execution_amendment_sha256") != EXECUTION_AMENDMENT_SHA256
    ):
        raise FineMappingError("materialization provenance drifted")
    for relative, expected in lock.get("script_identities", {}).items():
        if not identity_matches(stable_identity(root, relative), expected):
            raise FineMappingError(f"materialization helper/script drifted: {relative}")
    if family is None:
        family = validate_pre_materialization_family(root)
    if (
        lock.get("pre_family_lock") != family["lock_identity"]
        or lock.get("policy") != family["policy_identity"]
    ):
        raise FineMappingError("materialized bundle does not bind the current frozen pre-family")
    pre_matches = [
        row for row in family["rows"] if row["locus_entry_id"] == lock.get("locus_entry_id")
    ]
    if len(pre_matches) != 1:
        raise FineMappingError("materialized bundle locus is outside frozen pre-family")
    pre_row = pre_matches[0]
    for field in (
        "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
        "CHR", "START", "STOP", "ld_block_id", "inclusion_sources",
        "local_evidence_ids", "method_availability", "claim_limit",
    ):
        if manifest_row.get(field) != pre_row[field]:
            raise FineMappingError(f"materialized bundle changed frozen pre-family field: {field}")
    if manifest_row.get("pleiotropy_evidence_ids") != pre_row["pleiotropy_clump_evidence_ids"]:
        raise FineMappingError("materialized bundle lost the full retained clump evidence union")
    task_rel = bundle_rel / "task.tsv"
    if lock.get("task_sha256") != identities["task.tsv"]["sha256"]:
        raise FineMappingError("materialized task hash differs from bundle lock")
    validated = deep_validate_locus_task(root, task_rel, expected_manifest_row=manifest_row)
    if provenance.get("ld_validation_attestation") != validated["ld_validation_attestation"]:
        raise FineMappingError("materialization provenance differs from signed-LD attestation")
    if (
        provenance.get("ld_runtime_attestation") != validated["ld_runtime_attestation"]
        or provenance.get("ld_runtime_package_lock_sha256")
        != validated["task"]["ld_runtime_package_lock_sha256"]
        or lock.get("ld_runtime_package_lock_sha256")
        != validated["task"]["ld_runtime_package_lock_sha256"]
        or lock.get("ld_runtime_attestation_sha256")
        != validated["task"]["ld_runtime_attestation_sha256"]
    ):
        raise FineMappingError("materialized bundle differs from signed-LD runtime attestation")
    return {
        "bundle_rel": str(bundle_rel), "bundle_lock": lock,
        "bundle_lock_identity": lock_identity, "provenance": provenance,
        "manifest_row": dict(manifest_row), "task": validated,
    }


def _materialized_manifest_fields(policy: Mapping[str, Any]) -> list[str]:
    fields = policy.get("locus_entry_contract", {}).get("locus_manifest_schema")
    if not isinstance(fields, list) or any(not isinstance(field, str) for field in fields):
        raise FineMappingError("policy lacks exact materialized locus-manifest schema")
    return fields


def freeze_materialized_manifest(
    root: Path = ROOT, *, pre_manifest_rel: Path | str = PRE_MANIFEST_REL,
    pre_lock_rel: Path | str = PRE_LOCK_REL,
    destination_root_rel: Path | str = INPUT_ROOT_REL,
) -> dict[str, Any]:
    """Exclusively publish the full materialized family before result access."""

    family = validate_pre_materialization_family(root, pre_manifest_rel, pre_lock_rel)
    policy = family["policy"]
    canonical_paths = [
        Path(policy["required_science_outputs"][key]["path"])
        for key in ("output_10", "output_11", "output_12")
    ]
    if any(safe_path(root, path, "canonical science output").exists() for path in canonical_paths):
        raise FineMappingError("fine-mapping/coloc results exist before materialized family freeze")
    materialized_path = safe_path(root, MATERIALIZED_MANIFEST_REL, "materialized manifest")
    materialized_lock_path = safe_path(root, MATERIALIZED_LOCK_REL, "materialized lock")
    zero_path = safe_path(root, ZERO_PROVENANCE_REL, "zero-family provenance")
    if family["state"] == ZERO_STATE:
        if materialized_path.exists() or materialized_lock_path.exists():
            raise FineMappingError("zero family cannot coexist with a materialized locus manifest")
        payload = {
            "schema_version": ZERO_SCHEMA, "analysis_id": ANALYSIS_ID, "state": ZERO_STATE,
            "artifact_role": "SEALED_ZERO_FAMILY_NOT_SCIENTIFIC_RESULT",
            "pre_family_manifest": family["manifest_identity"],
            "pre_family_lock": family["lock_identity"], "policy": family["policy_identity"],
            "upstream_family_provenance": family["lock"]["upstream_family_provenance"],
            "upstream_family_counts": {"eligible_evidence": 0, "pair_loci": 0},
            "science_outputs_10_11_12_created": False,
            "header_only_science_outputs_created": False,
            "scalar_n_convention": SCALAR_N_CONVENTION,
            "execution_amendment": EXECUTION_AMENDMENT,
            "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
            "script_identities": _runtime_script_identities(root),
        }
        content = json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"
        if zero_path.exists() or zero_path.is_symlink():
            observed, identity = read_json(root, ZERO_PROVENANCE_REL)
            if observed != payload:
                raise FineMappingError("existing zero-family provenance differs")
            return {"state": ZERO_STATE, "zero_provenance": identity}
        identity = publish_bytes_no_replace(root, ZERO_PROVENANCE_REL, content)
        return {"state": ZERO_STATE, "zero_provenance": identity}
    if zero_path.exists() or zero_path.is_symlink():
        raise FineMappingError("nonzero family cannot coexist with zero-family provenance")
    manifest_already_present = materialized_path.exists() or materialized_path.is_symlink()
    lock_already_present = materialized_lock_path.exists() or materialized_lock_path.is_symlink()
    if manifest_already_present and lock_already_present:
        return validate_materialized_manifest(root, pre_family=family)
    if lock_already_present and not manifest_already_present:
        raise FineMappingError("materialized lock exists without its manifest")

    bundle_root = Path(destination_root_rel)
    rows: list[dict[str, str]] = []
    bundles: dict[str, dict[str, object]] = {}
    for pre_row in family["rows"]:
        locus_id = pre_row["locus_entry_id"]
        verified = verify_materialized_bundle(root, bundle_root / locus_id, family=family)
        manifest_row = verified["manifest_row"]
        rows.append(manifest_row)
        bundles[locus_id] = verified["bundle_lock_identity"]
    fields = _materialized_manifest_fields(policy)
    if any(list(row) != fields for row in rows):
        raise FineMappingError("materializer produced a row outside exact policy manifest schema")
    content = tsv_bytes(fields, rows)
    manifest_identity_expected = {
        "path": str(MATERIALIZED_MANIFEST_REL), "bytes": len(content),
        "sha256": _sha256_content(content),
    }
    lock = {
        "schema_version": MATERIALIZED_LOCK_SCHEMA, "analysis_id": ANALYSIS_ID,
        "artifact_role": "PRE_RESULT_COMPLETE_MATERIALIZED_LOCUS_FAMILY",
        "state": "READY_FRESH_PROCESS_SUSIE_COLOC",
        "results_accessed_before_lock": False, "manifest": manifest_identity_expected,
        "locus_count": len(rows), "locus_entry_ids_in_order": [row["locus_entry_id"] for row in rows],
        "pre_family_manifest": family["manifest_identity"],
        "pre_family_lock": family["lock_identity"], "policy": family["policy_identity"],
        "unavailable_manifest": family["unavailable_identity"],
        "unavailable_count": len(family["unavailable_rows"]),
        "unavailable_evidence_ids_sha256": digest_json(family["unavailable_evidence_ids"]),
        "materialized_bundles": bundles, "scalar_n_convention": SCALAR_N_CONVENTION,
        "execution_amendment": EXECUTION_AMENDMENT,
        "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
        "representative_order_rule": (
            "SMALLEST_THEN_LOWER_MEDIAN_THEN_LARGEST_DISTINCT_LOCUS_BY_FROZEN_"
            "VARIANT_COUNT;TIES_BY_FROZEN_MANIFEST_ORDER;THEN_ALL_REMAINING_IN_MANIFEST_ORDER"
        ),
        "no_subsetting": {
            "maximum_locus_variants": None, "locus_splitting": False,
            "variant_thinning": False, "lead_centered_truncation": False,
        },
        "script_identities": _runtime_script_identities(root),
    }
    lock_content = json.dumps(lock, indent=2, sort_keys=True).encode() + b"\n"
    # Manifest first, lock last.  A crash can leave only a visibly unsealed
    # manifest; a later invocation refuses to overwrite it.
    if manifest_already_present:
        observed_manifest = stable_identity(root, MATERIALIZED_MANIFEST_REL)
        if not identity_matches(observed_manifest, manifest_identity_expected):
            raise FineMappingError(
                "unsealed materialized manifest is not the exact deterministic crash prefix"
            )
    else:
        publish_bytes_no_replace(root, MATERIALIZED_MANIFEST_REL, content)
    # Close all input TOCTOU windows before the final seal.
    validate_pre_materialization_family(root, pre_manifest_rel, pre_lock_rel)
    for pre_row in family["rows"]:
        verify_materialized_bundle(root, bundle_root / pre_row["locus_entry_id"], family=family)
    publish_bytes_no_replace(root, MATERIALIZED_LOCK_REL, lock_content)
    return validate_materialized_manifest(root, pre_family=family)


def validate_materialized_manifest(
    root: Path = ROOT, *, pre_family: Mapping[str, Any] | None = None,
    manifest_rel: Path | str = MATERIALIZED_MANIFEST_REL,
    lock_rel: Path | str = MATERIALIZED_LOCK_REL,
) -> dict[str, Any]:
    """Side-effect-free, exhaustive validation of the full frozen task family."""

    if pre_family is None:
        pre_family = validate_pre_materialization_family(root)
    if pre_family["state"] == ZERO_STATE:
        if safe_path(root, manifest_rel, "materialized manifest").exists():
            raise FineMappingError("zero family unexpectedly has a materialized manifest")
        zero, identity = read_json(root, ZERO_PROVENANCE_REL)
        if (
            zero.get("schema_version") != ZERO_SCHEMA or zero.get("state") != ZERO_STATE
            or zero.get("science_outputs_10_11_12_created") is not False
            or zero.get("pre_family_lock") != pre_family["lock_identity"]
            or zero.get("execution_amendment") != EXECUTION_AMENDMENT
            or zero.get("execution_amendment_sha256") != EXECUTION_AMENDMENT_SHA256
        ):
            raise FineMappingError("zero-family production provenance drifted")
        return {
            "state": ZERO_STATE,
            "zero_provenance": identity,
            "rows": [],
            "policy": pre_family["policy"],
            "pre_family": pre_family,
        }
    policy = pre_family["policy"]
    fields = _materialized_manifest_fields(policy)
    _, rows, manifest_identity = read_tsv(
        root, manifest_rel, fields, allow_header_only=(len(pre_family["rows"]) == 0),
    )
    lock, lock_identity = read_json(root, lock_rel)
    if (
        set(lock) != {
            "schema_version", "analysis_id", "artifact_role", "state",
            "results_accessed_before_lock", "manifest", "locus_count",
            "locus_entry_ids_in_order", "pre_family_manifest", "pre_family_lock", "policy",
            "unavailable_manifest", "unavailable_count", "unavailable_evidence_ids_sha256",
            "materialized_bundles", "scalar_n_convention", "representative_order_rule",
            "execution_amendment", "execution_amendment_sha256", "no_subsetting",
            "script_identities",
        }
        or lock.get("schema_version") != MATERIALIZED_LOCK_SCHEMA
        or lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("artifact_role") != "PRE_RESULT_COMPLETE_MATERIALIZED_LOCUS_FAMILY"
        or lock.get("state") != "READY_FRESH_PROCESS_SUSIE_COLOC"
        or lock.get("results_accessed_before_lock") is not False
        or not identity_matches(lock.get("manifest"), manifest_identity)
        or lock.get("locus_count") != len(rows)
        or lock.get("locus_entry_ids_in_order") != [row["locus_entry_id"] for row in rows]
        or lock.get("pre_family_manifest") != pre_family["manifest_identity"]
        or lock.get("pre_family_lock") != pre_family["lock_identity"]
        or lock.get("policy") != pre_family["policy_identity"]
        or lock.get("unavailable_manifest") != pre_family["unavailable_identity"]
        or lock.get("unavailable_count") != len(pre_family["unavailable_rows"])
        or lock.get("unavailable_evidence_ids_sha256")
        != digest_json(pre_family["unavailable_evidence_ids"])
        or lock.get("scalar_n_convention") != SCALAR_N_CONVENTION
        or lock.get("execution_amendment") != EXECUTION_AMENDMENT
        or lock.get("execution_amendment_sha256") != EXECUTION_AMENDMENT_SHA256
        or lock.get("no_subsetting") != {
            "maximum_locus_variants": None, "locus_splitting": False,
            "variant_thinning": False, "lead_centered_truncation": False,
        }
    ):
        raise FineMappingError("materialized locus-manifest lock drifted")
    if len(rows) != len(pre_family["rows"]):
        raise FineMappingError("materialized family is not complete")
    bundles = lock.get("materialized_bundles")
    if not isinstance(bundles, Mapping) or set(bundles) != {row["locus_entry_id"] for row in rows}:
        raise FineMappingError("materialized lock omits a locus bundle")
    verified_bundles: dict[str, dict[str, Any]] = {}
    for row, pre_row in zip(rows, pre_family["rows"], strict=True):
        for field in (
            "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
            "CHR", "START", "STOP", "ld_block_id", "inclusion_sources",
            "local_evidence_ids", "method_availability", "claim_limit",
        ):
            if row[field] != pre_row[field]:
                raise FineMappingError(f"materialized manifest changed frozen family field: {field}")
        if row["pleiotropy_evidence_ids"] != pre_row["pleiotropy_clump_evidence_ids"]:
            raise FineMappingError("materialized manifest lost full clump-evidence union")
        if row["memory_admission_status"] != "ADMITTED_WHOLE_LOCUS":
            raise FineMappingError("blocked locus was incorrectly placed in runnable manifest")
        count = _integer(row["variant_count"], "materialized variant count", minimum=50)
        if int(row["estimated_peak_rss_bytes"]) != estimate_peak_rss_bytes(count, policy):
            raise FineMappingError("materialized manifest memory estimate drifted")
        task_path = Path(row["task_path"])
        bundle_rel = task_path.parent
        verified = verify_materialized_bundle(root, bundle_rel, family=pre_family)
        if verified["manifest_row"] != row:
            raise FineMappingError("materialized manifest row differs from sealed bundle")
        if bundles[row["locus_entry_id"]] != verified["bundle_lock_identity"]:
            raise FineMappingError("materialized bundle lock identity drifted")
        verified_bundles[row["locus_entry_id"]] = verified
    for relative, expected in lock.get("script_identities", {}).items():
        if not identity_matches(stable_identity(root, relative), expected):
            raise FineMappingError(f"materialized family helper/script drifted: {relative}")
    return {
        "state": "READY_FRESH_PROCESS_SUSIE_COLOC", "rows": rows, "lock": lock,
        "manifest_identity": manifest_identity, "lock_identity": lock_identity,
        "pre_family": pre_family, "bundles": verified_bundles, "policy": policy,
    }


def representative_execution_order(rows: Sequence[Mapping[str, str]]) -> list[str]:
    """Small/lower-median/large by frozen burden, then manifest order."""

    ids = [row["locus_entry_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise FineMappingError("cannot order a duplicated locus family")
    if len(rows) < 3:
        return ids
    indexed = [(int(row["variant_count"]), index, row["locus_entry_id"]) for index, row in enumerate(rows)]
    ranked = sorted(indexed)
    candidates = [ranked[0], ranked[(len(ranked) - 1) // 2], ranked[-1]]
    representatives: list[str] = []
    for _, _, locus_id in candidates:
        if locus_id not in representatives:
            representatives.append(locus_id)
    return representatives + [locus_id for locus_id in ids if locus_id not in representatives]


def _parse_optional_number(value: str, label: str) -> float | None:
    if value in {"", "NA"}:
        return None
    return _number(value, label)


def _true_false(value: str, label: str) -> bool:
    if value == "TRUE":
        return True
    if value == "FALSE":
        return False
    raise FineMappingError(f"{label} must be TRUE or FALSE")


def _row_identity(row: Mapping[str, str], task: Mapping[str, str]) -> None:
    expected = {
        "analysis_id": ANALYSIS_ID, "pair_id": task["pair_id"],
        "family_role": task["family_role"], "locus_entry_id": task["locus_entry_id"],
        "locus_id": task["ld_block_id"], "CHR": task["CHR"], "START": task["START"],
        "STOP": task["STOP"],
    }
    if any(row.get(field) != value for field, value in expected.items()):
        raise FineMappingError("engine row escaped its exact pair x official block")


def _n_claim_capped(task: Mapping[str, str], prefix: str | None = None) -> bool:
    prefixes = (prefix,) if prefix is not None else ("trait1", "trait2")
    return any(task[f"{item}_N_dispersion_status"] == N_DISPERSION_VARIABLE for item in prefixes)


def _validate_n_diagnostics(
    rows: Sequence[Mapping[str, str]], task_validation: Mapping[str, Any],
) -> dict[str, dict[str, str]]:
    task = task_validation["task"]
    roles = {"TRAIT1_SLEEP": "trait1", "TRAIT2_EXTERNAL": "trait2"}
    if len(rows) != 2 or {row["dataset_role"] for row in rows} != set(roles):
        raise FineMappingError("sample-size diagnostics do not retain exactly both trait roles")
    output: dict[str, dict[str, str]] = {}
    suffixes = (
        "count", "min", "q1", "median", "q3", "max", "max_to_min_ratio",
        "fraction_below_90pct_max", "fraction_below_50pct_max",
    )
    for raw in rows:
        row = dict(raw)
        role = row["dataset_role"]
        prefix = roles[role]
        expected = {
            "analysis_id": ANALYSIS_ID, "pair_id": task["pair_id"],
            "family_role": task["family_role"], "locus_entry_id": task["locus_entry_id"],
            "dataset_role": role, "trait_id": task[prefix],
            "per_snp_N_sha256": task[f"{prefix}_per_snp_N_sha256"],
            "scalar_N_rule": task[f"{prefix}_scalar_N_rule"],
            "dispersion_status": task[f"{prefix}_N_dispersion_status"],
            "claim_cap": "TRUE" if _n_claim_capped(task, prefix) else "FALSE",
        }
        if any(row[field] != value for field, value in expected.items()):
            raise FineMappingError("sample-size diagnostics differ from the full-vector task lock")
        if _integer(row["per_snp_N_count"], "sample-size diagnostic count", minimum=1) != int(
            task[f"{prefix}_per_snp_N_count"]
        ):
            raise FineMappingError("sample-size diagnostic count differs from task")
        for suffix in suffixes[1:]:
            if not math.isclose(
                _number(row[f"per_snp_N_{suffix}"], f"sample-size {suffix}"),
                _number(task[f"{prefix}_per_snp_N_{suffix}"], f"task sample-size {suffix}"),
                rel_tol=1e-12, abs_tol=1e-12,
            ):
                raise FineMappingError("sample-size diagnostic numeric summary differs from task")
        if not math.isclose(
            _number(row["scalar_N"], "sample-size scalar N"),
            _number(task[f"{prefix}_scalar_N"], "task scalar N"),
            rel_tol=1e-12, abs_tol=1e-12,
        ):
            raise FineMappingError("sample-size diagnostic scalar N differs from task")
        output[role] = row
    return output


def _validate_runtime_attestation(
    rows: Sequence[Mapping[str, str]], task_validation: Mapping[str, Any],
) -> dict[str, str]:
    if len(rows) != 1:
        raise FineMappingError("runtime attestation must contain exactly one row")
    task = task_validation["task"]
    expected = {
        "analysis_id": ANALYSIS_ID, "locus_entry_id": task["locus_entry_id"],
        "execution_amendment_sha256": EXECUTION_AMENDMENT_SHA256,
        "runtime_package_lock_sha256": task["runtime_package_lock_sha256"],
        "precheck_status": "PASS_EXACT_PATHS_AND_BYTES",
        "postcheck_status": "PASS_EXACT_LOADED_CLOSURE_PATHS_AND_BYTES",
        "loaded_namespace_count": str(len(PINNED_RUNTIME_NAMESPACE_CLOSURE)),
    }
    if dict(rows[0]) != expected:
        raise FineMappingError("R runtime pre/post path/byte attestation drifted")
    return dict(rows[0])


def _validate_diagnostics(
    rows: Sequence[Mapping[str, str]], task_validation: Mapping[str, Any],
) -> dict[str, dict[str, str]]:
    task = task_validation["task"]
    expected_roles = {
        "TRAIT1_SLEEP": task["trait1"], "TRAIT2_EXTERNAL": task["trait2"],
    }
    if len(rows) != 2 or {row["dataset_role"] for row in rows} != set(expected_roles):
        raise FineMappingError("diagnostics do not contain exactly both trait roles")
    output: dict[str, dict[str, str]] = {}
    allowed = {
        "PASS", "FAILED_ALLELE_OR_ORDER", "FAILED_LD_QC",
        "FAILED_MODEL_CONVERGENCE", "FAILED_KRIGING_QC", "FAILED_REQUIRED_DIAGNOSTIC",
    }
    for raw in rows:
        row = dict(raw)
        role = row["dataset_role"]
        if (
            row["analysis_id"] != ANALYSIS_ID or row["pair_id"] != task["pair_id"]
            or row["family_role"] != task["family_role"]
            or row["locus_entry_id"] != task["locus_entry_id"]
            or row["trait_id"] != expected_roles[role]
            or row["diagnostic_status"] not in allowed
            or _integer(row["variant_count"], "diagnostic variant count", minimum=1)
            != int(task["variant_count"])
            or _integer(row["ld_variant_count"], "diagnostic LD count", minimum=1)
            != int(task["variant_count"])
            or _integer(row["allele_match_count"], "diagnostic allele count", minimum=1)
            != int(task["variant_count"])
            or _integer(row["reference_sample_size"], "diagnostic reference N", minimum=1)
            != int(task["reference_sample_size"])
        ):
            raise FineMappingError("diagnostic identity/count/status drifted")
        converged = _true_false(row["model_converged"], "model_converged")
        niter = _parse_optional_number(row["niter"], "diagnostic niter")
        rss = _parse_optional_number(row["rss_ld_s"], "RSS-LD s")
        outlier_count = _integer(
            row["kriging_allele_switch_outlier_count"], "kriging outlier count", minimum=0,
        )
        outliers = [] if row["kriging_allele_switch_outliers"] == "NONE" else row[
            "kriging_allele_switch_outliers"
        ].split(";")
        if len(outliers) != outlier_count or len(outliers) != len(set(outliers)):
            raise FineMappingError("kriging outlier list/count drifted")
        if row["diagnostic_status"] == "PASS":
            if not converged or niter is None or niter > 1000 or rss is None or not 0 <= rss <= 1:
                raise FineMappingError("PASS diagnostic lacks valid convergence/RSS evidence")
            if outlier_count:
                raise FineMappingError("PASS diagnostic contains kriging outliers")
            prefix = "trait1" if role == "TRAIT1_SLEEP" else "trait2"
            expected_warnings: set[str] = set()
            if rss > 0.10:
                expected_warnings.add("WARNING_RSS_LD_S_GT_0.10_NONFATAL")
            if _n_claim_capped(task, prefix):
                expected_warnings.add(
                    "WARNING_VARIABLE_PER_SNP_N_MEDIAN_SCALAR_APPROXIMATION_CLAIM_CAPPED"
                )
            observed_warnings = set() if row["error"] == "NA" else set(row["error"].split(";"))
            if observed_warnings != expected_warnings:
                raise FineMappingError("PASS diagnostic warning family differs from locked RSS/N rules")
        else:
            if row["error"] in {"", "NA"}:
                raise FineMappingError("failed diagnostic lacks retained error")
        for field, expected in (
            ("ld_symmetry_max_abs", task_validation["ld_metrics"]["symmetry_max_abs"]),
            ("ld_diagonal_max_abs", task_validation["ld_metrics"]["diagonal_max_abs"]),
            ("ld_minimum_eigenvalue", task_validation["ld_metrics"]["minimum_eigenvalue"]),
        ):
            if not math.isclose(_number(row[field], field), float(expected), rel_tol=0, abs_tol=1e-8):
                raise FineMappingError("diagnostic LD metric differs from prelaunch validation")
        output[role] = row
    return output


def _validate_credible_sets(
    rows: Sequence[Mapping[str, str]], task_validation: Mapping[str, Any],
    diagnostics: Mapping[str, Mapping[str, str]],
) -> dict[tuple[str, str], dict[str, str]]:
    task = task_validation["task"]
    summary_snps = {row["SNP"] for row in task_validation["order_rows"]}
    output: dict[tuple[str, str], dict[str, str]] = {}
    counts = Counter()
    for raw in rows:
        row = dict(raw)
        role = row["trait_role"]
        expected_trait = task["trait1"] if role == "TRAIT1_SLEEP" else task["trait2"] if role == "TRAIT2_EXTERNAL" else None
        expected_analysis_status = (
            "COMPLETE" if role in diagnostics and diagnostics[role]["diagnostic_status"] == "PASS"
            else "FAILED_QC"
        )
        if (
            expected_trait is None or row["analysis_id"] != ANALYSIS_ID
            or row["pair_id"] != task["pair_id"] or row["family_role"] != task["family_role"]
            or row["locus_entry_id"] != task["locus_entry_id"] or row["trait_id"] != expected_trait
            or row["analysis_status"] != expected_analysis_status
            or row["claim_limit"] != task_validation["policy"]["claim_limits"]["fine_mapping"]
            or _true_false(row["model_converged"], "credible-set convergence")
            != _true_false(diagnostics[role]["model_converged"], "diagnostic convergence")
        ):
            raise FineMappingError("credible-set identity/status/claim drifted")
        component = _integer(row["component_index"], "credible-set component", minimum=1)
        if component > 10 or row["signal_id"] != f"L{component}":
            raise FineMappingError("credible-set component exceeds locked L=10")
        members = row["credible_set_snps"].split(";")
        size = _integer(row["credible_set_size"], "credible-set size", minimum=1)
        if (
            len(members) != size or len(members) != len(set(members))
            or not set(members) <= summary_snps or row["lead_snp"] not in members
        ):
            raise FineMappingError("credible-set member family is invalid")
        requested = _number(row["requested_coverage"], "requested coverage")
        achieved = _number(row["achieved_coverage"], "achieved coverage")
        minimum = _number(row["min_abs_corr"], "credible-set minimum correlation")
        mean = _number(row["mean_abs_corr"], "credible-set mean correlation")
        median = _number(row["median_abs_corr"], "credible-set median correlation")
        if (
            not math.isclose(requested, 0.95, rel_tol=0, abs_tol=1e-12)
            or achieved + 1e-8 < requested or achieved > 1 + 1e-8
            or minimum + 1e-12 < 0.5
            or any(value < 0 or value > 1 + 1e-8 for value in (minimum, mean, median))
        ):
            raise FineMappingError("credible-set coverage or purity violates locked model")
        _probability(row["lead_pip"], "credible-set lead PIP")
        if row["cs_log10bf"] != "NA":
            _number(row["cs_log10bf"], "credible-set log10 Bayes factor")
        key = (role, row["signal_id"])
        if key in output:
            raise FineMappingError("duplicate credible-set signal")
        output[key] = row
        counts[role] += 1
    for role, diagnostic in diagnostics.items():
        if counts[role] != int(diagnostic["credible_set_count"]):
            raise FineMappingError("credible-set rows differ from diagnostics count")
    return output


def _classify_variant_rows(
    rows: Sequence[Mapping[str, str]], role: str, task_validation: Mapping[str, Any],
    diagnostic: Mapping[str, str], canonical_fields: Sequence[str],
) -> list[dict[str, str]]:
    task = task_validation["task"]
    summary = task_validation["summary1_rows"] if role == "TRAIT1_SLEEP" else task_validation["summary2_rows"]
    trait = task["trait1"] if role == "TRAIT1_SLEEP" else task["trait2"]
    if len(rows) != len(summary) or len(rows) != int(task["variant_count"]):
        raise FineMappingError("engine variant output does not retain the complete locus")
    objective_pass = diagnostic["diagnostic_status"] == "PASS"
    converged = _true_false(diagnostic["model_converged"], "diagnostic convergence")
    rss = _parse_optional_number(diagnostic["rss_ld_s"], "RSS-LD s")
    rss_warning = objective_pass and rss is not None and rss > 0.10
    n_warning = _n_claim_capped(task, "trait1" if role == "TRAIT1_SLEEP" else "trait2")
    model_has_high = False
    model_has_small = False
    parsed: list[tuple[dict[str, str], bool, bool]] = []
    for index, (raw, source) in enumerate(zip(rows, summary, strict=True)):
        row = dict(raw)
        _row_identity(row, task)
        if row["trait_role"] != role or row["trait_id"] != trait:
            raise FineMappingError("engine variant row escaped trait role")
        for field in SUMMARY_FIELDS:
            if row[field] != source[field]:
                raise FineMappingError(f"engine variant output changed materialized {field}")
        if row["prior_role"] != "PRIMARY_FLAT":
            raise FineMappingError("engine used a non-flat or non-primary variant prior")
        expected_weight = 1.0 / len(rows)
        if not math.isclose(_number(row["prior_weight"], "flat prior weight"), expected_weight, rel_tol=0, abs_tol=1e-12):
            raise FineMappingError("engine prior weights are not exactly flat")
        if row["diagnostic_status"] != diagnostic["diagnostic_status"]:
            raise FineMappingError("variant diagnostic status differs from diagnostic row")
        if (
            _true_false(row["model_converged"], "variant convergence") != converged
            or row["niter"] != diagnostic["niter"]
            or row["rss_ld_s"] != diagnostic["rss_ld_s"]
            or row["error"] != diagnostic["error"]
        ):
            raise FineMappingError("variant model/diagnostic fields drifted")
        pip = _parse_optional_number(row["PIP"], "PIP")
        high = (
            objective_pass and not n_warning and converged and pip is not None
            and _probability(row["PIP"], "PIP") >= 0.5
        )
        cs_ids = [] if row["credible_set_ids"] == "NONE" else row["credible_set_ids"].split(";")
        if len(cs_ids) != len(set(cs_ids)) or any(not re.fullmatch(r"L(?:[1-9]|10)", item) for item in cs_ids):
            raise FineMappingError("variant credible-set IDs are duplicated or outside L=10")
        annotations = [
            [] if row[field] == "NA" else row[field].split(";")
            for field in (
                "credible_set_sizes", "credible_set_requested_coverage",
                "credible_set_achieved_coverage", "credible_set_min_abs_corr",
                "credible_set_mean_abs_corr",
            )
        ]
        sizes, requested, achieved, minimums, means = annotations
        del means
        if any(len(values) != len(cs_ids) for values in annotations):
            raise FineMappingError("variant credible-set annotations are not parallel")
        if any(
            not math.isclose(_number(value, "variant CS requested coverage"), 0.95, rel_tol=0, abs_tol=1e-12)
            for value in requested
        ) or any(
            _number(value, "variant CS achieved coverage") + 1e-8 < 0.95 for value in achieved
        ):
            raise FineMappingError("variant credible-set coverage annotation drifted")
        small = objective_pass and not n_warning and converged and any(
            int(size) <= 10 and _number(minimum, "variant CS purity") >= 0.5
            for size, minimum in zip(sizes, minimums, strict=True)
        )
        model_has_high |= high
        model_has_small |= small
        parsed.append((row, high, small))
        if row["SNP"] != task_validation["order_rows"][index]["SNP"]:
            raise FineMappingError("engine variant rows are out of locked order")
        expected_outlier = row["SNP"] in (
            [] if diagnostic["kriging_allele_switch_outliers"] == "NONE"
            else diagnostic["kriging_allele_switch_outliers"].split(";")
        )
        if _true_false(row["kriging_allele_switch_outlier"], "variant kriging flag") != expected_outlier:
            raise FineMappingError("variant kriging flag differs from retained diagnostic outliers")
        if row["max_alpha_component"] != "NA" and not re.fullmatch(
            r"L(?:[1-9]|10)", row["max_alpha_component"]
        ):
            raise FineMappingError("variant maximum-alpha component is outside L=10")
    diffuse = (
        objective_pass and not n_warning and converged
        and not model_has_high and not model_has_small
    )
    output: list[dict[str, str]] = []
    label_order = [
        "HIGH_PIP_VARIANT", "SMALL_CREDIBLE_SET", "DIFFUSE_SIGNAL",
        "SAMPLE_SIZE_UNCERTAINTY", "LD_UNCERTAINTY", "MODEL_INSTABILITY",
    ]
    for row, high, small in parsed:
        labels: set[str] = set()
        if high:
            labels.add("HIGH_PIP_VARIANT")
        if small:
            labels.add("SMALL_CREDIBLE_SET")
        if diffuse:
            labels.add("DIFFUSE_SIGNAL")
        if n_warning:
            labels.add("SAMPLE_SIZE_UNCERTAINTY")
        if rss_warning or diagnostic["diagnostic_status"] in {
            "FAILED_ALLELE_OR_ORDER", "FAILED_LD_QC", "FAILED_KRIGING_QC",
            "FAILED_REQUIRED_DIAGNOSTIC",
        }:
            labels.add("LD_UNCERTAINTY")
        if diagnostic["diagnostic_status"] == "FAILED_MODEL_CONVERGENCE" or not converged:
            labels.add("MODEL_INSTABILITY")
        values = dict(row)
        values["ld_reference_id"] = task["reference_id"]
        values["ld_variant_order_sha256"] = task["variant_order_sha256"]
        values["fine_mapping_classifications"] = (
            ";".join(label for label in label_order if label in labels) if labels else "NA"
        )
        values["analysis_status"] = "COMPLETE" if objective_pass else "FAILED_QC"
        values["claim_limit"] = task_validation["policy"]["claim_limits"]["fine_mapping"]
        canonical = {field: values[field] for field in canonical_fields}
        if list(canonical) != list(canonical_fields):
            raise FineMappingError("internal canonical fine-mapping schema assembly failed")
        output.append(canonical)
    return output


def _complete_probability_row(row: Mapping[str, str]) -> tuple[float, float, float, float, float] | None:
    values = [row[field] for field in ("PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4")]
    if all(value == "NA" for value in values):
        return None
    if any(value == "NA" for value in values):
        raise FineMappingError("coloc probability row is partially missing")
    parsed = tuple(_probability(value, field) for field, value in zip(
        ("PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4"), values, strict=True,
    ))
    if not math.isclose(sum(parsed), 1.0, rel_tol=0, abs_tol=1e-6):
        raise FineMappingError("coloc posterior probabilities do not sum to one")
    return parsed  # type: ignore[return-value]


def _ratio_value(value: str, label: str) -> float | None:
    if value in {"", "NA"}:
        return None
    if value == "Inf":
        return math.inf
    return _number(value, label)


def classify_coloc_rows(
    raw_rows: Sequence[Mapping[str, str]], diagnostics: Mapping[str, Mapping[str, str]],
    task_validation: Mapping[str, Any], canonical_fields: Sequence[str],
) -> list[dict[str, str]]:
    """Apply exact Track B primary/prior-robust classification semantics."""

    task = task_validation["task"]
    grid = (1e-6, 5e-6, 1e-5, 5e-5)
    if not raw_rows:
        raise FineMappingError("R engine omitted the complete coloc prior family")
    by_pair: dict[tuple[str, str], dict[float, dict[str, str]]] = {}
    ordered: list[dict[str, str]] = []
    for raw in raw_rows:
        row = dict(raw)
        _row_identity(row, task)
        if (
            row["trait1"] != task["trait1"] or row["trait2"] != task["trait2"]
            or row["coloc_method"] != "coloc.susie"
            or _number(row["p1"], "coloc p1") != 1e-4
            or _number(row["p2"], "coloc p2") != 1e-4
            or _integer(row["nsnps"], "coloc SNP count", minimum=1)
            != int(task["variant_count"])
            or row["claim_limit"] != task_validation["policy"]["claim_limits"]["trait_coloc"]
        ):
            raise FineMappingError("coloc row method/trait/prior identity drifted")
        p12 = _number(row["p12"], "coloc p12")
        if p12 not in grid or row["prior_role"] != ("PRIMARY" if p12 == 1e-5 else "SENSITIVITY"):
            raise FineMappingError("coloc row p12/prior role is outside frozen grid")
        pair = (row["signal1"], row["signal2"])
        if (row["signal1"] == "NA") != (row["signal2"] == "NA"):
            raise FineMappingError("coloc row has a half-missing signal pair")
        if row["signal1"] != "NA" and any(
            not re.fullmatch(r"L(?:[1-9]|10)", value)
            for value in (row["signal1"], row["signal2"])
        ):
            raise FineMappingError("coloc signal pair is outside the locked L=10 model")
        if row["engine_status"] not in {
            "COLOC_SUSIE_COMPLETE", "COLOC_SUSIE_FAILURE",
            "COLOC_SUSIE_NO_ESTIMABLE_SIGNAL_PAIR", "NO_VALID_CREDIBLE_SET_PAIR",
            "INVALID_FINE_MAPPING_DIAGNOSTICS",
        }:
            raise FineMappingError("coloc engine status is unsupported")
        family = by_pair.setdefault(pair, {})
        if p12 in family:
            raise FineMappingError("duplicate coloc signal-pair/prior row")
        family[p12] = row
        probability = _complete_probability_row(row)
        if probability is not None:
            h3, h4 = probability[3], probability[4]
            ratio = math.inf if h3 == 0 and h4 > 0 else (h4 / h3 if h3 > 0 else math.nan)
            observed_ratio = _ratio_value(row["PP_H4_over_PP_H3"], "H4/H3")
            ratio_matches = (
                math.isnan(ratio) and observed_ratio is None
            ) or (
                observed_ratio is not None and (
                    (math.isinf(ratio) and math.isinf(observed_ratio))
                    or math.isclose(observed_ratio, ratio, rel_tol=1e-8, abs_tol=1e-10)
                )
            )
            if not ratio_matches:
                raise FineMappingError("coloc H4/H3 ratio drifted")
            if row["top_shared_variant"] not in {
                order_row["SNP"] for order_row in task_validation["order_rows"]
            }:
                raise FineMappingError("coloc top shared variant is outside locked SNP family")
            _probability(row["top_shared_variant_PP_H4"], "top shared-variant posterior")
        ordered.append(row)
    # Every emitted signal pair must cover every locked prior. No best-pair
    # pooling or result-ranked omission is accepted.
    if any(set(family) != set(grid) for family in by_pair.values()):
        raise FineMappingError("coloc signal pair does not cover the complete p12 grid")
    fatal = any(row["diagnostic_status"] != "PASS" for row in diagnostics.values())
    rss_warning = any(
        (_parse_optional_number(row["rss_ld_s"], "RSS-LD s") or 0.0) > 0.10
        for row in diagnostics.values()
    )
    n_warning = _n_claim_capped(task)
    decisions: dict[tuple[str, str], tuple[str, str, str]] = {}
    for pair, family in by_pair.items():
        primary = family[1e-5]
        primary_probability = _complete_probability_row(primary)
        complete = {
            prior: _complete_probability_row(family[prior]) for prior in grid
        }
        primary_strong = False
        primary_moderate = False
        primary_distinct = False
        if primary_probability is not None:
            h3, h4 = primary_probability[3], primary_probability[4]
            ratio_h4_h3 = math.inf if h3 == 0 and h4 > 0 else h4 / h3 if h3 > 0 else 0.0
            ratio_h3_h4 = math.inf if h4 == 0 and h3 > 0 else h3 / h4 if h4 > 0 else 0.0
            primary_strong = h4 >= 0.80 and ratio_h4_h3 >= 5
            primary_moderate = h4 >= 0.50 and ratio_h4_h3 >= 2
            primary_distinct = h3 >= 0.80 and ratio_h3_h4 >= 5
        robust = pair != ("NA", "NA") and all(
            probability is not None
            and probability[4] >= 0.80
            and (
                math.inf if probability[3] == 0 and probability[4] > 0
                else probability[4] / probability[3] if probability[3] > 0 else 0.0
            ) >= 5
            for probability in complete.values()
        )
        if fatal or any(row["engine_status"] in {"COLOC_SUSIE_FAILURE", "INVALID_FINE_MAPPING_DIAGNOSTICS"} for row in family.values()):
            classification, status = "INVALID_INPUT", "INVALID_INPUT"
        elif rss_warning or n_warning:
            classification, status = "INCONCLUSIVE", "TESTED"
        elif primary_strong and robust:
            classification, status = "STRONG_SHARED_SIGNAL", "TESTED"
        elif primary_moderate:
            classification, status = "MODERATE_SHARED_SIGNAL", "TESTED"
        elif primary_distinct:
            classification, status = "DISTINCT_SIGNALS", "TESTED"
        else:
            classification, status = "INCONCLUSIVE", "TESTED"
        decisions[pair] = (
            "TRUE" if primary_strong else "FALSE", "TRUE" if robust else "FALSE",
            classification + "\t" + status,
        )
    output: list[dict[str, str]] = []
    for row in ordered:
        pair = (row["signal1"], row["signal2"])
        primary_pass, robust, combined = decisions[pair]
        classification, analysis_status = combined.split("\t")
        values = dict(row)
        values["primary_shared_signal_rule_pass"] = primary_pass
        values["same_signal_pair_prior_robust"] = robust
        if rss_warning and not fatal:
            values["fine_mapping_qc"] = "LD_UNCERTAINTY_WARNING_RSS_LD_S_GT_0.10"
            values["ld_qc"] = "LD_UNCERTAINTY_WARNING_RSS_LD_S_GT_0.10"
        if n_warning and not fatal:
            values["fine_mapping_qc"] = (
                "SAMPLE_SIZE_UNCERTAINTY_VARIABLE_PER_SNP_N_MEDIAN_SCALAR_CLAIM_CAPPED"
            )
            if not rss_warning:
                values["ld_qc"] = "PASS"
        values["classification"] = classification
        values["analysis_status"] = analysis_status
        values["claim_limit"] = task_validation["policy"]["claim_limits"]["trait_coloc"]
        canonical = {field: values[field] for field in canonical_fields}
        if list(canonical) != list(canonical_fields):
            raise FineMappingError("internal canonical coloc schema assembly failed")
        output.append(canonical)
    return output


def adapt_engine_outputs(
    root: Path, output_dir_rel: Path | str, task_validation: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate raw R serialization and add exact canonical per-locus fragments."""

    output_dir_rel = Path(output_dir_rel)
    policy = task_validation["policy"]
    verify_runtime_package_lock(
        root, task_validation["task"]["runtime_package_lock_path"],
        expected_sha256=task_validation["task"]["runtime_package_lock_sha256"],
    )
    output_schemas = policy["required_science_outputs"]
    variant_fields = output_schemas["output_10"]["schema"]
    coloc_fields = output_schemas["output_12"]["schema"]
    _, raw1, _ = read_tsv(root, output_dir_rel / "trait1.raw.tsv", RAW_VARIANT_FIELDS)
    _, raw2, _ = read_tsv(root, output_dir_rel / "trait2.raw.tsv", RAW_VARIANT_FIELDS)
    _, diagnostic_rows, _ = read_tsv(
        root, output_dir_rel / "diagnostics.tsv",
        policy["diagnostic_contract"]["required_schema"],
    )
    diagnostics = _validate_diagnostics(diagnostic_rows, task_validation)
    _, n_diagnostic_rows, _ = read_tsv(
        root, output_dir_rel / "sample_size_diagnostics.tsv", N_DIAGNOSTIC_FIELDS,
    )
    _validate_n_diagnostics(n_diagnostic_rows, task_validation)
    _, runtime_attestation_rows, _ = read_tsv(
        root, output_dir_rel / "runtime_attestation.tsv", RUNTIME_ATTESTATION_FIELDS,
    )
    _validate_runtime_attestation(runtime_attestation_rows, task_validation)
    _, credible_rows, _ = read_tsv(
        root, output_dir_rel / "credible_sets.tsv", CREDIBLE_SET_FIELDS,
        allow_header_only=True,
    )
    credible = _validate_credible_sets(credible_rows, task_validation, diagnostics)
    del credible
    _, qc_rows, _ = read_tsv(root, output_dir_rel / "finemap_qc.tsv", FINEMAP_QC_FIELDS)
    if len(qc_rows) != 2:
        raise FineMappingError("fine-map QC does not retain both traits")
    for row in qc_rows:
        diagnostic = diagnostics.get(row["trait_role"])
        if diagnostic is None or any(
            row[field] != diagnostic[field] for field in (
                "analysis_id", "pair_id", "family_role", "locus_entry_id", "trait_id",
                "variant_count", "model_converged", "niter", "credible_set_count",
                "max_pip", "rss_ld_s", "kriging_allele_switch_outlier_count",
                "ld_symmetry_max_abs", "ld_diagonal_max_abs", "ld_minimum_eigenvalue",
                "reference_sample_size", "diagnostic_status", "error",
            )
        ):
            raise FineMappingError("fine-map QC differs from diagnostics")
        expected_status = "COMPLETE" if diagnostic["diagnostic_status"] == "PASS" else "FAILED_QC"
        if row["analysis_status"] != expected_status:
            raise FineMappingError("fine-map QC analysis status drifted")
    _, raw_coloc, _ = read_tsv(root, output_dir_rel / "coloc.raw.tsv", RAW_COLOC_FIELDS)
    _, engine_rows, _ = read_tsv(root, output_dir_rel / "engine_status.tsv", ENGINE_STATUS_FIELDS)
    if len(engine_rows) != 1 or engine_rows[0]["engine_status"] != "COMPLETE_SERIALIZED":
        raise FineMappingError("R engine did not serialize its terminal model state")
    canonical1 = _classify_variant_rows(
        raw1, "TRAIT1_SLEEP", task_validation, diagnostics["TRAIT1_SLEEP"], variant_fields,
    )
    canonical2 = _classify_variant_rows(
        raw2, "TRAIT2_EXTERNAL", task_validation, diagnostics["TRAIT2_EXTERNAL"], variant_fields,
    )
    canonical_coloc = classify_coloc_rows(raw_coloc, diagnostics, task_validation, coloc_fields)
    contents = {
        "10_finemap_trait1.fragment.tsv": tsv_bytes(variant_fields, canonical1),
        "11_finemap_trait2.fragment.tsv": tsv_bytes(variant_fields, canonical2),
        "12_trait_trait_coloc.fragment.tsv": tsv_bytes(coloc_fields, canonical_coloc),
    }
    output_dir = safe_path(root, output_dir_rel, "engine output directory", must_exist=True)
    for name, content in contents.items():
        _write_new_file(output_dir / name, content)
    fsync_directory(output_dir)
    return {
        "trait1_rows": canonical1, "trait2_rows": canonical2,
        "coloc_rows": canonical_coloc, "diagnostics": diagnostic_rows,
        "credible_sets": credible_rows, "finemap_qc": qc_rows,
        "fragment_sha256": {name: _sha256_content(content) for name, content in contents.items()},
    }


def validate_engine_bundle(
    root: Path, output_dir_rel: Path | str, task_validation: Mapping[str, Any],
) -> dict[str, Any]:
    """Side-effect-free validation of a fully adapted per-locus engine directory."""

    output_dir_rel = Path(output_dir_rel)
    policy = task_validation["policy"]
    verify_runtime_package_lock(
        root, task_validation["task"]["runtime_package_lock_path"],
        expected_sha256=task_validation["task"]["runtime_package_lock_sha256"],
    )
    expected = {
        "trait1.raw.tsv": RAW_VARIANT_FIELDS, "trait2.raw.tsv": RAW_VARIANT_FIELDS,
        "credible_sets.tsv": CREDIBLE_SET_FIELDS,
        "diagnostics.tsv": policy["diagnostic_contract"]["required_schema"],
        "sample_size_diagnostics.tsv": N_DIAGNOSTIC_FIELDS,
        "runtime_attestation.tsv": RUNTIME_ATTESTATION_FIELDS,
        "finemap_qc.tsv": FINEMAP_QC_FIELDS, "coloc.raw.tsv": RAW_COLOC_FIELDS,
        "engine_status.tsv": ENGINE_STATUS_FIELDS,
        "10_finemap_trait1.fragment.tsv": policy["required_science_outputs"]["output_10"]["schema"],
        "11_finemap_trait2.fragment.tsv": policy["required_science_outputs"]["output_11"]["schema"],
        "12_trait_trait_coloc.fragment.tsv": policy["required_science_outputs"]["output_12"]["schema"],
    }
    rows_by_name: dict[str, list[dict[str, str]]] = {}
    identities: dict[str, dict[str, object]] = {}
    for name, fields in expected.items():
        _, rows, identity = read_tsv(
            root, output_dir_rel / name, fields,
            allow_header_only=(name == "credible_sets.tsv"),
        )
        rows_by_name[name] = rows
        identities[name] = identity
    diagnostics = _validate_diagnostics(rows_by_name["diagnostics.tsv"], task_validation)
    _validate_n_diagnostics(rows_by_name["sample_size_diagnostics.tsv"], task_validation)
    _validate_runtime_attestation(rows_by_name["runtime_attestation.tsv"], task_validation)
    _validate_credible_sets(rows_by_name["credible_sets.tsv"], task_validation, diagnostics)
    expected1 = _classify_variant_rows(
        rows_by_name["trait1.raw.tsv"], "TRAIT1_SLEEP", task_validation,
        diagnostics["TRAIT1_SLEEP"], expected["10_finemap_trait1.fragment.tsv"],
    )
    expected2 = _classify_variant_rows(
        rows_by_name["trait2.raw.tsv"], "TRAIT2_EXTERNAL", task_validation,
        diagnostics["TRAIT2_EXTERNAL"], expected["11_finemap_trait2.fragment.tsv"],
    )
    expected12 = classify_coloc_rows(
        rows_by_name["coloc.raw.tsv"], diagnostics, task_validation,
        expected["12_trait_trait_coloc.fragment.tsv"],
    )
    if (
        rows_by_name["10_finemap_trait1.fragment.tsv"] != expected1
        or rows_by_name["11_finemap_trait2.fragment.tsv"] != expected2
        or rows_by_name["12_trait_trait_coloc.fragment.tsv"] != expected12
    ):
        raise FineMappingError("canonical per-locus fragments differ from deterministic adapter")
    return {"rows": rows_by_name, "identities": identities, "diagnostics": diagnostics}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--pre-manifest", default=str(PRE_MANIFEST_REL))
    parser.add_argument("--pre-lock", default=str(PRE_LOCK_REL))
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--print-adapter-schema", action="store_true")
    action.add_argument("--validate-pre-family", action="store_true")
    action.add_argument("--materialize-one", metavar="LOCUS_ENTRY_ID")
    action.add_argument("--freeze-manifest", action="store_true")
    action.add_argument("--verify-manifest", action="store_true")
    action.add_argument("--verify-task", metavar="TASK_PATH")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument(
        "--memory-envelope-bytes", type=int, default=BASE_PHYSICAL_ENVELOPE_BYTES,
    )
    parser.add_argument("--larger-host-continuation", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve(strict=True)
    if args.print_adapter_schema:
        print(json.dumps({
            "schema_version": PRE_FAMILY_SCHEMA,
            "pre_family_fields": PRE_FAMILY_FIELDS,
            "pre_unavailable_fields": PRE_UNAVAILABLE_FIELDS,
            "required_lock_keys": sorted({
                "schema_version", "analysis_id", "artifact_role", "state",
                "results_accessed_before_lock", "pair_order", "manifest", "locus_count",
                "locus_entry_ids_in_order", "evidence_count", "evidence_ids_sha256",
                "policy", "official_blocks", "upstream_family_provenance",
                "generator_scripts", "selection_contract", "zero_family_provenance",
                "unavailable_manifest", "unavailable_count", "unavailable_entry_ids_in_order",
                "unavailable_evidence_count", "unavailable_evidence_ids_sha256",
                "upstream_evidence_count", "upstream_evidence_ids_sha256",
            }),
            "scalar_n_convention": SCALAR_N_CONVENTION,
        }, indent=2, sort_keys=True))
        return 0
    if args.validate_pre_family:
        family = validate_pre_materialization_family(root, args.pre_manifest, args.pre_lock)
        print(json.dumps({
            "status": "PRE_MATERIALIZATION_FAMILY_VERIFIED", "state": family["state"],
            "mapped_loci": len(family["rows"]),
            "unavailable_loci": len(family["unavailable_rows"]),
            "upstream_evidence": len(family["evidence_ids"]) + len(family["unavailable_evidence_ids"]),
        }, sort_keys=True))
        return 0
    if args.materialize_one:
        if not args.execute:
            raise FineMappingError("explicit --execute is required for materialization")
        result = materialize_one_locus(
            root, args.materialize_one, manifest_rel=args.pre_manifest, lock_rel=args.pre_lock,
            configured_envelope_bytes=args.memory_envelope_bytes,
            larger_host_continuation=args.larger_host_continuation,
        )
        print(json.dumps({
            "status": "WHOLE_LOCUS_MATERIALIZED_AND_DEEP_VALIDATED",
            "locus_entry_id": args.materialize_one,
            "variant_count": result["manifest_row"]["variant_count"],
            "bundle_lock": result["bundle_lock_identity"],
        }, sort_keys=True))
        return 0
    if args.freeze_manifest:
        if not args.execute:
            raise FineMappingError("explicit --execute is required to freeze materialized manifest")
        result = freeze_materialized_manifest(
            root, pre_manifest_rel=args.pre_manifest, pre_lock_rel=args.pre_lock,
        )
        print(json.dumps({
            "status": "MATERIALIZED_MANIFEST_FROZEN", "state": result["state"],
            "locus_count": len(result.get("rows", [])),
        }, sort_keys=True))
        return 0
    if args.verify_manifest:
        result = validate_materialized_manifest(root)
        print(json.dumps({
            "status": "MATERIALIZED_MANIFEST_VERIFIED", "state": result["state"],
            "locus_count": len(result.get("rows", [])),
        }, sort_keys=True))
        return 0
    validated = deep_validate_locus_task(root, args.verify_task)
    print(json.dumps({
        "status": "LOCUS_TASK_DEEP_VALIDATED",
        "locus_entry_id": validated["task"]["locus_entry_id"],
        "variant_count": len(validated["order_rows"]),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ComputeBlocked as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(75) from error
    except FineMappingError as error:
        raise SystemExit(f"ERROR: {error}") from error
