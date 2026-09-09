#!/usr/bin/env python3
"""Freeze or verify revision 3 of the result-blind Track B mechanism contract.

Revision 3 preserves the two earlier frozen families byte-for-byte.  It splits
the previously shared QTL readiness fact into exact tissue-eQTL, tissue-sQTL,
cell-eQTL and cell-sQTL families, exposes the missing human brain-pericyte
resource, checksum-binds the complete interpretation code family, and requires
explicit executors for broad-tissue MAGMA, Allen-MTG broad-class aggregation,
cross-atlas replication and Pair-A-versus-Pair-B comparison.  This program
never reads future Track B fine-mapping, pleiotropy or mechanism results.
"""
from __future__ import annotations

import argparse
import ast
import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
import stat
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath
from typing import Callable, Mapping, Sequence

sys.dont_write_bytecode = True

ROOT_AT_IMPORT = Path(__file__).resolve().parents[1]
BASE_REL = "scripts/161_build_track_b_mechanism_readiness_v2.py"
BASE_SHA256 = "2fc37fcb150bd5c0425ab0416b4fe01830f63e5270685dabddb598994ed106a4"
BASE_POLICY_REL = "config/track_b_mechanism_followup_policy_v2.json"
BASE_POLICY_SHA256 = "d6adcd5ec9863ddef2cd7c36d6fbff95cd243ebe973bb0e0fa726c1817e70e81"
BASE_LOCK_REL = "results/track_b/mechanism_followup/v2/PRE_RESULT_READINESS.lock.json"
BASE_LOCK_SHA256 = "6643561657f69d13e6e5c8e0352d9c095a401943e4028585562f81e3de3ac530"


def _bootstrap_hash(path: Path) -> str:
    info = os.lstat(path)
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise RuntimeError(f"base implementation is not a real regular file: {path}")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(path, flags)
    try:
        before = os.fstat(fd)
        digest = hashlib.sha256()
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_size,
        value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before) != identity(after) or identity(after) != identity(final):
        raise RuntimeError("base implementation drifted during bootstrap")
    return digest.hexdigest()


if _bootstrap_hash(ROOT_AT_IMPORT / BASE_REL) != BASE_SHA256:
    raise RuntimeError("revision-2 implementation differs from the revision-3 pin")
_BASE_SPEC = importlib.util.spec_from_file_location(
    "track_b_mechanism_readiness_v2_base", ROOT_AT_IMPORT / BASE_REL,
)
if _BASE_SPEC is None or _BASE_SPEC.loader is None:
    raise RuntimeError("cannot load revision-2 base implementation")
M2 = importlib.util.module_from_spec(_BASE_SPEC)
sys.modules[_BASE_SPEC.name] = M2
_BASE_SPEC.loader.exec_module(M2)
B = M2.B

SCHEMA = "track-b-mechanism-pre-result-readiness.3"
POLICY_REL = "config/track_b_mechanism_followup_policy_v3.json"
SCRIPT_REL = "scripts/162_build_track_b_mechanism_readiness_v3.py"
OUTPUT_DIR_REL = "results/track_b/mechanism_followup/v3"
OUTPUT_RELS = (
    f"{OUTPUT_DIR_REL}/PRE_RESULT_RESOURCE_EVIDENCE.tsv",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.tsv",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.json",
    f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS_REPORT.md",
)
SUPERSESSION_REL = "results/track_b/mechanism_followup/PRE_RESULT_READINESS_SUPERSESSION_V3.json"
LOCK_REL = f"{OUTPUT_DIR_REL}/PRE_RESULT_READINESS.lock.json"
EXPECTED_ANALYSIS_IDS = (
    "P07_CONJFDR", "P09_HDL_L", "P09_RHO_HESS", "P13_MAGMA_BROAD_TISSUE",
    "P13_SLDSC_GTEX", "P13_GTEX_TISSUE_EQTL", "P13_GTEX_TISSUE_SQTL",
    "P14_BRAIN_CELL_CLASS_DISCOVERY", "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY",
    "P15_BRAIN_CELL_SUBTYPE_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION",
    "P16_HUMAN_CELL_SUBTYPE_REPLICATION", "P17_CELL_EQTL", "P18_CELL_SQTL",
    "P19_THREE_WAY_EQTL", "P19_THREE_WAY_SQTL", "P20_CATLAS_SCATAC",
    "P20_SCREEN_CHROMATIN", "P20_HOCOMOCO_MOTIF", "P21_ABC_ENHANCER_GENE",
    "P21_PCHIC_ENHANCER_GENE", "P21_TWAS", "P22_REGULATORY_CHAIN_EQTL",
    "P22_REGULATORY_CHAIN_SQTL", "P23_SPATIAL", "P24_LOCKED_ATLAS_CROSS_SLEEP",
    "P25_PAIR_A_VS_B", "P26_POSITIVE_CONTROL", "P27_PATHWAYS", "P28_GRN",
    "P29_BIDIRECTIONAL_MR",
)
EXPECTED_OVERLAY_KEYS = {
    "schema_version", "contract_revision", "base_policy_path", "base_policy_sha256",
    "base_implementation_path", "base_implementation_sha256", "supersedes",
    "new_absence_probes", "expected_analysis_ids", "replacement_rows",
    "interpretation_code_sha256", "interpretation_code_analysis_ids",
    "brain_class_override", "broad_tissue_override", "sldsc_override", "cell_qtl_overrides",
    "cell_replication_overrides", "pair_comparison_override",
    "result_blind", "future_results_accessed",
}
NEW_PROBE_KEYS = {
    "brain_class_mapping_executor", "broad_tissue_magma_executor",
    "cell_replication_executor", "gtex_eqtl_payload_root",
    "gtex_sqtl_source_manifest", "human_brain_pericyte_reference",
    "pair_comparison_executor", "sldsc_broad_tissue_coverage",
}
FACT_TYPES = dict(B.FACT_TYPES)
del FACT_TYPES["qtl_metadata"]
del FACT_TYPES["qtl_payloads"]
FACT_TYPES.update({
    "brain_class_mapping_executor": "SOFTWARE",
    "broad_tissue_magma_executor": "SOFTWARE",
    "cell_eqtl_metadata": "DATA",
    "cell_eqtl_payloads": "DATA",
    "cell_sqtl_metadata": "DATA",
    "cell_sqtl_payloads": "DATA",
    "cell_replication_executor": "SOFTWARE",
    "gtex_eqtl_metadata": "DATA",
    "gtex_eqtl_payloads": "DATA",
    "gtex_sqtl_payloads": "DATA",
    "human_brain_pericyte_reference": "DATA",
    "interpretation_code_family": "SOFTWARE",
    "pair_comparison_executor": "SOFTWARE",
    "sldsc_broad_tissue_coverage": "DATA",
})
INTERPRETATION_CODE_PATHS = (
    "scripts/downstream_contract.py", "scripts/liftover_chain.py",
    "scripts/pathway_sources.py", "scripts/catlas_policy.py",
    "scripts/74_interpretation_preflight.py", "scripts/75_prepare_interpretation_tasks.py",
    "scripts/76_run_interpretation_task.py", "scripts/76_record_interpretation_task.py",
    "scripts/77_collate_interpretation.py", "scripts/78_build_atlas_edges.py",
    "scripts/82_run_regulatory_task.py", "scripts/83_run_motif_task.py",
    "scripts/84_prepare_abc_overlap_cache.py", "scripts/85_run_abc_task.py",
    "scripts/86_prepare_pchic_overlap_cache.py", "scripts/87_run_pchic_task.py",
    "scripts/88_materialize_fuma_resources.py", "scripts/89_prepare_magma_annotation.py",
    "scripts/90_prepare_magma_gene_results.py", "scripts/91_run_fuma_scrna_task.py",
    "scripts/92_run_pathway_task.py", "scripts/93_setup_causal_runtime.sh",
    "scripts/94_prepare_causal_reference.py", "scripts/95_prepare_catlas_reference.py",
    "scripts/96_prepare_catlas_trait_cache.py", "scripts/97_run_catlas_task.py",
    "scripts/98_prepare_ldsc_seg_gtex_source.py",
    "scripts/98_prepare_ldsc_seg_reference.py",
    "scripts/98_prepare_ldsc_seg_trait.py", "scripts/98_run_ldsc_seg_task.py",
)
INTERPRETATION_CODE_ANALYSIS_IDS = (
    "P13_MAGMA_BROAD_TISSUE", "P13_SLDSC_GTEX",
    "P14_BRAIN_CELL_CLASS_DISCOVERY", "P15_BRAIN_CELL_SUBTYPE_DISCOVERY",
    "P16_HUMAN_CELL_CLASS_REPLICATION", "P16_HUMAN_CELL_SUBTYPE_REPLICATION",
    "P20_CATLAS_SCATAC", "P20_SCREEN_CHROMATIN", "P20_HOCOMOCO_MOTIF",
    "P21_ABC_ENHANCER_GENE", "P21_PCHIC_ENHANCER_GENE", "P27_PATHWAYS",
)
EXPECTED_REQUIREMENTS = {
    "P07_CONJFDR": ["conjfdr_runtime"],
    "P09_HDL_L": ["hdl_l_reference"],
    "P09_RHO_HESS": ["rho_hess_reference"],
    "P13_MAGMA_BROAD_TISSUE": ["fuma_bundle", "broad_tissue_magma", "broad_tissue_magma_executor", "interpretation_code_family"],
    "P13_SLDSC_GTEX": ["ldsc_source_bundle", "sldsc_broad_tissue_coverage", "ldsc_runtime", "ldsc_static_reference", "interpretation_code_family"],
    "P13_GTEX_TISSUE_EQTL": ["gtex_eqtl_metadata", "gtex_eqtl_payloads", "tabix"],
    "P13_GTEX_TISSUE_SQTL": ["gtex_sqtl_payloads", "tabix"],
    "P14_BRAIN_CELL_CLASS_DISCOVERY": ["fuma_bundle", "brain_discovery_human", "brain_class_mapping_executor", "interpretation_code_family"],
    "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY": ["human_brain_pericyte_reference"],
    "P15_BRAIN_CELL_SUBTYPE_DISCOVERY": ["fuma_bundle", "brain_discovery_human", "interpretation_code_family"],
    "P16_HUMAN_CELL_CLASS_REPLICATION": ["fuma_bundle", "brain_replication_human_class", "brain_class_mapping_executor", "cell_replication_executor", "interpretation_code_family"],
    "P16_HUMAN_CELL_SUBTYPE_REPLICATION": ["human_subtype_replication", "cell_replication_executor", "interpretation_code_family"],
    "P17_CELL_EQTL": ["cell_eqtl_metadata", "cell_eqtl_payloads", "tabix"],
    "P18_CELL_SQTL": ["cell_sqtl_metadata", "cell_sqtl_payloads", "tabix"],
    "P19_THREE_WAY_EQTL": ["cell_eqtl_metadata", "cell_eqtl_payloads", "tabix", "three_way_coloc_executor"],
    "P19_THREE_WAY_SQTL": ["cell_sqtl_metadata", "cell_sqtl_payloads", "tabix", "three_way_coloc_executor"],
    "P20_CATLAS_SCATAC": ["catlas_bundle", "interpretation_code_family"],
    "P20_SCREEN_CHROMATIN": ["screen_bundle", "interpretation_code_family"],
    "P20_HOCOMOCO_MOTIF": ["hocomoco_bundle", "hocomoco_sequence", "interpretation_code_family"],
    "P21_ABC_ENHANCER_GENE": ["abc_bundle", "interpretation_code_family"],
    "P21_PCHIC_ENHANCER_GENE": ["pchic_bundle", "interpretation_code_family"],
    "P21_TWAS": ["twas_models", "twas_runtime"],
    "P22_REGULATORY_CHAIN_EQTL": ["screen_bundle", "catlas_bundle", "abc_bundle", "pchic_bundle", "cell_eqtl_metadata", "cell_eqtl_payloads", "regulatory_chain_executor"],
    "P22_REGULATORY_CHAIN_SQTL": ["screen_bundle", "catlas_bundle", "abc_bundle", "pchic_bundle", "cell_sqtl_metadata", "cell_sqtl_payloads", "regulatory_chain_executor"],
    "P23_SPATIAL": ["spatial"],
    "P24_LOCKED_ATLAS_CROSS_SLEEP": ["analysis_panel", "dense_sleep_family", "cross_sleep_executor"],
    "P25_PAIR_A_VS_B": ["pair_manifest", "pair_comparison_executor"],
    "P26_POSITIVE_CONTROL": ["pair_manifest"],
    "P27_PATHWAYS": ["pathway_bundle", "pathway_executor", "interpretation_code_family"],
    "P28_GRN": ["matched_multimodal_grn", "grn_executor"],
    "P29_BIDIRECTIONAL_MR": ["causal_bundle", "track_b_mr_executor"],
}
EXPECTED_UNBLOCK_TERMS = {
    "P07_CONJFDR": ("conjfdr", "runtime", "reference", "executor"),
    "P09_HDL_L": ("hdl-l", "implementation", "reference"),
    "P09_RHO_HESS": ("rho-hess", "contract", "executor", "reference"),
    "P13_MAGMA_BROAD_TISSUE": ("broad-tissue", "lung/airway", "skeletal muscle", "adapter", "correction family"),
    "P13_SLDSC_GTEX": ("lung index-36", "muscle_skeletal index-38", "18-tissue manifest", "static reference", "airway-specific"),
    "P13_GTEX_TISSUE_EQTL": ("all 49", "payloads", "tabix"),
    "P13_GTEX_TISSUE_SQTL": ("complete gtex v8 tissue-sqtl", "payloads", "indexes", "tabix"),
    "P14_BRAIN_CELL_CLASS_DISCOVERY": ("75-subtype-to-seven-class", "mapping/aggregation executor", "upstream tissue gate"),
    "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY": ("human brain pericyte", "atlas", "executor/adapter", "before upstream"),
    "P15_BRAIN_CELL_SUBTYPE_DISCOVERY": ("upstream discovery-class gate", "frozen subtype family"),
    "P16_HUMAN_CELL_CLASS_REPLICATION": ("allen-mtg broad-class mapping", "cross-atlas class-alignment/replication-classification executors", "upstream discovery"),
    "P16_HUMAN_CELL_SUBTYPE_REPLICATION": ("independent human subtype-resolved brain atlas", "dataset adapter", "cross-atlas subtype-alignment/replication classifier", "before discovery"),
    "P17_CELL_EQTL": ("qtd000559", "qtd000569", ".tbi indexes", "tabix"),
    "P18_CELL_SQTL": ("qtd000563", "qtd000573", "leafcutter", ".tbi indexes", "tabix"),
    "P19_THREE_WAY_EQTL": ("both cell-eqtl payloads", "multi-signal three-way executor", "before results"),
    "P19_THREE_WAY_SQTL": ("both leafcutter cell-sqtl payloads", "multi-signal three-way executor", "before results"),
    "P20_CATLAS_SCATAC": ("upstream variant/trait gates",),
    "P20_SCREEN_CHROMATIN": ("upstream fine-mapped/shared-signal variant gate",),
    "P20_HOCOMOCO_MOTIF": ("checksum-pinned hg38 fasta/index", "validated scorer"),
    "P21_ABC_ENHANCER_GENE": ("upstream eligible variants", "molecularly supported same-locus genes"),
    "P21_PCHIC_ENHANCER_GENE": ("upstream eligible variants", "molecularly supported same-locus genes"),
    "P21_TWAS": ("all 98", ".molecular-env/python runtime"),
    "P22_REGULATORY_CHAIN_EQTL": ("both cell-eqtl payloads/indexes", "regulatory-chain integrator", "upstream signals"),
    "P22_REGULATORY_CHAIN_SQTL": ("both leafcutter cell-sqtl payloads/indexes", "regulatory-chain integrator", "upstream signals"),
    "P23_SPATIAL": ("authorized immutable human spatial dataset", "metadata", "reference build", "validated executor", "before upstream"),
    "P24_LOCKED_ATLAS_CROSS_SLEEP": ("cross-sleep integration executor", "upstream mechanism targets"),
    "P25_PAIR_A_VS_B": ("both pairs", "same preregistered upstream/mechanism families", "comparison executor"),
    "P26_POSITIVE_CONTROL": ("control", "homologous upstream and mechanism analyses"),
    "P27_PATHWAYS": ("upstream high-confidence convergent-gene lists", "without changing the resource family"),
    "P28_GRN": ("matched multimodal human data", "grn executor", "before upstream"),
    "P29_BIDIRECTIONAL_MR": ("automatic executor", "instrument", "overlap", "steiger", "pleiotropy", "family-wide correction"),
}
INTERPRETATION_ENTRYPOINTS = tuple(
    relative for relative in INTERPRETATION_CODE_PATHS
    if relative not in {
        "scripts/downstream_contract.py", "scripts/liftover_chain.py",
        "scripts/pathway_sources.py", "scripts/catlas_policy.py",
    }
)
GTEX_FIELDS = [
    "study", "qtl_group", "tissue_ontology_id", "tissue_ontology_term",
    "tissue_label", "condition_label", "quant_method", "ftp_path",
]
CELL_METADATA_RELS = (
    "ref/molecular/eqtl_catalogue_r7/dataset_metadata_r7.tsv",
    "ref/molecular/eqtl_catalogue_r7/tabix_ftp_paths.tsv",
)
GTEX_METADATA_REL = "ref/molecular/eqtl_catalogue_r7/tabix_ftp_paths_imported.tsv"
QTL_ROOT_REL = "ref/molecular/eqtl_catalogue_r7/payloads"
CELL_EQTL_IDS = ("QTD000559", "QTD000569")
CELL_SQTL_IDS = ("QTD000563", "QTD000573")
SAFE_QTL_GROUP = re.compile(r"[A-Za-z0-9_-]+")


def _validate_old_lock(
    root: Path,
    relative: str,
    digest: str,
    schema: str,
    output_rels: Sequence[str],
    expected_policy_sha: str,
    expected_script_sha: str,
) -> dict[str, object]:
    payload = B.read_stable_bytes(root, relative)
    if hashlib.sha256(payload).hexdigest() != digest:
        raise B.ContractError(f"superseded lock bytes differ from revision-3 lineage: {relative}")
    try:
        lock = json.loads(payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise B.ContractError(f"superseded lock is invalid JSON: {relative}: {exc}") from exc
    if (
        not isinstance(lock, dict) or B.canonical_json(lock) != payload
        or lock.get("schema_version") != schema
        or lock.get("contract_kind") != "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK"
        or lock.get("result_blind") is not True
        or lock.get("future_results_accessed") is not False
        or lock.get("policy_sha256") != expected_policy_sha
        or lock.get("script_sha256") != expected_script_sha
    ):
        raise B.ContractError(f"superseded lock has invalid canonical lineage: {relative}")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != set(output_rels):
        raise B.ContractError(f"superseded lock output family is incomplete: {relative}")
    for output_rel in output_rels:
        expected = outputs[output_rel]
        if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
            raise B.ContractError(f"superseded output identity is malformed: {output_rel}")
        item = B.stable_file(root, output_rel)
        if item.observed_bytes != expected["bytes"] or item.sha256 != expected["sha256"]:
            raise B.ContractError(f"superseded output is not byte-exact: {output_rel}")
    return lock


def validate_superseded_families(root: Path, overlay: Mapping[str, object]) -> dict[str, object]:
    supersedes = overlay.get("supersedes")
    expected_keys = {"lock_path", "lock_sha256", "retention", "reason"}
    if not isinstance(supersedes, dict) or set(supersedes) != expected_keys:
        raise B.ContractError("revision-3 supersession record is malformed")
    if (
        supersedes.get("lock_path") != BASE_LOCK_REL
        or supersedes.get("lock_sha256") != BASE_LOCK_SHA256
        or supersedes.get("retention")
        != "PRESERVE_V1_AND_V2_BYTE_EXACT_AS_AUDIT_TRAILS_DO_NOT_USE_FOR_EXECUTION"
    ):
        raise B.ContractError("revision-3 supersession identity/retention differs from its exact pin")
    v2_overlay = B.read_json(root, BASE_POLICY_REL)
    M2.validate_superseded_family(root, v2_overlay)
    v2_outputs = tuple(M2.OUTPUT_RELS) + (M2.SUPERSESSION_REL,)
    lock = _validate_old_lock(
        root, BASE_LOCK_REL, BASE_LOCK_SHA256, M2.SCHEMA, v2_outputs,
        BASE_POLICY_SHA256, BASE_SHA256,
    )
    receipt_payload = B.read_stable_bytes(root, M2.SUPERSESSION_REL)
    try:
        receipt = json.loads(receipt_payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise B.ContractError(f"revision-2 supersession receipt is invalid: {exc}") from exc
    replacement = receipt.get("replacement_contract") if isinstance(receipt, dict) else None
    if (
        B.canonical_json(receipt) != receipt_payload or not isinstance(replacement, dict)
        or replacement.get("policy_path") != BASE_POLICY_REL
        or replacement.get("policy_sha256") != BASE_POLICY_SHA256
        or replacement.get("script_path") != BASE_REL
        or replacement.get("script_sha256") != BASE_SHA256
        or replacement.get("lock_path") != BASE_LOCK_REL
    ):
        raise B.ContractError("revision-2 supersession receipt has invalid replacement lineage")
    return lock


def _validate_analysis_row(row: object) -> dict[str, object]:
    if not isinstance(row, dict) or set(row) != B.ANALYSIS_REQUIRED_KEYS:
        raise B.ContractError("revision-3 replacement analysis has an unexpected schema")
    result = copy.deepcopy(row)
    requirements = result.get("requirements")
    if (
        not isinstance(requirements, list) or not requirements
        or len(requirements) != len(set(requirements))
        or any(value not in FACT_TYPES for value in requirements)
    ):
        raise B.ContractError(f"revision-3 analysis has invalid requirements: {result.get('analysis_id')}")
    for key in B.ANALYSIS_REQUIRED_KEYS - {"requirements"}:
        if not isinstance(result.get(key), str) or not str(result[key]).strip():
            raise B.ContractError(f"revision-3 analysis has invalid field {key}: {result.get('analysis_id')}")
    output = str(result["output"])
    pure = PurePosixPath(output)
    if (
        not output.startswith("results/track_b/mechanism_followup/science/")
        or pure.is_absolute() or "\\" in output or any(part in {"", ".", ".."} for part in pure.parts)
    ):
        raise B.ContractError(f"revision-3 planned output escapes the future science namespace: {result['analysis_id']}")
    return result


def validate_effective_policy(policy: Mapping[str, object]) -> None:
    if policy.get("analysis_build") != "GRCh37" or policy.get("analysis_ancestry") != "European":
        raise B.ContractError("effective revision-3 build/ancestry must remain GRCh37/European")
    if policy.get("forbidden_inputs") != list(B.FORBIDDEN_RESULT_PREFIXES):
        raise B.ContractError("effective revision-3 result-blind boundary changed")
    if policy.get("expected_analysis_ids") != list(EXPECTED_ANALYSIS_IDS):
        raise B.ContractError("effective revision-3 expected analysis family is incomplete")
    analyses = policy.get("analyses")
    if not isinstance(analyses, list):
        raise B.ContractError("effective revision-3 analyses are absent")
    validated = [_validate_analysis_row(row) for row in analyses]
    identifiers = [row["analysis_id"] for row in validated]
    if identifiers != list(EXPECTED_ANALYSIS_IDS) or len(set(identifiers)) != len(identifiers):
        raise B.ContractError("effective revision-3 analysis family is duplicated, reordered or incomplete")
    if tuple(EXPECTED_REQUIREMENTS) != EXPECTED_ANALYSIS_IDS or set(EXPECTED_UNBLOCK_TERMS) != set(EXPECTED_ANALYSIS_IDS):
        raise B.ContractError("hard-coded revision-3 requirement/unblock audit family is incomplete or reordered")
    for row in validated:
        identity = str(row["analysis_id"])
        if row["requirements"] != EXPECTED_REQUIREMENTS[identity]:
            raise B.ContractError(f"effective revision-3 exact requirement family drifted: {identity}")
        unblock = str(row["unblock"]).casefold()
        missing_terms = [term for term in EXPECTED_UNBLOCK_TERMS[identity] if term.casefold() not in unblock]
        if missing_terms:
            raise B.ContractError(
                f"effective revision-3 unblock predicate is not tied to its exact units: "
                f"{identity}: {missing_terms}"
            )
    outputs = [row["output"] for row in validated]
    if len(outputs) != len(set(outputs)):
        raise B.ContractError("effective revision-3 planned output family is not unique")
    if policy.get("required_qtl_datasets") != {
        "eQTL": list(CELL_EQTL_IDS), "sQTL": list(CELL_SQTL_IDS),
    }:
        raise B.ContractError("effective revision-3 cell-QTL family differs from the leafcutter policy")


def load_effective_policy(
    root: Path,
) -> tuple[dict[str, object], dict[str, object], str, dict[str, object]]:
    overlay = B.read_json(root, POLICY_REL)
    if (
        set(overlay) != EXPECTED_OVERLAY_KEYS or overlay.get("schema_version") != SCHEMA
        or overlay.get("contract_revision") != 3
        or overlay.get("result_blind") is not True
        or overlay.get("future_results_accessed") is not False
    ):
        raise B.ContractError("revision-3 overlay schema/result-blind boundary is invalid")
    if (
        overlay.get("base_policy_path") != BASE_POLICY_REL
        or overlay.get("base_policy_sha256") != BASE_POLICY_SHA256
        or overlay.get("base_implementation_path") != BASE_REL
        or overlay.get("base_implementation_sha256") != BASE_SHA256
    ):
        raise B.ContractError("revision-3 base lineage differs from its exact pins")
    if (
        B.stable_file(root, BASE_POLICY_REL).sha256 != BASE_POLICY_SHA256
        or B.stable_file(root, BASE_REL).sha256 != BASE_SHA256
    ):
        raise B.ContractError("revision-2 policy/implementation drifted")
    v2_lock = validate_superseded_families(root, overlay)
    _v2_overlay, effective, v2_effective_hash = M2.load_effective_policy(root)
    if v2_lock.get("effective_policy_sha256") != v2_effective_hash:
        raise B.ContractError("revision-2 effective-policy hash differs from its frozen lock")

    probes = overlay.get("new_absence_probes")
    if not isinstance(probes, dict) or set(probes) != NEW_PROBE_KEYS:
        raise B.ContractError("revision-3 absence-probe family is incomplete")
    for relative in probes.values():
        B.safe_relative(str(relative))
    code_hashes = overlay.get("interpretation_code_sha256")
    if not isinstance(code_hashes, dict) or set(code_hashes) != set(INTERPRETATION_CODE_PATHS):
        raise B.ContractError("revision-3 interpretation code family is incomplete")
    for relative in INTERPRETATION_CODE_PATHS:
        digest = code_hashes.get(relative)
        B.safe_relative(relative)
        if not isinstance(digest, str) or not B.SHA256_RE.fullmatch(digest):
            raise B.ContractError(f"revision-3 interpretation code hash is invalid: {relative}")
    if overlay.get("interpretation_code_analysis_ids") != list(INTERPRETATION_CODE_ANALYSIS_IDS):
        raise B.ContractError("revision-3 interpretation code consumers are incomplete or reordered")
    if overlay.get("expected_analysis_ids") != list(EXPECTED_ANALYSIS_IDS):
        raise B.ContractError("revision-3 overlay expected analysis family is incomplete")

    replacements = overlay.get("replacement_rows")
    if not isinstance(replacements, list):
        raise B.ContractError("revision-3 replacement rows are absent")
    replacements_by_id = {
        str(row.get("analysis_id")): _validate_analysis_row(row) for row in replacements
        if isinstance(row, dict)
    }
    replacement_order = [str(row.get("analysis_id")) for row in replacements if isinstance(row, dict)]
    expected_replacements = [
        "P13_GTEX_TISSUE_EQTL", "P13_GTEX_TISSUE_SQTL",
        "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY", "P19_THREE_WAY_EQTL",
        "P19_THREE_WAY_SQTL", "P22_REGULATORY_CHAIN_EQTL",
        "P22_REGULATORY_CHAIN_SQTL",
    ]
    if replacement_order != expected_replacements or len(replacements_by_id) != len(expected_replacements):
        raise B.ContractError("revision-3 replacement family is duplicated, reordered or incomplete")

    brain_override = overlay.get("brain_class_override")
    if (
        not isinstance(brain_override, dict)
        or set(brain_override) != {
            "analysis_id", "scope", "requirements", "required_inputs",
            "required_software", "unblock", "claim_limit",
        }
        or brain_override.get("analysis_id") != "P14_BRAIN_CELL_CLASS_DISCOVERY"
        or brain_override.get("requirements")
        != ["fuma_bundle", "brain_discovery_human", "brain_class_mapping_executor"]
        or any(
            not isinstance(brain_override.get(key), str) or not str(brain_override[key]).strip()
            for key in ("scope", "required_inputs", "required_software", "unblock", "claim_limit")
        )
    ):
        raise B.ContractError("revision-3 brain-class override is invalid")
    broad_override = overlay.get("broad_tissue_override")
    if (
        not isinstance(broad_override, dict)
        or set(broad_override) != {"analysis_id", "requirements", "unblock"}
        or broad_override.get("analysis_id") != "P13_MAGMA_BROAD_TISSUE"
        or broad_override.get("requirements") != [
            "fuma_bundle", "broad_tissue_magma", "broad_tissue_magma_executor",
            "interpretation_code_family",
        ]
        or not isinstance(broad_override.get("unblock"), str)
        or not str(broad_override["unblock"]).strip()
    ):
        raise B.ContractError("revision-3 broad-tissue override is invalid")
    sldsc_override = overlay.get("sldsc_override")
    if (
        not isinstance(sldsc_override, dict)
        or set(sldsc_override) != {
            "analysis_id", "scope", "requirements", "required_inputs",
            "ram_decomposition", "unblock", "claim_limit",
        }
        or sldsc_override.get("analysis_id") != "P13_SLDSC_GTEX"
        or sldsc_override.get("requirements") != [
            "ldsc_source_bundle", "sldsc_broad_tissue_coverage", "ldsc_runtime",
            "ldsc_static_reference", "interpretation_code_family",
        ]
        or any(
            not isinstance(sldsc_override.get(key), str) or not str(sldsc_override[key]).strip()
            for key in ("scope", "required_inputs", "ram_decomposition", "unblock", "claim_limit")
        )
    ):
        raise B.ContractError("revision-3 S-LDSC coverage override is invalid")
    cell_overrides = overlay.get("cell_qtl_overrides")
    if not isinstance(cell_overrides, list) or len(cell_overrides) != 2:
        raise B.ContractError("revision-3 cell-QTL overrides are invalid")
    cell_by_id: dict[str, dict[str, object]] = {}
    for row in cell_overrides:
        if not isinstance(row, dict) or set(row) != {"analysis_id", "requirements", "unblock"}:
            raise B.ContractError("revision-3 cell-QTL override schema is invalid")
        if not isinstance(row.get("unblock"), str) or not str(row["unblock"]).strip():
            raise B.ContractError("revision-3 cell-QTL unblock condition is invalid")
        cell_by_id[str(row["analysis_id"])] = copy.deepcopy(row)
    if set(cell_by_id) != {"P17_CELL_EQTL", "P18_CELL_SQTL"}:
        raise B.ContractError("revision-3 cell-QTL override family is incomplete")
    if cell_by_id["P17_CELL_EQTL"]["requirements"] != ["cell_eqtl_metadata", "cell_eqtl_payloads", "tabix"]:
        raise B.ContractError("revision-3 P17 requirements are not exactly cell-eQTL")
    if cell_by_id["P18_CELL_SQTL"]["requirements"] != ["cell_sqtl_metadata", "cell_sqtl_payloads", "tabix"]:
        raise B.ContractError("revision-3 P18 requirements are not exactly leafcutter cell-sQTL")
    replication_overrides = overlay.get("cell_replication_overrides")
    if not isinstance(replication_overrides, list) or len(replication_overrides) != 2:
        raise B.ContractError("revision-3 cell-replication overrides are invalid")
    replication_by_id: dict[str, dict[str, object]] = {}
    for row in replication_overrides:
        if (
            not isinstance(row, dict) or set(row) != {"analysis_id", "requirements", "unblock"}
            or not isinstance(row.get("unblock"), str) or not str(row["unblock"]).strip()
        ):
            raise B.ContractError("revision-3 cell-replication override schema is invalid")
        replication_by_id[str(row["analysis_id"])] = copy.deepcopy(row)
    expected_replication_requirements = {
        "P16_HUMAN_CELL_CLASS_REPLICATION": [
            "fuma_bundle", "brain_replication_human_class",
            "brain_class_mapping_executor", "cell_replication_executor",
            "interpretation_code_family",
        ],
        "P16_HUMAN_CELL_SUBTYPE_REPLICATION": [
            "human_subtype_replication", "cell_replication_executor",
            "interpretation_code_family",
        ],
    }
    if set(replication_by_id) != set(expected_replication_requirements) or any(
        replication_by_id[identity]["requirements"] != requirements
        for identity, requirements in expected_replication_requirements.items()
    ):
        raise B.ContractError("revision-3 cell-replication requirement family is invalid")
    pair_override = overlay.get("pair_comparison_override")
    if (
        not isinstance(pair_override, dict)
        or set(pair_override) != {"analysis_id", "requirements", "unblock"}
        or pair_override.get("analysis_id") != "P25_PAIR_A_VS_B"
        or pair_override.get("requirements") != ["pair_manifest", "pair_comparison_executor"]
        or not isinstance(pair_override.get("unblock"), str)
        or not str(pair_override["unblock"]).strip()
    ):
        raise B.ContractError("revision-3 pair-comparison override is invalid")

    transformed: list[dict[str, object]] = []
    for source in effective["analyses"]:
        row = copy.deepcopy(source)
        identity = str(row["analysis_id"])
        if identity == "P13_MAGMA_BROAD_TISSUE":
            row.update({key: copy.deepcopy(value) for key, value in broad_override.items() if key != "analysis_id"})
            transformed.append(row)
        elif identity == "P13_SLDSC_GTEX":
            row.update({key: copy.deepcopy(value) for key, value in sldsc_override.items() if key != "analysis_id"})
            transformed.append(row)
        elif identity == "P13_GTEX_TISSUE_QTL":
            transformed.extend(copy.deepcopy(replacements_by_id[name]) for name in expected_replacements[:2])
        elif identity == "P14_BRAIN_CELL_CLASS_DISCOVERY":
            row.update({key: copy.deepcopy(value) for key, value in brain_override.items() if key != "analysis_id"})
            transformed.append(row)
            transformed.append(copy.deepcopy(replacements_by_id["P14_HUMAN_BRAIN_PERICYTE_DISCOVERY"]))
        elif identity == "P19_THREE_WAY_COLOC":
            transformed.extend(copy.deepcopy(replacements_by_id[name]) for name in expected_replacements[3:5])
        elif identity == "P22_REGULATORY_CHAIN":
            transformed.extend(copy.deepcopy(replacements_by_id[name]) for name in expected_replacements[5:7])
        elif identity in cell_by_id:
            row.update({key: copy.deepcopy(value) for key, value in cell_by_id[identity].items() if key != "analysis_id"})
            transformed.append(row)
        elif identity in replication_by_id:
            row.update({key: copy.deepcopy(value) for key, value in replication_by_id[identity].items() if key != "analysis_id"})
            transformed.append(row)
        elif identity == "P25_PAIR_A_VS_B":
            row.update({key: copy.deepcopy(value) for key, value in pair_override.items() if key != "analysis_id"})
            transformed.append(row)
        else:
            transformed.append(row)
    for row in transformed:
        if row["analysis_id"] in INTERPRETATION_CODE_ANALYSIS_IDS:
            requirements = row["requirements"]
            assert isinstance(requirements, list)
            if "interpretation_code_family" not in requirements:
                requirements.append("interpretation_code_family")
            marker = "checksum-pinned complete interpretation dispatcher/runner/recorder/BH-collation code family"
            if marker not in str(row["required_software"]):
                row["required_software"] = f"{row['required_software']}; {marker}"
    effective["expected_analysis_ids"] = list(EXPECTED_ANALYSIS_IDS)
    effective["analyses"] = transformed
    validate_effective_policy(effective)
    return overlay, effective, B.bytes_sha256(B.canonical_json(effective)), v2_lock


def validate_gtex_contexts(root: Path) -> list[str]:
    fields, rows = B.read_tsv(root, GTEX_METADATA_REL)
    if fields != GTEX_FIELDS or len(rows) != 49:
        raise B.ContractError("GTEx imported eQTL metadata is not the exact 49-context schema/family")
    groups: list[str] = []
    for row in rows:
        group = row.get("qtl_group", "")
        expected_ftp = f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/{group}.tsv.gz"
        if (
            row.get("study") != "GTEx_V8" or row.get("quant_method") != "ge"
            or not SAFE_QTL_GROUP.fullmatch(group) or row.get("ftp_path") != expected_ftp
            or any(not row.get(key, "").strip() for key in GTEX_FIELDS[:-1])
        ):
            raise B.ContractError(f"GTEx imported context has wrong release/method/path/schema: {group!r}")
        groups.append(group)
    if len(groups) != len(set(groups)):
        raise B.ContractError("GTEx imported eQTL context family contains duplicates")
    return groups


def validate_interpretation_code_semantics(root: Path) -> dict[str, object]:
    parsed: dict[str, ast.Module] = {}
    for relative in INTERPRETATION_CODE_PATHS:
        payload = B.read_stable_bytes(root, relative)
        if relative.endswith(".sh"):
            if not payload.startswith(b"#!/usr/bin/env bash\nset -euo pipefail\n"):
                raise B.ContractError(f"interpretation shell entrypoint lacks fail-fast grammar: {relative}")
            continue
        try:
            parsed[relative] = ast.parse(payload.decode("utf-8"), filename=relative)
        except (UnicodeError, SyntaxError) as exc:
            raise B.ContractError(f"interpretation Python source is invalid: {relative}: {exc}") from exc
    for relative in INTERPRETATION_ENTRYPOINTS:
        if relative.endswith(".sh"):
            continue
        functions = {
            node.name for node in parsed[relative].body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        if "main" not in functions:
            raise B.ContractError(f"interpretation production entrypoint has no main function: {relative}")
    downstream = parsed["scripts/downstream_contract.py"]
    declared: object | None = None
    for node in downstream.body:
        if (
            isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "INTERPRETATION_SCRIPTS" for target in node.targets)
        ):
            try:
                declared = ast.literal_eval(node.value)
            except (ValueError, TypeError) as exc:
                raise B.ContractError("downstream interpretation code family is not a literal tuple") from exc
    if declared != INTERPRETATION_CODE_PATHS:
        raise B.ContractError("downstream contract interpretation script family differs from revision-3")
    collator_functions = {
        node.name for node in parsed["scripts/77_collate_interpretation.py"].body
        if isinstance(node, ast.FunctionDef)
    }
    if not {"bh", "corrected_rows", "classify_cells", "classify_causal", "main"}.issubset(collator_functions):
        raise B.ContractError("interpretation collator lacks required BH/classification entrypoints")
    return {
        "python_sources_ast_validated": len(parsed),
        "production_entrypoints_with_main": len(INTERPRETATION_ENTRYPOINTS) - 1,
        "fail_fast_shell_entrypoints": 1,
        "downstream_declared_code_family_matches": True,
        "bh_and_classification_functions_present": True,
    }


def _fuma_header(
    root: Path,
    interpretation: Mapping[str, object],
    manifest: Mapping[str, object],
    dataset_id: str,
) -> list[str]:
    spec = interpretation.get("fuma_scrna")
    if not isinstance(spec, dict):
        raise B.ContractError("FUMA policy block is malformed")
    datasets = manifest.get("datasets")
    components = manifest.get("components")
    if not isinstance(datasets, list) or not isinstance(components, list):
        raise B.ContractError("FUMA dataset/component family is malformed")
    selected = [row for row in datasets if isinstance(row, dict) and row.get("dataset_id") == dataset_id]
    archives = [
        row for row in components if isinstance(row, dict)
        and row.get("component_id") == "FUMA_SCRNA_COMMIT_ARCHIVE"
    ]
    if len(selected) != 1 or len(archives) != 1:
        raise B.ContractError(f"FUMA header source is ambiguous: {dataset_id}")
    dataset, archive_row = selected[0], archives[0]
    archive_rel = str(archive_row["path"])
    archive_path = B._check_components(root, archive_rel)
    member_name = str(spec["archive_member_prefix"]) + str(dataset["member"])
    try:
        with tarfile.open(archive_path, "r:gz") as archive:
            member = archive.getmember(member_name)
            if not member.isfile() or Path(member.name).is_absolute() or ".." in Path(member.name).parts:
                raise B.ContractError(f"unsafe/non-file FUMA member: {member_name}")
            compressed = archive.extractfile(member)
            if compressed is None:
                raise B.ContractError(f"unreadable FUMA member: {member_name}")
            with compressed, gzip.GzipFile(fileobj=compressed) as stream:
                raw_header = stream.readline()
    except (OSError, tarfile.TarError, KeyError, gzip.BadGzipFile) as exc:
        raise B.ContractError(f"cannot stream FUMA header {dataset_id}: {exc}") from exc
    if hashlib.sha256(raw_header).hexdigest() != dataset.get("header_sha256"):
        raise B.ContractError(f"FUMA header hash differs from its pin: {dataset_id}")
    try:
        header = raw_header.decode("utf-8").split()
    except UnicodeError as exc:
        raise B.ContractError(f"FUMA header is not UTF-8: {dataset_id}") from exc
    if (
        len(header) != int(dataset["cell_type_count"]) + 2
        or header[0] != "GENE" or header[-1] != "Average"
        or len(header) != len(set(header))
    ):
        raise B.ContractError(f"FUMA header shape differs from its pin: {dataset_id}")
    return header[1:-1]


def validate_brain_scope(
    root: Path, manifests: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    allen = _fuma_header(
        root, manifests["interpretation"], manifests["interpretation_fuma_scrna"],
        "Allen_Human_MTG_level2",
    )
    replication = _fuma_header(
        root, manifests["interpretation"], manifests["interpretation_fuma_scrna"],
        "GSE67835_Human_Cortex_woFetal",
    )
    prefixes = ("Astro_", "Endo_", "Exc_", "Inh_", "Micro_", "Oligo_", "OPC_")
    counts = {prefix.rstrip("_"): sum(value.startswith(prefix) for value in allen) for prefix in prefixes}
    if len(allen) != 75 or any(count < 1 for count in counts.values()):
        raise B.ContractError("Allen MTG does not contain all seven named class prefixes across 75 subtypes")
    expected_replication = [
        "astrocytes", "endothelial", "hybrid", "microglia", "neurons",
        "oligodendrocytes", "OPC",
    ]
    if replication != expected_replication:
        raise B.ContractError("independent human replication is not the exact seven-class GSE67835 family")
    if any("pericy" in value.lower() for value in allen + replication):
        raise B.ContractError("pericyte absence assumption changed; revision-3 requires a new policy")
    return {
        "allen_subtype_columns": len(allen),
        "allen_named_class_prefix_counts": counts,
        "gse67835_class_columns": len(replication),
        "human_pericyte_columns": 0,
    }


def _exact_payload_ids(registry: B.EvidenceRegistry, dataset_ids: Sequence[str]) -> tuple[str, ...]:
    result: list[str] = []
    for dataset_id in dataset_ids:
        for suffix in (".tsv.gz", ".tsv.gz.tbi"):
            result.append(registry.id_for(f"{QTL_ROOT_REL}/{dataset_id}{suffix}"))
    return tuple(result)


def _always_blocked_fact(
    name: str,
    blocker_type: str,
    evidence_ids: Sequence[str],
    registry: B.EvidenceRegistry,
    absent_reason: str,
    present_reason: str,
) -> B.Fact:
    evidence = [registry.get(value) for value in evidence_ids]
    reason = absent_reason if all(item.kind == "MISSING" for item in evidence) else present_reason
    return B.Fact(name, False, blocker_type, reason, tuple(sorted(set(evidence_ids))))


def build_facts_v3(
    root: Path,
    effective: Mapping[str, object],
    overlay: Mapping[str, object],
    registry: B.EvidenceRegistry,
    manifests: Mapping[str, Mapping[str, object]],
    gtex_groups: Sequence[str],
) -> tuple[dict[str, B.Fact], dict[str, object]]:
    facts = B.build_facts(root, effective, registry, manifests)
    del facts["qtl_metadata"]
    del facts["qtl_payloads"]
    brain_scope = validate_brain_scope(root, manifests)
    code_scope = validate_interpretation_code_semantics(root)
    metadata_ids = tuple(registry.id_for(relative) for relative in CELL_METADATA_RELS)
    facts["cell_eqtl_metadata"] = B.Fact(
        "cell_eqtl_metadata", True, "DATA",
        "eQTL Catalogue r7 metadata and payload paths validate QTD000559 Young-2019 microglia ge and QTD000569 Aygun-2021 neuron ge",
        metadata_ids,
    )
    facts["cell_sqtl_metadata"] = B.Fact(
        "cell_sqtl_metadata", True, "DATA",
        "eQTL Catalogue r7 metadata and payload paths validate QTD000563 Young-2019 microglia leafcutter and QTD000573 Aygun-2021 neuron leafcutter",
        metadata_ids,
    )
    eqtl_ids = _exact_payload_ids(registry, CELL_EQTL_IDS)
    sqtl_ids = _exact_payload_ids(registry, CELL_SQTL_IDS)
    facts["cell_eqtl_payloads"] = _always_blocked_fact(
        "cell_eqtl_payloads", "DATA", eqtl_ids, registry,
        "QTD000559/QTD000569 eQTL payloads and matching .tbi indexes are absent",
        "one or more cell-eQTL files appeared without a complete post-acquisition checksum/semantic lock and cannot confer readiness",
    )
    facts["cell_sqtl_payloads"] = _always_blocked_fact(
        "cell_sqtl_payloads", "DATA", sqtl_ids, registry,
        "QTD000563/QTD000573 leafcutter sQTL payloads and matching .tbi indexes are absent",
        "one or more leafcutter cell-sQTL files appeared without a complete post-acquisition checksum/semantic lock and cannot confer readiness",
    )
    facts["gtex_eqtl_metadata"] = B.Fact(
        "gtex_eqtl_metadata", True, "DATA",
        "GTEx v8 imported ge metadata validates the complete ordered 49-context family",
        (registry.id_for(GTEX_METADATA_REL),),
    )
    code_ids = tuple(registry.id_for(relative) for relative in INTERPRETATION_CODE_PATHS)
    facts["interpretation_code_family"] = B.Fact(
        "interpretation_code_family", True, "SOFTWARE",
        "all 30 pre-result interpretation preparation, dispatch, runner, recorder and family-wide BH-collation entrypoints match the revision-3 checksum family",
        code_ids,
    )
    probes = overlay["new_absence_probes"]
    assert isinstance(probes, dict)
    gtex_root = str(probes["gtex_eqtl_payload_root"])
    gtex_ids = [registry.id_for(gtex_root)]
    for group in gtex_groups:
        gtex_ids.append(registry.id_for(f"{gtex_root}/{group}.tsv.gz"))
        gtex_ids.append(registry.id_for(f"{gtex_root}/{group}.tsv.gz.tbi"))
    facts["gtex_eqtl_payloads"] = _always_blocked_fact(
        "gtex_eqtl_payloads", "DATA", gtex_ids, registry,
        "all 49 GTEx v8 ge payloads and their 49 tabix indexes are absent",
        "one or more GTEx eQTL files appeared without a complete 49-context post-acquisition checksum/format/index lock and cannot confer readiness",
    )
    single_probes = {
        "gtex_sqtl_payloads": (
            "DATA", "no preregistered GTEx v8 tissue-sQTL manifest/payload/index family is local",
            "a GTEx tissue-sQTL manifest appeared without a checksum-pinned complete semantic family and cannot confer readiness",
            "gtex_sqtl_source_manifest",
        ),
        "sldsc_broad_tissue_coverage": (
            "DATA", "the exact local S-LDSC family contains 16 selected GTEx tissues plus control; "
            "the required chromosome-wide Lung index-36 and Muscle_Skeletal index-38 annotations and "
            "an 18-tissue checksum manifest are absent",
            "an S-LDSC 18-tissue manifest appeared without checksum-pinned chromosome-wide Lung-36 and "
            "Muscle_Skeletal-38 validation and cannot confer readiness",
            "sldsc_broad_tissue_coverage",
        ),
        "human_brain_pericyte_reference": (
            "DATA", "no checksum-pinned human brain expression atlas containing pericytes/vascular-associated cells is local",
            "a pericyte reference manifest appeared without a checksum-pinned human/build/context validator and cannot confer readiness",
            "human_brain_pericyte_reference",
        ),
        "brain_class_mapping_executor": (
            "SOFTWARE", "no frozen Allen-MTG 75-subtype-to-seven-class mapping/aggregation executor is local",
            "an Allen-MTG class executor appeared without a checksum pin and semantic mapping validator and cannot confer readiness",
            "brain_class_mapping_executor",
        ),
        "broad_tissue_magma_executor": (
            "SOFTWARE", "no checksum-pinned complete broad-tissue MAGMA adapter is local",
            "a broad-tissue MAGMA adapter appeared without a checksum pin and full-family semantic validator and cannot confer readiness",
            "broad_tissue_magma_executor",
        ),
        "cell_replication_executor": (
            "SOFTWARE", "no checksum-pinned cross-atlas class/subtype alignment and replication-classification executor is local",
            "a cell-replication executor appeared without a checksum pin and semantic alignment/classification validator and cannot confer readiness",
            "cell_replication_executor",
        ),
        "pair_comparison_executor": (
            "SOFTWARE", "no checksum-pinned deterministic Pair-A-versus-Pair-B comparison executor is local",
            "a pair-comparison executor appeared without a checksum pin and semantic output validator and cannot confer readiness",
            "pair_comparison_executor",
        ),
    }
    for name, (blocker, absent, present, probe_name) in single_probes.items():
        evidence_id = registry.id_for(str(probes[probe_name]))
        facts[name] = _always_blocked_fact(
            name, blocker, (evidence_id,), registry, absent, present,
        )
    if set(facts) != set(FACT_TYPES):
        raise B.ContractError(
            f"revision-3 fact family mismatch; missing={sorted(set(FACT_TYPES) - set(facts))} "
            f"extra={sorted(set(facts) - set(FACT_TYPES))}"
        )
    for name, fact in facts.items():
        if fact.blocker_type != FACT_TYPES[name]:
            raise B.ContractError(f"revision-3 fact blocker type drifted: {name}")
    return facts, {**brain_scope, **code_scope}


def evaluate_analyses(
    effective: Mapping[str, object],
    facts: Mapping[str, B.Fact],
    registry: B.EvidenceRegistry,
) -> list[dict[str, object]]:
    analyses = effective.get("analyses")
    if not isinstance(analyses, list):
        raise B.ContractError("revision-3 analysis rows are absent")
    result: list[dict[str, object]] = []
    for declared in analyses:
        assert isinstance(declared, dict)
        requirements = [str(value) for value in declared["requirements"]]
        failed = [facts[name] for name in requirements if not facts[name].ready]
        blocker_types = sorted({fact.blocker_type for fact in failed}, key=lambda value: (value != "DATA", value))
        resource_status = (
            "BLOCKED_BY_DATA" if "DATA" in blocker_types
            else "BLOCKED_BY_SOFTWARE" if blocker_types else "READY"
        )
        status = resource_status if failed else "NOT_APPLICABLE_UNTIL_UPSTREAM"
        evidence_ids = sorted({
            evidence_id for name in requirements for evidence_id in facts[name].evidence_ids
        })
        if not evidence_ids:
            raise B.ContractError(f"revision-3 analysis has no exact evidence: {declared['analysis_id']}")
        evidence = [registry.get(value) for value in evidence_ids]
        result.append({
            "analysis_id": declared["analysis_id"],
            "phase": declared["phase"],
            "analysis": declared["analysis"],
            "scope": declared["scope"],
            "resource_status": resource_status,
            "status": status,
            "blocker_types": ";".join(blocker_types) if blocker_types else "NONE",
            "status_reason": " | ".join(fact.reason for fact in failed) if failed else "local resource/software contract is READY; scientific execution awaits its predeclared upstream gate",
            "exact_unblock_condition": declared["unblock"],
            "required_inputs": declared["required_inputs"],
            "required_software": declared["required_software"],
            "reference": declared["reference"],
            "build": declared["build"],
            "ancestry": declared["ancestry"],
            "requirements": requirements,
            "evidence_ids": evidence_ids,
            "evidence_paths": ";".join(item.path for item in evidence),
            "evidence_sha256": ";".join(f"{item.path}={item.sha256}" for item in evidence),
            "ram_decomposition_safety": declared["ram_decomposition"],
            "planned_output": declared["output"],
            "claim_limit": declared["claim_limit"],
        })
    if [row["analysis_id"] for row in result] != list(EXPECTED_ANALYSIS_IDS):
        raise B.ContractError("revision-3 derived readiness lost or reordered an analysis")
    return result


def _paths(row: Mapping[str, object]) -> set[str]:
    return set(str(row["evidence_paths"]).split(";"))


def audit_analysis_family(
    readiness: Sequence[Mapping[str, object]],
    gtex_groups: Sequence[str],
) -> dict[str, object]:
    by_id = {str(row["analysis_id"]): row for row in readiness}
    if list(by_id) != list(EXPECTED_ANALYSIS_IDS) or len(by_id) != len(EXPECTED_ANALYSIS_IDS):
        raise B.ContractError("revision-3 audit found an incomplete/duplicated analysis family")
    for identity, expected in EXPECTED_REQUIREMENTS.items():
        if by_id[identity]["requirements"] != expected:
            raise B.ContractError(f"revision-3 audit found cross-family requirement drift: {identity}")
        unblock = str(by_id[identity]["exact_unblock_condition"]).casefold()
        if any(term.casefold() not in unblock for term in EXPECTED_UNBLOCK_TERMS[identity]):
            raise B.ContractError(f"revision-3 audit found an untied conditional unblock predicate: {identity}")
    eqtl_tokens = set(CELL_EQTL_IDS)
    sqtl_tokens = set(CELL_SQTL_IDS)
    for identity in ("P17_CELL_EQTL", "P19_THREE_WAY_EQTL", "P22_REGULATORY_CHAIN_EQTL"):
        paths = _paths(by_id[identity])
        if not all(any(token in path for path in paths) for token in eqtl_tokens) or any(
            token in path for path in paths for token in sqtl_tokens
        ):
            raise B.ContractError(f"revision-3 audit found eQTL/sQTL evidence contamination: {identity}")
    for identity in ("P18_CELL_SQTL", "P19_THREE_WAY_SQTL", "P22_REGULATORY_CHAIN_SQTL"):
        paths = _paths(by_id[identity])
        if not all(any(token in path for path in paths) for token in sqtl_tokens) or any(
            token in path for path in paths for token in eqtl_tokens
        ):
            raise B.ContractError(f"revision-3 audit found sQTL/eQTL evidence contamination: {identity}")
    gtex_paths = _paths(by_id["P13_GTEX_TISSUE_EQTL"])
    gtex_root = "ref/molecular/eqtl_catalogue_r7/imported/GTEx_V8/ge"
    expected_gtex = {
        f"{gtex_root}/{group}{suffix}" for group in gtex_groups
        for suffix in (".tsv.gz", ".tsv.gz.tbi")
    }
    if (
        not expected_gtex.issubset(gtex_paths) or GTEX_METADATA_REL not in gtex_paths
        or any("QTD000" in path for path in gtex_paths)
    ):
        raise B.ContractError("revision-3 audit found incomplete/contaminated GTEx eQTL evidence")
    gtex_sqtl_paths = _paths(by_id["P13_GTEX_TISSUE_SQTL"])
    if any("QTD000" in path or "/GTEx_V8/ge/" in path for path in gtex_sqtl_paths):
        raise B.ContractError("revision-3 audit found contaminated GTEx sQTL evidence")
    for identity in INTERPRETATION_CODE_ANALYSIS_IDS:
        if not set(INTERPRETATION_CODE_PATHS).issubset(_paths(by_id[identity])):
            raise B.ContractError(f"revision-3 audit found incomplete interpretation code evidence: {identity}")
    ready_ids = {
        "P15_BRAIN_CELL_SUBTYPE_DISCOVERY", "P20_CATLAS_SCATAC",
        "P20_SCREEN_CHROMATIN", "P21_ABC_ENHANCER_GENE",
        "P21_PCHIC_ENHANCER_GENE", "P26_POSITIVE_CONTROL", "P27_PATHWAYS",
    }
    software_blocked_ids = {
        "P07_CONJFDR", "P14_BRAIN_CELL_CLASS_DISCOVERY",
        "P16_HUMAN_CELL_CLASS_REPLICATION",
        "P24_LOCKED_ATLAS_CROSS_SLEEP", "P25_PAIR_A_VS_B",
        "P29_BIDIRECTIONAL_MR",
    }
    data_blocked_ids = set(EXPECTED_ANALYSIS_IDS) - ready_ids - software_blocked_ids
    expected_statuses = {
        **{identity: ("READY", "NOT_APPLICABLE_UNTIL_UPSTREAM") for identity in ready_ids},
        **{identity: ("BLOCKED_BY_SOFTWARE", "BLOCKED_BY_SOFTWARE") for identity in software_blocked_ids},
        **{identity: ("BLOCKED_BY_DATA", "BLOCKED_BY_DATA") for identity in data_blocked_ids},
    }
    for identity, expected in expected_statuses.items():
        observed = (by_id[identity]["resource_status"], by_id[identity]["status"])
        if observed != expected:
            raise B.ContractError(f"revision-3 audit found unexpected status for {identity}: {observed}")
    planned = [str(row["planned_output"]) for row in readiness]
    if len(planned) != len(set(planned)):
        raise B.ContractError("revision-3 audit found duplicated planned outputs")
    for row in readiness:
        if row["resource_status"] == "READY" and row["status"] != "NOT_APPLICABLE_UNTIL_UPSTREAM":
            raise B.ContractError(f"pre-result READY escaped upstream gating: {row['analysis_id']}")
        if row["resource_status"] != "READY" and "null" in str(row["status_reason"]).lower():
            raise B.ContractError(f"blocked unit was described as null: {row['analysis_id']}")
        for path in _paths(row):
            if any(path == prefix or path.startswith(prefix + "/") for prefix in B.FORBIDDEN_RESULT_PREFIXES):
                raise B.ContractError(f"future result leaked into revision-3 evidence: {row['analysis_id']}")
    return {
        "analysis_ids_audited": len(readiness),
        "exact_requirement_families_audited": len(EXPECTED_REQUIREMENTS),
        "conditional_unblock_predicates_audited": len(EXPECTED_UNBLOCK_TERMS),
        "qtl_families_separated": True,
        "gtex_eqtl_contexts": len(gtex_groups),
        "broad_class_mapping_fail_closed": True,
        "broad_tissue_adapter_fail_closed": True,
        "sldsc_18_tissue_coverage_fail_closed": True,
        "cell_replication_classifier_fail_closed": True,
        "interpretation_code_files_pinned": len(INTERPRETATION_CODE_PATHS),
        "pericyte_gap_exposed": True,
        "pair_comparison_executor_fail_closed": True,
    }


def render_report(
    policy_hash: str,
    evidence_rows: Sequence[Mapping[str, object]],
    readiness: Sequence[Mapping[str, object]],
    overlay: Mapping[str, object],
) -> bytes:
    resource_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    for row in readiness:
        resource_counts[str(row["resource_status"])] = resource_counts.get(str(row["resource_status"]), 0) + 1
        status_counts[str(row["status"])] = status_counts.get(str(row["status"]), 0) + 1
    supersedes = overlay["supersedes"]
    assert isinstance(supersedes, dict)
    lines = [
        "# Track B mechanism follow-up: pre-result local readiness — revision 3",
        "",
        f"Revision 3 supersedes `{supersedes['lock_path']}` (SHA-256 `{supersedes['lock_sha256']}`). Revisions 1 and 2 remain byte-exact audit trails and must not be used for execution. Revision 3 is the authoritative pre-result readiness contract.",
        "",
        "This contract used only local resources and immutable upstream metadata. It did not read fine-mapping, PLACO/pleiotropy, or future mechanism results. Missing, blocked, skipped and unavailable resources are never null findings.",
        "",
        f"- Policy SHA-256: `{policy_hash}`",
        f"- Complete analysis family: {len(readiness)} records",
        f"- Exact evidence family: {len(evidence_rows)} filesystem records",
        f"- Resource readiness: {', '.join(f'{key}={resource_counts[key]}' for key in sorted(resource_counts))}",
        f"- Current pre-result state: {', '.join(f'{key}={status_counts[key]}' for key in sorted(status_counts))}",
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
        "## Corrected scientific boundaries",
        "",
        "- GTEx tissue eQTL, GTEx tissue sQTL, cell eQTL and leafcutter cell sQTL are four distinct families. GTEx eQTL metadata has 49 frozen `ge` contexts, but the 49 payloads and indexes are absent. No GTEx tissue-sQTL family is locally preregistered.",
        "- Cell eQTL evidence is limited to QTD000559 microglia and QTD000569 neuron. Cell sQTL evidence is limited to QTD000563 microglia and QTD000573 neuron leafcutter data. Payloads/indexes and tabix are absent. PsychENCODE remains controlled supplementary access, not human replication or local evidence.",
        "- Allen Human MTG contains 75 subtype columns spanning seven named class prefixes, but no frozen subtype-to-broad-class mapping/aggregation executor exists. Broad-class discovery is therefore software-blocked; the full subtype unit remains resource-ready but awaits its upstream class gate. Neither Allen MTG nor GSE67835 contains pericytes, so the human pericyte unit is data-blocked.",
        "- GSE67835 is an independent human seven-class resource, but no frozen cross-atlas class alignment/replication classifier exists; two unadjusted enrichment tables are not replication. The two vascular FUMA matrices are mouse mapped to human gene identifiers and cannot be used as human replication. No airway or vascular human expression reference was invented.",
        "- The complete 30-file interpretation code family, including task preparation, runners, recorder and family-wide BH collator, is checksum-pinned. S-LDSC has exactly its locally frozen 16 selected GTEx tissues plus the shared control (375 files) and exact runtime; this is not the complete intended 18-tissue family. Chromosome-wide GTEx Lung index 36 and Muscle_Skeletal index 38 annotations are absent, and the static EUR baselineLD/weights/HapMap3/genotype reference is also absent. Lung must not be described as an airway-specific reference. Broad-tissue MAGMA lacks both its complete tissue matrix and a dedicated adapter.",
        "- CATlas adult human scATAC, SCREEN, ABC, immune-only PCHiC, HOCOMOCO motifs and four pathway resources validate locally. HOCOMOCO lacks a pinned hg38 sequence family. The Pair-A-versus-Pair-B executor is absent.",
        "- The 98-file/49-context TWAS inventory is frozen, but model payloads and the exact MetaXcan runtime are absent. Spatial data, a matched human multimodal GRN resource/runtime, HDL-L, rho-HESS and conjFDR remain blocked.",
        "- MR package/reference resources validate, but the automatic Track B bidirectional executor is absent. A causal claim remains prohibited unless every predeclared robustness gate passes.",
        "",
        "## RAM-equivalent decomposition",
        "",
        "Locus-dependent QTL, colocalization, chromatin linking and local-correlation analyses may run one complete locus per fresh process; an LD-dependent locus must never be split. Trait, tissue, cell and pathway work may run one frozen unit at a time while preserving the complete preregistered family and correction denominator. Global nuisance estimation must remain genome-wide; chunking is allowed only after global parameters are frozen and mathematical identity is demonstrated. GRN batching remains blocked until method-specific equivalence is proved.",
        "",
        "## Claim boundary",
        "",
        "`READY` means only that local resources/software passed this result-blind contract. Every such unit remains `NOT_APPLICABLE_UNTIL_UPSTREAM` until its scientific gate is satisfied. A blocker, unavailable resource, skipped unit or failure must be preserved and must never be reported as a negative biological result.",
        "",
    ])
    return "\n".join(lines).encode("utf-8")


def construct_artifacts(root: Path) -> tuple[dict[str, bytes], dict[str, object]]:
    overlay, effective, effective_hash, v2_lock = load_effective_policy(root)
    policy_hash = B.stable_file(root, POLICY_REL).sha256
    script_hash = B.stable_file(root, SCRIPT_REL).sha256
    registry, manifests = B.collect_evidence(root, effective)
    registry.add(POLICY_REL, expected_sha256=policy_hash, detail="revision-3 corrective policy")
    registry.add(SCRIPT_REL, expected_sha256=script_hash, detail="revision-3 contract implementation")
    registry.add(BASE_POLICY_REL, expected_sha256=BASE_POLICY_SHA256, detail="retained revision-2 policy")
    registry.add(BASE_REL, expected_sha256=BASE_SHA256, detail="retained revision-2 implementation")
    registry.add(BASE_LOCK_REL, expected_sha256=BASE_LOCK_SHA256, detail="retained byte-exact revision-2 lock")
    code_hashes = overlay["interpretation_code_sha256"]
    assert isinstance(code_hashes, dict)
    for relative in INTERPRETATION_CODE_PATHS:
        registry.add(
            relative, expected_sha256=str(code_hashes[relative]),
            detail="revision-3 checksum-pinned complete interpretation code family",
        )
    for relative, expected in sorted(v2_lock["outputs"].items()):
        registry.add(
            relative, expected_bytes=int(expected["bytes"]), expected_sha256=str(expected["sha256"]),
            detail="retained byte-exact revision-2 output",
        )
    registry.add(
        M2.OLD_LOCK_REL,
        expected_sha256=str(B.read_json(root, BASE_POLICY_REL)["supersedes"]["lock_sha256"]),
        detail="retained byte-exact revision-1 lock",
    )
    gtex_groups = validate_gtex_contexts(root)
    probes = overlay["new_absence_probes"]
    assert isinstance(probes, dict)
    for name in sorted(probes):
        registry.add(str(probes[name]), required=False, detail=f"revision-3 fail-closed readiness probe: {name}")
    gtex_root = str(probes["gtex_eqtl_payload_root"])
    for group in gtex_groups:
        registry.add(f"{gtex_root}/{group}.tsv.gz", required=False, detail=f"required GTEx v8 ge payload: {group}")
        registry.add(f"{gtex_root}/{group}.tsv.gz.tbi", required=False, detail=f"required GTEx v8 ge tabix index: {group}")
    facts, brain_scope = build_facts_v3(
        root, effective, overlay, registry, manifests, gtex_groups,
    )
    registry.recheck()
    evidence_rows = registry.rows()
    readiness = evaluate_analyses(effective, facts, registry)
    audit = audit_analysis_family(readiness, gtex_groups)
    evidence_digest = B.bytes_sha256(B.canonical_json(evidence_rows))
    resource_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    for row in readiness:
        resource_counts[str(row["resource_status"])] = resource_counts.get(str(row["resource_status"]), 0) + 1
        status_counts[str(row["status"])] = status_counts.get(str(row["status"]), 0) + 1
    document = {
        "schema_version": SCHEMA,
        "contract_revision": 3,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS",
        "authoritative_for_execution": True,
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": POLICY_REL,
        "policy_sha256": policy_hash,
        "effective_policy_sha256": effective_hash,
        "script_path": SCRIPT_REL,
        "script_sha256": script_hash,
        "supersedes": overlay["supersedes"],
        "retained_audit_trails": {
            "revision_1_lock": {"path": M2.OLD_LOCK_REL, "sha256": B.stable_file(root, M2.OLD_LOCK_REL).sha256},
            "revision_2_lock": {"path": BASE_LOCK_REL, "sha256": BASE_LOCK_SHA256},
        },
        "analysis_family_count": len(readiness),
        "evidence_record_count": len(evidence_rows),
        "evidence_bundle_sha256": evidence_digest,
        "resource_status_counts": dict(sorted(resource_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "scientific_scope_audit": {**audit, **brain_scope},
        "facts": [facts[name].record() for name in sorted(facts)],
        "analyses": readiness,
    }
    supersedes = overlay["supersedes"]
    assert isinstance(supersedes, dict)
    receipt = {
        "schema_version": SCHEMA,
        "record_kind": "IMMUTABLE_SUPERSESSION_RECEIPT",
        "superseded_contract": {
            "lock_path": BASE_LOCK_REL,
            "lock_sha256": BASE_LOCK_SHA256,
            "status": "SUPERSEDED_DO_NOT_USE_FOR_EXECUTION",
            "retention": supersedes["retention"],
        },
        "retained_predecessor": {
            "lock_path": M2.OLD_LOCK_REL,
            "lock_sha256": B.stable_file(root, M2.OLD_LOCK_REL).sha256,
            "status": "SUPERSEDED_DO_NOT_USE_FOR_EXECUTION",
        },
        "replacement_contract": {
            "policy_path": POLICY_REL,
            "policy_sha256": policy_hash,
            "effective_policy_sha256": effective_hash,
            "script_path": SCRIPT_REL,
            "script_sha256": script_hash,
            "lock_path": LOCK_REL,
        },
        "scientific_correction": supersedes["reason"],
        "result_blind": True,
        "future_results_accessed": False,
    }
    artifacts: dict[str, bytes] = {
        OUTPUT_RELS[0]: B.tsv_bytes(B.EVIDENCE_FIELDS, evidence_rows),
        OUTPUT_RELS[1]: B.tsv_bytes(B.READINESS_FIELDS, readiness),
        OUTPUT_RELS[2]: B.canonical_json(document),
        OUTPUT_RELS[3]: render_report(policy_hash, evidence_rows, readiness, overlay),
        SUPERSESSION_REL: B.canonical_json(receipt),
    }
    lock = {
        "schema_version": SCHEMA,
        "contract_revision": 3,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK",
        "authoritative_for_execution": True,
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": POLICY_REL,
        "policy_sha256": policy_hash,
        "effective_policy_sha256": effective_hash,
        "script_path": SCRIPT_REL,
        "script_sha256": script_hash,
        "superseded_lock_path": BASE_LOCK_REL,
        "superseded_lock_sha256": BASE_LOCK_SHA256,
        "retained_revision_1_lock_path": M2.OLD_LOCK_REL,
        "retained_revision_1_lock_sha256": B.stable_file(root, M2.OLD_LOCK_REL).sha256,
        "analysis_family_count": len(readiness),
        "evidence_record_count": len(evidence_rows),
        "evidence_bundle_sha256": evidence_digest,
        "outputs": {
            relative: {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
            for relative, payload in sorted(artifacts.items())
        },
    }
    artifacts[LOCK_REL] = B.canonical_json(lock)
    return artifacts, lock


def _same_payload(root: Path, relative: str, payload: bytes) -> bool:
    try:
        item = B.stable_file(root, relative)
    except FileNotFoundError:
        return False
    return item.observed_bytes == len(payload) and item.sha256 == hashlib.sha256(payload).hexdigest()


def freeze_no_replace(
    root: Path,
    artifacts: Mapping[str, bytes],
    *,
    after_link: Callable[[str], None] | None = None,
) -> None:
    output_dir = B._secure_mkdirs(root, OUTPUT_DIR_REL)
    ordered = list(OUTPUT_RELS) + [SUPERSESSION_REL, LOCK_REL]
    if set(artifacts) != set(ordered):
        raise B.ContractError("revision-3 artifact family differs from the exact canonical family")
    for relative in ordered:
        try:
            os.lstat(root / relative)
        except FileNotFoundError:
            continue
        if not _same_payload(root, relative, artifacts[relative]):
            raise B.ContractError(f"no-replace revision-3 freeze refused nonidentical target: {relative}")
    stage = Path(tempfile.mkdtemp(prefix=".readiness-v3-freeze-", dir=output_dir))
    staged: dict[str, Path] = {}
    try:
        for index, relative in enumerate(ordered):
            if _same_payload(root, relative, artifacts[relative]):
                continue
            stage_path = stage / f"{index:02d}.stage"
            fd = os.open(
                stage_path,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            try:
                view = memoryview(artifacts[relative])
                while view:
                    written = os.write(fd, view)
                    if written <= 0:
                        raise B.ContractError(f"short write while staging revision-3 artifact: {relative}")
                    view = view[written:]
                os.fsync(fd)
            finally:
                os.close(fd)
            staged[relative] = stage_path
        for relative in ordered:
            if not _same_payload(root, relative, artifacts[relative]):
                parent_rel = (root / relative).parent.relative_to(root).as_posix()
                B._secure_mkdirs(root, parent_rel)
                try:
                    os.link(staged[relative], root / relative, follow_symlinks=False)
                except FileExistsError:
                    if not _same_payload(root, relative, artifacts[relative]):
                        raise B.ContractError(f"revision-3 publication race produced nonidentical target: {relative}")
                if not _same_payload(root, relative, artifacts[relative]):
                    raise B.ContractError(f"published revision-3 artifact mismatch: {relative}")
                if after_link is not None:
                    after_link(relative)
        for parent in sorted({(root / relative).parent for relative in ordered}, key=str):
            fd = os.open(parent, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
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
    lock_payload = B.read_stable_bytes(root, LOCK_REL)
    try:
        lock = json.loads(lock_payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise B.ContractError(f"invalid revision-3 lock: {exc}") from exc
    if (
        not isinstance(lock, dict) or B.canonical_json(lock) != lock_payload
        or lock.get("schema_version") != SCHEMA or lock.get("contract_revision") != 3
        or lock.get("contract_kind") != "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK"
        or lock.get("authoritative_for_execution") is not True
        or lock.get("result_blind") is not True
        or lock.get("future_results_accessed") is not False
        or lock.get("superseded_lock_path") != BASE_LOCK_REL
        or lock.get("superseded_lock_sha256") != BASE_LOCK_SHA256
    ):
        raise B.ContractError("revision-3 lock schema/lineage/result-blind invariants failed")
    expected_paths = set(OUTPUT_RELS) | {SUPERSESSION_REL}
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != expected_paths:
        raise B.ContractError("revision-3 lock output family is incomplete")
    for relative in sorted(expected_paths):
        expected = outputs[relative]
        item = B.stable_file(root, relative)
        if (
            not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}
            or item.observed_bytes != expected["bytes"] or item.sha256 != expected["sha256"]
        ):
            raise B.ContractError(f"revision-3 output differs from lock: {relative}")
    recomputed, expected_lock = construct_artifacts(root)
    for relative in sorted(expected_paths):
        if B.read_stable_bytes(root, relative) != recomputed[relative]:
            raise B.ContractError(f"current resources no longer reproduce revision-3 output: {relative}")
    if lock != expected_lock or lock_payload != recomputed[LOCK_REL]:
        raise B.ContractError("current resources no longer reproduce revision-3 lock")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true", help="publish the revision-3 canonical family without replacement")
    mode.add_argument("--verify", action="store_true", help="verify frozen artifacts and current resources without writes")
    parser.add_argument("--root", type=Path, default=ROOT_AT_IMPORT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        raise B.ContractError(f"repository root is not a directory: {root}")
    if root != ROOT_AT_IMPORT and _bootstrap_hash(root / BASE_REL) != BASE_SHA256:
        raise B.ContractError("selected repository root lacks the exact revision-2 base")
    if args.freeze:
        artifacts, _lock = construct_artifacts(root)
        freeze_no_replace(root, artifacts)
        print(f"FROZEN_REVISION_3 {len(EXPECTED_ANALYSIS_IDS)} analyses at {OUTPUT_DIR_REL}")
    else:
        verify_frozen(root)
        print(f"VERIFIED_REVISION_3 {len(EXPECTED_ANALYSIS_IDS)} analyses at {OUTPUT_DIR_REL}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except B.ContractError as exc:
        raise SystemExit(f"ERROR: {exc}")
