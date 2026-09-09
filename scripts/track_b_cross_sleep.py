#!/usr/bin/env python3
"""Preflight and publish the result-blind Track B Phase 24 evidence view.

Phase 24 performs no new hypothesis test.  It exposes the locked 396-pair
global-rg result as context and accepts locus, variant, and gene-mechanism
evidence only through complete, pre-result-locked upstream families.  The
positive control is validated and published separately.  Missing upstream
science stops preflight without creating a directory or file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import stat
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

sys.dont_write_bytecode = True
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import track_b_phase24_25_common as C  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
POLICY_REL = Path("config/track_b_cross_sleep_policy.json")
CONTRACT_LOCK_REL = Path(
    "results/track_b/mechanism_followup/phase24_25/P24_CROSS_SLEEP.contract.lock.json"
)
SCRIPT_REL = Path("scripts/track_b_cross_sleep.py")
COMMON_REL = Path("scripts/track_b_phase24_25_common.py")

ANALYSIS_ID = "P24_LOCKED_ATLAS_CROSS_SLEEP"
POLICY_SCHEMA = "sleep-atlas-track-b-cross-sleep-policy.1"
LAYER_LOCK_SCHEMA = "sleep-atlas-track-b-cross-sleep-layer-lock.1"
FAMILY_LOCK_SCHEMA = "sleep-atlas-track-b-cross-sleep-pre-result-family.1"
CONTRACT_LOCK_SCHEMA = "sleep-atlas-track-b-cross-sleep-contract-lock.1"

PAIR_ORDER = ("A", "B", "CONTROL")
PAIR_IDENTITIES = {
    "A": ("snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    "B": ("insomnia", "adhd", "PRIMARY_DISCOVERY"),
    "CONTROL": ("insomnia", "frailty", "POSITIVE_CONTROL_NON_NOVELTY"),
}
SLEEP_TRAITS = (
    "insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype",
    "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency",
    "accel_sleep_duration", "sleep_timing",
)
LAYERS = (
    "GLOBAL_RG_CONTEXT", "LOCUS_SHARED_EVIDENCE", "VARIANT_SHARED_SIGNAL",
    "GENE_MECHANISM_EVIDENCE",
)
SUPPLEMENT_LAYERS = LAYERS[1:]
DIMENSIONS = {
    "BREATHING": ("snoring", "sleep_apnea"),
    "DISRUPTION": ("insomnia", "sleepiness", "napping", "sleep_efficiency"),
    "DURATION": ("sleepdur", "shortsleep", "longsleep", "accel_sleep_duration"),
    "TIMING": ("chronotype", "sleep_timing"),
}
UPSTREAM_STATUSES = {
    "COMPLETE", "NOT_APPLICABLE_UPSTREAM_ZERO", "BLOCKED_UPSTREAM",
    "FAILED_UPSTREAM",
}
OUTPUT_STATUSES = {
    "SUPPORTED", "TESTED_NO_SUPPORT", "NOT_APPLICABLE", "BLOCKED", "FAILED",
}
FORBIDDEN_MARKERS = ("SYNTHETIC", "PLACEHOLDER", "FAKE_RESULT", "SMOKE_TEST")
SAFE_ENTITY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:@/+|-]{0,511}$")
NA = "NA"

ATLAS_TRAIT_FIELDS = [
    "trait_id", "label", "domain", "type", "phenotype_definition", "source_id",
    "dataset_version", "pmid", "doi", "ancestry", "build", "ncase", "ncontrol",
    "n_total", "population_prevalence", "source_status", "harmonized_path",
    "harmonized_sha256", "harmonized_rows", "munged_path", "munged_sha256",
    "h2_scale", "h2", "h2_se", "h2_z", "ldsc_intercept", "ldsc_ratio",
    "phase1_verdict", "phase1_qc_reason",
]
ATLAS_PAIR_FIELDS = [
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
    "interpretation_status", "global_rg", "global_rg_se", "global_rg_z",
    "global_rg_p", "global_rg_fdr_all_396", "global_rg_fdr_primary_372",
    "global_rg_primary_significant", "effect_direction", "ldsc_input_log",
]
PAIR_MANIFEST_FIELDS = [
    "pair_id", "sleep_trait", "external_trait", "exact_GWAS", "source",
    "publication", "ancestry", "build", "sample_size", "case_control_counts",
    "discovery_rg", "SE", "P", "FDR", "cross_trait_intercept", "SNP_overlap",
    "h2", "replication_availability", "dense_data_readiness",
    "primary_scientific_role",
]
SUPPLEMENT_FIELDS = [
    "analysis_id", "pair_id", "family_role", "anchor_sleep_trait",
    "external_trait", "queried_sleep_trait", "evidence_level", "record_kind",
    "entity_id", "upstream_status", "family_complete", "tested_hypothesis_count",
    "failed_hypothesis_count", "supported_entity_count", "supported_entity_ids",
    "estimate", "standard_error", "p_value", "adjusted_p_value",
    "effect_direction", "power_qc_status", "upstream_correction_method",
    "upstream_correction_denominator", "pre_result_family_lock_path",
    "pre_result_family_lock_sha256", "source_result_path", "source_result_sha256",
    "claim_limit",
]
OUTPUT_FIELDS = [
    "analysis_id", "evidence_row_id", "pair_id", "family_role",
    "anchor_sleep_trait", "external_trait", "queried_sleep_trait",
    "evidence_level", "record_kind", "entity_id", "evidence_status",
    "classification", "supported_entity_count", "supported_entity_ids",
    "estimate", "standard_error", "p_value", "adjusted_p_value",
    "effect_direction", "power_qc_status", "tested_hypothesis_count",
    "failed_hypothesis_count", "upstream_correction_method",
    "upstream_correction_denominator", "phase24_correction_denominator",
    "source_result_path", "source_result_sha256", "claim_limit",
]
GENERALIZATION_FIELDS = [
    "analysis_id", "pair_id", "family_role", "anchor_sleep_trait",
    "external_trait", "evidence_level", "entity_id",
    "supported_sleep_trait_count", "supported_sleep_traits",
    "supported_dimension_count", "supported_dimensions",
    "generalization_classification", "coverage_status", "source_evidence_row_ids",
    "claim_limit",
]

LAYER_PREFIX = {
    "LOCUS_SHARED_EVIDENCE": "LOCUS",
    "VARIANT_SHARED_SIGNAL": "VARIANT",
    "GENE_MECHANISM_EVIDENCE": "GENE_MECHANISM",
}


def _integer(value: str, label: str, *, minimum: int = 0) -> int:
    if not re.fullmatch(r"0|[1-9][0-9]*", value):
        raise C.ContractError(f"invalid integer {label}: {value!r}")
    result = int(value)
    if result < minimum:
        raise C.ContractError(f"{label} is below {minimum}")
    return result


def _probability(value: str, label: str, *, allow_na: bool = False) -> float | None:
    if value == NA and allow_na:
        return None
    try:
        result = float(value)
    except ValueError as error:
        raise C.ContractError(f"invalid probability {label}: {value!r}") from error
    if not math.isfinite(result) or not 0 <= result <= 1:
        raise C.ContractError(f"invalid probability {label}: {value!r}")
    return result


def _finite_or_na(value: str, label: str) -> float | None:
    if value == NA:
        return None
    try:
        result = float(value)
    except ValueError as error:
        raise C.ContractError(f"invalid numeric value {label}: {value!r}") from error
    if not math.isfinite(result):
        raise C.ContractError(f"non-finite numeric value {label}")
    return result


def _yes_no(value: str, label: str) -> bool:
    if value not in {"TRUE", "FALSE"}:
        raise C.ContractError(f"{label} must be TRUE or FALSE")
    return value == "TRUE"


def query_keys(layer: str | None = None) -> list[str]:
    layers: Sequence[str] = (layer,) if layer is not None else LAYERS
    return [
        f"{pair}\x1f{sleep}\x1f{evidence_layer}"
        for pair in PAIR_ORDER
        for sleep in SLEEP_TRAITS
        for evidence_layer in layers
    ]


def evidence_row_id(row: Mapping[str, str]) -> str:
    values = (
        row["pair_id"], row["queried_sleep_trait"], row["evidence_level"],
        row["record_kind"], row["entity_id"], row["source_result_sha256"],
    )
    return "P24ROW__" + hashlib.sha256("\x1f".join(values).encode()).hexdigest()[:24]


def validate_policy(root: Path) -> tuple[dict[str, Any], dict[str, object]]:
    policy, identity = C.stable_json(root, POLICY_REL)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("analysis_id") != ANALYSIS_ID
        or policy.get("selection_timing")
        != "BEFORE_ANY_PHASE24_CROSS_SLEEP_RESULT_ACCESS"
        or policy.get("result_blind") is not True
        or policy.get("future_results_accessed_before_freeze") is not False
    ):
        raise C.ContractError("Phase24 policy identity or result-blind timing drifted")
    if tuple(policy.get("sleep_traits_in_locked_panel_order", [])) != SLEEP_TRAITS:
        raise C.ContractError("Phase24 sleep-trait family/order drifted")
    if tuple(policy.get("evidence_layers_in_order", [])) != LAYERS:
        raise C.ContractError("Phase24 evidence-layer order drifted")
    scope = policy.get("pair_scope", {})
    observed_pairs = {
        row.get("pair_id"): (
            row.get("anchor_sleep_trait"), row.get("external_trait"),
            row.get("family_role"),
        )
        for row in scope.get("pairs", [])
    }
    if (
        tuple(scope.get("pair_order", [])) != PAIR_ORDER
        or scope.get("primary_pairs") != ["A", "B"]
        or scope.get("control_pairs") != ["CONTROL"]
        or observed_pairs != PAIR_IDENTITIES
    ):
        raise C.ContractError("Phase24 pair identities or CONTROL separation drifted")
    queries = policy.get("query_family", {})
    if (
        queries.get("primary_query_count") != 24
        or queries.get("control_query_count") != 12
        or queries.get("total_query_count") != 36
        or not all(queries.get(key) is True for key in (
            "unrestricted_phenome_scan_forbidden",
            "non_sleep_query_traits_forbidden_except_frozen_external_targets",
            "result_driven_sleep_trait_selection_forbidden",
            "result_ranked_locus_or_entity_selection_forbidden",
        ))
    ):
        raise C.ContractError("Phase24 fixed query scope or anti-selection rule drifted")
    dimensions = policy.get("sleep_dimensions", {})
    if {key: tuple(dimensions.get(key, [])) for key in DIMENSIONS} != DIMENSIONS:
        raise C.ContractError("Phase24 sleep-dimension partition drifted")
    flattened = [trait for name in DIMENSIONS for trait in DIMENSIONS[name]]
    if len(flattened) != len(set(flattened)) or set(flattened) != set(SLEEP_TRAITS):
        raise C.ContractError("sleep dimensions do not partition the twelve traits")
    mechanism = policy.get("mechanism_level_generalization", {})
    expected_classes = [
        "SNORING_SPECIFIC", "INSOMNIA_SPECIFIC", "SLEEP_BREATHING_SHARED",
        "SLEEP_DISRUPTION_SHARED", "MULTIPLE_SLEEP_DIMENSIONS", "NONSPECIFIC",
        "NO_DETECTED_SUPPORT", "INCONCLUSIVE", "BLOCKED",
    ]
    if (
        mechanism.get("applicable_layers") != list(SUPPLEMENT_LAYERS)
        or mechanism.get("allowed_classifications") != expected_classes
        or "never" not in mechanism.get("global_rg_exclusion", "")
    ):
        raise C.ContractError("Phase24 mechanism-level classification contract drifted")
    if policy.get("supplement_schema") != SUPPLEMENT_FIELDS:
        raise C.ContractError("Phase24 supplement schema drifted")
    output = policy.get("output_contract", {})
    if output.get("fields") != OUTPUT_FIELDS or output.get("generalization_fields") != GENERALIZATION_FIELDS:
        raise C.ContractError("Phase24 output schemas drifted")
    if (
        policy.get("phase24_multiple_testing", {}).get("new_hypothesis_count_primary") != 0
        or policy.get("phase24_multiple_testing", {}).get("new_hypothesis_count_control") != 0
        or policy.get("layer_contract", {}).get("GLOBAL_RG_CONTEXT", {}).get(
            "source_correction_denominator"
        ) != 396
    ):
        raise C.ContractError("Phase24 correction denominators drifted")
    if set(policy.get("upstream_statuses", [])) != UPSTREAM_STATUSES:
        raise C.ContractError("Phase24 upstream states drifted")
    if set(policy.get("output_statuses", [])) != OUTPUT_STATUSES:
        raise C.ContractError("Phase24 output states drifted")
    if policy.get("upstream_gate", {}).get("global_rg_not_sufficient") is not True:
        raise C.ContractError("global rg must not satisfy the Phase24 upstream gate")
    return policy, identity


def _validate_pinned_artifact(
    root: Path, definition: Mapping[str, Any], label: str,
) -> dict[str, object]:
    path = definition.get("path")
    digest = definition.get("sha256")
    if not isinstance(path, str):
        raise C.ContractError(f"missing immutable path for {label}")
    C.require_sha256(digest, label)
    _payload, identity = C.stable_bytes(root, path)
    if identity["sha256"] != digest:
        raise C.ContractError(f"immutable Phase24 artifact drifted: {path}")
    return identity


def validate_core(
    root: Path, policy: Mapping[str, Any], *, deep_hash_dense: bool = False,
) -> dict[str, Any]:
    immutable = policy["immutable_core"]
    identities: dict[str, dict[str, object]] = {}
    for label, definition in immutable.items():
        if isinstance(definition, dict) and "path" in definition:
            identities[label] = _validate_pinned_artifact(root, definition, label)

    _, panel, _ = C.stable_tsv(root, immutable["panel"]["path"])
    panel_ids = [row.get("trait_id", "") for row in panel]
    sleep_ids = [row.get("trait_id", "") for row in panel if row.get("domain") == "sleep"]
    non_sleep_ids = [row.get("trait_id", "") for row in panel if row.get("domain") != "sleep"]
    if len(panel_ids) != 45 or len(set(panel_ids)) != 45 or sleep_ids != list(SLEEP_TRAITS) or len(non_sleep_ids) != 33:
        raise C.ContractError("analysis panel is not the locked 45/12/33 family")
    panel_lock, _ = C.stable_json(root, immutable["panel_lock"]["path"])
    ordered = hashlib.sha256("".join(f"{trait}\n" for trait in panel_ids).encode()).hexdigest()
    if (
        panel_lock.get("trait_count") != 45
        or panel_lock.get("sleep_trait_count") != 12
        or panel_lock.get("non_sleep_trait_count") != 33
        or panel_lock.get("ordered_trait_ids_sha256") != ordered
    ):
        raise C.ContractError("analysis-panel lock no longer binds the ordered family")

    _, traits, _ = C.stable_tsv(
        root, immutable["atlas_traits"]["path"], ATLAS_TRAIT_FIELDS,
    )
    if len(traits) != 45 or [row["trait_id"] for row in traits] != panel_ids:
        raise C.ContractError("atlas traits table is not the ordered 45-trait family")
    trait_by_id = {row["trait_id"]: row for row in traits}
    _, atlas_pairs, pair_identity = C.stable_tsv(
        root, immutable["atlas_pairs"]["path"], ATLAS_PAIR_FIELDS,
    )
    expected_pairs = [(sleep, other) for sleep in SLEEP_TRAITS for other in non_sleep_ids]
    observed_pairs = [(row["sleep_trait"], row["non_sleep_trait"]) for row in atlas_pairs]
    if len(atlas_pairs) != 396 or observed_pairs != expected_pairs:
        raise C.ContractError("atlas pair table is not the exact ordered 396-pair family")
    if any(row["pair_id"] != f"{row['sleep_trait']}__{row['non_sleep_trait']}" for row in atlas_pairs):
        raise C.ContractError("atlas pair IDs are noncanonical")
    core, _ = C.stable_json(root, immutable["atlas_core_provenance"]["path"])
    if (
        core.get("trait_count") != 45
        or core.get("pair_count") != 396
        or core.get("panel_sha256") != immutable["panel"]["sha256"]
        or core.get("panel_lock_sha256") != immutable["panel_lock"]["sha256"]
        or core.get("outputs", {}).get("results/atlas/traits.tsv")
        != immutable["atlas_traits"]["sha256"]
        or core.get("outputs", {}).get("results/atlas/trait_pairs.tsv")
        != pair_identity["sha256"]
        or core.get("status") != "CORE_ATLAS_COMPLETE_DOWNSTREAM_LAYERS_PENDING"
    ):
        raise C.ContractError("atlas core provenance does not bind the locked core")

    _, pair_manifest, _ = C.stable_tsv(
        root, immutable["pair_manifest"]["path"], PAIR_MANIFEST_FIELDS,
    )
    if len(pair_manifest) != 3 or [row["pair_id"] for row in pair_manifest] != list(PAIR_ORDER):
        raise C.ContractError("Track B pair manifest is not ordered A/B/CONTROL")
    for row in pair_manifest:
        anchor, external, _role = PAIR_IDENTITIES[row["pair_id"]]
        if (
            row["sleep_trait"] != anchor
            or row["external_trait"] != external
            or row["dense_data_readiness"] != "READY_FULL_SUMSTATS_BOTH"
        ):
            raise C.ContractError(f"Track B pair identity/readiness drifted for {row['pair_id']}")
    pair_lock, _ = C.stable_json(root, immutable["pair_manifest_lock"]["path"])
    if (
        pair_lock.get("output_sha256", {}).get("results/track_b/pair_manifest.tsv")
        != immutable["pair_manifest"]["sha256"]
        or pair_lock.get("pair_replacement_policy", "").split(";", 1)[0] != "FORBIDDEN"
    ):
        raise C.ContractError("Track B pair-manifest lock drifted")

    dense_identities: list[dict[str, object]] = []
    dense_definitions = policy.get("dense_sleep_family", [])
    if [row.get("trait_id") for row in dense_definitions] != list(SLEEP_TRAITS):
        raise C.ContractError("dense sleep family is not in locked order")
    for definition in dense_definitions:
        trait = definition["trait_id"]
        atlas = trait_by_id[trait]
        if (
            definition.get("path") != atlas["harmonized_path"]
            or definition.get("sha256") != atlas["harmonized_sha256"]
            or str(definition.get("rows")) != atlas["harmonized_rows"]
            or definition.get("build") != atlas["build"]
        ):
            raise C.ContractError(f"dense input definition differs from atlas for {trait}")
        path = C.safe_path(root, definition["path"], "dense sleep input", must_exist=True)
        info = os.lstat(path)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_size <= 0:
            raise C.ContractError(f"dense sleep input is not a real non-empty file: {trait}")
        observed: dict[str, object] = {
            "path": definition["path"], "bytes": info.st_size,
            "expected_sha256": definition["sha256"], "hash_mode": "LOCK_CHAIN_ONLY",
        }
        if deep_hash_dense:
            checked = C.stable_identity(root, definition["path"])
            if checked["sha256"] != definition["sha256"]:
                raise C.ContractError(f"dense sleep summary differs from its lock: {trait}")
            observed.update(checked)
            observed["hash_mode"] = "FULL_REHASH"
        dense_identities.append(observed)

    atlas_schema, _ = C.stable_json(root, immutable["atlas_schema"]["path"])
    tables = atlas_schema.get("tables", {})
    if (
        atlas_schema.get("schema_version") != "atlas-v1.0.0"
        or tables.get("traits.tsv", {}).get("fields") != ATLAS_TRAIT_FIELDS
        or tables.get("trait_pairs.tsv", {}).get("fields") != ATLAS_PAIR_FIELDS
        or list(tables) != [
            "traits.tsv", "trait_pairs.tsv", "loci.tsv", "variants.tsv", "genes.tsv",
            "regulatory_elements.tsv", "cell_types.tsv", "pathways.tsv",
            "causal_tests.tsv", "edges.tsv",
        ]
    ):
        raise C.ContractError("downstream atlas schema no longer has the frozen ten-table family")
    return {
        "artifact_identities": identities,
        "dense_sleep": dense_identities,
        "atlas_pair_rows": atlas_pairs,
        "atlas_pair_identity": pair_identity,
        "pair_manifest_rows": pair_manifest,
    }


def validate_contract_lock(
    root: Path, policy_identity: Mapping[str, object],
) -> tuple[dict[str, Any], dict[str, object]]:
    lock, identity = C.stable_json(root, CONTRACT_LOCK_REL)
    if (
        lock.get("schema_version") != CONTRACT_LOCK_SCHEMA
        or lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("result_blind") is not True
        or lock.get("future_results_accessed") is not False
        or lock.get("policy", {}).get("path") != str(POLICY_REL)
        or lock.get("policy", {}).get("sha256") != policy_identity["sha256"]
    ):
        raise C.ContractError("Phase24 contract lock identity drifted")
    for relative, key in ((SCRIPT_REL, "executor"), (COMMON_REL, "common_io")):
        observed = C.stable_identity(root, relative)
        expected = lock.get("implementation", {}).get(key, {})
        if expected.get("path") != str(relative) or expected.get("sha256") != observed["sha256"]:
            raise C.ContractError(f"Phase24 locked implementation drifted: {relative}")
    if lock.get("output_contract_sha256") is None:
        raise C.ContractError("Phase24 lock lacks an output-contract digest")
    C.require_sha256(lock["output_contract_sha256"], "Phase24 output contract")
    return lock, identity


def _entity_list(value: str, *, expected_count: int, label: str) -> list[str]:
    if expected_count == 0:
        if value != NA:
            raise C.ContractError(f"{label} must be NA when no entity is supported")
        return []
    if value == NA:
        raise C.ContractError(f"{label} is missing despite supported entities")
    values = value.split(";")
    if (
        len(values) != expected_count
        or values != sorted(values)
        or len(values) != len(set(values))
        or any(SAFE_ENTITY.fullmatch(item) is None for item in values)
    ):
        raise C.ContractError(f"{label} is not an exact sorted unique entity family")
    return values


def validate_family_lock(
    root: Path,
    relative: str,
    expected_sha256: str,
    layer: str,
) -> tuple[dict[str, Any], dict[str, object]]:
    lock, identity = C.stable_json(root, relative)
    if identity["sha256"] != C.require_sha256(expected_sha256, "pre-result family lock"):
        raise C.ContractError(f"pre-result family lock hash drifted for {layer}")
    expected_queries = query_keys(layer)
    if (
        lock.get("schema_version") != FAMILY_LOCK_SCHEMA
        or lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("evidence_level") != layer
        or lock.get("query_keys_in_order") != expected_queries
        or lock.get("query_count") != 36
        or lock.get("results_accessed_before_family_lock") is not False
        or lock.get("complete_hypothesis_family") is not True
        or lock.get("result_driven_trait_or_locus_selection") is not False
        or lock.get("all_hypotheses_retained_for_correction") is not True
    ):
        raise C.ContractError(f"invalid pre-result family lock for {layer}")
    primary_denominator = lock.get("primary_correction_denominator")
    control_denominator = lock.get("control_correction_denominator")
    if (
        type(primary_denominator) is not int or primary_denominator <= 0
        or type(control_denominator) is not int or control_denominator <= 0
    ):
        raise C.ContractError(f"pre-result correction denominators are not exact for {layer}")
    if layer == "LOCUS_SHARED_EVIDENCE" and (
        primary_denominator != 59_880 or control_denominator != 29_940
    ):
        raise C.ContractError("locus family must account for 2495 blocks x fixed queries")
    C.require_sha256(lock.get("family_definition_sha256"), f"{layer} family definition")
    return lock, identity


def _row_base_validation(row: Mapping[str, str], layer: str) -> None:
    if row["analysis_id"] != ANALYSIS_ID or row["evidence_level"] != layer:
        raise C.ContractError(f"supplement row has wrong analysis/layer for {layer}")
    pair = row["pair_id"]
    if pair not in PAIR_IDENTITIES:
        raise C.ContractError(f"out-of-scope pair in {layer}: {pair}")
    anchor, external, role = PAIR_IDENTITIES[pair]
    if (
        row["anchor_sleep_trait"] != anchor
        or row["external_trait"] != external
        or row["family_role"] != role
        or row["queried_sleep_trait"] not in SLEEP_TRAITS
    ):
        raise C.ContractError(f"trait or role substitution in {layer}/{pair}")
    if any(marker in "\t".join(row.values()).upper() for marker in FORBIDDEN_MARKERS):
        raise C.ContractError(f"forbidden synthetic/placeholder marker in {layer}")
    if row["upstream_status"] not in UPSTREAM_STATUSES:
        raise C.ContractError(f"invalid upstream state in {layer}")
    if row["record_kind"] not in {"COVERAGE", "EVIDENCE"}:
        raise C.ContractError(f"invalid record kind in {layer}")
    C.require_sha256(row["pre_result_family_lock_sha256"], "pre-result family lock")
    C.require_sha256(row["source_result_sha256"], "source result")
    if not row["claim_limit"] or row["claim_limit"] == NA:
        raise C.ContractError(f"missing claim limit in {layer}")


def validate_supplement_rows(
    rows: Sequence[dict[str, str]],
    layer: str,
    *,
    primary_denominator: int,
    control_denominator: int,
) -> list[dict[str, str]]:
    """Validate one complete 36-query layer; exposed for adversarial tests."""

    if layer not in SUPPLEMENT_LAYERS:
        raise C.ContractError(f"invalid supplement layer: {layer}")
    expected_coverage = [
        (pair, sleep, layer) for pair in PAIR_ORDER for sleep in SLEEP_TRAITS
    ]
    coverage: dict[tuple[str, str, str], dict[str, str]] = {}
    evidence: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    seen_evidence: set[tuple[str, str, str, str]] = set()
    for row in rows:
        if list(row) != SUPPLEMENT_FIELDS:
            raise C.ContractError("supplement row field order differs from the contract")
        _row_base_validation(row, layer)
        key = (row["pair_id"], row["queried_sleep_trait"], layer)
        if row["record_kind"] == "COVERAGE":
            if row["entity_id"] != NA or key in coverage:
                raise C.ContractError(f"duplicate or entity-bearing coverage row: {key}")
            coverage[key] = row
        else:
            identity = (*key, row["entity_id"])
            if row["entity_id"] == NA or SAFE_ENTITY.fullmatch(row["entity_id"]) is None:
                raise C.ContractError(f"invalid supported entity ID in {layer}")
            if identity in seen_evidence:
                raise C.ContractError(f"duplicate supported entity in {layer}: {identity}")
            seen_evidence.add(identity)
            evidence[key].append(row)
    if list(coverage) != expected_coverage:
        raise C.ContractError(f"{layer} lacks exact A/B/CONTROL x 12 ordered coverage")

    for key in expected_coverage:
        row = coverage[key]
        pair = key[0]
        status = row["upstream_status"]
        complete = _yes_no(row["family_complete"], "family_complete")
        supported_count = _integer(row["supported_entity_count"], "supported_entity_count")
        supported = _entity_list(
            row["supported_entity_ids"], expected_count=supported_count,
            label="supported_entity_ids",
        )
        observed_evidence = sorted(item["entity_id"] for item in evidence.get(key, []))
        if observed_evidence != supported:
            raise C.ContractError(f"coverage/evidence entity mismatch for {key}")
        expected_denominator = primary_denominator if pair != "CONTROL" else control_denominator
        if status == "COMPLETE":
            tested = _integer(row["tested_hypothesis_count"], "tested_hypothesis_count", minimum=1)
            failed = _integer(row["failed_hypothesis_count"], "failed_hypothesis_count")
            denominator = _integer(
                row["upstream_correction_denominator"],
                "upstream_correction_denominator", minimum=1,
            )
            if not complete or failed > tested or denominator != expected_denominator:
                raise C.ContractError(f"invalid complete-family accounting for {key}")
            if row["power_qc_status"] in {"", NA}:
                raise C.ContractError(f"complete row lacks power/QC status for {key}")
        elif status == "NOT_APPLICABLE_UPSTREAM_ZERO":
            if (
                not complete or supported_count != 0
                or row["tested_hypothesis_count"] != "0"
                or row["failed_hypothesis_count"] != "0"
                or row["upstream_correction_denominator"] != str(expected_denominator)
                or row["power_qc_status"] != "NOT_APPLICABLE_UPSTREAM_ZERO"
            ):
                raise C.ContractError(f"invalid upstream-zero accounting for {key}")
        elif status == "BLOCKED_UPSTREAM":
            if (
                complete or supported_count != 0
                or any(row[field] != NA for field in (
                    "tested_hypothesis_count", "failed_hypothesis_count",
                    "upstream_correction_denominator",
                ))
                or not row["power_qc_status"].startswith("BLOCKED_")
            ):
                raise C.ContractError(f"blocked row is disguised as tested for {key}")
        else:
            if (
                complete or supported_count != 0
                or row["tested_hypothesis_count"] == NA
                or row["failed_hypothesis_count"] == NA
                or row["upstream_correction_denominator"] != str(expected_denominator)
                or not row["power_qc_status"].startswith("FAILED_")
            ):
                raise C.ContractError(f"failed row lacks explicit accounting for {key}")
            tested = _integer(row["tested_hypothesis_count"], "tested_hypothesis_count")
            failed = _integer(row["failed_hypothesis_count"], "failed_hypothesis_count")
            if failed <= 0 or failed > tested:
                raise C.ContractError(f"failed row has invalid failure count for {key}")
        if row["record_kind"] != "COVERAGE":
            raise AssertionError("coverage map contains non-coverage row")
        for statistic in ("estimate", "standard_error", "p_value", "adjusted_p_value", "effect_direction"):
            if row[statistic] != NA:
                raise C.ContractError(f"coverage summary carries selected statistic for {key}")

        for evidence_row in evidence.get(key, []):
            if status != "COMPLETE" or supported_count == 0:
                raise C.ContractError(f"evidence row exists for non-complete family {key}")
            if (
                evidence_row["upstream_status"] != "COMPLETE"
                or evidence_row["family_complete"] != "TRUE"
                or evidence_row["supported_entity_count"] != "1"
                or evidence_row["supported_entity_ids"] != evidence_row["entity_id"]
                or evidence_row["tested_hypothesis_count"] != row["tested_hypothesis_count"]
                or evidence_row["failed_hypothesis_count"] != row["failed_hypothesis_count"]
                or evidence_row["upstream_correction_method"] != row["upstream_correction_method"]
                or evidence_row["upstream_correction_denominator"]
                != row["upstream_correction_denominator"]
                or evidence_row["pre_result_family_lock_path"]
                != row["pre_result_family_lock_path"]
                or evidence_row["pre_result_family_lock_sha256"]
                != row["pre_result_family_lock_sha256"]
                or evidence_row["source_result_path"] != row["source_result_path"]
                or evidence_row["source_result_sha256"] != row["source_result_sha256"]
            ):
                raise C.ContractError(f"evidence row escapes coverage identity for {key}")
            p_value = _probability(evidence_row["p_value"], "evidence p", allow_na=True)
            adjusted = _probability(
                evidence_row["adjusted_p_value"], "evidence adjusted p", allow_na=True,
            )
            _finite_or_na(evidence_row["estimate"], "evidence estimate")
            se = _finite_or_na(evidence_row["standard_error"], "evidence standard error")
            if se is not None and se < 0:
                raise C.ContractError(f"negative standard error for {key}")
            if p_value is None or adjusted is None:
                if "NO_SINGLE_P_VALUE_BY_DESIGN" not in evidence_row["power_qc_status"]:
                    raise C.ContractError(f"missing evidence statistics without method reason for {key}")
            if evidence_row["power_qc_status"] in {"", NA}:
                raise C.ContractError(f"evidence row lacks exact power/QC status for {key}")
    return [dict(row) for row in rows]


def validate_supplement(
    root: Path,
    policy: Mapping[str, Any],
    policy_identity: Mapping[str, object],
    layer: str,
) -> tuple[list[dict[str, str]], dict[str, object], dict[str, object]]:
    definition = policy["layer_contract"][layer]
    table_relative = definition["future_input"]
    lock_relative = definition["future_input_lock"]
    fields, rows, table_identity = C.stable_tsv(root, table_relative, SUPPLEMENT_FIELDS)
    del fields
    lock, lock_identity = C.stable_json(root, lock_relative)
    required = policy["supplement_lock_requirements"]
    if (
        lock.get("schema_version") != LAYER_LOCK_SCHEMA
        or lock.get("analysis_id") != ANALYSIS_ID
        or lock.get("evidence_level") != layer
        or lock.get("policy_path") != str(POLICY_REL)
        or lock.get("policy_sha256") != policy_identity["sha256"]
        or lock.get("coverage_table") != table_relative
        or lock.get("coverage_table_sha256") != table_identity["sha256"]
        or lock.get("coverage_row_count") != 36
        or lock.get("evidence_row_count") != sum(
            row["record_kind"] == "EVIDENCE" for row in rows
        )
        or lock.get("query_keys_sha256") != C.digest_json(query_keys(layer))
        or lock.get("results_accessed_before_family_lock")
        is not required["results_accessed_before_family_lock"]
        or lock.get("complete_fixed_query_family")
        is not required["complete_fixed_query_family"]
        or lock.get("result_ranked_selection") is not required["result_ranked_selection"]
        or lock.get("all_supported_entities_retained")
        is not required["all_supported_entities_retained"]
        or lock.get("exact_statistics_and_power_qc_retained")
        is not required["exact_statistics_and_power_qc_retained"]
        or lock.get("raw_or_uncorrected_p_value_substitution_forbidden")
        is not required["raw_or_uncorrected_p_value_substitution_forbidden"]
    ):
        raise C.ContractError(f"invalid Phase24 layer result lock for {layer}")
    family_relative = lock.get("pre_result_family_lock_path")
    family_sha = lock.get("pre_result_family_lock_sha256")
    if not isinstance(family_relative, str) or not isinstance(family_sha, str):
        raise C.ContractError(f"missing pre-result family reference for {layer}")
    family, _family_identity = validate_family_lock(
        root, family_relative, family_sha, layer,
    )
    if any(
        row["pre_result_family_lock_path"] != family_relative
        or row["pre_result_family_lock_sha256"] != family_sha
        for row in rows
    ):
        raise C.ContractError(f"rows do not bind one pre-result family for {layer}")
    validated = validate_supplement_rows(
        rows,
        layer,
        primary_denominator=family["primary_correction_denominator"],
        control_denominator=family["control_correction_denominator"],
    )
    source_results = lock.get("source_results")
    if not isinstance(source_results, list) or not source_results:
        raise C.ContractError(f"{layer} lock has no source result identities")
    source_by_path: dict[str, str] = {}
    for item in source_results:
        if not isinstance(item, dict) or set(item) != {"path", "sha256"}:
            raise C.ContractError(f"invalid source-result identity in {layer}")
        path, expected = item["path"], C.require_sha256(item["sha256"], "source result")
        observed = C.stable_identity(root, path)
        if observed["sha256"] != expected or path in source_by_path:
            raise C.ContractError(f"source-result hash drift or duplication in {layer}")
        source_by_path[path] = expected
    for row in rows:
        if source_by_path.get(row["source_result_path"]) != row["source_result_sha256"]:
            raise C.ContractError(f"unbound row source result in {layer}")
    return validated, table_identity, lock_identity


def _output_classification(layer: str, status: str, supported_count: int) -> tuple[str, str]:
    if layer == "GLOBAL_RG_CONTEXT":
        if status == "SUPPORTED":
            return status, "GLOBAL_RG_CONTEXT_SUPPORTED"
        if status == "TESTED_NO_SUPPORT":
            return status, "GLOBAL_RG_CONTEXT_NO_DETECTED_SUPPORT"
        raise C.ContractError("global rg cannot have a non-tested output state")
    prefix = LAYER_PREFIX[layer]
    if status == "SUPPORTED":
        suffix = "SUPPORTED_HYPOTHESIS" if prefix == "GENE_MECHANISM" else "SUPPORTED"
    elif status == "TESTED_NO_SUPPORT":
        suffix = "TESTED_NO_SUPPORT"
    elif status == "NOT_APPLICABLE":
        suffix = "NOT_APPLICABLE"
    elif status == "BLOCKED":
        suffix = "BLOCKED"
    elif status == "FAILED":
        suffix = "FAILED"
    else:
        raise C.ContractError(f"invalid output status: {status}")
    if status == "SUPPORTED" and supported_count <= 0:
        raise C.ContractError("SUPPORTED classification has no supported entity")
    return status, f"{prefix}_{suffix}"


def _output_status(upstream: str, supported_count: int) -> str:
    if upstream == "COMPLETE":
        return "SUPPORTED" if supported_count else "TESTED_NO_SUPPORT"
    return {
        "NOT_APPLICABLE_UPSTREAM_ZERO": "NOT_APPLICABLE",
        "BLOCKED_UPSTREAM": "BLOCKED",
        "FAILED_UPSTREAM": "FAILED",
    }[upstream]


def _global_rows(
    policy: Mapping[str, Any], core: Mapping[str, Any],
) -> list[dict[str, str]]:
    by_pair = {
        (row["sleep_trait"], row["non_sleep_trait"]): row
        for row in core["atlas_pair_rows"]
    }
    source_path = policy["immutable_core"]["atlas_pairs"]["path"]
    source_sha = core["atlas_pair_identity"]["sha256"]
    layer = policy["layer_contract"]["GLOBAL_RG_CONTEXT"]
    result: list[dict[str, str]] = []
    for pair in PAIR_ORDER:
        anchor, external, role = PAIR_IDENTITIES[pair]
        for sleep in SLEEP_TRAITS:
            source = by_pair.get((sleep, external))
            if source is None:
                raise C.ContractError(f"locked global-rg row is absent: {sleep}/{external}")
            q_value = _probability(source["global_rg_fdr_all_396"], "global rg FDR")
            _probability(source["global_rg_p"], "global rg P")
            _finite_or_na(source["global_rg"], "global rg")
            standard_error = _finite_or_na(source["global_rg_se"], "global rg SE")
            if standard_error is None or standard_error <= 0:
                raise C.ContractError(f"invalid global-rg standard error: {sleep}/{external}")
            supported = bool(q_value is not None and q_value <= 0.05)
            status = "SUPPORTED" if supported else "TESTED_NO_SUPPORT"
            classification = (
                "GLOBAL_RG_CONTEXT_SUPPORTED"
                if supported else "GLOBAL_RG_CONTEXT_NO_DETECTED_SUPPORT"
            )
            base = {
                "analysis_id": ANALYSIS_ID,
                "pair_id": pair,
                "family_role": role,
                "anchor_sleep_trait": anchor,
                "external_trait": external,
                "queried_sleep_trait": sleep,
                "evidence_level": "GLOBAL_RG_CONTEXT",
                "record_kind": "COVERAGE",
                "entity_id": "GENOME_WIDE",
                "evidence_status": status,
                "classification": classification,
                "supported_entity_count": "1" if supported else "0",
                "supported_entity_ids": "GENOME_WIDE" if supported else NA,
                "estimate": source["global_rg"],
                "standard_error": source["global_rg_se"],
                "p_value": source["global_rg_p"],
                "adjusted_p_value": source["global_rg_fdr_all_396"],
                "effect_direction": source["effect_direction"],
                "power_qc_status": (
                    f"{source['analysis_tier']};{source['interpretation_status']}"
                ),
                "tested_hypothesis_count": "1",
                "failed_hypothesis_count": "0",
                "upstream_correction_method": layer["source_correction"],
                "upstream_correction_denominator": "396",
                "phase24_correction_denominator": "0",
                "source_result_path": source_path,
                "source_result_sha256": str(source_sha),
                "claim_limit": layer["claim_limit"],
            }
            output = {"analysis_id": base["analysis_id"], "evidence_row_id": ""}
            output.update({key: base[key] for key in OUTPUT_FIELDS if key not in output})
            output["evidence_row_id"] = evidence_row_id(output)
            result.append(output)
    return result


def _supplement_output_rows(
    policy: Mapping[str, Any], rows: Sequence[dict[str, str]], layer: str,
) -> list[dict[str, str]]:
    claim = policy["layer_contract"][layer]["claim_limit"]
    result: list[dict[str, str]] = []
    for row in rows:
        supported_count = _integer(row["supported_entity_count"], "supported entity count")
        status = _output_status(row["upstream_status"], supported_count)
        if row["record_kind"] == "EVIDENCE":
            status = "SUPPORTED"
            supported_count = 1
        status, classification = _output_classification(layer, status, supported_count)
        base = {
            "analysis_id": ANALYSIS_ID,
            "pair_id": row["pair_id"],
            "family_role": row["family_role"],
            "anchor_sleep_trait": row["anchor_sleep_trait"],
            "external_trait": row["external_trait"],
            "queried_sleep_trait": row["queried_sleep_trait"],
            "evidence_level": layer,
            "record_kind": row["record_kind"],
            "entity_id": row["entity_id"],
            "evidence_status": status,
            "classification": classification,
            "supported_entity_count": str(supported_count),
            "supported_entity_ids": row["supported_entity_ids"],
            "estimate": row["estimate"],
            "standard_error": row["standard_error"],
            "p_value": row["p_value"],
            "adjusted_p_value": row["adjusted_p_value"],
            "effect_direction": row["effect_direction"],
            "power_qc_status": row["power_qc_status"],
            "tested_hypothesis_count": row["tested_hypothesis_count"],
            "failed_hypothesis_count": row["failed_hypothesis_count"],
            "upstream_correction_method": row["upstream_correction_method"],
            "upstream_correction_denominator": row["upstream_correction_denominator"],
            "phase24_correction_denominator": "0",
            "source_result_path": row["source_result_path"],
            "source_result_sha256": row["source_result_sha256"],
            "claim_limit": claim,
        }
        output = {"analysis_id": base["analysis_id"], "evidence_row_id": ""}
        output.update({key: base[key] for key in OUTPUT_FIELDS if key not in output})
        output["evidence_row_id"] = evidence_row_id(output)
        result.append(output)
    return result


def classify_sleep_pattern(supported: set[str]) -> tuple[str, list[str]]:
    """Classify one complete non-global entity over all twelve sleep traits."""

    if not supported or not supported.issubset(SLEEP_TRAITS):
        raise C.ContractError("sleep-pattern classification requires supported locked traits")
    if supported == {"snoring"}:
        classification = "SNORING_SPECIFIC"
    elif supported == {"insomnia"}:
        classification = "INSOMNIA_SPECIFIC"
    elif supported == {"snoring", "sleep_apnea"}:
        classification = "SLEEP_BREATHING_SHARED"
    elif (
        supported.issubset(DIMENSIONS["DISRUPTION"])
        and "insomnia" in supported
        and len(supported) >= 2
    ):
        classification = "SLEEP_DISRUPTION_SHARED"
    else:
        dimensions = [
            name for name, members in DIMENSIONS.items() if supported.intersection(members)
        ]
        classification = "MULTIPLE_SLEEP_DIMENSIONS" if len(dimensions) >= 2 else "NONSPECIFIC"
    dimensions = [
        name for name, members in DIMENSIONS.items() if supported.intersection(members)
    ]
    return classification, dimensions


def build_generalization_rows(
    policy: Mapping[str, Any], rows: Sequence[dict[str, str]], pair: str,
) -> list[dict[str, str]]:
    anchor, external, role = PAIR_IDENTITIES[pair]
    claim = policy["mechanism_level_generalization"]["claim_limit"]
    result: list[dict[str, str]] = []
    for layer in SUPPLEMENT_LAYERS:
        relevant = [
            row for row in rows
            if row["pair_id"] == pair and row["evidence_level"] == layer
        ]
        coverage = {
            row["queried_sleep_trait"]: row
            for row in relevant if row["record_kind"] == "COVERAGE"
        }
        if list(coverage) != list(SLEEP_TRAITS):
            raise C.ContractError(f"generalization lacks ordered twelve-trait coverage: {pair}/{layer}")
        states = [coverage[trait]["evidence_status"] for trait in SLEEP_TRAITS]
        if any(state in {"BLOCKED", "FAILED"} for state in states):
            classification, coverage_status = "BLOCKED", "BLOCKED_OR_FAILED_QUERY_PRESENT"
            entities: dict[str, set[str]] = {}
        elif any(state == "NOT_APPLICABLE" for state in states):
            classification, coverage_status = "INCONCLUSIVE", "INCOMPLETE_COMPARABLE_COVERAGE"
            entities = {}
        elif any(state not in {"SUPPORTED", "TESTED_NO_SUPPORT"} for state in states):
            raise C.ContractError(f"unexpected coverage state: {pair}/{layer}")
        else:
            classification, coverage_status = "", "COMPLETE_12_OF_12"
            entities = defaultdict(set)
            for row in relevant:
                if row["record_kind"] == "EVIDENCE" and row["evidence_status"] == "SUPPORTED":
                    entities[row["entity_id"]].add(row["queried_sleep_trait"])
        if not entities:
            if not classification:
                classification = "NO_DETECTED_SUPPORT"
            result.append({
                "analysis_id": ANALYSIS_ID,
                "pair_id": pair,
                "family_role": role,
                "anchor_sleep_trait": anchor,
                "external_trait": external,
                "evidence_level": layer,
                "entity_id": NA,
                "supported_sleep_trait_count": "0",
                "supported_sleep_traits": NA,
                "supported_dimension_count": "0",
                "supported_dimensions": NA,
                "generalization_classification": classification,
                "coverage_status": coverage_status,
                "source_evidence_row_ids": NA,
                "claim_limit": claim,
            })
            continue
        for entity in sorted(entities):
            supported = entities[entity]
            pattern, dimensions = classify_sleep_pattern(supported)
            source_ids = sorted(
                row["evidence_row_id"] for row in relevant
                if row["record_kind"] == "EVIDENCE" and row["entity_id"] == entity
            )
            result.append({
                "analysis_id": ANALYSIS_ID,
                "pair_id": pair,
                "family_role": role,
                "anchor_sleep_trait": anchor,
                "external_trait": external,
                "evidence_level": layer,
                "entity_id": entity,
                "supported_sleep_trait_count": str(len(supported)),
                "supported_sleep_traits": ";".join(
                    trait for trait in SLEEP_TRAITS if trait in supported
                ),
                "supported_dimension_count": str(len(dimensions)),
                "supported_dimensions": ";".join(dimensions),
                "generalization_classification": pattern,
                "coverage_status": coverage_status,
                "source_evidence_row_ids": ";".join(source_ids),
                "claim_limit": claim,
            })
    return result


def build_phase24_rows(
    policy: Mapping[str, Any],
    core: Mapping[str, Any],
    supplements: Mapping[str, Sequence[dict[str, str]]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    if set(supplements) != set(SUPPLEMENT_LAYERS):
        raise C.ContractError("all three non-global layers are required; global rg alone is insufficient")
    rows = _global_rows(policy, core)
    for layer in SUPPLEMENT_LAYERS:
        rows.extend(_supplement_output_rows(policy, supplements[layer], layer))
    pair_rank = {pair: index for index, pair in enumerate(PAIR_ORDER)}
    sleep_rank = {trait: index for index, trait in enumerate(SLEEP_TRAITS)}
    layer_rank = {layer: index for index, layer in enumerate(LAYERS)}
    kind_rank = {"COVERAGE": 0, "EVIDENCE": 1}
    rows.sort(key=lambda row: (
        pair_rank[row["pair_id"]], sleep_rank[row["queried_sleep_trait"]],
        layer_rank[row["evidence_level"]], kind_rank[row["record_kind"]], row["entity_id"],
    ))
    ids = [row["evidence_row_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise C.ContractError("Phase24 evidence row IDs collide")
    primary = [row for row in rows if row["pair_id"] != "CONTROL"]
    control = [row for row in rows if row["pair_id"] == "CONTROL"]
    primary_coverage = sum(row["record_kind"] == "COVERAGE" for row in primary)
    control_coverage = sum(row["record_kind"] == "COVERAGE" for row in control)
    if primary_coverage != 96 or control_coverage != 48:
        raise C.ContractError("Phase24 fixed coverage family is not 96 primary + 48 control")
    general_primary = []
    for pair in ("A", "B"):
        general_primary.extend(build_generalization_rows(policy, rows, pair))
    general_control = build_generalization_rows(policy, rows, "CONTROL")
    return primary, control, general_primary, general_control


def preflight(
    root: Path = ROOT, *, deep_hash_dense: bool = False,
) -> dict[str, Any]:
    """Side-effect-free production preflight; raises when upstream is absent."""

    policy, policy_identity = validate_policy(root)
    _contract, contract_identity = validate_contract_lock(root, policy_identity)
    core = validate_core(root, policy, deep_hash_dense=deep_hash_dense)
    missing: list[str] = []
    for layer in SUPPLEMENT_LAYERS:
        definition = policy["layer_contract"][layer]
        for key in ("future_input", "future_input_lock"):
            relative = definition[key]
            try:
                path = C.safe_path(root, relative, "future upstream artifact")
                info = os.lstat(path)
            except FileNotFoundError:
                missing.append(relative)
                continue
            if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_size <= 0:
                raise C.ContractError(f"upstream artifact is not a real non-empty file: {relative}")
    if missing:
        raise C.UpstreamBlocked(
            "Phase24 requires all non-global upstream families; missing: " + ", ".join(missing)
        )
    supplements: dict[str, list[dict[str, str]]] = {}
    inputs: dict[str, dict[str, object]] = {}
    for layer in SUPPLEMENT_LAYERS:
        layer_rows, table_identity, lock_identity = validate_supplement(
            root, policy, policy_identity, layer,
        )
        supplements[layer] = layer_rows
        inputs[f"{layer}:table"] = table_identity
        inputs[f"{layer}:lock"] = lock_identity
    primary, control, general_primary, general_control = build_phase24_rows(
        policy, core, supplements,
    )
    return {
        "schema_version": "sleep-atlas-track-b-cross-sleep-preflight.1",
        "analysis_id": ANALYSIS_ID,
        "status": "READY_TO_EXECUTE",
        "side_effect_free": True,
        "policy": policy,
        "policy_identity": policy_identity,
        "contract_identity": contract_identity,
        "core": core,
        "upstream_inputs": inputs,
        "supplements": supplements,
        "primary_rows": primary,
        "control_rows": control,
        "primary_generalization_rows": general_primary,
        "control_generalization_rows": general_control,
    }


def execute(
    root: Path = ROOT, *, deep_hash_dense: bool = True,
) -> dict[str, Any]:
    state = preflight(root, deep_hash_dense=deep_hash_dense)
    policy = state["policy"]
    output = policy["output_contract"]
    primary_bytes = C.table_bytes(OUTPUT_FIELDS, state["primary_rows"])
    control_bytes = C.table_bytes(OUTPUT_FIELDS, state["control_rows"])
    general_primary_bytes = C.table_bytes(
        GENERALIZATION_FIELDS, state["primary_generalization_rows"],
    )
    general_control_bytes = C.table_bytes(
        GENERALIZATION_FIELDS, state["control_generalization_rows"],
    )
    input_identities = {
        key: value for key, value in sorted(state["upstream_inputs"].items())
    }
    fingerprint = C.digest_json({
        "policy": state["policy_identity"]["sha256"],
        "contract": state["contract_identity"]["sha256"],
        "inputs": input_identities,
    })
    output_payloads = {
        output["primary_path"]: primary_bytes,
        output["control_path"]: control_bytes,
        output["primary_generalization_path"]: general_primary_bytes,
        output["control_generalization_path"]: general_control_bytes,
    }
    provenance = {
        "schema_version": "sleep-atlas-track-b-cross-sleep-provenance.1",
        "analysis_id": ANALYSIS_ID,
        "execution_fingerprint": fingerprint,
        "contract_frozen_before_results": True,
        "phase24_new_hypothesis_count_primary": 0,
        "phase24_new_hypothesis_count_control": 0,
        "global_rg_source_correction_denominator": 396,
        "global_rg_is_context_only": True,
        "global_rg_alone_cannot_satisfy_upstream_gate": True,
        "all_twelve_sleep_traits_queried": True,
        "unrestricted_phenome_scan_performed": False,
        "control_published_separately": True,
        "policy": state["policy_identity"],
        "contract": state["contract_identity"],
        "implementation": {
            str(SCRIPT_REL): C.stable_identity(root, SCRIPT_REL),
            str(COMMON_REL): C.stable_identity(root, COMMON_REL),
        },
        "upstream_inputs": input_identities,
        "output_rows": {
            output["primary_path"]: len(state["primary_rows"]),
            output["control_path"]: len(state["control_rows"]),
            output["primary_generalization_path"]: len(state["primary_generalization_rows"]),
            output["control_generalization_path"]: len(state["control_generalization_rows"]),
        },
        "outputs": {
            path: {"bytes": len(payload), "sha256": C.digest_bytes(payload)}
            for path, payload in output_payloads.items()
        },
        "claim_limit": policy["claim_limits"]["phase"],
    }
    provenance_bytes = json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n"
    payloads = [(path, payload) for path, payload in output_payloads.items()]
    payloads.append((output["provenance_path"], provenance_bytes))
    publication = C.no_replace_publish(root, payloads)
    return {
        "status": "COMPLETE_SEALED",
        "execution_fingerprint": fingerprint,
        "publication": publication,
        "outputs": provenance["outputs"],
        "provenance": output["provenance_path"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--preflight", action="store_true")
    modes.add_argument("--execute", action="store_true")
    parser.add_argument("--deep-hash-dense", action="store_true")
    parser.add_argument(
        "--acknowledge-production-science",
        choices=["PHASE24_PRODUCTION"],
        help="required only with --execute",
    )
    args = parser.parse_args()
    root = Path(args.root).resolve()
    try:
        if args.execute:
            if args.acknowledge_production_science != "PHASE24_PRODUCTION":
                raise C.ContractError(
                    "--execute requires --acknowledge-production-science PHASE24_PRODUCTION"
                )
            result = execute(root, deep_hash_dense=True)
        else:
            result = preflight(root, deep_hash_dense=args.deep_hash_dense)
            result = {
                "schema_version": result["schema_version"],
                "analysis_id": ANALYSIS_ID,
                "status": result["status"],
                "side_effect_free": True,
                "primary_rows_if_executed": len(result["primary_rows"]),
                "control_rows_if_executed": len(result["control_rows"]),
            }
    except C.UpstreamBlocked as error:
        print(json.dumps({
            "analysis_id": ANALYSIS_ID, "status": "BLOCKED_UPSTREAM",
            "side_effect_free": True, "reason": str(error),
        }, sort_keys=True))
        return 3
    except C.ContractError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
