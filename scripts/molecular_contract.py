#!/usr/bin/env python3
"""Shared immutable-contract checks for molecular-QTL and TWAS production."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import fine_mapping_contract


SEARCH_FIELDS = [
    "search_task_id", "locus_id", "pair_id", "sleep_trait", "non_sleep_trait",
    "chromosome_grch37", "start_grch37", "end_grch37", "chromosome_grch38",
    "start_grch38", "end_grch38", "source_family_id", "study_id", "dataset_id",
    "modality", "quant_method", "context", "sample_size", "source_build", "query_mode",
    "source_url", "source_index_url", "exact_release", "planned_outcome",
    "results_accessed_before_lock",
]
SEARCH_COVERAGE_FIELDS = [
    "search_task_id", "locus_id", "source_family_id", "study_id", "dataset_id",
    "modality", "context", "search_outcome", "normalized_row_count", "feature_count",
    "analyzable_feature_count", "query_provenance_path", "query_provenance_sha256",
]
FEATURE_FIELDS = [
    "comparison_id", "search_task_id", "locus_id", "pair_id", "trait_id", "trait_role",
    "source_family_id", "study_id", "dataset_id", "exact_release", "modality",
    "feature_id", "feature_name", "gene_id", "gene_symbol", "context", "molecular_N",
    "shared_variant_count", "normalized_qtl_path", "normalized_qtl_sha256",
    "trait_task_path", "trait_task_sha256", "trait_task_lock_path", "trait_task_lock_sha256",
    "trait_molecular_results_accessed_before_feature_lock",
]
FEATURE_EXCLUSION_FIELDS = [
    "search_task_id", "locus_id", "source_family_id", "dataset_id", "modality",
    "feature_id", "shared_variant_count", "exclusion_reason", "query_provenance_path",
]
MODEL_REGISTRY_FIELDS = [
    "model_id", "model_family", "modality", "context", "model_db_path", "model_db_sha256",
    "covariance_path", "covariance_sha256", "model_snp_key", "gene_count", "weight_count",
    "phi_gene_count", "phi_missing_gene_count", "phi_min", "phi_max", "source_release",
    "model_file_id", "covariance_file_id",
]
PHI_EXCLUSION_FIELDS = ["model_id", "model_family", "context", "gene_id", "exclusion_reason"]
MODEL_VARIANT_FIELDS = ["model_variant_id", "chromosome_grch38", "position_grch38", "ref", "alt"]
TWAS_MAPPING_FIELDS = [
    "mapping_id", "trait_id", "model_family", "modality", "full_gwas_path", "full_gwas_sha256",
    "model_variant_path", "model_variant_sha256", "mapped_gwas_path", "mapped_gwas_lock_path",
    "results_accessed_before_lock",
]
TWAS_ELIGIBILITY_FIELDS = [
    "trait_id", "trait_domain", "trait_type", "gwas_N", "gwas_h2", "gwas_h2_scale",
    "analysis_status", "reason", "results_accessed_before_lock",
]
TWAS_RUN_FIELDS = [
    "run_id", "trait_id", "trait_domain", "model_id", "model_family", "modality", "context",
    "mapping_id", "mapped_gwas_path", "mapped_gwas_lock_path", "model_db_path", "model_db_sha256",
    "covariance_path", "covariance_sha256", "model_snp_key", "gwas_N", "gwas_h2",
    "gwas_h2_scale", "variance_control_status", "output_path", "results_accessed_before_lock",
]
TWAS_RESULT_FIELDS = [
    "twas_id", "trait_id", "model_family", "modality", "context", "gene_id", "gene_symbol",
    "zscore", "uncalibrated_zscore", "effect_size", "p_value", "uncalibrated_p_value", "fdr",
    "n_snps_used", "n_snps_in_model", "gwas_N", "gwas_h2", "gwas_h2_scale",
    "coverage_fraction", "status", "model_id", "provenance_id",
]
TWAS_COVERAGE_FIELDS = [
    "run_id", "trait_id", "model_id", "model_family", "context", "analysis_status", "reason",
    "result_rows", "provenance_path", "provenance_sha256",
]
COMMON_SCRIPTS = (
    "scripts/61_molecular_preflight.py", "scripts/73_collate_molecular.py",
    "scripts/liftover_chain.py", "scripts/molecular_contract.py",
    "scripts/fine_mapping_contract.py",
)
QTL_SCRIPTS = COMMON_SCRIPTS + (
    "scripts/61_fetch_molecular_metadata.sh", "scripts/62_prepare_molecular_search_plan.py",
    "scripts/63_query_eqtl_catalogue.py", "scripts/63_record_molecular_search.py",
    "scripts/63_run_molecular_search.py", "scripts/64_lock_molecular_features.py",
    "scripts/65_materialize_molecular_coloc.py", "scripts/66_run_molecular_coloc.py",
    "discovery_extension/scripts/35_run_susie_coloc.R",
)
TWAS_SCRIPTS = COMMON_SCRIPTS + (
    "scripts/67_fetch_twas_resources.sh", "scripts/67_lock_twas_model_inventory.py",
    "scripts/68_index_twas_models.py", "scripts/69_prepare_twas_manifest.py",
    "scripts/70_materialize_twas_gwas.py", "scripts/71_run_twas.py",
    "scripts/72_collate_twas.py",
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty checksum-bound artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, expected_fields: list[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty tabular artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = list(reader.fieldnames or []), list(reader)
    if expected_fields is not None and fields != expected_fields:
        fail(f"wrong schema for {path}: {fields}")
    return fields, rows


def safe_relative(value: str, label: str) -> Path:
    path = Path(value)
    if not value or path.is_absolute() or ".." in path.parts:
        fail(f"unsafe {label} path: {value}")
    return path


def script_hashes(root: Path, family: str) -> dict[str, str]:
    paths = QTL_SCRIPTS if family == "qtl" else TWAS_SCRIPTS if family == "twas" else None
    if paths is None:
        fail(f"unknown molecular script family: {family}")
    return {relative: sha256(root / relative) for relative in paths}


def load_policy(root: Path, relative: str = "config/molecular_analysis_policy.json") -> tuple[Path, dict[str, object]]:
    path = root / relative
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"molecular policy is unreadable: {exc}")
    if (
        policy.get("analysis_id") != "atlas-v1.0-molecular"
        or policy.get("schema_version") != "1.0.0"
        or policy.get("required_modalities") != ["eQTL", "sQTL", "pQTL_or_PWAS", "TWAS"]
        or policy.get("twas", {}).get("trait_family") != "all_45_locked_traits"
        or policy.get("twas", {}).get("phi_model_source", {}).get("expected_context_count") != 49
    ):
        fail("molecular policy differs from the frozen atlas scope")
    return path, policy


def validate_fine_mapping(root: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], dict[str, object]]:
    loci_path = root / "results/atlas/loci.tsv"
    variants_path = root / "results/atlas/variants.tsv"
    coloc_path = root / "results/tables/trait_trait_colocalization.tsv"
    provenance_path = root / "results/atlas/fine_mapping.provenance.json"
    _, loci = read_tsv(loci_path, fine_mapping_contract.CANONICAL_LOCUS_FIELDS)
    _, variants = read_tsv(variants_path, fine_mapping_contract.CANONICAL_VARIANT_FIELDS)
    _, coloc = read_tsv(coloc_path, fine_mapping_contract.ENGINE_OUTPUT_FIELDS["colocalization.tsv"])
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"fine-mapping provenance is unreadable: {exc}")
    outputs = provenance.get("outputs", {})
    if (
        provenance.get("schema_version") != "atlas-v1.0-finemapping-canonical.2"
        or outputs.get("results/atlas/loci.tsv") != sha256(loci_path)
        or outputs.get("results/atlas/variants.tsv") != sha256(variants_path)
        or outputs.get("results/tables/trait_trait_colocalization.tsv") != sha256(coloc_path)
        or provenance.get("locus_count") != len(loci)
        or provenance.get("variant_count") != len(variants)
        or provenance.get("colocalization_row_count") != len(coloc)
        or provenance.get("zero_family_not_applicable") != (len(loci) == 0)
        or provenance.get("script_sha256") != fine_mapping_contract.script_hashes(root)
    ):
        fail("fine-mapping canonical family differs from provenance")
    locus_ids = [row["locus_id"] for row in loci]
    if len(set(locus_ids)) != len(locus_ids) or any(row["analysis_tier"] != "PRIMARY_PHASE1" for row in loci):
        fail("fine-mapping locus family is duplicated or outside PRIMARY_PHASE1")
    variant_loci = {row["locus_id"] for row in variants}
    if not set(locus_ids).issubset(variant_loci) or variant_loci - set(locus_ids):
        fail("fine-mapping variants differ from the canonical locus family")
    return loci, variants, coloc, provenance


def validate_preflight(
    root: Path, policy_path: Path, path: Path, *, purpose: str,
) -> dict[str, object]:
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"molecular preflight is unreadable: {exc}")
    sources_path = root / "config/molecular_source_registry.tsv"
    readiness_path = root / safe_relative(str(report.get("readiness_path", "")), "molecular readiness")
    loci, _, _, fine = validate_fine_mapping(root)
    ready_field = "qtl_planning_ready" if purpose == "qtl_plan" else "twas_production_ready"
    if (
        report.get("policy_sha256") != sha256(policy_path)
        or report.get("sources_sha256") != sha256(sources_path)
        or report.get("readiness_sha256") != sha256(readiness_path)
        or report.get("fine_mapping_provenance_sha256")
        != sha256(root / "results/atlas/fine_mapping.provenance.json")
        or report.get("fine_mapping_locus_count") != len(loci)
        or report.get("qtl_analysis_required") != (len(loci) > 0)
        or report.get("script_sha256", {}).get("qtl") != script_hashes(root, "qtl")
        or report.get("script_sha256", {}).get("twas") != script_hashes(root, "twas")
        or report.get(ready_field) is not True
        or fine.get("locus_count") != len(loci)
    ):
        fail(f"molecular preflight is not ready for {purpose} or has drifted")
    return report


def validate_search_plan(
    root: Path, plan_path: Path, lock_path: Path, policy_path: Path,
    preflight_path: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    _, plan = read_tsv(plan_path, SEARCH_FIELDS)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"molecular search-plan lock is unreadable: {exc}")
    loci, variants, _, _ = validate_fine_mapping(root)
    sources_path = root / "config/molecular_source_registry.tsv"
    ids = [row["search_task_id"] for row in plan]
    if (
        lock.get("schema_version") != "atlas-v1.0-molecular-search.2"
        or lock.get("results_accessed_before_lock") is not False
        or lock.get("plan_sha256") != sha256(plan_path)
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("sources_sha256") != sha256(sources_path)
        or lock.get("preflight_sha256") != sha256(preflight_path)
        or lock.get("loci_sha256") != sha256(root / "results/atlas/loci.tsv")
        or lock.get("variants_sha256") != sha256(root / "results/atlas/variants.tsv")
        or lock.get("fine_mapping_provenance_sha256")
        != sha256(root / "results/atlas/fine_mapping.provenance.json")
        or lock.get("locus_count") != len(loci)
        or lock.get("search_task_count") != len(plan)
        or lock.get("search_task_ids_in_locked_order") != ids
        or len(ids) != len(set(ids))
        or lock.get("script_sha256") != script_hashes(root, "qtl")
    ):
        fail("molecular search plan differs from its immutable lock")
    locus_ids = {row["locus_id"] for row in loci}
    if any(row["locus_id"] not in locus_ids or row["results_accessed_before_lock"] != "NO" for row in plan):
        fail("molecular search plan contains an out-of-family or post-result task")
    if not loci and (plan or lock.get("dataset_count_per_locus") != 0):
        fail("zero-locus molecular search plan is not exactly empty")
    del variants
    return plan, lock


def validate_feature_family(
    root: Path, manifest_path: Path, lock_path: Path, policy_path: Path,
    preflight_path: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    plan_path = root / "results/tables/molecular_search_plan.tsv"
    plan_lock_path = root / "results/tables/molecular_search_plan.lock.json"
    plan, _ = validate_search_plan(
        root, plan_path, plan_lock_path, policy_path, preflight_path,
    )
    coverage_path = root / "results/tables/molecular_search_coverage.tsv"
    exclusions_path = root / "results/tables/molecular_feature_exclusions.tsv"
    _, coverage = read_tsv(coverage_path, SEARCH_COVERAGE_FIELDS)
    _, manifest = read_tsv(manifest_path, FEATURE_FIELDS)
    read_tsv(exclusions_path, FEATURE_EXCLUSION_FIELDS)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"molecular feature lock is unreadable: {exc}")
    ids = [row["comparison_id"] for row in manifest]
    if (
        lock.get("schema_version") != "atlas-v1.0-molecular-features.2"
        or lock.get("trait_molecular_results_accessed_before_feature_lock") is not False
        or lock.get("coverage_sha256") != sha256(coverage_path)
        or lock.get("manifest_sha256") != sha256(manifest_path)
        or lock.get("exclusions_sha256") != sha256(exclusions_path)
        or lock.get("search_plan_sha256") != sha256(plan_path)
        or lock.get("search_plan_lock_sha256") != sha256(plan_lock_path)
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("search_task_count") != len(plan)
        or lock.get("analyzable_feature_comparison_count") != len(manifest)
        or lock.get("comparison_ids_in_locked_order") != ids
        or len(ids) != len(set(ids))
        or lock.get("zero_locus_qtl_not_applicable") != (len(plan) == 0)
        or lock.get("script_sha256") != script_hashes(root, "qtl")
        or [row["search_task_id"] for row in coverage]
        != [row["search_task_id"] for row in plan]
    ):
        fail("molecular feature family differs from its immutable pre-result lock")
    if not plan and (coverage or manifest or lock.get("search_input_hashes") != {}):
        fail("zero-locus molecular feature family is not exactly empty")
    return manifest, lock


def validate_model_registry(
    root: Path, registry_path: Path, lock_path: Path, policy_path: Path,
    phi_exclusions_path: Path | None = None,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    _, models = read_tsv(registry_path, MODEL_REGISTRY_FIELDS)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"TWAS model-registry lock is unreadable: {exc}")
    _, policy = load_policy(root, str(policy_path.relative_to(root)))
    source = policy["twas"]["phi_model_source"]
    inventory_path = root / source["inventory_path"]
    inventory_lock_path = root / source["inventory_lock_path"]
    download_lock_path = root / source["download_lock_path"]
    model_ids = [row["model_id"] for row in models]
    if (
        lock.get("schema_version") != "atlas-v1.0-twas-models.2"
        or lock.get("results_accessed") is not False
        or lock.get("registry_sha256") != sha256(registry_path)
        or lock.get("model_count") != len(models)
        or lock.get("model_ids_in_locked_order") != model_ids
        or len(model_ids) != len(set(model_ids))
        or lock.get("inventory_sha256") != sha256(inventory_path)
        or lock.get("inventory_lock_sha256") != sha256(inventory_lock_path)
        or lock.get("download_lock_sha256") != sha256(download_lock_path)
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("script_sha256") != script_hashes(root, "twas")
    ):
        fail("TWAS model registry differs from its immutable lock")
    expected_families = set(policy["twas"]["model_families"])
    contexts: dict[str, set[str]] = {}
    for row in models:
        if row["modality"] != "TWAS":
            fail(f"invalid TWAS model registry row: {row['model_id']}")
        contexts.setdefault(row["model_family"], set()).add(row["context"])
        for path_field, hash_field in (("model_db_path", "model_db_sha256"), ("covariance_path", "covariance_sha256")):
            path = root / safe_relative(row[path_field], path_field)
            if sha256(path) != row[hash_field]:
                fail(f"TWAS model file differs from registry: {row['model_id']}/{path_field}")
    expected_context_count = int(source["expected_context_count"])
    if set(contexts) != expected_families or any(len(values) != expected_context_count for values in contexts.values()):
        fail("TWAS model registry does not cover the exact locked family/context grid")
    for family in sorted(expected_families):
        variant_path = root / "results/tables/twas_model_variants" / f"{family}.tsv.gz"
        if lock.get("variant_registry_sha256", {}).get(family) != sha256(variant_path):
            fail(f"TWAS model-variant registry differs from lock: {family}")
    if phi_exclusions_path is not None:
        _, exclusions = read_tsv(phi_exclusions_path, PHI_EXCLUSION_FIELDS)
        if (
            lock.get("phi_exclusions_sha256") != sha256(phi_exclusions_path)
            or lock.get("phi_exclusion_count") != len(exclusions)
        ):
            fail("TWAS phi-exclusion family differs from model lock")
    return models, lock


def validate_twas_manifest(
    root: Path, mapping_path: Path, eligibility_path: Path, run_path: Path,
    lock_path: Path, policy_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], dict[str, object]]:
    _, mappings = read_tsv(mapping_path, TWAS_MAPPING_FIELDS)
    _, eligibility = read_tsv(eligibility_path, TWAS_ELIGIBILITY_FIELDS)
    _, runs = read_tsv(run_path, TWAS_RUN_FIELDS)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"TWAS run-manifest lock is unreadable: {exc}")
    run_ids = [row["run_id"] for row in runs]
    mapping_ids = [row["mapping_id"] for row in mappings]
    trait_ids = [row["trait_id"] for row in eligibility]
    preflight_path = root / "results/tables/twas_preflight.json"
    atlas_traits_path = root / "results/atlas/traits.tsv"
    model_registry_path = root / "results/tables/twas_model_registry.tsv"
    model_lock_path = root / "results/tables/twas_model_registry.lock.json"
    panel_ids = [row["trait_id"] for row in read_tsv(root / "config/analysis_panel.tsv")[1]]
    eligible_count = sum(row["analysis_status"] == "ELIGIBLE" for row in eligibility)
    if (
        lock.get("schema_version") != "atlas-v1.0-twas-runs.2"
        or lock.get("results_accessed_before_lock") is not False
        or lock.get("trait_count") != 45
        or len(eligibility) != 45 or len(trait_ids) != len(set(trait_ids))
        or lock.get("mapping_count") != len(mappings)
        or lock.get("run_count") != len(runs)
        or lock.get("run_ids_in_locked_order") != run_ids
        or len(run_ids) != len(set(run_ids)) or len(mapping_ids) != len(set(mapping_ids))
        or lock.get("mapping_sha256") != sha256(mapping_path)
        or lock.get("eligibility_sha256") != sha256(eligibility_path)
        or lock.get("run_manifest_sha256") != sha256(run_path)
        or lock.get("eligible_trait_count") != eligible_count
        or lock.get("not_applicable_trait_count") != 45 - eligible_count
        or lock.get("model_registry_sha256") != sha256(model_registry_path)
        or lock.get("model_registry_lock_sha256") != sha256(model_lock_path)
        or lock.get("atlas_traits_sha256") != sha256(atlas_traits_path)
        or lock.get("preflight_sha256") != sha256(preflight_path)
        or lock.get("fine_mapping_provenance_sha256")
        != sha256(root / "results/atlas/fine_mapping.provenance.json")
        or trait_ids != panel_ids
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("script_sha256") != script_hashes(root, "twas")
    ):
        fail("TWAS run family differs from its immutable pre-result lock")
    validate_preflight(root, policy_path, preflight_path, purpose="twas")
    models, _ = validate_model_registry(
        root, model_registry_path, model_lock_path, policy_path,
        root / "results/tables/twas_model_phi_exclusions.tsv",
    )
    if lock.get("model_count") != len(models) or len(runs) != eligible_count * len(models):
        fail("TWAS run family does not cover the exact eligible trait-by-model grid")
    observed_full = lock.get("full_gwas_sha256_by_trait", {})
    if set(observed_full) != set(panel_ids):
        fail("TWAS lock does not bind full dense inputs for the exact 45-trait panel")
    for trait_id, expected in observed_full.items():
        source, _, ready = fine_mapping_contract.choose_full_input(root, trait_id)
        if not ready or source is None or sha256(source) != expected:
            fail(f"full dense TWAS input differs from lock: {trait_id}")
    mapping_by_id = {row["mapping_id"]: row for row in mappings}
    if any(
        row["mapping_id"] not in mapping_by_id
        or row["trait_id"] != mapping_by_id[row["mapping_id"]]["trait_id"]
        or row["model_family"] != mapping_by_id[row["mapping_id"]]["model_family"]
        or row["results_accessed_before_lock"] != "NO"
        for row in runs
    ) or any(row["results_accessed_before_lock"] != "NO" for row in mappings + eligibility):
        fail("TWAS run family contains an out-of-lock or post-result row")
    return mappings, eligibility, runs, lock
