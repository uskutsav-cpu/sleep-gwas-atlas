#!/usr/bin/env python3
"""Shared immutable-contract checks for atlas fine-mapping and trait coloc."""
from __future__ import annotations

import csv
from functools import lru_cache
import hashlib
import json
import re
from pathlib import Path

import lava_contract
import variant_map


SHARED_FIELDS = [
    "shared_locus_id", "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
    "locus_id", "CHR", "START", "STOP", "placo_lead_snp", "placo_lead_p",
    "conjfdr_lead_snp", "conjfdr_lead_fdr", "same_lead_snp", "effect_direction",
    "placo_significant_variant_count", "conjfdr_significant_variant_count",
    "evidence_status", "claim_limit",
]
MANIFEST_FIELDS = [
    "comparison_id", "shared_locus_id", "pair_id", "sleep_trait", "non_sleep_trait",
    "analysis_tier", "lava_block_id", "chromosome", "start_bp", "end_bp",
    "placo_lead_snp", "placo_lead_p", "conjfdr_lead_snp", "conjfdr_lead_fdr",
    "effect_direction", "sleep_full_input", "non_sleep_full_input", "reference_prefix",
    "summary1_path", "summary2_path", "variant_order_path", "ld_path", "task_path",
]
ENGINE_OUTPUT_FIELDS = {
    "variants.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role", "dataset_type", "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO", "prior_method", "normalized_prior_weight", "PIP", "credible_set_ids", "max_alpha_component", "model_converged", "rss_ld_s", "kriging_allele_switch_outlier"],
    "credible_sets.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role", "signal_id", "component_index", "lead_snp", "lead_pip", "credible_set_size", "credible_set_snps", "requested_coverage", "achieved_coverage", "min_abs_corr", "mean_abs_corr", "median_abs_corr", "cs_log10bf", "model_converged"],
    "colocalization.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset1_id", "dataset2_id", "molecular_feature_id", "tissue_cell_context", "coloc_method", "p1", "p2", "p12", "prior_role", "signal1", "signal2", "hit1", "hit2", "nsnps", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant", "top_shared_variant_PP_H4", "fine_mapping_qc", "analysis_status", "single_signal_fallback_justification", "claim_limit"],
    "shared_variant_posteriors.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "coloc_method", "p12", "signal1", "signal2", "SNP", "SNP_PP_H4"],
    "diagnostics.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role", "variant_count", "model_converged", "niter", "credible_set_count", "max_pip", "rss_ld_s", "kriging_allele_switch_outlier_count", "kriging_allele_switch_outliers", "diagnostic_status"],
}
CANONICAL_LOCUS_FIELDS = [
    "locus_id", "chromosome", "start_bp", "end_bp", "lead_snp", "sleep_trait",
    "non_sleep_trait", "analysis_tier", "placo_lead_p", "conjfdr_min_q",
    "effect_direction", "lava_block_id", "method_support", "evidence_level", "provenance_id",
]
CANONICAL_VARIANT_FIELDS = [
    "variant_id", "locus_id", "rsid", "chromosome", "position_bp", "effect_allele",
    "other_allele", "minor_allele_frequency", "sleep_beta", "sleep_se", "non_sleep_beta",
    "non_sleep_se", "pip_sleep", "pip_non_sleep", "credible_set_sleep",
    "credible_set_non_sleep", "shared_signal_posterior", "fine_mapping_method",
    "ld_reference", "qc_status", "provenance_id",
]
SCRIPT_PATHS = (
    "scripts/56_finemapping_preflight.py",
    "scripts/57_prepare_finemapping_loci.py",
    "scripts/58_materialize_finemapping_locus.py",
    "scripts/58_extract_lava_ld.R",
    "scripts/59_run_finemapping_locus.py",
    "scripts/60_collate_finemapping.py",
    "scripts/fine_mapping_contract.py",
    "scripts/lava_contract.py",
)
HEX64 = re.compile(r"[0-9a-f]{64}")


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


def script_hashes(root: Path) -> dict[str, str]:
    return {relative: sha256(root / relative) for relative in SCRIPT_PATHS}


def qc_value(path: Path, key: str) -> str:
    if not path.is_file():
        return ""
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0] == key:
            return fields[1]
    return ""


def qc_is_full_resolution(path: Path) -> bool:
    variant_map_strategy = qc_value(path, "variant_map_strategy")
    variant_map_scope = qc_value(path, "variant_map_scope")
    if qc_value(path, "prefilter_strategy") not in {"", "not supplied"}:
        return False
    if variant_map_strategy == "not supplied":
        return True
    if (
        variant_map_strategy not in {"BY_COORD_ALLELES", "BY_RSID_ALLELES"}
        or variant_map_scope != "GENOME_WIDE_IMPUTED_VARIANT_IDENTITY"
    ):
        return False
    map_path = qc_value(path, "variant_map")
    map_bytes = qc_value(path, "variant_map_bytes")
    map_sha256 = qc_value(path, "variant_map_sha256")
    provenance_sha256 = qc_value(path, "variant_map_provenance_sha256")
    schema_version = qc_value(path, "variant_map_schema_version")
    if (
        not map_bytes.isdigit() or not HEX64.fullmatch(map_sha256)
        or not HEX64.fullmatch(provenance_sha256)
    ):
        return False
    return sealed_genome_wide_map(
        map_path, int(map_bytes), map_sha256, provenance_sha256, schema_version,
    )


@lru_cache(maxsize=4)
def sealed_genome_wide_map(
    map_path_value: str,
    map_bytes: int,
    map_sha256: str,
    provenance_sha256: str,
    schema_version: str,
) -> bool:
    map_path = Path(map_path_value)
    provenance_path = Path(f"{map_path}.provenance.json")
    try:
        if (
            not map_path.is_absolute() or not map_path.is_file()
            or map_path.stat().st_size != map_bytes or not provenance_path.is_file()
            or sha256(provenance_path) != provenance_sha256
        ):
            return False
        provenance = variant_map.validate_provenance(
            map_path, expected_sha256=map_sha256, expected_bytes=map_bytes,
        )
    except (OSError, UnicodeError, ValueError, SystemExit):
        return False
    return (
        provenance.get("schema_version") == schema_version
        and provenance.get("map_scope") == "GENOME_WIDE_IMPUTED_VARIANT_IDENTITY"
    )


def choose_full_input(root: Path, trait: str) -> tuple[Path, Path, bool]:
    candidates = [
        (
            root / "data/harmonized_mixer_full" / f"{trait}.harmonized.tsv.gz",
            root / "data/harmonized_mixer_full" / f"{trait}.qc.txt",
        ),
        (
            root / "data/harmonized" / f"{trait}.harmonized.tsv.gz",
            root / "data/harmonized" / f"{trait}.qc.txt",
        ),
    ]
    observed: tuple[Path, Path] | None = None
    for harmonized, qc in candidates:
        if harmonized.is_file() and harmonized.stat().st_size and qc.is_file():
            if qc_is_full_resolution(qc):
                return harmonized, qc, True
            if observed is None:
                observed = (harmonized, qc)
    harmonized, qc = observed or candidates[0]
    return harmonized, qc, False


def load_policy(root: Path, relative: str = "config/fine_mapping_analysis_policy.json") -> tuple[Path, dict[str, object]]:
    path = root / relative
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"fine-mapping policy is unreadable: {exc}")
    if (
        policy.get("analysis_id") != "atlas-v1.0-fine-mapping-colocalization"
        or policy.get("analysis_panel") != "atlas-v1.0"
        or policy.get("entry_source") != "results/atlas/shared_loci.tsv"
        or policy.get("entry_provenance") != "results/atlas/shared_loci.provenance.json"
        or policy.get("lava_policy") != "config/lava_analysis_policy.json"
    ):
        fail("fine-mapping policy differs from the frozen atlas scope")
    return path, policy


def validate_shared_upstream(root: Path, policy: dict[str, object]) -> tuple[Path, Path, list[dict[str, str]], dict[str, object]]:
    shared_path = root / str(policy["entry_source"])
    provenance_path = root / str(policy["entry_provenance"])
    _, shared = read_tsv(shared_path, SHARED_FIELDS)
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"shared-locus provenance is unreadable: {exc}")
    pleio_policy_path = root / "config/pleiotropy_analysis_policy.json"
    pleio_policy = json.loads(pleio_policy_path.read_text(encoding="utf-8"))
    manifest_path = root / "results/tables/pleiotropy_pair_manifest.tsv"
    manifest_lock_path = root / "results/tables/pleiotropy_pair_manifest.lock.json"
    locus_path = root / str(pleio_policy["locus_definition"])
    if (
        provenance.get("schema_version") != "sleep-atlas-pleiotropic-loci.1"
        or provenance.get("completed_pair_scans") != 396
        or provenance.get("canonical_shared_locus_count") != len(shared)
        or provenance.get("canonical_output_sha256", {}).get(str(policy["entry_source"])) != sha256(shared_path)
        or provenance.get("policy_sha256") != sha256(pleio_policy_path)
        or provenance.get("manifest_sha256") != sha256(manifest_path)
        or provenance.get("manifest_lock_sha256") != sha256(manifest_lock_path)
        or provenance.get("locus_definition_sha256") != sha256(locus_path)
    ):
        fail("shared-locus table differs from its complete upstream provenance")
    runtime_path = root / str(pleio_policy["runtime_provenance"])
    if provenance.get("runtime_provenance_sha256") != sha256(runtime_path):
        fail("shared-locus runtime provenance drifted")
    upstream_scripts = (
        "scripts/43_prepare_pleiotropy_pairs.py", "scripts/44_materialize_pleiotropy_pair.py",
        "scripts/45_prepare_placo_task.py", "scripts/46_run_placo_pair.R",
        "scripts/47_prepare_pleiofdr_trait.py", "scripts/48_prepare_conjfdr_task.py",
        "scripts/49_run_conjfdr_pair.py", "scripts/50_collate_pleiotropy.py",
        "scripts/pleiotropy_contract.py",
    )
    if set(provenance.get("script_sha256", {})) != set(upstream_scripts):
        fail("shared-locus upstream script family is incomplete")
    for relative, expected in provenance["script_sha256"].items():
        if not HEX64.fullmatch(str(expected)) or sha256(root / safe_relative(relative, "upstream script")) != expected:
            fail(f"shared-locus upstream script drifted: {relative}")
    identities = [row["shared_locus_id"] for row in shared]
    if len(set(identities)) != len(identities):
        fail("duplicate shared-locus identifier")
    if any(row["evidence_status"] != "PLACO_PLUS_AND_CONJFDR_SAME_LOCKED_LD_BLOCK" for row in shared):
        fail("shared-locus table contains a non-consensus locus")
    return shared_path, provenance_path, shared, provenance


def validate_lava_reference(root: Path, policy: dict[str, object], *, rehash_payloads: bool) -> tuple[Path, dict[str, object]]:
    lava_policy_path, lava_policy = lava_contract.load_policy(root)
    reference = policy["reference"]
    if (
        str(reference["prefix"]) != lava_policy["reference_prefix"]
        or reference["id"] != "LAVA_UKB_v1.1_EUR_GRCh37"
        or reference["build"] != "GRCh37"
        or reference["ancestry"] != "EUR"
        or str(policy["lava_policy"]) != str(lava_policy_path.relative_to(root))
    ):
        fail("fine-mapping signed-LD policy differs from the pinned LAVA reference")
    provenance = lava_contract.validate_reference(root, lava_policy, rehash=False)
    if rehash_payloads:
        for record in provenance.get("extracted_files", []):
            path = root / safe_relative(str(record.get("path", "")), "LAVA reference payload")
            if sha256(path) != record.get("sha256"):
                fail(f"LAVA extracted reference SHA-256 mismatch: {path}")
    provenance_path = root / str(lava_policy["reference_provenance"])
    return provenance_path, provenance


def validate_preflight(root: Path, policy_path: Path, policy: dict[str, object], path: Path) -> dict[str, object]:
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"fine-mapping preflight is unreadable: {exc}")
    shared_path, provenance_path, shared, _ = validate_shared_upstream(root, policy)
    readiness_path = root / safe_relative(str(report.get("trait_readiness_path", "")), "trait readiness")
    _, readiness = read_tsv(readiness_path)
    primary_count = sum(row["analysis_tier"] == "PRIMARY_PHASE1" for row in shared)
    if (
        report.get("ready") is not True
        or report.get("policy_sha256") != sha256(policy_path)
        or report.get("shared_loci_sha256") != sha256(shared_path)
        or report.get("shared_loci_provenance_sha256") != sha256(provenance_path)
        or report.get("primary_shared_locus_count") != primary_count
        or report.get("analysis_required") != (primary_count > 0)
        or report.get("trait_readiness_sha256") != sha256(readiness_path)
        or len(readiness) != 45
    ):
        fail("fine-mapping preflight is absent, blocked, or provenance-drifted")
    if primary_count > 0:
        reference_path, _ = validate_lava_reference(root, policy, rehash_payloads=False)
        if (
            report.get("lava_reference_provenance_sha256") != sha256(reference_path)
            or report.get("lava_reference_content_verified") is not True
        ):
            fail("fine-mapping preflight does not bind a fully verified LAVA reference")
    return report


def validate_manifest(
    root: Path, manifest_path: Path, lock_path: Path, policy_path: Path,
    policy: dict[str, object], preflight_path: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    _, rows = read_tsv(manifest_path, MANIFEST_FIELDS)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"fine-mapping manifest lock is unreadable: {exc}")
    shared_path, shared_provenance_path, shared, _ = validate_shared_upstream(root, policy)
    primary = [row for row in shared if row["analysis_tier"] == "PRIMARY_PHASE1"]
    expected_ids = [f"{row['shared_locus_id']}__TRAIT_TRAIT" for row in primary]
    expected_contract_inputs = {
        relative: sha256(root / relative) for relative in (
            "config/analysis_panel.tsv", "config/downstream_analysis_policy.json",
            "config/fine_mapping_sources.tsv",
        )
    }
    if (
        lock.get("schema_version") != "atlas-v1.0-finemapping-locus-manifest.2"
        or lock.get("results_accessed_before_lock") is not False
        or lock.get("manifest_sha256") != sha256(manifest_path)
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("preflight_sha256") != sha256(preflight_path)
        or lock.get("shared_loci_sha256") != sha256(shared_path)
        or lock.get("shared_loci_provenance_sha256") != sha256(shared_provenance_path)
        or lock.get("locus_count") != len(rows)
        or lock.get("comparison_ids_in_locked_order") != expected_ids
        or [row["comparison_id"] for row in rows] != expected_ids
        or lock.get("selection_rule") != policy["entry_rule"]
        or lock.get("zero_family_rule") != policy["zero_family_rule"]
        or lock.get("claim_limit") != policy["claim_limit"]
        or lock.get("contract_input_sha256") != expected_contract_inputs
        or lock.get("script_sha256") != script_hashes(root)
    ):
        fail("fine-mapping locus manifest differs from its immutable family lock")
    for observed, source in zip(rows, primary):
        base = f"data/fine_mapping/{source['shared_locus_id']}"
        expected = {
            "comparison_id": f"{source['shared_locus_id']}__TRAIT_TRAIT",
            "shared_locus_id": source["shared_locus_id"],
            "pair_id": source["pair_id"],
            "sleep_trait": source["sleep_trait"],
            "non_sleep_trait": source["non_sleep_trait"],
            "analysis_tier": source["analysis_tier"],
            "lava_block_id": source["locus_id"],
            "chromosome": source["CHR"], "start_bp": source["START"], "end_bp": source["STOP"],
            "placo_lead_snp": source["placo_lead_snp"], "placo_lead_p": source["placo_lead_p"],
            "conjfdr_lead_snp": source["conjfdr_lead_snp"],
            "conjfdr_lead_fdr": source["conjfdr_lead_fdr"],
            "effect_direction": source["effect_direction"],
            "reference_prefix": str(policy["reference"]["prefix"]),
            "summary1_path": f"{base}/sleep.tsv.gz",
            "summary2_path": f"{base}/non_sleep.tsv.gz",
            "variant_order_path": f"{base}/variants.tsv",
            "ld_path": f"{base}/ld.tsv.gz",
            "task_path": f"results/fine_mapping/tasks/{source['shared_locus_id']}.tsv",
        }
        if any(observed.get(field) != value for field, value in expected.items()):
            fail(f"fine-mapping manifest row drifted from upstream locus: {source['shared_locus_id']}")
        for field in ("sleep_full_input", "non_sleep_full_input"):
            safe_relative(observed[field], field)
        for field, trait_field in (
            ("sleep_full_input", "sleep_trait"), ("non_sleep_full_input", "non_sleep_trait"),
        ):
            expected_input, _, ready = choose_full_input(root, observed[trait_field])
            if not ready or observed[field] != str(expected_input.relative_to(root)):
                fail(f"fine-mapping manifest does not use the current dense input: {observed[trait_field]}")
    return rows, lock
