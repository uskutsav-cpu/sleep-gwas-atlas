#!/usr/bin/env python3
"""Freeze the pre-result Track B fine-mapping and trait-coloc contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "config/track_b_finemapping_policy.json"
PAIR_MANIFEST = ROOT / "results/track_b/pair_manifest.tsv"
PAIR_LOCK = ROOT / "results/track_b/pair_manifest.lock.json"
DENSE_QC = ROOT / "results/track_b/03_dense_input_qc.tsv"
DENSE_LOCK = ROOT / "results/track_b/03_dense_input_qc.lock.json"
LOCAL_POLICY = ROOT / "config/track_b_local_analysis_policy.json"
LOCAL_INPUT_LOCK = ROOT / "results/track_b/local_analysis_input.lock.json"
SOFTWARE_SOURCES = ROOT / "config/fine_mapping_sources.tsv"
METHOD_REFERENCES = ROOT / "discovery_extension/config/fine_mapping_method_references.tsv"
SHARED_ENGINE = ROOT / "discovery_extension/scripts/35_run_susie_coloc.R"
SHARED_CONTRACT = ROOT / "scripts/fine_mapping_contract.py"
PLEIOTROPY_CONTRACT_LOCK = ROOT / "results/track_b/pleiotropy/contract.lock.json"
PLEIOTROPY_INPUT_GATE_LOCK = ROOT / "results/track_b/pleiotropy/input_gate.lock.json"
OUT_DIR = ROOT / "results/track_b/finemapping"
READINESS = OUT_DIR / "readiness.tsv"
LOCK = OUT_DIR / "contract.lock.json"

READINESS_FIELDS = [
    "gate_order", "gate_id", "status", "blocker_class", "observed_state",
    "required_to_unblock", "authoritative_evidence", "result_accessed", "claim_limit",
]
PAIR_SCOPE = [
    ("A", "snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    ("B", "insomnia", "adhd", "PRIMARY_DISCOVERY"),
    ("CONTROL", "insomnia", "frailty", "POSITIVE_CONTROL_NON_NOVELTY"),
]
DENSE_TRAITS = ["snoring", "parental_lifespan", "insomnia", "adhd", "frailty"]
FINE_MAPPING_LABELS = [
    "HIGH_PIP_VARIANT", "SMALL_CREDIBLE_SET", "DIFFUSE_SIGNAL", "LD_UNCERTAINTY",
    "MODEL_INSTABILITY",
]
COLOC_LABELS = [
    "STRONG_SHARED_SIGNAL", "MODERATE_SHARED_SIGNAL", "DISTINCT_SIGNALS",
    "INCONCLUSIVE", "INVALID_INPUT",
]
PLEIOTROPY_LABELS = ["PLACO_AND_CONJFDR", "PLACO_ONLY", "CONJFDR_ONLY"]
LOCAL_RESULT_SCHEMA = [
    "pair_id", "chr", "start", "end", "locus_id", "trait1", "trait2",
    "local_h2_trait1", "local_h2_trait2", "local_h2_scale", "local_covariance",
    "local_rg", "SE", "SE_method", "local_rg_CI_lower", "local_rg_CI_upper",
    "P", "FDR", "direction", "SNP_count", "QC", "interpretation_status",
]
PLEIOTROPY_UNION_SCHEMA_PREFIX = [
    "pair_id", "lead_variant", "chr", "position", "PLACO_P", "FDR", "trait1_P",
    "trait2_P", "locus_start", "locus_end", "independent_signal", "annotations",
]
RESOURCE_METRICS_SCHEMA = [
    "analysis_id", "run_fingerprint", "locus_entry_id", "pair_id", "attempt_id",
    "worker_environment_id", "process_isolation", "variant_count", "ld_dimension",
    "estimated_peak_rss_bytes", "host_physical_memory_bytes", "available_memory_before_bytes",
    "admission_status", "started_utc", "finished_utc", "elapsed_monotonic_seconds",
    "peak_rss_bytes", "rss_measurement_backend", "exit_code", "terminal_status",
    "task_sha256", "trait1_summary_sha256", "trait2_summary_sha256",
    "variant_order_sha256", "signed_ld_sha256", "engine_sha256", "policy_sha256",
    "output_manifest_sha256", "error",
]
EXPECTED_METHOD_REFERENCES = {
    "SUSIE_MODEL_2020", "SUSIE_RSS_2022", "COLOC_ORIGINAL_2014",
    "COLOC_PRIORS_2020", "COLOC_SUSIE_2021", "POLYFUN_2020", "FINEMAP_2016",
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {rel(path)}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty tabular artifact: {rel(path)}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"unreadable JSON artifact {rel(path)}: {error}")
    if not isinstance(value, dict):
        fail(f"JSON artifact must contain an object: {rel(path)}")
    return value


def load_policy() -> dict[str, Any]:
    policy = load_json(POLICY)
    if (
        policy.get("schema_version") != "sleep-atlas-track-b-finemapping-policy.1"
        or policy.get("analysis_id") != "track-b-v1.0-finemapping-trait-coloc"
        or policy.get("analysis_build") != "hg19"
        or policy.get("ancestry") != "EUR"
        or policy.get("selection_timing")
        != "BEFORE_TRACK_B_FINE_MAPPING_OR_TRAIT_COLOC_RESULT_ACCESS"
        or policy.get("results_accessed_before_contract_freeze") is not False
    ):
        fail("Track B fine-mapping policy identity or pre-result timing drifted")

    immutable = policy.get("immutable_pre_result_inputs", {})
    if immutable != {
        "pair_manifest": rel(PAIR_MANIFEST),
        "pair_manifest_lock": rel(PAIR_LOCK),
        "dense_qc": rel(DENSE_QC),
        "dense_qc_lock": rel(DENSE_LOCK),
        "local_analysis_policy": rel(LOCAL_POLICY),
        "local_analysis_input_lock": rel(LOCAL_INPUT_LOCK),
        "software_sources": rel(SOFTWARE_SOURCES),
        "method_references": rel(METHOD_REFERENCES),
        "shared_r_engine": rel(SHARED_ENGINE),
        "shared_python_contract": rel(SHARED_CONTRACT),
        "pleiotropy_contract_lock": rel(PLEIOTROPY_CONTRACT_LOCK),
        "pleiotropy_input_gate_lock": rel(PLEIOTROPY_INPUT_GATE_LOCK),
    }:
        fail("Track B fine-mapping immutable pre-result dependency map drifted")

    pair_scope = policy.get("pair_scope", {})
    observed_pairs = [
        (row.get("pair_id"), row.get("trait1"), row.get("trait2"), row.get("family_role"))
        for row in pair_scope.get("pairs", [])
    ]
    if (
        observed_pairs != PAIR_SCOPE
        or pair_scope.get("pair_order") != [row[0] for row in PAIR_SCOPE]
        or pair_scope.get("primary_pairs") != ["A", "B"]
        or pair_scope.get("control_pairs") != ["CONTROL"]
    ):
        fail("Track B fine-mapping pair roles drifted")

    locus = policy.get("locus_entry_contract", {})
    if (
        locus.get("intersection_required") is not False
        or locus.get("consensus_only_selection_forbidden") is not True
        or locus.get("result_ranked_maximum_forbidden") is not True
        or locus.get("global_rg_filter_forbidden") is not True
        or locus.get("zero_family_terminal_status")
        != "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY"
    ):
        fail("Track B locus-entry union or zero-family semantics drifted")

    dense = policy.get("dense_input_contract", {})
    if (
        dense.get("unique_traits") != DENSE_TRAITS
        or dense.get("minimum_variant_rows_per_trait") != 5_000_000
        or dense.get("hapmap3_only_forbidden") is not True
        or dense.get("minimum_locus_variants") != 50
        or dense.get("maximum_locus_variants") is not None
        or "no count cap" not in dense.get("locus_size_rule", "")
    ):
        fail("Track B full-resolution dense-input scope drifted")

    signed_ld = policy.get("signed_ld_contract", {})
    if (
        signed_ld.get("reference_id") != "LAVA_UKB_v1.1_EUR_GRCh37"
        or signed_ld.get("ancestry") != "EUR"
        or signed_ld.get("build") != "GRCh37/hg19"
        or signed_ld.get("exact_variant_order_sha256_required") is not True
        or signed_ld.get("reference_payload_sha256_required") is not True
        or signed_ld.get("minimum_reference_sample_size") != 10000
    ):
        fail("Track B signed-LD contract drifted")

    fine_mapping = policy.get("fine_mapping_contract", {})
    coloc = policy.get("colocalization_contract", {})
    if (
        fine_mapping.get("primary_method") != "SuSiE-RSS"
        or fine_mapping.get("maximum_causal_signals") != 10
        or fine_mapping.get("credible_set_coverage") != 0.95
        or fine_mapping.get("primary_prior") != "FLAT"
        or coloc.get("primary_method") != "coloc.susie"
        or coloc.get("p1") != 0.0001
        or coloc.get("p2") != 0.0001
        or coloc.get("p12_primary") != 0.00001
        or coloc.get("p12_sensitivity_grid") != [0.000001, 0.000005, 0.00001, 0.00005]
        or "same signal1/signal2 pair" not in coloc.get("prior_robust_rule", "")
    ):
        fail("Track B SuSiE/coloc prior or sensitivity contract drifted")
    software = policy.get("software_contract", {})
    if "generic atlas consensus-only" not in software.get("shared_python_contract_use_rule", ""):
        fail("generic atlas consensus-only fine-mapping selection is not explicitly forbidden")

    classes = policy.get("classification_contract", {})
    if (
        classes.get("fine_mapping_labels_in_order") != FINE_MAPPING_LABELS
        or classes.get("trait_coloc_allowed") != COLOC_LABELS
    ):
        fail("Track B fine-mapping or trait-coloc classifications drifted")

    future = policy.get("future_upstream_results", {})
    pleiotropy = future.get("pleiotropy", {})
    if (
        future.get("local_sharing", {}).get("canonical_result")
        != "results/track_b/04_lava_local_results.tsv"
        or future.get("local_sharing", {}).get("required_schema") != LOCAL_RESULT_SCHEMA
        or pleiotropy.get("placo_variants") != "results/track_b/07_placo_plus_variants.tsv"
        or pleiotropy.get("method_union_loci") != "results/track_b/08_shared_loci.tsv"
        or pleiotropy.get("method_comparison") != "results/track_b/09_pleiotropy_comparison.tsv"
        or pleiotropy.get("required_method_labels") != PLEIOTROPY_LABELS
        or pleiotropy.get("method_union_required_schema_prefix")
        != PLEIOTROPY_UNION_SCHEMA_PREFIX
    ):
        fail("Track B future local/pleiotropy dependency paths drifted")

    outputs = policy.get("required_science_outputs", {})
    if [outputs.get(key, {}).get("path") for key in ("output_10", "output_11", "output_12")] != [
        "results/track_b/10_finemap_trait1.tsv",
        "results/track_b/11_finemap_trait2.tsv",
        "results/track_b/12_trait_trait_coloc.tsv",
    ]:
        fail("Track B canonical 10/11/12 paths drifted")
    if (
        outputs.get("output_10", {}).get("schema")
        != outputs.get("output_11", {}).get("schema")
        or outputs.get("output_10", {}).get("trait_role") != "TRAIT1_SLEEP"
        or outputs.get("output_11", {}).get("trait_role") != "TRAIT2_EXTERNAL"
        or outputs.get("output_12", {}).get("schema", [None])[-3:]
        != ["classification", "analysis_status", "claim_limit"]
    ):
        fail("Track B exact 10/11/12 schemas or trait roles drifted")

    execution = policy.get("execution_and_publication", {})
    if not all(execution.get(key) is True for key in (
        "pre_result_locus_manifest_lock_required",
        "results_access_before_manifest_lock_forbidden",
        "complete_locus_family_required",
        "failed_loci_retained",
        "parameter_changes_after_result_access_forbidden",
        "canonical_output_overwrite_forbidden",
        "staged_validation_required",
        "exclusive_publication_required",
        "result_provenance_required",
        "generic_atlas_science_output_substitution_forbidden",
    )):
        fail("Track B no-result-access or no-overwrite publication contract drifted")
    statuses = policy.get("result_status_semantics", {})
    if (
        statuses.get("fine_mapping_allowed") != ["COMPLETE", "FAILED_QC"]
        or statuses.get("trait_coloc_allowed") != ["TESTED", "INVALID_INPUT"]
        or "forbidden inside 10/11/12" not in statuses.get("blocked_rule", "")
    ):
        fail("Track B complete-family result status semantics drifted")
    if policy.get("current_readiness", {}).get("status") != "BLOCKED_UPSTREAM":
        fail("pre-result Track B fine-mapping readiness must remain BLOCKED_UPSTREAM")

    resources = policy.get("ram_aware_execution_contract", {})
    estimate = resources.get("admission_estimate", {})
    atomic_resume = resources.get("atomic_resume", {})
    if (
        resources.get("execution_unit") != "ONE_VALIDATED_PAIR_LOCUS_PER_FRESH_OS_PROCESS"
        or resources.get("maximum_concurrent_loci_per_worker") != 1
        or resources.get("process_reuse_across_loci_forbidden") is not True
        or resources.get("full_locus_variant_universe_required") is not True
        or resources.get("full_ancestry_matched_signed_ld_required") is not True
        or resources.get("locus_splitting_forbidden") is not True
        or resources.get("variant_thinning_for_compute_forbidden") is not True
        or estimate.get("ld_scalar_bytes") != 8
        or estimate.get("ld_live_copy_multiplier") != 8
        or estimate.get("fixed_process_overhead_bytes") != 2 * 1024**3
        or estimate.get("absolute_minimum_memory_bytes") != 4 * 1024**3
        or estimate.get("minimum_unallocated_host_reserve_bytes") != 1024**3
        or resources.get("resource_metrics_schema") != RESOURCE_METRICS_SCHEMA
        or "exclusive no-replace" not in atomic_resume.get("staging_rule", "")
        or "never overwritten" not in atomic_resume.get("resume_rule", "")
    ):
        fail("Track B per-locus RAM, measurement, or atomic-resume contract drifted")
    return policy


def validate_pair_and_dense_inputs(policy: dict[str, Any]) -> dict[str, str]:
    _, pairs = read_tsv(PAIR_MANIFEST)
    if [
        (row["pair_id"], row["sleep_trait"], row["external_trait"])
        for row in pairs
    ] != [row[:3] for row in PAIR_SCOPE]:
        fail("live Track B pair manifest differs from the frozen fine-mapping family")
    pair_lock = load_json(PAIR_LOCK)
    if (
        pair_lock.get("output_sha256", {}).get(rel(PAIR_MANIFEST)) != sha256(PAIR_MANIFEST)
        or pair_lock.get("downstream_results_accessed_before_pair_freeze") is not False
        or pair_lock.get("pair_replacement_policy", "").split(";", 1)[0] != "FORBIDDEN"
    ):
        fail("Track B pair manifest differs from its pre-result lock")

    _, dense_rows = read_tsv(DENSE_QC)
    dense_lock = load_json(DENSE_LOCK)
    dense_hash = sha256(DENSE_QC)
    if (
        dense_lock.get("dense_qc_sha256") != dense_hash
        or dense_lock.get("pair_manifest_sha256") != sha256(PAIR_MANIFEST)
        or dense_lock.get("all_are_dense_not_hapmap3_only") is not True
        or dense_lock.get("unique_dense_trait_family") != DENSE_TRAITS
        or dense_lock.get("trait_count") != len(DENSE_TRAITS)
    ):
        fail("Track B dense QC differs from its frozen lock")
    if [row["trait_id"] for row in dense_rows] != DENSE_TRAITS:
        fail("Track B dense trait order or uniqueness drifted")

    required_schema = ",".join(policy["dense_input_contract"]["required_columns"])
    live_dense_hashes: dict[str, str] = {}
    for row in dense_rows:
        trait = row["trait_id"]
        path = ROOT / row["dense_file"]
        try:
            expected_bytes = int(row["compressed_size_bytes"])
            rows_out = int(row["rows_out"])
        except ValueError:
            fail(f"invalid dense input counts for {trait}")
        if (
            row["analysis_build"] != "hg19"
            or row["ancestry"] != "EUR"
            or row["observed_schema"] != required_schema
            or rows_out < policy["dense_input_contract"]["minimum_variant_rows_per_trait"]
            or not row["fine_mapping_readiness"].startswith("PASS_DENSE")
            or not row["overall_QC"].startswith("PASS_DENSE")
            or not path.is_file()
            or path.stat().st_size != expected_bytes
        ):
            fail(f"full non-HapMap3 dense-input gate failed for {trait}")
        observed = sha256(path)
        if observed != row["dense_file_sha256"]:
            fail(f"live dense input differs from its QC hash for {trait}")
        live_dense_hashes[rel(path)] = observed
    return live_dense_hashes


def validate_software(policy: dict[str, Any]) -> None:
    _, sources = read_tsv(SOFTWARE_SOURCES)
    by_method = {row["method_id"]: row for row in sources}
    if set(by_method) != {"SUSIE_RSS_PRIMARY", "COLOC_SUSIE_PRIMARY"} or len(sources) != 2:
        fail("fine-mapping source registry is not the exact pinned two-method family")

    software = policy["software_contract"]
    expected = {
        "SUSIE_RSS_PRIMARY": ("susieR", software["susieR"]),
        "COLOC_SUSIE_PRIMARY": ("coloc", software["coloc"]),
    }
    for method_id, (name, specification) in expected.items():
        row = by_method[method_id]
        archive = ROOT / specification["source_archive"]
        description = ROOT / specification["installed_description"]
        if (
            row["software"] != name
            or row["software_version"] != specification["version"]
            or row["source_archive_sha256"] != specification["source_archive_sha256"]
            or row["installed_source_archive_path"] != specification["source_archive"]
            or row["installed_description_path"] != specification["installed_description"]
            or sha256(archive) != specification["source_archive_sha256"]
            or not description.is_file()
            or f"Version: {specification['version']}" not in description.read_text(
                encoding="utf-8", errors="strict",
            )
        ):
            fail(f"pinned {name} source/runtime identity drifted")

    _, references = read_tsv(METHOD_REFERENCES)
    if (
        len(references) != len(EXPECTED_METHOD_REFERENCES)
        or {row["reference_id"] for row in references} != EXPECTED_METHOD_REFERENCES
    ):
        fail("fine-mapping method reference registry is incomplete or duplicated")
    if sha256(SHARED_ENGINE) != software["shared_r_engine_sha256"]:
        fail("shared SuSiE/coloc R engine differs from its policy pin")
    sha256(SHARED_CONTRACT)

    rscript = ROOT / ".r-env/bin/Rscript"
    if not rscript.is_file():
        fail("pinned fine-mapping R runtime is absent")
    expression = (
        "suppressPackageStartupMessages(library(susieR));"
        "suppressPackageStartupMessages(library(coloc));"
        "ok<-as.character(packageVersion('susieR'))=='0.14.2' && "
        "as.character(packageVersion('coloc'))=='5.2.3' && "
        "all(vapply(c('susie_rss','susie_get_cs','estimate_s_rss','kriging_rss'),"
        "exists,logical(1),where=asNamespace('susieR'),inherits=FALSE)) && "
        "all(vapply(c('coloc.susie','runsusie','sensitivity'),exists,logical(1),"
        "where=asNamespace('coloc'),inherits=FALSE));cat(ok)"
    )
    smoke = subprocess.run(
        [str(rscript), "-e", expression], cwd=ROOT, capture_output=True, text=True, check=False,
    )
    if smoke.returncode != 0 or smoke.stdout.strip() != "TRUE":
        fail(f"pinned fine-mapping R runtime smoke test failed: {(smoke.stdout + smoke.stderr).strip()}")


def validate_local_pre_result_contract() -> None:
    local_policy = load_json(LOCAL_POLICY)
    local_lock = load_json(LOCAL_INPUT_LOCK)
    ram = local_policy.get("ram_aware_execution", {})
    if (
        local_policy.get("analysis_id") != "track-b-v1.0-local"
        or local_policy.get("lava_version") != "0.1.5"
        or local_policy.get("expected_loci") != 2495
        or local_policy.get("expected_analysis_traits") != 8
        or local_policy.get("expected_discovery_pairs") != 3
        or local_policy.get("planned_bivariate_pair_locus_family_max") != 2495 * 3
        or ram.get("execution_unit") != "WHOLE_PREDECLARED_LAVA_LOCUS"
        or ram.get("maximum_loci_per_process") != 1
        or ram.get("worker_process_rule")
        != "FRESH_R_PROCESS_PER_LOCUS_WITH_EXPLICIT_OBJECT_REMOVAL_AND_FULL_GARBAGE_COLLECTION"
    ):
        fail("settled Track B local policy is not the exact per-locus family contract")
    if (
        local_lock.get("analysis_id") != "track-b-v1.0-local"
        or local_lock.get("policy_sha256") != sha256(LOCAL_POLICY)
        or local_lock.get("pair_manifest_sha256") != sha256(PAIR_MANIFEST)
        or local_lock.get("local_results_accessed_before_input_freeze") is not False
        or local_lock.get("result_substitution_policy")
        != "FORBIDDEN_SYNTHETIC_OR_PARTIAL_OUTPUTS_CANNOT_BE_PROMOTED_TO_REAL_LOCAL_RESULTS"
    ):
        fail("settled Track B local input lock does not bind the pre-result local family")


def validate_pleiotropy_pre_result_locks() -> None:
    contract = load_json(PLEIOTROPY_CONTRACT_LOCK)
    input_gate = load_json(PLEIOTROPY_INPUT_GATE_LOCK)
    contract_sha256 = sha256(PLEIOTROPY_CONTRACT_LOCK)
    if (
        contract.get("schema_version") != "sleep-atlas-track-b-pleiotropy-contract.1"
        or contract.get("analysis_id") != "track-b-v1.0-pleiotropy"
        or contract.get("method_union_not_intersection") is not True
        or contract.get("pleiotropy_results_accessed_before_contract_freeze") is not False
        or contract.get("upstream", {}).get("pair_manifest_lock_sha256") != sha256(PAIR_LOCK)
        or contract.get("upstream", {}).get("dense_qc_lock_sha256") != sha256(DENSE_LOCK)
    ):
        fail("settled Track B pleiotropy contract lock is not an immutable method-union contract")
    if (
        input_gate.get("schema_version") != "sleep-atlas-track-b-pleiotropy-input-gate.1"
        or input_gate.get("analysis_id") != "track-b-v1.0-pleiotropy"
        or input_gate.get("contract_lock_sha256") != contract_sha256
        or input_gate.get("all_pair_dense_input_gates_pass") is not True
        or input_gate.get("pleiotropy_results_accessed_before_input_gate_freeze") is not False
        or input_gate.get("scientific_result_count") != 0
    ):
        fail("settled Track B pleiotropy input-gate lock is not a valid pre-result dependency")


def table_text(rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=READINESS_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def build() -> dict[Path, str]:
    policy = load_policy()
    dense_inputs = validate_pair_and_dense_inputs(policy)
    validate_software(policy)
    validate_local_pre_result_contract()
    validate_pleiotropy_pre_result_locks()
    future = policy["future_upstream_results"]
    outputs = policy["required_science_outputs"]
    claims = policy["claim_limits"]
    readiness_rows = [
        {
            "gate_order": "1", "gate_id": "FROZEN_PAIR_SCOPE", "status": "PASS",
            "blocker_class": "NONE", "observed_state": "A_AND_B_PRIMARY;CONTROL_NON_NOVELTY",
            "required_to_unblock": "NONE",
            "authoritative_evidence": f"{rel(PAIR_MANIFEST)};{rel(PAIR_LOCK)}",
            "result_accessed": "FALSE", "claim_limit": claims["blocked"],
        },
        {
            "gate_order": "2", "gate_id": "FULL_DENSE_INPUTS", "status": "PASS",
            "blocker_class": "NONE",
            "observed_state": "5/5_FULL_AUTOSOMAL_POST_QC_NON_HAPMAP3_INPUTS_CHECKSUM_VERIFIED",
            "required_to_unblock": "NONE",
            "authoritative_evidence": f"{rel(DENSE_QC)};{rel(DENSE_LOCK)}",
            "result_accessed": "FALSE", "claim_limit": claims["fine_mapping"],
        },
        {
            "gate_order": "3", "gate_id": "PINNED_METHODS_AND_ENGINE", "status": "PASS",
            "blocker_class": "NONE",
            "observed_state": "susieR_0.14.2;coloc_5.2.3;ENTRYPOINT_SMOKE_PASS;SHARED_ENGINE_HASH_PASS",
            "required_to_unblock": "NONE",
            "authoritative_evidence": (
                f"{rel(SOFTWARE_SOURCES)};{rel(METHOD_REFERENCES)};"
                f"{rel(SHARED_ENGINE)};{rel(SHARED_CONTRACT)}"
            ),
            "result_accessed": "FALSE", "claim_limit": claims["fine_mapping"],
        },
        {
            "gate_order": "4", "gate_id": "VALIDATED_LOCAL_SHARING_FAMILY",
            "status": "BLOCKED_UPSTREAM", "blocker_class": "UPSTREAM_LOCAL_RESULT",
            "observed_state": "FUTURE_COMPLETE_SEALED_LOCAL_FAMILY_REQUIRED",
            "required_to_unblock": (
                f"{future['local_sharing']['canonical_result']};"
                f"{future['local_sharing']['result_provenance']}"
            ),
            "authoritative_evidence": (
                f"{rel(LOCAL_POLICY)};{rel(LOCAL_INPUT_LOCK)};{rel(POLICY)}"
            ), "result_accessed": "FALSE",
            "claim_limit": claims["blocked"],
        },
        {
            "gate_order": "5", "gate_id": "VALIDATED_PLEIOTROPY_METHOD_UNION",
            "status": "BLOCKED_UPSTREAM", "blocker_class": "UPSTREAM_PLEIOTROPY_RESULT",
            "observed_state": "FUTURE_COMPLETE_PLACO_OR_CONJFDR_METHOD_UNION_REQUIRED;INTERSECTION_NOT_REQUIRED",
            "required_to_unblock": (
                f"{future['pleiotropy']['method_union_loci']};"
                f"{future['pleiotropy']['method_comparison']};"
                f"{future['pleiotropy']['result_provenance']}"
            ),
            "authoritative_evidence": (
                f"{future['pleiotropy']['contract_lock']};"
                f"{future['pleiotropy']['input_gate_lock']}"
            ),
            "result_accessed": "FALSE", "claim_limit": claims["blocked"],
        },
        {
            "gate_order": "6", "gate_id": "NONDUPLICATED_LOCUS_ENTRY_UNION",
            "status": "BLOCKED_UPSTREAM", "blocker_class": "UPSTREAM_LOCUS_FAMILIES",
            "observed_state": "NO_LOCUS_MANIFEST_CREATED;NO_CONSENSUS_ONLY_OR_RESULT_RANKED_SELECTION",
            "required_to_unblock": (
                f"{policy['locus_entry_contract']['future_locus_manifest']};"
                f"{policy['locus_entry_contract']['future_locus_manifest_lock']}"
            ),
            "authoritative_evidence": rel(POLICY), "result_accessed": "FALSE",
            "claim_limit": claims["blocked"],
        },
        {
            "gate_order": "7", "gate_id": "ANCESTRY_MATCHED_SIGNED_LD",
            "status": "DEFERRED_UPSTREAM", "blocker_class": "UPSTREAM_LOCUS_FAMILY",
            "observed_state": "NO_LOCUS_SPECIFIC_SIGNED_LD_REQUESTED_BEFORE_LOCUS_UNION_FREEZE",
            "required_to_unblock": (
                f"{policy['signed_ld_contract']['reference_provenance']};"
                "LOCUS_SPECIFIC_EXACT_VARIANT_ORDER_AND_ALLELE_LOCK"
            ),
            "authoritative_evidence": rel(POLICY), "result_accessed": "FALSE",
            "claim_limit": claims["fine_mapping"],
        },
        {
            "gate_order": "8", "gate_id": "FRESH_PROCESS_RESOURCE_EXECUTION",
            "status": "DEFERRED_UPSTREAM", "blocker_class": "UPSTREAM_LOCUS_FAMILY",
            "observed_state": "ONE_FULL_PAIR_LOCUS_PER_FRESH_PROCESS;NO_SPLITTING;NO_ATTEMPTS_STARTED",
            "required_to_unblock": (
                f"{policy['ram_aware_execution_contract']['future_resource_metrics']};"
                f"{policy['ram_aware_execution_contract']['atomic_resume']['future_run_index']}"
            ),
            "authoritative_evidence": rel(POLICY), "result_accessed": "FALSE",
            "claim_limit": claims["blocked"],
        },
        {
            "gate_order": "9", "gate_id": "CANONICAL_10_11_12_PUBLICATION",
            "status": "NOT_STARTED", "blocker_class": "UPSTREAM_LOCUS_AND_LD",
            "observed_state": "NO_10_11_12_SCIENCE_OUTPUT_OR_HEADER_ONLY_PLACEHOLDER_CREATED",
            "required_to_unblock": ";".join(outputs[key]["path"] for key in (
                "output_10", "output_11", "output_12",
            )),
            "authoritative_evidence": rel(POLICY), "result_accessed": "FALSE",
            "claim_limit": claims["trait_coloc"],
        },
        {
            "gate_order": "10", "gate_id": "OVERALL", "status": "BLOCKED_UPSTREAM",
            "blocker_class": "UPSTREAM_LOCAL_AND_PLEIOTROPY_RESULTS",
            "observed_state": "CONTRACT_READY;SCIENCE_NOT_RUN;ZERO_FAMILY_NOT_YET_EVALUABLE",
            "required_to_unblock": "COMPLETE_AND_SEAL_BOTH_UPSTREAM_FAMILIES_THEN_FREEZE_THEIR_NONDUPLICATED_UNION",
            "authoritative_evidence": f"{rel(POLICY)};{rel(LOCK)}",
            "result_accessed": "FALSE", "claim_limit": claims["blocked"],
        },
    ]
    readiness_text = table_text(readiness_rows)

    immutable_paths = [
        POLICY, PAIR_MANIFEST, PAIR_LOCK, DENSE_QC, DENSE_LOCK, LOCAL_POLICY,
        LOCAL_INPUT_LOCK, SOFTWARE_SOURCES,
        METHOD_REFERENCES, SHARED_ENGINE, SHARED_CONTRACT, PLEIOTROPY_CONTRACT_LOCK,
        PLEIOTROPY_INPUT_GATE_LOCK,
        ROOT / policy["software_contract"]["susieR"]["source_archive"],
        ROOT / policy["software_contract"]["susieR"]["installed_description"],
        ROOT / policy["software_contract"]["coloc"]["source_archive"],
        ROOT / policy["software_contract"]["coloc"]["installed_description"],
    ]
    schema_sha256 = {
        outputs[key]["path"]: hashlib.sha256(
            json.dumps(outputs[key]["schema"], separators=(",", ":")).encode()
        ).hexdigest()
        for key in ("output_10", "output_11", "output_12")
    }
    immutable_relatives = {rel(path) for path in immutable_paths}
    future_required_paths = sorted({
        future["local_sharing"]["canonical_result"],
        future["local_sharing"]["result_provenance"],
        future["pleiotropy"]["contract_lock"],
        future["pleiotropy"]["input_gate_lock"],
        future["pleiotropy"]["placo_variants"],
        future["pleiotropy"]["method_union_loci"],
        future["pleiotropy"]["method_comparison"],
        future["pleiotropy"]["result_provenance"],
        future["pleiotropy"]["control_directory"],
        policy["signed_ld_contract"]["reference_provenance"],
        policy["signed_ld_contract"]["reference_manifest"],
        policy["locus_entry_contract"]["future_locus_manifest"],
        policy["locus_entry_contract"]["future_locus_manifest_lock"],
        policy["ram_aware_execution_contract"]["future_resource_metrics"],
        policy["ram_aware_execution_contract"]["atomic_resume"]["future_run_index"],
    } - immutable_relatives)
    lock = {
        "schema_version": "sleep-atlas-track-b-finemapping-contract-lock.1",
        "analysis_id": policy["analysis_id"],
        "artifact_role": "PRE_RESULT_FINE_MAPPING_AND_TRAIT_COLOC_CONTRACT_NOT_SCIENTIFIC_RESULT",
        "status": "BLOCKED_UPSTREAM",
        "policy_sha256": sha256(POLICY),
        "script_sha256": sha256(Path(__file__)),
        "immutable_input_sha256": {rel(path): sha256(path) for path in immutable_paths},
        "dense_input_sha256": dense_inputs,
        "pair_order": [row[0] for row in PAIR_SCOPE],
        "primary_pairs": ["A", "B"],
        "control_pairs": ["CONTROL"],
        "locus_family_rule": "NONDUPLICATED_UNION_OF_VALID_LOCAL_OR_COMPLETED_PLACO_OR_COMPLETED_CONJFDR_SIGNALS",
        "consensus_only_selection": "FORBIDDEN",
        "result_ranked_maximum": "FORBIDDEN",
        "future_required_paths_not_hashed": future_required_paths,
        "future_science_outputs": [outputs[key]["path"] for key in (
            "output_10", "output_11", "output_12",
        )],
        "future_science_output_schema_sha256": schema_sha256,
        "readiness_sha256": hashlib.sha256(readiness_text.encode()).hexdigest(),
        "results_accessed_before_contract_freeze": False,
        "science_outputs_created_by_contract": False,
        "zero_family_terminal_status": policy["locus_entry_contract"]["zero_family_terminal_status"],
        "zero_family_creates_10_11_12": False,
        "canonical_output_overwrite": "FORBIDDEN",
        "execution_unit": policy["ram_aware_execution_contract"]["execution_unit"],
        "maximum_concurrent_loci_per_worker": 1,
        "full_locus_and_signed_ld_required": True,
        "locus_splitting": "FORBIDDEN",
        "resource_metrics_schema_sha256": hashlib.sha256(
            json.dumps(RESOURCE_METRICS_SCHEMA, separators=(",", ":")).encode()
        ).hexdigest(),
        "atomic_resume": "FINGERPRINTED_STAGE_VALIDATE_FSYNC_EXCLUSIVE_NO_REPLACE",
        "claim_limit": claims["blocked"],
    }
    return {
        READINESS: readiness_text,
        LOCK: json.dumps(lock, indent=2, sort_keys=True) + "\n",
    }


def publish_exclusive(payloads: dict[Path, str]) -> None:
    existing = [path for path in payloads if path.exists()]
    if existing:
        fail(
            "pre-result Track B fine-mapping contract already exists; overwrite is forbidden: "
            + ",".join(rel(path) for path in existing)
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    try:
        for path, text in payloads.items():
            with path.open("x", encoding="utf-8") as handle:
                handle.write(text)
            created.append(path)
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build()
    if args.verify:
        for path, text in expected.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                fail(f"Track B fine-mapping contract drifted: {rel(path)}")
        observed_lock = load_json(LOCK)
        if observed_lock.get("readiness_sha256") != sha256(READINESS):
            fail("Track B fine-mapping readiness hash drifted")
        print("TRACK_B_FINEMAPPING_CONTRACT_VERIFIED status=BLOCKED_UPSTREAM science_outputs=0")
        return

    policy = load_policy()
    science_outputs = [
        ROOT / policy["required_science_outputs"][key]["path"]
        for key in ("output_10", "output_11", "output_12")
    ]
    if any(path.exists() for path in science_outputs):
        fail("10/11/12 result access or publication preceded the required pre-result contract")
    publish_exclusive(expected)
    print("TRACK_B_FINEMAPPING_CONTRACT_WRITTEN status=BLOCKED_UPSTREAM science_outputs=0")


if __name__ == "__main__":
    main()
