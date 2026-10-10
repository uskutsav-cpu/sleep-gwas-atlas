#!/usr/bin/env python3
"""Seal a prepare-only LD-source/generator plan; contains no fit/acquisition path."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "sleep_unified_research_v4"
PLAN = V4 / "manifests/canonical_realistic_ld_source_admission_plan_v4.json"
GATES = V4 / "statistical_validation/canonical_realistic_ld_feasibility_gates_v4.tsv"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    if PLAN.exists() or GATES.exists():
        raise SystemExit("preserve existing files; choose a distinct version")
    probe_path = V4 / "statistical_validation/canonical_realistic_ld_feasibility_probe_receipt_v4.json"
    probe = json.loads(probe_path.read_text())
    fixture_path = V4 / "statistical_validation/canonical_signed_genotype_generator_fixture_receipt_v4.json"
    fixture = json.loads(fixture_path.read_text())
    blocks_path = V4 / "manifests/canonical_200_interval_preparation_v4.json"
    blocks = json.loads(blocks_path.read_text())
    native_path = V4 / "manifests/ssd_native_execution_plan_v4_3.json"
    native = json.loads(native_path.read_text())
    mfiles = {p: h for p, h in native["dependencies_sha256"].items() if p.endswith(".l2.M_5_50")}
    assert len(mfiles) == 22
    assert all(sha(p) == h for p, h in mfiles.items())
    M = sum(float(Path(p).read_text().strip()) for p in mfiles)
    assert M == 1173569.0
    counts = Counter()
    for r in probe["raw_genotype_per_chr_compatibility"]:
        for k, v in r.items():
            if k not in {"CHR", "raw_variants"}:
                counts[k] += v
    zenodo_path = V4 / "statistical_validation/canonical_ld_feasibility_sources_v4/zenodo_8292725_metadata.json"
    zenodo = json.loads(zenodo_path.read_text())
    archive = next(f for f in zenodo["files"] if f["key"] == "1000G_Phase3_plinkfiles.tgz")
    assert archive["size"] == 288277344 and archive["checksum"] == "md5:a7773ab485827b533cb300c76356d76b"
    gates = [
        ("L01", "consumed reference and coordinate identities", "PASS_SOURCE_HASHES_ONLY", "probe before/after; neither actual fit nor scientific calibration"),
        ("L02", "503 EUR FAM population/sample metadata", "PASS_METADATA_ONLY", "exact authoritative panel sample set; BED genotype content unauthenticated"),
        ("L03", "fresh BED identity, upstream derivation and dosage QC", "PENDING", "prior expected hash, size/header and BIM/FAM do not prove genotype payload"),
        ("L04", "actual surviving native SNP source coverage", "PENDING_NATIVE_CAPTURE", "reference-ID missingness does not establish any pair's surviving coverage"),
        ("L05", "restricted source-only SNP universe and orientation seal", "PROPOSED_NOT_FROZEN", "complete exact intersection recipe; no favorable GWAS-mask/outcome selection"),
        ("L06", "actual-source decoder, signed allele and operator verification", "PENDING", "independent dense-window/factor/noise controls under same source and double precision"),
        ("L07", "small genetic/noise covariance PSD and factor algebra", "PASS_INVENTED_FIXTURE_ONLY", "41 checks; no assets, random replicates or native fits"),
        ("L08", "unchanged historical score/weight/M science", "PINNED_NOT_CALIBRATED", "retain M1173569 and all original input hashes; generator/reference mismatch explicit"),
        ("L09", "integrated-vs-realized-effect truth and sample-overlap law", "PROPOSED_REVIEW_REQUIRED", "Gaussian source-factor definition; physical-cohort interpretation not established"),
        ("L10", "100/400/shift partitions and generalized native adapter", "PENDING", "source-only proposed rules; sealed control adapter still requires200"),
        ("L11", "historical completion and native 16 controls", "PENDING_SEPARATE_EXECUTOR", "v4_3 operational preparation passed; native admission remainsfalse"),
        ("L12", "one-worker total RSS/time/storage pilot", "PENDING_ROOT_ADMISSION", "2GiB, BLAS1, shared heavy lock, unchanged3/5GiB floors; no quota reduction"),
        ("L13", "full realistic signed LD and causal/regression universe", "NOT_ADMITTED", "complete provenance/coverage/population-LD and score/M match unresolved"),
        ("L14", "8000 full primary calibrations and sensitivities", "NOT_RUN", "S qualifications cannot count as full C primary replicates"),
        ("L15", "marginal/covariance/contrast coverage and MC intervals", "NOT_RUN", "all frozen protocol thresholds and failed-replicate accounting retained"),
        ("L16", "independent final scientific release decision", "CLOSED", "no 41Cov outcomes/calibrated P values admitted"),
    ]
    with GATES.open("x", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["gate_id", "requirement", "current_status", "qualification"])
        writer.writerows(gates)
    plan = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "schema": "canonical-realistic-ld-source-plan-v4.1",
        "prepare_only": True, "any_acquisition_or_estimator_path": False,
        "admitted": False, "restricted_S_qualification_admitted": False,
        "full_C_calibration_admitted": False, "biological_41_covariance_admitted": False,
        "preserved_protocol": {"path": str(V4 / "statistical_validation/COMMON_BOUNDARY_METHOD_VALIDATION_PROTOCOL_v4.md"),
                               "sha256": sha(V4 / "statistical_validation/COMMON_BOUNDARY_METHOD_VALIDATION_PROTOCOL_v4.md"),
                               "scenarios": 4, "replicates_each": 2000, "total": 8000,
                               "base_seed": 20261009, "bootstrap_replicates": 5000, "bootstrap_seed": 202610095,
                               "unchanged_thresholds_and_failed_replicate_accounting": True},
        "bound_small_artifacts": {str(p): sha(p) for p in [probe_path, fixture_path, blocks_path, native_path, zenodo_path,
                                    V4 / "reviews/genomicsem_realistic_ld_calibration_feasibility_v4.md", GATES,
                                    V4 / "reviews/independent_canonical_executor_correction_review_v4_3.md"]},
        "coordinate_map_path": blocks["coordinate_map_path"], "coordinate_map_sha256": blocks["coordinate_map_sha256"],
        "primary_interval_path": blocks["interval_tsv_path"], "primary_interval_sha256": blocks["interval_tsv_sha256"],
        "reference_sha256": blocks["reference_sha256"], "M_5_50_sha256": mfiles, "M_5_50_sum": M,
        "existing_raw_reference_metadata_counts": dict(counts),
        "reference_counts_are_not_actual_pair_retained_sets": True,
        "all200_reference_match_counts_are_not_full_or_two_step_mask_proof": True,
        "required_raw_BED_fresh_SHA": probe["raw_bed_metadata_only"]["historical_expected_sha256_not_freshly_verified"],
        "raw_BED_fresh_hash_achieved": False,
        "restricted_generator_definition": {
            "name": "503-EUR signed Gram-factor Gaussian summary-statistic qualification",
            "labels": ["S01", "S02", "S03", "S04"], "C_quota_credit": 0,
            "m_total": "pending sealed exact reference/allele/coordinate intersection after raw BED QC",
            "matrix": "R=G'G/503 within chromosomes; zero between chromosomes; no jackknife-boundary factorization",
            "dosage_rule": "reference mean imputation, centering and unit empirical variance; QC/exclusion criteria must be frozen first",
            "genetic_draw": "beta=U chol(Sigma_g)'/sqrt(m_total)",
            "noise_draw": "epsilon=G'H chol(Sigma_e)'/sqrt(503), Sigma_e unit-diagonal GWAS-noise correlation",
            "GWAS": "Z=sqrt(N)*R*beta+epsilon; shared X is exact same vector; independent keyed role substreams",
            "signal_covariance": "R^2 tensor Sigma_g/m_total", "noise_covariance": "R tensor Sigma_e",
            "C03_rank": "original whole historical master reference rank, never subset reindexing",
            "target_semantics": "integrated random-effect parameter ratios; root independent scientific review required before native simulation",
            "source_mismatch": "historical scores/weights/M retained exactly; finite source R need not agree with unbiased windowed historical L2",
            "not_established": ["population LD/long-range correlation length", "actual native retained-set coverage", "physical variant-N cohort overlap",
                                "binary/liability transformations", "conditional Schur models", "selection/transport validity"]},
        "proposed_sensitivity_rule_review_required": {
            "B100_B400": "same source/static MHC filter; Hamilton one perchr plus B-22 proportional count, fractional ties numericchr",
            "shift200": "same allocation; cut at BP index floor((j+0.5)*N/B_chr), start1/endmaxBP+1; reject repeatedBP, never split ties",
            "manifests_achieved": False, "generalized_native_executor_review_achieved": False},
        "future_resource_contract_review_required": {
            "worker_count": 1, "blas_threads": 1, "aggregate_RSS_limit_bytes": 2 * 1024**3,
            "internal_floor_bytes": 3 * 1024**3, "SSD_floor_bytes": 5 * 1024**3,
            "required_heavy_lock": "/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/native_heavy_worker.lock",
            "lock": "fcntl LOCK_EX|LOCK_NB for full owned worker/group lifetime; reviewed deferred catchable signals and cleanup guard",
            "decoded_chunk_max_SNPs": 2048, "503sample_chunk_float64_bytes": 8241152,
            "pair_fits_if_both_methods_and_four_partitions": 128000,
            "measured_native_throughput_or_total_RSS": None,
            "existing16control_capture_cap_bytes": 8 * 1024**3,
            "new_empirical_storage_cap": None, "reuse16control_cap_without_new_plan": False,
            "four_array_delete_result_bytes_lower_bound": 921600000,
            "every_scheduled_seed_status_delete_result_retained": True,
            "full_input_design_preservation_or_regeneration_rule": "pending explicit future admission"},
        "bounded_external_candidate_no_acquisition": {
            "source_type": "author-deposited S-LDSC EUR genotype archive candidate",
            "record_doi": zenodo["doi"], "metadata_path": str(zenodo_path), "metadata_sha256": sha(zenodo_path),
            "file_name": archive["key"], "source_bytes": archive["size"], "upstream_checksum": archive["checksum"],
            "content_url": archive["links"]["self"], "local_payload_verified": False,
            "planned_archive_cap_bytes": 512 * 1024**2, "planned_extracted_cap_bytes": 2 * 1024**3,
            "separate_existing_mirror_uncompressed_claim_bytes": 1593837597,
            "mirror_source_metadata_is_not_proof_of_Zenodo_member_identity": True,
            "actual_archive_members_samples_alleles_build_and_coverage_verified": False,
            "requires_pre_acquisition_source_method_resource_review": True,
            "historical_scores_weights_M_substitution_allowed": False},
        "fixture_checks": fixture["checks"], "empirical_calibration_checks": 0,
        "source_assessment_estimator_calls": 0, "calibrated_P_values_computed": 0,
        "failure_or_uncalibrated_result_may_not_be_dropped_or_replaced": True,
    }
    assert all(sha(p) == h for p, h in mfiles.items())
    PLAN.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"plan": str(PLAN), "sha256": sha(PLAN), "admitted": False, "M": M, "estimator_calls": 0}))


if __name__ == "__main__":
    main()
