#!/usr/bin/env python3
"""Seal and verify the pre-result Track B three-pair pleiotropy contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import Any


POLICY = Path("config/track_b_pleiotropy_policy.json")
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
PAIR_LOCK = Path("results/track_b/pair_manifest.lock.json")
DENSE_QC = Path("results/track_b/03_dense_input_qc.tsv")
DENSE_LOCK = Path("results/track_b/03_dense_input_qc.lock.json")
CONTRACT_LOCK = Path("results/track_b/pleiotropy/contract.lock.json")
INPUT_GATE_SCRIPT = Path("scripts/124_build_track_b_pleiotropy_input_gate.py")
LD_MATERIALIZER_SCRIPT = Path("scripts/127_prepare_track_b_pleiotropy_ld.py")
CONTRACT_SCRIPT = Path("scripts/123_track_b_pleiotropy_contract.py")

SHA256 = re.compile(r"[0-9a-f]{64}")
EXPECTED_PAIRS = [
    ("A", "snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    ("B", "insomnia", "adhd", "PRIMARY_DISCOVERY"),
    ("CONTROL", "insomnia", "frailty", "POSITIVE_CONTROL"),
]
EXPECTED_TRAITS = ["snoring", "parental_lifespan", "insomnia", "adhd", "frailty"]
EXPECTED_DENSE_SCHEMA = "SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N"
EXPECTED_PRIMARY_LABELS = [
    "PLACO_PRIMARY_HEADLINE",
    "PLACO_PAIRWISE_GWS",
    "PLACO_WITHIN_PAIR_FDR",
    "CONJFDR_INDEPENDENT",
    "PLACO_AND_CONJFDR_SAME_LD_LOCUS",
]
EXPECTED_CONTROL_LABELS = [
    "CONTROL_PLACO_GWS_RECOVERED",
    "CONTROL_PLACO_WITHIN_PAIR_FDR",
    "CONTROL_CONJFDR_RECOVERED",
    "CONTROL_BOTH_METHODS_SAME_LD_LOCUS",
]
EXPECTED_COMPARISON_LABELS = ["PLACO_AND_CONJFDR", "PLACO_ONLY", "CONJFDR_ONLY"]
EXPECTED_GATE_STATES = [
    "READY",
    "BLOCKED_UPSTREAM",
    "BLOCKED_BY_DATA",
    "BLOCKED_BY_SOFTWARE",
    "BLOCKED_BY_COMPUTE",
    "BLOCKED_BY_IMPLEMENTATION",
    "BLOCKED_BY_MULTIPLE_GATES",
]
EXPECTED_TERMINAL_STATES = ["COMPLETE_WITH_HITS", "TESTED_NO_HIT", "FAILED_QC"]
EXPECTED_SHARED_LOCI_PREFIX = [
    "pair_id", "lead_variant", "chr", "position", "PLACO_P", "FDR",
    "trait1_P", "trait2_P", "locus_start", "locus_end",
    "independent_signal", "annotations",
]
EXPECTED_LD_MEMBERS = {
    "bed": {
        "name": "g1000_eur.bed", "expected_bytes": 2_855_798_067,
        "sha256": "fadce7321e622ed20e3cb7798b3e54699c7f342752c642a14ddb982d93d8e08c",
        "records": 22_665_064,
    },
    "bim": {
        "name": "g1000_eur.bim", "expected_bytes": 658_724_714,
        "sha256": "b7dfeb30cef0573059ae43dd85a6924e46e418d6fe03df6b0ff1166c79cfa131",
        "records": 22_665_064,
    },
    "fam": {
        "name": "g1000_eur.fam", "expected_bytes": 12_575,
        "sha256": "79dee21226fd04b5a889cbbeacd08079e1740e7d63b0a4dea115527e8fea4dae",
        "records": 503,
    },
}
EXPECTED_LD_CHROMOSOME_COUNTS = {
    "1": 1_762_677, "2": 1_888_262, "3": 1_580_662, "4": 1_576_069,
    "5": 1_410_778, "6": 1_417_881, "7": 1_299_428, "8": 1_232_855,
    "9": 968_848, "10": 1_107_431, "11": 1_095_598, "12": 1_053_520,
    "13": 787_058, "14": 717_296, "15": 648_961, "16": 726_794,
    "17": 629_662, "18": 614_823, "19": 511_554, "20": 487_406,
    "21": 306_951, "22": 308_143, "23": 532_407,
}
EXPECTED_GRCH37_CHROMOSOME_LENGTHS = {
    "1": 249_250_621, "2": 243_199_373, "3": 198_022_430, "4": 191_154_276,
    "5": 180_915_260, "6": 171_115_067, "7": 159_138_663, "8": 146_364_022,
    "9": 141_213_431, "10": 135_534_747, "11": 135_006_516, "12": 133_851_895,
    "13": 115_169_878, "14": 107_349_540, "15": 102_531_392, "16": 90_354_753,
    "17": 81_195_210, "18": 78_077_248, "19": 59_128_983, "20": 63_025_520,
    "21": 48_129_895, "22": 51_304_566, "23": 155_270_560,
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def root_path(root: Path, relative: Path | str) -> Path:
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        fail(f"unsafe repository-relative path: {relative}")
    return root / path


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        fail(f"unreadable JSON artifact {path}: {error}")
    if not isinstance(value, dict):
        fail(f"JSON artifact must contain one object: {path}")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        fail(f"unreadable TSV artifact {path}: {error}")
    if not fields or len(fields) != len(set(fields)):
        fail(f"TSV artifact has an empty or duplicate schema: {path}")
    return fields, rows


def close(observed: Any, expected: float, label: str) -> None:
    if isinstance(observed, bool) or not isinstance(observed, (int, float)):
        fail(f"{label} must be numeric")
    if not math.isclose(float(observed), expected, rel_tol=1e-12, abs_tol=0.0):
        fail(f"{label} drifted: observed={observed!r} expected={expected!r}")


def validate_policy(root: Path) -> dict[str, Any]:
    policy = read_json(root_path(root, POLICY))
    if (
        policy.get("schema_version") != "sleep-atlas-track-b-pleiotropy-policy.1"
        or policy.get("analysis_id") != "track-b-v1.0-pleiotropy"
        or policy.get("analysis_build") != "hg19"
        or policy.get("ancestry") != "EUR"
        or policy.get("selection_timing") != "BEFORE_TRACK_B_PLEIOTROPY_RESULT_ACCESS"
        or policy.get("expected_pairs") != 3
        or policy.get("expected_primary_pairs") != 2
        or policy.get("expected_control_pairs") != 1
        or policy.get("expected_unique_traits") != 5
        or policy.get("pair_order") != ["A", "B", "CONTROL"]
        or policy.get("primary_pair_family") != ["A", "B"]
        or policy.get("control_pair_family") != ["CONTROL"]
    ):
        fail("Track B pleiotropy policy scope drifted")
    pairs = policy.get("pairs")
    observed_pairs = [
        (row.get("pair_id"), row.get("trait1"), row.get("trait2"), row.get("family_role"))
        for row in pairs
    ] if isinstance(pairs, list) and all(isinstance(row, dict) for row in pairs) else []
    if observed_pairs != EXPECTED_PAIRS:
        fail("Track B pleiotropy pair identities or family roles drifted")
    if not str(policy.get("pair_replacement_policy", "")).startswith("FORBIDDEN"):
        fail("Track B pleiotropy pair replacement must remain forbidden")

    upstream = policy.get("upstream_contracts", {})
    if upstream != {
        "pair_manifest": str(PAIR_MANIFEST),
        "pair_manifest_lock": str(PAIR_LOCK),
        "dense_qc": str(DENSE_QC),
        "dense_qc_lock": str(DENSE_LOCK),
    }:
        fail("Track B pleiotropy upstream-contract paths drifted")

    dense = policy.get("dense_input_contract", {})
    if (
        dense.get("required_columns") != EXPECTED_DENSE_SCHEMA.split(",")
        or dense.get("minimum_variant_rows_per_trait") != 5_000_000
        or dense.get("hapmap3_only_inputs_forbidden") is not True
        or dense.get("full_gzip_integrity_scan_required") is not True
        or dense.get("full_compressed_sha256_required") is not True
        or dense.get("full_row_count_required") is not True
        or dense.get("exact_schema_required") is not True
        or len(dense.get("alignment_provenance_counts", [])) != 10
        or "position-only" not in str(dense.get("alignment_rule", ""))
    ):
        fail("Track B genome-wide dense-input policy drifted")
    plausibility = dense.get("pair_materialization_plausibility", {})
    if (
        plausibility.get("minimum_aligned_eligible_variants") != 1_000_000
        or plausibility.get("minimum_fraction_of_smaller_dense_input") != 0.8
        or plausibility.get("required_autosomes") != list(range(1, 23))
        or plausibility.get("nonzero_eligible_variants_on_every_required_autosome") is not True
        or plausibility.get("per_autosome_provenance_required") is not True
        or plausibility.get("per_autosome_provenance_fields") != [
            "CHR", "trait1_rows", "trait2_rows", "exact_rsid_coordinate_matches",
            "allele_matches", "eligible_written",
        ]
        or "FAILED_QC" not in str(plausibility.get("failure_rule", ""))
        or "no PLACO task or scientific result" not in str(plausibility.get("failure_rule", ""))
    ):
        fail("Track B pair-materialization plausibility floors drifted")

    placo = policy.get("placo_plus", {})
    if (
        placo.get("version") != "0.2.0"
        or placo.get("source_path") != ".r-env/share/placo/PLACO_v0.2.0.R"
        or placo.get("source_sha256") != "fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124"
        or placo.get("upstream_commit") != "3ba3cae1d323ad117fb4540e620bcefa79f70663"
        or placo.get("full_p_ledger_required") is not True
        or placo.get("full_family_size_required_in_every_row") is not True
        or "every aligned eligible genome-wide variant" not in str(placo.get("within_pair_bh_family", ""))
        or "P=1" not in str(placo.get("within_pair_bh_family", ""))
        or "never supports novelty" not in str(placo.get("control_role", ""))
    ):
        fail("Track B PLACO+ family-accounting policy drifted")
    close(placo.get("primary_pair_family_headline_threshold"), 2.5e-8, "A/B headline threshold")
    close(placo.get("per_pair_genome_wide_threshold"), 5e-8, "per-pair threshold")
    close(placo.get("control_genome_wide_threshold"), 5e-8, "CONTROL threshold")
    close(placo.get("within_pair_bh_alpha"), 0.05, "within-pair BH alpha")
    if not math.isclose(
        2 * float(placo["primary_pair_family_headline_threshold"]),
        float(placo["per_pair_genome_wide_threshold"]), rel_tol=1e-12,
    ):
        fail("A/B headline threshold is not exactly 5e-8 divided by two primary pairs")

    conjfdr = policy.get("conjfdr", {})
    close(conjfdr.get("threshold"), 0.05, "conjFDR threshold")
    if (
        conjfdr.get("mode") != "conjunction_FDR"
        or conjfdr.get("independent_inclusion") is not True
        or "intersection is not required" not in str(conjfdr.get("independence_rule", ""))
        or conjfdr.get("code_commit") != "0da963cac22fe8de9030166d7aea974bcbb0a367"
        or conjfdr.get("runtime_provenance_path") != "ref/pleiofdr/runtime.provenance.json"
        or conjfdr.get("variant_template_sha256") != "06268420a0ec04e4529e832e1d4f4a53231b078cc5a3741a3eb215a5a4e1a9d5"
        or conjfdr.get("overlap_patch_sha256") != "a5ea51cafd08b783903733f6485272fd8d3563664b060ad31ae8a9266dced200"
        or conjfdr.get("overlap_correction_required") is not True
        or conjfdr.get("matlab_required") is not True
        or conjfdr.get("octave_primary_run_forbidden") is not True
    ):
        fail("Track B conjFDR independent-evidence policy drifted")

    ld = policy.get("ld_reference_policy", {})
    extracted_bytes = sum(member["expected_bytes"] for member in EXPECTED_LD_MEMBERS.values())
    if (
        ld.get("ancestry") != "EUR"
        or ld.get("build") != "GRCh37/hg19"
        or ld.get("source_release") != "1000_Genomes_Phase3_EUR_503_GRCh37_2018-09-17"
        or ld.get("source_archive_path")
        != "ref/interpretation/cell_type/magma_v1.10/g1000_eur.zip"
        or ld.get("source_archive_expected_bytes") != 511_626_945
        or ld.get("source_archive_sha256")
        != "83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954"
        or ld.get("source_members") != EXPECTED_LD_MEMBERS
        or ld.get("archive_allowed_members") != [
            "README", "g1000_eur.bed", "g1000_eur.bim", "g1000_eur.fam",
            "g1000_eur.synonyms",
        ]
        or ld.get("materialized_directory") != "ref/track_b/pleiotropy/g1000_eur"
        or ld.get("materialized_prefix")
        != "ref/track_b/pleiotropy/g1000_eur/g1000_eur"
        or ld.get("materialized_manifest")
        != "ref/track_b/pleiotropy/g1000_eur/g1000_eur.manifest.tsv"
        or ld.get("materialized_provenance")
        != "ref/track_b/pleiotropy/g1000_eur/g1000_eur.provenance.json"
        or ld.get("staging_root") != "ref/track_b/pleiotropy/.staging"
        or ld.get("materializer_script") != str(LD_MATERIALIZER_SCRIPT)
        or ld.get("staging_safety_bytes") != 1_073_741_824
        or ld.get("minimum_free_storage_bytes") != extracted_bytes + 1_073_741_824
        or ld.get("conflicting_activity_markers")
        != ["ref/lava/ukb_v1.1/download_manifest.tsv.tmp"]
        or ld.get("sample_count") != 503
        or ld.get("variant_count") != 22_665_064
        or ld.get("reference_chromosomes") != list(range(1, 24))
        or ld.get("analysis_autosomes") != list(range(1, 23))
        or ld.get("nonanalysis_chromosomes_retained_for_source_fidelity") != [23]
        or ld.get("expected_autosomal_variant_count") != 22_132_657
        or ld.get("expected_reference_chromosome_variant_counts")
        != EXPECTED_LD_CHROMOSOME_COUNTS
        or ld.get("grch37_chromosome_lengths") != EXPECTED_GRCH37_CHROMOSOME_LENGTHS
        or ld.get("full_genome_wide_reference_required") is not True
        or ld.get("hapmap3_subset_forbidden") is not True
        or "untouched chr1-23" not in str(ld.get("raw_reference_fidelity_rule", ""))
        or "do not derive or relabel an autosome-only BED/BIM/FAM"
        not in str(ld.get("raw_reference_fidelity_rule", ""))
        or "22,132,657" not in str(ld.get("analysis_chromosome_rule", ""))
        or "chr23/X is always forbidden" not in str(ld.get("analysis_chromosome_rule", ""))
        or "Byte-identical extraction" not in str(ld.get("materialization_transform_rule", ""))
        or "no LD-block splitting" not in str(ld.get("materialization_transform_rule", ""))
        or "SNP reduction" not in str(ld.get("materialization_transform_rule", ""))
        or ld.get("materialization_resource_metrics_required") != [
            "staging_extraction_validation_wall_seconds", "peak_rss_bytes",
            "streaming_chunk_bytes", "extracted_uncompressed_bytes",
        ]
        or ld.get("plink_path")
        != "ref/interpretation/causal/plink_1.9_b7.11/plink"
        or ld.get("plink_version") != "PLINK v1.9.0-b.7.11 64-bit (19 Aug 2025)"
        or ld.get("plink_sha256")
        != "5e49297d82c680f2ccfe45c62c974773f6c6e360e5d8a453f831ad95a480b642"
        or ld.get("clump_r2_strict_upper_bound") != 0.1
        or ld.get("clump_window_kb") != 1000
        or "r2<0.1 within 1 Mb" not in str(ld.get("lead_independence_rule", ""))
        or "r2>=0.6" not in str(ld.get("cross_method_locus_rule", ""))
        or "distance-only" not in str(ld.get("cross_method_locus_rule", ""))
        or "broad LAVA-block" not in str(ld.get("cross_method_locus_rule", ""))
    ):
        fail("Track B full EUR hg19 LD-reference policy drifted")
    if (
        sum(EXPECTED_LD_CHROMOSOME_COUNTS.values()) != ld["variant_count"]
        or sum(EXPECTED_LD_CHROMOSOME_COUNTS[str(value)] for value in range(1, 23))
        != ld["expected_autosomal_variant_count"]
        or EXPECTED_LD_MEMBERS["bed"]["expected_bytes"]
        != 3 + math.ceil(ld["sample_count"] / 4) * ld["variant_count"]
    ):
        fail("Track B raw PLINK dimensions or autosomal counts do not reconcile")

    evidence = policy.get("evidence_union", {})
    if (
        evidence.get("intersection_required") is not False
        or evidence.get("primary_labels_in_order") != EXPECTED_PRIMARY_LABELS
        or evidence.get("control_labels_in_order") != EXPECTED_CONTROL_LABELS
        or evidence.get("comparison_labels_in_order") != EXPECTED_COMPARISON_LABELS
        or "only when both methods completed" not in str(evidence.get("comparison_only_label_rule", ""))
        or "do not emit an ONLY label" not in str(evidence.get("comparison_only_label_rule", ""))
        or "union" not in str(evidence.get("primary_pair_inclusion_rule", "")).lower()
        or "never merge CONTROL" not in str(evidence.get("control_inclusion_rule", ""))
    ):
        fail("Track B method-union evidence policy drifted")

    states = policy.get("method_state_semantics", {})
    if (
        states.get("allowed_gate_states") != EXPECTED_GATE_STATES
        or states.get("allowed_terminal_result_states") != EXPECTED_TERMINAL_STATES
        or "not P=1" not in str(states.get("blocked_rule", ""))
        or "complete checksum-sealed genome-wide scan" not in str(states.get("zero_hit_rule", ""))
        or "Never publish an empty scientific result table" not in str(states.get("partial_completion_rule", ""))
    ):
        fail("Track B blocked-method semantics drifted")

    execution = policy.get("execution_order", {})
    if (
        execution.get("pre_result_freeze_required") is not True
        or execution.get("replication_terminal_before_method_execution") is not True
        or execution.get("local_analysis_terminal_before_method_execution") is not True
        or not str(execution.get("placo_and_conjfdr_ordering", "")).startswith("INDEPENDENT")
        or execution.get("generic_396_pair_artifacts_as_substitute") != "FORBIDDEN"
    ):
        fail("Track B pleiotropy execution-order policy drifted")

    result_schemas = policy.get("required_future_result_schemas", {})
    future_paths = policy.get("future_result_paths", {})
    if (
        result_schemas.get("shared_loci_08_required_prefix") != EXPECTED_SHARED_LOCI_PREFIX
        or future_paths.get("placo_plus_variants_07") != "results/track_b/07_placo_plus_variants.tsv"
        or future_paths.get("shared_loci_08") != "results/track_b/08_shared_loci.tsv"
        or future_paths.get("pleiotropy_comparison_09") != "results/track_b/09_pleiotropy_comparison.tsv"
        or future_paths.get("control_directory") != "results/track_b/control_insomnia_frailty/"
        or future_paths.get("result_provenance")
        != "results/track_b/pleiotropy/results/results.provenance.json"
    ):
        fail("Track B canonical 07/08/09/control paths or 08 schema drifted")

    claims = policy.get("claim_limits", {})
    if set(claims) != {"primary", "placo", "conjfdr", "cross_method", "control", "blocked"}:
        fail("Track B pleiotropy claim-limit family is incomplete")
    joined_claims = " ".join(str(value) for value in claims.values())
    if (
        "does not establish a shared causal variant" not in joined_claims
        or "never supports novelty" in joined_claims
        or "excluded from novelty" not in joined_claims
        or "cannot be described as null" not in joined_claims
    ):
        fail("Track B pleiotropy claim limits drifted")
    return policy


def validate_upstream(root: Path, policy: dict[str, Any]) -> dict[str, Any]:
    pair_path = root_path(root, policy["upstream_contracts"]["pair_manifest"])
    pair_lock_path = root_path(root, policy["upstream_contracts"]["pair_manifest_lock"])
    dense_path = root_path(root, policy["upstream_contracts"]["dense_qc"])
    dense_lock_path = root_path(root, policy["upstream_contracts"]["dense_qc_lock"])
    pair_hash = sha256(pair_path)
    pair_lock = read_json(pair_lock_path)
    if (
        pair_lock.get("schema_version") != "1.0.0"
        or pair_lock.get("output_sha256", {}).get(str(PAIR_MANIFEST)) != pair_hash
        or not SHA256.fullmatch(str(pair_lock.get("pair_b_identity_sha256", "")))
        or not str(pair_lock.get("pair_replacement_policy", "")).startswith("FORBIDDEN")
    ):
        fail("frozen Track B pair lock does not bind the live pair manifest")
    pair_fields, pair_rows = read_tsv(pair_path)
    required_pair_fields = {"pair_id", "sleep_trait", "external_trait", "ancestry", "build", "dense_data_readiness"}
    if not required_pair_fields.issubset(pair_fields):
        fail("frozen Track B pair manifest schema is incomplete")
    observed_pairs = [
        (row["pair_id"], row["sleep_trait"], row["external_trait"])
        for row in pair_rows
    ]
    if observed_pairs != [(row[0], row[1], row[2]) for row in EXPECTED_PAIRS]:
        fail("live Track B pair manifest differs from the three-pair policy")
    for row in pair_rows:
        if (
            row["ancestry"] != "sleep=EUR;external=EUR"
            or row["build"] != "sleep=hg19;external=hg19"
            or row["dense_data_readiness"] != "READY_FULL_SUMSTATS_BOTH"
        ):
            fail(f"frozen pair is not EUR hg19 dense-ready: {row['pair_id']}")

    dense_hash = sha256(dense_path)
    dense_lock = read_json(dense_lock_path)
    if (
        dense_lock.get("schema_version") != 1
        or dense_lock.get("pair_manifest_sha256") != pair_hash
        or dense_lock.get("dense_qc_sha256") != dense_hash
        or dense_lock.get("trait_count") != 5
        or dense_lock.get("unique_dense_trait_family") != EXPECTED_TRAITS
        or dense_lock.get("all_are_dense_not_hapmap3_only") is not True
        or dense_lock.get("all_analysis_builds") != ["hg19"]
    ):
        fail("Track B dense-QC lock differs from the live five-trait family")
    dense_fields, dense_rows = read_tsv(dense_path)
    required_dense_fields = {
        "pair_ids", "trait_id", "ancestry", "analysis_build", "dense_file",
        "dense_file_sha256", "compressed_size_bytes", "rows_out", "observed_schema",
        "INFO_status", "overall_QC",
    }
    if not required_dense_fields.issubset(dense_fields):
        fail("Track B dense-QC schema is incomplete")
    if [row["trait_id"] for row in dense_rows] != EXPECTED_TRAITS:
        fail("Track B dense-QC trait order drifted")
    expected_pair_ids = ["A", "A", "B;CONTROL", "B", "CONTROL"]
    for row, pair_ids in zip(dense_rows, expected_pair_ids, strict=True):
        try:
            size = int(row["compressed_size_bytes"])
            count = int(row["rows_out"])
        except ValueError:
            fail(f"invalid frozen dense size/count for {row['trait_id']}")
        if (
            row["pair_ids"] != pair_ids
            or row["ancestry"] != "EUR"
            or row["analysis_build"] != "hg19"
            or row["observed_schema"] != EXPECTED_DENSE_SCHEMA
            or not SHA256.fullmatch(row["dense_file_sha256"])
            or size <= 0
            or count < int(policy["dense_input_contract"]["minimum_variant_rows_per_trait"])
            or not row["overall_QC"].startswith("PASS_DENSE")
        ):
            fail(f"Track B dense-QC identity/readiness drifted for {row['trait_id']}")
    return {
        "pair_manifest_sha256": pair_hash,
        "pair_manifest_lock_sha256": sha256(pair_lock_path),
        "dense_qc_sha256": dense_hash,
        "dense_qc_lock_sha256": sha256(dense_lock_path),
        "pairs": [
            {"pair_id": pair_id, "trait1": trait1, "trait2": trait2, "family_role": role}
            for pair_id, trait1, trait2, role in EXPECTED_PAIRS
        ],
        "dense_traits": [
            {
                "trait_id": row["trait_id"], "path": row["dense_file"],
                "bytes": int(row["compressed_size_bytes"]),
                "rows": int(row["rows_out"]), "sha256": row["dense_file_sha256"],
                "schema": row["observed_schema"],
            }
            for row in dense_rows
        ],
    }


def result_evidence(root: Path, policy: dict[str, Any]) -> list[Path]:
    result_paths = policy["future_result_paths"]
    candidates = [
        root_path(root, result_paths["placo_plus_variants_07"]),
        root_path(root, result_paths["shared_loci_08"]),
        root_path(root, result_paths["pleiotropy_comparison_09"]),
        root_path(root, result_paths["result_provenance"]),
    ]
    evidence = [path for path in candidates if path.exists()]
    for key in ("placo_full_ledger_directory", "conjfdr_full_ledger_directory", "control_directory"):
        directory = root_path(root, result_paths[key])
        if directory.is_dir() and any(directory.iterdir()):
            evidence.append(directory)
    return evidence


def contract_payload(root: Path) -> dict[str, Any]:
    policy = validate_policy(root)
    upstream = validate_upstream(root, policy)
    return {
        "schema_version": "sleep-atlas-track-b-pleiotropy-contract.1",
        "analysis_id": policy["analysis_id"],
        "selection_timing": policy["selection_timing"],
        "pleiotropy_results_accessed_before_contract_freeze": False,
        "policy": str(POLICY),
        "policy_sha256": sha256(root_path(root, POLICY)),
        "upstream": upstream,
        "script_sha256": {
            str(CONTRACT_SCRIPT): sha256(root_path(root, CONTRACT_SCRIPT)),
            str(INPUT_GATE_SCRIPT): sha256(root_path(root, INPUT_GATE_SCRIPT)),
            str(LD_MATERIALIZER_SCRIPT): sha256(root_path(root, LD_MATERIALIZER_SCRIPT)),
        },
        "primary_pair_family": ["A", "B"],
        "control_pair_family": ["CONTROL"],
        "thresholds": {
            "primary_pair_family_headline": 2.5e-8,
            "per_pair_genome_wide": 5e-8,
            "control_genome_wide": 5e-8,
            "within_pair_bh_alpha": 0.05,
            "conjfdr": 0.05,
        },
        "method_union_not_intersection": True,
        "generic_396_pair_substitution_forbidden": True,
        "scientific_result_substitution_policy": "FORBIDDEN; gates and locks are not scientific results",
    }


def exclusive_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            fail(f"immutable Track B pleiotropy contract already exists: {path}")
    finally:
        temporary.unlink(missing_ok=True)


def seal_contract(root: Path) -> None:
    policy = validate_policy(root)
    evidence = result_evidence(root, policy)
    if evidence:
        fail(
            "refusing to freeze the Track B pleiotropy contract after result access: "
            + ", ".join(str(path.relative_to(root)) for path in evidence[:3])
        )
    path = root_path(root, CONTRACT_LOCK)
    if path.exists():
        verify_contract(root)
        return
    payload = contract_payload(root)
    exclusive_json(path, payload)
    print("TRACK_B_PLEIOTROPY_CONTRACT_SEALED pairs=3 primary=2 control=1 results_accessed=false")


def verify_contract(root: Path) -> dict[str, Any]:
    path = root_path(root, CONTRACT_LOCK)
    observed = read_json(path)
    expected = contract_payload(root)
    if observed != expected:
        fail("Track B pleiotropy contract lock differs from policy, scripts, or upstream locks")
    print("TRACK_B_PLEIOTROPY_CONTRACT_VERIFIED pairs=A,B control=CONTROL")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--seal-contract", action="store_true")
    action.add_argument("--verify-contract", action="store_true")
    action.add_argument("--validate-policy", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    if args.seal_contract:
        seal_contract(root)
    elif args.verify_contract:
        verify_contract(root)
    else:
        policy = validate_policy(root)
        validate_upstream(root, policy)
        print("TRACK_B_PLEIOTROPY_POLICY_VALIDATED")


if __name__ == "__main__":
    main()
