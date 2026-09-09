from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    specification = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        specification.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


RUNTIME = load_module("track_b_finemap_runtime_148_test", "scripts/148_track_b_finemapping_runtime.py")
COORDINATOR = load_module(
    "track_b_finemap_coordinator_151_test", "scripts/151_run_track_b_finemapping_sequential.py"
)
POLICY = json.loads((ROOT / "config/track_b_finemapping_policy.json").read_text())


def write_ld(path: Path, names: list[str], matrix: np.ndarray, row_names: list[str] | None = None):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", *names])
        for name, row in zip(row_names or names, matrix, strict=True):
            writer.writerow([name, *[format(float(value), ".17g") for value in row]])


def order_rows(names: list[str]):
    return [
        {"SNP": name, "CHR": "1", "BP": str(index + 1), "A1": "A", "A2": "C"}
        for index, name in enumerate(names)
    ]


@pytest.fixture
def measured_ld_context(tmp_path, monkeypatch):
    lease = (tmp_path / "worker.lease").open("w+b")
    monkeypatch.setenv("TRACK_B_FINEMAP_MEASURED_PHASE", "MATERIALIZATION")
    monkeypatch.setenv("TRACK_B_FINEMAP_WORKER_LEASE_FD", str(lease.fileno()))
    monkeypatch.setattr(RUNTIME.os, "getsid", lambda _: RUNTIME.os.getpid())
    try:
        yield
    finally:
        lease.close()


def coloc_task_validation(
    *, rss1="0.01", status1="PASS", rss2="0.02", status2="PASS",
    n_status1=None, n_status2=None,
):
    task = {
        "analysis_id": RUNTIME.ANALYSIS_ID, "pair_id": "A",
        "family_role": "PRIMARY_DISCOVERY", "locus_entry_id": "A_LOC1",
        "ld_block_id": "LOC1", "CHR": "1", "START": "1", "STOP": "10",
        "trait1": "snoring", "trait2": "parental_lifespan",
        "reference_id": "LAVA_UKB_v1.1_EUR_GRCh37",
        "variant_order_sha256": "a" * 64, "variant_count": "2",
        "trait1_N_dispersion_status": n_status1 or RUNTIME.N_DISPERSION_CONSTANT,
        "trait2_N_dispersion_status": n_status2 or RUNTIME.N_DISPERSION_CONSTANT,
    }
    diagnostics = {
        "TRAIT1_SLEEP": {
            "diagnostic_status": status1, "rss_ld_s": rss1, "model_converged": "TRUE",
        },
        "TRAIT2_EXTERNAL": {
            "diagnostic_status": status2, "rss_ld_s": rss2, "model_converged": "TRUE",
        },
    }
    validation = {
        "task": task, "order_rows": [{"SNP": "rs1"}, {"SNP": "rs2"}],
        "policy": POLICY,
    }
    return validation, diagnostics


def probability_rows(values_by_prior: dict[float, tuple[float, float]]):
    rows = []
    for p12 in (1e-6, 5e-6, 1e-5, 5e-5):
        h3, h4 = values_by_prior[p12]
        remainder = 1 - h3 - h4
        assert remainder >= 0
        h0, h1, h2 = remainder / 3, remainder / 3, remainder / 3
        ratio = "Inf" if h3 == 0 and h4 > 0 else ("NA" if h3 == h4 == 0 else format(h4 / h3, ".17g"))
        rows.append({
            "analysis_id": RUNTIME.ANALYSIS_ID, "pair_id": "A",
            "family_role": "PRIMARY_DISCOVERY", "locus_entry_id": "A_LOC1",
            "locus_id": "LOC1", "CHR": "1", "START": "1", "STOP": "10",
            "trait1": "snoring", "trait2": "parental_lifespan",
            "signal1": "L1", "signal2": "L2", "coloc_method": "coloc.susie",
            "p1": "0.0001", "p2": "0.0001", "p12": format(p12, ".17g"),
            "prior_role": "PRIMARY" if p12 == 1e-5 else "SENSITIVITY",
            "nsnps": "2", "PP_H0": format(h0, ".17g"), "PP_H1": format(h1, ".17g"),
            "PP_H2": format(h2, ".17g"), "PP_H3": format(h3, ".17g"),
            "PP_H4": format(h4, ".17g"), "PP_H4_over_PP_H3": ratio,
            "top_shared_variant": "rs1", "top_shared_variant_PP_H4": "0.7",
            "fine_mapping_qc": "PASS", "ld_qc": "PASS",
            "engine_status": "COLOC_SUSIE_COMPLETE", "error": "NA",
            "claim_limit": POLICY["claim_limits"]["trait_coloc"],
        })
    return rows


def test_scalar_n_exact_when_constant():
    assert RUNTIME.scalar_n([1000, 1000, 1000]) == (
        1000.0, RUNTIME.SCALAR_N_RULE_CONSTANT
    )


def test_scalar_n_is_unrounded_complete_vector_median_and_claim_capped_when_varying():
    value, rule = RUNTIME.scalar_n([900.25, 10, 700.5, 2000])
    assert value == (700.5 + 900.25) / 2
    assert rule == RUNTIME.SCALAR_N_RULE_MEDIAN
    diagnostics = RUNTIME.per_snp_n_diagnostics([900.25, 10, 700.5, 2000])
    assert diagnostics == {
        "count": 4, "min": 10.0, "q1": 527.875, "median": 800.375,
        "q3": 1175.1875, "max": 2000.0, "max_to_min_ratio": 200.0,
        "fraction_below_90pct_max": 0.75, "fraction_below_50pct_max": 0.75,
        "dispersion_status": RUNTIME.N_DISPERSION_VARIABLE,
    }
    assert RUNTIME.per_snp_n_sha256([1, 2]) != RUNTIME.per_snp_n_sha256([2, 1])


@pytest.mark.parametrize("values", [[], [1, float("nan")], [0, 1], [-1, 1]])
def test_scalar_n_rejects_incomplete_or_invalid_vector(values):
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.scalar_n(values)


def test_exact_allele_orientation_and_eaf_swap():
    assert RUNTIME.oriented_alleles("A", "C", 0.2, "A", "C") == (1, 0.2, "DIRECT")
    sign, eaf, orientation = RUNTIME.oriented_alleles("C", "A", 0.2, "A", "C")
    assert sign == -1
    assert eaf == pytest.approx(0.8)
    assert orientation == "SWAP"


def test_unambiguous_complement_orientation_preserves_or_swaps_eaf_exactly():
    assert RUNTIME.oriented_alleles("T", "G", 0.2, "A", "C") == (
        1, 0.2, "COMPLEMENT",
    )
    sign, eaf, orientation = RUNTIME.oriented_alleles("G", "T", 0.2, "A", "C")
    assert sign == -1
    assert eaf == pytest.approx(0.8)
    assert orientation == "COMPLEMENT_SWAP"


@pytest.mark.parametrize(
    "source,reference",
    [("A/T", "A/T"), ("C/G", "G/C"), ("T/C", "A/C"), ("AA/C", "A/C"), ("A/A", "A/C")],
)
def test_alleles_reject_palindrome_complement_multibase_and_mismatch(source, reference):
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.oriented_alleles(*source.split("/"), 0.2, *reference.split("/"))


def test_dense_materialization_flips_beta_and_eaf(tmp_path):
    path = tmp_path / "dense.tsv"
    path.write_text(
        "SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\n"
        "rs1\t1\t100\tC\tA\t0.2\t0.4\t0.1\t0.01\t1000\n"
    )
    reference = {"rs1": {"BP": "100", "A1": "A", "A2": "C"}}
    rows, counts = RUNTIME.extract_dense_block(
        path, reference, chromosome=1, minimum_maf=0.01,
        info_status="ABSENT_IN_RELEASE_SOURCE_LEVEL_ONLY",
    )
    assert float(rows["rs1"]["BETA"]) == -0.4
    assert float(rows["rs1"]["EAF"]) == pytest.approx(0.8)
    assert rows["rs1"]["INFO"] == "NA"
    assert counts["allele_swaps"] == 1


def test_dense_materialization_tracks_complement_without_false_eaf_flip(tmp_path):
    path = tmp_path / "dense.tsv"
    path.write_text(
        "SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\n"
        "rs1\t1\t100\tT\tG\t0.2\t0.4\t0.1\t0.01\t1000\n"
    )
    reference = {"rs1": {"BP": "100", "A1": "A", "A2": "C"}}
    rows, counts = RUNTIME.extract_dense_block(
        path, reference, chromosome=1, minimum_maf=0.01,
        info_status="ABSENT_IN_RELEASE_SOURCE_LEVEL_ONLY",
    )
    assert float(rows["rs1"]["BETA"]) == 0.4
    assert float(rows["rs1"]["EAF"]) == 0.2
    assert counts["allele_complements"] == 1
    assert counts.get("allele_swaps", 0) == 0


@pytest.mark.parametrize(
    "second",
    [
        "rs1\t1\t100\tA\tC\t0.2\t0.4\t0.1\t0.01\t1000\n",
        "rs1\t1\t101\tA\tC\t0.2\t0.4\t0.1\t0.01\t1000\n",
        "rs1\t1\t100\tA\tC\t0.001\t0.4\t0.1\t0.01\t1000\n",
        "rs1\t1\t100\tA\tC\t0.2\tNaN\t0.1\t0.01\t1000\n",
    ],
)
def test_dense_materialization_rejects_duplicate_coordinate_maf_and_nonfinite(tmp_path, second):
    path = tmp_path / "dense.tsv"
    first = "rs1\t1\t100\tA\tC\t0.2\t0.4\t0.1\t0.01\t1000\n"
    # A duplicate case needs two rows; other cases replace the valid first row.
    rows = first + second if second == first else second
    path.write_text("SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\n" + rows)
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.extract_dense_block(
            path, {"rs1": {"BP": "100", "A1": "A", "A2": "C"}},
            chromosome=1, minimum_maf=0.01, info_status="ROW_WISE_FILTERED",
        )


def test_signed_ld_full_stream_accepts_negative_correlation(tmp_path, measured_ld_context):
    names = ["rs1", "rs2", "rs3"]
    matrix = np.array([[1, -0.2, 0.1], [-0.2, 1, 0.3], [0.1, 0.3, 1]])
    path = tmp_path / "ld.tsv.gz"
    write_ld(path, names, matrix)
    result = RUNTIME.deep_validate_signed_ld(
        path, order_rows(names), policy=POLICY, scratch_directory=tmp_path,
    )
    assert result["variant_count"] == 3
    assert result["matrix_kind"] == "SIGNED_PEARSON_CORRELATION"
    assert not list(tmp_path.glob(".ld-validate-*.bin"))


def test_signed_ld_spectral_validation_requires_explicit_measured_process(tmp_path):
    path = tmp_path / "ld.tsv.gz"
    write_ld(path, ["rs1"], np.eye(1))
    with pytest.raises(RUNTIME.FineMappingError, match="measured-phase"):
        RUNTIME.deep_validate_signed_ld(
            path, order_rows(["rs1"]), policy=POLICY, scratch_directory=tmp_path,
        )
    assert not list(tmp_path.glob(".ld-validate-*.bin"))


@pytest.mark.parametrize("defect", ["truncated", "order", "symmetry", "diagonal", "psd"])
def test_signed_ld_rejects_truncation_order_symmetry_diagonal_and_non_psd(
    tmp_path, measured_ld_context, defect,
):
    names = ["rs1", "rs2", "rs3"]
    matrix = np.eye(3)
    row_names = names.copy()
    if defect == "symmetry":
        matrix[0, 1] = 0.2
    elif defect == "diagonal":
        matrix[0, 0] = 0.9
    elif defect == "psd":
        matrix = np.array([[1, 0.9, 0.9], [0.9, 1, -0.9], [0.9, -0.9, 1]])
    elif defect == "order":
        row_names[0], row_names[1] = row_names[1], row_names[0]
    path = tmp_path / "ld.tsv.gz"
    if defect == "truncated":
        with gzip.open(path, "wt") as handle:
            handle.write("SNP\trs1\trs2\trs3\nrs1\t1\t0\t0\nrs2\t0\t1\t0\n")
    else:
        write_ld(path, names, matrix, row_names)
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.deep_validate_signed_ld(
            path, order_rows(names), policy=POLICY, scratch_directory=tmp_path,
        )


def test_no_5000_cap_or_thinning_and_exact_admission_formula(monkeypatch):
    assert POLICY["dense_input_contract"]["maximum_locus_variants"] is None
    assert POLICY["ram_aware_execution_contract"]["locus_splitting_forbidden"] is True
    assert POLICY["ram_aware_execution_contract"]["variant_thinning_for_compute_forbidden"] is True
    expected = max(4 * 1024**3, 64 * 6001 * 6001 + 2 * 1024**3)
    assert RUNTIME.estimate_peak_rss_bytes(6001, POLICY) == expected
    monkeypatch.setattr(RUNTIME, "host_physical_memory_bytes", lambda: 8 * 1024**3)
    blocked = RUNTIME.memory_admission(10000, POLICY, available_bytes=8 * 1024**3)
    assert blocked["admission_status"] == "BLOCKED_BY_COMPUTE"


def test_materialization_admission_uses_real_explicit_larger_envelope(monkeypatch):
    monkeypatch.setattr(RUNTIME, "host_physical_memory_bytes", lambda: 32 * 1024**3)
    with pytest.raises(RUNTIME.FineMappingError, match="explicit continuation"):
        RUNTIME.memory_admission(
            10000, POLICY, available_bytes=32 * 1024**3,
            configured_envelope_bytes=16 * 1024**3,
        )
    admitted = RUNTIME.memory_admission(
        10000, POLICY, available_bytes=32 * 1024**3,
        configured_envelope_bytes=16 * 1024**3, larger_host_continuation=True,
    )
    assert admitted["physical_envelope_bytes"] == 16 * 1024**3
    assert admitted["configured_envelope_bytes"] == 16 * 1024**3
    assert admitted["admission_status"] == "ADMITTED_WHOLE_LOCUS"


def test_materialization_resource_continuation_requires_strictly_larger_host_and_envelope():
    prior = [{
        "terminal_status": "FAILED_RESOURCE_OOM",
        "receipt": {
            "host_physical_memory_bytes": 8 * 1024**3,
            "effective_envelope_bytes": 8 * 1024**3,
        },
    }]
    with pytest.raises(COORDINATOR.CoordinatorError, match="strictly larger"):
        COORDINATOR.validate_larger_host_materialization_continuation(
            prior, current_physical_bytes=8 * 1024**3,
            configured_envelope_bytes=16 * 1024**3,
        )
    effective = COORDINATOR.validate_larger_host_materialization_continuation(
        prior, current_physical_bytes=32 * 1024**3,
        configured_envelope_bytes=16 * 1024**3,
    )
    assert effective == 16 * 1024**3
    with pytest.raises(COORDINATOR.CoordinatorError, match="only retained resource"):
        COORDINATOR.validate_larger_host_materialization_continuation(
            [{"terminal_status": "FAILED_VALIDATION", "receipt": {}}],
            current_physical_bytes=32 * 1024**3,
            configured_envelope_bytes=16 * 1024**3,
        )


def test_r_engine_continuation_requires_resource_predecessor_and_strict_effective_growth(
    tmp_path, monkeypatch,
):
    with pytest.raises(COORDINATOR.CoordinatorError, match="retained resource failure"):
        COORDINATOR.validate_larger_host_engine_continuation(
            [], current_physical_bytes=32 * 1024**3,
            configured_envelope_bytes=16 * 1024**3,
        )
    prior = [{
        "terminal_status": "FAILED_RESOURCE_OOM",
        "seal": {"attempt_id": "attempt-1"},
        "seal_identity": {"path": "run/attempt-1/run.seal.json", "bytes": 1, "sha256": "a" * 64},
        "provenance": {"memory_envelope_bytes": 16 * 1024**3},
        "provenance_identity": {
            "path": "run/attempt-1/run.provenance.json", "bytes": 1, "sha256": "b" * 64,
        },
        "resource_identity": {
            "path": "run/attempt-1/resource_metrics.tsv", "bytes": 1, "sha256": "c" * 64,
        },
        "resource_row": {"host_physical_memory_bytes": str(16 * 1024**3)},
        "run_rel": "results/track_b/finemapping/runs/A_LOC1/attempt-1",
    }]
    with pytest.raises(COORDINATOR.CoordinatorError, match="strictly larger"):
        COORDINATOR.validate_larger_host_engine_continuation(
            prior, current_physical_bytes=32 * 1024**3,
            configured_envelope_bytes=9 * 1024**3,
        )
    continuation = COORDINATOR.validate_larger_host_engine_continuation(
        prior, current_physical_bytes=32 * 1024**3,
        configured_envelope_bytes=24 * 1024**3,
    )
    assert continuation["current_effective_envelope_bytes"] == 24 * 1024**3
    assert continuation["predecessors"][0]["run_seal"]["sha256"] == "a" * 64
    assert continuation["predecessors"][0]["run_provenance"]["sha256"] == "b" * 64

    monkeypatch.setattr(COORDINATOR, "_existing_runs", lambda *args, **kwargs: [])
    with pytest.raises(COORDINATOR.CoordinatorError, match="retained resource failure"):
        COORDINATOR.execute_locus(
            object(), tmp_path, {}, {"locus_entry_id": "A_LOC1"},
            configured_envelope_bytes=16 * 1024**3,
            larger_host_continuation=True,
        )


def test_r_engine_continuation_verifier_binds_exact_predecessor_evidence(tmp_path):
    locus = "A_LOC1"
    fingerprint = "f" * 64
    prior_attempt = "20260101T000000000000Z.aaaa"
    current_attempt = "20260102T000000000000Z.bbbb"
    parent = tmp_path / COORDINATOR.RUN_ROOT_REL / locus
    prior_directory = parent / prior_attempt
    current_directory = parent / current_attempt
    prior_directory.mkdir(parents=True)
    current_directory.mkdir()
    prior_provenance = {
        "run_fingerprint": fingerprint, "attempt_id": prior_attempt,
        "locus_entry_id": locus, "terminal_status": "FAILED_RESOURCE_OOM",
        "memory_envelope_bytes": 16 * 1024**3,
        "r_engine_measurement": {"configured_envelope_bytes": 16 * 1024**3},
    }
    provenance_path = prior_directory / "run.provenance.json"
    provenance_path.write_text(json.dumps(prior_provenance, sort_keys=True) + "\n")
    fields = POLICY["ram_aware_execution_contract"]["resource_metrics_schema"]
    metric = {field: "NA" for field in fields}
    metric.update({
        "analysis_id": COORDINATOR.ANALYSIS_ID, "run_fingerprint": fingerprint,
        "locus_entry_id": locus, "pair_id": "A", "attempt_id": prior_attempt,
        "host_physical_memory_bytes": str(16 * 1024**3),
        "terminal_status": "FAILED_RESOURCE_OOM",
    })
    metrics_path = prior_directory / "resource_metrics.tsv"
    metrics_path.write_bytes(RUNTIME.tsv_bytes(fields, [metric]))
    provenance_identity = RUNTIME.stable_identity(tmp_path, provenance_path)
    resource_identity = RUNTIME.stable_identity(tmp_path, metrics_path)
    seal = {
        "schema_version": COORDINATOR.RUN_SCHEMA, "analysis_id": COORDINATOR.ANALYSIS_ID,
        "run_fingerprint": fingerprint, "attempt_id": prior_attempt,
        "locus_entry_id": locus, "terminal_status": "FAILED_RESOURCE_OOM",
        "records": {
            "run.provenance.json": {
                "bytes": provenance_identity["bytes"], "sha256": provenance_identity["sha256"],
            },
            "resource_metrics.tsv": {
                "bytes": resource_identity["bytes"], "sha256": resource_identity["sha256"],
            },
        },
    }
    (prior_directory / "run.seal.json").write_text(json.dumps(seal, sort_keys=True) + "\n")
    materialized = {"policy": POLICY}
    predecessors = COORDINATOR._bound_engine_predecessors(
        RUNTIME, tmp_path, materialized, current_directory.relative_to(tmp_path),
        locus=locus, fingerprint=fingerprint,
    )
    current_provenance = {
        "larger_host_continuation": True,
        "engine_continuation_predecessors": predecessors,
        "memory_envelope_bytes": 24 * 1024**3,
    }
    COORDINATOR._verify_engine_continuation_provenance(
        RUNTIME, tmp_path, materialized, current_directory.relative_to(tmp_path),
        locus=locus, fingerprint=fingerprint, provenance=current_provenance,
        measurement={"configured_envelope_bytes": 24 * 1024**3},
        metric={"host_physical_memory_bytes": str(32 * 1024**3)},
    )
    current_provenance["engine_continuation_predecessors"][0]["run_seal"]["sha256"] = "0" * 64
    with pytest.raises(COORDINATOR.CoordinatorError, match="lineage"):
        COORDINATOR._verify_engine_continuation_provenance(
            RUNTIME, tmp_path, materialized, current_directory.relative_to(tmp_path),
            locus=locus, fingerprint=fingerprint, provenance=current_provenance,
            measurement={"configured_envelope_bytes": 24 * 1024**3},
            metric={"host_physical_memory_bytes": str(32 * 1024**3)},
        )


def test_representative_order_is_burden_only_and_keeps_every_locus():
    rows = [
        {"locus_entry_id": "primary", "variant_count": "300", "priority_tier": "PRIMARY_TIER"},
        {"locus_entry_id": "secondary_small", "variant_count": "50", "priority_tier": "SECONDARY_TIER"},
        {"locus_entry_id": "secondary_large", "variant_count": "900", "priority_tier": "SECONDARY_TIER"},
        {"locus_entry_id": "control", "variant_count": "200", "priority_tier": "CONTROL_TIER"},
    ]
    ordered = RUNTIME.representative_execution_order(rows)
    assert ordered[:3] == ["secondary_small", "control", "secondary_large"]
    assert set(ordered) == {row["locus_entry_id"] for row in rows}


def test_coloc_strong_requires_same_pair_across_every_prior():
    validation, diagnostics = coloc_task_validation()
    values = {prior: (0.1, 0.84) for prior in (1e-6, 5e-6, 1e-5, 5e-5)}
    result = RUNTIME.classify_coloc_rows(
        probability_rows(values), diagnostics, validation,
        POLICY["required_science_outputs"]["output_12"]["schema"],
    )
    assert {row["classification"] for row in result} == {"STRONG_SHARED_SIGNAL"}
    assert {row["primary_shared_signal_rule_pass"] for row in result} == {"TRUE"}
    assert {row["same_signal_pair_prior_robust"] for row in result} == {"TRUE"}


def test_nonconstant_per_snp_n_forces_result_blind_coloc_claim_cap():
    validation, diagnostics = coloc_task_validation(n_status2=RUNTIME.N_DISPERSION_VARIABLE)
    values = {prior: (0.1, 0.84) for prior in (1e-6, 5e-6, 1e-5, 5e-5)}
    result = RUNTIME.classify_coloc_rows(
        probability_rows(values), diagnostics, validation,
        POLICY["required_science_outputs"]["output_12"]["schema"],
    )
    assert {row["classification"] for row in result} == {"INCONCLUSIVE"}
    assert {row["analysis_status"] for row in result} == {"TESTED"}
    assert {row["primary_shared_signal_rule_pass"] for row in result} == {"TRUE"}
    assert all("SAMPLE_SIZE_UNCERTAINTY" in row["fine_mapping_qc"] for row in result)
    assert {row["ld_qc"] for row in result} == {"PASS"}


def test_primary_strong_prior_sensitive_downgrades_to_moderate():
    validation, diagnostics = coloc_task_validation()
    values = {prior: (0.3, 0.6) for prior in (1e-6, 5e-6, 5e-5)}
    values[1e-5] = (0.1, 0.84)
    result = RUNTIME.classify_coloc_rows(
        probability_rows(values), diagnostics, validation,
        POLICY["required_science_outputs"]["output_12"]["schema"],
    )
    assert {row["classification"] for row in result} == {"MODERATE_SHARED_SIGNAL"}
    assert {row["primary_shared_signal_rule_pass"] for row in result} == {"TRUE"}
    assert {row["same_signal_pair_prior_robust"] for row in result} == {"FALSE"}


def test_rss_s_above_point_one_is_warning_not_invalid():
    validation, diagnostics = coloc_task_validation(rss1="0.1000000001")
    values = {prior: (0.1, 0.84) for prior in (1e-6, 5e-6, 1e-5, 5e-5)}
    result = RUNTIME.classify_coloc_rows(
        probability_rows(values), diagnostics, validation,
        POLICY["required_science_outputs"]["output_12"]["schema"],
    )
    assert {row["classification"] for row in result} == {"INCONCLUSIVE"}
    assert {row["analysis_status"] for row in result} == {"TESTED"}
    assert {row["primary_shared_signal_rule_pass"] for row in result} == {"TRUE"}
    assert all("LD_UNCERTAINTY_WARNING" in row["ld_qc"] for row in result)


def test_invalid_diagnostic_has_highest_classification_priority():
    validation, diagnostics = coloc_task_validation(status1="FAILED_KRIGING_QC")
    values = {prior: (0.1, 0.84) for prior in (1e-6, 5e-6, 1e-5, 5e-5)}
    result = RUNTIME.classify_coloc_rows(
        probability_rows(values), diagnostics, validation,
        POLICY["required_science_outputs"]["output_12"]["schema"],
    )
    assert {row["classification"] for row in result} == {"INVALID_INPUT"}
    assert {row["analysis_status"] for row in result} == {"INVALID_INPUT"}


def test_distinct_signals_rule_is_evaluated_after_shared_rules():
    validation, diagnostics = coloc_task_validation()
    values = {prior: (0.85, 0.1) for prior in (1e-6, 5e-6, 1e-5, 5e-5)}
    result = RUNTIME.classify_coloc_rows(
        probability_rows(values), diagnostics, validation,
        POLICY["required_science_outputs"]["output_12"]["schema"],
    )
    assert {row["classification"] for row in result} == {"DISTINCT_SIGNALS"}


def test_coloc_rejects_missing_prior_and_probability_drift():
    validation, diagnostics = coloc_task_validation()
    values = {prior: (0.1, 0.84) for prior in (1e-6, 5e-6, 1e-5, 5e-5)}
    rows = probability_rows(values)
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.classify_coloc_rows(
            rows[:-1], diagnostics, validation,
            POLICY["required_science_outputs"]["output_12"]["schema"],
        )
    rows = probability_rows(values)
    rows[0]["PP_H0"] = "0.5"
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.classify_coloc_rows(
            rows, diagnostics, validation,
            POLICY["required_science_outputs"]["output_12"]["schema"],
        )


def test_exact_r_model_prior_grid_and_no_fallback_are_source_locked():
    source = (ROOT / "scripts/150_run_track_b_susie_coloc.R").read_text()
    assert 'packageVersion("susieR")' in source and '!= "0.14.2"' in source
    assert 'packageVersion("coloc")' in source and '!= "5.2.3"' in source
    assert "coloc::runsusie(" in source and "coloc::coloc.susie(" in source
    assert "L = 10L" in source and "coverage = 0.95" in source
    assert "min_abs_corr = 0.5" in source and "maxit = 1000L" in source
    assert "repeat_until_convergence = FALSE" in source
    assert "estimate_residual_variance = FALSE" in source
    assert "prior_weights = prior" in source
    assert "c(1e-6, 5e-6, 1e-5, 5e-5)" in source
    assert "coloc::coloc.abf(" not in source and "coloc.abf(" not in source
    assert "conditional$logLR > 2 & abs(conditional$z) > 2" in source
    assert "diagnostic$s > 0.10" in source and "WARNING_RSS_LD_S_GT_0.10_NONFATAL" in source
    assert "verify_runtime_lock(\"PRE\"" in source
    assert "verify_runtime_lock(\"POST\", require_loaded_closure = TRUE)" in source
    assert "row$name %in% loadedNamespaces()" in source
    assert "getNamespaceInfo(asNamespace(row$name), \"path\")" in source
    assert "WARNING_VARIABLE_PER_SNP_N_MEDIAN_SCALAR_APPROXIMATION_CLAIM_CAPPED" in source
    status_function = source.split("diagnostic_status <- function", 1)[1].split("make_diagnostic_row <- function", 1)[0]
    assert "0.10" not in status_function


def test_python_and_r_engine_task_schemas_include_exact_ld_runtime_binding():
    source = (ROOT / "scripts/150_run_track_b_susie_coloc.R").read_text()
    task_block = source.split("task_fields <- c(", 1)[1].split("\n)", 1)[0]
    import re
    assert re.findall(r'"([A-Za-z0-9_]+)"', task_block) == RUNTIME.TASK_FIELDS
    for field in (
        "ld_runtime_package_lock_path", "ld_runtime_package_lock_sha256",
        "ld_runtime_attestation_path", "ld_runtime_attestation_sha256",
    ):
        assert field in RUNTIME.TASK_FIELDS


def test_runtime_lock_binds_rscript_and_exact_loaded_dependency_closure():
    content = RUNTIME.runtime_package_lock_bytes(ROOT).decode()
    rows = list(csv.DictReader(content.splitlines(), delimiter="\t"))
    assert list(rows[0]) == RUNTIME.RUNTIME_PACKAGE_LOCK_FIELDS
    assert rows[0]["record_type"] == "R_EXECUTABLE"
    assert rows[0]["path"] == ".r-env/bin/Rscript"
    package_rows = rows[1:]
    assert [row["name"] for row in package_rows] == list(
        RUNTIME.PINNED_RUNTIME_NAMESPACE_CLOSURE
    )
    assert {"susieR", "coloc", "Rcpp", "Matrix", "data.table"} <= {
        row["name"] for row in package_rows
    }
    assert all(int(row["file_count"]) > 0 and int(row["bytes"]) > 0 for row in rows)
    assert all(len(row["sha256"]) == 64 for row in rows)


def test_signed_ld_runtime_lock_binds_exact_lava_loaded_dependency_closure():
    content = RUNTIME.ld_runtime_package_lock_bytes(ROOT).decode()
    rows = list(csv.DictReader(content.splitlines(), delimiter="\t"))
    assert rows[0]["record_type"] == "R_EXECUTABLE"
    assert rows[0]["path"] == ".r-env/bin/Rscript"
    assert [row["name"] for row in rows[1:]] == list(
        RUNTIME.PINNED_LD_RUNTIME_NAMESPACE_CLOSURE
    )
    lava = [row for row in rows if row["name"] == "LAVA"]
    assert len(lava) == 1 and lava[0]["version"] == "0.1.5"
    assert lava[0]["path"] == ".r-env/lib/R/library/LAVA"
    assert int(lava[0]["file_count"]) > 0 and len(lava[0]["sha256"]) == 64


def test_signed_ld_extractor_ignores_adversarial_same_version_inherited_lava_tree(tmp_path):
    alternate_library = tmp_path / "alternate-library"
    shutil.copytree(ROOT / ".r-env/lib/R/library/LAVA", alternate_library / "LAVA")
    (alternate_library / "LAVA" / "SAME_VERSION_ALTERNATE_TREE").write_text("drift\n")
    runtime_lock = tmp_path / "ld-runtime.tsv"
    runtime_lock.write_bytes(RUNTIME.ld_runtime_package_lock_bytes(ROOT))
    runtime_lock_sha256 = hashlib.sha256(runtime_lock.read_bytes()).hexdigest()
    order = tmp_path / "order.tsv"
    order.write_text("SNP\tCHR\tBP\tA1\tA2\nrs1\t1\t1\tA\tC\n")
    environment = os.environ.copy()
    environment["R_LIBS"] = str(alternate_library)
    environment["R_LIBS_USER"] = str(alternate_library)
    completed = subprocess.run(
        [
            str(ROOT / ".r-env/bin/Rscript"), "--vanilla",
            str(ROOT / "scripts/149_extract_track_b_signed_ld.R"),
            str(tmp_path / "unused-reference"), "1", str(order),
            str(tmp_path / "ld.tsv.gz"), str(tmp_path / "ld.qc.tsv"), "10000",
            "LAVA_UKB_v1.1_EUR_GRCh37", "SIGNED_PEARSON_FROM_LAVA_READ_LD",
            str(runtime_lock), runtime_lock_sha256,
            str(tmp_path / "ld.runtime.attestation.tsv"), RUNTIME.ANALYSIS_ID,
            "A_LOC1", RUNTIME.EXECUTION_AMENDMENT_SHA256,
        ],
        cwd=ROOT, env=environment, capture_output=True, text=True, check=False,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    # Exact pinned-runtime PRE verification completed; the deliberately short
    # order then stops execution before any reference or scientific output use.
    assert "variant-order schema/count drifted" in output
    assert "runtime lock" not in output and "loaded-package path/version mismatch" not in output
    assert not (tmp_path / "ld.tsv.gz").exists()


def test_read_only_task_verifier_uses_attestation_not_spectral_recomputation():
    source = (ROOT / "scripts/148_track_b_finemapping_runtime.py").read_text()
    verifier = source.split("def deep_validate_locus_task", 1)[1].split(
        "def _description_version", 1
    )[0]
    assert "deep_validate_signed_ld(" not in verifier
    assert "_validate_ld_attestation(" in verifier
    assert "_validate_ld_runtime_attestation(" in verifier
    assert "runtime_package_lock" in verifier


def test_current_frozen_v1_lineage_files_remain_byte_exact():
    expected = {
        "config/track_b_finemapping_policy.json": "3f236d7427d6b5483dc62f5fb62a2041a679af5a00aaf9d4b4c6327b3864f374",
        "scripts/130_track_b_finemapping_contract.py": "473fddd7299ef87ac171c4c26923142df42840221eb5c77367a33175ab3799da",
        "results/track_b/finemapping/contract.lock.json": "273f8d5ab358e66e50dbd6619f0585b262efe4be52ee75eacd6de1b4aaa80366",
    }
    for relative, digest in expected.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest


def test_signed_ld_extractor_never_squares_or_emits_r2():
    source = (ROOT / "scripts/149_extract_track_b_signed_ld.R").read_text()
    assert "R <- loaded$ld" in source
    assert "SIGNED_PEARSON_CORRELATION" in source
    assert "R <- loaded$ld^2" not in source
    assert "r2" not in source.lower()


def test_safe_identity_rejects_symlink_and_hardlink(tmp_path):
    target = tmp_path / "target"
    target.write_text("data")
    (tmp_path / "link").symlink_to(target)
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.stable_identity(tmp_path, "link")
    os.link(target, tmp_path / "hard")
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.stable_identity(tmp_path, "target")


def test_stable_identity_detects_named_file_toctou_replacement(tmp_path, monkeypatch):
    target = tmp_path / "target"
    target.write_bytes(b"a" * (4 * 1024 * 1024 + 1))
    replacement = tmp_path / "replacement"
    replacement.write_bytes(b"b" * target.stat().st_size)
    original_read = RUNTIME.os.read
    replaced = False

    def replacing_read(descriptor, count):
        nonlocal replaced
        block = original_read(descriptor, count)
        if block and not replaced:
            replaced = True
            os.replace(replacement, target)
        return block

    monkeypatch.setattr(RUNTIME.os, "read", replacing_read)
    with pytest.raises(RUNTIME.FineMappingError, match="changed while hashing"):
        RUNTIME.stable_identity(tmp_path, "target")


def test_no_replace_publication_refuses_second_content(tmp_path):
    first = RUNTIME.publish_bytes_no_replace(tmp_path, "nested/value.tsv", b"a\n")
    assert first["sha256"] == hashlib.sha256(b"a\n").hexdigest()
    with pytest.raises(RUNTIME.FineMappingError):
        RUNTIME.publish_bytes_no_replace(tmp_path, "nested/value.tsv", b"b\n")


def test_parent_plus_pgid_sampler_deduplicates_and_keeps_reparented_descendants(monkeypatch):
    snapshot = subprocess.CompletedProcess(
        args=[], returncode=0,
        stdout="100 50 11\n200 200 20\n201 200 30\n202 999 999\n",
        stderr="",
    )
    monkeypatch.setattr(COORDINATOR.subprocess, "run", lambda *args, **kwargs: snapshot)
    # PID 201 has no PPID relationship encoded here; PGID membership alone
    # retains it, and the parent PID 100 is outside the worker PGID.
    assert COORDINATOR.process_group_rss_bytes(200, coordinator_process_id=100) == 61 * 1024


def test_parent_in_owned_group_is_not_double_counted(monkeypatch):
    snapshot = subprocess.CompletedProcess(
        args=[], returncode=0, stdout="100 200 11\n201 200 30\n", stderr="",
    )
    monkeypatch.setattr(COORDINATOR.subprocess, "run", lambda *args, **kwargs: snapshot)
    assert COORDINATOR.process_group_rss_bytes(200, coordinator_process_id=100) == 41 * 1024


def test_rusage_failure_fallback_units_are_platform_normalized(monkeypatch):
    monkeypatch.setattr(COORDINATOR.sys, "platform", "linux")
    assert COORDINATOR._normalize_ru_maxrss(123) == 123 * 1024
    monkeypatch.setattr(COORDINATOR.sys, "platform", "darwin")
    assert COORDINATOR._normalize_ru_maxrss(123) == 123


def test_sigkill_process_group_cleanup_is_confirmed():
    process = subprocess.Popen(
        [
            sys.executable, "-u", "-c",
            "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
            "print('READY', flush=True); time.sleep(30)",
        ],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,
    )
    try:
        assert process.stdout is not None and process.stdout.readline().strip() == "READY"
        existed, confirmed = COORDINATOR._terminate_process_group(process)
        assert existed is True and confirmed is True
        assert process.wait(timeout=2) == -9
        with pytest.raises(ProcessLookupError):
            os.killpg(process.pid, 0)
    finally:
        if process.poll() is None:
            os.killpg(process.pid, 9)
            process.wait(timeout=2)
        if process.stdout is not None:
            process.stdout.close()
        if process.stderr is not None:
            process.stderr.close()


def test_interrupted_stage_is_retained_with_reconciliation_receipt(tmp_path):
    stage = tmp_path / COORDINATOR.STAGE_ROOT_REL / "A_LOC1" / "attempt.building"
    stage.mkdir(parents=True)
    (stage / "partial.tsv").write_text("partial\n")
    receipts = COORDINATOR._quarantine_interrupted_run_stages(
        RUNTIME, tmp_path, "A_LOC1"
    )
    assert len(receipts) == 1
    assert not stage.exists()
    quarantines = list((tmp_path / COORDINATOR.QUARANTINE_ROOT_REL / "interrupted_r_engine_attempt").iterdir())
    assert len(quarantines) == 1
    assert (quarantines[0] / "partial.tsv").read_text() == "partial\n"
    assert (quarantines[0] / "quarantine.receipt.json").is_file()
    assert (quarantines[0] / "quarantine.evidence.json").is_file()
    assert (quarantines[0] / "quarantine.seal.json").is_file()
    verified = COORDINATOR.verify_quarantine_directory(
        RUNTIME, tmp_path, quarantines[0],
        category="interrupted_r_engine_attempt", locus_entry_id="A_LOC1",
    )
    assert verified["evidence_identity"]["sha256"]


def test_failed_and_orphan_materialization_input_stages_are_sealed_without_deletion(tmp_path):
    input_root = tmp_path / RUNTIME.INPUT_ROOT_REL
    stage = input_root / ".A_LOC1.1234.failed.building"
    stage.mkdir(parents=True)
    (stage / "signed_ld.qc.tsv").write_text("diagnostic\n")
    quarantines = COORDINATOR._quarantine_materialization_input_stages(
        RUNTIME, tmp_path, "A_LOC1", reason="TEST_CAUGHT_MATERIALIZATION_FAILURE",
    )
    assert not stage.exists()
    assert len(quarantines) == 1
    retained = Path(quarantines[0]["receipt"]["quarantine_path"])
    retained_directory = tmp_path / retained
    assert (retained_directory / "signed_ld.qc.tsv").read_text() == "diagnostic\n"
    assert quarantines[0]["receipt"]["scientific_output_deleted"] is False
    assert quarantines[0]["receipt"]["recovered_after_interruption"] is False

    # Simulate a process death immediately after the atomic move and before its
    # metadata was sealed; the next pair-lease reconciliation finishes it.
    crash_directory = (
        tmp_path / COORDINATOR.QUARANTINE_ROOT_REL /
        "materialization_input_attempt" / "A_LOC1.crash"
    )
    crash_directory.mkdir(parents=True)
    (crash_directory / "partial.tsv").write_text("crash evidence\n")
    reconciled = COORDINATOR._quarantine_materialization_input_stages(
        RUNTIME, tmp_path, "A_LOC1", reason="TEST_RESTART",
    )
    assert len(reconciled) == 2
    recovered = [item for item in reconciled if item["receipt"]["recovered_after_interruption"]]
    assert len(recovered) == 1
    assert (crash_directory / "partial.tsv").read_text() == "crash evidence\n"
    assert (crash_directory / "quarantine.seal.json").is_file()


def test_materializer_source_never_unlinks_failed_attempt_stage():
    source = (ROOT / "scripts/148_track_b_finemapping_runtime.py").read_text()
    failure_boundary = source.split("return verify_materialized_bundle", 1)[1].split(
        "def _file_sha256", 1
    )[0]
    assert "_remove_owned_empty_tree" not in source
    assert ".unlink(" not in failure_boundary and ".rmdir(" not in failure_boundary
    assert "Never erase a failed materialization" in failure_boundary


def test_stream_publication_is_no_replace_and_detects_conflict(tmp_path):
    source = tmp_path / "source"
    source.write_bytes(b"complete\n")
    expected = {"bytes": source.stat().st_size, "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    first = COORDINATOR._stream_publish_no_replace(
        RUNTIME, tmp_path, source, Path("out/value"), expected,
    )
    second = COORDINATOR._stream_publish_no_replace(
        RUNTIME, tmp_path, source, Path("out/value"), expected,
    )
    assert first["sha256"] == second["sha256"]
    source2 = tmp_path / "source2"
    source2.write_bytes(b"different\n")
    expected2 = {"bytes": source2.stat().st_size, "sha256": hashlib.sha256(source2.read_bytes()).hexdigest()}
    with pytest.raises(COORDINATOR.CoordinatorError):
        COORDINATOR._stream_publish_no_replace(
            RUNTIME, tmp_path, source2, Path("out/value"), expected2,
        )


def test_component_ram_schema_collision_key_and_unavailable_retention():
    family = {
        "rows": [{"locus_entry_id": "A_LOC1", "pair_id": "A", "CHR": "1"}],
        "unavailable_rows": [{"unavailable_entry_id": "B_UNMAPPED_1", "pair_id": "B", "CHR": "2"}],
        "unavailable_identity": {"sha256": "b" * 64},
    }
    outcome = {
        "terminal_status": "COMPLETE_MATERIALIZED",
        "receipt": {
            "peak_rss_bytes": 1024**3, "elapsed_monotonic_seconds": 2.0,
            "bundle_lock": {"sha256": "a" * 64},
        },
        "bundle": {"manifest_row": {"variant_count": "70"}},
        "seal_identity": {"sha256": "c" * 64},
    }
    rows = COORDINATOR.build_component_ram_rows(RUNTIME, family, {"A_LOC1": outcome})
    assert list(rows[0]) == COORDINATOR.COMPONENT_RAM_FIELDS
    assert rows[0]["peak_ram_gb"] == "1"
    assert rows[1]["exit_status"] == "BLOCKED_BY_DATA:NO_PREDECLARED_SIGNED_LD_BLOCK"
    assert len({(r["analysis"], r["pair"], r["locus"], r["chromosome"]) for r in rows}) == 2


def test_canonical_policy_schemas_are_used_exactly():
    assert RUNTIME._materialized_manifest_fields(POLICY) == POLICY["locus_entry_contract"]["locus_manifest_schema"]
    assert COORDINATOR.COMPONENT_RAM_FIELDS == [
        "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
        "runtime_sec", "exit_status", "output_hash",
    ]
    assert RUNTIME.RAW_VARIANT_FIELDS != POLICY["required_science_outputs"]["output_10"]["schema"]
    assert POLICY["required_science_outputs"]["output_10"]["schema"][-3:] == [
        "fine_mapping_classifications", "analysis_status", "claim_limit"
    ]


def test_zero_family_verifier_forbids_science_placeholders(tmp_path):
    class Stub:
        ZERO_STATE = RUNTIME.ZERO_STATE

        @staticmethod
        def safe_path(root, value, label):
            del label
            return root / value

    materialized = {"state": RUNTIME.ZERO_STATE, "policy": POLICY}
    result = COORDINATOR.validate_canonical_family(
        tmp_path, runtime=Stub(), materialized=materialized,
    )
    assert result["state"] == RUNTIME.ZERO_STATE
    forbidden = tmp_path / POLICY["required_science_outputs"]["output_10"]["path"]
    forbidden.parent.mkdir(parents=True)
    forbidden.write_text("header\n")
    with pytest.raises(COORDINATOR.CoordinatorError):
        COORDINATOR.validate_canonical_family(tmp_path, runtime=Stub(), materialized=materialized)


def test_source_binds_failure_retention_fresh_process_and_no_parameter_retry():
    py = (ROOT / "scripts/151_run_track_b_finemapping_sequential.py").read_text()
    r = (ROOT / "scripts/150_run_track_b_susie_coloc.R").read_text()
    assert "start_new_session=True" in py
    assert "PARENT_PLUS_PGID_AGGREGATE_PS_RSS_BYTES" in py
    assert '"FAILED_RESOURCE_OOM"' in py and '"BLOCKED_BY_COMPUTE"' in py
    assert "model_or_prior_retry\": False" in py
    assert "any_coloc_failure" in r and "COLOC_SUSIE_FAILURE" in r
    assert "repeat_until_convergence = FALSE" in r
    assert '"TRACK_B_FINEMAP_MEASURED_PHASE": "MATERIALIZATION"' in py
    assert "verify_failed_materialization_attempt" in py
    assert "current_physical <= max(" in py
    assert "current_effective <= max(" in py
    assert "engine_continuation_predecessors" in py
    assert "_bound_engine_predecessors(" in py
    assert "R-engine larger-host continuation requires a retained resource failure" in py
    assert "prelaunch_quarantine_seals" in py and "attempt_quarantine_seals" in py
    execute_family = py.split("def execute_family", 1)[1]
    assert "configured_envelope_bytes=configured_envelope_bytes" in execute_family
    assert "larger_host_continuation=larger_host_continuation" in execute_family


def test_scripts_parse_and_python_compile():
    subprocess.run(
        [sys.executable, "-m", "py_compile", str(ROOT / "scripts/148_track_b_finemapping_runtime.py"),
         str(ROOT / "scripts/151_run_track_b_finemapping_sequential.py")],
        check=True,
    )
    subprocess.run(
        ["Rscript", "-e", (
            "invisible(parse(file='scripts/149_extract_track_b_signed_ld.R'));"
            "invisible(parse(file='scripts/150_run_track_b_susie_coloc.R'))"
        )], cwd=ROOT, check=True,
    )
