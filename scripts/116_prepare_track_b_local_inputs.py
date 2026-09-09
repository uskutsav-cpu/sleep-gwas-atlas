#!/usr/bin/env python3
"""Prepare and freeze the three-pair Track B local-analysis input family."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np


POLICY = Path("config/track_b_local_analysis_policy.json")
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
PANEL = Path("config/analysis_panel.tsv")
INTERCEPTS = Path("results/tables/ldsc_intercept_45x45.tsv")
LAVA_LOCUS = Path("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
INPUT_INFO = Path("results/track_b/lava_input_info.tsv")
PAIR_INPUT = Path("results/track_b/lava_pair_manifest.tsv")
OVERLAP = Path("results/track_b/lava_sample_overlap.txt")
PROVENANCE = Path("results/track_b/lava_input_provenance.tsv")
CONDITIONAL = Path("results/track_b/local_conditional_manifest.tsv")
ROBUSTNESS = Path("results/track_b/local_method_robustness_plan.tsv")
RUNTIME_POLICY = Path("results/track_b/lava_runtime_policy.tsv")
HANDOFF = Path("results/track_b/LOCAL_ANALYSIS_HANDOFF.md")
LOCK = Path("results/track_b/local_analysis_input.lock.json")
REFERENCE_PROVENANCE = Path("ref/lava/ukb_v1.1/reference.provenance.json")
RAM_BENCHMARK = Path("results/track_b/RAM_BENCHMARK.tsv")
RAM_BENCHMARK_PROVENANCE = Path("results/track_b/RAM_BENCHMARK.provenance.json")
CHROMOSOME_INPUT_LOCK = Path("results/track_b/lava_chromosome_inputs.provenance.json")
LAVA_CONTRACT = Path("scripts/119_track_b_lava_contract.py")
CHECKPOINT_DIR = Path("results/track_b/checkpoints/lava")
STAGING_ROOT = Path("results/track_b/staging")
RESULT_EVIDENCE = [
    Path("results/track_b/local/lava_results.provenance.json"),
    Path("results/track_b/local/lava_locus_status.tsv"),
    Path("results/track_b/local/lava_univariate.tsv"),
    Path("results/track_b/local/lava_bivariate.tsv"),
    Path("results/track_b/local/lava_conditional.tsv"),
    Path("results/track_b/04_lava_local_results.tsv"),
    Path("results/track_b/05_local_conditional_results.tsv"),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def load_policy() -> dict[str, object]:
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    if (
        policy["lava_version"] != "0.1.5"
        or policy["expected_loci"] != 2495
        or policy["expected_discovery_pairs"] != 3
        or policy["expected_analysis_traits"] != 8
        or policy["planned_univariate_tests"] != 2495 * 8
        or policy["planned_bivariate_pair_locus_family_max"] != 2495 * 3
        or policy.get("ram_aware_execution", {}).get("execution_unit") != "WHOLE_PREDECLARED_LAVA_LOCUS"
        or policy.get("ram_aware_execution", {}).get("maximum_loci_per_process") != 1
        or policy.get("locus_definition", {}).get("path") != str(LAVA_LOCUS)
        or policy.get("locus_definition", {}).get("sha256") != "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882"
    ):
        raise SystemExit("ERROR: Track B local policy scope drifted")
    expected = policy["univariate_alpha"] / policy["planned_univariate_tests"]
    if not math.isclose(expected, policy["univariate_p_threshold"], rel_tol=1e-12):
        raise SystemExit("ERROR: Track B local univariate threshold is inconsistent")
    model_ids: list[str] = []
    for pair_id in policy["pair_order"]:
        models = policy["conditional_models"][pair_id]
        observed_covariates = [
            covariate for model in models for covariate in model["covariates"]
        ]
        if observed_covariates != policy["conditional_covariates"][pair_id]:
            raise SystemExit(f"ERROR: conditional model/covariate drift for {pair_id}")
        for model in models:
            model_ids.append(str(model["conditional_model_id"]))
            if not model["covariates"] or not str(model["rationale"]).strip():
                raise SystemExit(f"ERROR: incomplete conditional model for {pair_id}")
    if len(model_ids) != len(set(model_ids)):
        raise SystemExit("ERROR: duplicate Track B conditional model ID")
    return policy


def validate_locus_file(policy: dict[str, object]) -> None:
    if not LAVA_LOCUS.is_file():
        raise SystemExit("ERROR: LAVA 2,495-locus definition is missing")
    with LAVA_LOCUS.open(encoding="utf-8") as handle:
        lines = sum(1 for _ in handle)
    if lines != int(policy["expected_loci"]) + 1:
        raise SystemExit(f"ERROR: expected 2,495 LAVA loci, found {lines-1}")
    if sha256(LAVA_LOCUS) != policy["locus_definition"]["sha256"]:
        raise SystemExit("ERROR: LAVA locus definition differs from the policy-pinned SHA-256")


def refuse_post_result_refreeze() -> None:
    evidence = [path for path in RESULT_EVIDENCE if path.exists()]
    if CHECKPOINT_DIR.is_dir():
        evidence.extend(sorted(CHECKPOINT_DIR.glob("*/discovery/locus_*"))[:1])
    if STAGING_ROOT.is_dir():
        evidence.extend(sorted(STAGING_ROOT.glob("lava_*"))[:1])
    if evidence:
        raise SystemExit(
            "ERROR: refusing to refreeze Track B local inputs after checkpoint/staging/result evidence exists: "
            + ", ".join(str(path) for path in evidence[:3])
        )


def validate_sumstats(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing LAVA input: {path}")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    required = {"SNP", "A1", "A2", "Z", "N"}
    if not required.issubset(fields):
        raise SystemExit(f"ERROR: {path} lacks LAVA fields: {sorted(required-set(fields))}")
    return ",".join(fields)


def input_rows(policy: dict[str, object], panel: dict[str, dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    info: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    for trait in policy["trait_order"]:
        row = panel[str(trait)]
        path = Path("data/munged") / f"{trait}.sumstats.gz"
        schema = validate_sumstats(path)
        if row["type"] == "binary":
            cases, controls, prevalence = row["ncase"], row["ncontrol"], row["pop_prev"]
            if not cases.isdigit() or not controls.isdigit() or not 0 < float(prevalence) < 1:
                raise SystemExit(f"ERROR: unresolved binary metadata for {trait}")
        else:
            cases = controls = prevalence = "NA"
        info.append({
            "phenotype": trait, "cases": cases, "controls": controls,
            "prevalence": prevalence, "filename": str(path),
        })
        provenance.append({
            "phenotype": trait, "filename": str(path), "bytes": path.stat().st_size,
            "sha256": sha256(path), "observed_columns": schema,
            "validation_status": "VALIDATED_FOR_LAVA_LOCAL_NOT_FINE_MAPPING",
        })
    return info, provenance


def pair_rows(policy: dict[str, object]) -> list[dict[str, object]]:
    manifest = read_tsv(PAIR_MANIFEST)
    by_id = {row["pair_id"]: row for row in manifest}
    if set(by_id) != {"A", "B", "CONTROL"}:
        raise SystemExit("ERROR: frozen Track B pair family drifted")
    rows = []
    for pair_id in policy["pair_order"]:
        row = by_id[str(pair_id)]
        rows.append({
            "pair_id": pair_id,
            "trait1": row["sleep_trait"],
            "trait2": row["external_trait"],
            "discovery_rg": row["discovery_rg"],
            "discovery_SE": row["SE"],
            "discovery_P": row["P"],
            "discovery_FDR": row["FDR"],
            "planned_loci": policy["expected_loci"],
            "bivariate_family": policy["bivariate_multiple_testing"],
            "pair_replacement": "FORBIDDEN",
        })
    return rows


def overlap_text(policy: dict[str, object]) -> tuple[str, dict[str, float]]:
    all_rows = read_tsv(INTERCEPTS)
    all_traits = [row["trait_id"] for row in all_rows]
    if len(all_rows) != 45 or list(all_rows[0].keys())[1:] != all_traits:
        raise SystemExit("ERROR: canonical 45-trait intercept matrix drifted")
    index = {trait: row for trait, row in zip(all_traits, all_rows)}
    traits = [str(value) for value in policy["trait_order"]]
    matrix = np.array([[float(index[first][second]) for second in traits] for first in traits])
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, atol=1e-10):
        raise SystemExit("ERROR: selected LDSC intercept matrix is invalid")
    diagonal = np.diag(matrix)
    correlation = matrix / np.sqrt(np.outer(diagonal, diagonal))
    correlation = np.round((correlation + correlation.T) / 2, 5)
    np.fill_diagonal(correlation, 1.0)
    eigenvalues = np.linalg.eigvalsh(correlation)
    if eigenvalues[0] <= 0:
        raise SystemExit("ERROR: Track B overlap matrix is not positive definite")
    buffer = io.StringIO()
    buffer.write(" " + " ".join(traits) + "\n")
    for trait, values in zip(traits, correlation):
        buffer.write(trait + " " + " ".join(f"{value:.5f}" for value in values) + "\n")
    return buffer.getvalue(), {
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "condition_number": float(eigenvalues[-1] / eigenvalues[0]),
    }


def conditional_rows(policy: dict[str, object]) -> list[dict[str, object]]:
    pairs = {row["pair_id"]: row for row in read_tsv(PAIR_MANIFEST)}
    output = []
    for pair_id in policy["pair_order"]:
        models = policy["conditional_models"][pair_id]
        if not models:
            output.append({
                "pair_id": pair_id, "trait1": pairs[pair_id]["sleep_trait"],
                "trait2": pairs[pair_id]["external_trait"], "conditional_model_id": "NONE",
                "covariates": "NONE",
                "rationale": policy["conditional_covariate_rationale"][pair_id],
                "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS", "execution_status": "NOT_PLANNED_FOR_CONTROL",
            })
            continue
        for model in models:
            output.append({
                "pair_id": pair_id, "trait1": pairs[pair_id]["sleep_trait"],
                "trait2": pairs[pair_id]["external_trait"],
                "conditional_model_id": model["conditional_model_id"],
                "covariates": ";".join(model["covariates"]), "rationale": model["rationale"],
                "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS", "execution_status": "PENDING_LAVA_REFERENCE",
            })
    return output


def robustness_rows() -> list[dict[str, object]]:
    return [
        {
            "method": "HDL-L", "status": "BLOCKED_BY_SOFTWARE_AND_REFERENCE",
            "reference": "NO_PINNED_ANCESTRY_MATCHED_HDL_L_REFERENCE_PRESENT",
            "comparison_rule": "Compare locus overlap, sign, approximate magnitude and corrected significance; preserve discordance.",
            "allowed_labels": "CONCORDANT;PARTIAL_CONCORDANCE;DISCORDANT",
        },
        {
            "method": "rho-HESS", "status": "FALLBACK_REQUIRES_SEPARATE_PINNED_REFERENCE_AND_CONTRACT",
            "reference": "NO_TRACK_B_SPECIFIC_SEALED_RHO_HESS_REFERENCE_PRESENT",
            "comparison_rule": "May be used only as an independently implemented local-correlation sensitivity, not as a proxy for LAVA.",
            "allowed_labels": "CONCORDANT;PARTIAL_CONCORDANCE;DISCORDANT",
        },
    ]


def runtime_rows(policy: dict[str, object]) -> list[dict[str, object]]:
    runtime = policy["runtime"]
    ram = policy["ram_aware_execution"]
    values = {
        "analysis_id": policy["analysis_id"], "lava_version": policy["lava_version"],
        "reference_prefix": policy["reference_prefix"], "expected_loci": policy["expected_loci"],
        "expected_traits": policy["expected_analysis_traits"], "expected_pairs": policy["expected_discovery_pairs"],
        "planned_univariate_tests": policy["planned_univariate_tests"],
        "univariate_p_threshold": policy["univariate_p_threshold"],
        "bivariate_fdr_alpha": policy["bivariate_fdr_alpha"],
        "conditional_fdr_alpha": policy["bivariate_fdr_alpha"],
        "maximum_locus_failure_fraction": policy["maximum_locus_failure_fraction"],
        "maximum_univariate_untested_fraction": policy["maximum_univariate_untested_fraction"],
        "maximum_bivariate_failure_fraction": policy["maximum_bivariate_failure_fraction"],
        "maximum_conditional_failure_fraction": policy["maximum_conditional_failure_fraction"],
        "ram_execution_unit": ram["execution_unit"],
        "maximum_loci_per_process": ram["maximum_loci_per_process"],
        "worker_process_rule": ram["worker_process_rule"],
        "reference_loading_rule": ram["reference_loading_rule"],
        "family_correction_rule": ram["family_correction_rule"],
        "memory_decision_rule": ram["memory_decision_rule"],
        "memory_safety_reserve_bytes": ram["memory_safety_reserve_bytes"],
        "memory_admission_rule": ram["memory_admission_rule"],
        "resume_rule": ram["resume_rule"],
        "random_seed": policy["random_seed"], "min_K": runtime["min_K"],
        "prune_threshold": runtime["prune_threshold"], "max_proportion_K": runtime["max_proportion_K"],
        "max_block_size": runtime["max_block_size"], "cap_estimates": str(runtime["cap_estimates"]).lower(),
        "conditional_max_r2": runtime["conditional_max_r2"],
        "conditional_attenuation_fraction_partial": policy["interpretation_policy"]["conditional_attenuation_fraction_partial"],
        "conditional_attenuation_fraction_full": policy["interpretation_policy"]["conditional_attenuation_fraction_full"],
        "conditional_execution_gate": policy["conditional_execution_gate"],
        "conditional_multiple_testing": policy["conditional_multiple_testing"],
    }
    return [{"key": key, "value": value} for key, value in values.items()]


@lru_cache(maxsize=1)
def current_lava_fingerprint() -> str | None:
    if not CHROMOSOME_INPUT_LOCK.is_file() or not REFERENCE_PROVENANCE.is_file():
        return None
    result = subprocess.run(
        [sys.executable, str(LAVA_CONTRACT), "--execution-fingerprint"],
        text=True, capture_output=True, check=False,
    )
    value = result.stdout.strip()
    return value if result.returncode == 0 and len(value) == 64 else None


def ram_benchmark_state(
    physical_memory_bytes: int, safety_reserve_bytes: int,
) -> tuple[str, float | None, int]:
    """Return an operational state without inventing an a-priori RAM minimum."""
    if (
        not RAM_BENCHMARK.is_file() or RAM_BENCHMARK.stat().st_size == 0
        or not RAM_BENCHMARK_PROVENANCE.is_file()
    ):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    try:
        provenance = json.loads(RAM_BENCHMARK_PROVENANCE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    fingerprint = current_lava_fingerprint()
    if (
        provenance.get("schema_version") != 2
        or fingerprint is None
        or provenance.get("execution_fingerprint") != fingerprint
        or provenance.get("benchmark_sha256") != sha256(RAM_BENCHMARK)
        or provenance.get("memory_safety_reserve_bytes") != safety_reserve_bytes
    ):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    selected = provenance.get("selected_representatives", {})
    discovery = selected.get("discovery", {}) if isinstance(selected, dict) else {}
    conditional = selected.get("conditional", {}) if isinstance(selected, dict) else {}
    if set(discovery) != {"SMALL", "MEDIAN", "LARGE"}:
        return "RAM_BENCHMARK_REQUIRED", None, len(discovery)
    evidence = provenance.get("measured_admission_units", [])
    attempt_rows = provenance.get("all_lava_attempt_rows", [])
    if not isinstance(evidence, list) or not isinstance(attempt_rows, list):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    try:
        exact_peaks = [int(item["peak_rss_bytes"]) for item in evidence]
        exact_peaks.extend(int(round(float(item["peak_ram_gb"]) * 1024**3)) for item in attempt_rows)
        peaks = [value / 1024**3 for value in exact_peaks]
    except (KeyError, TypeError, ValueError):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    if not peaks or any(not math.isfinite(peak) or peak <= 0 for peak in peaks):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    if provenance.get("maximum_measured_admission_peak_rss_bytes") != max(exact_peaks):
        return "RAM_BENCHMARK_REQUIRED", None, 0
    if provenance.get("state") == "DISCOVERY_REPRESENTATIVES_COMPLETE_CONDITIONAL_PENDING":
        return "CONDITIONAL_RAM_BENCHMARK_PENDING", max(peaks), len(attempt_rows)
    if provenance.get("state") != "ALL_REQUIRED_PHASE_REPRESENTATIVES_COMPLETE":
        return "RAM_BENCHMARK_REQUIRED", max(peaks), len(attempt_rows)
    maximum = max(peaks)
    state = (
        "MEASURED_SEQUENTIAL_RAM_PASS"
        if maximum * 1024**3 + safety_reserve_bytes <= physical_memory_bytes
        else "BLOCKED_BY_MEASURED_PER_LOCUS_RAM"
    )
    return state, maximum, len(attempt_rows)


def handoff_text(
    policy: dict[str, object], free_bytes: int, physical_memory_bytes: int,
    reference_ready: bool,
) -> str:
    acquisition_storage_ready = free_bytes >= int(policy["reference_minimum_free_bytes"])
    safety_reserve = int(policy["ram_aware_execution"]["memory_safety_reserve_bytes"])
    benchmark_state, measured_peak_gib, measured_loci = ram_benchmark_state(
        physical_memory_bytes, safety_reserve,
    )
    data_status = "READY" if reference_ready else "BLOCKED_BY_DATA"
    if not reference_ready and not acquisition_storage_ready:
        compute_status = "BLOCKED_BY_COMPUTE_STORAGE"
    elif benchmark_state == "BLOCKED_BY_MEASURED_PER_LOCUS_RAM":
        compute_status = benchmark_state
    elif benchmark_state == "MEASURED_SEQUENTIAL_RAM_PASS":
        compute_status = "READY_FOR_SEQUENTIAL_EXECUTION"
    elif benchmark_state == "CONDITIONAL_RAM_BENCHMARK_PENDING":
        compute_status = "READY_FOR_SEQUENTIAL_DISCOVERY_CONDITIONAL_BENCHMARK_PENDING"
    else:
        compute_status = "READY_FOR_RAM_BENCHMARK"
    reference_sentence = (
        "The official LAVA UK Biobank European v1.1 LD payload is present and sealed."
        if reference_ready else
        "The official LAVA UK Biobank European v1.1 LD payload is absent."
    )
    return f"""# Track B local-analysis production handoff

Current data status: `{data_status}`. Current compute status: `{compute_status}`.

The three frozen pairs, eight required analysis traits, 2,495 loci, LDSC overlap submatrix, conditional covariates, correction families, and input hashes are prepared and locked. {reference_sentence}

Current free space at preparation: {free_bytes} bytes ({free_bytes/1024**3:.3f} GiB).
Policy minimum for checksum-ledgered acquisition and extraction: {policy['reference_minimum_free_bytes']} bytes (35 GiB).
Physical memory at preparation: {physical_memory_bytes} bytes ({physical_memory_bytes/1024**3:.3f} GiB).
Fingerprint-bound RAM evidence: `{benchmark_state}` ({measured_loci} selected intact phase representatives; maximum observed peak {measured_peak_gib if measured_peak_gib is not None else 'NA'} GiB).
No assumed monolithic RAM minimum is used for the production decision; admission preserves a {safety_reserve} byte (1 GiB) non-worker reserve.

## Minimum production environment

- x86_64 Linux recommended
- R 4.3.x with LAVA 0.1.5
- enough RAM for the largest measured whole-locus worker; execute one fresh R process per locus
- at least 35 GiB free for reference acquisition; 60 GiB recommended for outputs/checkpoints
- the exact repository commit and ignored dense/munged inputs whose hashes are in `local_analysis_input.lock.json`

## Exact execution order

```bash
python3 scripts/110_track_b_checkpoint.py --verify
python3 scripts/111_freeze_track_b_pairs.py --verify
python3 scripts/113_freeze_track_b_replication_sources.py --verify
python3 scripts/115_build_track_b_dense_qc.py --verify
python3 scripts/116_prepare_track_b_local_inputs.py --verify
bash scripts/32_download_lava_reference.sh --download
python3 scripts/lava_contract.py --verify-reference
python3 scripts/129_prepare_track_b_lava_chromosome_inputs.py
python3 scripts/129_prepare_track_b_lava_chromosome_inputs.py --verify
python3 scripts/119_track_b_lava_contract.py --preflight
python3 scripts/131_run_track_b_lava_sequential.py --benchmark
python3 scripts/131_run_track_b_lava_sequential.py --run
python3 scripts/121_validate_track_b_lava.py
```

Do not run the existing 45-trait/396-pair LAVA result family as a substitute for Track B. The dedicated Track B supervisor consumes `results/track_b/lava_pair_manifest.tsv`, uses all 2,495 complete loci, runs each locus in a fresh R process, and applies BH FDR only after collating the full frozen family. Conditional tests run in a second per-locus pass after the discovery-family BH gate and use only the separately frozen BMI-only, sleep-apnea-only, and MDD-only models in `local_conditional_manifest.tsv`.

No local result exists yet. No empty table or synthetic output is presented as science.
"""


def build() -> dict[Path, str]:
    policy = load_policy()
    validate_locus_file(policy)
    panel = {row["trait_id"]: row for row in read_tsv(PANEL)}
    if [trait for trait in policy["trait_order"] if trait not in panel]:
        raise SystemExit("ERROR: Track B local trait absent from locked panel")
    info, provenance = input_rows(policy, panel)
    pairs = pair_rows(policy)
    overlap, overlap_diagnostics = overlap_text(policy)
    conditional = conditional_rows(policy)
    robustness = robustness_rows()
    free_bytes = shutil.disk_usage(Path(".")).free
    physical_memory_bytes = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    reference_ready = REFERENCE_PROVENANCE.is_file()
    safety_reserve = int(policy["ram_aware_execution"]["memory_safety_reserve_bytes"])
    benchmark_state, measured_peak_gib, measured_loci = ram_benchmark_state(
        physical_memory_bytes, safety_reserve,
    )
    outputs = {
        INPUT_INFO: tsv_text(["phenotype", "cases", "controls", "prevalence", "filename"], info),
        PAIR_INPUT: tsv_text([
            "pair_id", "trait1", "trait2", "discovery_rg", "discovery_SE", "discovery_P",
            "discovery_FDR", "planned_loci", "bivariate_family", "pair_replacement",
        ], pairs),
        OVERLAP: overlap,
        PROVENANCE: tsv_text([
            "phenotype", "filename", "bytes", "sha256", "observed_columns", "validation_status",
        ], provenance),
        CONDITIONAL: tsv_text([
            "pair_id", "trait1", "trait2", "conditional_model_id", "covariates",
            "rationale", "selection_timing", "execution_status",
        ], conditional),
        ROBUSTNESS: tsv_text(["method", "status", "reference", "comparison_rule", "allowed_labels"], robustness),
        RUNTIME_POLICY: tsv_text(["key", "value"], runtime_rows(policy)),
        HANDOFF: handoff_text(policy, free_bytes, physical_memory_bytes, reference_ready),
    }
    immutable_hashes = {str(path): hashlib.sha256(value.encode()).hexdigest() for path, value in outputs.items() if path != HANDOFF}
    lock = {
        "schema_version": 1,
        "analysis_id": policy["analysis_id"],
        "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS",
        "local_results_accessed_before_input_freeze": False,
        "policy_sha256": sha256(POLICY),
        "pair_manifest_sha256": sha256(PAIR_MANIFEST),
        "locus_file_sha256": sha256(LAVA_LOCUS),
        "input_artifact_sha256": immutable_hashes,
        "overlap_diagnostics": overlap_diagnostics,
        "reference_status": "READY" if reference_ready else "BLOCKED_BY_DATA",
        "required_reference_provenance": str(REFERENCE_PROVENANCE),
        "reference_minimum_free_bytes": policy["reference_minimum_free_bytes"],
        "current_compute_snapshot": {
            "free_bytes": free_bytes,
            "physical_memory_bytes": physical_memory_bytes,
            "ram_decision_rule": policy["ram_aware_execution"]["memory_decision_rule"],
            "memory_safety_reserve_bytes": safety_reserve,
            "memory_admission_rule": policy["ram_aware_execution"]["memory_admission_rule"],
            "benchmark_state": benchmark_state,
            "measured_discovery_loci": measured_loci,
            "maximum_observed_peak_ram_gib": measured_peak_gib,
            "reference_acquisition_storage_gate_applies": not reference_ready,
            "status": (
                "BLOCKED_BY_COMPUTE_STORAGE"
                if not reference_ready and free_bytes < int(policy["reference_minimum_free_bytes"])
                else benchmark_state
            ),
            "note": (
                "Operational disk/RAM snapshot excluded from immutable artifact hashes; "
                "the 35 GiB storage gate applies only while the reference must be acquired, "
                "and no compute-memory blocker is asserted before per-locus measurement."
            ),
        },
        "result_substitution_policy": "FORBIDDEN_SYNTHETIC_OR_PARTIAL_OUTPUTS_CANNOT_BE_PROMOTED_TO_REAL_LOCAL_RESULTS",
    }
    outputs[LOCK] = json.dumps(lock, indent=2, sort_keys=True) + "\n"
    return outputs


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build()
    if args.verify:
        # Operational free-space and handoff snapshots may change. Refresh them, then
        # require every immutable table to match the current scientific inputs.
        for path, value in expected.items():
            if path in {HANDOFF, LOCK}:
                continue
            if not path.is_file() or path.read_text(encoding="utf-8") != value:
                raise SystemExit(f"ERROR: Track B local input artifact missing or drifted: {path}")
        observed = json.loads(LOCK.read_text(encoding="utf-8"))
        for key in (
            "analysis_id", "selection_timing", "local_results_accessed_before_input_freeze",
            "policy_sha256", "pair_manifest_sha256", "locus_file_sha256",
            "input_artifact_sha256", "overlap_diagnostics", "result_substitution_policy",
        ):
            expected_lock = json.loads(expected[LOCK])
            if observed.get(key) != expected_lock.get(key):
                raise SystemExit(f"ERROR: Track B local input lock drifted: {key}")
        print("verified Track B local inputs: 3 pairs, 8 traits, 2495 loci")
        return
    refuse_post_result_refreeze()
    for path, value in expected.items():
        atomic_text(path, value)
    reference_status = json.loads(expected[LOCK])["reference_status"]
    print(
        "wrote Track B local inputs: 3 pairs, 8 traits, 2495 loci; "
        f"reference status {reference_status.lower()}"
    )


if __name__ == "__main__":
    main()
