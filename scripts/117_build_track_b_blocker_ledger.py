#!/usr/bin/env python3
"""Build a fail-closed ledger for every unexecuted Track B analysis stage."""

from __future__ import annotations

import argparse
import csv
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


ROOT = Path(__file__).resolve().parents[1]
TRACK_B = ROOT / "results/track_b"
OUT = TRACK_B / "04_compute_data_blockers.tsv"
REPORT = TRACK_B / "04_BLOCKER_LEDGER.md"
LOCK = TRACK_B / "04_compute_data_blockers.lock.json"

PAIR_MANIFEST = TRACK_B / "pair_manifest.tsv"
SOURCE_AUDIT = TRACK_B / "replication_source_audit.tsv"
SOURCE_LOCK = TRACK_B / "replication_source.lock.json"
DENSE_QC = TRACK_B / "03_dense_input_qc.tsv"
LOCAL_POLICY = ROOT / "config/track_b_local_analysis_policy.json"
LOCAL_LOCK = TRACK_B / "local_analysis_input.lock.json"
PLEIOTROPY_POLICY = ROOT / "config/pleiotropy_analysis_policy.json"
MOLECULAR_POLICY = ROOT / "config/molecular_analysis_policy.json"
INTERPRETATION_READINESS = ROOT / "results/tables/interpretation_source_readiness.tsv"

PAIR_B_RECEIPT = TRACK_B / "replication/finngen_r13_F5_ADHD_ingest_receipt.json"
PAIR_B_QC = TRACK_B / "replication/finngen_r13_F5_ADHD_ingest_qc.tsv"
PAIR_B_LDSC = TRACK_B / "replication/pair_b_ldsc.tsv"
PAIR_B_LDSC_RECEIPT = TRACK_B / "replication/pair_b_ldsc.provenance.json"
REPLICATION_SUMMARY = TRACK_B / "02_independent_global_replication.tsv"
REPLICATION_SUMMARY_RECEIPT = TRACK_B / "02_independent_global_replication.provenance.json"
LAVA_REFERENCE = ROOT / "ref/lava/ukb_v1.1/reference.provenance.json"
RAM_BENCHMARK = TRACK_B / "RAM_BENCHMARK.tsv"
RAM_BENCHMARK_PROVENANCE = TRACK_B / "RAM_BENCHMARK.provenance.json"
CHROMOSOME_INPUT_LOCK = TRACK_B / "lava_chromosome_inputs.provenance.json"
LAVA_CONTRACT = ROOT / "scripts/119_track_b_lava_contract.py"
LAVA_RESULT_LOCK = TRACK_B / "local/lava_results.provenance.json"
LAVA_EXECUTION_EVIDENCE = TRACK_B / "local/lava_execution_evidence.provenance.json"
LAVA_VALIDATOR = ROOT / "scripts/121_validate_track_b_lava.py"
LAVA_SUPERVISOR = ROOT / "scripts/131_run_track_b_lava_sequential.py"
HDL_REFERENCE = ROOT / "discovery_extension/ref/hdl_ukb_336k/reference.lock.json"
PLEIOFDR_REFERENCE = ROOT / "ref/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat"
PLEIOFDR_TEMPLATE = ROOT / "ref/pleiofdr/9545380.ref"

FIELDS = [
    "stage_order", "stage_id", "analysis", "scope", "status", "blocker_class",
    "observed_state", "required_to_unblock", "authoritative_evidence",
    "claim_status",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def memory_bytes() -> int:
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        return 0


def required_file(path: Path) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: required Track B evidence is missing: {rel(path)}")


def lava_completed_and_validated() -> bool:
    if not LAVA_RESULT_LOCK.is_file() or not LAVA_EXECUTION_EVIDENCE.is_file():
        return False
    semantic = subprocess.run(
        [sys.executable, str(LAVA_VALIDATOR)], cwd=ROOT,
        text=True, capture_output=True, check=False,
    )
    execution = subprocess.run(
        [sys.executable, str(LAVA_SUPERVISOR), "--verify-execution-evidence"], cwd=ROOT,
        text=True, capture_output=True, check=False,
    )
    if semantic.returncode != 0 or execution.returncode != 0:
        raise SystemExit(
            "ERROR: sealed LAVA completion evidence failed validation: "
            + semantic.stdout + semantic.stderr + execution.stdout + execution.stderr
        )
    return True


def classify_lava_readiness(
    reference_ok: bool, acquisition_storage_ok: bool, ram_state: str,
) -> tuple[str, list[str]]:
    classes = []
    if not reference_ok:
        classes.append("DATA_REFERENCE")
    if not reference_ok and not acquisition_storage_ok:
        classes.append("COMPUTE_STORAGE")
    if reference_ok and ram_state == "RAM_BENCHMARK_REQUIRED":
        classes.append("RAM_MEASUREMENT_REQUIRED")
    elif reference_ok and ram_state == "CONDITIONAL_RAM_BENCHMARK_PENDING":
        classes.append("CONDITIONAL_RAM_MEASUREMENT_PENDING")
    elif reference_ok and ram_state == "BLOCKED_BY_MEASURED_PER_LOCUS_RAM":
        classes.append("MEASURED_PER_LOCUS_MEMORY")
    status = "READY_TO_RUN" if not classes else (
        "READY_FOR_RAM_BENCHMARK" if classes == ["RAM_MEASUREMENT_REQUIRED"] else
        "READY_FOR_SEQUENTIAL_DISCOVERY" if classes == ["CONDITIONAL_RAM_MEASUREMENT_PENDING"] else
        "BLOCKED_BY_DATA_AND_COMPUTE" if "DATA_REFERENCE" in classes and len(classes) > 1 else
        "BLOCKED_BY_DATA" if classes == ["DATA_REFERENCE"] else "BLOCKED_BY_COMPUTE"
    )
    return status, classes


@lru_cache(maxsize=1)
def current_lava_fingerprint() -> str | None:
    if not CHROMOSOME_INPUT_LOCK.is_file() or not LAVA_REFERENCE.is_file():
        return None
    result = subprocess.run(
        [sys.executable, str(LAVA_CONTRACT), "--execution-fingerprint"],
        cwd=ROOT, text=True, capture_output=True, check=False,
    )
    value = result.stdout.strip()
    return value if result.returncode == 0 and len(value) == 64 else None


def lava_ram_benchmark_state(
    physical_memory_bytes: int, safety_reserve_bytes: int,
) -> tuple[str, float | None, int]:
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
    return (
        "MEASURED_SEQUENTIAL_RAM_PASS"
        if maximum * 1024**3 + safety_reserve_bytes <= physical_memory_bytes
        else "BLOCKED_BY_MEASURED_PER_LOCUS_RAM",
        maximum,
        len(attempt_rows),
    )


def row(order: int, stage_id: str, analysis: str, scope: str, status: str,
        blocker_class: str, observed: str, required: str, evidence: str,
        claim_status: str = "NO_SCIENTIFIC_RESULT") -> dict[str, object]:
    return {
        "stage_order": order,
        "stage_id": stage_id,
        "analysis": analysis,
        "scope": scope,
        "status": status,
        "blocker_class": blocker_class,
        "observed_state": observed,
        "required_to_unblock": required,
        "authoritative_evidence": evidence,
        "claim_status": claim_status,
    }


def build_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    for path in (
        PAIR_MANIFEST, SOURCE_AUDIT, SOURCE_LOCK, DENSE_QC, LOCAL_POLICY,
        LOCAL_LOCK, PLEIOTROPY_POLICY, MOLECULAR_POLICY, INTERPRETATION_READINESS,
        PAIR_B_QC, PAIR_B_RECEIPT, PAIR_B_LDSC, PAIR_B_LDSC_RECEIPT,
        REPLICATION_SUMMARY, REPLICATION_SUMMARY_RECEIPT,
    ):
        required_file(path)

    pairs = read_tsv(PAIR_MANIFEST)
    if [(r["pair_id"], r["sleep_trait"], r["external_trait"]) for r in pairs] != [
        ("A", "snoring", "parental_lifespan"),
        ("B", "insomnia", "adhd"),
        ("CONTROL", "insomnia", "frailty"),
    ]:
        raise SystemExit("ERROR: frozen Track B pair family drifted")

    sources = read_tsv(SOURCE_AUDIT)
    no_pair_a = [r for r in sources if r["pair_id"] == "A" and r["selection_decision"] == "NO_VALID_REPLICATION"]
    selected_b = [r for r in sources if r["pair_id"] == "B" and r["selection_decision"] == "SELECT_PRIMARY_REPLICATION"]
    if len(no_pair_a) != 1 or no_pair_a[0]["independence_status"] != "NO_VALID_REPLICATION":
        raise SystemExit("ERROR: Pair A no-valid-replication determination drifted")
    if len(selected_b) != 1 or selected_b[0]["source_id"] != "finngen_r13_F5_ADHD":
        raise SystemExit("ERROR: Pair B selected replication source drifted")

    replication = read_tsv(REPLICATION_SUMMARY)
    if [r["pair_id"] for r in replication] != ["A", "B"]:
        raise SystemExit("ERROR: canonical Track B replication pair order drifted")
    replication_by_pair = {r["pair_id"]: r for r in replication}
    pair_a = replication_by_pair["A"]
    pair_b = replication_by_pair["B"]
    if pair_a["replication_status"] != "NO_VALID_REPLICATION":
        raise SystemExit("ERROR: canonical Pair A terminal status drifted")
    if (
        pair_b["replication_status"] != "DIRECTIONAL_REPLICATION"
        or pair_b["direction"] != "CONCORDANT"
        or pair_b["h2_gate_status"] != "PASS"
    ):
        raise SystemExit("ERROR: canonical Pair B directional-replication status drifted")

    pair_b_ldsc = read_tsv(PAIR_B_LDSC)
    if len(pair_b_ldsc) != 1:
        raise SystemExit("ERROR: Pair B LDSC artifact must contain exactly one row")
    ldsc = pair_b_ldsc[0]
    ldsc_summary_fields = {
        "h2": "h2",
        "h2_se": "h2_SE",
        "rg": "replication_rg",
        "rg_se": "replication_SE",
        "rg_p": "replication_P",
        "cross_trait_intercept": "cross_trait_intercept",
        "cross_trait_intercept_se": "cross_trait_intercept_SE",
        "replication_class": "replication_status",
    }
    if ldsc["pair_id"] != "B" or any(
        ldsc[ldsc_field] != pair_b[summary_field]
        for ldsc_field, summary_field in ldsc_summary_fields.items()
    ):
        raise SystemExit("ERROR: Pair B LDSC and canonical replication summary disagree")

    ingest_receipt = json.loads(PAIR_B_RECEIPT.read_text(encoding="utf-8"))
    ldsc_receipt = json.loads(PAIR_B_LDSC_RECEIPT.read_text(encoding="utf-8"))
    replication_receipt = json.loads(REPLICATION_SUMMARY_RECEIPT.read_text(encoding="utf-8"))
    if ingest_receipt.get("qc_sha256") != sha256(PAIR_B_QC):
        raise SystemExit("ERROR: Pair B ingest QC hash does not match its receipt")
    if ldsc_receipt.get("output_sha256") != sha256(PAIR_B_LDSC):
        raise SystemExit("ERROR: Pair B LDSC hash does not match its receipt")
    if (
        replication_receipt.get("output_sha256") != sha256(REPLICATION_SUMMARY)
        or replication_receipt.get("pair_a_terminal_status") != "NO_VALID_REPLICATION"
        or replication_receipt.get("pair_b_terminal_status") != "DIRECTIONAL_REPLICATION"
    ):
        raise SystemExit("ERROR: canonical replication summary does not match its receipt")

    dense = read_tsv(DENSE_QC)
    if len(dense) != 5 or any(not r["overall_QC"].startswith("PASS_DENSE") for r in dense):
        raise SystemExit("ERROR: Track B dense-input gate is not intact")

    local_policy = json.loads(LOCAL_POLICY.read_text(encoding="utf-8"))
    local_lock = json.loads(LOCAL_LOCK.read_text(encoding="utf-8"))
    pleiotropy = json.loads(PLEIOTROPY_POLICY.read_text(encoding="utf-8"))
    molecular = json.loads(MOLECULAR_POLICY.read_text(encoding="utf-8"))
    interpretation = {r["source_id"]: r for r in read_tsv(INTERPRETATION_READINESS)}

    free = shutil.disk_usage(ROOT).free
    memory = memory_bytes()
    matlab = shutil.which("matlab")
    lava_storage_ok = free >= int(local_policy["reference_minimum_free_bytes"])
    lava_reference_ok = LAVA_REFERENCE.is_file()
    lava_safety_reserve = int(local_policy["ram_aware_execution"]["memory_safety_reserve_bytes"])
    lava_ram_state, lava_peak_gib, lava_measured_loci = lava_ram_benchmark_state(
        memory, lava_safety_reserve,
    )
    pleio_storage_ok = free >= int(pleiotropy["minimum_free_storage_bytes"])
    pleio_memory_ok = memory >= int(pleiotropy["pleiofdr_minimum_memory_bytes"])
    pleio_reference_ok = (
        PLEIOFDR_REFERENCE.is_file()
        and PLEIOFDR_REFERENCE.stat().st_size == int(pleiotropy["pleiofdr_reference_bytes"])
        and PLEIOFDR_TEMPLATE.is_file()
        and PLEIOFDR_TEMPLATE.stat().st_size == int(pleiotropy["pleiofdr_variant_template_bytes"])
    )
    molecular_storage_ok = free >= int(molecular["runtime"]["minimum_free_disk_gib_for_production"]) * 1024**3
    molecular_memory_ok = memory >= int(molecular["runtime"]["minimum_memory_gib"]) * 1024**3

    pair_b_observed = (
        f"{pair_b['replication_status']};rg={pair_b['replication_rg']};"
        f"SE={pair_b['replication_SE']};P={pair_b['replication_P']};CHECKSUM_LOCKED"
    )

    lava_status, lava_classes = classify_lava_readiness(
        lava_reference_ok, lava_storage_ok, lava_ram_state,
    )
    lava_complete = lava_completed_and_validated()
    if lava_complete:
        lava_status, lava_classes = "COMPLETE", []

    pleio_classes = []
    if not pleio_reference_ok or not matlab:
        pleio_classes.append("DATA_RUNTIME_REFERENCE")
    if not pleio_memory_ok or not pleio_storage_ok:
        pleio_classes.append("COMPUTE_MEMORY_STORAGE")
    pleio_status = "READY_TO_RUN" if not pleio_classes else (
        "BLOCKED_BY_DATA_AND_COMPUTE" if len(pleio_classes) == 2 else
        "BLOCKED_BY_DATA" if pleio_classes[0] == "DATA_RUNTIME_REFERENCE" else "BLOCKED_BY_COMPUTE"
    )

    rows = [
        row(1, "REPLICATION_A", "Independent exact replication", "A: snoring x parental lifespan",
            "NO_VALID_REPLICATION", "DATA_SOURCE",
            "No well-powered exact independent public analysis-ready parental-lifespan GWAS was identified",
            "Report NO_VALID_REPLICATION; do not substitute discovery, overlapping, or nonexact aging data",
            f"{rel(SOURCE_AUDIT)};{rel(SOURCE_LOCK)};{rel(REPLICATION_SUMMARY)}",
            "NO_VALID_REPLICATION"),
        row(2, "REPLICATION_B", "Independent external-trait LDSC", "B: insomnia x ADHD",
            "COMPLETE", "NONE", pair_b_observed,
            "No action required; preserve the frozen FinnGen result and its Finnish-founder/phenotype caveats",
            f"{rel(SOURCE_AUDIT)};{rel(SOURCE_LOCK)};{rel(PAIR_B_QC)};{rel(PAIR_B_RECEIPT)};{rel(PAIR_B_LDSC)};{rel(PAIR_B_LDSC_RECEIPT)};{rel(REPLICATION_SUMMARY)};{rel(REPLICATION_SUMMARY_RECEIPT)};scripts/114_stream_track_b_pair_b_replication.py;scripts/118_run_track_b_pair_b_ldsc.py;scripts/122_build_track_b_replication_summary.py",
            "DIRECTIONAL_REPLICATION"),
        row(3, "LOCAL_LAVA", "LAVA univariate and bivariate local sharing", "A;B;CONTROL across 2,495 loci",
            lava_status, ";".join(lava_classes) or "NONE",
            (f"SEALED_VALIDATED;result_lock_sha256={sha256(LAVA_RESULT_LOCK)};execution_evidence_sha256={sha256(LAVA_EXECUTION_EVIDENCE)}" if lava_complete else f"reference_present={str(lava_reference_ok).upper()};free_bytes={free};acquisition_minimum_free_bytes={local_policy['reference_minimum_free_bytes']};memory_bytes={memory};memory_safety_reserve_bytes={lava_safety_reserve};ram_benchmark_state={lava_ram_state};measured_loci={lava_measured_loci};maximum_observed_peak_ram_gib={lava_peak_gib if lava_peak_gib is not None else 'NA'}"),
            "Seal the official LAVA UKB EUR v1.1 reference, benchmark whole small/median/large loci in fresh processes, and proceed sequentially unless measured peak RSS exceeds this machine",
            f"{rel(LOCAL_POLICY)};{rel(LOCAL_LOCK)};{rel(TRACK_B / 'lava_pair_manifest.tsv')};scripts/119_track_b_lava_contract.py;scripts/120_run_track_b_lava.R;scripts/121_validate_track_b_lava.py;scripts/131_run_track_b_lava_sequential.py;{rel(RAM_BENCHMARK)}" + (f";{rel(LAVA_RESULT_LOCK)};{rel(LAVA_EXECUTION_EVIDENCE)}" if lava_complete else ""),
            "VALIDATED_EXECUTION_RESULT" if lava_complete else "NO_SCIENTIFIC_RESULT"),
        row(4, "LOCAL_CONDITIONAL", "Conditional local sharing", "A: BMI-only and sleep-apnea-only; B: MDD-only",
            "COMPLETE" if lava_complete else "BLOCKED_UPSTREAM", "NONE" if lava_complete else "UPSTREAM_LOCAL_RESULT",
            "SEALED_VALIDATED_CONDITIONAL_FAMILY" if lava_complete else "Covariates and timing are frozen; no valid LAVA local result exists",
            "No action required; preserve the full conditional family including failures and ineligible models" if lava_complete else "Complete primary LAVA local analysis, then test only predeclared covariates in supported loci",
            f"{rel(TRACK_B / 'local_conditional_manifest.tsv')};{rel(LOCAL_LOCK)}" + (f";{rel(LAVA_RESULT_LOCK)};{rel(LAVA_EXECUTION_EVIDENCE)}" if lava_complete else ""),
            "VALIDATED_EXECUTION_RESULT" if lava_complete else "NO_SCIENTIFIC_RESULT"),
        row(5, "LOCAL_ROBUSTNESS", "HDL-L or separately contracted rho-HESS", "Supported local loci",
            "BLOCKED_BY_DATA", "DATA_REFERENCE",
            f"hdl_reference_present={str(HDL_REFERENCE.is_file()).upper()};no Track B rho-HESS contract/reference",
            "Pin an ancestry-matched independent local-method reference and contract before result access",
            f"{rel(TRACK_B / 'local_method_robustness_plan.tsv')};discovery_extension/config/local_reference_sources.tsv"),
        row(6, "PLEIOTROPY_PLACO", "PLACO+ genome-wide pleiotropy", "A;B;CONTROL",
            "BLOCKED_UPSTREAM", "UPSTREAM_REPLICATION_AND_LOCAL",
            "Replication is terminally classified and five dense pair traits pass; primary local sharing remains incomplete",
            "Complete or terminally classify local sharing, then freeze the distinct Track B three-pair PLACO family",
            f"{rel(DENSE_QC)};{rel(PLEIOTROPY_POLICY)}"),
        row(7, "PLEIOTROPY_CONJFDR", "conjFDR supporting evidence", "A;B;CONTROL",
            pleio_status, ";".join(pleio_classes) or "NONE",
            f"reference_present={str(pleio_reference_ok).upper()};matlab_present={str(bool(matlab)).upper()};memory_bytes={memory};minimum_memory_bytes={pleiotropy['pleiofdr_minimum_memory_bytes']};free_bytes={free};minimum_free_bytes={pleiotropy['minimum_free_storage_bytes']}",
            "Pinned MATLAB runtime, sealed reference/template, >=16 GiB RAM, and >=20 GiB free; then run after PLACO+",
            f"{rel(PLEIOTROPY_POLICY)};scripts/42_pleiotropy_preflight.py"),
        row(8, "FINE_MAPPING_COLOC", "SuSiE-RSS fine-mapping and trait-trait coloc", "Validated pleiotropic/local loci",
            "BLOCKED_UPSTREAM", "UPSTREAM_LOCUS_AND_LD",
            "SuSiE/coloc code is pinned; no Track B validated locus family or signed locus LD exists",
            "Freeze the primary locus family, extract signed ancestry-matched LD, then run SuSiE-RSS and coloc.susie",
            "discovery_extension/results/fine_mapping/fine_mapping_readiness.tsv;discovery_extension/config/fine_mapping_colocalization_contract.json"),
        row(9, "TISSUE_CELL", "Tissue and cell-type enrichment", "Bulk and single-cell references",
            "BLOCKED_UPSTREAM", "UPSTREAM_VARIANTS_AND_GENES",
            "LDSC-SEG, FUMA scRNA, and CATlas scATAC reference families are ready; Track B variant/gene families are absent",
            "Complete validated loci/fine-mapping and apply prespecified multiple-testing families",
            rel(INTERPRETATION_READINESS)),
        row(10, "MOLECULAR_QTL_TWAS", "eQTL/sQTL/pQTL colocalization and TWAS", "Result-dependent tissues and genes",
            "BLOCKED_BY_DATA_AND_COMPUTE", "UPSTREAM_FINE_MAPPING;DATA_QTL_PAYLOAD;COMPUTE_STORAGE",
            f"Exact QTL payloads are result-dependent and absent; upstream fine-mapping absent;free_bytes={free};minimum_free_bytes={int(molecular['runtime']['minimum_free_disk_gib_for_production']) * 1024**3};memory_gate_pass={str(molecular_memory_ok).upper()};storage_gate_pass={str(molecular_storage_ok).upper()}",
            "Validated locus/variant family, checksum-locked exact QTL payloads, and molecular production resource gates",
            f"{rel(MOLECULAR_POLICY)};scripts/61_molecular_preflight.py"),
        row(11, "REGULATORY_SCATAC", "Regulatory overlap and scATAC enrichment", "Fine-mapped variants and nominated cells",
            "BLOCKED_UPSTREAM", "UPSTREAM_VARIANTS_AND_CELLS",
            "GENCODE, SCREEN, HOCOMOCO, ABC, PCHi-C, and CATlas sources are ready; Track B fine-mapped variants/cells absent",
            "Complete variant fine-mapping and cell prioritization before overlap/enrichment interpretation",
            rel(INTERPRETATION_READINESS)),
        row(12, "SPATIAL", "Spatial transcriptomic support", "Prioritized tissue-cell-gene chain",
            "BLOCKED_BY_DATA", "DATA_AND_AUTHORIZATION",
            "No authorized, checksum-locked spatial transcriptomic dataset is present",
            "Obtain an authorized context-matched spatial dataset and freeze its release before querying results",
            "results/track_b/02_LITERATURE_AUDIT.md"),
        row(13, "PATHWAYS", "Competitive pathway enrichment", "High-confidence Track B genes",
            "BLOCKED_UPSTREAM", "UPSTREAM_GENE_FAMILY",
            "Reactome, GO, MSigDB, and MAGMA gene-set resources are ready; no high-confidence Track B gene family exists",
            "Complete molecular/regulatory triangulation and freeze the supported gene family",
            rel(INTERPRETATION_READINESS)),
        row(14, "MODEL_FUNCTION", "Human perturbation and model-system validation", "Top genes/pathways",
            "BLOCKED_UPSTREAM", "UPSTREAM_GENE_AND_CELL_HYPOTHESIS",
            "No result-derived gene-cell hypothesis exists; a broad untargeted search would not be valid functional evidence",
            "Prioritize genes/cells first, then lock exact perturbation/model search routes and accession evidence",
            "discovery_extension/config/mechanistic_sources.tsv"),
        row(15, "CAUSAL_MR", "Bidirectional MR and pleiotropy-robust sensitivity", "A;B;CONTROL where instruments pass",
            "BLOCKED_UPSTREAM", "UPSTREAM_VALIDATED_INSTRUMENTS",
            "MR software sources are ready; Track B validated instruments and locus-level mechanism family are absent",
            "Freeze independent instruments after upstream validation; run IVW plus prespecified robust sensitivities",
            rel(INTERPRETATION_READINESS)),
        row(16, "FINAL_SYNTHESIS", "Evidence grading, mechanistic models, figures, and reproducible report", "All three pairs",
            "BLOCKED_UPSTREAM", "UPSTREAM_ANALYSIS_RESULTS",
            "Pair B directional replication and Pair A no-valid-replication are complete; no Track B local, pleiotropic-locus, fine-map, or mechanism result exists; no synthetic/empty substitution is permitted",
            "Complete each prior gate or retain its terminal no-valid-data outcome; never fill missing stages with synthetic/empty science tables",
            "results/track_b/02_LITERATURE_AUDIT.md;results/track_b/LOCAL_ANALYSIS_HANDOFF.md"),
    ]

    allowed_complete_claims = {
        ("REPLICATION_B", "DIRECTIONAL_REPLICATION"),
        ("LOCAL_LAVA", "VALIDATED_EXECUTION_RESULT"),
        ("LOCAL_CONDITIONAL", "VALIDATED_EXECUTION_RESULT"),
    }
    if any(
        r["status"] == "COMPLETE"
        and (r["stage_id"], r["claim_status"]) not in allowed_complete_claims
        for r in rows
    ):
        raise SystemExit("ERROR: blocker ledger contains an unapproved completed scientific claim")
    snapshot = {
        "free_bytes": free,
        "memory_bytes": memory,
        "matlab_path": matlab or "NONE",
        "lava_minimum_free_bytes": int(local_policy["reference_minimum_free_bytes"]),
        "lava_ram_benchmark_state": lava_ram_state,
        "lava_measured_loci": lava_measured_loci,
        "lava_maximum_observed_peak_ram_gib": lava_peak_gib,
        "lava_memory_safety_reserve_bytes": lava_safety_reserve,
        "pleiofdr_minimum_memory_bytes": int(pleiotropy["pleiofdr_minimum_memory_bytes"]),
        "pleiofdr_minimum_free_bytes": int(pleiotropy["minimum_free_storage_bytes"]),
        "molecular_minimum_free_bytes": int(molecular["runtime"]["minimum_free_disk_gib_for_production"]) * 1024**3,
    }
    return rows, snapshot


def tsv_text(rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def report_text(rows: list[dict[str, object]], snapshot: dict[str, object]) -> str:
    counts: dict[str, int] = {}
    for item in rows:
        counts[str(item["status"])] = counts.get(str(item["status"]), 0) + 1
    status_lines = "\n".join(f"- `{key}`: {value}" for key, value in sorted(counts.items()))
    return f"""# Track B compute and data blocker ledger

This is an execution-readiness artifact that binds terminal global-replication classifications. It is not a local or mechanistic analysis result. It records terminal data limitations and the exact gates that prevent downstream stages from being run or interpreted.

Current machine snapshot:

- free storage: {snapshot['free_bytes']} bytes ({int(snapshot['free_bytes']) / 1024**3:.3f} GiB)
- physical memory visible to the process: {snapshot['memory_bytes']} bytes ({int(snapshot['memory_bytes']) / 1024**3:.3f} GiB)
- LAVA RAM assessment: `{snapshot['lava_ram_benchmark_state']}` from {snapshot['lava_measured_loci']} measured discovery loci; maximum observed peak {snapshot['lava_maximum_observed_peak_ram_gib'] if snapshot['lava_maximum_observed_peak_ram_gib'] is not None else 'NA'} GiB plus a {snapshot['lava_memory_safety_reserve_bytes']} byte non-worker reserve
- MATLAB: `{snapshot['matlab_path']}`

Status counts:

{status_lines}

The Pair A replication outcome is terminal `NO_VALID_REPLICATION`; it is not a negative replication result. Pair B completed external-cohort LDSC with a significant concordant global genetic correlation (`rg={rows[1]['observed_state'].split(';')[1].split('=')[1]}`, `SE={rows[1]['observed_state'].split(';')[2].split('=')[1]}`, `P={rows[1]['observed_state'].split(';')[3].split('=')[1]}`), earning `DIRECTIONAL_REPLICATION` with the frozen Finnish-founder, phenotype-equivalence, and unverified individual-overlap caveats. This clears only the global replication gate. No local correlation, pleiotropic locus, colocalization, tissue/cell mechanism, molecular-QTL, pathway, functional, causal, or final mechanistic result is present or promoted by this ledger.

After storage or runtime conditions change, rerun `python3 scripts/117_build_track_b_blocker_ledger.py` before executing the next stage. Empty or synthetic downstream tables are forbidden as substitutes for real analyses.
"""


def build() -> dict[Path, str]:
    rows, snapshot = build_rows()
    table = tsv_text(rows)
    report = report_text(rows, snapshot)
    inputs = [
        PAIR_MANIFEST, SOURCE_AUDIT, SOURCE_LOCK, DENSE_QC, LOCAL_POLICY,
        LOCAL_LOCK, PLEIOTROPY_POLICY, MOLECULAR_POLICY, INTERPRETATION_READINESS,
        PAIR_B_QC, PAIR_B_RECEIPT, PAIR_B_LDSC, PAIR_B_LDSC_RECEIPT,
        REPLICATION_SUMMARY, REPLICATION_SUMMARY_RECEIPT,
    ]
    lock = {
        "schema_version": 1,
        "analysis_id": "track-b-v1.0-blocker-ledger",
        "artifact_role": "EXECUTION_READINESS_WITH_TERMINAL_REPLICATION_STATUS",
        "input_sha256": {rel(path): sha256(path) for path in inputs},
        "ledger_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "report_sha256": hashlib.sha256(report.encode()).hexdigest(),
        "row_count": len(rows),
        "stage_order": [str(r["stage_id"]) for r in rows],
        "current_compute_snapshot": snapshot,
        "substitution_policy": "FORBID_EMPTY_SYNTHETIC_OR_PLACEHOLDER_SCIENCE_OUTPUTS",
    }
    return {
        OUT: table,
        REPORT: report,
        LOCK: json.dumps(lock, indent=2, sort_keys=True) + "\n",
    }


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build()
    if args.verify:
        for path in expected:
            if not path.is_file():
                raise SystemExit(f"ERROR: Track B blocker ledger is missing: {rel(path)}")
        observed_rows = read_tsv(OUT)
        expected_rows, _ = build_rows()
        if len(observed_rows) != len(expected_rows):
            raise SystemExit("ERROR: Track B blocker ledger row count drifted")
        dynamic_observation_stages = {
            "LOCAL_LAVA", "PLEIOTROPY_CONJFDR", "MOLECULAR_QTL_TWAS",
        }
        for observed, wanted in zip(observed_rows, expected_rows):
            for field in FIELDS:
                if field == "observed_state" and wanted["stage_id"] in dynamic_observation_stages:
                    continue
                if observed[field] != str(wanted[field]):
                    raise SystemExit(
                        f"ERROR: Track B blocker ledger drifted: {wanted['stage_id']}:{field}"
                    )
        observed_lock = json.loads(LOCK.read_text(encoding="utf-8"))
        expected_lock = json.loads(expected[LOCK])
        for key in (
            "schema_version", "analysis_id", "artifact_role", "input_sha256",
            "row_count", "stage_order", "substitution_policy",
        ):
            if observed_lock.get(key) != expected_lock.get(key):
                raise SystemExit(f"ERROR: Track B blocker lock drifted: {key}")
        if observed_lock.get("ledger_sha256") != sha256(OUT):
            raise SystemExit("ERROR: Track B blocker ledger hash drifted")
        if observed_lock.get("report_sha256") != sha256(REPORT):
            raise SystemExit("ERROR: Track B blocker report hash drifted")
        print("verified Track B blocker ledger: terminal replication status bound; no mechanistic result promoted")
        return
    for path, value in expected.items():
        atomic_text(path, value)
    print("wrote Track B blocker ledger: terminal replication status bound; no mechanistic result promoted")


if __name__ == "__main__":
    main()
