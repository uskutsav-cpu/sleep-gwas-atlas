#!/usr/bin/env python3
"""Freeze or verify the pre-analysis Track B repository checkpoint."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


CANONICAL_INPUTS = [
    Path("config/analysis_panel.tsv"),
    Path("config/analysis_panel.lock.json"),
    Path("results/tables/rg_matrix.tsv"),
    Path("results/tables/h2_summary.tsv"),
    Path("results/tables/trait_readiness.tsv"),
    Path("results/analysis/phase1_master_analysis.tsv"),
    Path("results/analysis/literature_novelty_audit.tsv"),
    Path("results/analysis/local_followup_priorities.tsv"),
    Path("discovery_extension/core_checkpoint.json"),
    Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
    Path("discovery_extension/results/replication/replication_results.tsv"),
    Path("results/tables/ldsc_covariance_45x45.tsv"),
    Path("results/tables/ldsc_sampling_covariance_1035x1035.tsv.gz"),
    Path("results/tables/dense_harmonization_readiness.tsv"),
    Path("results/tables/pleiotropy_input_readiness.tsv"),
    Path("results/tables/interpretation_source_readiness.tsv"),
    Path("environment/tool_versions.tsv"),
]

OPERATIONAL_SNAPSHOTS = [
    Path("results/tables/lava_input_diagnostics.json"),
    Path("results/tables/fine_mapping_preflight.json"),
    Path("results/tables/molecular_preflight.json"),
]

OUTPUT = Path("results/track_b/00_repository_checkpoint.json")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty canonical input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def git(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], check=True, text=True, capture_output=True,
    )
    return result.stdout.strip()


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def runtime_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "UNAVAILABLE"
    lines = (result.stdout or result.stderr).splitlines()
    return lines[0].strip() if lines else "UNAVAILABLE"


def validate_core_counts() -> dict[str, object]:
    panel = read_tsv(Path("config/analysis_panel.tsv"))
    rg = read_tsv(Path("results/tables/rg_matrix.tsv"))
    h2 = read_tsv(Path("results/tables/h2_summary.tsv"))
    extension = read_tsv(Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"))
    replication = read_tsv(Path("discovery_extension/results/replication/replication_results.tsv"))
    sleep = [row for row in panel if row["domain"] == "sleep"]
    external = [row for row in panel if row["domain"] != "sleep"]
    expected_pairs = {(a["trait_id"], b["trait_id"]) for a in sleep for b in external}
    observed_pairs = {(row["sleep_trait"], row["disease_trait"]) for row in rg}
    if len(panel) != 45 or len(sleep) != 12 or len(external) != 33:
        fail("locked analysis panel is not exactly 12 sleep + 33 external traits")
    if len(rg) != 396 or observed_pairs != expected_pairs:
        fail("canonical rg table is not the complete locked 12 x 33 family")
    if len(h2) != 45 or {row["trait"] for row in h2} != {row["trait_id"] for row in panel}:
        fail("canonical h2 table is not trait-complete")
    if len(extension) != 1200:
        fail("extension rg table is not the locked 1,200-pair family")
    return {
        "locked_panel_traits": len(panel),
        "sleep_traits": len(sleep),
        "external_traits": len(external),
        "canonical_rg_pairs": len(rg),
        "canonical_h2_traits": len(h2),
        "extension_rg_pairs": len(extension),
        "extension_replication_candidates": len(replication),
        "extension_replicated_pairs": sum(row.get("replication_class") == "REPLICATED" for row in replication),
    }


def dense_readiness() -> dict[str, object]:
    rows = read_tsv(Path("results/tables/pleiotropy_input_readiness.tsv"))
    ready = [row["trait_id"] for row in rows if row["input_status"] == "READY_FULL_SUMSTATS"]
    blocked = [row["trait_id"] for row in rows if row["input_status"] != "READY_FULL_SUMSTATS"]
    return {
        "source": "results/tables/pleiotropy_input_readiness.tsv",
        "available_count": len(ready),
        "available_traits": ready,
        "missing_or_hapmap3_only_count": len(blocked),
        "missing_or_hapmap3_only_traits": blocked,
        "schema_caveat": (
            "Repository-ready dense files carry SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N; "
            "INFO and case/control counts are retained in QC/manifests rather than row-wise. "
            "Track B dense-input QC is therefore still required."
        ),
    }


def resource_readiness() -> dict[str, object]:
    interpretation = read_tsv(Path("results/tables/interpretation_source_readiness.tsv"))
    cell_rows = [row for row in interpretation if row["analysis_family"] == "cell_type"]
    qtl_metadata = Path("ref/molecular/eqtl_catalogue_r7/dataset_metadata_r7.tsv")
    lava_reference = Path("ref/lava/ukb_v1.1/reference.provenance.json")
    lava_loci = Path("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
    ldsc_dir = Path("ref/eur_w_ld_chr")
    return {
        "ld_references": [
            {
                "resource": "1000G_EUR_LDSC_scores",
                "path": str(ldsc_dir),
                "status": "AVAILABLE" if len(list(ldsc_dir.glob("*.l2.ldscore.gz"))) >= 22 else "INCOMPLETE",
            },
            {
                "resource": "LAVA_2495_locus_definition",
                "path": str(lava_loci),
                "status": "AVAILABLE" if lava_loci.is_file() else "MISSING",
            },
            {
                "resource": "LAVA_UKB_EUR_LD_v1.1",
                "path": "ref/lava/ukb_v1.1",
                "status": "AVAILABLE" if lava_reference.is_file() else "MISSING_BLOCKED_BY_DISK",
            },
            {
                "resource": "pleioFDR_1KG_Phase3_EUR_reference",
                "path": "ref/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat",
                "status": "MISSING",
            },
        ],
        "qtl_resources": [
            {
                "resource": "eQTL_Catalogue_release_7_metadata",
                "path": str(qtl_metadata),
                "status": "AVAILABLE_METADATA_ONLY" if qtl_metadata.is_file() else "MISSING",
            },
            {
                "resource": "cell_specific_QTL_results",
                "path": "results/tables/molecular_evidence.tsv",
                "status": "MISSING_BLOCKED_UPSTREAM",
            },
        ],
        "single_cell_and_chromatin_resources": [resource_record(row) for row in cell_rows],
        "authorized_spatial_resource": "NOT_FOUND",
    }


def resource_record(row: dict[str, str]) -> dict[str, str]:
    path = Path(row["local_path"])
    status = row["readiness_status"]
    blocker = row["blocker"]
    if status == "READY" and not path.exists():
        status = "INCONSISTENT_REPORTED_READY_OUTPUT_MISSING"
        blocker = "readiness ledger reports READY but the referenced local artifact is absent"
    return {
        "source_id": row["source_id"],
        "method": row["method"],
        "release": row["exact_release"],
        "path": row["local_path"],
        "status": status,
        "blocker": blocker,
    }


def baseline_dirty_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for line in git("status", "--porcelain=v1").splitlines():
        if not line:
            continue
        status, raw_path = line[:2], line[3:]
        path = Path(raw_path)
        record: dict[str, object] = {"status": status, "path": raw_path}
        if path.is_file():
            record.update({"bytes": path.stat().st_size, "sha256": sha256(path)})
        records.append(record)
    return records


def build_checkpoint() -> dict[str, object]:
    missing = [
        str(path) for path in [*CANONICAL_INPUTS, *OPERATIONAL_SNAPSHOTS]
        if not path.is_file()
    ]
    if missing:
        fail(f"canonical inputs missing: {missing}")
    canonical = {str(path): file_record(path) for path in CANONICAL_INPUTS}
    disk = shutil.disk_usage(Path.cwd())
    preflight = json.loads(Path("results/tables/fine_mapping_preflight.json").read_text(encoding="utf-8"))
    tool_rows = read_tsv(Path("environment/tool_versions.tsv"))
    local_result_candidates = [
        Path("results/tables/lava_local_results.tsv"),
        Path("results/atlas/shared_loci.tsv"),
        Path("results/atlas/finemapping_variants.tsv"),
        Path("results/atlas/trait_colocalization.tsv"),
    ]
    present_local_results = [str(path) for path in local_result_candidates if path.is_file() and path.stat().st_size]
    if present_local_results:
        fail(f"unexpected pre-existing local/mechanistic results: {present_local_results}")
    return {
        "schema_version": "1.0.0",
        "analysis_id": "sleep-circadian-genetic-atlas-track-b",
        "checkpoint_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "checkpoint_timing": "before_Track_B_pair_selection_or_local_result_access",
        "git": {
            "commit": git("rev-parse", "HEAD"),
            "branch": git("branch", "--show-current"),
            "remote": git("remote", "get-url", "origin"),
            "baseline_dirty_entries": baseline_dirty_records(),
        },
        "integrity": validate_core_counts(),
        "analysis_panel_hash": canonical["config/analysis_panel.tsv"]["sha256"],
        "rg_table_hash": canonical["results/tables/rg_matrix.tsv"]["sha256"],
        "input_hashes": canonical,
        "operational_snapshots": {
            str(path): file_record(path) for path in OPERATIONAL_SNAPSHOTS
        },
        "operational_snapshot_policy": (
            "Observed preflight snapshots are recorded but not immutable inputs because tests and "
            "read-only preflight reruns legitimately refresh timestamps and free-space observations."
        ),
        "dense_full_resolution_gwas": dense_readiness(),
        "resources": resource_readiness(),
        "software_versions": {
            "pinned": {row["component"]: row["version_or_commit"] for row in tool_rows},
            "observed_python": sys.version.split()[0],
            "observed_R": runtime_version(["R", "--version"]),
            "observed_git": runtime_version(["git", "--version"]),
        },
        "compute_environment": {
            "system": platform.system(),
            "release": platform.release(),
            "architecture": platform.machine(),
            "cpu_count": os.cpu_count(),
            "memory_bytes": preflight["machine"]["memory_bytes"],
            "disk_total_bytes": disk.total,
            "disk_free_bytes": disk.free,
            "disk_free_gib": round(disk.free / 1024**3, 3),
            "assessment": "BLOCKED_BY_COMPUTE_FOR_LAVA_PLEIOFDR_AND_LARGE_FINE_MAPPING",
        },
        "real_analysis_outputs": {
            "phase1_global_rg": "COMPLETE_396_OF_396",
            "phase1_h2": "COMPLETE_45_OF_45",
            "extension_global_rg": "COMPLETE_1200_OF_1200",
            "extension_replication": "COMPLETE_LOCKED_FAMILY",
            "genomic_sem_covariance": "COMPLETE_WITH_RECORDED_NUMERICAL_WARNINGS",
            "track_b_local_or_mechanistic_results": "NONE_AT_CHECKPOINT",
        },
        "hard_constraints": {
            "locked_45_trait_core_unchanged": True,
            "canonical_phase1_rg_unchanged": True,
            "track_b_local_results_accessed_before_checkpoint": False,
            "synthetic_results_counted_as_scientific": False,
        },
    }


def verify_checkpoint() -> None:
    if not OUTPUT.is_file() or OUTPUT.stat().st_size == 0:
        fail(f"checkpoint is missing: {OUTPUT}")
    checkpoint = json.loads(OUTPUT.read_text(encoding="utf-8"))
    observed = {str(path): file_record(path) for path in CANONICAL_INPUTS}
    if checkpoint.get("input_hashes") != observed:
        fail("one or more canonical inputs differ from the Track B checkpoint")
    if checkpoint.get("analysis_panel_hash") != observed["config/analysis_panel.tsv"]["sha256"]:
        fail("analysis-panel hash mismatch")
    if checkpoint.get("rg_table_hash") != observed["results/tables/rg_matrix.tsv"]["sha256"]:
        fail("canonical rg-table hash mismatch")
    validate_core_counts()
    print(
        "TRACK_B_REPOSITORY_CHECKPOINT_OK "
        f"commit={checkpoint['git']['commit']} panel={checkpoint['analysis_panel_hash']} "
        f"rg={checkpoint['rg_table_hash']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_checkpoint()
        return
    if OUTPUT.exists():
        fail(f"checkpoint already exists; verify it with --verify: {OUTPUT}")
    payload = build_checkpoint()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "TRACK_B_REPOSITORY_CHECKPOINT_FROZEN "
        f"commit={payload['git']['commit']} dense={payload['dense_full_resolution_gwas']['available_count']}/45 "
        f"sha256={sha256(OUTPUT)}"
    )


if __name__ == "__main__":
    main()
