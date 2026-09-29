#!/usr/bin/env python3
"""Combine locked FI, sleep-panel, and latent-factor LDSC h2 tables.

This is an aggregation only: it does not refit LDSC or modify the source tables.
All three inputs and their 20 expected trait IDs are fixed here so missing,
duplicated, or unexpected rows fail closed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path


SOURCES = {
    "frailty_index": ("frailty_paper/results/frailty_v1/h2_frailty_verification.tsv", {"frailty"}),
    "sleep_panel": (
        "frailty_paper/results/frailty_v1/h2_sleep_panel.tsv",
        {"chronotype", "napping", "snoring", "shortsleep", "sleepdur", "insomnia", "sleepiness", "sleep_apnea", "accel_sleep_duration", "sleep_efficiency", "longsleep", "sleep_timing"},
    ),
    "latent_factors": (
        "frailty_paper/results/frailty_v1/h2_latent_factors.tsv",
        {"frailty_general", "frailty_factor_1", "frailty_factor_2", "frailty_factor_3", "frailty_factor_4", "frailty_factor_5", "frailty_factor_6"},
    ),
}
EXPECTED = set().union(*(ids for _, ids in SOURCES.values()))
H2_COLUMNS = "scale h2 se z intercept intercept_se lambda_gc mean_chi2 ratio ldsc_n_eff_h2 ldsc_n_eff_h2_gt_12000 mixer_univariate_status pass_z pass_intercept verdict qc_reason".split()
FI_VERIFICATION_SOURCE = {
    "source_table": "/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/h2_verification.tsv",
    "source_table_sha256": "1b03c306816a49ad1c2737bf08264acd55c7eae9398214b47bc8dba7831062eb",
    "source_log": "/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/verification_logs/h2_frailty.log",
    "source_log_sha256": "43b888392487ee88081bac85281b0a02f304cc972761453d93994bde1fbbb6d2",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path, key: str = "trait") -> list[dict[str, str]]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or key not in reader.fieldnames:
            raise ValueError(f"Missing {key} header: {path}")
        rows = list(reader)
    return rows


def index_unique(rows: list[dict[str, str]], label: str, key_column: str = "trait") -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        key = row[key_column].strip()
        if not key or key in result:
            raise ValueError(f"Missing or duplicate {label} key: {key!r}")
        result[key] = row
    return result


def load_registry(path: Path, wanted: set[str]) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    # The acquired-resource file has one row per source resource. Keep the
    # resource identifier so traits are joined explicitly, never by display name.
    # Some collection-log rows (e.g. search batches) are descriptive and do
    # not have a stable resource_id; they are outside this h2 aggregation.
    selected = [r for r in rows if r.get("resource_id", "").strip() in wanted]
    return index_unique(selected, "resource_id", "resource_id")


def assemble(repo: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    registry_path = repo / "frailty_paper/manifests/all_acquired_resources.tsv"
    factor_path = repo / "frailty_paper/config/secondary_frailty_factors_v1.tsv"
    plan_path = repo / "frailty_paper/config/analysis_plan_v1.yaml"
    plan_lock_path = repo / "frailty_paper/config/analysis_plan_v1.lock.json"
    factors = index_unique(read_tsv(factor_path, "trait_id"), "factor trait_id", "trait_id")
    wanted_resources = {"atkins_2021_frailty_index"} | {f"sleep_panel_{t}" for t in SOURCES["sleep_panel"][1]} | {f"gwas_catalog_{r.get('source_id', '')}" for r in factors.values()}
    registry = load_registry(registry_path, wanted_resources)
    source_hashes: dict[str, str] = {
        str(registry_path.relative_to(repo)): sha256(registry_path),
        str(factor_path.relative_to(repo)): sha256(factor_path),
        str(plan_path.relative_to(repo)): sha256(plan_path),
        str(plan_lock_path.relative_to(repo)): sha256(plan_lock_path),
    }

    records: dict[str, dict[str, str]] = {}
    provenance = {}
    for family, (relative_path, expected_ids) in SOURCES.items():
        source_path = repo / relative_path
        if not source_path.is_file():
            raise FileNotFoundError(f"Required h2 source table missing: {relative_path}")
        source_hashes[relative_path] = sha256(source_path)
        if family == "frailty_index" and source_hashes[relative_path] != FI_VERIFICATION_SOURCE["source_table_sha256"]:
            raise ValueError("FI verification table checksum differs from the captured external verification")
        rows = index_unique(read_tsv(source_path), family)
        actual = set(rows)
        if actual != expected_ids:
            raise ValueError(f"{family} trait IDs mismatch; missing={sorted(expected_ids-actual)} unexpected={sorted(actual-expected_ids)}")
        for trait_id, h2row in rows.items():
            if h2row.get("verdict") != "PASS" or h2row.get("pass_z") != "True" or h2row.get("pass_intercept") != "True":
                raise ValueError(f"{family}/{trait_id} did not pass the recorded h2 gates")
            required_h2 = set(H2_COLUMNS) - {"ratio"}
            if any(not h2row.get(col, "").strip() for col in required_h2):
                raise ValueError(f"{family}/{trait_id} has missing required h2 fields")

            if family == "frailty_index":
                rid = "atkins_2021_frailty_index"
                meta = registry.get(rid)
            elif family == "sleep_panel":
                rid = f"sleep_panel_{trait_id}"
                meta = registry.get(rid)
            else:
                factor = factors.get(trait_id)
                if not factor:
                    raise ValueError(f"Missing factor metadata: {trait_id}")
                rid = f"gwas_catalog_{factor.get('source_id', '')}"
                meta = registry.get(rid)
            if not meta:
                raise ValueError(f"Missing resource metadata for {trait_id} ({rid})")
            raw_file = meta.get("file", "")
            raw_hash = meta.get("sha256", "")
            if not raw_file or len(raw_hash) != 64:
                raise ValueError(f"Missing registered raw file/hash for {trait_id}")

            log_name = Path(h2row.get("input_log", "")).name
            if family == "frailty_index":
                log_path = repo / "frailty_paper/results/frailty_v1/logs_fi_verification/h2_frailty.log"
                if sha256(log_path) != FI_VERIFICATION_SOURCE["source_log_sha256"]:
                    raise ValueError("FI verification log checksum differs from the captured external verification")
            else:
                log_path = None
            log_candidates = [repo / "frailty_paper/results/frailty_v1" / log_name,
                              repo / "frailty_paper/results/frailty_v1/logs" / log_name,
                              repo / "frailty_paper/results/frailty_v1/logs_sleep_panel" / log_name,
                              repo / "frailty_paper/results/frailty_v1/logs_latent_factors" / log_name]
            if log_path is None:
                log_path = next((p for p in log_candidates if p.is_file()), None)
            if log_path is None:
                raise FileNotFoundError(f"No local LDSC log found for {family}/{trait_id}: {h2row['input_log']}")
            log_relative = str(log_path.relative_to(repo))
            source_hashes[log_relative] = sha256(log_path)
            if raw_file not in source_hashes:
                source_hashes[raw_file] = raw_hash  # SHA-256 verified by the 347-resource manifest audit.

            row = {
                "trait_id": trait_id,
                "trait": meta.get("trait", ""),
                "analysis_family": family,
                "resource_id": rid,
                "accession": meta.get("accession", ""),
                "study": meta.get("study", ""),
                "genome_build": meta.get("genome_build", ""),
                "ancestry": meta.get("ancestry", ""),
                "sample_size": meta.get("sample_size", ""),
                "cases": meta.get("cases", ""),
                "controls": meta.get("controls", ""),
                "cohorts": meta.get("cohorts", ""),
                "ukb_overlap": meta.get("ukb_overlap", ""),
                "finngen_overlap": meta.get("finngen_overlap", ""),
                "raw_file": raw_file,
                "raw_sha256": raw_hash,
                **{c: h2row[c] for c in H2_COLUMNS},
                "input_table": relative_path,
                "input_log": log_relative,
            }
            records[trait_id] = row
            provenance[trait_id] = {"resource_id": rid, "raw_file": raw_file, "raw_sha256": raw_hash, "input_table": relative_path, "input_log": log_relative}

    if set(records) != EXPECTED or len(records) != 20:
        raise ValueError(f"Master must contain exactly 20 distinct IDs; got {len(records)}")
    return [records[k] for k in sorted(records)], source_hashes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = (args.output or repo / "frailty_paper/results/frailty_v1/ldsc/h2_master.tsv").resolve()
    rows, input_hashes = assemble(repo)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "artifact": str(output),
        "rows": len(rows),
        "trait_ids": [r["trait_id"] for r in rows],
        "gate_source": "frailty_paper/config/analysis_plan_v1.yaml: univariate_heritability (h2 Z >= 4; intercept <= 1.20)",
        "operation": "aggregation only; no LDSC estimates recalculated",
        "frailty_index_verification_origin": FI_VERIFICATION_SOURCE,
        "command": " ".join(sys.argv),
        "python_version": platform.python_version(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "input_sha256": dict(sorted(input_hashes.items())),
        "output_sha256": sha256(output),
    }
    manifest_path = output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {len(rows)} rows: {output}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
