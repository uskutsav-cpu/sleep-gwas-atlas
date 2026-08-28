#!/usr/bin/env python3
"""Execute one locked S-PrediXcan trait-by-context analysis without overwrite."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/twas_run_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/twas_run_manifest.lock.json")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest_path, lock_path, policy_path = root / args.manifest, root / args.manifest_lock, root / args.policy
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if lock.get("run_manifest_sha256") != sha256(manifest_path) or lock.get("policy_sha256") != sha256(policy_path):
        fail("TWAS run manifest differs from lock")
    manifest = read_tsv(manifest_path)
    if lock.get("run_ids_in_locked_order") != [row["run_id"] for row in manifest]:
        fail("TWAS run family/order differs from lock")
    selected = [row for row in manifest if row["run_id"] == args.run_id]
    if len(selected) != 1:
        fail("run_id must identify exactly one locked TWAS run")
    row = selected[0]
    try:
        gwas_n, gwas_h2 = int(row["gwas_N"]), float(row["gwas_h2"])
    except ValueError as exc:
        fail("locked TWAS variance-control inputs are nonnumeric")
        raise AssertionError from exc
    if row.get("variance_control_status") != "REQUIRED" or gwas_n <= 0 or not math.isfinite(gwas_h2) or not 0 < gwas_h2 <= 1:
        fail("locked TWAS variance-control inputs are inadmissible")
    mapped, mapped_lock = root / row["mapped_gwas_path"], root / row["mapped_gwas_lock_path"]
    if not mapped.is_file() or not mapped_lock.is_file():
        fail("mapped TWAS GWAS and lock are absent")
    mapping_lock = json.loads(mapped_lock.read_text(encoding="utf-8"))
    if mapping_lock.get("output_sha256") != sha256(mapped) or mapping_lock.get("run_manifest_lock_sha256") != sha256(lock_path):
        fail("mapped TWAS GWAS differs from lock")
    model, covariance = root / row["model_db_path"], root / row["covariance_path"]
    if not model.is_file() or sha256(model) != row["model_db_sha256"] or not covariance.is_file() or sha256(covariance) != row["covariance_sha256"]:
        fail("TWAS model database or covariance differs from lock")
    twas = policy["twas"]
    python = root / twas["runtime_environment"] / "bin/python"
    entrypoint = root / twas["entrypoint_path"]
    if not python.is_file() or not entrypoint.is_file() or sha256(entrypoint) != twas["entrypoint_sha256"]:
        fail("pinned MetaXcan runtime or entrypoint is absent")
    final = root / row["output_path"]
    final_dir = final.parent
    if final_dir.exists():
        fail("TWAS run directory already exists; refusing overwrite")
    if not args.execute:
        print("TWAS task is locked. Re-run with --execute to perform the real S-PrediXcan analysis.")
        return 0
    work_root = root / "work/twas_runs"
    work_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=args.run_id + ".", dir=work_root))
    try:
        output = staging / "spredixcan.csv"
        command = [
            str(python), str(entrypoint), "--model_db_path", str(model),
            "--model_db_snp_key", row["model_snp_key"], "--covariance", str(covariance),
            "--stream_covariance", "--gwas_file", str(mapped), "--snp_column", "SNP",
            "--effect_allele_column", "A1", "--non_effect_allele_column", "A2",
            "--beta_column", "BETA", "--se_column", "SE", "--separator", "\t",
            "--gwas_N", str(gwas_n), "--gwas_h2", format(gwas_h2, ".15g"),
            "--output_file", str(output), "--remove_ens_version", "--additional_output",
            "--throw", "--verbosity", "10",
        ]
        result = subprocess.run(command, cwd=root, capture_output=True, text=True)
        if result.returncode or not output.is_file() or output.stat().st_size == 0:
            fail(f"S-PrediXcan failed or omitted output: {(result.stdout + result.stderr).strip()}")
        with output.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames or [], list(reader)
        required = {
            "gene", "zscore", "pvalue", "n_snps_used", "n_snps_in_model",
            "uncalibrated_zscore", "uncalibrated_pvalue", "phi",
        }
        if not required.issubset(fields):
            fail(f"S-PrediXcan output lacks required fields: {sorted(required - set(fields))}")
        seen = set()
        for record in rows:
            gene = record["gene"]
            if not gene or gene in seen:
                fail("S-PrediXcan output has missing or duplicate gene")
            seen.add(gene)
            try:
                zscore, pvalue = float(record["zscore"]), float(record["pvalue"])
                raw_z, raw_p, phi = (
                    float(record["uncalibrated_zscore"]), float(record["uncalibrated_pvalue"]),
                    float(record["phi"]),
                )
                used, model_n = int(float(record["n_snps_used"])), int(float(record["n_snps_in_model"]))
            except ValueError as exc:
                fail(f"invalid S-PrediXcan result for {gene}")
                raise AssertionError from exc
            if (
                not all(math.isfinite(value) for value in (zscore, pvalue, raw_z, raw_p, phi))
                or not 0 <= pvalue <= 1 or not 0 <= raw_p <= 1 or phi < 0
                or used < 1 or model_n < used
            ):
                fail(f"S-PrediXcan result fails numeric QC for {gene}")
        provenance = {
            "schema_version": "atlas-v1.0-twas-run.1", "run_id": args.run_id,
            "manifest_sha256": sha256(manifest_path), "manifest_lock_sha256": sha256(lock_path),
            "policy_sha256": sha256(policy_path), "mapped_gwas_sha256": sha256(mapped),
            "mapped_gwas_lock_sha256": sha256(mapped_lock), "model_db_sha256": sha256(model),
            "covariance_sha256": sha256(covariance), "entrypoint_sha256": sha256(entrypoint),
            "variance_control": {"status": "APPLIED", "gwas_N": gwas_n, "gwas_h2": gwas_h2,
                                 "gwas_h2_scale": row["gwas_h2_scale"], "phi_required": True},
            "command": command, "stdout": result.stdout.strip(), "stderr": result.stderr.strip(),
            "row_count": len(rows), "output_sha256": sha256(output), "claim_limit": policy["claim_limit"],
        }
        (staging / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        final_dir.parent.mkdir(parents=True, exist_ok=True); os.replace(staging, final_dir)
        print(f"TWAS_RUN_OK run={args.run_id} genes={len(rows)}")
        return 0
    finally:
        if staging.exists(): shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
