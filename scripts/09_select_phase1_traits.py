#!/usr/bin/env python3
"""Create the auditable Phase 1 inclusion table for the locked panel."""
import argparse
import hashlib
import os
import subprocess
import sys

import pandas as pd


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def boolean_series(series, column):
    normalized = series.astype(str).str.strip().str.lower()
    invalid = ~normalized.isin({"true", "false"})
    if invalid.any():
        raise SystemExit(
            f"ERROR: readiness column {column} contains non-boolean values: "
            f"{sorted(normalized[invalid].unique())}"
        )
    return normalized.eq("true")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--lock", default="config/analysis_panel.lock.json")
    parser.add_argument("--readiness", required=True)
    parser.add_argument("--h2", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    subprocess.run(
        [sys.executable, "scripts/00_validate_panel.py", "--manifest", args.config,
         "--lock", args.lock],
        check=True,
    )
    config = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
    readiness = pd.read_csv(args.readiness, sep="\t", dtype=str).fillna("")
    h2 = pd.read_csv(args.h2, sep="\t", dtype=str).fillna("")
    for frame, key, name in [
        (config, "trait_id", "analysis panel"),
        (readiness, "trait_id", "readiness ledger"),
        (h2, "trait", "h2 summary"),
    ]:
        if key not in frame.columns or frame[key].duplicated().any():
            raise SystemExit(f"ERROR: {name} must contain unique {key} values")
    required_readiness = {
        "manifest_sha256", "source_verified", "ldsc_ready",
        "liability_h2_ready", "phase1_pass", "readiness_stage",
    }
    missing = required_readiness.difference(readiness.columns)
    if missing:
        raise SystemExit(f"ERROR: readiness ledger missing columns: {sorted(missing)}")
    if not {"trait", "verdict", "qc_reason"}.issubset(h2.columns):
        raise SystemExit("ERROR: h2 summary must include trait, verdict, and qc_reason")

    panel_hash = file_sha256(args.config)
    ledger_hashes = set(readiness["manifest_sha256"])
    if ledger_hashes != {panel_hash}:
        raise SystemExit(
            "ERROR: readiness ledger was not generated from the current analysis_panel.tsv"
        )
    panel_ids = set(config["trait_id"])
    readiness_ids = set(readiness["trait_id"])
    if readiness_ids != panel_ids:
        raise SystemExit("ERROR: readiness ledger trait set differs from the locked panel")
    unknown_h2 = set(h2["trait"]).difference(panel_ids)
    if unknown_h2:
        raise SystemExit(f"ERROR: h2 summary contains traits outside the locked panel: {sorted(unknown_h2)}")

    readiness = readiness.copy()
    for column in ["source_verified", "ldsc_ready", "liability_h2_ready", "phase1_pass"]:
        readiness[column] = boolean_series(readiness[column], column)
    merged = config.merge(readiness, on="trait_id", how="left", suffixes=("", "_readiness"))
    merged = merged.merge(
        h2[["trait", "verdict", "qc_reason"]],
        how="left",
        left_on="trait_id",
        right_on="trait",
    )
    h2_pass = merged["verdict"].eq("PASS")
    derived_pass = (
        merged["source_verified"]
        & merged["ldsc_ready"]
        & merged["liability_h2_ready"]
        & h2_pass
    )
    if not derived_pass.equals(merged["phase1_pass"]):
        bad = merged.loc[derived_pass.ne(merged["phase1_pass"]), "trait_id"].tolist()
        raise SystemExit(f"ERROR: readiness phase1_pass disagrees with h2 QC for: {bad}")
    merged["include_phase1"] = derived_pass
    merged["inclusion_reason"] = "included"
    merged.loc[~merged["source_verified"], "inclusion_reason"] = "source_not_verified"
    merged.loc[merged["source_verified"] & ~merged["ldsc_ready"], "inclusion_reason"] = "ldsc_not_ready"
    merged.loc[
        merged["source_verified"] & merged["ldsc_ready"] & ~merged["liability_h2_ready"],
        "inclusion_reason",
    ] = "liability_h2_not_ready"
    merged.loc[
        merged["source_verified"] & merged["ldsc_ready"]
        & merged["liability_h2_ready"] & ~h2_pass,
        "inclusion_reason",
    ] = "h2_qc_not_pass"

    output = merged[[
        "panel_version", "manifest_sha256", "trait_id", "label", "domain", "type",
        "declared_source_status", "source_verified", "harmonization_ready",
        "ldsc_ready", "liability_h2_ready", "readiness_stage", "verdict",
        "qc_reason", "include_phase1", "inclusion_reason",
    ]]
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    output.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote Phase 1 inclusion table: {args.out}")
    print(f"Included: {int(output['include_phase1'].sum())} / {len(output)}")


if __name__ == "__main__":
    main()
