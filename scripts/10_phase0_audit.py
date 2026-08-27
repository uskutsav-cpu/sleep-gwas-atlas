#!/usr/bin/env python3
"""Derive evidence-backed readiness gates for every locked panel trait.

The manifest declares source selection only. This audit derives later stages
from source-registry evidence, metadata, local files, and (optionally) the h2
QC table. A missing artifact is a visible blocker, never an implicit status.
"""
import argparse
import hashlib
import os
import subprocess
import sys

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "panel_version", "trait_id", "label", "domain", "type",
    "phenotype_definition", "source_note", "dataset_version", "source_id",
    "raw_file", "ncase", "ncontrol", "n_total", "pop_prev", "build",
    "source_status", "pmid", "doi", "ancestry", "pop_prev_citation",
}
SOURCE_COLUMNS = {
    "source_id", "trait_ids", "source_page_url", "download_url", "access",
    "archive_sha256", "raw_files", "pmid", "ancestry_reported",
    "build_status", "acquisition_status",
}
SCHEMA_COLUMNS = {
    "source_id", "trait_id", "schema_status", "file_format", "variant_id", "chromosome",
    "position", "effect_allele", "other_allele", "effect",
    "effect_convention", "standard_error", "p_value", "eaf", "info",
    "sample_size", "notes",
}
MISSING_TEXT = {"", "na", "nan", "none", "null", "unresolved", "unregistered"}


def populated(value):
    if value is None:
        return False
    normalized = str(value).strip().lower()
    return normalized not in MISSING_TEXT and not normalized.startswith("unresolved")


def numeric(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def manifest_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_public_sources(path, trait_ids):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: source registry not found: {path}")
    sources = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = SOURCE_COLUMNS.difference(sources.columns)
    if missing:
        raise SystemExit(f"ERROR: source registry missing columns: {sorted(missing)}")
    by_trait = {}
    for _, source in sources.iterrows():
        mapped_ids = [item.strip() for item in source["trait_ids"].split(",") if item.strip()]
        for trait_id in mapped_ids:
            if trait_id not in trait_ids:
                continue  # source registry may retain archived candidate sources
            if trait_id in by_trait:
                raise SystemExit(f"ERROR: multiple source rows claim trait_id {trait_id}")
            by_trait[trait_id] = source
    return by_trait


def load_source_schemas(path, selected_sources):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: GWAS schema registry not found: {path}")
    schemas = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = SCHEMA_COLUMNS.difference(schemas.columns)
    if missing:
        raise SystemExit(f"ERROR: GWAS schema registry missing columns: {sorted(missing)}")
    duplicate_pairs = schemas.duplicated(subset=["source_id", "trait_id"], keep=False)
    if duplicate_pairs.any():
        duplicates = schemas.loc[duplicate_pairs, ["source_id", "trait_id"]].to_dict("records")
        raise SystemExit(f"ERROR: duplicate source/trait schema rows: {duplicates}")
    invalid_pairs = []
    for _, schema in schemas.iterrows():
        selected_source = selected_sources.get(schema["trait_id"])
        if selected_source != schema["source_id"]:
            invalid_pairs.append({
                "source_id": schema["source_id"],
                "trait_id": schema["trait_id"],
                "selected_source_id": selected_source or "UNSELECTED_TRAIT",
            })
    if invalid_pairs:
        raise SystemExit(f"ERROR: source schemas do not match selected trait sources: {invalid_pairs}")
    allowed = {"SCHEMA_PENDING", "SCHEMA_VERIFIED"}
    invalid = schemas.loc[~schemas["schema_status"].isin(allowed), ["source_id", "schema_status"]]
    if len(invalid):
        raise SystemExit(f"ERROR: invalid source schema status: {invalid.to_dict('records')}")
    return {(row["source_id"], row["trait_id"]): row for _, row in schemas.iterrows()}


def load_h2(path):
    if not path:
        return {}
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: h2 table not found: {path}")
    h2 = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if not {"trait", "verdict"}.issubset(h2.columns):
        raise SystemExit("ERROR: h2 table must contain trait and verdict columns")
    if h2["trait"].duplicated().any():
        raise SystemExit("ERROR: h2 table contains duplicate trait rows")
    return {row["trait"]: row for _, row in h2.iterrows()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--lock", default="config/analysis_panel.lock.json")
    parser.add_argument("--sources", default="config/public_gwas_sources.tsv")
    parser.add_argument("--schemas", default="config/gwas_schemas.tsv")
    parser.add_argument("--h2")
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--harmonized-dir", default="data/harmonized")
    parser.add_argument("--munged-dir", default="data/munged")
    parser.add_argument("--out", default="results/tables/trait_readiness.tsv")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    subprocess.run(
        [sys.executable, "scripts/00_validate_panel.py", "--manifest", args.config,
         "--lock", args.lock],
        check=True,
    )
    config = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
    missing = REQUIRED_COLUMNS.difference(config.columns)
    if missing:
        raise SystemExit(f"ERROR: analysis panel missing columns: {sorted(missing)}")
    sources_by_trait = load_public_sources(args.sources, set(config["trait_id"]))
    selected_sources = dict(zip(config["trait_id"], config["source_id"]))
    schemas_by_source_trait = load_source_schemas(args.schemas, selected_sources)
    h2_by_trait = load_h2(args.h2)
    panel_hash = manifest_sha256(args.config)

    rows = []
    for _, trait in config.iterrows():
        trait_id = trait["trait_id"]
        source = sources_by_trait.get(trait_id)
        schema = schemas_by_source_trait.get((trait["source_id"], trait_id))

        source_issues = []
        if trait["source_status"] != "SOURCE_VERIFIED":
            source_issues.append("source_status_pending")
        if source is None:
            source_issues.append("source_not_registered")
        else:
            if trait["source_id"] != source["source_id"]:
                source_issues.append("manifest_source_id_mismatch")
            if source["access"].upper() != "PUBLIC":
                source_issues.append("source_not_public")
            if source["acquisition_status"].upper() != "SOURCE_VERIFIED_PUBLIC":
                source_issues.append("source_acquisition_not_verified")
            if not populated(source["source_page_url"]) or not populated(source["download_url"]):
                source_issues.append("source_url_missing")
            if not populated(source["archive_sha256"]):
                source_issues.append("source_checksum_missing")
            raw_files = {item.strip() for item in source["raw_files"].split(",")}
            if trait["raw_file"] not in raw_files:
                source_issues.append("source_registry_raw_file_mismatch")
        source_verified = not source_issues

        harmonization_issues = []
        if not source_verified:
            harmonization_issues.append("source_not_verified")
        if not populated(trait["phenotype_definition"]):
            harmonization_issues.append("phenotype_definition_unresolved")
        if not populated(trait["dataset_version"]):
            harmonization_issues.append("dataset_version_unresolved")
        if not populated(trait["pmid"]):
            harmonization_issues.append("primary_publication_unresolved")
        if trait["ancestry"].upper() != "EUR":
            harmonization_issues.append("eur_subset_unresolved")
        if trait["build"] not in {"hg19", "GRCh37"}:
            harmonization_issues.append("not_hg19_requires_explicit_build_decision")
        if source is not None and source["build_status"].upper() != "HEADER_VALIDATED_HG19":
            harmonization_issues.append("source_build_not_header_validated_hg19")
        if schema is None:
            harmonization_issues.append("source_schema_not_registered")
            schema_status = "UNREGISTERED"
        else:
            schema_status = schema["schema_status"]
            if schema_status != "SCHEMA_VERIFIED":
                harmonization_issues.append("source_schema_not_verified")
        if trait["type"] == "binary":
            if not numeric(trait["ncase"]) or not numeric(trait["ncontrol"]):
                harmonization_issues.append("binary_ncase_or_ncontrol_unresolved")
        elif not numeric(trait["n_total"]):
            harmonization_issues.append("continuous_n_total_unresolved")
        raw_path = os.path.join(args.raw_dir, trait["raw_file"])
        if not os.path.isfile(raw_path):
            harmonization_issues.append("registered_raw_file_not_materialized")
        harmonization_ready = not harmonization_issues

        ldsc_issues = []
        if not harmonization_ready:
            ldsc_issues.append("harmonization_not_ready")
        qc_path = os.path.join(args.harmonized_dir, f"{trait_id}.qc.txt")
        munged_path = os.path.join(args.munged_dir, f"{trait_id}.sumstats.gz")
        if not os.path.isfile(qc_path):
            ldsc_issues.append("harmonization_qc_ledger_missing")
        if not os.path.isfile(munged_path):
            ldsc_issues.append("hapmap3_sumstats_missing")
        ldsc_ready = not ldsc_issues

        liability_issues = []
        if not ldsc_ready:
            liability_issues.append("ldsc_not_ready")
        if trait["type"] == "binary":
            prevalence = numeric(trait["pop_prev"])
            if prevalence is None or not 0 < prevalence < 1:
                liability_issues.append("population_prevalence_unresolved")
            if not populated(trait["pop_prev_citation"]):
                liability_issues.append("population_prevalence_citation_missing")
        liability_h2_ready = not liability_issues

        phase1_issues = []
        if not liability_h2_ready:
            phase1_issues.append("liability_h2_not_ready")
        h2_row = h2_by_trait.get(trait_id)
        if h2_row is None:
            phase1_issues.append("h2_qc_not_run" if args.h2 else "h2_table_not_supplied")
            h2_verdict = "UNASSESSED"
            h2_reason = "UNASSESSED"
        else:
            h2_verdict = h2_row["verdict"]
            h2_reason = h2_row.get("qc_reason", "") or "UNSPECIFIED"
            if h2_verdict != "PASS":
                phase1_issues.append("h2_qc_not_pass")
        phase1_pass = not phase1_issues

        if phase1_pass:
            stage = "PHASE1_PASS"
        elif liability_h2_ready:
            stage = "LIABILITY_H2_READY"
        elif ldsc_ready:
            stage = "LDSC_READY"
        elif harmonization_ready:
            stage = "HARMONIZATION_READY"
        elif source_verified:
            stage = "SOURCE_VERIFIED"
        else:
            stage = "SOURCE_PENDING"

        rows.append({
            "panel_version": trait["panel_version"],
            "manifest_sha256": panel_hash,
            "trait_id": trait_id,
            "label": trait["label"],
            "domain": trait["domain"],
            "type": trait["type"],
            "declared_source_status": trait["source_status"],
            "schema_status": schema_status,
            "source_verified": source_verified,
            "harmonization_ready": harmonization_ready,
            "ldsc_ready": ldsc_ready,
            "liability_h2_ready": liability_h2_ready,
            "phase1_pass": phase1_pass,
            "readiness_stage": stage,
            "h2_verdict": h2_verdict,
            "h2_qc_reason": h2_reason,
            "source_issues": ";".join(source_issues) if source_issues else "none",
            "harmonization_issues": ";".join(harmonization_issues) if harmonization_issues else "none",
            "ldsc_issues": ";".join(ldsc_issues) if ldsc_issues else "none",
            "liability_h2_issues": ";".join(liability_issues) if liability_issues else "none",
            "phase1_issues": ";".join(phase1_issues) if phase1_issues else "none",
        })

    audit = pd.DataFrame(rows)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    audit.to_csv(args.out, sep="\t", index=False)
    declaration_mismatches = audit[
        (audit["declared_source_status"] == "SOURCE_VERIFIED") & ~audit["source_verified"]
    ]
    print(f"Wrote readiness ledger: {args.out}")
    for stage in ["source_verified", "harmonization_ready", "ldsc_ready", "liability_h2_ready", "phase1_pass"]:
        print(f"{stage}: {int(audit[stage].sum())} / {len(audit)}")
    print(f"SOURCE_VERIFIED declarations without evidence: {len(declaration_mismatches)}")
    if args.strict and len(declaration_mismatches):
        bad = ", ".join(declaration_mismatches["trait_id"].tolist())
        raise SystemExit(f"ERROR: source-verification declaration mismatch: {bad}")


if __name__ == "__main__":
    main()
