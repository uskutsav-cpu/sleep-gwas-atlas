#!/usr/bin/env python3
"""Export the Phase 0 metadata table without inventing provenance fields.

The proposal requires Trait, Domain, Source, PMID, Sample Size, Ancestry,
Build, h2 SNP, and File. ``config/analysis_panel.tsv`` is the locked panel, not
proof that a PMID, ancestry subset, or population prevalence has been checked.
Consequently missing curation values are emitted as ``UNRESOLVED`` rather than
the previous, incorrect blanket statement that every PMID was verified.
"""
import argparse
import os

import numpy as np
import pandas as pd


def safe_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def format_build(value):
    build = str(value).strip()
    return {"hg19": "GRCh37", "GRCh37": "GRCh37", "hg38": "GRCh38", "GRCh38": "GRCh38"}.get(build, build or "UNRESOLVED")


def format_h2(row):
    value = row.get("h2")
    standard_error = row.get("se")
    scale = row.get("scale")
    if pd.isna(value):
        return "Pending LDSC"
    suffix = f" ({scale} scale)" if pd.notna(scale) else ""
    if pd.notna(standard_error):
        return f"{float(value):.4f} ({float(standard_error):.4f}){suffix}"
    return f"{float(value):.4f}{suffix}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--h2", default="results/tables/h2_summary.tsv")
    parser.add_argument("--out", default="results/tables/gwas_metadata_table.tsv")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        raise SystemExit(f"ERROR: config file not found: {args.config}")
    config = pd.read_csv(args.config, sep="\t", dtype=str)
    if config["trait_id"].duplicated().any():
        raise SystemExit("ERROR: config contains duplicate trait_id values")

    h2_by_trait = {}
    if os.path.exists(args.h2):
        h2 = pd.read_csv(args.h2, sep="\t")
        if "trait" in h2:
            h2_by_trait = {row["trait"]: row for _, row in h2.iterrows()}

    rows = []
    for _, metadata in config.iterrows():
        trait_id = metadata["trait_id"]
        n_total = safe_float(metadata.get("n_total"))
        n_case = safe_float(metadata.get("ncase"))
        n_control = safe_float(metadata.get("ncontrol"))
        if n_case is not None and n_control is not None and n_total is not None:
            sample_size = f"N={int(n_total):,} ({int(n_case):,} cases / {int(n_control):,} controls)"
        elif n_total is not None:
            sample_size = f"N={int(n_total):,}"
        else:
            sample_size = "UNRESOLVED"

        pmid = metadata.get("pmid", "")
        pmid_status = metadata.get("pmid_status", "")
        if not pmid or str(pmid).strip().upper() in {"NA", "UNKNOWN", "NAN"}:
            pmid_display = "UNRESOLVED - see docs/phase0_trait_manifest.csv"
        elif pmid_status:
            pmid_display = f"{pmid} ({pmid_status})"
        else:
            pmid_display = str(pmid)

        ancestry = metadata.get("ancestry", "")
        if not ancestry or str(ancestry).strip().upper() in {"NA", "UNKNOWN", "NAN"}:
            ancestry = "EUR target - source subset verification pending"

        h2_row = h2_by_trait.get(trait_id)
        h2_display = format_h2(h2_row) if h2_row is not None else "Pending LDSC"
        rows.append({
            "Trait": metadata.get("label", trait_id),
            "Domain": metadata.get("domain", ""),
            "Source": metadata.get("source_note", ""),
            "PMID": pmid_display,
            "Sample Size": sample_size,
            "Ancestry": ancestry,
            "Build": format_build(metadata.get("build", "")),
            "h2 SNP": h2_display,
            "File": f"data/raw/{metadata.get('raw_file', trait_id + '.txt.gz')}",
            "Trait_ID": trait_id,
            "Source_Status": metadata.get("source_status", "SOURCE_PENDING"),
        })

    output = pd.DataFrame(rows)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    proposal_columns = ["Trait", "Domain", "Source", "PMID", "Sample Size", "Ancestry", "Build", "h2 SNP", "File"]
    output[proposal_columns].to_csv(args.out, sep="\t", index=False)
    print(f"Exported provenance-aware Phase 0 metadata table: {args.out}")
    print(f"Registered traits: {len(output)}; rows with unsupplied PMID: {(output['PMID'].str.startswith('UNRESOLVED')).sum()}")


if __name__ == "__main__":
    main()
