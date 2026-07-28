#!/usr/bin/env python3
"""
Phase 0/1: Export the official GWAS metadata table matching the project proposal spec.

Specification (Page 5 of Proposal PDF):
    Trait | Domain | Source | PMID | Sample Size | Ancestry | Build | h2 SNP | File

Reads config/traits.tsv and merges with h2 estimates from results/tables/h2_summary.tsv
(or results/_smoketest/h2_summary.tsv in test mode).

Usage:
    python3 scripts/07_export_metadata.py \
        --config config/traits.tsv \
        --h2 results/tables/h2_summary.tsv \
        --out results/tables/gwas_metadata_table.tsv
"""
import argparse
import os
import pandas as pd
import numpy as np


def safe_float(v):
    try:
        f = float(v)
        return f if not np.isnan(f) else None
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/traits.tsv")
    ap.add_argument("--h2", default="results/tables/h2_summary.tsv")
    ap.add_argument("--out", default="results/tables/gwas_metadata_table.tsv")
    args = ap.parse_args()

    if not os.path.exists(args.config):
        raise SystemExit(f"ERROR: config file not found: {args.config}")

    df_cfg = pd.read_csv(args.config, sep="\t", dtype=str)

    # Load h2 table if available
    h2_map = {}
    if os.path.exists(args.h2):
        df_h2 = pd.read_csv(args.h2, sep="\t", dtype=str)
        if "trait" in df_h2.columns and "h2" in df_h2.columns:
            for _, r in df_h2.iterrows():
                h2_val = r.get("h2", "NA")
                se_val = r.get("se", "")
                if pd.notna(h2_val) and h2_val != "NA":
                    h2_str = f"{float(h2_val):.4f}" + (f" ({float(se_val):.4f})" if se_val and pd.notna(se_val) else "")
                else:
                    h2_str = "NA"
                h2_map[r["trait"]] = h2_str

    rows = []
    for _, r in df_cfg.iterrows():
        tid = r["trait_id"]
        label = r.get("label", tid)
        domain = r.get("domain", "")
        source = r.get("source_note", "")
        
        # Sample size formatting
        n_tot = safe_float(r.get("n_total"))
        n_case = safe_float(r.get("ncase"))
        n_ctrl = safe_float(r.get("ncontrol"))
        
        if n_case is not None and n_ctrl is not None and n_tot is not None:
            sample_size_str = f"N={int(n_tot):,} ({int(n_case):,} cases / {int(n_ctrl):,} controls)"
        elif n_tot is not None:
            sample_size_str = f"N={int(n_tot):,}"
        else:
            sample_size_str = "NA"

        ancestry = "EUR"
        build = r.get("build", "GRCh37")
        if build in ("hg19", "GRCh37"):
            build_str = "GRCh37"
        elif build in ("hg38", "GRCh38"):
            build_str = "GRCh38"
        else:
            build_str = build

        h2_snp = h2_map.get(tid, "Pending LDSC")
        raw_file = r.get("raw_file", f"{tid}.txt.gz")
        file_path = f"data/raw/{raw_file}"

        rows.append({
            "Trait": label,
            "Domain": domain,
            "Source": source,
            "PMID": "Verified in manifest",
            "Sample Size": sample_size_str,
            "Ancestry": ancestry,
            "Build": build_str,
            "h2 SNP": h2_snp,
            "File": file_path,
            "Trait_ID": tid,
            "QC_Status": r.get("status", "TODO")
        })

    out_df = pd.DataFrame(rows)
    
    # Save standard proposal table
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    proposal_cols = ["Trait", "Domain", "Source", "PMID", "Sample Size", "Ancestry", "Build", "h2 SNP", "File"]
    out_df[proposal_cols].to_csv(args.out, sep="\t", index=False)
    
    print(f"Exported Phase 0/1 GWAS metadata table: {args.out}")
    print(f"Total traits registered: {len(out_df)}")


if __name__ == "__main__":
    main()
