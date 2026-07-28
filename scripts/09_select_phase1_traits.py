#!/usr/bin/env python3
"""Create the auditable trait-inclusion table used by Phase 1 rg.

A trait must be curated in the registry *and* pass the h2 QC gate. This avoids
the old manual instruction to overwrite ``status`` with ``PASS`` and preserves
the distinction between source curation and an LDSC result.
"""
import argparse
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/traits.tsv")
    parser.add_argument("--h2", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    config = pd.read_csv(args.config, sep="\t", dtype=str)
    h2 = pd.read_csv(args.h2, sep="\t", dtype=str)
    if config["trait_id"].duplicated().any() or h2["trait"].duplicated().any():
        raise SystemExit("ERROR: trait IDs must be unique in config and h2 summary")
    if not {"trait", "verdict"}.issubset(h2.columns):
        raise SystemExit("ERROR: h2 summary must include trait and verdict columns")

    merged = config.merge(h2[["trait", "verdict", "qc_reason"]], how="left", left_on="trait_id", right_on="trait")
    curated = merged["status"].fillna("").str.upper().eq("CURATED")
    h2_pass = merged["verdict"].eq("PASS")
    merged["include_phase1"] = curated & h2_pass
    merged["inclusion_reason"] = "included"
    merged.loc[~curated, "inclusion_reason"] = "source_not_curated"
    merged.loc[curated & ~h2_pass, "inclusion_reason"] = "h2_qc_not_pass"
    output = merged[["trait_id", "label", "domain", "type", "status", "verdict", "qc_reason", "include_phase1", "inclusion_reason"]]

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    output.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote Phase 1 inclusion table: {args.out}")
    print(f"Included: {int(output['include_phase1'].sum())} / {len(output)}")


if __name__ == "__main__":
    main()
