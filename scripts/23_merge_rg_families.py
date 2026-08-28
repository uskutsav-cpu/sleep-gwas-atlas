#!/usr/bin/env python3
"""Merge primary and QC-failed LDSC rg families into the locked 396 pairs.

QC-failed traits remain sensitivity-only. Their estimates are retained for the
complete atlas family but are never relabelled as Phase 1 passes.
"""
import argparse
import os

import numpy as np
import pandas as pd


PAIR = ["sleep_trait", "disease_trait"]


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def bh(values):
    values = pd.to_numeric(values, errors="coerce")
    if values.isna().any() or (~values.between(0, 1)).any():
        fail("P values must all be numeric and in [0, 1]")
    order = np.argsort(values.to_numpy(), kind="stable")
    ranked = values.to_numpy()[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1].clip(max=1.0)
    result = np.empty(len(values), dtype=float)
    result[order] = adjusted
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--h2", default="results/tables/h2_summary.tsv")
    parser.add_argument("--primary", default="results/tables/rg_primary_phase1.tsv")
    parser.add_argument(
        "--sensitivity", default="results/tables/rg_qc_failed_sensitivity.tsv"
    )
    parser.add_argument("--out", default="results/tables/rg_matrix.tsv")
    args = parser.parse_args()

    panel = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
    h2 = pd.read_csv(args.h2, sep="\t", dtype=str).fillna("")
    primary = pd.read_csv(args.primary, sep="\t")
    sensitivity = pd.read_csv(args.sensitivity, sep="\t")
    if panel["trait_id"].duplicated().any() or len(panel) != 45:
        fail("panel must contain exactly 45 unique trait IDs")
    if h2["trait"].duplicated().any() or set(h2["trait"]) != set(panel["trait_id"]):
        fail("h2 table must cover exactly the locked 45 traits")
    if not {"trait", "verdict", "qc_reason"}.issubset(h2.columns):
        fail("h2 table lacks verdict or QC-reason columns")

    for name, frame in [("primary", primary), ("sensitivity", sensitivity)]:
        missing = {"sleep_trait", "disease_trait", "rg", "se", "p"}.difference(
            frame.columns
        )
        if missing:
            fail(f"{name} rg table lacks columns: {sorted(missing)}")
        if frame.duplicated(PAIR).any():
            fail(f"{name} rg table contains duplicate pairs")

    verdict = h2.set_index("trait")["verdict"].to_dict()
    reason = h2.set_index("trait")["qc_reason"].to_dict()
    sleep = panel.loc[panel["domain"].eq("sleep"), "trait_id"].tolist()
    disease = panel.loc[~panel["domain"].eq("sleep"), "trait_id"].tolist()
    expected = {(a, b) for a in sleep for b in disease}
    expected_primary = {
        pair for pair in expected
        if verdict[pair[0]] == "PASS" and verdict[pair[1]] == "PASS"
    }
    expected_sensitivity = expected.difference(expected_primary)
    observed_primary = set(primary[PAIR].itertuples(index=False, name=None))
    observed_sensitivity = set(sensitivity[PAIR].itertuples(index=False, name=None))
    if observed_primary != expected_primary:
        fail(
            f"primary family mismatch: expected {len(expected_primary)}, "
            f"found {len(observed_primary)}"
        )
    if observed_sensitivity != expected_sensitivity:
        fail(
            f"sensitivity family mismatch: expected {len(expected_sensitivity)}, "
            f"found {len(observed_sensitivity)}"
        )

    combined = pd.concat([primary, sensitivity], ignore_index=True)
    if len(combined) != 396 or set(
        combined[PAIR].itertuples(index=False, name=None)
    ) != expected:
        fail("combined rg table is not the exact locked 396-pair family")
    combined["sleep_h2_verdict"] = combined["sleep_trait"].map(verdict)
    combined["disease_h2_verdict"] = combined["disease_trait"].map(verdict)
    combined["disease_h2_qc_reason"] = combined["disease_trait"].map(reason)
    primary_mask = (
        combined["sleep_h2_verdict"].eq("PASS")
        & combined["disease_h2_verdict"].eq("PASS")
    )
    combined["analysis_tier"] = np.where(
        primary_mask, "PRIMARY_PHASE1", "QC_FAILED_SENSITIVITY"
    )
    combined["interpretation_status"] = np.where(
        primary_mask, "PRIMARY", "EXCLUDED_FROM_PRIMARY_INFERENCE"
    )
    combined["fdr"] = bh(combined["p"])
    combined["fdr_primary_phase1"] = np.nan
    combined.loc[primary_mask, "fdr_primary_phase1"] = bh(
        combined.loc[primary_mask, "p"]
    )

    panel_order = {trait: i for i, trait in enumerate(panel["trait_id"])}
    combined["_sleep_order"] = combined["sleep_trait"].map(panel_order)
    combined["_disease_order"] = combined["disease_trait"].map(panel_order)
    combined.sort_values(["_sleep_order", "_disease_order"], inplace=True)
    combined.drop(columns=["_sleep_order", "_disease_order"], inplace=True)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    combined.to_csv(args.out, sep="\t", index=False)
    print(
        f"Wrote locked rg family: {args.out} "
        f"({int(primary_mask.sum())} primary + "
        f"{int((~primary_mask).sum())} QC-failed sensitivity = {len(combined)})"
    )


if __name__ == "__main__":
    main()
