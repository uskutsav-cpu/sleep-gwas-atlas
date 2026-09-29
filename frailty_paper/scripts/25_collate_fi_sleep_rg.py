#!/usr/bin/env python3
"""Collate the 12 fresh FI rg estimates with the locked all-396 family q-values.

The frozen atlas extract is accepted only when all twelve fresh raw p-values
match the corresponding frozen values. This prevents accidentally carrying an
old q-value onto a changed analysis result or adjusting only the 12 displayed
pairs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_SLEEP = {
    "insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype",
    "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency",
    "accel_sleep_duration", "sleep_timing",
}


def collate(pair_dir: Path, frozen_path: Path, overlap_path: Path) -> pd.DataFrame:
    paths = sorted(pair_dir.glob("*.tsv"))
    if len(paths) != 12:
        raise ValueError(f"expected 12 pair tables, found {len(paths)} in {pair_dir}")
    fresh = pd.concat((pd.read_csv(path, sep="\t") for path in paths), ignore_index=True)
    if set(fresh["sleep_trait"]) != EXPECTED_SLEEP or set(fresh["disease_trait"]) != {"frailty"}:
        raise ValueError("fresh rg tables do not contain exactly the locked 12 sleep x FI pairs")
    if fresh.duplicated(["sleep_trait", "disease_trait"]).any():
        raise ValueError("duplicate fresh sleep x FI rg rows")

    frozen = pd.read_csv(frozen_path, sep="\t").rename(columns={"non_sleep_trait": "disease_trait"})
    needed_frozen = {"sleep_trait", "disease_trait", "global_rg_p", "global_rg_fdr_all_396"}
    if not needed_frozen.issubset(frozen.columns):
        raise ValueError(f"frozen atlas extract is missing {sorted(needed_frozen - set(frozen.columns))}")
    joined = fresh.merge(frozen, on=["sleep_trait", "disease_trait"], how="outer", validate="one_to_one", indicator=True)
    if not joined["_merge"].eq("both").all():
        raise ValueError("fresh and frozen tables do not contain the same pair family")
    if not np.allclose(joined["p"], joined["global_rg_p"], rtol=1e-10, atol=0.0):
        bad = joined.loc[~np.isclose(joined["p"], joined["global_rg_p"], rtol=1e-10, atol=0.0), "sleep_trait"].tolist()
        raise ValueError(f"fresh raw p-values differ from the frozen family for: {bad}")

    overlap = pd.read_csv(overlap_path, sep="\t")
    overlap = overlap.loc[overlap["frailty_endpoint"].eq("frailty"), [
        "sleep_trait", "cohort_overlap_status", "exact_participant_overlap"
    ]]
    if overlap["sleep_trait"].duplicated().any() or set(overlap["sleep_trait"]) != EXPECTED_SLEEP:
        raise ValueError("overlap ledger does not have one row for every primary sleep x FI pair")
    joined = joined.merge(overlap, on="sleep_trait", validate="one_to_one")
    joined["analysis_family"] = "primary_sleep_x_frailty_index"
    joined["fdr_family_denominator"] = 396
    joined["fdr_source"] = "frozen_atlas_all_396; fresh raw p-values verified equal"
    joined["interpretation"] = "global genetic correlation; not causal"
    columns = [
        "sleep_trait", "disease_trait", "analysis_family", "rg", "se", "z", "p",
        "global_rg_fdr_all_396", "fdr_family_denominator", "h2_obs", "h2_obs_se",
        "h2_int", "h2_int_se", "gcov_int", "gcov_int_se", "cohort_overlap_status",
        "exact_participant_overlap", "fdr_source", "input_log", "interpretation",
    ]
    result = joined[columns].sort_values("sleep_trait").reset_index(drop=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-dir", type=Path, default=Path("frailty_paper/results/frailty_v1/rg_sleep_fi"))
    parser.add_argument("--frozen", type=Path, default=Path("frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv"))
    parser.add_argument("--overlap", type=Path, default=Path("frailty_paper/analysis/sample_overlap_assessment.tsv"))
    parser.add_argument("--out", type=Path, default=Path("frailty_paper/results/frailty_v1/global_rg_sleep_frailty.tsv"))
    args = parser.parse_args()
    result = collate(args.pair_dir, args.frozen, args.overlap)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False, float_format="%.12g")
    print(f"PASS: {len(result)} FI-primary rg pairs; all fresh p-values match frozen all-396 family")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
