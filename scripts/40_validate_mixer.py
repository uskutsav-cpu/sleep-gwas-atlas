#!/usr/bin/env python3
"""Validate complete univariate-first and eligible-pair MiXeR outputs."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def finite(row: dict[str, str], columns: set[str]) -> None:
    for column in columns:
        try:
            value = float(row[column])
        except (KeyError, ValueError):
            fail(f"missing/non-numeric {column}")
        if not math.isfinite(value):
            fail(f"non-finite {column}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads((root / "config/mixer_analysis_policy.json").read_text(encoding="utf-8"))
    _, panel = read_tsv(root / "config/analysis_panel.tsv")
    traits = [row["trait_id"] for row in panel]
    if len(traits) != 45 or len(set(traits)) != 45:
        fail("panel is not the exact 45-trait set")
    univ_fields, univ = read_tsv(root / "results/tables/mixer_univariate.tsv")
    bivar_fields, bivar = read_tsv(root / "results/tables/mixer_bivariate.tsv")
    required_univ = {
        "trait_id", "mixer_version", "container_digest", "reference_commit", "input_scope",
        "replicates", "pi_mean", "pi_std", "sig2_beta_mean", "sig2_beta_std",
        "sig2_zero_mean", "sig2_zero_std", "h2_mean", "h2_std", "n_causal_p9_mean",
        "n_causal_p9_std", "AIC", "BIC", "power_interpretation", "bivariate_eligibility",
    }
    required_bivar = {
        "sleep_trait", "non_sleep_trait", "mixer_version", "container_digest",
        "reference_commit", "replicates", "dice_mean", "dice_std", "pi1_mean", "pi1_std",
        "pi2_mean", "pi2_std", "pi12_mean", "pi12_std", "n_shared_p9_mean",
        "n_shared_p9_std", "rho_beta_mean", "rho_beta_std", "rg_mean", "rg_std",
        "fraction_concordant_mean", "fraction_concordant_std", "best_vs_min_AIC",
        "best_vs_min_BIC", "best_vs_max_AIC", "best_vs_max_BIC",
    }
    if missing := required_univ.difference(univ_fields):
        fail(f"univariate table lacks columns: {sorted(missing)}")
    if missing := required_bivar.difference(bivar_fields):
        fail(f"bivariate table lacks columns: {sorted(missing)}")
    if len(univ) != 45 or [row["trait_id"] for row in univ] != traits:
        fail("univariate table is not the ordered 45-trait family")
    common_expected = {
        "mixer_version": policy["mixer_release"],
        "container_digest": policy["container_amd64_manifest_digest"],
        "reference_commit": policy["mixer_reference_commit_at_lock"],
        "replicates": str(policy["fit_replicates"]),
    }
    numeric_univ = {
        "pi_mean", "pi_std", "sig2_beta_mean", "sig2_beta_std", "sig2_zero_mean",
        "sig2_zero_std", "h2_mean", "h2_std", "n_causal_p9_mean", "n_causal_p9_std",
        "AIC", "BIC",
    }
    eligible = set()
    for row in univ:
        for key, expected in common_expected.items():
            if row[key] != expected:
                fail(f"wrong {key} for {row['trait_id']}")
        if row["input_scope"] != "FULL_POST_QC_AUTOSOMAL":
            fail(f"non-full MiXeR input for {row['trait_id']}")
        finite(row, numeric_univ)
        expected_eligibility = "ELIGIBLE" if float(row["AIC"]) > policy["univariate_aic_threshold"] else "INELIGIBLE_LOW_POWER"
        if row["bivariate_eligibility"] != expected_eligibility:
            fail(f"incorrect eligibility for {row['trait_id']}")
        expected_power = "SUPPORTED_AIC_BIC" if float(row["AIC"]) > 0 and float(row["BIC"]) > 0 else (
            "BORDERLINE_POWER_AIC_ONLY" if float(row["AIC"]) > 0 else "INSUFFICIENT_POWER"
        )
        if row["power_interpretation"] != expected_power:
            fail(f"incorrect power interpretation for {row['trait_id']}")
        if not 0 < float(row["pi_mean"]) <= 1:
            fail(f"polygenicity outside (0,1] for {row['trait_id']}")
        if float(row["sig2_beta_mean"]) <= 0 or float(row["sig2_zero_mean"]) < 0:
            fail(f"invalid variance component for {row['trait_id']}")
        if float(row["h2_mean"]) <= 0 or float(row["n_causal_p9_mean"]) <= 0:
            fail(f"non-positive h2/causal count for {row['trait_id']}")
        for column in ("pi_std", "sig2_beta_std", "sig2_zero_std", "h2_std", "n_causal_p9_std"):
            if float(row[column]) < 0:
                fail(f"negative uncertainty {column} for {row['trait_id']}")
        if expected_eligibility == "ELIGIBLE":
            eligible.add(row["trait_id"])

    expected_pairs = {
        (first["trait_id"], second["trait_id"])
        for first in panel if first["domain"] == "sleep" and first["trait_id"] in eligible
        for second in panel if second["domain"] != "sleep" and second["trait_id"] in eligible
    }
    observed_pairs = {(row["sleep_trait"], row["non_sleep_trait"]) for row in bivar}
    if len(observed_pairs) != len(bivar) or observed_pairs != expected_pairs:
        fail("bivariate table differs from the univariate-eligible locked pair family")
    numeric_bivar = required_bivar.difference({
        "sleep_trait", "non_sleep_trait", "mixer_version", "container_digest",
        "reference_commit", "replicates",
    })
    for row in bivar:
        for key, expected in common_expected.items():
            if row[key] != expected:
                fail(f"wrong {key} for {row['sleep_trait']}/{row['non_sleep_trait']}")
        finite(row, numeric_bivar)
        if not 0 <= float(row["dice_mean"]) <= 1:
            fail("Dice overlap outside [0,1]")
        if not 0 <= float(row["fraction_concordant_mean"]) <= 1:
            fail("concordance fraction outside [0,1]")
        for column in ("pi1_mean", "pi2_mean", "pi12_mean", "n_shared_p9_mean"):
            if float(row[column]) < 0:
                fail(f"negative overlap parameter {column}")
        for column in ("rho_beta_mean", "rg_mean"):
            if not -1 <= float(row[column]) <= 1:
                fail(f"correlation {column} outside [-1,1]")
        for column in (
            "dice_std", "pi1_std", "pi2_std", "pi12_std", "n_shared_p9_std",
            "rho_beta_std", "rg_std", "fraction_concordant_std",
        ):
            if float(row[column]) < 0:
                fail(f"negative uncertainty {column}")
    if not args.quiet:
        print(f"Validated MiXeR: 45 univariate models, {len(eligible)} eligible traits, {len(bivar)} eligible sleep-by-non-sleep pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
