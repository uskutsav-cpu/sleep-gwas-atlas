#!/usr/bin/env python3
"""Validate that canonical LAVA outputs cover the preregistered test families."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import lava_contract


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def finite_probability(value: str) -> float:
    try:
        number = float(value)
    except ValueError:
        fail(f"invalid probability: {value}")
    if not math.isfinite(number) or not 0 <= number <= 1:
        fail(f"probability outside [0,1]: {value}")
    return number


def bh(values: list[float]) -> list[float]:
    count = len(values)
    ranked = sorted(range(count), key=values.__getitem__)
    adjusted = [0.0] * count
    running = 1.0
    for zero_rank in range(count - 1, -1, -1):
        index = ranked[zero_rank]
        candidate = min(1.0, values[index] * count / (zero_rank + 1))
        running = min(running, candidate)
        adjusted[index] = running
    return adjusted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--seal-results", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads((root / "config/lava_analysis_policy.json").read_text(encoding="utf-8"))
    _, panel = read_tsv(root / "config/analysis_panel.tsv")
    traits = [row["trait_id"] for row in panel]
    if len(traits) != 45 or len(set(traits)) != 45:
        fail("panel is not the exact 45-trait set")

    locus_path = root / policy["locus_definition"]
    if not locus_path.is_file():
        fail(f"missing locus definition: {policy['locus_definition']}")
    if hashlib.sha256(locus_path.read_bytes()).hexdigest() != policy["locus_definition_sha256"]:
        fail("locus-definition checksum mismatch")
    locus_lines = [line.split() for line in locus_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not locus_lines or locus_lines[0][:4] != ["LOC", "CHR", "START", "STOP"]:
        fail("locus definition has an invalid header")
    locus_ids = [row[0] for row in locus_lines[1:]]
    if len(locus_ids) != policy["expected_loci"] or len(set(locus_ids)) != len(locus_ids):
        fail("locus definition is not the exact 2,495-locus family")

    status_fields, status = read_tsv(root / "results/tables/lava_locus_status.tsv")
    univ_fields, univ = read_tsv(root / "results/tables/lava_univariate.tsv")
    bivar_fields, bivar = read_tsv(root / "results/tables/lava_bivariate.tsv")
    _, pairs = read_tsv(root / "results/tables/lava_pair_manifest.tsv")
    required_status = {"LOC", "status", "univariate_tested", "eligible_bivariate_pairs", "analysis_fingerprint", "lava_version", "reference_prefix"}
    required_univ = {"LOC", "phen", "h2.obs", "h2.latent", "p", "analysis_status", "univariate_test_family_n", "univariate_p_threshold", "p_bonferroni", "p_fdr"}
    required_bivar = {"LOC", "sleep_trait", "non_sleep_trait", "analysis_tier", "rho", "p", "analysis_status", "bivariate_test_family_n", "p_fdr", "fdr_significant"}
    for label, fields, required in [
        ("locus status", status_fields, required_status),
        ("univariate", univ_fields, required_univ),
        ("bivariate", bivar_fields, required_bivar),
    ]:
        missing = required.difference(fields)
        if missing:
            fail(f"{label} output lacks columns: {sorted(missing)}")

    if len(status) != len(locus_ids) or [row["LOC"] for row in status] != locus_ids:
        fail("locus status does not cover the ordered 2,495-locus family")
    if {row["lava_version"] for row in status} != {policy["lava_version"]}:
        fail("locus status has the wrong LAVA version")
    if {row["reference_prefix"] for row in status} != {policy["reference_prefix"]}:
        fail("locus status has the wrong reference prefix")
    fingerprints = {row["analysis_fingerprint"] for row in status}
    expected_fingerprint = lava_contract.run_fingerprint(root)
    if fingerprints != {expected_fingerprint}:
        fail("locus checkpoints do not share the current immutable run fingerprint")
    failed_loci = sum(row["status"] != "PROCESSED" for row in status)
    if failed_loci / len(status) > policy["max_locus_failure_fraction"]:
        fail(f"locus failure fraction exceeds policy: {failed_loci}/{len(status)}")

    expected_univ = {(locus, trait) for locus in locus_ids for trait in traits}
    observed_univ = {(row["LOC"], row["phen"]) for row in univ}
    if len(univ) != policy["planned_univariate_tests"] or observed_univ != expected_univ:
        fail("univariate output is not the exact 2,495 x 45 planned family")
    if len(observed_univ) != len(univ):
        fail("duplicate univariate locus-trait rows")
    if {int(row["univariate_test_family_n"]) for row in univ} != {policy["planned_univariate_tests"]}:
        fail("univariate family-size provenance is inconsistent")
    threshold = policy["univariate_p_threshold"]
    if any(not math.isclose(float(row["univariate_p_threshold"]), threshold, rel_tol=1e-12) for row in univ):
        fail("univariate threshold differs from policy")
    tested_univ = [row for row in univ if row["analysis_status"] == "TESTED"]
    if (len(univ) - len(tested_univ)) / len(univ) > policy["max_univariate_untested_fraction"]:
        fail("univariate untested fraction exceeds policy")
    for row in tested_univ:
        p_value = finite_probability(row["p"])
        if not math.isclose(float(row["p_bonferroni"]), min(1.0, p_value * len(univ)), rel_tol=1e-8, abs_tol=1e-12):
            fail(f"incorrect Bonferroni value at locus {row['LOC']}/{row['phen']}")

    pair_set = {(row["sleep_trait"], row["non_sleep_trait"]) for row in pairs}
    if len(pairs) != 396 or len(pair_set) != 396:
        fail("pair manifest is not the exact 396-pair family")
    significant_by_locus: dict[str, set[str]] = {locus: set() for locus in locus_ids}
    for row in tested_univ:
        if finite_probability(row["p"]) <= threshold:
            significant_by_locus[row["LOC"]].add(row["phen"])
    expected_bivar = {
        (locus, sleep, non_sleep)
        for locus, significant in significant_by_locus.items()
        for sleep, non_sleep in pair_set
        if sleep in significant and non_sleep in significant
    }
    observed_bivar = {(row["LOC"], row["sleep_trait"], row["non_sleep_trait"]) for row in bivar}
    if observed_bivar != expected_bivar or len(observed_bivar) != len(bivar):
        fail("bivariate output differs from the locally eligible locked pair family")
    if bivar and {int(row["bivariate_test_family_n"]) for row in bivar} != {len(bivar)}:
        fail("bivariate test-family provenance is inconsistent")
    tested_bivar = [row for row in bivar if row["analysis_status"] == "TESTED"]
    if bivar and (len(bivar) - len(tested_bivar)) / len(bivar) > policy["max_bivariate_failure_fraction"]:
        fail("bivariate failure fraction exceeds policy")
    p_values = [finite_probability(row["p"]) for row in tested_bivar]
    expected_fdr = bh(p_values)
    for row, adjusted in zip(tested_bivar, expected_fdr):
        if not math.isclose(float(row["p_fdr"]), adjusted, rel_tol=1e-8, abs_tol=1e-12):
            fail(f"incorrect bivariate FDR at locus {row['LOC']}")

    if args.seal_results:
        lava_contract.seal_results(root, policy)
    else:
        lava_contract.validate_results(root, policy)

    if not args.quiet:
        print(f"Validated LAVA: {len(status)} loci, {len(tested_univ)}/{len(univ)} univariate tests, {len(tested_bivar)}/{len(bivar)} eligible bivariate tests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
