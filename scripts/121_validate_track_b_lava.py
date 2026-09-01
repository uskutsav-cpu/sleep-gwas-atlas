#!/usr/bin/env python3
"""Validate the exact Track B LAVA univariate, bivariate, and conditional families."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("track_b_lava_contract", ROOT / "scripts/119_track_b_lava_contract.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load Track B LAVA contract")
contract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contract)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty Track B LAVA result: {path.relative_to(ROOT)}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def probability(value: str) -> float:
    try:
        number = float(value)
    except ValueError:
        fail(f"invalid probability: {value}")
    if not math.isfinite(number) or not 0 <= number <= 1:
        fail(f"probability outside [0,1]: {value}")
    return number


def bh(values: list[float]) -> list[float]:
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [0.0] * count
    running = 1.0
    for rank0 in range(count - 1, -1, -1):
        index = order[rank0]
        running = min(running, min(1.0, values[index] * count / (rank0 + 1)))
        adjusted[index] = running
    return adjusted


def boolean(value: str) -> bool:
    if value.lower() not in {"true", "false"}:
        fail(f"invalid logical value: {value}")
    return value.lower() == "true"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seal-results", action="store_true")
    args = parser.parse_args()
    policy = contract.load_policy()
    _, input_info = read_tsv(ROOT / "results/track_b/lava_input_info.tsv")
    _, pairs = read_tsv(ROOT / "results/track_b/lava_pair_manifest.tsv")
    _, conditioners = read_tsv(ROOT / "results/track_b/local_conditional_manifest.tsv")
    traits = [r["phenotype"] for r in input_info]
    if traits != policy["trait_order"] or [r["pair_id"] for r in pairs] != policy["pair_order"]:
        fail("Track B LAVA result inputs differ from policy order")

    locus_rows = [line.split() for line in (ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile").read_text(encoding="utf-8").splitlines() if line]
    loci = [r[0] for r in locus_rows[1:]]
    if len(loci) != policy["expected_loci"] or len(set(loci)) != len(loci):
        fail("LAVA locus family drifted")

    status_fields, status = read_tsv(contract.RESULTS[0])
    univ_fields, univ = read_tsv(contract.RESULTS[1])
    bivar_fields, bivar = read_tsv(contract.RESULTS[2])
    conditional_fields, conditional = read_tsv(contract.RESULTS[3])
    requirements = [
        (status_fields, {"LOC", "status", "analysis_fingerprint", "lava_version", "reference_prefix"}, "status"),
        (univ_fields, {"LOC", "phen", "p", "analysis_status", "univariate_test_family_n", "univariate_p_threshold", "p_bonferroni", "p_fdr"}, "univariate"),
        (bivar_fields, {"LOC", "pair_id", "trait1", "trait2", "rho", "p", "analysis_status", "bivariate_test_family_n", "p_fdr", "fdr_significant"}, "bivariate"),
        (conditional_fields, {"LOC", "pair_id", "trait1", "trait2", "covariates", "pcor", "p", "analysis_status", "conditional_test_family_n", "p_fdr", "fdr_significant"}, "conditional"),
    ]
    for fields, required, label in requirements:
        if missing := required.difference(fields):
            fail(f"{label} result lacks fields: {sorted(missing)}")

    fingerprint = contract.run_fingerprint()
    if len(status) != len(loci) or [r["LOC"] for r in status] != loci:
        fail("locus status does not cover the ordered 2,495-locus family")
    if {r["analysis_fingerprint"] for r in status} != {fingerprint}:
        fail("locus status fingerprint differs from current Track B contract")
    if {r["lava_version"] for r in status} != {policy["lava_version"]}:
        fail("locus status has wrong LAVA version")
    failed_loci = sum(r["status"] != "PROCESSED" for r in status)
    if failed_loci / len(status) > policy["maximum_locus_failure_fraction"]:
        fail("locus failure fraction exceeds policy")

    expected_univ = {(locus, trait) for locus in loci for trait in traits}
    observed_univ = {(r["LOC"], r["phen"]) for r in univ}
    if len(univ) != policy["planned_univariate_tests"] or observed_univ != expected_univ or len(observed_univ) != len(univ):
        fail("univariate result is not the exact 2,495 x 8 family")
    if {int(r["univariate_test_family_n"]) for r in univ} != {policy["planned_univariate_tests"]}:
        fail("univariate family-size provenance drifted")
    threshold = float(policy["univariate_p_threshold"])
    tested_univ = [r for r in univ if r["analysis_status"] == "TESTED"]
    if (len(univ) - len(tested_univ)) / len(univ) > policy["maximum_univariate_untested_fraction"]:
        fail("univariate untested fraction exceeds policy")
    for item in tested_univ:
        p = probability(item["p"])
        if not math.isclose(float(item["p_bonferroni"]), min(1.0, p * len(univ)), rel_tol=1e-8, abs_tol=1e-12):
            fail(f"incorrect univariate Bonferroni value: {item['LOC']}/{item['phen']}")

    significant_traits: dict[str, set[str]] = {locus: set() for locus in loci}
    for item in tested_univ:
        if probability(item["p"]) <= threshold:
            significant_traits[item["LOC"]].add(item["phen"])
    pair_defs = {(r["pair_id"], r["trait1"], r["trait2"]) for r in pairs}
    expected_bivar = {
        (locus, pair_id, first, second)
        for locus, significant in significant_traits.items()
        for pair_id, first, second in pair_defs if first in significant and second in significant
    }
    observed_bivar = {(r["LOC"], r["pair_id"], r["trait1"], r["trait2"]) for r in bivar}
    if observed_bivar != expected_bivar or len(observed_bivar) != len(bivar):
        fail("bivariate result differs from the locally eligible frozen three-pair family")
    if bivar and {int(r["bivariate_test_family_n"]) for r in bivar} != {len(bivar)}:
        fail("bivariate test-family size drifted")
    tested_bivar = [r for r in bivar if r["analysis_status"] == "TESTED"]
    if bivar and (len(bivar) - len(tested_bivar)) / len(bivar) > policy["maximum_bivariate_failure_fraction"]:
        fail("bivariate failure fraction exceeds policy")
    adjusted = bh([probability(r["p"]) for r in tested_bivar])
    for item, wanted in zip(tested_bivar, adjusted):
        if not math.isclose(float(item["p_fdr"]), wanted, rel_tol=1e-8, abs_tol=1e-12):
            fail(f"incorrect bivariate FDR: {item['LOC']}/{item['pair_id']}")
        if boolean(item["fdr_significant"]) != (wanted <= policy["bivariate_fdr_alpha"]):
            fail("bivariate significance label differs from FDR")

    expected_conditional = {
        (r["LOC"], r["pair_id"], r["trait1"], r["trait2"])
        for r in tested_bivar if boolean(r["fdr_significant"]) and r["pair_id"] in {"A", "B"}
    }
    observed_conditional = {(r["LOC"], r["pair_id"], r["trait1"], r["trait2"]) for r in conditional}
    if observed_conditional != expected_conditional or len(observed_conditional) != len(conditional):
        fail("conditional result does not exactly cover FDR-supported A/B loci")
    covariates = {
        pair_id: [r["covariate"] for r in conditioners if r["pair_id"] == pair_id and r["covariate"] != "NONE"]
        for pair_id in ("A", "B")
    }
    univ_index = {(r["LOC"], r["phen"]): r for r in univ}
    for item in conditional:
        expected_covariates = covariates[item["pair_id"]]
        if item["covariates"].split(";") != expected_covariates:
            fail("conditional covariate set drifted")
        eligible = True
        for covariate in expected_covariates:
            source = univ_index[(item["LOC"], covariate)]
            eligible &= source["analysis_status"] == "TESTED" and probability(source["p"]) <= threshold
        if not eligible and item["analysis_status"] != "CONDITIONER_LOCAL_H2_INELIGIBLE":
            fail("ineligible conditioner locus was analyzed")
        if eligible and item["analysis_status"] not in {"TESTED", "CONDITIONAL_FAILED"}:
            fail("eligible conditional locus has invalid terminal status")
    tested_conditional = [r for r in conditional if r["analysis_status"] == "TESTED"]
    failed_conditional = [r for r in conditional if r["analysis_status"] == "CONDITIONAL_FAILED"]
    denominator = len(tested_conditional) + len(failed_conditional)
    if denominator and len(failed_conditional) / denominator > policy["maximum_conditional_failure_fraction"]:
        fail("conditional failure fraction exceeds policy")
    conditional_adjusted = bh([probability(r["p"]) for r in tested_conditional])
    for item, wanted in zip(tested_conditional, conditional_adjusted):
        if int(item["conditional_test_family_n"]) != len(tested_conditional):
            fail("conditional test-family size drifted")
        if not math.isclose(float(item["p_fdr"]), wanted, rel_tol=1e-8, abs_tol=1e-12):
            fail(f"incorrect conditional FDR: {item['LOC']}/{item['pair_id']}")

    if args.seal_results:
        contract.seal_results()
    else:
        contract.validate_results()
    print(f"TRACK_B_LAVA_VALIDATED loci={len(status)} univariate={len(univ)} bivariate={len(bivar)} conditional={len(conditional)}")


if __name__ == "__main__":
    main()
