#!/usr/bin/env python3
"""Freeze a source-independent LAVA failure matrix before new source selection."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

from build_lava_canonical_v3_not_run_causes import validate_source


ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
RUN = ROOT / "work/lava-canonical-v3-production" / RUN_ID
CONFIG = ROOT / "brain6/config/lava_family_canonical_v3.json"
CELL_CAUSES = ROOT / "brain6/results/lava/canonical_v3_not_run_cells.tsv"
CELL_RECEIPT = ROOT / "brain6/results/lava/canonical_v3_not_run_cells.provenance.json"
OUT = ROOT / "brain6/results/lava_confirmatory_rescue_plan_v1"
CAUSES = (
    "LOW_LOCAL_H2_UNDERPOWERED",
    "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS",
    "FEWER_THAN_MIN_K",
)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path: Path, rows: list[dict]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def strict_pass(row: dict[str, str], threshold: float) -> bool:
    return row["status"] == "TESTED" and row["p"] not in ("", "NA") and float(row["p"]) < threshold


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned rescue-planning output already exists")
    decision, _latest, aggregate_path, expected_causes = validate_source(RUN)
    config = json.loads(CONFIG.read_text())
    if config["analysis_id"] != decision["analysis_id"] or config["pairs"] != [
        "insomnia__adhd", "insomnia__mdd", "longsleep__scz", "longsleep__bipolar", "longsleep__parkinson"
    ]:
        raise ValueError("Frozen family/pair identities changed")
    if sha(CELL_CAUSES) != json.loads(CELL_RECEIPT.read_text())["output_sha256"]:
        raise ValueError("Existing cell-level cause receipt failed")
    if len(read(CELL_CAUSES)) != len(expected_causes) or len(expected_causes) != 3720:
        raise ValueError("Cell-level cause coverage differs from canonical aggregate")
    aggregate = read(aggregate_path)
    if len(aggregate) != 17465:
        raise ValueError("Incomplete seven-trait family")
    by_trait = defaultdict(dict)
    for row in aggregate:
        trait, locus = row["phen"], row["locus_id"]
        if locus in by_trait[trait]:
            raise ValueError(f"Duplicate canonical trait/locus cell: {trait}/{locus}")
        by_trait[trait][locus] = row
    traits = sorted(by_trait)
    loci = sorted({r["locus_id"] for r in aggregate}, key=int)
    if traits != sorted(config["trait_ids"]) or len(loci) != 2495 or any(len(by_trait[t]) != 2495 for t in traits):
        raise ValueError("Frozen 7 x 2,495 matrix is incomplete")
    threshold = config["univariate_gate"]["p_threshold_strictly_less_than"]
    if threshold != decision["multiple_testing"]["strict_p_threshold"]:
        raise ValueError("Frozen Bonferroni threshold mismatch")

    trait_rows = []
    cause_ranking = []
    for trait in traits:
        rows = list(by_trait[trait].values())
        c = Counter(r["status"] for r in rows)
        reasons = Counter(r["reason"] for r in rows if r["status"] == "NOT_RUN")
        if c["TESTED"] + c["NOT_RUN"] + c["FAILED"] != 2495 or sum(reasons.values()) != c["NOT_RUN"]:
            raise ValueError(f"Status/reason accounting failed for {trait}")
        item = {"trait_id": trait, "planned_cells": 2495, "TESTED": c["TESTED"],
                "NOT_RUN": c["NOT_RUN"], "FAILED": c["FAILED"],
                "fraction_testable_numeric": f"{c['TESTED']/2495:.8f}",
                "strict_univariate_gate_pass": sum(strict_pass(r, threshold) for r in rows),
                "low_local_h2_underpowered": reasons[CAUSES[0]],
                "fewer_than_min_k_shared_reference_variants": reasons[CAUSES[1]],
                "fewer_than_min_k_components": reasons[CAUSES[2]],
                "observed_primary_cause": "LOCAL_H2_SUPPORT_GATE" if reasons[CAUSES[0]] >= max(
                    reasons[CAUSES[1]], reasons[CAUSES[2]]) else "REFERENCE_OR_COMPONENT_SPARSE",
                "upstream_mechanism": "UNRESOLVED_FROM_CANONICAL_OUTPUT_ALONE",
                "analysis_label": "RESCUE_PLANNING_ONLY"}
        trait_rows.append(item)
        cause_ranking.append({"rank": "", "trait_id": trait, "NOT_RUN": c["NOT_RUN"],
                              "share_of_all_3720_not_run": f"{c['NOT_RUN']/3720:.8f}",
                              "best_possible_residual_if_this_trait_perfect": 3720-c["NOT_RUN"],
                              "priority_basis": "OBSERVED_NOT_RUN_COUNT_ONLY"})
    cause_ranking.sort(key=lambda r: (-r["NOT_RUN"], r["trait_id"]))
    for index, row in enumerate(cause_ranking, 1):
        row["rank"] = index
    if sum(r["NOT_RUN"] for r in trait_rows) != 3720:
        raise ValueError("Frozen NOT_RUN family count failed")

    pair_rows = []
    combinations = []
    for pair_id in config["pairs"]:
        left, right = pair_id.split("__")
        pair_counts = Counter()
        details = Counter()
        for locus in loci:
            a, b = by_trait[left][locus], by_trait[right][locus]
            a_state = "TESTED" if a["status"] == "TESTED" else a["reason"]
            b_state = "TESTED" if b["status"] == "TESTED" else b["reason"]
            details[(a_state, b_state)] += 1
            pair_counts["both_numeric"] += a["status"] == "TESTED" and b["status"] == "TESTED"
            pair_counts["at_least_one_not_run"] += a["status"] == "NOT_RUN" or b["status"] == "NOT_RUN"
            pair_counts["both_not_run"] += a["status"] == "NOT_RUN" and b["status"] == "NOT_RUN"
            pair_counts["strict_eligible"] += strict_pass(a, threshold) and strict_pass(b, threshold)
        if pair_counts["both_numeric"] + pair_counts["at_least_one_not_run"] != 2495:
            raise ValueError(f"Pair-locus eligibility accounting failed: {pair_id}")
        pair_rows.append({"pair_id": pair_id, "planned_pair_locus_slots": 2495,
                          "bivariate_TESTED_observed": 0, "bivariate_NOT_RUN_family_gate": 2495,
                          "bivariate_FAILED_observed": 0,
                          "bivariate_status_origin": "DERIVED_NOT_EXECUTED_FAMILY_GATE",
                          "both_univariate_TESTED_numeric": pair_counts["both_numeric"],
                          "fraction_both_numeric": f"{pair_counts['both_numeric']/2495:.8f}",
                          "at_least_one_univariate_NOT_RUN": pair_counts["at_least_one_not_run"],
                          "both_univariate_NOT_RUN": pair_counts["both_not_run"],
                          "both_strict_univariate_gates_pass": pair_counts["strict_eligible"],
                          "fraction_strict_gate_eligible": f"{pair_counts['strict_eligible']/2495:.8f}",
                          "left_trait_not_run": next(r["NOT_RUN"] for r in trait_rows if r["trait_id"] == left),
                          "right_trait_not_run": next(r["NOT_RUN"] for r in trait_rows if r["trait_id"] == right),
                          "pairwise_sample_overlap_status": "LDSC_INTERCEPT_BOUND_FOR_FROZEN_SOURCES_ONLY",
                          "analysis_label": "RESCUE_PLANNING_ONLY"})
        for (a_state, b_state), count in sorted(details.items()):
            combinations.append({"pair_id": pair_id, "left_trait": left, "right_trait": right,
                                 "left_observed_status_or_reason": a_state,
                                 "right_observed_status_or_reason": b_state,
                                 "n_loci": count, "bivariate_status": "NOT_EXECUTED_FAMILY_GATE"})
    if len(pair_rows) != 5 or sum(int(r["planned_pair_locus_slots"]) for r in pair_rows) != 12475:
        raise ValueError("Frozen pair family failed")
    for pair_id in config["pairs"]:
        if sum(r["n_loci"] for r in combinations if r["pair_id"] == pair_id) != 2495:
            raise ValueError(f"Detailed pair reason combination failed: {pair_id}")

    ceiling = decision["maximum_allowed_untested_cells"]
    valid_sets = []
    for n in range(1, len(traits)+1):
        for selected in itertools.combinations(traits, n):
            residual = 3720 - sum(next(r["NOT_RUN"] for r in trait_rows if r["trait_id"] == t) for t in selected)
            if residual <= ceiling:
                valid_sets.append({"traits": list(selected), "perfect_replacement_count": n,
                                   "best_possible_residual_NOT_RUN": residual})
        if valid_sets:
            break
    if len(valid_sets) != 1 or valid_sets[0]["perfect_replacement_count"] != 4 or valid_sets[0]["traits"] != [
        "insomnia", "longsleep", "mdd", "parkinson"
    ] or valid_sets[0]["best_possible_residual_NOT_RUN"] != 533:
        raise ValueError("Frozen theoretical minimal-replacement arithmetic changed")

    OUT.mkdir(parents=True, exist_ok=False)
    files = {
        "trait_failure_matrix.tsv": trait_rows,
        "pair_failure_matrix.tsv": pair_rows,
        "pair_reason_combinations.tsv": combinations,
        "trait_priority_by_not_run.tsv": cause_ranking,
    }
    for name, data in files.items():
        write(OUT / name, data)
    planning = {
        "analysis_label": "RESCUE_PLANNING_ONLY", "source_selection_performed": False,
        "canonical_status": decision["overall_status"], "canonical_NOT_RUN": 3720,
        "frozen_maximum_NOT_RUN": ceiling, "required_reduction": 3720-ceiling,
        "observed_reason_counts": dict(Counter(r["reason"] for r in aggregate if r["status"] == "NOT_RUN")),
        "smallest_perfect_replacement_count": 4,
        "smallest_perfect_replacement_sets": valid_sets,
        "no_single_trait_can_pass": True,
        "binary_N_semantics_status": "SOURCE_SPECIFIC_REVIEW_REQUIRED_BEFORE_RESCUE_INPUT",
        "pairwise_overlap_status": "FROZEN_MATRICES_APPLY_ONLY_TO_ORIGINAL_SOURCES; REESTIMATE_FOR_REPLACEMENTS",
        "unknown_upstream_causes": ["GWAS_coverage", "sample_size_or_power", "harmonization", "allele_ambiguity",
                                    "LD_reference_compatibility", "binary_N_semantics"],
        "interpretation": "Observed LAVA reasons are exact; causal upstream mechanisms behind low local h2 are not proven by the aggregate. Pairwise tests were never executed. Candidate biology was not read or used.",
    }
    (OUT / "planning_summary.json").write_text(json.dumps(planning, indent=2, sort_keys=True) + "\n")
    lines = ["# Confirmatory LAVA rescue planning, before additional source selection", "",
             "**RESCUE_PLANNING_ONLY. Frozen canonical v3 remains FAILED_QC_NOT_PROMOTED.**", "",
             "The immutable seven-trait × 2,495-locus aggregate has 13,745 TESTED, 3,720",
             "NOT_RUN, and zero FAILED trait-locus cells. The frozen limit is 873 NOT_RUN;",
             "a valid rescue must reduce the count by at least 2,847 without altering the",
             "family, thresholds, or criteria. No bivariate v3 tests were executed.", "",
             "## Trait failure matrix", "",
             "| Trait | Planned | TESTED | NOT_RUN | FAILED | Numeric testable | Low local h² | Shared-reference K | Component K |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in sorted(trait_rows, key=lambda x: -x["NOT_RUN"]):
        lines.append(f"| {r['trait_id']} | 2,495 | {r['TESTED']} | {r['NOT_RUN']} | 0 | "
                     f"{float(r['fraction_testable_numeric']):.2%} | {r['low_local_h2_underpowered']} | "
                     f"{r['fewer_than_min_k_shared_reference_variants']} | {r['fewer_than_min_k_components']} |")
    lines.extend(["", "Observed reasons total 3,564 low-local-h², 154 fewer shared-reference variants,",
                  "and two insufficient components. These LAVA outcomes do not identify whether",
                  "underlying GWAS power, variant coverage, harmonization, allele ambiguity,",
                  "LD-reference compatibility, or binary N handling caused a given low-h² cell.",
                  "The sparse-reference group is a direct overlap/representation problem; the",
                  "source-side gap at loci 950–969 is already documented separately. No canonical",
                  "numerical failure was observed. Sample overlap cannot explain this trait-only",
                  "univariate missingness, but any replacement changes the pairwise overlap",
                  "contract and requires new review.", "", "## Pair eligibility and unexecuted tests", "",
                  "| Frozen pair | Planned | Bivariate TESTED | Bivariate NOT_RUN by family gate | Bivariate FAILED observed | Both numeric | At least one NOT_RUN | Both strict gates passed |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|",])
    for r in pair_rows:
        lines.append(f"| {r['pair_id']} | 2,495 | 0 | 2,495 | 0 | {r['both_univariate_TESTED_numeric']} | "
                     f"{r['at_least_one_univariate_NOT_RUN']} | {r['both_strict_univariate_gates_pass']} |")
    lines.extend(["", "The 12,475 pair slots were not executed because the family failed first.",
                  "Zero observed bivariate failures therefore does not establish bivariate",
                  "numerical stability. The detailed pair-reason combination TSV accounts for",
                  "every locus and both observed trait reasons without double-counting pair",
                  "eligibility. Numeric testability and Bonferroni strict-gate eligibility are",
                  "reported separately.", "", "## Source-independent rescue bound", "",
                  "Ranking by NOT_RUN contribution alone gives longsleep, insomnia, parkinson,",
                  "mdd, adhd, bipolar, then scz. Perfectly repairing the three largest leaves",
                  "1,122 NOT_RUN, still above 873. The unique smallest perfect-repair set has",
                  "four traits: insomnia, longsleep, mdd, and parkinson; it would leave 533.",
                  "This is an arithmetic lower bound, not evidence that any such replacement",
                  "exists or will pass. MDD alone leaves 3,131. Engineering the 156 sparse-K",
                  "cells alone leaves 3,564. Do not run a full rescue family until source",
                  "semantics, multi-trait pilots, and a pre-run receipt justify it.", ""])
    (OUT / "rescue_planning_report.md").write_text("\n".join(lines))
    provenance = {"analysis_label": "RESCUE_PLANNING_ONLY", "run_id": RUN_ID,
                  "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                                    (CONFIG, RUN / "canonical_family_decision.json", aggregate_path, CELL_CAUSES, CELL_RECEIPT)},
                  "output_sha256": {name: sha(OUT / name) for name in [*files, "planning_summary.json", "rescue_planning_report.md"]}}
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"traits": len(trait_rows), "pairs": len(pair_rows), "not_run": 3720,
                      "minimal_perfect_repair_traits": valid_sets[0]["traits"]}, sort_keys=True))


if __name__ == "__main__":
    main()
