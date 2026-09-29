#!/usr/bin/env python3
"""Freeze the canonical Brain6 repair burden before new source outcomes."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/confirmatory_source_rescue_20260927/phase1"
INPUT = {
    "canonical": ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv",
    "source_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/canonical_family_trait_power_audit_v1.tsv",
    "pair": ROOT / "brain6/results/lava_confirmatory_rescue_plan_v1/pair_failure_matrix.tsv",
    "pair_reasons": ROOT / "brain6/results/lava_confirmatory_rescue_plan_v1/pair_reason_combinations.tsv",
    "source_triage": ROOT / "brain6/results/lava_longsleep_source_rescue_v1/OTHER_BINARY_SOURCE_TRIAGE.md",
    "long_sleep_semantics": ROOT / "brain6/results/lava_longsleep_source_rescue_v1/LAVA_SAMPLE_SIZE_SEMANTICS.md",
    "historical_harmonizer": ROOT / "scripts/01_harmonize.py",
    "insomnia_historical_qc": Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/insomnia.qc.txt"),
    "decision": ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json",
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(name: str) -> list[dict[str, str]]:
    with INPUT[name].open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(path)
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


SOURCE_SEMANTICS = {
    "longsleep": ("BINARY", "NO_NATIVE_VARIANT_N", "STUDY_LEVEL_34184_CASES_305742_CONTROLS;HISTORICAL_CONSTANT_NEFF", "EXTERNAL_SOURCE_N_AND_MODEL_REQUIRED"),
    "insomnia": ("BINARY", "NATIVE_PER_VARIANT_N_PRESENT_BUT_DISCARDED_CANONICAL", "STUDY_LEVEL_109402_CASES_277131_CONTROLS;HISTORICAL_CONSTANT_NEFF", "SOURCE_MODEL_AND_VARIANT_COVERAGE_REVIEW_POSSIBLE"),
    "parkinson": ("BINARY_PROXY_MIXED", "NATIVE_PER_VARIANT_CASE_CONTROL_FIELDS", "VARYING_CASE_CONTROL_COUNTS;HISTORICAL_NEFF_DERIVED", "SOURCE_FAITHFUL_REPRESENTATION_AND_COVERAGE_REVIEW_POSSIBLE"),
    "mdd": ("BINARY", "NO_NATIVE_VARIANT_N_IN_PRESERVED_RELEASE", "STUDY_LEVEL_CASE_CONTROL_ONLY;HISTORICAL_CONSTANT_N", "EXTERNAL_SOURCE_N_OR_EXACT_REPLACEMENT_REQUIRED"),
    "adhd": ("BINARY", "NATIVE_PER_VARIANT_CASE_CONTROL_FIELDS", "VARYING_CASE_CONTROL_COUNTS;HISTORICAL_NEFF_DERIVED", "SOURCE_FAITHFUL_REPRESENTATION_REVIEW_POSSIBLE"),
    "bipolar": ("BINARY", "PER_VARIANT_EFFECTIVE_COUNT_FIELDS_ONLY", "NCAS_NCON_HEADER_DEFINES_EFFECTIVE_COUNTS", "LITERAL_N_OR_METHOD_MAPPING_REQUIRED"),
    "scz": ("BINARY", "PER_VARIANT_EFFECTIVE_COUNT_FIELDS_ONLY", "NCAS_NCON_HEADER_DEFINES_EFFECTIVE_COUNTS", "LITERAL_N_OR_METHOD_MAPPING_REQUIRED"),
}


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Frozen output already exists: {OUT}")
    frozen = json.loads(INPUT["decision"].read_text())
    if frozen.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or frozen.get("maximum_allowed_untested_cells") != 873:
        raise ValueError("Canonical decision changed")
    source = {r["trait_id"]: r for r in read("source_audit")}
    canonical = {r["trait_id"]: r for r in read("canonical")}
    if len(source) != len(canonical) or set(source) != set(SOURCE_SEMANTICS) or len(source) != 7:
        raise ValueError("Seven-trait source/canonical join changed")
    trait_rows = []
    for trait, c in canonical.items():
        s = source[trait]
        kind, per_snp_n, counts, feasibility = SOURCE_SEMANTICS[trait]
        if any(int(c[x]) != int(s[y]) for x, y in (
            ("tested_cells", "tested"), ("not_run_cells", "NOT_RUN"), ("failed_cells", "failed"),
            ("low_local_h2_underpowered", "low_local_h2_NOT_RUN"),
            ("fewer_than_min_k_shared_reference_variants", "shared_reference_minK_NOT_RUN"),
            ("fewer_than_min_k", "other_minK_NOT_RUN"))):
            raise ValueError(f"Frozen source/canonical mismatch: {trait}")
        trait_rows.append({
            "rank_by_not_run": 0, "trait_id": trait, "implicated_source": s["study"],
            "PMID": s["PMID"], "DOI": s["DOI"],
            "planned_cells": c["planned_cells"], "TESTED": c["tested_cells"],
            "NOT_RUN": c["not_run_cells"], "FAILED": c["failed_cells"],
            "testability_pct": s["tested_pct"],
            "low_local_h2_not_run": c["low_local_h2_underpowered"],
            "shared_reference_min_k_not_run": c["fewer_than_min_k_shared_reference_variants"],
            "component_min_k_not_run": c["fewer_than_min_k"],
            "source_sample_size": s["source_sample_size"], "cases": s["cases"], "controls": s["controls"],
            "available_dense_unique_snps": s["dense_unique_variants_after_frozen_QC"],
            "exact_snp_id_reference_overlap": s["exact_SNP_ID_reference_overlap"],
            "overlap_scope": "EXACT_SNP_ID_AND_CHROMOSOME;NOT_ALLELE_HARMONIZATION",
            "ancestry": s["ancestry"], "genome_build": s["genome_build"],
            "phenotype_type": kind, "variant_n_status": per_snp_n,
            "case_control_n_semantics": counts,
            "maximum_theoretical_recovery_cells": c["not_run_cells"],
            "source_first_repair_feasibility": feasibility,
            "classification": "FROZEN_PRE_OUTCOME_RESCUE_PLANNING",
        })
    trait_rows.sort(key=lambda r: (-int(r["NOT_RUN"]), r["trait_id"]))
    for i, row in enumerate(trait_rows, 1):
        row["rank_by_not_run"] = i
    if sum(int(r["NOT_RUN"]) for r in trait_rows) != 3720:
        raise ValueError("Canonical NOT_RUN total changed")
    pair_rows = []
    reasons = read("pair_reasons")
    for p in read("pair"):
        pair = p["pair_id"]
        pair_reasons = [r for r in reasons if r["pair_id"] == pair]
        if sum(int(r["n_loci"]) for r in pair_reasons) != 2495:
            raise ValueError(f"Pair reason coverage changed: {pair}")
        left, right = pair.split("__")
        pair_rows.append({
            "pair_id": pair, "planned_pair_locus_slots": p["planned_pair_locus_slots"],
            "bivariate_TESTED_observed": p["bivariate_TESTED_observed"],
            "bivariate_NOT_RUN_family_gate": p["bivariate_NOT_RUN_family_gate"],
            "bivariate_FAILED_observed": p["bivariate_FAILED_observed"],
            "both_univariate_TESTED_numeric": p["both_univariate_TESTED_numeric"],
            "both_numeric_testability_pct": f"{100*float(p['fraction_both_numeric']):.4f}",
            "at_least_one_univariate_NOT_RUN": p["at_least_one_univariate_NOT_RUN"],
            "both_univariate_NOT_RUN": p["both_univariate_NOT_RUN"],
            "both_strict_univariate_gates_pass": p["both_strict_univariate_gates_pass"],
            "left_trait_not_run": p["left_trait_not_run"],
            "right_trait_not_run": p["right_trait_not_run"],
            "left_source": source[left]["study"], "right_source": source[right]["study"],
            "exact_reason_combinations": len(pair_reasons),
            "reason_detail_source": str(INPUT["pair_reasons"].relative_to(ROOT)),
            "pair_overlap_status": p["pairwise_sample_overlap_status"],
            "interpretation": "PAIRWISE_NOT_EXECUTED_FROZEN_FAMILY_GATE",
        })
    not_run = [int(r["NOT_RUN"]) for r in trait_rows]
    minimal = []
    for n in range(1, 8):
        options = [tuple(trait_rows[i]["trait_id"] for i in indices)
                   for indices in itertools.combinations(range(7), n)
                   if 3720 - sum(not_run[i] for i in indices) <= 873]
        if options:
            minimal = options
            break
    if n != 4 or minimal != [("longsleep", "insomnia", "parkinson", "mdd")]:
        raise ValueError("Frozen minimum theoretical repair set changed")
    OUT.mkdir(parents=True, exist_ok=False)
    write(OUT / "trait_failure_burden.tsv", trait_rows)
    write(OUT / "pair_failure_burden.tsv", pair_rows)
    summary = {
        "analysis_label": "FROZEN_PRE_OUTCOME_RESCUE_PLANNING",
        "canonical_status": frozen["overall_status"],
        "canonical_tested": 13745, "canonical_not_run": 3720, "canonical_failed": 0,
        "frozen_maximum_not_run": 873, "minimum_required_recovery": 2847,
        "smallest_perfect_repair_set": list(minimal[0]),
        "best_possible_residual_for_that_set": 533,
        "long_sleep_perfect_repair_alone_residual": 2429,
        "exact_reason_counts": {"LOW_LOCAL_H2_UNDERPOWERED": 3564,
                                "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS": 154,
                                "FEWER_THAN_MIN_K": 2},
        "caveat": "Observed LAVA reasons are exact. Upstream causes and actual recoverability are unproven; pairwise tests were not executed.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    provenance = {
        "schema_version": 1, "source_sha256": {
            str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path): sha(path)
            for path in INPUT.values()},
        "builder_sha256": sha(Path(__file__)),
        "output_sha256": {p.name: sha(p) for p in OUT.iterdir() if p.is_file()},
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
