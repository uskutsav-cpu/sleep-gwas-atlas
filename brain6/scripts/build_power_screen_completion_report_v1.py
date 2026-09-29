#!/usr/bin/env python3
"""Build a source-bound comparison for the completed long-sleep power screen."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/power_optimized_sensitivity_v1"
SCREEN = BASE / "lava_trait_screen_v1"
OUTPUTS = {
    "table": BASE / "completed_trait_only_screen_comparison_v1.tsv",
    "report": BASE / "completed_trait_only_screen_comparison_v1.md",
    "provenance": BASE / "completed_trait_only_screen_comparison_v1.provenance.json",
}
FIELDS = ("metric", "canonical_longsleep", "continuous_duration_candidate", "difference", "interpretation")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"Expected a JSON object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build(root: Path = ROOT) -> tuple[list[dict[str, str]], dict, dict[str, Path]]:
    base = root / "brain6/results/power_optimized_sensitivity_v1"
    screen = base / "lava_trait_screen_v1"
    audit_prov_path = base / "provenance.json"
    aggregate_prov_path = screen / "aggregate.provenance.json"
    aggregate_prov = read_json(aggregate_prov_path)
    summary = read_json(screen / "summary.json")
    audit_prov = read_json(audit_prov_path)
    inputs = {
        "current_trait_audit": base / "current_sleep_trait_audit.tsv",
        "canonical_cause_summary": root / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv",
        "upstream_power_audit_provenance": audit_prov_path,
        "candidate_aggregate": screen / "sleep_duration_continuous_lava_univariate.tsv",
        "candidate_aggregate_summary": screen / "summary.json",
        "candidate_aggregate_provenance": aggregate_prov_path,
        "aggregate_builder": root / "brain6/scripts/aggregate_power_optimized_sleep_screen_v1.py",
    }
    for rel, expected in audit_prov.get("inputs", {}).items():
        path = root / rel
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Power-audit input hash mismatch: {rel}")
    for rel, expected in audit_prov.get("outputs", {}).items():
        path = base / rel
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Power-audit output hash mismatch: {rel}")
    if aggregate_prov.get("status") != "COMPLETE_SCREEN_NO_FAILED_ROWS":
        raise ValueError("Candidate screen does not have a complete aggregate")
    for rel, expected in aggregate_prov.get("inputs", {}).items():
        path = root / rel
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Candidate aggregate input hash mismatch: {rel}")
    for rel, expected in aggregate_prov.get("outputs", {}).items():
        path = root / rel
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Candidate aggregate output hash mismatch: {rel}")
    if aggregate_prov.get("upstream_audit_provenance_sha256") != sha256(audit_prov_path):
        raise ValueError("Candidate aggregate no longer binds to the frozen audit provenance")

    baseline = next(row for row in read_tsv(inputs["current_trait_audit"]) if row["trait"] == "longsleep")
    cause = next(row for row in read_tsv(inputs["canonical_cause_summary"]) if row["trait_id"] == "longsleep")
    candidate_rows = read_tsv(inputs["candidate_aggregate"])
    if len(candidate_rows) != 2495 or len({r["locus_id"] for r in candidate_rows}) != 2495:
        raise ValueError("Candidate aggregate is incomplete or has duplicate loci")
    count = {
        "tested": sum(r["status"] == "TESTED" for r in candidate_rows),
        "not_run": sum(r["status"] == "NOT_RUN" for r in candidate_rows),
        "failed": sum(r["status"] == "FAILED" for r in candidate_rows),
        "strict_gate_pass": sum(r["strict_gate_pass"].lower() == "true" for r in candidate_rows),
        "low_h2": sum(r["reason"] == "LOW_LOCAL_H2_UNDERPOWERED" for r in candidate_rows),
        "min_k": sum(r["reason"] == "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS" for r in candidate_rows),
        "other_min_k": sum(r["reason"] == "FEWER_THAN_MIN_K" for r in candidate_rows),
    }
    if (count["failed"] or count["tested"] + count["not_run"] != 2495 or
            count["tested"] != summary["tested"] or count["not_run"] != summary["not_run"] or
            count["strict_gate_pass"] != summary["strict_gate_pass"]):
        raise ValueError("Candidate rows disagree with their receipt-bound aggregate summary")
    old = {
        "tested": int(baseline["processable_tested"]),
        "not_run": int(baseline["not_run"]),
        "low_h2": int(cause["low_local_h2_underpowered"]),
        "min_k": int(cause["fewer_than_min_k_shared_reference_variants"]),
        "other_min_k": int(cause["fewer_than_min_k"]),
        "strict_gate_pass": int(baseline["strict_univariate_gate_pass"]),
    }
    new = count
    additional = new["strict_gate_pass"] - old["strict_gate_pass"]
    minimum = int(audit_prov["decision_rule"]["minimum_additional_loci"])
    if additional != summary["additional_strict_gate_pass"] or additional >= minimum:
        raise ValueError("Outcome-blinded promotion decision does not match the frozen criterion")
    denominator = 2495
    rows = [
        ("planned_loci", 2495, 2495, 0, "Same frozen locus set."),
        ("processable_tested_loci", old["tested"], new["tested"], new["tested"]-old["tested"], "Finite local-h2 estimates; not evidence of sharing."),
        ("processable_percent", 100*old["tested"]/denominator, 100*new["tested"]/denominator, 100*(new["tested"]-old["tested"])/denominator, "Absolute processability change in percentage points."),
        ("not_run_loci", old["not_run"], new["not_run"], new["not_run"]-old["not_run"], "Total NOT_RUN burden."),
        ("low_local_h2_not_run", old["low_h2"], new["low_h2"], new["low_h2"]-old["low_h2"], "Low-local-h2 failure count."),
        ("low_local_h2_relative_reduction_percent", 0, 100*(old["low_h2"]-new["low_h2"])/old["low_h2"], 100*(old["low_h2"]-new["low_h2"])/old["low_h2"], "Relative decrease from canonical long sleep."),
        ("shared_reference_min_k_failures", old["min_k"], new["min_k"], new["min_k"]-old["min_k"], "Per-locus minimum-K failures; canonical total-panel variant overlap is not reported."),
        ("other_min_k_failures", old["other_min_k"], new["other_min_k"], new["other_min_k"]-old["other_min_k"], "Other frozen minimum-K failure class."),
        ("strict_local_h2_gate_pass", old["strict_gate_pass"], new["strict_gate_pass"], additional, "Frozen p < 0.05/17,465 criterion."),
        ("strict_gate_improvement_percentage_points", 0, 100*additional/denominator, 100*additional/denominator, f"Predeclared minimum: {minimum} additional loci (10 percentage points); not met."),
        ("candidate_reference_overlap", "NOT_REPORTED_COMPARABLY", "6093758/6344850 (96.04%)", "NOT_COMPARABLE", "Candidate versus sealed reference; canonical full-panel overlap unavailable."),
        ("newly_testable_local_rg_cells", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_AVAILABLE", "No bivariate candidate family was launched."),
        ("directional_or_numeric_rg_agreement", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_AVAILABLE", "No candidate local-rg estimates were computed."),
    ]
    table = [dict(zip(FIELDS, map(str, row), strict=True)) for row in rows]
    result = {
        "status": "COMPLETE_SCREEN_NO_PROMOTION",
        "loci": denominator,
        "canonical": old,
        "candidate": new,
        "additional_strict_gate_pass": additional,
        "minimum_additional_for_full_sensitivity": minimum,
        "full_sensitivity_eligible": False,
        "canonical_decision_unchanged": "FAILED_QC_NOT_PROMOTED",
        "downstream_association_results_consulted": False,
    }
    return table, result, inputs


def main() -> None:
    output_paths = {
        "table": ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.tsv",
        "report": ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.md",
        "provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.provenance.json",
    }
    if any(path.exists() for path in output_paths.values()):
        raise FileExistsError("Refusing to overwrite a completion report artifact")
    table_rows, summary, inputs = build()
    from io import StringIO
    stream = StringIO()
    writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(table_rows)
    table_text = stream.getvalue()
    old, new = summary["canonical"], summary["candidate"]
    report_text = f"""# Continuous sleep-duration power screen: completion and decision

Snapshot: 2026-09-26. Brain6-only; the canonical long-sleep result remains immutable and FAILED_QC_NOT_PROMOTED.

The four checksum-bound, disjoint shards completed all 2,495 loci with no failed rows. The continuous Dashti 2019 candidate yielded {new['tested']:,} finite local-h2 estimates ({100*new['tested']/2495:.2f}%) and {new['not_run']:,} NOT_RUN, all from low local h2. Canonical long sleep yielded {old['tested']:,} processable loci ({100*old['tested']/2495:.2f}%) and {old['not_run']:,} NOT_RUN ({old['low_h2']:,} low-local-h2; {old['min_k']:,} shared-reference minimum-K).

Processability improved by {new['tested']-old['tested']:,} loci ({100*(new['tested']-old['tested'])/2495:.2f} percentage points). Low-local-h2 NOT_RUN fell by {old['low_h2']-new['low_h2']:,} ({100*(old['low_h2']-new['low_h2'])/old['low_h2']:.2f}% relative). Strict frozen local-h2 gate passes increased from {old['strict_gate_pass']} to {new['strict_gate_pass']}: {summary['additional_strict_gate_pass']} additional loci ({100*summary['additional_strict_gate_pass']/2495:.2f} percentage points).

The predeclared threshold is {summary['minimum_additional_for_full_sensitivity']} additional strict-gate loci (10 percentage points). This candidate adds only {summary['additional_strict_gate_pass']}, so it is not eligible for a full bivariate sensitivity family and does not replace canonical long sleep. No new local-rg cells or estimates exist; direction/numeric agreement is not applicable.

The candidate is a related but distinct continuous self-reported duration estimand from the same study/UK Biobank source; substantial participant overlap is expected. It has 6,093,758/6,344,850 QC-retained variants in the frozen reference (96.04%). Canonical long-sleep total variant overlap with the full panel is not reported; per-locus minimum-K failures are {old['min_k']} canonical versus {new['min_k']} candidate.

No canonical thresholds, receipts, run directories, or downstream analyses were changed. Detailed machine-readable comparison: completed_trait_only_screen_comparison_v1.tsv.
"""
    output_paths["table"].parent.mkdir(parents=True, exist_ok=True)
    contents = {"table": table_text, "report": report_text}
    for name, content in contents.items():
        with output_paths[name].open("x", encoding="utf-8", newline="") as stream:
            stream.write(content)
    provenance = {
        "schema_version": 1,
        "status": summary["status"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "decision_rule": {"minimum_additional_strict_gate_loci": summary["minimum_additional_for_full_sensitivity"],
                          "downstream_association_results_consulted": False},
        "inputs": {name: sha256(path) for name, path in inputs.items()},
        "outputs": {str(output_paths[name].relative_to(ROOT)): hashlib.sha256(content.encode()).hexdigest()
                    for name, content in contents.items()},
    }
    with output_paths["provenance"].open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
