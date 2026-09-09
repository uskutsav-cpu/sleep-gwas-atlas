#!/usr/bin/env python3
"""Build the specification-complete Pair A/B independent-replication summary."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAIR_MANIFEST = ROOT / "results/track_b/pair_manifest.tsv"
SOURCE_AUDIT = ROOT / "results/track_b/replication_source_audit.tsv"
SOURCE_LOCK = ROOT / "results/track_b/replication_source.lock.json"
PAIR_B_RESULT = ROOT / "results/track_b/replication/pair_b_ldsc.tsv"
PAIR_B_PROVENANCE = ROOT / "results/track_b/replication/pair_b_ldsc.provenance.json"
OUTPUT = ROOT / "results/track_b/02_independent_global_replication.tsv"
PROVENANCE = ROOT / "results/track_b/02_independent_global_replication.provenance.json"
Z_975 = 1.959963984540054

FIELDS = [
    "pair_id", "sleep_trait", "external_trait", "discovery_rg", "discovery_SE", "discovery_P",
    "discovery_CI_lower", "discovery_CI_upper", "replication_source_id", "replication_rg",
    "replication_SE", "replication_P", "replication_CI_lower", "replication_CI_upper",
    "direction", "heterogeneity_z", "heterogeneity_P", "heterogeneity_assumption",
    "sample_overlap", "phenotype_equivalence", "replication_ancestry", "replication_build",
    "replication_sample_size", "replication_case_control", "independence_status", "power_status",
    "h2_gate_status", "h2", "h2_SE", "h2_z", "cross_trait_intercept",
    "cross_trait_intercept_SE", "replication_status", "go_no_go", "classification_rule",
    "claim_limit",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty replication artifact: {path.relative_to(ROOT)}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite(value: str, label: str) -> float:
    try:
        observed = float(value)
    except ValueError:
        fail(f"invalid {label}: {value}")
    if not math.isfinite(observed):
        fail(f"non-finite {label}: {value}")
    return observed


def fmt(value: float) -> str:
    return f"{value:.12g}"


def ci(estimate: float, standard_error: float) -> tuple[str, str]:
    return fmt(estimate - Z_975 * standard_error), fmt(estimate + Z_975 * standard_error)


def tsv_text(rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def build() -> tuple[str, dict[str, object]]:
    pairs = {row["pair_id"]: row for row in read_tsv(PAIR_MANIFEST)}
    if set(pairs) != {"A", "B", "CONTROL"}:
        fail("frozen Track B pair manifest drifted")
    audit = read_tsv(SOURCE_AUDIT)
    if any(row["results_accessed_before_selection"] != "FALSE" for row in audit):
        fail("replication source selection did not precede result access")
    pair_a_terminal = [row for row in audit if row["pair_id"] == "A" and row["selection_decision"] == "NO_VALID_REPLICATION"]
    pair_b_sources = [row for row in audit if row["pair_id"] == "B" and row["selection_decision"] == "SELECT_PRIMARY_REPLICATION"]
    if len(pair_a_terminal) != 1 or len(pair_b_sources) != 1 or pair_b_sources[0]["source_id"] != "finngen_r13_F5_ADHD":
        fail("frozen Pair A/B replication decisions drifted")
    pair_b_result = read_tsv(PAIR_B_RESULT)
    if len(pair_b_result) != 1 or pair_b_result[0]["pair_id"] != "B":
        fail("Pair B LDSC result is not the exact frozen pair")
    raw = pair_b_result[0]
    raw_provenance = json.loads(PAIR_B_PROVENANCE.read_text(encoding="utf-8"))
    if raw_provenance.get("analysis_status") != "PAIR_B_REPLICATION_LDSC_COMPLETE" or raw_provenance.get("output_sha256") != sha256(PAIR_B_RESULT):
        fail("Pair B LDSC provenance is incomplete or drifted")
    source = pair_b_sources[0]

    rows: list[dict[str, str]] = []
    pair_a = pairs["A"]
    discovery_a, discovery_a_se = finite(pair_a["discovery_rg"], "Pair A discovery rg"), finite(pair_a["SE"], "Pair A discovery SE")
    discovery_a_ci = ci(discovery_a, discovery_a_se)
    rows.append({
        "pair_id": "A", "sleep_trait": pair_a["sleep_trait"], "external_trait": pair_a["external_trait"],
        "discovery_rg": pair_a["discovery_rg"], "discovery_SE": pair_a["SE"], "discovery_P": pair_a["P"],
        "discovery_CI_lower": discovery_a_ci[0], "discovery_CI_upper": discovery_a_ci[1],
        "replication_source_id": pair_a_terminal[0]["source_id"], "replication_rg": "NA",
        "replication_SE": "NA", "replication_P": "NA", "replication_CI_lower": "NA",
        "replication_CI_upper": "NA", "direction": "NOT_APPLICABLE", "heterogeneity_z": "NA",
        "heterogeneity_P": "NA", "heterogeneity_assumption": "NOT_COMPUTABLE_NO_VALID_REPLICATION",
        "sample_overlap": pair_a_terminal[0]["cohort_overlap_assessment"],
        "phenotype_equivalence": pair_a_terminal[0]["phenotype_relation"],
        "replication_ancestry": pair_a_terminal[0]["ancestry"], "replication_build": pair_a_terminal[0]["build"],
        "replication_sample_size": pair_a_terminal[0]["sample_size"], "replication_case_control": "NA/NA",
        "independence_status": pair_a_terminal[0]["independence_status"], "power_status": pair_a_terminal[0]["power_status"],
        "h2_gate_status": "NOT_RUN", "h2": "NA", "h2_SE": "NA", "h2_z": "NA",
        "cross_trait_intercept": "NA", "cross_trait_intercept_SE": "NA",
        "replication_status": "NO_VALID_REPLICATION", "go_no_go": "PROCEED_AS_UNREPLICATED_TIER1_DISCOVERY",
        "classification_rule": "Track B Section 5 terminal source audit: no exact independent public analysis-ready GWAS",
        "claim_limit": "Pair A remains unreplicated and cannot exceed Tier 1 without a valid independent dataset",
    })

    pair_b = pairs["B"]
    discovery_b, discovery_b_se = finite(pair_b["discovery_rg"], "Pair B discovery rg"), finite(pair_b["SE"], "Pair B discovery SE")
    replication_b, replication_b_se = finite(raw["rg"], "Pair B replication rg"), finite(raw["rg_se"], "Pair B replication SE")
    discovery_b_ci, replication_b_ci = ci(discovery_b, discovery_b_se), ci(replication_b, replication_b_se)
    heterogeneity_z = (discovery_b - replication_b) / math.sqrt(discovery_b_se**2 + replication_b_se**2)
    heterogeneity_p = math.erfc(abs(heterogeneity_z) / math.sqrt(2))
    concordant = raw["direction_vs_discovery"] == "CONCORDANT"
    significant = finite(raw["rg_p"], "Pair B replication P") < 0.05
    if raw["h2_gate_status"] != "PASS":
        status, gate = "UNDERPOWERED", "PROCEED_WITHOUT_REPLICATION_CREDIT"
    elif significant and not concordant:
        status, gate = "OPPOSITE_DIRECTION", "NO_GO_STOP_MECHANISTIC_ESCALATION"
    elif significant and concordant:
        # Individual-level overlap is not verifiable and the Finnish register phenotype/founder LD differ.
        status, gate = "DIRECTIONAL_REPLICATION", "GO_WITH_FINNISH_FOUNDER_AND_PHENOTYPE_CAVEATS"
    else:
        status, gate = "UNDERPOWERED", "PROCEED_WITHOUT_REPLICATION_CREDIT"
    if raw["replication_class"] != status:
        fail("Pair B raw LDSC classification differs from the specification-complete summary rule")
    rows.append({
        "pair_id": "B", "sleep_trait": pair_b["sleep_trait"], "external_trait": pair_b["external_trait"],
        "discovery_rg": pair_b["discovery_rg"], "discovery_SE": pair_b["SE"], "discovery_P": pair_b["P"],
        "discovery_CI_lower": discovery_b_ci[0], "discovery_CI_upper": discovery_b_ci[1],
        "replication_source_id": source["source_id"], "replication_rg": raw["rg"], "replication_SE": raw["rg_se"],
        "replication_P": raw["rg_p"], "replication_CI_lower": replication_b_ci[0], "replication_CI_upper": replication_b_ci[1],
        "direction": raw["direction_vs_discovery"], "heterogeneity_z": fmt(heterogeneity_z),
        "heterogeneity_P": fmt(heterogeneity_p),
        "heterogeneity_assumption": "INDEPENDENT_ESTIMATE_APPROXIMATION_WITH_ZERO_COVARIANCE;INDIVIDUAL_OVERLAP_NOT_VERIFIED",
        "sample_overlap": source["cohort_overlap_assessment"], "phenotype_equivalence": source["phenotype_relation"],
        "replication_ancestry": source["ancestry"], "replication_build": source["build"],
        "replication_sample_size": source["sample_size"],
        "replication_case_control": f"{source['cases']}/{source['controls']}",
        "independence_status": source["independence_status"], "power_status": source["power_status"],
        "h2_gate_status": raw["h2_gate_status"], "h2": raw["h2"], "h2_SE": raw["h2_se"], "h2_z": raw["h2_z"],
        "cross_trait_intercept": raw["cross_trait_intercept"],
        "cross_trait_intercept_SE": raw["cross_trait_intercept_se"], "replication_status": status,
        "go_no_go": gate,
        "classification_rule": "Track B Section 5: significant concordant external-cohort rg earns directional credit; STRONG requires verified nonoverlap and closer population/phenotype equivalence",
        "claim_limit": "Replicates the global direction in a Finnish register cohort; does not establish a local shared locus, mechanism, mediation, or causality",
    })

    table = tsv_text(rows)
    inputs = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (PAIR_MANIFEST, SOURCE_AUDIT, SOURCE_LOCK, PAIR_B_RESULT, PAIR_B_PROVENANCE)
    }
    provenance = {
        "schema_version": 1,
        "analysis_status": "TRACK_B_INDEPENDENT_GLOBAL_REPLICATION_COMPLETE",
        "pair_order": ["A", "B"],
        "classification_family": [
            "STRONG_REPLICATION", "DIRECTIONAL_REPLICATION", "UNDERPOWERED",
            "FAILED_REPLICATION", "OPPOSITE_DIRECTION", "NO_VALID_REPLICATION",
        ],
        "pair_a_terminal_status": "NO_VALID_REPLICATION",
        "pair_b_terminal_status": rows[1]["replication_status"],
        "input_sha256": inputs,
        "output_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "script_sha256": sha256(Path(__file__)),
    }
    return table, provenance


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    table, provenance = build()
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    if args.verify:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != table:
            fail("canonical Track B independent-replication table is absent or drifted")
        if not PROVENANCE.is_file() or PROVENANCE.read_text(encoding="utf-8") != provenance_text:
            fail("canonical Track B independent-replication provenance is absent or drifted")
        print("verified Track B independent global replication: A=NO_VALID_REPLICATION B=DIRECTIONAL_REPLICATION")
        return
    atomic_text(OUTPUT, table)
    atomic_text(PROVENANCE, provenance_text)
    print("wrote Track B independent global replication: A=NO_VALID_REPLICATION B=DIRECTIONAL_REPLICATION")


if __name__ == "__main__":
    main()
