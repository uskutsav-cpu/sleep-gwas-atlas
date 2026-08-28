#!/usr/bin/env python3
"""Build transparent Phase-1 novelty, local-analysis, and final ranking tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path


DIRECT = {
    "DIRECT_RG_PREVIOUSLY_REPORTED",
    "DIRECT_RG_REPLICATION_DIFFERENT_DATASET",
}
NO_DIRECT = {
    "APPARENTLY_NOVEL", "RELATED_GENETIC_EVIDENCE_ONLY", "MR_ONLY",
    "OBSERVATIONAL_ONLY", "NO_DIRECT_RG_FOUND", "UNCERTAIN",
}
TIER_ORDER = {"TIER_A": 0, "TIER_B": 1, "TIER_C": 2, "TIER_D": 3}
CONFIDENCE_ORDER = {"HIGH_CONFIDENCE": 0, "MODERATE_CONFIDENCE": 1, "QC_CAUTION": 2}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "NA") for field in fields})
    return buffer.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def pair_key(row: dict[str, str]) -> tuple[str, str]:
    return row["sleep_trait"], row["external_trait"]


def pair_id(row: dict[str, str]) -> str:
    return "__".join(pair_key(row))


def h2_strength(row: dict[str, str]) -> str:
    sleep_z = abs(float(row["sleep_h2"]) / float(row["sleep_h2_se"]))
    external_z = abs(float(row["external_h2"]) / float(row["external_h2_se"]))
    if min(sleep_z, external_z) >= 10:
        return "STRONG_BOTH_H2_Z_GE_10"
    if min(sleep_z, external_z) >= 4:
        return "ADEQUATE_BOTH_H2_Z_GE_4"
    return "WEAK_AT_LEAST_ONE_H2_Z_LT_4"


def biological_interpretability(row: dict[str, str]) -> str:
    if row["external_trait"] in {"parental_lifespan", "healthspan", "frailty", "telomere_length"}:
        return "AGING_AXIS_WITH_CLEAR_DIRECTIONAL_HEALTH_INTERPRETATION"
    if row["external_domain"] in {"metabolic", "cardio"}:
        return "CARDIOMETABOLIC_AXIS_WITH_TESTABLE_SHARED_PATHWAYS"
    if row["external_domain"] == "psychiatric":
        return "NEUROBEHAVIORAL_AXIS_WITH_TESTABLE_SHARED_PATHWAYS"
    if row["external_domain"] == "immune":
        return "IMMUNE_AXIS_REQUIRING_DIRECTION_AND_PHENOTYPE_REVIEW"
    return "BIOLOGICALLY_PLAUSIBLE_BUT_MECHANISM_UNRESOLVED"


def novelty_tier(row: dict[str, str], audit: dict[str, str]) -> tuple[str, str]:
    classification = audit["novelty_classification"]
    clean = row["confidence_class"] != "QC_CAUTION"
    dense = "DENSE" in audit["resource_context"]
    if classification == "APPARENTLY_NOVEL" and float(row["abs_rg"]) >= .20 and clean and dense:
        return (
            "TIER_A",
            "Apparently unreported; |rg|>=0.20; no QC_CAUTION; dense summary statistics support an alternative-dataset replication design, although exact independent datasets may be limited.",
        )
    if classification in NO_DIRECT and clean and float(row["abs_rg"]) >= .10:
        return (
            "TIER_B",
            "No explicit direct rg located; clean result with |rg|>=0.10, but exact replication or literature status remains incomplete.",
        )
    if classification in NO_DIRECT:
        return (
            "TIER_C",
            "No explicit direct rg located, but effect size, QC, phenotype definition, or resource availability reduces immediate priority.",
        )
    return (
        "TIER_D",
        "Direct prior rg evidence makes this a known-replication or pipeline-positive-control finding.",
    )


def select_novelty_rows(
    master: list[dict[str, str]], audits: dict[tuple[str, str], dict[str, str]],
) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    for row in master:
        audit = audits[pair_key(row)]
        tier, tier_reason = novelty_tier(row, audit)
        candidates.append({
            "pair_id": pair_id(row),
            "sleep_trait": row["sleep_trait"],
            "external_trait": row["external_trait"],
            "external_domain": row["external_domain"],
            "rg": row["rg"], "se": row["se"], "p": row["p"], "fdr": row["fdr"],
            "abs_rg": row["abs_rg"],
            "sleep_h2": row["sleep_h2"], "sleep_h2_se": row["sleep_h2_se"],
            "external_h2": row["external_h2"], "external_h2_se": row["external_h2_se"],
            "h2_strength": h2_strength(row),
            "qc_confidence": row["confidence_class"],
            "literature_novelty": audit["novelty_classification"],
            "phenotype_distinctiveness": audit["phenotype_definition_match"],
            "independent_gwas_availability": audit["replication_feasibility"],
            "dense_summary_statistics_availability": audit["resource_context"],
            "molecular_qtl_availability": (
                "LIKELY_GENERAL_QTL_RESOURCES_BUT_PAIR_SPECIFIC_COLOCALIZATION_NOT_AUDITED"
                if row["external_domain"] != "aging"
                else "AGING_QTL_RESOURCES_EXIST_BUT_ENDPOINT_MATCHING_REQUIRES_REVIEW"
            ),
            "biological_interpretability": biological_interpretability(row),
            "priority_tier": tier,
            "tier_rule_applied": tier_reason,
            "ranking_rule": "Within tier: descending |rg|, then ascending locked FDR, then pair_id; no composite score.",
        })
    selected: list[dict[str, object]] = []
    quotas = {"TIER_A": 5, "TIER_B": 10, "TIER_C": 5, "TIER_D": 5}
    for tier in TIER_ORDER:
        tier_rows = [row for row in candidates if row["priority_tier"] == tier]
        tier_rows.sort(key=lambda row: (-float(row["abs_rg"]), float(row["fdr"]), row["pair_id"]))
        selected.extend(tier_rows[:quotas[tier]])
    selected.sort(key=lambda row: (TIER_ORDER[str(row["priority_tier"])], -float(row["abs_rg"]), float(row["fdr"]), row["pair_id"]))
    for rank, row in enumerate(selected, 1):
        row["priority_rank"] = rank
    fields = ["priority_rank"] + [field for field in selected[0] if field != "priority_rank"]
    return [{field: row[field] for field in fields} for row in selected]


def build_local_priorities(
    primary_all: list[dict[str, str]],
    significant: list[dict[str, str]],
    audits: dict[tuple[str, str], dict[str, str]],
    cancellation: list[dict[str, str]],
    measurement: list[dict[str, str]],
    bridge_sleep: set[str],
    external_degrees: dict[str, int],
) -> list[dict[str, object]]:
    master = {pair_key(row): row for row in primary_all}
    sig = {pair_key(row): row for row in significant}
    reasons: defaultdict[tuple[str, str], list[tuple[str, str, str]]] = defaultdict(list)

    for key, audit in audits.items():
        if audit["novelty_classification"] == "APPARENTLY_NOVEL":
            reasons[key].append(("A", "LAVA;HDL-L;PLACO;COLOC", "Strongest apparently unreported global connection; test whether sharing is localized and replicable."))

    known = [row for row in significant if audits[pair_key(row)]["novelty_classification"] in DIRECT]
    known.sort(key=lambda row: (-float(row["abs_rg"]), float(row["fdr"]), pair_id(row)))
    for row in known[:6]:
        reasons[pair_key(row)].append(("B", "LAVA;HDL-L", "Strong known direct-rg positive control for local-analysis calibration."))

    negatives = [row for row in significant if float(row["rg"]) < 0]
    negatives.sort(key=lambda row: (-float(row["abs_rg"]), float(row["fdr"]), pair_id(row)))
    for row in negatives[:12]:
        reasons[pair_key(row)].append(("C", "LAVA;PLACO;COLOC", "Unusual negative global rg; test signed local sharing without calling antagonistic pleiotropy in advance."))

    discordant_external = {
        row["external_trait"] for row in measurement
        if row["significant_pair_direction_discordance"] == "True"
    }
    discordant = [row for row in significant if row["external_trait"] in discordant_external]
    discordant.sort(key=lambda row: (-float(row["abs_rg"]), float(row["fdr"]), pair_id(row)))
    for row in discordant[:12]:
        reasons[pair_key(row)].append(("D", "LAVA;GENOMIC_SEM", "External trait shows significant objective-versus-subjective direction discordance; test measurement-specific local architecture."))

    for row in cancellation:
        key = pair_key(row)
        reasons[key].append(("E", "LAVA;HDL-L;MiXeR", "Globally weak/null candidate with strong related-trait effects; test for cancelling signed local components. Cancellation is not yet demonstrated."))

    bridges = [
        row for row in significant
        if row["sleep_trait"] in bridge_sleep and external_degrees.get(row["external_trait"], 0) >= 3
    ]
    bridges.sort(key=lambda row: (-float(row["abs_rg"]), float(row["fdr"]), pair_id(row)))
    for row in bridges[:12]:
        reasons[pair_key(row)].append(("F", "PLACO;LAVA;GENOMIC_SEM", "High-weight edge linking a broad sleep-domain bridge to an external multi-sleep hub."))

    rows: list[dict[str, object]] = []
    for key, items in reasons.items():
        atlas = master[key]
        categories = sorted({item[0] for item in items})
        methods = sorted({method for item in items for method in item[1].split(";")})
        reason_text = " | ".join(dict.fromkeys(item[2] for item in items))
        audit = audits.get(key)
        rows.append({
            "priority_rank": 0,
            "sleep_trait": key[0], "external_trait": key[1],
            "external_domain": atlas["external_domain"],
            "rg": atlas["rg"], "se": atlas["se"], "p": atlas["p"], "fdr": atlas["fdr"],
            "abs_rg": atlas["abs_rg"],
            "locked_fdr_significant": str(key in sig),
            "priority_categories": ";".join(categories),
            "recommended_methods": ";".join(methods),
            "confidence_class": atlas["confidence_class"],
            "literature_status": audit["novelty_classification"] if audit else "NOT_APPLICABLE_GLOBAL_NULL",
            "main_reason_for_priority": reason_text,
            "claim_limit": "Global rg is not a local mechanism; all local sharing, cancellation, and antagonism labels remain hypotheses until tested.",
            "ranking_rule": "Lexicographic by earliest category A-F, then descending |rg|, then ascending locked FDR; no composite score.",
        })
    rows.sort(key=lambda row: (str(row["priority_categories"])[0], -float(row["abs_rg"]), float(row["fdr"]), f"{row['sleep_trait']}__{row['external_trait']}"))
    for rank, row in enumerate(rows, 1):
        row["priority_rank"] = rank
    return rows


def build_top_findings(
    significant: list[dict[str, str]],
    audits: dict[tuple[str, str], dict[str, str]],
    local: list[dict[str, object]],
    bridge_sleep: set[str],
    external_degrees: dict[str, int],
) -> list[dict[str, object]]:
    local_by_key = {(str(row["sleep_trait"]), str(row["external_trait"])): row for row in local}
    novelty_order = {
        "APPARENTLY_NOVEL": 0, "NO_DIRECT_RG_FOUND": 1, "UNCERTAIN": 2,
        "RELATED_GENETIC_EVIDENCE_ONLY": 3, "MR_ONLY": 4, "OBSERVATIONAL_ONLY": 5,
        "DIRECT_RG_REPLICATION_DIFFERENT_DATASET": 6, "DIRECT_RG_PREVIOUSLY_REPORTED": 7,
    }
    views = {
        "EFFECT_MAGNITUDE": (
            lambda row: (-float(row["abs_rg"]), float(row["fdr"]), pair_id(row)),
            "Descending |rg|, ascending locked FDR, pair_id; no composite score.",
        ),
        "STATISTICAL_EVIDENCE": (
            lambda row: (float(row["fdr"]), -float(row["abs_rg"]), pair_id(row)),
            "Ascending locked-family BH FDR, descending |rg|, pair_id; no composite score.",
        ),
        "NOVELTY": (
            lambda row: (novelty_order[audits[pair_key(row)]["novelty_classification"]], -float(row["abs_rg"]), float(row["fdr"]), pair_id(row)),
            "Literature class (apparently novel first, direct prior rg last), descending |rg|, ascending FDR, pair_id; no composite score.",
        ),
        "CONFIDENCE": (
            lambda row: (CONFIDENCE_ORDER[row["confidence_class"]], 0 if h2_strength(row).startswith("STRONG") else 1, float(row["fdr"]), -float(row["abs_rg"]), pair_id(row)),
            "QC confidence class, both-trait h2 strength, ascending FDR, descending |rg|, pair_id; no composite score.",
        ),
        "MECHANISTIC_FOLLOWUP": (
            lambda row: (
                str(local_by_key.get(pair_key(row), {}).get("priority_categories", "Z"))[0],
                novelty_order[audits[pair_key(row)]["novelty_classification"]],
                -float(row["abs_rg"]), float(row["fdr"]), pair_id(row),
            ),
            "Earliest explicit local-follow-up category A-F, literature class, descending |rg|, ascending FDR, pair_id; no composite score.",
        ),
    }
    output: list[dict[str, object]] = []
    for view, (key_function, rule) in views.items():
        ordered = sorted(significant, key=key_function)[:25]
        for rank, row in enumerate(ordered, 1):
            key = pair_key(row)
            audit = audits[key]
            local_row = local_by_key.get(key)
            bridge = row["sleep_trait"] in bridge_sleep and external_degrees.get(row["external_trait"], 0) >= 3
            output.append({
                "ranking_view": view,
                "rank": rank,
                "sleep_trait": row["sleep_trait"],
                "external_trait": row["external_trait"],
                "external_domain": row["external_domain"],
                "rg": row["rg"], "se": row["se"], "p": row["p"], "fdr": row["fdr"],
                "abs_rg": row["abs_rg"],
                "sleep_h2": row["sleep_h2"], "external_h2": row["external_h2"],
                "gcov_intercept": row["gcov_intercept"],
                "objective_or_subjective": row["sleep_measurement_category"],
                "known_or_novel": audit["novelty_classification"],
                "published_rg_if_known": audit["representative_published_rg"],
                "confidence_class": row["confidence_class"],
                "cross_domain_bridge": str(bridge),
                "negative_relationship": str(float(row["rg"]) < 0),
                "local_followup_priority": local_row["priority_categories"] if local_row else "NO",
                "replication_dataset_available": audit["replication_feasibility"],
                "main_reason_for_priority": (
                    str(local_row["main_reason_for_priority"]) if local_row
                    else f"Ranks highly in {view.lower().replace('_', ' ')} under the stated lexicographic rule."
                ),
                "ranking_rule": rule,
            })
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    analysis = root / "results/analysis"

    master_all = read_tsv(analysis / "phase1_master_analysis.tsv")
    primary_all = [row for row in master_all if row["primary_or_sensitivity"] == "PRIMARY_PHASE1"]
    significant = [row for row in primary_all if row["locked_primary_significant"] == "True"]
    if (len(primary_all), len(significant)) != (372, 153):
        fail("expected 372 primary pairs and 153 locked discoveries")
    confidence = {pair_key(row): row for row in read_tsv(analysis / "result_confidence.tsv")}
    for row in primary_all:
        row["confidence_class"] = confidence[pair_key(row)]["confidence_class"]
    audits = {pair_key(row): row for row in read_tsv(analysis / "literature_novelty_audit.tsv")}
    if set(audits) != {pair_key(row) for row in significant}:
        fail("literature audit does not exactly cover significant primary pairs")
    bridge_rows = read_tsv(analysis / "cross_domain_bridges.tsv")
    bridge_sleep = {
        row["trait_id"] for row in bridge_rows
        if row["bridge_type"] == "SLEEP_DOMAIN_BRIDGE" and row["bridge_tier"] == "BROAD_5PLUS_DOMAINS"
    }
    external_hubs = read_tsv(analysis / "external_trait_hubs.tsv")
    external_degrees = {row["external_trait"]: int(row["significant_sleep_connections"]) for row in external_hubs}

    novelty = select_novelty_rows(significant, audits)
    local = build_local_priorities(
        primary_all, significant, audits,
        read_tsv(analysis / "global_cancellation_candidates.tsv"),
        read_tsv(analysis / "measurement_mode_comparison.tsv"),
        bridge_sleep, external_degrees,
    )
    top = build_top_findings(significant, audits, local, bridge_sleep, external_degrees)
    if len(novelty) != 25 or len(top) != 125:
        fail(f"priority output invariant failed: novelty={len(novelty)} top={len(top)}")
    if any(row["external_trait"] in {"t2d", "melanoma"} for row in novelty + local + top):
        fail("sensitivity-only endpoint leaked into primary priority output")

    novelty_payload = tsv_text(novelty, list(novelty[0]))
    local_payload = tsv_text(local, list(local[0]))
    top_payload = tsv_text(top, list(top[0]))
    summary = {
        "schema_version": "sleep-atlas-phase1-priorities.1",
        "novel_connection_priority_rows": len(novelty),
        "novel_priority_tier_counts": {
            tier: sum(row["priority_tier"] == tier for row in novelty) for tier in TIER_ORDER
        },
        "local_followup_pair_count": len(local),
        "local_followup_category_counts": {
            letter: sum(letter in str(row["priority_categories"]).split(";") for row in local)
            for letter in "ABCDEF"
        },
        "top_findings_ranked_views": 5,
        "top_findings_rows_per_view": 25,
        "ranking_policy": "Every ranking is lexicographic and printed in its output row; no hidden composite score is used.",
        "primary_only": True,
        "output_content_hashes": {
            "novel_connection_priorities.tsv": hashlib.sha256(novelty_payload.encode()).hexdigest(),
            "local_followup_priorities.tsv": hashlib.sha256(local_payload.encode()).hexdigest(),
            "top_findings.tsv": hashlib.sha256(top_payload.encode()).hexdigest(),
        },
    }
    summary_payload = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    outputs = {
        analysis / "novel_connection_priorities.tsv": novelty_payload,
        analysis / "local_followup_priorities.tsv": local_payload,
        analysis / "top_findings.tsv": top_payload,
        analysis / "priority_summary.json": summary_payload,
    }
    if args.validate_only:
        for path, payload in outputs.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"priority output drifted: {path.relative_to(root)}")
        print(f"PHASE1_PRIORITIES_VALID novel={len(novelty)} local={len(local)} ranked={len(top)}")
        return 0
    for path, payload in outputs.items():
        atomic_text(path, payload)
    print(f"PHASE1_PRIORITIES_BUILT novel={len(novelty)} local={len(local)} ranked={len(top)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
