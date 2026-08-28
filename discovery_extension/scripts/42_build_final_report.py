#!/usr/bin/env python3
"""Build the Stage-17 scientific report without converting missing analyses to zero."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("discovery_extension")
BLOCKED = "NA_BLOCKED_UPSTREAM"
TOP_FIELDS = [
    "rank", "pair_id", "sleep_trait", "external_phenotype", "discovery_rg",
    "discovery_fdr", "replication_rg", "replication_p", "prior_literature_status",
    "local_loci", "pleiotropic_loci", "colocalized_genes",
    "strongest_mechanistic_evidence", "confidence_level", "analysis_status",
    "claim_limit",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "status": "PRESENT" if path.is_file() else "ABSENT_BLOCKED_UPSTREAM",
        "sha256": sha256(path) if path.is_file() else None,
    }


def pair_id(row: dict[str, str]) -> str:
    if row.get("pair_id"):
        return row["pair_id"]
    return f"{row['sleep_trait']}__{row['extension_trait_id']}"


def available_count(path: Path, count: int) -> tuple[str, str, str]:
    if not path.is_file():
        return BLOCKED, "NOT_ESTIMABLE", "Required result artifact is absent because acquisition is blocked."
    return str(count), "OBSERVED_COUNT", "Computed from the completed canonical result artifact."


def metric(
    order: int, name: str, value: str | int, value_type: str, status: str,
    definition: str, evidence: Path,
) -> dict[str, object]:
    return {
        "report_order": order, "metric": name, "value": value, "value_type": value_type,
        "analysis_status": status, "definition": definition, "evidence_path": str(evidence),
    }


def result_metric(order: int, name: str, path: Path, count: int, definition: str) -> dict[str, object]:
    value, value_type, note = available_count(path, count)
    return metric(
        order, name, value, value_type,
        "COMPLETE" if path.is_file() else "BLOCKED_UPSTREAM",
        note + " " + definition, path,
    )


def join_values(values: set[str]) -> str:
    cleaned = sorted(value for value in values if value and value not in {"NA", "NONE"})
    return ";".join(cleaned) if cleaned else "NONE"


def build_top_rows(paths: dict[str, Path]) -> list[dict[str, object]]:
    """Build a ranked table only from completed canonical result artifacts."""
    required = ("priority", "novelty")
    if any(not paths[name].is_file() for name in required):
        return []
    priority = {
        row["pair_id"]: row for row in read_tsv(paths["priority"])
        if row.get("priority_tier") in {"A", "B"}
        and row.get("pair_novelty_class") in {"APPARENTLY_NOVEL", "NO_DIRECT_RG_FOUND"}
    }
    novelty = {row["pair_id"]: row for row in read_tsv(paths["novelty"])}
    replication = {row["pair_id"]: row for row in read_tsv(paths["replication"])}
    local_loci: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(paths["local"]):
        if row.get("local_significance_status") == "FDR_SIGNIFICANT":
            local_loci[row["pair_id"]].add(row["locus_id"])
    pleiotropic_loci: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(paths["pleiotropy"]):
        pleiotropic_loci[row["pair_id"]].add(row["locus_id"])
    coloc_genes: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(paths["coloc"]):
        if (
            row.get("colocalization_interpretation") in {
                "SHARED_SIGNAL_MODEL_SUPPORTED", "SINGLE_SIGNAL_FALLBACK_SHARED_MODEL_SUPPORTED",
            }
            and row.get("prior_robust") == "True"
            and row.get("molecular_feature_id") not in {None, "", "NA"}
        ):
            coloc_genes[row["pair_id"]].add(row["molecular_feature_id"])
    mechanism: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(paths["mechanism"]):
        mechanism[row["pair_id"]].append(row)
    confidence_order = {
        "MULTI_LAYER_HUMAN_GENETIC_SUPPORT": 0, "COLOCALIZED_MOLECULAR_SUPPORT": 1,
        "ANNOTATION_ONLY_SUPPORT": 2, "INSUFFICIENT_OR_CONFLICTING": 3,
    }
    ranked = sorted(
        priority.values(),
        key=lambda row: (
            0 if row["priority_tier"] == "B" else 1,
            float(row["extension_fdr"]), -abs(float(row["rg"])), row["pair_id"],
        ),
    )
    output: list[dict[str, object]] = []
    for rank, row in enumerate(ranked, start=1):
        identity = row["pair_id"]
        rep = replication.get(identity, {})
        audit = novelty.get(identity, {})
        mech = sorted(
            mechanism.get(identity, []),
            key=lambda item: (
                confidence_order.get(item.get("confidence_class", ""), 99),
                -int(item.get("supported_required_edge_count", "0")),
                item.get("gene_id", ""),
            ),
        )
        output.append({
            "rank": rank, "pair_id": identity, "sleep_trait": row["sleep_trait"],
            "external_phenotype": row["phenotype_name"], "discovery_rg": row["rg"],
            "discovery_fdr": row["extension_fdr"],
            "replication_rg": rep.get("replication_rg", "NA_NOT_AVAILABLE"),
            "replication_p": rep.get("replication_p", "NA_NOT_AVAILABLE"),
            "prior_literature_status": audit.get("novelty_class", row["pair_novelty_class"]),
            "local_loci": join_values(local_loci[identity]),
            "pleiotropic_loci": join_values(pleiotropic_loci[identity]),
            "colocalized_genes": join_values(coloc_genes[identity]),
            "strongest_mechanistic_evidence": (
                mech[0]["strongest_mechanistic_evidence"] if mech else "NONE"
            ),
            "confidence_level": (
                mech[0]["confidence_class"] if mech else row["priority_tier"]
            ),
            "analysis_status": "REPORTABLE_WITH_LOCKED_CLAIM_LIMITS",
            "claim_limit": "DISCOVERY_ASSOCIATION_AND_SUPPORTING_EVIDENCE_NOT_CAUSAL_PROOF",
        })
    return output


def main() -> None:
    paths = {
        "universe": ROOT / "results/panukbb_eur_eligible_universe.tsv",
        "pool": ROOT / "results/candidate_pool.tsv",
        "panel": ROOT / "config/candidate_traits.tsv",
        "h2": ROOT / "results/ldsc/extension_trait_readiness.tsv",
        "rg": ROOT / "results/ldsc/extension_rg_matrix.tsv",
        "novelty": ROOT / "results/novelty/extension_novelty_audit.tsv",
        "priority": ROOT / "results/prioritization/novel_hit_priority.tsv",
        "replication": ROOT / "results/replication/replication_results.tsv",
        "local": ROOT / "results/local/local_rg_results.tsv",
        "pleiotropy": ROOT / "results/pleiotropy/novel_shared_loci.tsv",
        "coloc": ROOT / "results/fine_mapping/fine_mapping_colocalization.tsv",
        "mechanism": ROOT / "results/mechanism/mechanistic_synthesis.tsv",
    }
    lock = json.loads((ROOT / "config/extension_panel.lock.json").read_text(encoding="utf-8"))
    preflight = json.loads((ROOT / "provenance/acquisition_preflight.json").read_text(encoding="utf-8"))
    universe, pool, panel = (read_tsv(paths[name]) for name in ("universe", "pool", "panel"))
    if len(panel) < 75 or len(panel) > 125:
        raise SystemExit(f"ERROR: locked panel size is outside the brief's 75-125 range: {len(panel)}")
    h2, rg, novelty = (read_tsv(paths[name]) for name in ("h2", "rg", "novelty"))
    replication, local, pleiotropy = (read_tsv(paths[name]) for name in ("replication", "local", "pleiotropy"))
    coloc, mechanism = (read_tsv(paths[name]) for name in ("coloc", "mechanism"))

    h2_value, h2_type, h2_note = available_count(
        paths["h2"], sum(row.get("primary_rg_eligibility") == "PRIMARY_PASS" for row in h2),
    )
    test_value, test_type, test_note = available_count(paths["rg"], len(rg))
    fdr_value, fdr_type, fdr_note = available_count(
        paths["rg"], sum(float(row["extension_fdr"]) < 0.05 for row in rg),
    )
    novelty_counts = {
        label: sum(row.get("novelty_class") == label for row in novelty)
        for label in ("KNOWN_REPLICATION", "PARTIAL_EXTENSION", "NO_DIRECT_RG_FOUND", "APPARENTLY_NOVEL")
    }
    coloc_signals = {
        (row["comparison_id"], row.get("signal1", "NA"), row.get("signal2", "NA"))
        for row in coloc
        if row.get("colocalization_interpretation") in {
            "SHARED_SIGNAL_MODEL_SUPPORTED", "SINGLE_SIGNAL_FALLBACK_SHARED_MODEL_SUPPORTED",
        } and row.get("prior_robust") == "True"
    }
    mechanism_pairs = {
        row["pair_id"] for row in mechanism
        if row.get("confidence_class") != "INSUFFICIENT_OR_CONFLICTING"
        and row.get("strongest_mechanistic_evidence") not in {None, "", "NA", "NONE"}
    }
    rows = [
        metric(1, "candidates_considered", len(pool), "OBSERVED_PREANALYSIS_COUNT", "COMPLETE",
               f"Candidate pool after prospective screening; eligible source universe={len(universe)}.", paths["pool"]),
        metric(2, "selected_before_analysis", len(panel), "OBSERVED_PREANALYSIS_COUNT", "COMPLETE",
               "Independently locked extension traits; no extension rg result existed at lock time.", paths["panel"]),
        metric(3, "passing_rerun_h2_qc", h2_value, h2_type, "COMPLETE" if paths["h2"].is_file() else "BLOCKED_UPSTREAM",
               h2_note + " Source-screen h2 is not substituted for the required rerun.", paths["h2"]),
        metric(4, "sleep_x_trait_tests_executed", test_value, test_type, "COMPLETE" if paths["rg"].is_file() else "BLOCKED_UPSTREAM",
               test_note + f" Locked planned maximum={lock['planned_raw_rg_test_count']} (12 x {len(panel)}).", paths["rg"]),
        metric(5, "extension_fdr_significant_pairs", fdr_value, fdr_type, "COMPLETE" if paths["rg"].is_file() else "BLOCKED_UPSTREAM",
               fdr_note + " BH FDR is confined to the executed primary extension family.", paths["rg"]),
        result_metric(6, "already_established_pairs", paths["novelty"], novelty_counts["KNOWN_REPLICATION"],
                      "Pair-level novelty class KNOWN_REPLICATION."),
        result_metric(7, "partial_extension_pairs", paths["novelty"], novelty_counts["PARTIAL_EXTENSION"],
                      "Pair-level novelty class PARTIAL_EXTENSION."),
        result_metric(8, "no_direct_prior_rg_found_pairs", paths["novelty"], novelty_counts["NO_DIRECT_RG_FOUND"],
                      "Pair-level novelty class NO_DIRECT_RG_FOUND after completed targeted review."),
        result_metric(9, "apparently_novel_pairs", paths["novelty"], novelty_counts["APPARENTLY_NOVEL"],
                      "Pair-level novelty class APPARENTLY_NOVEL after several targeted searches."),
        result_metric(10, "independently_replicated_pairs", paths["replication"],
                      len({row["pair_id"] for row in replication if row.get("replication_class") == "REPLICATED"}),
                      "Unique pairs classified REPLICATED in the locked non-overlapping replication family."),
        result_metric(11, "pairs_with_significant_local_correlations", paths["local"],
                      len({row["pair_id"] for row in local if row.get("local_significance_status") == "FDR_SIGNIFICANT"}),
                      "Unique pairs with at least one FDR_SIGNIFICANT canonical local-rg row."),
        result_metric(12, "pairs_with_pleiotropic_loci", paths["pleiotropy"],
                      len({row["pair_id"] for row in pleiotropy}),
                      "Unique pairs with at least one canonical PLACO+ locus; not proof of a shared causal variant."),
        result_metric(13, "prior_robust_colocalized_signals", paths["coloc"], len(coloc_signals),
                      "Unique signal comparisons supporting a shared model across the locked prior grid."),
        result_metric(14, "pairs_with_candidate_genes_or_mechanisms", paths["mechanism"], len(mechanism_pairs),
                      "Unique pairs with non-insufficient evidence-linked synthesis; annotation is noncausal."),
    ]
    counts_path = ROOT / "results/final_extension_counts.tsv"
    counts_path.parent.mkdir(parents=True, exist_ok=True)
    with counts_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    top_rows = build_top_rows(paths)
    top_path = ROOT / "results/top_novel_discoveries.tsv"
    with top_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=TOP_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(top_rows)

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    report_path = ROOT / "final_report.md"
    table_lines = [
        "| # | Required count | Value | Status |",
        "|---:|---|---:|---|",
        *[
            f"| {row['report_order']} | {str(row['metric']).replace('_', ' ')} | {row['value']} | {row['analysis_status']} |"
            for row in rows
        ],
    ]
    top_markdown = (
        "No rows are reported. The header-only canonical table is "
        "`results/top_novel_discoveries.tsv`; its empty body means results are unavailable, not that zero discoveries were observed."
        if not top_rows else
        f"{len(top_rows)} ranked rows are available in `results/top_novel_discoveries.tsv`."
    )
    report_path.write_text(
        "# Novelty-Enriched Phenome Discovery Extension — final scientific report\n\n"
        f"Generated: {generated}\n\n"
        "## Publication-readiness verdict\n\n"
        "**BLOCKED BEFORE FULL-RESOLUTION ACQUISITION; NO EXTENSION FINDINGS EXIST.**\n\n"
        "The immutable 45-trait/396-pair core remains separate and unchanged. The extension completed its "
        "eligible-universe review, 242-candidate screen, independent 100-trait lock, novelty-prescreen inventory, "
        "source/schema contracts, analytic software locks, synthetic verification, acceptance gates, and pre-result "
        "adversarial review. It did not acquire or analyze the real dense extension statistics.\n\n"
        f"The exact compressed sources total {preflight['compressed_source_gib']} GiB. The locked {preflight['safety_factor']} "
        f"safety factor requires {preflight['required_free_gib']} GiB free; the recorded preflight found "
        f"{preflight['available_free_gib']} GiB ({preflight['shortfall_gib']} GiB short). No bulk download was started.\n\n"
        "Consequently, all result-dependent counts below are `NA_BLOCKED_UPSTREAM`, not zero. The 1,200 figure is "
        "the locked planned maximum (12 sleep traits × 100 extension traits), not an executed-test count.\n\n"
        "## Required scientific summary\n\n"
        + "\n".join(table_lines)
        + "\n\n## Top Novel Discoveries\n\n"
        + top_markdown
        + "\n\n## Claim boundary\n\n"
        "There is no extension h² rerun, global genetic-correlation result, extension-FDR hit, pair-level novelty "
        "claim, independent replication, local correlation, pleiotropic locus, fine-mapped signal, colocalization, "
        "candidate gene, mechanistic inference, or findings-level adversarial review to publish. Synthetic outputs "
        "validate software contracts only and are excluded from every count above.\n\n"
        "Clearing the storage gate requires a different volume with at least the locked required free space (practically, "
        "approximately 260 GiB free to retain a modest operational margin), followed by a fresh preflight. Thresholds, "
        "panel membership, and the isolated extension FDR family must not be weakened to fit the current volume.\n",
        encoding="utf-8",
    )

    provenance_path = ROOT / "provenance/final_report.json"
    provenance = {
        "schema_version": "1.0.0", "generated_utc": generated,
        "publication_readiness": "BLOCKED_BEFORE_FULL_RESOLUTION_ACQUISITION",
        "real_result_rows_reported": len(top_rows),
        "missing_result_policy": "NA_BLOCKED_UPSTREAM_NEVER_ZERO",
        "planned_test_count": lock["planned_raw_rg_test_count"],
        "executed_test_count": len(rg) if paths["rg"].is_file() else None,
        "inputs": {name: artifact(path) for name, path in paths.items()},
        "outputs": {
            "counts": artifact(counts_path), "top_novel_discoveries": artifact(top_path),
            "report": artifact(report_path),
        },
        "warning": "A header-only top-discoveries table denotes unavailable results, not an observed count of zero.",
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"FINAL_EXTENSION_REPORT_OK selected={len(panel)} planned_tests={lock['planned_raw_rg_test_count']} "
        f"executed_tests={len(rg) if paths['rg'].is_file() else BLOCKED} top_rows={len(top_rows)}"
    )


if __name__ == "__main__":
    main()
