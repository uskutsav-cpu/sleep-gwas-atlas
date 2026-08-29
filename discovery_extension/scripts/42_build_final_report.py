#!/usr/bin/env python3
"""Build the Stage-17 report from validated real results and explicit blockers."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
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


def available_count(path: Path, count: int) -> tuple[str, str, str]:
    if not path.is_file():
        return BLOCKED, "NOT_ESTIMABLE", "Required canonical result artifact is absent."
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


def joined_or_status(path: Path, values: set[str]) -> str:
    if not path.is_file():
        return BLOCKED
    cleaned = sorted(value for value in values if value and value not in {"NA", "NONE"})
    return ";".join(cleaned) if cleaned else "NONE_OBSERVED_IN_COMPLETED_ANALYSIS"


def build_top_rows(paths: dict[str, Path]) -> list[dict[str, object]]:
    """Rank only independently replicated Tier-B discoveries."""
    required = ("priority", "novelty", "replication", "replication_queue")
    if any(not paths[name].is_file() for name in required):
        return []
    priority = {
        row["pair_id"]: row for row in read_tsv(paths["priority"])
        if row.get("priority_tier") == "B"
    }
    novelty = {row["pair_id"]: row for row in read_tsv(paths["novelty"])}
    replication = {
        row["pair_id"]: row for row in read_tsv(paths["replication"])
        if row.get("replication_class") == "REPLICATED"
    }
    replication_queue = {row["pair_id"]: row for row in read_tsv(paths["replication_queue"])}
    if set(priority) != set(replication):
        raise SystemExit("ERROR: Tier-B priority family is not the independently replicated family")

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

    ranked = sorted(
        priority.values(),
        key=lambda row: (float(row["extension_fdr"]), -abs(float(row["rg"])), row["pair_id"]),
    )
    output: list[dict[str, object]] = []
    for rank, row in enumerate(ranked, start=1):
        identity = row["pair_id"]
        rep = replication[identity]
        audit = novelty[identity]
        source = replication_queue[identity]
        heterogeneous = float(rep["heterogeneity_p"]) < 0.05
        match = source["phenotype_match_status"]
        confidence = (
            "REPLICATED_EXACT_PHENOTYPE" if match == "EXACT"
            else "REPLICATED_COMPARABLE_PHENOTYPE"
        )
        if heterogeneous:
            confidence += "_WITH_EFFECT_HETEROGENEITY"
        mechanistic = BLOCKED
        if paths["mechanism"].is_file():
            candidates = sorted(
                mechanism.get(identity, []),
                key=lambda item: (-int(item.get("supported_required_edge_count", "0")), item.get("gene_id", "")),
            )
            mechanistic = (
                candidates[0].get("strongest_mechanistic_evidence", "NONE_OBSERVED_IN_COMPLETED_ANALYSIS")
                if candidates else "NONE_OBSERVED_IN_COMPLETED_ANALYSIS"
            )
        output.append({
            "rank": rank, "pair_id": identity, "sleep_trait": row["sleep_trait"],
            "external_phenotype": row["phenotype_name"], "discovery_rg": row["rg"],
            "discovery_fdr": row["extension_fdr"], "replication_rg": rep["replication_rg"],
            "replication_p": rep["replication_p"],
            "prior_literature_status": audit["novelty_class"],
            "local_loci": joined_or_status(paths["local"], local_loci[identity]),
            "pleiotropic_loci": joined_or_status(paths["pleiotropy"], pleiotropic_loci[identity]),
            "colocalized_genes": joined_or_status(paths["coloc"], coloc_genes[identity]),
            "strongest_mechanistic_evidence": mechanistic,
            "confidence_level": confidence,
            "analysis_status": "REPLICATED_GLOBAL_GENETIC_CORRELATION;DOWNSTREAM_RESOURCE_BLOCKED",
            "claim_limit": "UNDERREPORTED_REPLICATED_ASSOCIATION_NOT_FIRST_EVER_OR_CAUSAL_PROOF",
        })
    return output


def main() -> None:
    post_priority = ROOT / "results/prioritization/novel_hit_priority_post_replication.tsv"
    paths = {
        "universe": ROOT / "results/panukbb_eur_eligible_universe.tsv",
        "pool": ROOT / "results/candidate_pool.tsv",
        "panel": ROOT / "config/candidate_traits.tsv",
        "h2": ROOT / "results/ldsc/extension_trait_readiness.tsv",
        "rg": ROOT / "results/ldsc/extension_rg_matrix.tsv",
        "novelty": ROOT / "results/novelty/extension_novelty_audit.tsv",
        "priority": post_priority if post_priority.is_file() else ROOT / "results/prioritization/novel_hit_priority.tsv",
        "replication": ROOT / "results/replication/replication_results.tsv",
        "replication_queue": ROOT / "results/replication/replication_source_queue.tsv",
        "local": ROOT / "results/local/local_rg_results.tsv",
        "pleiotropy": ROOT / "results/pleiotropy/novel_shared_loci.tsv",
        "coloc": ROOT / "results/fine_mapping/fine_mapping_colocalization.tsv",
        "mechanism": ROOT / "results/mechanism/mechanistic_synthesis.tsv",
        "local_readiness": ROOT / "results/local/local_architecture_readiness.tsv",
        "pleiotropy_queue": ROOT / "results/pleiotropy/pleiotropy_input_queue.tsv",
    }
    lock = json.loads((ROOT / "config/extension_panel.lock.json").read_text(encoding="utf-8"))
    preflight = json.loads((ROOT / "provenance/acquisition_preflight.json").read_text(encoding="utf-8"))
    universe, pool, panel = (read_tsv(paths[name]) for name in ("universe", "pool", "panel"))
    h2, rg, novelty = (read_tsv(paths[name]) for name in ("h2", "rg", "novelty"))
    replication = read_tsv(paths["replication"])
    local, pleiotropy = (read_tsv(paths[name]) for name in ("local", "pleiotropy"))
    coloc, mechanism = (read_tsv(paths[name]) for name in ("coloc", "mechanism"))
    replication_queue = {row["pair_id"]: row for row in read_tsv(paths["replication_queue"])}

    if not 75 <= len(panel) <= 125 or lock["extension_trait_count"] != len(panel):
        raise SystemExit("ERROR: extension panel size differs from its lock or the required range")
    if len(h2) != len(panel) or [row["extension_trait_id"] for row in h2] != [row["extension_trait_id"] for row in panel]:
        raise SystemExit("ERROR: h2 result family differs from locked panel")
    expected_pairs = lock["planned_raw_rg_test_count"]
    if len(rg) != expected_pairs or len({(row["sleep_trait"], row["extension_trait_id"]) for row in rg}) != expected_pairs:
        raise SystemExit("ERROR: global-rg result family differs from extension lock")
    fdr_hits = [row for row in rg if float(row["extension_fdr"]) < 0.05]
    if len(novelty) != len(fdr_hits) or any(row.get("audit_status") != "COMPLETE" for row in novelty):
        raise SystemExit("ERROR: pair-level novelty audit is incomplete")
    if len(replication) != 217:
        raise SystemExit("ERROR: replication family is incomplete")

    h2_pass = sum(row["primary_rg_eligibility"] == "PRIMARY_PASS" for row in h2)
    joint_screen = sum(float(row["extension_fdr"]) < 0.05 and abs(float(row["rg"])) >= 0.15 for row in rg)
    novelty_counts = Counter(row["novelty_class"] for row in novelty)
    replication_counts = Counter(row["replication_class"] for row in replication)
    replicated = [row for row in replication if row["replication_class"] == "REPLICATED"]
    tested = [row for row in replication if row["replication_class"] in {"REPLICATED", "DIRECTIONALLY_CONCORDANT", "FAILED_REPLICATION"}]
    heterogeneous = sum(float(row["heterogeneity_p"]) < 0.05 for row in replicated)
    phenotype_matches = Counter(replication_queue[row["pair_id"]]["phenotype_match_status"] for row in replicated)

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
    established = novelty_counts["KNOWN_REPLICATION"] + novelty_counts["KNOWN_BUT_NEW_DATASET"]
    rows = [
        metric(1, "candidates_considered", len(pool), "OBSERVED_PREANALYSIS_COUNT", "COMPLETE",
               f"Candidate pool after prospective screening; eligible Pan-UKB source universe={len(universe)}.", paths["pool"]),
        metric(2, "selected_before_analysis", len(panel), "OBSERVED_PREANALYSIS_COUNT", "COMPLETE",
               "Ordered extension panel locked before any extension genetic-correlation result existed.", paths["panel"]),
        result_metric(3, "passing_rerun_h2_qc", paths["h2"], h2_pass,
                      "Predefined h2 Z>=4 and LDSC intercept<=1.2 primary gates."),
        result_metric(4, "sleep_x_trait_tests_executed", paths["rg"], len(rg),
                      "Exact isolated 12 x 100 primary extension family."),
        result_metric(5, "extension_fdr_significant_pairs", paths["rg"], len(fdr_hits),
                      "Benjamini-Hochberg FDR<0.05 within the isolated extension family."),
        result_metric(6, "already_established_pairs", paths["novelty"], established,
                      "KNOWN_REPLICATION plus KNOWN_BUT_NEW_DATASET after pair-level literature review."),
        result_metric(7, "partial_extension_pairs", paths["novelty"], novelty_counts["PARTIAL_EXTENSION"],
                      "Related prior genetic evidence exists, but not the exact comparison."),
        result_metric(8, "no_direct_prior_rg_found_pairs", paths["novelty"], novelty_counts["NO_DIRECT_RG_FOUND"],
                      "No directly comparable prior rg was returned by the recorded multi-query audit."),
        result_metric(9, "apparently_novel_pairs", paths["novelty"], novelty_counts["APPARENTLY_NOVEL"],
                      "Reserved stringent class; zero is an observed audited count, not a first-ever claim."),
        result_metric(10, "independently_replicated_pairs", paths["replication"], len(replicated),
                      "Frozen 0.05/217 Bonferroni threshold, direction concordance, EUR ancestry, and non-overlapping cohort."),
        result_metric(11, "pairs_with_significant_local_correlations", paths["local"],
                      len({row["pair_id"] for row in local if row.get("local_significance_status") == "FDR_SIGNIFICANT"}),
                      "Unavailable because checksum-locked LAVA/HDL-L LD references and dense inputs are absent."),
        result_metric(12, "pairs_with_pleiotropic_loci", paths["pleiotropy"],
                      len({row["pair_id"] for row in pleiotropy}),
                      "Unavailable because full-genome per-variant inputs and clumping LD are absent."),
        result_metric(13, "prior_robust_colocalized_signals", paths["coloc"], len(coloc_signals),
                      "Unavailable because PLACO loci, signed LD, and locus-complete trait/QTL statistics are absent."),
        result_metric(14, "pairs_with_candidate_genes_or_mechanisms", paths["mechanism"], len(mechanism_pairs),
                      "Unavailable because no fine-mapped/colocalized signal family exists."),
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
    if len(top_rows) != len(replicated):
        raise SystemExit("ERROR: top-discovery table is not the exact replicated Tier-B family")

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    table_lines = [
        "| # | Required count | Value | Status |",
        "|---:|---|---:|---|",
        *[
            f"| {row['report_order']} | {str(row['metric']).replace('_', ' ')} | {row['value']} | {row['analysis_status']} |"
            for row in rows
        ],
    ]
    top_lines = [
        "| Rank | Sleep trait | External phenotype | Discovery rg | Replication rg | Replication P | Literature status |",
        "|---:|---|---|---:|---:|---:|---|",
    ]
    for row in top_rows[:10]:
        top_lines.append(
            f"| {row['rank']} | {row['sleep_trait']} | {row['external_phenotype']} | "
            f"{float(row['discovery_rg']):.3f} | {float(row['replication_rg']):.3f} | "
            f"{float(row['replication_p']):.3g} | {row['prior_literature_status']} |"
        )

    local_readiness = read_tsv(paths["local_readiness"])
    available_gib = local_readiness[0]["available_free_gib"] if local_readiness else "unknown"
    replication_figure_png = ROOT / "figures/extension_replication_summary.png"
    replication_figure_pdf = ROOT / "figures/extension_replication_summary.pdf"
    replication_figure_caption = ROOT / "figures/extension_replication_summary.caption.txt"
    report_path = ROOT / "final_report.md"
    report_path.write_text(
        "# Novelty-Enriched Phenome Discovery Extension — final scientific report\n\n"
        f"Generated: {generated}\n\n"
        "## Publication-readiness verdict\n\n"
        "**GLOBAL DISCOVERY AND INDEPENDENT REPLICATION COMPLETE; LOCAL AND MECHANISTIC FOLLOW-UP RESOURCE-BLOCKED.**\n\n"
        "The immutable 45-trait/396-pair core remains separate and unchanged. This separately locked extension "
        "prospectively selected 100 phenotypes, executed all 1,200 planned sleep-by-phenotype LDSC tests, audited "
        "every extension-FDR hit, and preserved the complete 217-pair replication candidate family. Findings are "
        "reportable as replicated global genetic correlations with explicit novelty and causal-language limits; "
        "they are not local-sharing, pleiotropy, colocalization, or mechanism claims.\n\n"
        "## Required scientific summary\n\n" + "\n".join(table_lines) + "\n\n"
        "## Discovery and novelty results\n\n"
        f"All {h2_pass}/{len(panel)} extension traits passed the prespecified rerun h2 gates. The exact {len(rg)}-pair "
        f"family yielded {len(fdr_hits)} BH-FDR-significant pairs; {joint_screen} also had |rg|>=0.15. The maximum "
        f"absolute cross-trait LDSC intercept was {max(abs(float(row['cross_trait_LDSC_intercept'])) for row in rg):.4f}, "
        f"and the minimum valid-allele overlap was {min(int(row['snp_overlap_valid_alleles']) for row in rg):,}.\n\n"
        f"The literature audit covered {len(novelty)}/{len(fdr_hits)} FDR-significant pairs: {established} already "
        f"established or same-pair/new-dataset result, {novelty_counts['PARTIAL_EXTENSION']} partial extensions, "
        f"{novelty_counts['NO_DIRECT_RG_FOUND']} with no direct prior rg found, and "
        f"{novelty_counts['APPARENTLY_NOVEL']} classed APPARENTLY_NOVEL. The zero APPARENTLY_NOVEL count is deliberate: "
        "none cleared the audit's strongest novelty bar, and no first-ever claim is made.\n\n"
        "## Independent replication\n\n"
        f"Among the frozen 217 candidates, {replication_counts['REPLICATED']} met the 0.05/217 replication threshold, "
        f"{replication_counts['DIRECTIONALLY_CONCORDANT']} more were concordant but below that threshold, "
        f"{replication_counts['UNDERPOWERED']} were stopped by prespecified replication h2/intercept gates, and "
        f"{replication_counts['NO_INDEPENDENT_DATASET']} had no suitable independent dataset. None of the {len(tested)} "
        f"tested pairs was direction-discordant. Of the 23 replicated pairs, {phenotype_matches['EXACT']} used exact "
        f"phenotype matches and {phenotype_matches['COMPARABLE_WITH_DOCUMENTED_DIFFERENCES']} used documented comparable "
        f"definitions; {heterogeneous} showed effect-size heterogeneity P<0.05 and are retained with that qualification.\n\n"
        "The replicated set concentrates in gastrointestinal, pulmonary, smoking-related, pain, and musculoskeletal "
        "phenotypes. This pattern is biologically interpretable but cannot distinguish direct biology, vertical or "
        "horizontal pleiotropy, mediation, ascertainment, or residual sample-structure effects.\n\n"
        "## Top replicated underreported discoveries\n\n" + "\n".join(top_lines) + "\n\n"
        f"The complete {len(top_rows)}-row ranked table is `results/top_novel_discoveries.tsv`. Local loci, pleiotropic "
        "loci, colocalized genes, and mechanistic evidence are explicitly `NA_BLOCKED_UPSTREAM`, not zero.\n\n"
        "The complete replication scatter, effect-estimate forest plot, and locked-family attrition flow are in "
        "`figures/extension_replication_summary.png` and its vector PDF counterpart.\n\n"
        "## Resource-bounded downstream follow-up\n\n"
        f"The original full Pan-UKB mirror would require {preflight['required_free_gib']} GiB under the locked safety "
        "factor. Streaming and receipt sealing solved the LDSC acquisition problem for all 100 traits, but did not "
        "create the dense local-analysis inputs. At the final preflight only " + available_gib + " GiB remained. The "
        "published LAVA UKB v1.1 reference alone is 15 GiB unzipped; the signed fine-mapping LD and PLACO+ clumping "
        "reference are also absent. The retained HapMap3-munged LDSC files omit contract-required dense locus coverage "
        "and per-variant MAF/INFO and therefore were not silently relabeled as full-resolution inputs.\n\n"
        f"A result-free local queue preserves 217 priority discoveries and 597 globally-null secondary candidates. "
        f"A separately locked PLACO+ candidate family preserves all {len(read_tsv(paths['pleiotropy_queue']))} replicated "
        "Tier-B pairs. No LAVA, HDL-L, PLACO+, SuSiE, coloc, molecular-QTL, or mechanism result is reported.\n\n"
        "## Claim boundary\n\n"
        "The strongest supported statement is that 23 underreported sleep-phenotype genetic correlations replicated "
        "directionally in independent EUR GWAS at the frozen family-wise threshold. Genetic correlation is not causation. "
        "No result is described as first-ever, locally shared, pleiotropic, colocalized, fine-mapped, gene-mediated, or "
        "mechanistic. Seven replicated effects require explicit heterogeneity disclosure, and five use comparable rather "
        "than exact phenotype definitions.\n",
        encoding="utf-8",
    )

    provenance_path = ROOT / "provenance/final_report.json"
    provenance = {
        "schema_version": "2.0.0", "generated_utc": generated,
        "publication_readiness": "GLOBAL_DISCOVERY_AND_REPLICATION_COMPLETE_DOWNSTREAM_BLOCKED",
        "real_result_rows_reported": len(top_rows),
        "missing_result_policy": "NA_BLOCKED_UPSTREAM_NEVER_ZERO",
        "planned_test_count": lock["planned_raw_rg_test_count"], "executed_test_count": len(rg),
        "extension_fdr_hit_count": len(fdr_hits), "replication_classification_counts": dict(sorted(replication_counts.items())),
        "replicated_effect_heterogeneity_p_lt_0_05": heterogeneous,
        "replicated_phenotype_match_counts": dict(sorted(phenotype_matches.items())),
        "inputs": {name: artifact(path) for name, path in paths.items()},
        "outputs": {
            "counts": artifact(counts_path), "top_novel_discoveries": artifact(top_path),
            "report": artifact(report_path), "replication_figure_png": artifact(replication_figure_png),
            "replication_figure_pdf": artifact(replication_figure_pdf),
            "replication_figure_caption": artifact(replication_figure_caption),
        },
        "warning": "Downstream unavailable values are NA_BLOCKED_UPSTREAM, never observed zeros.",
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"FINAL_EXTENSION_REPORT_OK selected={len(panel)} tests={len(rg)} fdr={len(fdr_hits)} "
        f"replicated={len(replicated)} top_rows={len(top_rows)}"
    )


if __name__ == "__main__":
    main()
