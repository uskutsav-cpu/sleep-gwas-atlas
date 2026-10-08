#!/usr/bin/env python3
"""Exercise full 100-trait h2 and 12x98 rg collation in isolation."""

from __future__ import annotations

import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MARKER = "SYNTHETIC EXTENSION LDSC OUTPUT - NOT REAL RESULTS"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    panel = read_tsv(ROOT / "discovery_extension/config/candidate_traits.tsv")
    sleeps = [
        row["trait_id"] for row in read_tsv(ROOT / "config/analysis_panel.tsv")
        if row["domain"] == "sleep"
    ]
    ids = [row["extension_trait_id"] for row in panel]
    with tempfile.TemporaryDirectory(prefix="synthetic_extension_ldsc_") as temporary:
        work = Path(temporary)
        h2_logs, rg_logs = work / "h2", work / "rg"
        h2_logs.mkdir()
        rg_logs.mkdir()
        # Explicitly synthetic, disposable gate fixture. No original h2 result
        # is required or written. Trait identifiers come from the public panel.
        core_h2 = work / "synthetic_core_sleep_h2.tsv"
        core_h2_fields = ["trait", "h2", "se", "z", "intercept", "verdict", "analysis_status"]
        core_h2_rows = [
            {"trait": trait, "h2": 0.2, "se": 0.02, "z": 10,
             "intercept": 1.0, "verdict": "PASS", "analysis_status": "SYNTHETIC_TEST_ONLY"}
            for trait in sleeps
        ]
        if len(core_h2_rows) != 12 or len({row["trait"] for row in core_h2_rows}) != 12:
            raise SystemExit("ERROR: synthetic core h2 fixture is not the exact 12-sleep family")

        def write_core_h2(rows: list[dict[str, object]]) -> None:
            with core_h2.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=core_h2_fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

        write_core_h2(core_h2_rows)
        h2_out = work / "h2.tsv"
        h2_failed_out = work / "h2_failed.tsv"
        rg_out = work / "rg.tsv"
        universe_out = work / "universe.tsv"

        for index, trait_id in enumerate(ids):
            h2, se, intercept = 0.2, 0.02, 1.0
            if index == 0:
                se = 0.1
            if index == 1:
                intercept = 1.3
            (h2_logs / f"h2_{trait_id}.log").write_text(
                f"{MARKER}\nRead summary statistics for 1200000 SNPs.\n"
                f"After merging with regression SNP LD, 1100000 SNPs remain.\n"
                f"Total Observed scale h2: {h2} ({se})\n"
                f"Intercept: {intercept} (0.01)\nMean Chi^2: 1.2\nLambda GC: 1.1\nRatio: 0.1\n"
            )
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/13_collate_extension_ldsc.py"),
                "--mode", "h2", "--logdir", str(h2_logs), "--out", str(h2_out),
                "--h2-failed-out", str(h2_failed_out),
                "--provenance-out", str(work / "h2.json"),
                "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
                "--core", str(ROOT / "config/analysis_panel.tsv"),
            ], check=True,
        )
        h2_rows = read_tsv(h2_out)
        if len(h2_rows) != 100:
            raise SystemExit(f"ERROR: synthetic h2 rows={len(h2_rows)} expected 100")
        if sum(row["primary_rg_eligibility"] == "PRIMARY_PASS" for row in h2_rows) != 98:
            raise SystemExit("ERROR: synthetic h2 gate did not yield exactly 98 primary traits")
        if len(read_tsv(h2_failed_out)) != 2:
            raise SystemExit("ERROR: synthetic separate h2-failed table is not exactly two rows")

        passing = ids[2:]
        header = "p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se"
        for sleep_index, sleep in enumerate(sleeps):
            p_values = [1e-8 if sleep_index == 0 and index == 0 else 0.5 for index in range(len(passing))]
            lines = [MARKER]
            for index, trait_id in enumerate(passing):
                lines.extend([
                    f"Computing rg for phenotype {index + 2}/{len(passing) + 1}",
                    f"Reading summary statistics from discovery_extension/data/munged/{trait_id}.sumstats.gz ...",
                    "Read summary statistics for 1150000 SNPs.",
                    "After merging with summary statistics, 1100000 SNPs remain.",
                    "1090000 SNPs with valid alleles.",
                ])
            lines.extend(f"P: {p_value}" for p_value in p_values)
            lines.extend(["Summary of Genetic Correlation Results", header])
            for index, (trait_id, p_value) in enumerate(zip(passing, p_values)):
                rg = 0.2 if sleep_index == 0 and index == 0 else (-0.25 if sleep_index == 1 and index == 0 else 0.01)
                lines.append(
                    f"data/munged/{sleep}.sumstats.gz "
                    f"discovery_extension/data/munged/{trait_id}.sumstats.gz "
                    f"{rg} 0.03 {rg/0.03:.5g} {p_value} 0.2 0.02 1.0 0.01 0.01 0.001"
                )
            lines.append("")
            (rg_logs / f"rg_{sleep}.log").write_text("\n".join(lines))

        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/13_collate_extension_ldsc.py"),
                "--mode", "rg", "--logdir", str(rg_logs), "--h2", str(h2_out),
                "--out", str(rg_out), "--pair-universe-out", str(universe_out),
                "--provenance-out", str(work / "rg.json"),
                "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
                "--core", str(ROOT / "config/analysis_panel.tsv"),
            ], check=True,
        )
        rg_rows, universe = read_tsv(rg_out), read_tsv(universe_out)
        if len(rg_rows) != 1176 or len(universe) != 1200:
            raise SystemExit(f"ERROR: rg rows={len(rg_rows)} universe={len(universe)}")
        if sum(row["pair_status"] == "H2_FAILED_PRIMARY_EXCLUSION" for row in universe) != 24:
            raise SystemExit("ERROR: synthetic pair universe did not preserve 24 h2 exclusions")
        top = next(row for row in rg_rows if row["sleep_trait"] == sleeps[0] and row["extension_trait_id"] == passing[0])
        if float(top["extension_fdr"]) >= 0.05 or top["initial_screen_status"] != "PAIR_NOVELTY_AUDIT_REQUIRED":
            raise SystemExit(f"ERROR: synthetic isolated FDR/priority failed: {top}")
        if top["cross_trait_LDSC_intercept"] != "0.01" or top["snp_overlap_valid_alleles"] != "1090000" or top["ancestry"] != "EUR":
            raise SystemExit(f"ERROR: synthetic pair diagnostics are incomplete: {top}")
        positive_out, negative_out = work / "positive.tsv", work / "negative.tsv"
        domain_pairs_out, domain_summary_out = work / "domain_pairs.tsv", work / "domain_summary.tsv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/19_build_extension_views.py"),
                "--rg", str(rg_out), "--positive-out", str(positive_out),
                "--negative-out", str(negative_out), "--domain-pairs-out", str(domain_pairs_out),
                "--domain-summary-out", str(domain_summary_out),
                "--provenance-out", str(work / "views.json"),
            ], check=True,
        )
        if len(read_tsv(positive_out)) != 1175 or len(read_tsv(negative_out)) != 1:
            raise SystemExit("ERROR: synthetic ranked positive/negative families are incomplete")
        expected_domains = len({row["phenotype_domain"] for row in panel[2:]})
        if len(read_tsv(domain_pairs_out)) != 1176 or len(read_tsv(domain_summary_out)) != expected_domains:
            raise SystemExit("ERROR: synthetic domain views are incomplete")
        figure = work / "screen.png"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "discovery_extension/scripts/15_plot_extension.py"),
                "--rg", str(rg_out), "--out", str(figure),
            ], check=True,
        )
        network = work / "extension_discovery_network.png"
        domain_views = work / "extension_domain_views.pdf"
        if not figure.is_file() or figure.stat().st_size < 10000 or not figure.with_suffix(".pdf").is_file():
            raise SystemExit("ERROR: synthetic discovery figure was not created")
        if not network.is_file() or network.stat().st_size < 10000 or not domain_views.is_file() or domain_views.stat().st_size < 10000:
            raise SystemExit("ERROR: synthetic network/domain figures were not created")
        audit = work / "novelty.tsv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/17_prepare_pair_novelty_audit.py"),
                "--rg", str(rg_out), "--out", str(audit),
                "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
                "--provenance-out", str(work / "novelty.json"),
            ], check=True,
        )
        audit_rows = read_tsv(audit)
        if len(audit_rows) != 1 or audit_rows[0]["novelty_class"] != "PENDING":
            raise SystemExit(f"ERROR: unexpected synthetic novelty template: {audit_rows}")
        incomplete = subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/18_validate_pair_novelty_audit.py"),
                "--rg", str(rg_out), "--audit", str(audit),
            ], check=False, capture_output=True, text=True,
        )
        if incomplete.returncode == 0:
            raise SystemExit("ERROR: incomplete novelty audit incorrectly passed")
        completed = audit_rows[0]
        completed.update({
            "audit_status": "COMPLETE", "direct_prior_same_pair": "NO",
            "same_sleep_trait_context": "NO", "same_or_equivalent_phenotype": "NO",
            "same_direction": "NOT_APPLICABLE", "broad_phenome_screen_overlap": "NO",
            "near_neighbor_evidence": "NO", "discovery_vs_replication_in_prior_work": "NONE",
            "exact_prior_rg_found": "NO", "closest_prior_result": "NONE_FOUND_AFTER_TARGETED_SEARCHES",
            "prior_method": "NONE_FOUND", "prior_effect": "NONE_FOUND",
            "prior_publication": "NONE_FOUND", "prior_DOI": "NONE_FOUND", "prior_PMID": "NONE_FOUND",
            "search_databases": "PubMed;GWAS Catalog",
            "search_queries_used": "query one || query two || query three",
            "search_date": "2099-01-01", "evidence_PMIDs_DOIs_URLs": "PMID:00000000",
            "biological_plausibility": "MODERATE", "connection_obviousness": "NON_OBVIOUS",
            "independent_replication_dataset_availability": "AVAILABLE",
            "dense_summary_statistics_available": "YES", "molecular_qtl_data_available": "YES",
            "independent_replication_class": "NOT_YET_ATTEMPTED",
            "novelty_class": "APPARENTLY_NOVEL", "novelty_strength": "MODERATE",
            "decision_rationale": "synthetic validator exercise", "reviewer_notes": "synthetic only",
            "reviewer": "synthetic_test",
        })
        with audit.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(completed), lineterminator="\n")
            writer.writeheader()
            writer.writerow(completed)
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/18_validate_pair_novelty_audit.py"),
                "--rg", str(rg_out), "--audit", str(audit),
            ], check=True,
        )
        priorities = work / "priorities.tsv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/20_prioritize_extension_pairs.py"),
                "--rg", str(rg_out), "--extension-h2", str(h2_out),
                "--core-h2", str(core_h2),
                "--novelty-audit", str(audit), "--out", str(priorities),
                "--provenance-out", str(work / "priorities.json"),
            ], check=True,
        )
        priority_rows = read_tsv(priorities)
        if len(priority_rows) != 1176 or sum(row["priority_tier"] == "A" for row in priority_rows) != 1:
            raise SystemExit("ERROR: synthetic pre-replication A/B/C prioritization family is incorrect")

        # The sleep-side eligibility gate must matter: removing it should
        # demote the sole otherwise-eligible discovery while other gates stay fixed.
        failed_core_rows = [dict(row) for row in core_h2_rows]
        failed_core_rows[0].update({"h2": 0.02, "se": 0.01, "z": 2, "verdict": "FAIL"})
        write_core_h2(failed_core_rows)
        failed_priorities = work / "synthetic_sleep_qc_failed_priorities.tsv"
        subprocess.run(
            [sys.executable, str(ROOT / "discovery_extension/scripts/20_prioritize_extension_pairs.py"),
             "--rg", str(rg_out), "--extension-h2", str(h2_out), "--core-h2", str(core_h2),
             "--novelty-audit", str(audit), "--out", str(failed_priorities),
             "--provenance-out", str(work / "synthetic_sleep_qc_failed_priorities.json")], check=True,
        )
        demoted = next(row for row in read_tsv(failed_priorities) if row["pair_id"] == completed["pair_id"])
        if demoted["priority_tier"] != "C" or "sleep_h2_failed" not in demoted["priority_rationale"]:
            raise SystemExit("ERROR: failed synthetic sleep h2 gate did not demote the eligible discovery")
        write_core_h2(core_h2_rows)

        local_queue = work / "local_queue.tsv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/25_prepare_local_analysis_queue.py"),
                "--priority", str(priorities), "--out", str(local_queue),
                "--provenance-out", str(work / "local_queue.json"),
            ], check=True,
        )
        local_queue_rows = read_tsv(local_queue)
        if len(local_queue_rows) != 1176 or sum(
            row["selection_status"] == "SELECTED_PRIORITY_DISCOVERY" for row in local_queue_rows
        ) != 1 or sum(
            row["selection_status"] == "PENDING_GLOBAL_NULL_CURATION" for row in local_queue_rows
        ) != 1175:
            raise SystemExit("ERROR: synthetic priority/globally-null local queue is incorrect")

        replication_queue = work / "replication_queue.tsv"
        replication_candidate_lock = work / "replication_candidate_family.lock.json"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/21_prepare_replication_queue.py"),
                "--priority", str(priorities), "--out", str(replication_queue),
                "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
                "--core", str(ROOT / "config/analysis_panel.tsv"),
                "--provenance-out", str(work / "replication_queue.json"),
                "--candidate-lock-out", str(replication_candidate_lock),
            ], check=True,
        )
        queue_rows = read_tsv(replication_queue)
        if len(queue_rows) != 1 or queue_rows[0]["source_curation_status"] != "PENDING_INDEPENDENT_SOURCE_CURATION":
            raise SystemExit("ERROR: synthetic replication queue is incorrect")
        curated = queue_rows[0]
        replication_source = work / "synthetic_replication_sumstats.tsv.gz"
        replication_source.write_bytes(b"SYNTHETIC REPLICATION SOURCE - NOT REAL RESULTS\n")
        replication_source_sha = hashlib.sha256(replication_source.read_bytes()).hexdigest()
        curated.update({
            "replication_source_id": "synthetic_non_ukb", "replication_study_accession": "SYNTHETIC001",
            "replication_publication": "synthetic test only", "replication_PMID": "00000000",
            "replication_DOI": "10.0000/synthetic", "replication_source_url": "https://example.org/synthetic.gz",
            "replication_checksum": f"sha256:{replication_source_sha}",
            "replication_source_generation": "NA_LOCAL_SYNTHETIC",
            "replication_etag": "NA_LOCAL_SYNTHETIC",
            "replication_content_length_bytes": str(replication_source.stat().st_size),
            "replication_storage_mode": "LOCAL_FULL_RESOLUTION",
            "replication_local_path": str(replication_source),
            "replication_receipt_path": "NA_LOCAL_SYNTHETIC",
            "replication_munged_path": str(replication_source),
            "replication_phenotype_definition": "synthetic exact phenotype",
            "phenotype_match_status": "EXACT", "ancestry": "EUR", "build": "GRCh37",
            "sample_size": "50000", "cases": "NA", "controls": "NA",
            "discovery_cohort_relation": "NON_UKB", "participant_overlap_status": "NON_OVERLAPPING_CONFIRMED",
            "participant_overlap_evidence": "synthetic disjoint cohorts", "source_identity_status": "VERIFIED",
            "schema_status": "VERIFIED", "effect_allele_status": "UNAMBIGUOUS",
            "full_resolution_availability": "YES", "source_curation_status": "COMPLETE_BEFORE_RESULTS",
            "replication_search_databases": "GWAS Catalog;PubMed;consortium repository",
            "replication_search_queries": "synthetic exact phenotype GWAS || synthetic consortium summary statistics",
            "replication_search_date": "2099-01-01", "replication_search_evidence": "synthetic test evidence only",
            "unavailable_reason": "NA",
        })
        with replication_queue.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(curated), lineterminator="\n")
            writer.writeheader()
            writer.writerow(curated)
        replication_manifest, replication_lock = work / "replication_manifest.tsv", work / "replication_manifest.lock.json"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/22_lock_replication_manifest.py"),
                "--queue", str(replication_queue), "--out", str(replication_manifest),
                "--candidate-lock", str(replication_candidate_lock), "--lock", str(replication_lock),
            ], check=True,
        )
        if not replication_manifest.is_file() or not replication_lock.is_file():
            raise SystemExit("ERROR: synthetic replication manifest did not lock")

        replication_h2 = work / "replication_h2.tsv"
        with replication_h2.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, delimiter="\t",
                fieldnames=["replication_source_id", "h2", "h2_se", "h2_z", "LDSC_intercept", "primary_status"],
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow({
                "replication_source_id": "synthetic_non_ukb", "h2": 0.2, "h2_se": 0.02,
                "h2_z": 10, "LDSC_intercept": 1.0, "primary_status": "PASS",
            })
        replication_rg = work / "replication_rg.tsv"
        with replication_rg.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, delimiter="\t",
                fieldnames=[
                    "pair_id", "sleep_trait", "replication_source_id", "rg", "se", "z", "p",
                    "cross_trait_LDSC_intercept", "cross_trait_LDSC_intercept_se",
                    "snp_overlap_valid_alleles", "ancestry", "analysis_status",
                ], lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow({
                "pair_id": curated["pair_id"], "sleep_trait": curated["sleep_trait"],
                "replication_source_id": "synthetic_non_ukb", "rg": 0.18, "se": 0.04,
                "z": 4.5, "p": 0.01, "cross_trait_LDSC_intercept": 0.0,
                "cross_trait_LDSC_intercept_se": 0.005, "snp_overlap_valid_alleles": 900000,
                "ancestry": "EUR", "analysis_status": "REPLICATION_RG_COMPLETE",
            })
        replication_results = work / "replication_results.tsv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/23_collate_replication.py"),
                "--manifest", str(replication_manifest), "--lock", str(replication_lock),
                "--h2", str(replication_h2), "--rg", str(replication_rg),
                "--out", str(replication_results), "--provenance-out", str(work / "replication_results.json"),
            ], check=True,
        )
        if read_tsv(replication_results)[0]["replication_class"] != "REPLICATED":
            raise SystemExit("ERROR: synthetic independent replication was not classified REPLICATED")
        replicated_priorities = work / "replicated_priorities.tsv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "discovery_extension/scripts/20_prioritize_extension_pairs.py"),
                "--rg", str(rg_out), "--extension-h2", str(h2_out),
                "--core-h2", str(core_h2),
                "--novelty-audit", str(audit), "--replication", str(replication_results),
                "--out", str(replicated_priorities), "--provenance-out", str(work / "replicated_priorities.json"),
            ], check=True,
        )
        if sum(row["priority_tier"] == "B" for row in read_tsv(replicated_priorities)) != 1:
            raise SystemExit("ERROR: synthetic replicated pair did not advance to Tier B")
    print("EXTENSION_LDSC_COLLATION_SYNTHETIC_OK h2=100 primary_pairs=1176 universe=1200 ranked_views=true figures=true novelty_gate=true prioritization=true replication_lock=true local_queue=true sleep_qc_fixture=true sleep_qc_failure_gate=true isolated=true")


if __name__ == "__main__":
    main()
