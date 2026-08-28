#!/usr/bin/env python3
"""Exercise full 100-trait h2 and 12x98 rg collation in isolation."""

from __future__ import annotations

import csv
import subprocess
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
                "python3", str(ROOT / "discovery_extension/scripts/13_collate_extension_ldsc.py"),
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
                "python3", str(ROOT / "discovery_extension/scripts/13_collate_extension_ldsc.py"),
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
                "python3", str(ROOT / "discovery_extension/scripts/19_build_extension_views.py"),
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
                str(ROOT / ".venv/bin/python"),
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
                "python3", str(ROOT / "discovery_extension/scripts/17_prepare_pair_novelty_audit.py"),
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
                "python3", str(ROOT / "discovery_extension/scripts/18_validate_pair_novelty_audit.py"),
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
            "independent_replication_status": "INDEPENDENT_REPLICATION_PASS",
            "novelty_class": "APPARENTLY_NOVEL", "novelty_strength": "STRONG",
            "decision_rationale": "synthetic validator exercise", "reviewer_notes": "synthetic only",
            "reviewer": "synthetic_test",
        })
        with audit.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(completed), lineterminator="\n")
            writer.writeheader()
            writer.writerow(completed)
        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/18_validate_pair_novelty_audit.py"),
                "--rg", str(rg_out), "--audit", str(audit),
            ], check=True,
        )
        priorities = work / "priorities.tsv"
        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/20_prioritize_extension_pairs.py"),
                "--rg", str(rg_out), "--extension-h2", str(h2_out),
                "--core-h2", str(ROOT / "results/tables/h2_summary.tsv"),
                "--novelty-audit", str(audit), "--out", str(priorities),
                "--provenance-out", str(work / "priorities.json"),
            ], check=True,
        )
        priority_rows = read_tsv(priorities)
        if len(priority_rows) != 1176 or sum(row["priority_tier"] == "B" for row in priority_rows) != 1:
            raise SystemExit("ERROR: synthetic A/B/C prioritization family is incorrect")
    print("EXTENSION_LDSC_COLLATION_SYNTHETIC_OK h2=100 primary_pairs=1176 universe=1200 ranked_views=true figures=true novelty_gate=true prioritization=true isolated=true")


if __name__ == "__main__":
    main()
