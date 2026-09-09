#!/usr/bin/env python3
"""Exercise the result-free local lock and canonical local-result flags."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="synthetic_local_architecture_") as temporary:
        work = Path(temporary)
        trait_a, trait_b1, trait_b2, loci = (work / name for name in ("sleep.tsv", "trait1.tsv", "trait2.tsv", "loci.tsv"))
        trait_a.write_text("SYNTHETIC_NOT_REAL\n", encoding="utf-8")
        trait_b1.write_text("SYNTHETIC_NOT_REAL_TRAIT_1\n", encoding="utf-8")
        trait_b2.write_text("SYNTHETIC_NOT_REAL_TRAIT_2\n", encoding="utf-8")
        loci.write_text("locus_id\tchr\tstart\tend\nL1\t1\t1\t100\nL2\t1\t101\t200\n", encoding="utf-8")

        readiness = work / "readiness.tsv"
        write_tsv(readiness, ["method_id", "readiness"], [
            {"method_id": "LAVA_PRIMARY", "readiness": "READY"},
            {"method_id": "HDL_L_ROBUSTNESS", "readiness": "READY"},
        ])
        queue = work / "queue.tsv"
        common = {
            "sleep_trait": "SYNTHETIC_SLEEP", "phenotype_domain": "Synthetic",
            "global_rg_se": "0.03", "global_rg_p": "0.001",
            "global_analysis_status": "PRIMARY_EXTENSION_RG_COMPLETE",
            "replication_class": "REPLICATED", "pair_novelty_class": "APPARENTLY_NOVEL",
            "pair_novelty_strength": "MODERATE", "snp_overlap_valid_alleles": "1000000",
            "biological_relevance": "HIGH", "polygenic_overlap_evidence": "synthetic documented overlap",
            "dense_summary_statistics_available": "YES", "trait_a_dense_path": str(trait_a),
            "trait_a_dense_sha256": sha256(trait_a), "locus_definition_path": str(loci),
            "locus_definition_sha256": sha256(loci), "locus_count": "2",
            "hdl_feasibility_rationale": "synthetic reference and inputs are compatible",
            "lava_simulation_random_seed": "1234567",
            "result_access_status": "NOT_ACCESSED", "curator": "synthetic_test", "curation_date": "2099-01-01",
        }
        queue_rows = [
            {
                **common, "pair_id": "SYNTHETIC_SLEEP__TRAIT1", "extension_trait_id": "TRAIT1",
                "phenotype_name": "Synthetic trait 1", "global_rg": "0.20", "global_rg_fdr": "0.01",
                "discovery_priority_tier": "B", "selection_stratum": "PRIORITY_DISCOVERY",
                "selection_status": "SELECTED_PRIORITY_DISCOVERY",
                "secondary_selection_rationale": "synthetic priority pair",
                "trait_b_dense_path": str(trait_b1), "trait_b_dense_sha256": sha256(trait_b1),
                "planned_methods": "LAVA_PRIMARY;HDL_L_ROBUSTNESS",
            },
            {
                **common, "pair_id": "SYNTHETIC_SLEEP__TRAIT2", "extension_trait_id": "TRAIT2",
                "phenotype_name": "Synthetic trait 2", "global_rg": "0.01", "global_rg_fdr": "0.50",
                "discovery_priority_tier": "C", "selection_stratum": "GLOBAL_NULL_SECONDARY",
                "selection_status": "SELECTED_GLOBAL_NULL_SECONDARY",
                "secondary_selection_rationale": "synthetic biologically relevant globally-null pair",
                "trait_b_dense_path": str(trait_b2), "trait_b_dense_sha256": sha256(trait_b2),
                "planned_methods": "LAVA_PRIMARY",
            },
        ]
        queue_fields = [
            "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
            "global_rg", "global_rg_se", "global_rg_p", "global_rg_fdr", "global_analysis_status",
            "discovery_priority_tier", "replication_class", "pair_novelty_class", "pair_novelty_strength",
            "snp_overlap_valid_alleles", "selection_stratum", "selection_status",
            "biological_relevance", "polygenic_overlap_evidence", "secondary_selection_rationale",
            "dense_summary_statistics_available", "trait_a_dense_path", "trait_a_dense_sha256",
            "trait_b_dense_path", "trait_b_dense_sha256", "locus_definition_path",
            "locus_definition_sha256", "locus_count", "planned_methods", "hdl_feasibility_rationale",
            "lava_simulation_random_seed", "result_access_status", "curator", "curation_date",
        ]
        write_tsv(queue, queue_fields, queue_rows)
        manifest, lock = work / "manifest.tsv", work / "manifest.lock.json"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/26_lock_local_analysis_manifest.py"),
            "--queue", str(queue), "--readiness", str(readiness), "--out", str(manifest), "--lock", str(lock),
        ], check=True)
        locked = json.loads(lock.read_text(encoding="utf-8"))
        if locked["pair_count"] != 2 or locked["lava_pair_locus_family_size"] != 4 or locked["results_accessed_before_lock"]:
            raise SystemExit("ERROR: synthetic local family lock is incorrect")

        raw = work / "normalized_local.tsv"
        raw_fields = [
            "pair_id", "method_id", "locus_id", "chr", "start", "end",
            "local_h2_trait_a", "local_h2_trait_a_se", "local_h2_trait_a_p",
            "local_h2_trait_b", "local_h2_trait_b_se", "local_h2_trait_b_p",
            "local_rg", "local_rg_se", "local_rg_p", "bivariate_test_status",
            "effect_direction", "convergence_status", "qc_status",
        ]
        raw_rows = []
        specifications = [
            ("SYNTHETIC_SLEEP__TRAIT1", "LAVA_PRIMARY", "L1", 0.30, 1e-6),
            ("SYNTHETIC_SLEEP__TRAIT1", "LAVA_PRIMARY", "L2", -0.25, 1e-5),
            ("SYNTHETIC_SLEEP__TRAIT1", "HDL_L_ROBUSTNESS", "L1", 0.28, 0.01),
            ("SYNTHETIC_SLEEP__TRAIT1", "HDL_L_ROBUSTNESS", "L2", 0.01, 0.80),
            ("SYNTHETIC_SLEEP__TRAIT2", "LAVA_PRIMARY", "L1", 0.22, 2e-6),
            ("SYNTHETIC_SLEEP__TRAIT2", "LAVA_PRIMARY", "L2", 0.01, 0.50),
        ]
        for pair_id, method, locus_id, rg, p_value in specifications:
            raw_rows.append({
                "pair_id": pair_id, "method_id": method, "locus_id": locus_id, "chr": "1",
                "start": "1" if locus_id == "L1" else "101", "end": "100" if locus_id == "L1" else "200",
                "local_h2_trait_a": "0.02", "local_h2_trait_a_se": "0.003", "local_h2_trait_a_p": "0.001",
                "local_h2_trait_b": "0.03", "local_h2_trait_b_se": "0.004", "local_h2_trait_b_p": "0.001",
                "local_rg": str(rg), "local_rg_se": "0.05", "local_rg_p": str(p_value),
                "bivariate_test_status": "TESTED", "effect_direction": "POSITIVE" if rg > 0 else "NEGATIVE",
                "convergence_status": "PASS", "qc_status": "PASS",
            })
        write_tsv(raw, raw_fields, raw_rows)
        results, summary = work / "local.tsv", work / "summary.tsv"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/27_collate_local_architecture.py"),
            "--manifest", str(manifest), "--lock", str(lock), "--raw", str(raw),
            "--out", str(results), "--pair-summary-out", str(summary),
            "--provenance-out", str(work / "local.json"),
        ], check=True)
        summaries = {row["pair_id"]: row for row in read_tsv(summary)}
        first = summaries["SYNTHETIC_SLEEP__TRAIT1"]
        second = summaries["SYNTHETIC_SLEEP__TRAIT2"]
        if first["local_support_status"] != "STRONG_LOCAL_SUPPORT" or set(first["architecture_flags"].split(";")) != {"OPPOSING_LOCAL_EFFECTS", "MULTI_LOCUS_SHARED_ARCHITECTURE"}:
            raise SystemExit(f"ERROR: synthetic opposing/multi-locus flags failed: {first}")
        if second["architecture_flags"] != "GLOBAL_NULL_LOCAL_POSITIVE" or second["local_support_status"] != "LOCAL_SUPPORT_LAVA_ONLY":
            raise SystemExit(f"ERROR: synthetic globally-null local flag failed: {second}")
        if len(read_tsv(results)) != 6:
            raise SystemExit("ERROR: synthetic local result family did not retain every method/locus row")
    print("LOCAL_ARCHITECTURE_SYNTHETIC_OK pairs=2 rows=6 global_null=true opposing=true multi_locus=true robustness=true isolated=true")


if __name__ == "__main__":
    main()
