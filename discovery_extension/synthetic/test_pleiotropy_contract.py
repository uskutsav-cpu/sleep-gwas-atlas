#!/usr/bin/env python3
"""Exercise PLACO+ pair/input locking and canonical shared-locus collation."""

from __future__ import annotations

import csv
import hashlib
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    pair_id = "SYNTHETIC_SLEEP__SYNTHETIC_TRAIT"
    with tempfile.TemporaryDirectory(prefix="synthetic_pleiotropy_") as temporary:
        work = Path(temporary)
        priority = work / "priority.tsv"
        write_tsv(priority, [
            "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name", "rg", "extension_fdr", "priority_tier",
        ], [{
            "pair_id": pair_id, "sleep_trait": "SYNTHETIC_SLEEP", "extension_trait_id": "SYNTHETIC_TRAIT",
            "phenotype_name": "Synthetic trait", "rg": "0.2", "extension_fdr": "0.01", "priority_tier": "B",
        }])
        replication = work / "replication.tsv"
        write_tsv(replication, ["pair_id", "replication_rg", "replication_p", "replication_class"], [{
            "pair_id": pair_id, "replication_rg": "0.18", "replication_p": "0.01", "replication_class": "REPLICATED",
        }])
        queue, candidate_lock = work / "queue.tsv", work / "candidate.lock.json"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/29_prepare_pleiotropy_queue.py"),
            "--priority", str(priority), "--replication", str(replication),
            "--sources", str(ROOT / "discovery_extension/config/pleiotropy_sources.tsv"),
            "--out", str(queue), "--candidate-lock-out", str(candidate_lock),
            "--provenance-out", str(work / "queue.json"),
        ], check=True)
        rows = read_tsv(queue)
        merged, ld_reference = work / "merged.tsv", work / "ld_reference.bin"
        merged_rows = []
        synthetic_z2 = ("2.7", "-1.5", "0.5", "-2.2")
        for index in range(1, 5):
            merged_rows.append({
                "SNP": f"rs{index}", "CHR": "1", "BP": str(index * 100), "A1": "A", "A2": "G",
                "Z1": str(2 + index / 10), "Z2": synthetic_z2[index - 1], "P1": "0.001", "P2": "0.002",
                "MAF": "0.2", "INFO1": "0.99", "INFO2": "0.98",
            })
        write_tsv(merged, list(merged_rows[0]), merged_rows)
        ld_reference.write_bytes(b"SYNTHETIC EUR LD REFERENCE - NOT REAL DATA\n")
        placo_source = ROOT / ".r-env/share/placo/PLACO_v0.2.0.R"
        rows[0].update({
            "merged_genomewide_path": str(merged), "merged_genomewide_sha256": sha256(merged),
            "input_variant_count": "5", "eligible_variant_count": "4", "z2_excluded_count": "1",
            "input_scope": "FULL_GENOME", "schema_status": "VERIFIED",
            "effect_allele_alignment_status": "VERIFIED", "maf_info_filter_status": "MAF_GE_0.01_INFO_GE_0.9",
            "z_squared_filter_status": "PRIMARY_Z2_LE_80_WITH_SEPARATE_LARGE_EFFECT_LEDGER",
            "ld_reference_id": "SYNTHETIC_EUR_LD", "ld_reference_path": str(ld_reference),
            "ld_reference_sha256": sha256(ld_reference), "clumping_parameters": "r2<0.1 within 1Mb synthetic",
            "placo_source_path": str(placo_source), "placo_source_sha256": sha256(placo_source),
            "curator": "synthetic_test", "curation_date": "2099-01-01",
        })
        write_tsv(queue, list(rows[0]), rows)
        readiness = work / "readiness.tsv"
        write_tsv(readiness, ["method_id", "code_status", "readiness"], [{
            "method_id": "PLACO_PLUS_PRIMARY", "code_status": "PASS", "readiness": "READY_FOR_PAIR_MANIFEST_CURATION",
        }])
        manifest, lock = work / "manifest.tsv", work / "manifest.lock.json"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/30_lock_pleiotropy_manifest.py"),
            "--queue", str(queue), "--candidate-lock", str(candidate_lock), "--readiness", str(readiness),
            "--out", str(manifest), "--lock", str(lock), "--minimum-genomewide-variants", "3",
        ], check=True)
        runner_variants, runner_summary = work / "runner_variants.tsv", work / "runner_summary.tsv"
        subprocess.run([
            str(ROOT / ".r-env/bin/Rscript"), str(ROOT / "discovery_extension/scripts/31_run_placo_plus.R"),
            str(manifest), pair_id, str(runner_variants), str(runner_summary), "1e-8",
        ], check=True)
        if len(read_tsv(runner_variants)) != 4 or read_tsv(runner_summary)[0]["analysis_status"] != "PLACO_PLUS_COMPLETE":
            raise SystemExit("ERROR: pinned PLACO+ R runner smoke failed")
        run_summary = work / "run_summary.tsv"
        write_tsv(run_summary, [
            "pair_id", "input_variant_count", "eligible_variant_count", "z2_excluded_count",
            "VarZ1", "VarZ2", "CorZ", "marginal_p_threshold", "analysis_status",
        ], [{
            "pair_id": pair_id, "input_variant_count": "5", "eligible_variant_count": "4",
            "z2_excluded_count": "1", "VarZ1": "1.0", "VarZ2": "1.1", "CorZ": "0.2",
            "marginal_p_threshold": "0.0001", "analysis_status": "PLACO_PLUS_COMPLETE",
        }])
        variants = work / "variants.tsv"
        variant_rows = []
        for index, placo_p in enumerate(("1e-9", "2e-9", "0.1", "0.5"), start=1):
            variant_rows.append({
                "pair_id": pair_id, "SNP": f"rs{index}", "CHR": "1", "BP": str(index * 100),
                "A1": "A", "A2": "G", "Z1": str(2 + index / 10), "Z2": str(2.5 + index / 10),
                "P1": "0.001", "P2": "0.002", "T_PLACO_PLUS": str(5 + index / 10),
                "P_PLACO_PLUS": placo_p, "analysis_status": "PLACO_PLUS_COMPLETE",
            })
        write_tsv(variants, list(variant_rows[0]), variant_rows)
        clumps = work / "clumps.tsv"
        write_tsv(clumps, [
            "pair_id", "SNP", "locus_id", "lead_snp", "CHR", "start", "end", "r2_to_lead", "clumping_reference_id",
        ], [
            {"pair_id": pair_id, "SNP": "rs1", "locus_id": "L1", "lead_snp": "rs1", "CHR": "1", "start": "50", "end": "250", "r2_to_lead": "1", "clumping_reference_id": "SYNTHETIC_EUR_LD"},
            {"pair_id": pair_id, "SNP": "rs2", "locus_id": "L1", "lead_snp": "rs1", "CHR": "1", "start": "50", "end": "250", "r2_to_lead": "0.2", "clumping_reference_id": "SYNTHETIC_EUR_LD"},
        ])
        output = work / "novel_shared_loci.tsv"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/31_collate_pleiotropy.py"),
            "--manifest", str(manifest), "--lock", str(lock), "--run-summary", str(run_summary),
            "--variant-results", str(variants), "--clump-map", str(clumps), "--out", str(output),
            "--provenance-out", str(work / "results.json"),
        ], check=True)
        loci = read_tsv(output)
        if len(loci) != 1 or loci[0]["lead_snp"] != "rs1" or loci[0]["significant_variant_count"] != "2":
            raise SystemExit(f"ERROR: synthetic PLACO+ locus collation failed: {loci}")
        if loci[0]["claim_limit"] != "NOT_SHARED_CAUSAL_VARIANT_OR_CAUSAL_MECHANISM":
            raise SystemExit("ERROR: synthetic PLACO+ claim guard is missing")
    print("PLEIOTROPY_SYNTHETIC_OK pairs=1 significant_variants=2 loci=1 family_locked=true isolated=true")


if __name__ == "__main__":
    main()
