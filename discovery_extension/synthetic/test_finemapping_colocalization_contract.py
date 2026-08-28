#!/usr/bin/env python3
"""Exercise result-free selection, locking, SuSiE-RSS, coloc, and unavailable-QTL retention."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "discovery_extension"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def run(command: list[str]) -> str:
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    if result.returncode:
        raise AssertionError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}{result.stderr}")
    return result.stdout + result.stderr


class FineMappingColocalizationContractTest(unittest.TestCase):
    def test_end_to_end_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="sleep_gwas_finemap_synthetic_") as temporary:
            temp = Path(temporary)
            loci = temp / "novel_shared_loci.tsv"
            priority = temp / "priority.tsv"
            replication = temp / "replication.tsv"
            queue = temp / "queue.tsv"
            family_lock = temp / "candidate.lock.json"
            queue_provenance = temp / "queue.json"
            locus_row = {
                "pair_id": "sleep_test__external_test", "sleep_trait": "sleep_test",
                "extension_trait_id": "external_test", "external_phenotype_name": "Synthetic external trait",
                "locus_id": "locus_1", "CHR": 1, "start": 1, "end": 1000000,
                "lead_snp": "rs61", "lead_P_PLACO_PLUS": "1e-12",
            }
            write_tsv(loci, list(locus_row), [locus_row])
            priority_row = {
                "pair_id": locus_row["pair_id"], "priority_tier": "B", "rg": "0.31", "extension_fdr": "0.001",
            }
            write_tsv(priority, list(priority_row), [priority_row])
            replication_row = {"pair_id": locus_row["pair_id"], "replication_class": "REPLICATED"}
            write_tsv(replication, list(replication_row), [replication_row])
            output = run([
                "python3", str(EXTENSION / "scripts/33_prepare_finemapping_queue.py"),
                "--loci", str(loci), "--priority", str(priority), "--replication", str(replication),
                "--contract", str(EXTENSION / "config/fine_mapping_colocalization_contract.json"),
                "--out", str(queue), "--candidate-lock-out", str(family_lock),
                "--provenance-out", str(queue_provenance),
            ])
            self.assertIn("FINEMAPPING_QUEUE_OK loci=1 rows=4", output)

            variant_count = 120
            snps = [f"rs{index + 1}" for index in range(variant_count)]
            positions = [1000 + 100 * index for index in range(variant_count)]
            rho = 0.75
            indices = np.arange(variant_count)
            ld = rho ** np.abs(indices[:, None] - indices[None, :])
            ld_path = temp / "ld.tsv"
            with ld_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["SNP", *snps])
                for snp, values in zip(snps, ld, strict=True):
                    writer.writerow([snp, *(f"{value:.12g}" for value in values)])
            order_path = temp / "variant_order.tsv"
            write_tsv(order_path, ["SNP"], [{"SNP": snp} for snp in snps])

            def summary(path: Path, amplitude: float) -> None:
                z = amplitude * ld[:, 60]
                rows = []
                for index, (snp, position) in enumerate(zip(snps, positions, strict=True)):
                    rows.append({
                        "SNP": snp, "CHR": 1, "BP": position,
                        "A1": "A" if index % 2 == 0 else "C", "A2": "G" if index % 2 == 0 else "T",
                        "BETA": f"{z[index] * 0.01:.12g}", "SE": "0.01", "MAF": "0.25",
                        "INFO": "0.99", "PRIOR_WEIGHT": "1",
                    })
                write_tsv(path, list(rows[0]), rows)

            sleep_summary = temp / "sleep.tsv"
            external_summary = temp / "external.tsv"
            eqtl_summary = temp / "eqtl.tsv"
            summary(sleep_summary, 8.0)
            summary(external_summary, 7.5)
            summary(eqtl_summary, 9.0)
            eqtl_log = temp / "eqtl_search.tsv"
            sqtl_log = temp / "sqtl_search.tsv"
            pqtl_log = temp / "pqtl_search.tsv"
            for path, modality, outcome in (
                (eqtl_log, "EQTL", "SYNTHETIC_COMPATIBLE_DATASET"),
                (sqtl_log, "SQTL", "NO_SUITABLE_DATASET"),
                (pqtl_log, "PQTL", "NO_SUITABLE_DATASET"),
            ):
                write_tsv(path, ["modality", "outcome"], [{"modality": modality, "outcome": outcome}])

            fields, rows = read_tsv(queue)
            for row in rows:
                modality = row["qtl_modality"]
                row.update({
                    "dataset1_type": "quant", "dataset1_N": "100000", "dataset1_effective_N": "100000",
                    "dataset1_case_fraction": "NA", "dataset1_sdY": "1", "dataset1_prior_method": "FLAT",
                    "dataset1_prior_source_path": "NA", "dataset1_prior_source_sha256": "NA",
                    "dataset2_type": "quant", "dataset2_N": "90000", "dataset2_effective_N": "90000",
                    "dataset2_case_fraction": "NA", "dataset2_sdY": "1", "dataset2_prior_method": "FLAT",
                    "dataset2_prior_source_path": "NA", "dataset2_prior_source_sha256": "NA",
                    "source_id": "SYNTHETIC_FIXTURE", "source_url": "synthetic://fine-mapping-fixture",
                    "source_access_date": "2026-08-28", "source_license": "SYNTHETIC_ONLY",
                    "phenotype_compatibility": "VERIFIED", "sample_overlap_status": "OVERLAP_DISCLOSED",
                    "summary1_path": str(sleep_summary), "summary1_sha256": sha256(sleep_summary),
                    "ld_reference_id": "SYNTHETIC_AR1", "ld_path": str(ld_path), "ld_sha256": sha256(ld_path),
                    "ld_variant_order_path": str(order_path), "ld_variant_order_sha256": sha256(order_path),
                    "ld_source_type": "IN_SAMPLE", "ld_ancestry": "EUR", "ld_build": "GRCh37",
                    "ld_sample_size": "20000", "schema_status": "VERIFIED",
                    "effect_allele_alignment_status": "VERIFIED", "dense_variant_status": "FULL_RESOLUTION",
                    "single_signal_fallback_justification": "NOT_JUSTIFIED", "results_accessed_before_lock": "NO",
                    "curator": "synthetic_test", "curation_date": "2026-08-28", "notes": "SYNTHETIC_TEST_ONLY",
                })
                if modality == "NONE":
                    row.update({
                        "source_search_status": "VERIFIED_ANALYSIS", "summary2_path": str(external_summary),
                        "summary2_sha256": sha256(external_summary), "molecular_feature_id": "NA",
                        "molecular_feature_name": "NA", "tissue_cell_context": "NA",
                        "source_search_log_path": "NA", "source_search_log_sha256": "NA",
                    })
                elif modality == "EQTL":
                    row.update({
                        "queue_row_id": f"{locus_row['pair_id']}__locus_1__SLEEP_EQTL__GENE1",
                        "comparison_type": "SLEEP_EQTL", "comparison_role": "MOLECULAR_QTL_FOLLOWUP",
                        "dataset2_id": "GENE1_eqtl", "dataset2_role": "MOLECULAR_QTL", "dataset2_N": "1000",
                        "dataset2_effective_N": "1000", "molecular_feature_id": "GENE1",
                        "molecular_feature_name": "Synthetic gene 1", "tissue_cell_context": "SYNTHETIC_TISSUE",
                        "source_search_status": "VERIFIED_ANALYSIS", "source_search_log_path": str(eqtl_log),
                        "source_search_log_sha256": sha256(eqtl_log), "summary2_path": str(eqtl_summary),
                        "summary2_sha256": sha256(eqtl_summary), "sample_overlap_status": "NON_OVERLAPPING",
                    })
                else:
                    log = sqtl_log if modality == "SQTL" else pqtl_log
                    row.update({
                        "source_search_status": "NO_SUITABLE_DATASET", "source_search_log_path": str(log),
                        "source_search_log_sha256": sha256(log), "source_id": f"SYNTHETIC_{modality}_SEARCH",
                        "source_url": f"synthetic://{modality.lower()}-search", "notes": f"SYNTHETIC_{modality}_UNAVAILABLE",
                    })
            write_tsv(queue, fields, rows)

            manifest = temp / "manifest.tsv"
            manifest_lock = temp / "manifest.lock.json"
            lock_output = run([
                "python3", str(EXTENSION / "scripts/34_lock_finemapping_manifest.py"),
                "--queue", str(queue), "--candidate-lock", str(family_lock),
                "--readiness", str(EXTENSION / "results/fine_mapping/fine_mapping_readiness.tsv"),
                "--contract", str(EXTENSION / "config/fine_mapping_colocalization_contract.json"),
                "--out", str(manifest), "--lock", str(manifest_lock),
            ])
            self.assertIn("analyses=2 unavailable=2", lock_output)
            lock = json.loads(manifest_lock.read_text(encoding="utf-8"))

            raw_paths: dict[str, list[Path]] = {key: [] for key in ("variant", "credible", "coloc", "shared", "diagnostic")}
            for index, comparison_id in enumerate(lock["analysis_comparison_ids_in_locked_order"]):
                paths = {key: temp / f"raw_{index}_{key}.tsv" for key in raw_paths}
                for key, path in paths.items():
                    raw_paths[key].append(path)
                run_output = run([
                    str(ROOT / ".r-env/bin/Rscript"), str(EXTENSION / "scripts/35_run_susie_coloc.R"),
                    str(manifest), comparison_id, str(paths["variant"]), str(paths["credible"]),
                    str(paths["coloc"]), str(paths["shared"]), str(paths["diagnostic"]),
                    str(lock["susie_parameters"]["L"]), str(lock["susie_parameters"]["coverage"]),
                    str(lock["susie_parameters"]["min_abs_corr"]), str(lock["susie_parameters"]["maxit"]),
                    str(lock["coloc_primary_priors"]["p1"]), str(lock["coloc_primary_priors"]["p2"]),
                    ",".join(str(value) for value in lock["coloc_p12_sensitivity_grid"]),
                ])
                self.assertIn("SUSIE_COLOC_COMPARISON_OK", run_output)

            variant_out = temp / "fine_mapping_variants.tsv"
            credible_out = temp / "credible_sets.tsv"
            shared_out = temp / "shared_variants.tsv"
            coloc_out = temp / "fine_mapping_colocalization.tsv"
            provenance = temp / "results.json"
            collate_output = run([
                "python3", str(EXTENSION / "scripts/36_collate_finemapping_coloc.py"),
                "--manifest", str(manifest), "--lock", str(manifest_lock),
                "--variant-results", *(str(path) for path in raw_paths["variant"]),
                "--credible-set-results", *(str(path) for path in raw_paths["credible"]),
                "--coloc-results", *(str(path) for path in raw_paths["coloc"]),
                "--shared-variant-results", *(str(path) for path in raw_paths["shared"]),
                "--diagnostic-results", *(str(path) for path in raw_paths["diagnostic"]),
                "--variant-out", str(variant_out), "--credible-set-out", str(credible_out),
                "--shared-variant-out", str(shared_out), "--out", str(coloc_out),
                "--provenance-out", str(provenance),
            ])
            self.assertIn("FINEMAPPING_COLOC_RESULTS_OK analyses=2 unavailable=2", collate_output)
            _, variant_rows = read_tsv(variant_out)
            _, credible_rows = read_tsv(credible_out)
            _, coloc_rows = read_tsv(coloc_out)
            _, shared_rows = read_tsv(shared_out)
            self.assertEqual(len(variant_rows), 2 * 2 * variant_count)
            self.assertGreaterEqual(len(credible_rows), 4)
            self.assertEqual(sum(row["analysis_status"] == "NO_SUITABLE_DATASET" for row in coloc_rows), 2)
            primary = [row for row in coloc_rows if row["prior_role"] == "PRIMARY"]
            self.assertEqual(len(primary), 2)
            self.assertTrue(all(row["colocalization_interpretation"] == "SHARED_SIGNAL_MODEL_SUPPORTED" for row in primary))
            self.assertTrue(all(row["prior_robust"] == "True" for row in primary))
            self.assertGreater(len(shared_rows), 0)
            self.assertTrue(all(row["claim_limit"].startswith("SHARED_SIGNAL_MODEL_NOT_") for row in primary))

            bad_ld = temp / "bad_asymmetric_ld.tsv"
            asymmetric = ld.copy()
            asymmetric[0, 1] = 0.1
            with bad_ld.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["SNP", *snps])
                for snp, values in zip(snps, asymmetric, strict=True):
                    writer.writerow([snp, *(f"{value:.12g}" for value in values)])
            bad_queue = temp / "bad_queue.tsv"
            for row in rows:
                if row["source_search_status"] == "VERIFIED_ANALYSIS":
                    row["ld_path"] = str(bad_ld)
                    row["ld_sha256"] = sha256(bad_ld)
            write_tsv(bad_queue, fields, rows)
            rejected = subprocess.run([
                "python3", str(EXTENSION / "scripts/34_lock_finemapping_manifest.py"),
                "--queue", str(bad_queue), "--candidate-lock", str(family_lock),
                "--readiness", str(EXTENSION / "results/fine_mapping/fine_mapping_readiness.tsv"),
                "--contract", str(EXTENSION / "config/fine_mapping_colocalization_contract.json"),
                "--out", str(temp / "bad_manifest.tsv"), "--lock", str(temp / "bad_manifest.lock.json"),
            ], cwd=ROOT, text=True, capture_output=True, check=False)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("LD matrix fails symmetry/unit-diagonal contract", rejected.stdout + rejected.stderr)


if __name__ == "__main__":
    unittest.main()
