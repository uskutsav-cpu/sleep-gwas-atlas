import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/131_run_track_b_lava_sequential.py"
SPEC = importlib.util.spec_from_file_location("track_b_ram_aware_execution", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
RAM = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RAM)


class TrackBRAMAwareExecutionTests(unittest.TestCase):
    def test_input_only_selection_is_min_median_max_reference_snp_count(self) -> None:
        loci = [
            {"LOC": value, "CHR": 1, "START": value * 100, "STOP": value * 100 + 99}
            for value in range(1, 6)
        ]
        rows = RAM.benchmark_selection(loci, [10, 30, 20, 50, 40], "a" * 64)
        self.assertEqual([row["selection_role"] for row in rows], ["SMALL", "MEDIAN", "LARGE"])
        self.assertEqual([row["locus_index"] for row in rows], [1, 2, 4])
        self.assertEqual([row["reference_n_snps"] for row in rows], [10, 30, 50])
        self.assertTrue(all("INPUT_ONLY" in str(row["selection_basis"]) for row in rows))

    def test_report_and_benchmark_have_the_exact_required_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            benchmark = root / "RAM_BENCHMARK.tsv"
            report = root / "RAM_AWARE_EXECUTION_REPORT.md"
            with (
                patch.object(RAM, "BENCHMARK", benchmark),
                patch.object(RAM, "REPORT", report),
                patch.object(RAM, "physical_memory_bytes", return_value=8 * 1024**3),
            ):
                RAM.ensure_benchmark_file()
                header = benchmark.read_text(encoding="utf-8").splitlines()[0].split("\t")
                self.assertEqual(header, RAM.BENCHMARK_FIELDS)
                text = report.read_text(encoding="utf-8")
                self.assertIn("PENDING", text)
                self.assertIn("global nuisance/correlation fit is not decomposable", text)
                self.assertIn("never within an LD block", text)

    def test_unsealed_historical_rows_cannot_enable_a_ram_pass(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            benchmark = root / "RAM_BENCHMARK.tsv"
            report = root / "RAM_AWARE_EXECUTION_REPORT.md"
            rows = [
                {
                    "analysis": "LAVA_DISCOVERY", "pair": "NONE", "locus": index,
                    "chromosome": 1, "n_snps": n_snps, "peak_ram_gb": peak,
                    "runtime_sec": runtime, "exit_status": 0, "output_hash": str(index) * 64,
                }
                for index, n_snps, peak, runtime in (
                    (1, 10, 0.5, 2.0), (2, 100, 1.5, 3.0), (3, 1000, 2.5, 4.0)
                )
            ]
            benchmark.write_text(RAM.tsv_text(RAM.BENCHMARK_FIELDS, rows), encoding="utf-8")
            with (
                patch.object(RAM, "BENCHMARK", benchmark),
                patch.object(RAM, "REPORT", report),
                patch.object(RAM, "physical_memory_bytes", return_value=8 * 1024**3),
            ):
                RAM.build_report()
            text = report.read_text(encoding="utf-8")
            self.assertIn("PENDING — fingerprint-bound representative measurements are incomplete", text)
            self.assertNotIn("Current LAVA assessment: **PASS", text)

    def test_orphan_attempt_is_ignored_but_published_partial_bundle_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            checkpoint_root = Path(directory) / "checkpoints"
            fingerprint = "b" * 64
            attempts = checkpoint_root / fingerprint / "discovery" / ".attempts"
            attempt = attempts / "locus_0001.crashed"
            attempt.mkdir(parents=True)
            (attempt / "result.rds").write_bytes(b"uncommitted result")
            with patch.object(RAM, "CHECKPOINT_ROOT", checkpoint_root):
                self.assertIsNone(RAM.validate_bundle("discovery", 1, fingerprint))
                bundle = checkpoint_root / fingerprint / "discovery" / "locus_0001"
                os.symlink(os.path.relpath(attempt, bundle.parent), bundle, target_is_directory=True)
                with self.assertRaisesRegex(RAM.ExecutionError, "incomplete"):
                    RAM.validate_bundle("discovery", 1, fingerprint)

    def test_malformed_rds_cannot_receive_a_semantic_attestation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            malformed = Path(directory) / "malformed.rds"
            malformed.write_bytes(b"not an RDS object")
            with self.assertRaisesRegex(RAM.ExecutionError, "semantic validation"):
                RAM._semantic_validation("discovery", 1, "c" * 64, malformed)

    def test_tested_univariate_na_is_rejected_as_an_explicit_checkpoint_failure(self) -> None:
        fingerprint = "d" * 64
        with tempfile.TemporaryDirectory() as directory:
            rds = Path(directory) / "invalid-univariate.rds"
            code = r'''
args <- commandArgs(trailingOnly=TRUE); path <- args[[1L]]; fingerprint <- args[[2L]]
traits <- read.delim("results/track_b/lava_input_info.tsv", stringsAsFactors=FALSE)$phenotype
loc <- read.table("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile", header=TRUE, nrows=1, stringsAsFactors=FALSE)
univ <- data.frame(LOC=loc$LOC, CHR=loc$CHR, START=loc$START, STOP=loc$STOP, phen=traits,
  h2.obs=rep(0.1,8), h2.latent=rep(NA_real_,8), ascertained=rep(FALSE,8), p=c(NA_real_,rep(0.5,7)),
  analysis_status=rep("TESTED",8), error=rep(NA_character_,8), n_snps=rep(10L,8), K=rep(2L,8), stringsAsFactors=FALSE)
bivar <- data.frame(LOC=numeric(0), CHR=numeric(0), START=numeric(0), STOP=numeric(0), pair_id=character(0),
  trait1=character(0), trait2=character(0), discovery_rg=numeric(0), discovery_SE=numeric(0), discovery_P=numeric(0),
  discovery_FDR=numeric(0), local_covariance=numeric(0), rho=numeric(0), rho.lower=numeric(0), rho.upper=numeric(0),
  r2=numeric(0), r2.lower=numeric(0), r2.upper=numeric(0), p=numeric(0), analysis_status=character(0), error=character(0))
status <- data.frame(LOC=loc$LOC, CHR=loc$CHR, START=loc$START, STOP=loc$STOP, status="PROCESSED", n_snps=10L,
  K=2L, univariate_tested=8L, eligible_bivariate_pairs=0L, elapsed_seconds=0.1, stringsAsFactors=FALSE)
saveRDS(list(fingerprint=fingerprint, phase="discovery", locus_index=1L, LOC=as.character(loc$LOC), CHR=as.integer(loc$CHR),
  status=status, univ=univ, bivar=bivar), path)
'''
            built = subprocess.run(
                [str(RAM.R_SCRIPT), "-e", code, str(rds), fingerprint], cwd=ROOT,
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            with self.assertRaisesRegex(RAM.ExecutionError, "numerically invalid"):
                RAM._semantic_validation("discovery", 1, fingerprint, rds)

    def test_complete_ready_bundle_reconciles_immutable_semantic_and_input_attestations(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            checkpoint_root = root / "checkpoints"
            fingerprint = "c" * 64
            attempts = checkpoint_root / fingerprint / "discovery" / ".attempts"
            attempt = attempts / "locus_0001.complete"
            attempt.mkdir(parents=True)
            result_path = attempt / "result.rds"
            result_path.write_bytes(b"synthetic RDS bytes")
            marker = {
                "phase": "discovery", "index": "1", "locus": "101", "chromosome": "1",
                "n_snps": "10", "pair": "NONE", "qc": "PROCESSED",
                "construction_complete": "TRUE",
                "output": str(result_path),
            }
            log = "TRACK_B_LAVA_WORKER\t" + "\t".join(f"{key}={value}" for key, value in marker.items()) + "\n"
            log_path = attempt / "worker.log"
            log_path.write_text(log, encoding="utf-8")
            bundle = checkpoint_root / fingerprint / "discovery" / "locus_0001"
            semantic_text = "synthetic immutable semantic attestation\n"
            semantic_path = attempt / "semantic_validation.txt"
            semantic_path.write_text(semantic_text, encoding="utf-8")
            active_path = root / "active-input.dat"
            active_path.write_bytes(b"frozen input")
            active_stat = RAM.stable_stat(active_path)
            active_spec = {
                "path": "active-input.dat", "bytes": active_path.stat().st_size,
                "sha256": RAM.sha256(active_path), "role": "COMMON_INPUT", "chromosome": None,
            }
            active_receipt = active_spec | {"pre_stat": active_stat, "post_stat": active_stat}
            baseline = checkpoint_root / fingerprint / "active_input_identity_baseline.json"
            baseline.parent.mkdir(parents=True, exist_ok=True)
            baseline.write_text("synthetic baseline\n", encoding="utf-8")
            result_size, result_hash = RAM.file_identity(result_path)
            receipt = {
                "schema_version": 3, "analysis": "LAVA_DISCOVERY", "phase": "discovery",
                "locus_index": 1, "execution_fingerprint": fingerprint, "marker": marker,
                "command": ["Rscript"], "runtime_sec": 1.0, "peak_rss_bytes": 1024,
                "exit_status": 0, "bundle_path": str(bundle.relative_to(root)),
                "log_bytes": log_path.stat().st_size, "log_sha256": RAM.sha256(log_path),
                "semantic_validator_sha256": RAM.sha256(RAM.CHECKPOINT_VALIDATOR),
                "input_baseline_sha256": RAM.sha256(baseline),
                "active_inputs": [active_receipt],
                "artifacts": [
                    {"path": "result.rds", "bytes": result_size, "sha256": result_hash},
                    {"path": "semantic_validation.txt", "bytes": semantic_path.stat().st_size, "sha256": RAM.sha256(semantic_path)},
                ],
            }
            receipt_path = attempt / "receipt.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            (attempt / "READY").write_text(json.dumps({
                "schema_version": 1, "receipt_sha256": RAM.sha256(receipt_path), "state": "READY",
            }), encoding="utf-8")
            os.symlink(os.path.relpath(attempt, bundle.parent), bundle, target_is_directory=True)
            loci = [{"LOC": 101, "CHR": 1, "START": 1, "STOP": 2}]
            semantic_marker = {
                key: value for key, value in marker.items() if key != "output"
            } | {"fingerprint": fingerprint, "status": "PASS"}
            with (
                patch.object(RAM, "ROOT", root),
                patch.object(RAM, "CHECKPOINT_ROOT", checkpoint_root),
                patch.object(RAM, "load_loci", return_value=loci),
                patch.object(RAM, "_semantic_validation", return_value=(semantic_marker, semantic_text)),
                patch.object(RAM, "_active_specs", return_value=[active_spec]),
            ):
                observed = RAM.validate_bundle("discovery", 1, fingerprint)
            self.assertEqual(observed["marker"]["qc"], "PROCESSED")

    def test_ram_gate_includes_high_peak_scientific_fallback_attempt(self) -> None:
        selected = {
            role: {"peak_rss_bytes": 2 * 1024**3}
            for role in ("SMALL", "MEDIAN", "LARGE")
        }
        high_failed_candidate = {
            "peak_rss_bytes": 8 * 1024**3 - 512 * 1024**2,
            "marker": {"construction_complete": "FALSE"},
        }
        with (
            patch.object(RAM, "_representative_receipts", return_value=selected),
            patch.object(RAM, "_expected_roles", return_value=list(selected)),
            patch.object(RAM, "measured_admission_receipts", return_value=list(selected.values()) + [high_failed_candidate]),
            patch.object(RAM, "physical_memory_bytes", return_value=8 * 1024**3),
            patch.object(RAM, "memory_safety_reserve_bytes", return_value=1024**3),
        ):
            with self.assertRaisesRegex(RAM.ExecutionError, "BLOCKED_BY_MEASURED_PER_LOCUS_RAM"):
                RAM.assert_representatives_fit_machine("e" * 64, "discovery")

    def test_execution_evidence_verification_does_not_repair_a_tampered_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = root / "evidence.json"
            report = root / "report.md"
            fingerprint = "f" * 64
            observed = {"execution_fingerprint": fingerprint, "ram_report": {"sha256": "sealed"}}
            evidence.write_text(json.dumps(observed), encoding="utf-8")
            report.write_text("tampered report\n", encoding="utf-8")
            before = (report.read_bytes(), report.stat().st_ino, report.stat().st_mtime_ns)
            with (
                patch.object(RAM, "EXECUTION_EVIDENCE", evidence),
                patch.object(RAM, "quick_execution_fingerprint", return_value=None),
                patch.object(RAM, "execution_evidence_payload", return_value={"execution_fingerprint": fingerprint, "ram_report": {"sha256": "different"}}),
            ):
                with self.assertRaisesRegex(RAM.ExecutionError, "differs"):
                    RAM.verify_execution_evidence()
            after = (report.read_bytes(), report.stat().st_ino, report.stat().st_mtime_ns)
            self.assertEqual(after, before)

    def test_runner_source_preserves_full_family_and_fresh_process_invariants(self) -> None:
        r_source = (ROOT / "scripts/120_run_track_b_lava.R").read_text(encoding="utf-8")
        python_source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("rep.int(1, nrow(univ))", r_source)
        self.assertIn("rep.int(1, nrow(bivar))", r_source)
        self.assertIn('ref.prefix = chromosome_prefix', r_source)
        self.assertIn("CONDITIONER_LOCAL_H2_INELIGIBLE", r_source)
        self.assertIn("gc(full = TRUE)", r_source)
        self.assertIn("RUSAGE_CHILDREN", python_source)
        self.assertIn('for index in range(1, 2496)', python_source)
        self.assertIn("os.symlink", python_source)
        self.assertIn("RAM_BENCHMARK.provenance.json", python_source)


if __name__ == "__main__":
    unittest.main()
