from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/53_audit_lava_full_receipts.py"
SPEC = importlib.util.spec_from_file_location("lava_full_receipt_audit", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FullLavaReceiptAuditTests(unittest.TestCase):
    def _write_receipt(self, root: Path, trait: str, locus: int, value: dict, log: str) -> None:
        pair = root / "pairs" / trait
        (pair / "loci").mkdir(parents=True, exist_ok=True)
        (pair / "logs").mkdir(parents=True, exist_ok=True)
        (pair / "loci" / f"locus_{locus:04d}.json").write_text(json.dumps(value))
        (pair / "logs" / f"locus_{locus:04d}.log").write_text(log)

    def test_audit_reconciles_valid_receipts_and_failure_categories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ok = {
                "schema_version": 1, "trait": "accel_sleep_duration", "pair_order": 1,
                "locus_index": 1, "process_status": "PROCESSED", "univ_status": "TESTED",
                "fi_status": "TESTED", "sleep_status": "TESTED", "bivar_status": "NOT_ELIGIBLE",
                "fi_h2_obs": 0.01, "fi_p": 0.2, "sleep_h2_obs": 0.02, "sleep_p": 0.1,
            }
            failed = {
                "schema_version": 1, "trait": "chronotype", "pair_order": 2,
                "process_status": "PROCESS_FAILED", "univ_status": "UNIVARIATE_FAILED",
                "fi_status": "PROCESS_FAILED", "sleep_status": "PROCESS_FAILED",
                "bivar_status": "NOT_ELIGIBLE", "error": "process.locus returned NULL",
            }
            self._write_receipt(root, "accel_sleep_duration", 1, ok, "completed")
            negative_variance = {
                **failed, "trait": "accel_sleep_duration", "pair_order": 1,
                "locus_index": 2,
            }
            self._write_receipt(
                root, "accel_sleep_duration", 2, negative_variance,
                'Error: Negative variance estimate for all phenotypes in locus 2. This locus cannot be analysed',
            )
            for locus in range(1, 26):
                self._write_receipt(
                    root, "chronotype", locus, {**failed, "locus_index": locus},
                    "Error: none of specified SNP IDs are present in reference data",
                )
            result = MODULE.audit_run(root)
            self.assertEqual(result["receipt_files_validated"], 27)
            self.assertEqual(result["receipt_issues"], [])
            self.assertEqual(result["process_failure_categories"], {
                "all_phenotypes_negative_variance": 1,
                "no_specified_snps_in_reference": 25,
            })
            self.assertEqual(result["process_failures_by_trait"]["chronotype"][
                "process_or_univariate_failure_lower_bound"], 25)
            self.assertEqual(result["process_failures_by_trait"]["chronotype"]["process_failed"], 25)
            self.assertEqual(result["process_failures_by_trait"]["chronotype"][
                "no_specified_snps_in_reference"], 25)
            self.assertNotIn("unexpected_process_failure", result["process_failures_by_trait"]["chronotype"])
            self.assertEqual(result["process_failures_by_trait"]["accel_sleep_duration"][
                "all_phenotypes_negative_variance"], 1)
            self.assertFalse(result["all_traits_pass_1pct_gate"])
            self.assertFalse(result["mutates_run_or_receipts"])

    def test_missing_or_empty_log_is_reported_as_audit_issue(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            receipt = {
                "schema_version": 1, "trait": "accel_sleep_duration", "pair_order": 1,
                "locus_index": 1, "process_status": "PROCESSED", "univ_status": "TESTED",
                "fi_status": "TESTED", "sleep_status": "TESTED", "bivar_status": "NOT_ELIGIBLE",
                "fi_h2_obs": 0.01, "fi_p": 0.2, "sleep_h2_obs": 0.02, "sleep_p": 0.1,
            }
            self._write_receipt(root, "accel_sleep_duration", 1, receipt, "")
            result = MODULE.audit_run(root)
            self.assertEqual(result["receipt_files_validated"], 0)
            self.assertEqual(len(result["receipt_issues"]), 1)
            self.assertIn("missing or empty log", result["receipt_issues"][0])


if __name__ == "__main__":
    unittest.main()
