from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/47_audit_lava_interim_receipts.py"
SPEC = importlib.util.spec_from_file_location("lava_interim_receipt_audit", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class LavaInterimReceiptAuditTests(unittest.TestCase):
    def receipt(self, **changes):
        value = {"process_status": "PROCESSED", "fi_status": "TESTED", "sleep_status": "TESTED"}
        value.update(changes)
        return value

    def test_fully_tested_receipt_is_not_flagged(self) -> None:
        flagged, counts, message = MODULE.classify_receipt(self.receipt(), "LOCUS_RECEIPT", "accel_sleep_duration")
        self.assertFalse(flagged)
        self.assertEqual(counts, {})
        self.assertIsNone(message)

    def test_negative_frailty_variance_is_classified_and_flagged(self) -> None:
        value = self.receipt(fi_status="PHENOTYPE_DROPPED")
        log = "Warning: Negative variance estimate for phenotype(s) 'frailty' in locus 8"
        flagged, counts, message = MODULE.classify_receipt(value, log, "accel_sleep_duration")
        self.assertTrue(flagged)
        self.assertEqual(counts, {"frailty_univariate_not_tested": 1, "frailty_negative_variance": 1})
        self.assertIsNone(message)

    def test_negative_sleep_variance_is_classified_and_flagged(self) -> None:
        value = self.receipt(sleep_status="PHENOTYPE_DROPPED")
        log = "Warning: Negative variance estimate for phenotype(s) 'accel_sleep_duration' in locus 8"
        flagged, counts, message = MODULE.classify_receipt(value, log, "accel_sleep_duration")
        self.assertTrue(flagged)
        self.assertEqual(counts, {"sleep_univariate_not_tested": 1, "sleep_negative_variance": 1})
        self.assertIsNone(message)

    def test_null_process_requires_all_phenotype_negative_variance_evidence(self) -> None:
        value = self.receipt(process_status="PROCESS_FAILED", fi_status="PROCESS_FAILED",
                             sleep_status="PROCESS_FAILED", error="process.locus returned NULL")
        log = "Error: Negative variance estimate for all phenotypes in locus 8"
        flagged, counts, message = MODULE.classify_receipt(value, log, "accel_sleep_duration")
        self.assertTrue(flagged)
        self.assertEqual(counts, {"process_failed": 1, "all_phenotypes_negative_variance": 1})
        self.assertIsNone(message)

    def test_known_no_reference_snp_failure_is_classified(self) -> None:
        value = self.receipt(process_status="PROCESS_FAILED", fi_status="PROCESS_FAILED",
                             sleep_status="PROCESS_FAILED", error="process.locus returned NULL")
        flagged, counts, message = MODULE.classify_receipt(
            value, 'Error: none of specified SNP IDs are present in reference data', "accel_sleep_duration")
        self.assertTrue(flagged)
        self.assertEqual(counts, {"process_failed": 1, "no_specified_snps_in_reference": 1})
        self.assertIsNone(message)

    def test_unexpected_process_failure_is_counted_and_retains_error_evidence(self) -> None:
        value = self.receipt(process_status="PROCESS_FAILED", fi_status="PROCESS_FAILED",
                             sleep_status="PROCESS_FAILED", error="process.locus returned NULL")
        flagged, counts, message = MODULE.classify_receipt(
            value, "Error: unexpected numerical failure", "accel_sleep_duration")
        self.assertTrue(flagged)
        self.assertEqual(counts, {"process_failed": 1, "unexpected_process_failure": 1})
        self.assertEqual(message, "receipt: process.locus returned NULL; log: Error: unexpected numerical failure")

    def test_unknown_process_status_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown process status"):
            MODULE.classify_receipt(self.receipt(process_status="UNKNOWN"), "", "accel_sleep_duration")

    def test_latest_prefix_stops_before_missing_or_empty_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            loci = root / "loci"
            logs = root / "logs"
            loci.mkdir()
            logs.mkdir()
            for index in (1, 2, 3, 4):
                value = self.receipt(trait="accel_sleep_duration", locus_index=index)
                (loci / f"locus_{index:04d}.json").write_text(json.dumps(value))
            (logs / "locus_0001.log").write_text("completed")
            (logs / "locus_0002.log").write_text("completed")
            (logs / "locus_0003.log").write_text("")
            (logs / "locus_0004.log").write_text("completed")
            self.assertEqual(MODULE.latest_complete_prefix(root, "accel_sleep_duration", 4), 2)

    def test_latest_prefix_fails_closed_on_identity_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "loci").mkdir()
            (root / "logs").mkdir()
            value = self.receipt(trait="wrong_trait", locus_index=1)
            (root / "loci/locus_0001.json").write_text(json.dumps(value))
            (root / "logs/locus_0001.log").write_text("completed")
            with self.assertRaisesRegex(ValueError, "receipt identity mismatch"):
                MODULE.latest_complete_prefix(root, "accel_sleep_duration", 4)

    def test_latest_prefix_returns_zero_when_no_completed_slot_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(MODULE.latest_complete_prefix(Path(tmp), "accel_sleep_duration", 4), 0)


if __name__ == "__main__":
    unittest.main()
