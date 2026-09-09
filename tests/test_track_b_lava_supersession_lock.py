import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
FINGERPRINT = "ceb36dad6534e5f3256c32e53e2a5f0d8df3a56ad68cab93ebf96843a5fa1461"
DIRECT_FINGERPRINT = "09cd41888ba405ae2b54b8584c5bc15c7d0213e81e7c1cd8f661136bca34a152"
SOURCE_FINGERPRINT = "caeb6b5a1188560f27609cad77801715e0bedfb998660a531c5491afad1177c7"
STOP_REASON = "PERFORMANCE_ONLY_REDUNDANT_SUPERVISOR_VALIDATION_OVERHEAD"
PERFORMANCE_RATIONALE = (
    "OBSERVED_APPROXIMATELY_20_SECOND_CONDITIONAL_START_CADENCE_WITH_2_TO_3_SECOND_"
    "WORKERS;OPTIMIZE_DUPLICATE_CONTRACT_LINEAGE_AGGREGATE_SOURCE_LOOKUP_AND_POST_"
    "PUBLICATION_VALIDATION_WITHOUT_CHANGING_SCIENTIFIC_GATES"
)
STOP_BOUNDARY = (
    "CLEAN_INTERRUPT_BETWEEN_LOCI_DURING_NEXT_AGGREGATE_ACTIVE_INPUT_HASH_VERIFICATION"
)


def load_module():
    path = ROOT / "scripts/144_freeze_track_b_lava_supersession.py"
    specification = importlib.util.spec_from_file_location("lava_supersession_lock_test", path)
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load supersession lock module")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


LOCK = load_module()


class TrackBLavaSupersessionLockTests(unittest.TestCase):
    def test_direct_snapshot_is_exactly_contiguous_one_through_twenty_two(self):
        payload = LOCK.build_direct_lock_payload(
            DIRECT_FINGERPRINT, SOURCE_FINGERPRINT, 22,
            stop_reason=STOP_REASON,
            performance_rationale=PERFORMANCE_RATIONALE,
            stop_boundary=STOP_BOUNDARY,
        )
        family = payload["verified_snapshot"]["conditional_partial_family"]
        self.assertEqual(payload["schema_version"], LOCK.DIRECT_LOCK_SCHEMA)
        self.assertEqual(family["ready_count"], 22)
        self.assertEqual(family["ready_indices"], list(range(1, 23)))
        self.assertTrue(family["no_other_ready_attempts"])
        self.assertEqual(payload["verified_snapshot"]["failed_attempt_count"], 0)
        self.assertFalse(payload["verified_snapshot"]["finalize_present"])
        self.assertFalse(payload["verified_snapshot"]["terminal_qc_present"])
        self.assertEqual(payload["verified_snapshot"]["benchmark"]["row_count"], 23)
        self.assertEqual(
            payload["successor_binding"],
            {
                "mode": "SUCCESSOR_CONTRACT_BINDS_THIS_LOCK",
                "successor_fingerprint": None,
            },
        )
        self.assertEqual(
            family["receipt_result_family_sha256"],
            hashlib.sha256(
                json.dumps(
                    family["members"], sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
            ).hexdigest(),
        )
        with self.assertRaises(LOCK.SupersessionError):
            LOCK.build_direct_lock_payload(
                DIRECT_FINGERPRINT, SOURCE_FINGERPRINT, 21,
                stop_reason=STOP_REASON,
                performance_rationale=PERFORMANCE_RATIONALE,
                stop_boundary=STOP_BOUNDARY,
            )

    def test_real_snapshot_is_exactly_ten_and_has_honest_disclosure(self):
        payload = LOCK.build_lock_payload(FINGERPRINT)
        family = payload["verified_snapshot"]["conditional_partial_family"]
        self.assertEqual(payload["state"], "FINAL_LOCKED")
        self.assertTrue(
            payload["historical_disclosure"][
                "supersession_json_was_edited_during_drafting"
            ]
        )
        self.assertFalse(
            payload["historical_disclosure"]["pre_freeze_draft_versions_attested"]
        )
        self.assertEqual(family["ready_count"], 10)
        self.assertEqual(family["ready_indices"], list(range(1, 11)))
        self.assertTrue(family["no_other_ready_attempts"])
        self.assertEqual(
            family["receipt_result_family_sha256"],
            hashlib.sha256(
                json.dumps(
                    family["members"], sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
            ).hexdigest(),
        )
        legacy = payload["legacy_pre_freeze_digest"]
        self.assertFalse(legacy["recipe_available"])
        self.assertEqual(
            legacy["verification_status"],
            "UNVERIFIABLE_PRE_FREEZE_DRAFT_VALUE_NOT_USED_AS_EVIDENCE",
        )
        self.assertNotEqual(legacy["value_preserved_verbatim"], family["receipt_result_family_sha256"])

    def test_final_lock_verifies_against_current_frozen_state(self):
        payload = LOCK.verify(FINGERPRINT)
        self.assertEqual(payload["superseded_continuation_fingerprint"], FINGERPRINT)

    def test_direct_final_lock_verifies_exact_twenty_two_member_state(self):
        path = LOCK.lock_path(DIRECT_FINGERPRINT)
        before = path.read_bytes()
        payload = LOCK.verify(DIRECT_FINGERPRINT)
        family = payload["verified_snapshot"]["conditional_partial_family"]
        self.assertEqual(payload["schema_version"], LOCK.DIRECT_LOCK_SCHEMA)
        self.assertEqual(payload["superseded_continuation_fingerprint"], DIRECT_FINGERPRINT)
        self.assertEqual(family["ready_indices"], list(range(1, 23)))
        self.assertEqual(
            family["receipt_result_family_sha256"],
            "98af285deb9ec5a7fc45a17ff1424d1b1cbc8cf90eefdf87a31d278bbd2878a6",
        )
        self.assertEqual(path.read_bytes(), before)

    def test_no_replace_refuses_even_identical_final_lock(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "final.lock.json"
            LOCK.publish_no_replace(path, b"sealed\n")
            with self.assertRaises(LOCK.SupersessionError):
                LOCK.publish_no_replace(path, b"sealed\n")
            self.assertEqual(path.read_bytes(), b"sealed\n")

    def test_historical_bundle_rejects_receipt_and_result_tampering(self):
        source = "caeb6b5a1188560f27609cad77801715e0bedfb998660a531c5491afad1177c7"
        run = LOCK.RUN_ROOT / FINGERPRINT
        original = (run / "conditional/locus_0001").resolve(strict=True)
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            receipt_target = temporary_root / "receipt-tamper"
            result_target = temporary_root / "result-tamper"
            shutil.copytree(original, receipt_target)
            shutil.copytree(original, result_target)
            receipt = json.loads((receipt_target / "receipt.json").read_text())
            receipt["runtime_sec"] = float(receipt["runtime_sec"]) + 1
            (receipt_target / "receipt.json").write_text(
                json.dumps(receipt, indent=2, sort_keys=True) + "\n"
            )
            with patch.object(LOCK, "ROOT", temporary_root):
                with self.assertRaisesRegex(LOCK.SupersessionError, "READY marker differs"):
                    LOCK._validate_bundle(
                        receipt_target, "conditional", 1, source, FINGERPRINT
                    )

            with (result_target / "result.rds").open("ab") as handle:
                handle.write(b"tamper")
            with patch.object(LOCK, "ROOT", temporary_root):
                with self.assertRaisesRegex(LOCK.SupersessionError, "artifact differs"):
                    LOCK._validate_bundle(
                        result_target, "conditional", 1, source, FINGERPRINT
                    )


if __name__ == "__main__":
    unittest.main()
