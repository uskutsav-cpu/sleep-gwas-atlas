"""Regression tests for the Brain6 sleep power-screen receipt auditor."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_brain6_sleep_power_pilot.py"
SPEC = importlib.util.spec_from_file_location("audit_power_optimized_sleep_pilot", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load sleep power auditor: {SCRIPT}")
AUDITOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDITOR)


def _write_trait(root: Path, trait: str, tested: int, gate_passed: int | None = None) -> None:
    pair = root / "pairs" / trait
    (pair / "loci").mkdir(parents=True)
    (pair / "logs").mkdir()
    for index in range(1, 101):
        status = "TESTED" if index <= tested else "PHENOTYPE_DROPPED"
        receipt = {
            "trait": trait,
            "locus_index": index,
            "process_status": "PROCESSED",
            "sleep_status": status,
            "sleep_p": (
                AUDITOR.CANONICAL_LOCAL_H2_P_THRESHOLD / 2
                if index <= (tested if gate_passed is None else gate_passed)
                else AUDITOR.CANONICAL_LOCAL_H2_P_THRESHOLD
            ),
            "K": 10,
        }
        (pair / "loci" / f"locus_{index:04d}.json").write_text(json.dumps(receipt))
        (pair / "logs" / f"locus_{index:04d}.log").write_text("completed\n")


def _write_locus_file(path: Path) -> None:
    rows = ["LOC CHR START STOP"]
    rows.extend(f"{index} {(index - 1) % 22 + 1} 1 2" for index in range(1, 2496))
    path.write_text("\n".join(rows) + "\n")


class TestAuditBrain6SleepPowerPilot(unittest.TestCase):
    def test_audit_applies_strict_local_h2_threshold_separately_from_computability(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _write_trait(root, "longsleep", 60, gate_passed=12)
            _write_trait(root, "sleepdur", 70, gate_passed=30)
            locus_file = root / "loci.locfile"
            _write_locus_file(locus_file)

            result = AUDITOR.audit(root, locus_file)

        self.assertEqual(result["matched_locus_count"], 100)
        self.assertEqual(result["chromosomes_covered"], 22)
        self.assertEqual(result["traits"]["longsleep"]["computable_univariate"], 60)
        self.assertEqual(result["traits"]["sleepdur"]["computable_univariate"], 70)
        self.assertEqual(
            result["traits"]["longsleep"]["canonical_threshold_gate_pass_pair_context"], 12
        )
        self.assertEqual(
            result["traits"]["sleepdur"]["canonical_threshold_gate_pass_pair_context"], 30
        )
        self.assertAlmostEqual(result["pair_context_gate_rate_difference_percentage_points"], 18)
        self.assertAlmostEqual(result["relative_reduction_in_pair_context_gate_failures"], 18 / 88)
        self.assertEqual(result["traits"]["sleepdur"]["shared_reference_K"]["median"], 10)
        self.assertIs(result["mutates_run_or_receipts"], False)

    def test_audit_gate_is_strict_at_canonical_threshold(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _write_trait(root, "longsleep", 100, gate_passed=1)
            _write_trait(root, "sleepdur", 100, gate_passed=1)
            # Put exactly-threshold p-values in one otherwise computable receipt.
            receipt_path = root / "pairs" / "sleepdur" / "loci" / "locus_0001.json"
            value = json.loads(receipt_path.read_text())
            value["sleep_p"] = AUDITOR.CANONICAL_LOCAL_H2_P_THRESHOLD
            receipt_path.write_text(json.dumps(value))
            locus_file = root / "loci.locfile"
            _write_locus_file(locus_file)

            result = AUDITOR.audit(root, locus_file)

        self.assertEqual(
            result["traits"]["sleepdur"]["canonical_threshold_gate_pass_pair_context"], 0
        )

    def test_audit_rejects_missing_receipt_log(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _write_trait(root, "longsleep", 60)
            _write_trait(root, "sleepdur", 70)
            locus_file = root / "loci.locfile"
            _write_locus_file(locus_file)
            (root / "pairs" / "sleepdur" / "logs" / "locus_0001.log").unlink()

            with self.assertRaisesRegex(ValueError, "missing/empty paired log"):
                AUDITOR.audit(root, locus_file)


if __name__ == "__main__":
    unittest.main()
