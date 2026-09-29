from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "41_build_phase16_sensitivity_summary.py"
SPEC = importlib.util.spec_from_file_location("phase16_sensitivity_summary", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Phase16SensitivitySummaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo = Path(__file__).resolve().parents[2]
        cls.rows, cls.inputs = MODULE.build(cls.repo)
        cls.by_claim_domain = {
            (row["conclusion_id"], row["sensitivity_domain"]): row for row in cls.rows
        }

    def test_matrix_covers_all_headlines_and_sensitivity_domains(self) -> None:
        self.assertEqual(len(self.rows), 119)
        self.assertEqual(len({row["conclusion_id"] for row in self.rows}), 7)
        self.assertEqual(len({row["sensitivity_domain"] for row in self.rows}), 17)
        self.assertEqual(len(self.inputs), 9)

    def test_completed_correction_and_workflow_checks_are_evidence_labeled(self) -> None:
        self.assertEqual(
            self.by_claim_domain[("FI_INSOMNIA", "alternate_multiple_testing")]["conclusion_survival"],
            "SURVIVES",
        )
        workflow = self.by_claim_domain[("FI_PANEL", "workflow_reproducibility")]
        self.assertEqual(workflow["conclusion_survival"], "WORKFLOW_CONCORDANT")
        self.assertIn("11/12", workflow["observed_result"])
        self.assertIn("not independent replication", workflow["interpretation"])

    def test_unavailable_tests_are_not_reported_as_negative_results(self) -> None:
        for domain in ("alternate_frailty_definitions", "physical_frailty_components", "LD_reference", "coloc_prior_sensitivity"):
            row = self.by_claim_domain[("FI_INSOMNIA", domain)]
            self.assertEqual(row["conclusion_survival"], "NOT_TESTED")
            self.assertIn("not", row["interpretation"].lower())

    def test_overlap_and_replication_remain_explicitly_unresolved(self) -> None:
        self.assertEqual(self.by_claim_domain[("FI_INSOMNIA", "UKB_overlap")]["conclusion_survival"], "UNRESOLVED")
        self.assertEqual(self.by_claim_domain[("FI_INSOMNIA", "independent_replication")]["conclusion_survival"], "NOT_ESTABLISHED")
        self.assertEqual(self.by_claim_domain[("FI_INSOMNIA", "exact_vs_comparable_replication")]["conclusion_survival"], "NOT_ESTABLISHED")


if __name__ == "__main__":
    unittest.main()
