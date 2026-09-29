import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "reassess_candidate_replacement_decisions_v1.py"
SPEC = importlib.util.spec_from_file_location("candidate_decision_reassessment", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CandidateDecisionReassessmentTests(unittest.TestCase):
    def test_technical_advancement_is_distinct_from_family_qc(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "candidates.tsv"
            provenance_path = Path(temporary) / "provenance.json"
            provenance = MODULE.build(output, provenance_path)

        self.assertEqual(provenance["status"], "PASS_TECHNICAL_ADVANCEMENT_FAMILY_QC_BLOCKED")
        self.assertEqual(provenance["technical_advancement_candidates"], 1)
        self.assertEqual(provenance["substituted_family_not_run"], 3305)
        self.assertEqual(provenance["substituted_family_maximum_not_run"], 873)
        self.assertEqual(provenance["full_family_qc_passing_candidates"], 0)
        self.assertFalse(provenance["bivariate_lava_started"])
        self.assertEqual(provenance["sleep_only_best_case_both_traits_excess"], 885)

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "candidates.tsv"
            MODULE.build(output, Path(temporary) / "provenance.json")
            text = output.read_text()
        self.assertIn("TECHNICAL ADVANCEMENT CRITERIA MET", text)
        self.assertIn("Bivariate LAVA remains blocked", text)


if __name__ == "__main__":
    unittest.main()
