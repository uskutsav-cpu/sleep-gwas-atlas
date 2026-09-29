import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "audit_power_optimized_family_feasibility_v1.py"
SPEC = importlib.util.spec_from_file_location("family_feasibility_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FamilyFeasibilityAuditTests(unittest.TestCase):
    def test_sleep_replacements_cannot_pass_frozen_family_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = MODULE.audit(Path(temporary) / "feasibility.json")

        self.assertEqual(result["status"], "PASS_SLEEP_ONLY_REPLACEMENTS_CANNOT_PASS_FROZEN_FAMILY_QC")
        self.assertFalse(result["association_results_consulted"])
        lower = result["best_case_lower_bounds"]
        self.assertEqual(lower["one_longsleep_replacement_zero_not_run"]["not_run"], 2429)
        self.assertEqual(lower["one_longsleep_replacement_zero_not_run"]["excess_over_frozen_limit"], 1556)
        self.assertEqual(lower["both_sleep_traits_zero_not_run"]["retained_five_disorder_traits_not_run"], 1758)
        self.assertEqual(lower["both_sleep_traits_zero_not_run"]["excess_over_frozen_limit"], 885)
        self.assertFalse(lower["both_sleep_traits_zero_not_run"]["passes"])


if __name__ == "__main__":
    unittest.main()
