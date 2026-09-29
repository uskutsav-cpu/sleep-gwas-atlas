import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/30_build_aging_context_overlap.py"
SPEC = importlib.util.spec_from_file_location("aging_context_overlap", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AgingContextOverlapTests(unittest.TestCase):
    def test_cohort_status_never_calls_unknown_zero(self):
        self.assertEqual(MODULE.pair_status("YES", "YES"), "SHARED_COHORT_EXPECTED")
        self.assertEqual(MODULE.pair_status("NO", "NO"), "NO_SHARED_COHORT_REPORTED")
        self.assertEqual(MODULE.pair_status("NO", "UNKNOWN"), "POSSIBLE_OR_UNKNOWN")

    def test_builds_one_pair_annotation_with_exact_overlap_unknown(self):
        results = [{"pair_id": "insomnia__healthspan", "sleep_trait": "insomnia",
                    "non_sleep_trait": "healthspan"}]
        ledger = [{"trait_id": t, "source_id": t, "uk_biobank": flag,
                   "finngen": "UNKNOWN", "23andme": "NO", "charge": "UNKNOWN"}
                  for t, flag in (("insomnia", "YES"), ("healthspan", "YES"))]
        rows = MODULE.build_rows(results, ledger, expected_count=1)
        self.assertEqual(rows[0]["uk_biobank_pair_status"], "SHARED_COHORT_EXPECTED")
        self.assertEqual(rows[0]["declared_shared_cohorts"], "uk_biobank")
        self.assertEqual(rows[0]["exact_participant_intersection"], "UNKNOWN")
        self.assertEqual(rows[0]["overall_cohort_status"], "SHARED_COHORT_EXPECTED")


if __name__ == "__main__":
    unittest.main()
