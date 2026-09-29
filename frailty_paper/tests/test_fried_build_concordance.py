import unittest
import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/56_audit_fried_build_concordance.py"
SPEC = importlib.util.spec_from_file_location("fried_build_concordance", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
merge_concordance = MODULE.merge_concordance


class FriedBuildConcordanceTests(unittest.TestCase):
    def test_coordinate_intersection_and_unordered_alleles(self):
        summary = [
            (1, 10, "rs10", "A", "G"),
            (1, 30, "rs30", "C", "T"),
        ]
        reference = [
            (1, 10, "rs10", "G", "A"),
            (1, 20, "rs20", "A", "C"),
            (1, 30, "rs30", "C", "T"),
        ]
        result = merge_concordance(summary, reference)
        self.assertEqual(result["coordinate_intersection"], 2)
        self.assertEqual(result["same_coordinate_same_rsid"], 2)
        self.assertEqual(result["same_rsid_swapped_alleles"], 1)
        self.assertEqual(result["same_rsid_same_allele_order"], 1)
        self.assertEqual(result["summary_sort_violations"], 0)
        self.assertEqual(result["reference_sort_violations"], 0)

    def test_unsorted_input_is_detected(self):
        result = merge_concordance(
            [(1, 30, "rs30", "C", "T"), (1, 10, "rs10", "A", "G")],
            [(1, 10, "rs10", "A", "G"), (1, 30, "rs30", "C", "T")],
        )
        self.assertGreater(result["summary_sort_violations"], 0)

    def test_same_coordinate_multiallelic_rows_match_by_rsid(self):
        result = merge_concordance(
            [(1, 10, "rs1", "A", "G"), (1, 10, "rs2", "A", "C")],
            [(1, 10, "rs2", "A", "C"), (1, 10, "rs1", "G", "A")],
        )
        self.assertEqual(result["coordinate_intersection"], 1)
        self.assertEqual(result["same_coordinate_same_rsid"], 2)
        self.assertEqual(result["unmatched_variant_rows_at_shared_coordinates"], 0)
        self.assertEqual(result["same_rsid_same_allele_order"], 1)
        self.assertEqual(result["same_rsid_swapped_alleles"], 1)


if __name__ == "__main__":
    unittest.main()
