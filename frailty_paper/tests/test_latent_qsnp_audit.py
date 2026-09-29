import unittest
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "latent_qsnp_audit",
    Path(__file__).resolve().parents[1] / "scripts" / "43_audit_latent_qsnp_pruning.py",
)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(audit)


class LatentQsnpAuditTests(unittest.TestCase):
    def test_threshold_matches_seven_factor_bonferroni(self):
        self.assertAlmostEqual(audit.Q_THRESHOLD, 5e-8 / 7, places=20)

    def test_author_inclusive_q_threshold_includes_exact_boundary(self):
        self.assertTrue(audit.is_q_significant(audit.Q_THRESHOLD))
        self.assertTrue(audit.is_q_significant(audit.Q_THRESHOLD - 1e-20))
        self.assertFalse(audit.is_q_significant(audit.Q_THRESHOLD + 1e-20))

    def test_windows_merge_overlap_and_adjacency_per_chromosome(self):
        windows = audit.merge_windows([("1", 1_000_000), ("1", 3_000_001), ("2", 50)])
        self.assertEqual(windows["1"], [(1, 4_000_001)])
        self.assertEqual(windows["2"], [(1, 1_000_050)])
        self.assertTrue(audit.in_windows("1", 4_000_001, windows))
        self.assertFalse(audit.in_windows("2", 1_000_051, windows))

    def test_chromosome_prefix_is_normalized(self):
        windows = audit.merge_windows([("chr1", 1_500_000)])
        self.assertTrue(audit.in_windows("1", 2_500_000, windows))
        self.assertTrue(audit.in_windows("chr1", 2_500_000, windows))


if __name__ == "__main__":
    unittest.main()
