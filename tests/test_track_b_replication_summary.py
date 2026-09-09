import csv
import math
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/track_b/02_independent_global_replication.tsv"


def rows() -> list[dict[str, str]]:
    with OUTPUT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBReplicationSummaryTests(unittest.TestCase):
    def test_canonical_summary_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/122_build_track_b_replication_summary.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_pair_a_terminal_and_pair_b_directional_evidence(self) -> None:
        observed = rows()
        self.assertEqual([row["pair_id"] for row in observed], ["A", "B"])
        self.assertEqual(observed[0]["replication_status"], "NO_VALID_REPLICATION")
        self.assertEqual(observed[1]["replication_status"], "DIRECTIONAL_REPLICATION")
        self.assertEqual(observed[1]["direction"], "CONCORDANT")
        self.assertLess(float(observed[1]["replication_P"]), 0.05)
        self.assertLess(float(observed[1]["replication_CI_lower"]), float(observed[1]["replication_rg"]))
        self.assertGreater(float(observed[1]["replication_CI_upper"]), float(observed[1]["replication_rg"]))
        self.assertGreater(float(observed[1]["heterogeneity_P"]), 0.05)
        self.assertTrue(math.isfinite(float(observed[1]["heterogeneity_z"])))
        self.assertIn("NOT_INDIVIDUAL_LEVEL_VERIFIED", observed[1]["sample_overlap"])

    def test_summary_never_promotes_global_replication_to_mechanism(self) -> None:
        pair_b = rows()[1]
        self.assertIn("does not establish a local shared locus", pair_b["claim_limit"])
        self.assertNotEqual(pair_b["replication_status"], "STRONG_REPLICATION")


if __name__ == "__main__":
    unittest.main()
