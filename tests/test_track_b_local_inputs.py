import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRACK_B = ROOT / "results/track_b"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBLocalInputTests(unittest.TestCase):
    def test_local_input_family_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/116_prepare_track_b_local_inputs.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_exact_three_pair_and_eight_trait_scope(self) -> None:
        pairs = read_tsv(TRACK_B / "lava_pair_manifest.tsv")
        self.assertEqual([(r["pair_id"], r["trait1"], r["trait2"]) for r in pairs], [
            ("A", "snoring", "parental_lifespan"),
            ("B", "insomnia", "adhd"),
            ("CONTROL", "insomnia", "frailty"),
        ])
        self.assertTrue(all(r["planned_loci"] == "2495" for r in pairs))
        info = read_tsv(TRACK_B / "lava_input_info.tsv")
        self.assertEqual([r["phenotype"] for r in info], [
            "snoring", "parental_lifespan", "insomnia", "adhd", "frailty",
            "bmi", "sleep_apnea", "mdd",
        ])

    def test_conditioners_are_minimal_and_pre_result(self) -> None:
        rows = read_tsv(TRACK_B / "local_conditional_manifest.tsv")
        by_pair: dict[str, list[str]] = {}
        for row in rows:
            by_pair.setdefault(row["pair_id"], []).append(row["covariate"])
            self.assertEqual(row["selection_timing"], "BEFORE_LOCAL_RESULT_ACCESS")
        self.assertEqual(by_pair, {"A": ["bmi", "sleep_apnea"], "B": ["mdd"], "CONTROL": ["NONE"]})

    def test_missing_reference_is_preserved_as_blocker(self) -> None:
        lock = json.loads((TRACK_B / "local_analysis_input.lock.json").read_text(encoding="utf-8"))
        self.assertFalse(lock["local_results_accessed_before_input_freeze"])
        self.assertEqual(lock["reference_status"], "BLOCKED_BY_DATA")
        self.assertEqual(lock["current_compute_snapshot"]["status"], "BLOCKED_BY_COMPUTE")
        self.assertEqual(lock["reference_minimum_free_bytes"], 37580963840)


if __name__ == "__main__":
    unittest.main()
