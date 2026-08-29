import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results/validation/candidate_validation_manifest.tsv"
LOCK = ROOT / "results/validation/candidate_validation_manifest.lock.json"


class CandidateValidationManifestTests(unittest.TestCase):
    def test_exact_six_candidate_family_is_prevalidation_locked(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/100_freeze_candidate_validation_manifest.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with MANIFEST.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        lock = json.loads(LOCK.read_text(encoding="utf-8"))
        self.assertEqual(len(rows), 6)
        self.assertEqual([row["candidate_id"] for row in rows], lock["candidate_ids_in_locked_order"])
        self.assertFalse(lock["validation_results_accessed_before_lock"])
        self.assertTrue(all(
            row["frozen_candidate_family_status"] == "FROZEN_BEFORE_VALIDATION_NO_REPLACEMENT"
            for row in rows
        ))


if __name__ == "__main__":
    unittest.main()
