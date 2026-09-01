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


class TrackBBlockerLedgerTests(unittest.TestCase):
    def test_ledger_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/117_build_track_b_blocker_ledger.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_all_required_stages_are_retained_in_order(self) -> None:
        rows = read_tsv(TRACK_B / "04_compute_data_blockers.tsv")
        self.assertEqual([r["stage_order"] for r in rows], [str(i) for i in range(1, 17)])
        self.assertEqual([r["stage_id"] for r in rows], [
            "REPLICATION_A", "REPLICATION_B", "LOCAL_LAVA", "LOCAL_CONDITIONAL",
            "LOCAL_ROBUSTNESS", "PLEIOTROPY_PLACO", "PLEIOTROPY_CONJFDR",
            "FINE_MAPPING_COLOC", "TISSUE_CELL", "MOLECULAR_QTL_TWAS",
            "REGULATORY_SCATAC", "SPATIAL", "PATHWAYS", "MODEL_FUNCTION",
            "CAUSAL_MR", "FINAL_SYNTHESIS",
        ])
        self.assertTrue(all(r["claim_status"] == "NO_SCIENTIFIC_RESULT" for r in rows))

    def test_terminal_data_and_current_compute_gates_are_explicit(self) -> None:
        rows = {r["stage_id"]: r for r in read_tsv(TRACK_B / "04_compute_data_blockers.tsv")}
        self.assertEqual(rows["REPLICATION_A"]["status"], "NO_VALID_REPLICATION")
        self.assertEqual(rows["REPLICATION_B"]["status"], "BLOCKED_BY_COMPUTE")
        self.assertEqual(rows["LOCAL_LAVA"]["status"], "BLOCKED_BY_DATA_AND_COMPUTE")
        self.assertEqual(rows["SPATIAL"]["status"], "BLOCKED_BY_DATA")
        self.assertIn("synthetic/empty", rows["FINAL_SYNTHESIS"]["observed_state"].lower())

    def test_lock_is_readiness_only_and_forbids_substitution(self) -> None:
        lock = json.loads((TRACK_B / "04_compute_data_blockers.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(lock["artifact_role"], "EXECUTION_READINESS_NOT_SCIENTIFIC_RESULT")
        self.assertEqual(lock["row_count"], 16)
        self.assertIn("FORBID", lock["substitution_policy"])
        self.assertLess(
            lock["current_compute_snapshot"]["free_bytes"],
            lock["current_compute_snapshot"]["pair_b_minimum_free_bytes"],
        )


if __name__ == "__main__":
    unittest.main()
