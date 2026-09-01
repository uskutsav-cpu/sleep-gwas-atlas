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


class TrackBReplicationFreezeTests(unittest.TestCase):
    def test_freeze_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/113_freeze_track_b_replication_sources.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_pair_a_has_no_invalid_replication_credit(self) -> None:
        rows = read_tsv(TRACK_B / "replication_source_audit.tsv")
        pair_a = [row for row in rows if row["pair_id"] == "A"]
        self.assertTrue(any(row["selection_decision"] == "NO_VALID_REPLICATION" for row in pair_a))
        self.assertFalse(any(row["selection_decision"] == "SELECT_PRIMARY_REPLICATION" for row in pair_a))
        self.assertTrue(any(row["discovery_external_relation"] == "ANCESTOR_COMPONENT_OF_TIMMERS_META_ANALYSIS" for row in pair_a))

    def test_pair_b_source_is_exactly_version_pinned_finngen(self) -> None:
        rows = read_tsv(TRACK_B / "replication_source_audit.tsv")
        selected = [row for row in rows if row["selection_decision"] == "SELECT_PRIMARY_REPLICATION"]
        self.assertEqual(len(selected), 1)
        row = selected[0]
        self.assertEqual((row["pair_id"], row["source_id"]), ("B", "finngen_r13_F5_ADHD"))
        self.assertEqual(row["source_generation"], "1777989549884404")
        self.assertEqual(row["etag"], "e7aedc4602071840bbd5abe62d6fcdcc")
        self.assertEqual(row["content_length_bytes"], "802195177")
        self.assertEqual((row["cases"], row["controls"]), ("5559", "489493"))

    def test_source_selection_precedes_results(self) -> None:
        lock = json.loads((TRACK_B / "replication_source.lock.json").read_text(encoding="utf-8"))
        self.assertFalse(lock["replication_results_accessed_before_selection"])
        self.assertEqual(lock["selection_timing"], "BEFORE_REPLICATION_RESULT_ACCESS")
        self.assertTrue(lock["post_result_source_replacement"].startswith("FORBIDDEN"))


if __name__ == "__main__":
    unittest.main()
