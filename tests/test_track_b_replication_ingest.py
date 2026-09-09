import csv
import gzip
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz"
QC = ROOT / "results/track_b/replication/finngen_r13_F5_ADHD_ingest_qc.tsv"
RECEIPT = ROOT / "results/track_b/replication/finngen_r13_F5_ADHD_ingest_receipt.json"


class TrackBReplicationIngestTests(unittest.TestCase):
    @unittest.skipUnless(OUTPUT.is_file() and QC.is_file() and RECEIPT.is_file(), "real frozen replication ingest not materialized")
    def test_real_ingest_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/114_stream_track_b_pair_b_replication.py", "--verify", "--acknowledge-network-gib", "0.75"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    @unittest.skipUnless(OUTPUT.is_file(), "real frozen replication ingest not materialized")
    def test_munged_schema_and_minimum_rows(self) -> None:
        with gzip.open(OUTPUT, "rt", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["SNP", "A1", "A2", "Z", "N"])
            count = sum(1 for _ in reader)
        self.assertGreaterEqual(count, 700_000)

    @unittest.skipUnless(RECEIPT.is_file(), "real frozen replication ingest not materialized")
    def test_receipt_preserves_pre_result_state(self) -> None:
        receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
        self.assertFalse(receipt["replication_results_accessed"])
        self.assertEqual(receipt["source_verification"]["verification_status"], "PASS")
        self.assertFalse(receipt["remote_full_resolution_retained_locally"])


if __name__ == "__main__":
    unittest.main()
