import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/152_freeze_track_b_ram_v1_history.py"
SPEC = importlib.util.spec_from_file_location("track_b_ram_v1_history", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
HISTORY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HISTORY)


class TrackBRAMV1HistoryTests(unittest.TestCase):
    def fixture(self, root: Path, rows: int = 2495) -> str:
        fingerprint = "a" * 64
        result = root / "results/track_b"
        result.mkdir(parents=True)
        header = "\t".join(HISTORY.BENCHMARK_FIELDS) + "\n"
        body = "".join(
            f"LAVA_DISCOVERY\tNONE\t{i}\t1\t10\t0.1\t1.0\t0\t{'b' * 64}\n"
            for i in range(1, rows + 1)
        )
        benchmark = (header + body).encode()
        (result / "RAM_BENCHMARK.tsv").write_bytes(benchmark)
        (result / "RAM_BENCHMARK.namespace.json").write_text(json.dumps({
            "schema_version": 1,
            "analysis_id": "track-b-v1.0-local",
            "execution_fingerprint": fingerprint,
        }))
        provenance = {
            "schema_version": 2,
            "analysis_id": "track-b-v1.0-local",
            "execution_fingerprint": fingerprint,
            "benchmark_path": "results/track_b/RAM_BENCHMARK.tsv",
            "benchmark_bytes": len(benchmark),
            "benchmark_sha256": HISTORY.hashlib.sha256(benchmark).hexdigest(),
            "all_lava_attempt_rows": [{} for _ in range(rows)],
        }
        (result / "RAM_BENCHMARK.provenance.json").write_text(json.dumps(provenance))
        (result / "RAM_AWARE_EXECUTION_REPORT.md").write_text(
            "# Track B RAM-aware execution report\n"
        )
        return fingerprint

    def run_with_root(self, root: Path):
        with patch.object(HISTORY, "ROOT", root):
            return HISTORY.freeze(root)

    def test_freeze_copies_exact_bytes_and_receipt_is_last_verifiable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = self.fixture(root)
            receipt = self.run_with_root(root)
            self.assertEqual(receipt["benchmark_row_count"], 2495)
            snapshot = root / HISTORY.HISTORY_ROOT / fingerprint / "RAM_BENCHMARK.tsv"
            source = root / HISTORY.SOURCE_FILES["benchmark"]
            self.assertEqual(snapshot.read_bytes(), source.read_bytes())
            self.assertNotEqual(snapshot.stat().st_ino, source.stat().st_ino)
            source.write_text("future canonical replacement\n")
            with patch.object(HISTORY, "ROOT", root):
                verified = HISTORY.verify_snapshot(root, fingerprint)
            self.assertEqual(verified["execution_fingerprint"], fingerprint)

    def test_rerun_is_idempotent_only_for_identical_history(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = self.fixture(root)
            self.run_with_root(root)
            self.run_with_root(root)
            path = root / HISTORY.HISTORY_ROOT / fingerprint / "RAM_BENCHMARK.tsv"
            os.chmod(path, 0o644)
            path.write_text("tampered\n")
            with patch.object(HISTORY, "ROOT", root):
                with self.assertRaises(HISTORY.HistoryError):
                    HISTORY.verify_snapshot(root, fingerprint)

    def test_truncated_benchmark_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, rows=2494)
            with patch.object(HISTORY, "ROOT", root):
                with self.assertRaisesRegex(HISTORY.HistoryError, "incomplete"):
                    HISTORY.load_source(root)

    def test_retained_process_failure_encoding_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = self.fixture(root)
            path = root / HISTORY.SOURCE_FILES["benchmark"]
            lines = path.read_text().splitlines()
            lines[1] = (
                f"LAVA_DISCOVERY\tNONE\t1\t1\tNA\t0.1\t1.0\t"
                f"0:PROCESS_FAILED\t{'b' * 64}"
            )
            content = ("\n".join(lines) + "\n").encode()
            path.write_bytes(content)
            provenance_path = root / HISTORY.SOURCE_FILES["provenance"]
            provenance = json.loads(provenance_path.read_text())
            provenance["benchmark_bytes"] = len(content)
            provenance["benchmark_sha256"] = HISTORY.hashlib.sha256(content).hexdigest()
            provenance_path.write_text(json.dumps(provenance))
            with patch.object(HISTORY, "ROOT", root):
                observed = HISTORY.load_source(root)
            self.assertEqual(observed["fingerprint"], fingerprint)

    def test_symlinked_source_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            source = root / HISTORY.SOURCE_FILES["report"]
            real = root / "real-report.md"
            source.rename(real)
            source.symlink_to(real)
            with patch.object(HISTORY, "ROOT", root):
                with self.assertRaisesRegex(HISTORY.HistoryError, "symbolic link"):
                    HISTORY.load_source(root)

    def test_stale_but_well_formed_provenance_is_preserved_and_disclosed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            path = root / HISTORY.SOURCE_FILES["provenance"]
            value = json.loads(path.read_text())
            value["benchmark_sha256"] = "c" * 64
            path.write_text(json.dumps(value))
            with patch.object(HISTORY, "ROOT", root):
                observed = HISTORY.load_source(root)
            self.assertFalse(observed["provenance_audit"]["fully_consistent"])
            self.assertFalse(observed["provenance_audit"]["benchmark_sha256_match"])

    def test_invalid_provenance_hash_shape_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            path = root / HISTORY.SOURCE_FILES["provenance"]
            value = json.loads(path.read_text())
            value["benchmark_sha256"] = "invalid"
            path.write_text(json.dumps(value))
            with patch.object(HISTORY, "ROOT", root):
                with self.assertRaisesRegex(HISTORY.HistoryError, "invalid benchmark hash"):
                    HISTORY.load_source(root)


if __name__ == "__main__":
    unittest.main()
