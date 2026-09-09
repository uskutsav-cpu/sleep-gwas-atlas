import csv
import importlib.util
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRACK_B = ROOT / "results/track_b"
LAVA_REFERENCE = ROOT / "ref/lava/ukb_v1.1/reference.provenance.json"
SCRIPT = ROOT / "scripts/117_build_track_b_blocker_ledger.py"

SPEC = importlib.util.spec_from_file_location("track_b_blocker_ledger", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
LEDGER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LEDGER)


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
        claims = {r["stage_id"]: r["claim_status"] for r in rows}
        self.assertEqual(claims.pop("REPLICATION_A"), "NO_VALID_REPLICATION")
        self.assertEqual(claims.pop("REPLICATION_B"), "DIRECTIONAL_REPLICATION")
        self.assertEqual(set(claims.values()), {"NO_SCIENTIFIC_RESULT"})

    def test_terminal_data_and_current_compute_gates_are_explicit(self) -> None:
        rows = {r["stage_id"]: r for r in read_tsv(TRACK_B / "04_compute_data_blockers.tsv")}
        lock = json.loads((TRACK_B / "04_compute_data_blockers.lock.json").read_text(encoding="utf-8"))
        snapshot = lock["current_compute_snapshot"]
        lava_reference_ok = LAVA_REFERENCE.is_file()
        lava_storage_ok = snapshot["free_bytes"] >= snapshot["lava_minimum_free_bytes"]
        lava_status, lava_classes = LEDGER.classify_lava_readiness(
            lava_reference_ok,
            lava_storage_ok,
            snapshot["lava_ram_benchmark_state"],
        )

        self.assertEqual(rows["REPLICATION_A"]["status"], "NO_VALID_REPLICATION")
        self.assertEqual(rows["REPLICATION_A"]["claim_status"], "NO_VALID_REPLICATION")
        self.assertEqual(rows["REPLICATION_B"]["status"], "COMPLETE")
        self.assertEqual(rows["REPLICATION_B"]["blocker_class"], "NONE")
        self.assertEqual(rows["REPLICATION_B"]["claim_status"], "DIRECTIONAL_REPLICATION")
        self.assertIn("rg=0.3817", rows["REPLICATION_B"]["observed_state"])
        self.assertEqual(rows["LOCAL_LAVA"]["status"], lava_status)
        self.assertEqual(rows["LOCAL_LAVA"]["blocker_class"], ";".join(lava_classes) or "NONE")
        self.assertEqual(rows["SPATIAL"]["status"], "BLOCKED_BY_DATA")
        self.assertIn("synthetic/empty", rows["FINAL_SYNTHESIS"]["observed_state"].lower())

    def test_lava_gate_transitions_after_reference_download(self) -> None:
        self.assertEqual(
            LEDGER.classify_lava_readiness(False, False, "RAM_BENCHMARK_REQUIRED"),
            (
                "BLOCKED_BY_DATA_AND_COMPUTE",
                ["DATA_REFERENCE", "COMPUTE_STORAGE"],
            ),
        )
        self.assertEqual(
            LEDGER.classify_lava_readiness(False, True, "RAM_BENCHMARK_REQUIRED"),
            ("BLOCKED_BY_DATA", ["DATA_REFERENCE"]),
        )
        self.assertEqual(
            LEDGER.classify_lava_readiness(True, False, "RAM_BENCHMARK_REQUIRED"),
            ("READY_FOR_RAM_BENCHMARK", ["RAM_MEASUREMENT_REQUIRED"]),
        )
        self.assertEqual(
            LEDGER.classify_lava_readiness(True, True, "MEASURED_SEQUENTIAL_RAM_PASS"),
            ("READY_TO_RUN", []),
        )
        self.assertEqual(
            LEDGER.classify_lava_readiness(True, True, "BLOCKED_BY_MEASURED_PER_LOCUS_RAM"),
            ("BLOCKED_BY_COMPUTE", ["MEASURED_PER_LOCUS_MEMORY"]),
        )

    def test_lock_binds_replication_evidence_and_forbids_substitution(self) -> None:
        lock = json.loads((TRACK_B / "04_compute_data_blockers.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(
            lock["artifact_role"],
            "EXECUTION_READINESS_WITH_TERMINAL_REPLICATION_STATUS",
        )
        self.assertEqual(lock["row_count"], 16)
        self.assertIn("FORBID", lock["substitution_policy"])
        replication_inputs = {
            "results/track_b/replication/finngen_r13_F5_ADHD_ingest_qc.tsv",
            "results/track_b/replication/finngen_r13_F5_ADHD_ingest_receipt.json",
            "results/track_b/replication/pair_b_ldsc.tsv",
            "results/track_b/replication/pair_b_ldsc.provenance.json",
            "results/track_b/02_independent_global_replication.tsv",
            "results/track_b/02_independent_global_replication.provenance.json",
        }
        self.assertTrue(replication_inputs.issubset(lock["input_sha256"]))

    def test_report_states_terminal_replication_without_mechanistic_promotion(self) -> None:
        report = (TRACK_B / "04_BLOCKER_LEDGER.md").read_text(encoding="utf-8")
        self.assertIn("Pair A replication outcome is terminal `NO_VALID_REPLICATION`", report)
        self.assertIn("earning `DIRECTIONAL_REPLICATION`", report)
        self.assertIn("No local correlation", report)
        self.assertIn("final mechanistic result is present or promoted", report)
        self.assertNotIn("Pair B has a frozen independent-cohort source but no result", report)


if __name__ == "__main__":
    unittest.main()
