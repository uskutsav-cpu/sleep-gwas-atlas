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
        by_pair: dict[str, list[tuple[str, str]]] = {}
        for row in rows:
            by_pair.setdefault(row["pair_id"], []).append((row["conditional_model_id"], row["covariates"]))
            self.assertEqual(row["selection_timing"], "BEFORE_LOCAL_RESULT_ACCESS")
        self.assertEqual(by_pair, {
            "A": [("A_BMI_ONLY", "bmi"), ("A_SLEEP_APNEA_ONLY", "sleep_apnea")],
            "B": [("B_MDD_ONLY", "mdd")],
            "CONTROL": [("NONE", "NONE")],
        })

    def test_reference_and_compute_snapshot_are_state_aware(self) -> None:
        lock = json.loads((TRACK_B / "local_analysis_input.lock.json").read_text(encoding="utf-8"))
        self.assertFalse(lock["local_results_accessed_before_input_freeze"])
        reference_ready = (ROOT / lock["required_reference_provenance"]).is_file()
        self.assertEqual(lock["reference_status"], "READY" if reference_ready else "BLOCKED_BY_DATA")
        self.assertEqual(lock["reference_minimum_free_bytes"], 37580963840)
        self.assertEqual(
            lock["current_compute_snapshot"]["reference_acquisition_storage_gate_applies"],
            not reference_ready,
        )
        expected_compute = (
            "BLOCKED_BY_COMPUTE_STORAGE"
            if not reference_ready
            and lock["current_compute_snapshot"]["free_bytes"] < lock["reference_minimum_free_bytes"]
            else lock["current_compute_snapshot"]["benchmark_state"]
        )
        self.assertEqual(lock["current_compute_snapshot"]["status"], expected_compute)
        self.assertNotIn("minimum_ram_bytes", lock["current_compute_snapshot"])
        self.assertIn("MEASURE", lock["current_compute_snapshot"]["ram_decision_rule"])

    def test_runtime_and_conditional_family_are_frozen(self) -> None:
        runtime = {row["key"]: row["value"] for row in read_tsv(TRACK_B / "lava_runtime_policy.tsv")}
        self.assertEqual(runtime["expected_traits"], "8")
        self.assertEqual(runtime["expected_pairs"], "3")
        self.assertEqual(runtime["planned_univariate_tests"], str(2495 * 8))
        self.assertEqual(runtime["conditional_attenuation_fraction_partial"], "0.25")
        self.assertEqual(runtime["conditional_attenuation_fraction_full"], "0.75")
        self.assertEqual(runtime["ram_execution_unit"], "WHOLE_PREDECLARED_LAVA_LOCUS")
        self.assertEqual(runtime["maximum_loci_per_process"], "1")
        self.assertIn("FRESH_R_PROCESS_PER_LOCUS", runtime["worker_process_rule"])
        self.assertIn("DO_NOT_APPLY_AN_ASSUMED", runtime["memory_decision_rule"])
        self.assertNotIn("minimum_ram_bytes", runtime)
        self.assertIn("FDR<=0.05", runtime["conditional_execution_gate"])


if __name__ == "__main__":
    unittest.main()
