import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBLAVARunnerTests(unittest.TestCase):
    def test_runtime_policy_is_pre_result_and_exact_scope(self) -> None:
        policy = json.loads((ROOT / "config/track_b_local_analysis_policy.json").read_text(encoding="utf-8"))
        self.assertEqual(policy["planned_univariate_tests"], 2495 * 8)
        self.assertEqual(policy["planned_bivariate_pair_locus_family_max"], 2495 * 3)
        self.assertIn("FDR<=0.05", policy["conditional_execution_gate"])
        runtime = {r["key"]: r["value"] for r in read_tsv(ROOT / "results/track_b/lava_runtime_policy.tsv")}
        self.assertEqual(runtime["expected_traits"], "8")
        self.assertEqual(runtime["expected_pairs"], "3")
        self.assertEqual(runtime["conditional_max_r2"], "0.95")

    def test_r_runner_parses_and_uses_exact_track_b_inputs(self) -> None:
        result = subprocess.run(
            [str(ROOT / ".r-env/bin/Rscript"), "-e", "parse(file='scripts/120_run_track_b_lava.R')"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        source = (ROOT / "scripts/120_run_track_b_lava.R").read_text(encoding="utf-8")
        self.assertIn('!identical(pairs$pair_id, c("A", "B", "CONTROL"))', source)
        self.assertIn("run.pcor", source)
        self.assertIn("CONDITIONER_LOCAL_H2_INELIGIBLE", source)
        self.assertIn("scripts/121_validate_track_b_lava.py", source)

    def test_contract_fails_closed_while_reference_is_absent(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/119_track_b_lava_contract.py", "--preflight"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("LAVA reference provenance is unreadable", result.stdout + result.stderr)

    def test_validator_enforces_exact_families_and_conditional_gate(self) -> None:
        result = subprocess.run(
            ["python3", "-m", "py_compile", "scripts/119_track_b_lava_contract.py", "scripts/121_validate_track_b_lava.py"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        source = (ROOT / "scripts/121_validate_track_b_lava.py").read_text(encoding="utf-8")
        self.assertIn("univariate result is not the exact 2,495 x 8 family", source)
        self.assertIn("bivariate result differs from the locally eligible frozen three-pair family", source)
        self.assertIn("conditional result does not exactly cover FDR-supported A/B loci", source)


if __name__ == "__main__":
    unittest.main()
