from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/28_run_secondary_factor_rg.py"
SPEC = importlib.util.spec_from_file_location("run_secondary_factor_rg", SCRIPT)
RUNNER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(RUNNER)


class SecondaryFactorRgRunTests(unittest.TestCase):
    def test_pair_family_is_exactly_twelve_by_seven(self) -> None:
        self.assertEqual(len(RUNNER.SLEEP_TRAITS), 12)
        self.assertEqual(len(RUNNER.FACTOR_TRAITS), 7)
        self.assertEqual(len(RUNNER.SLEEP_TRAITS) * len(RUNNER.FACTOR_TRAITS), 84)

    def test_resume_requires_seven_rows_not_just_seven_input_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "rg_sleep.log"
            factors = RUNNER.FACTOR_TRAITS
            header = "p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se"
            rows = [
                f"/data/sleep.sumstats.gz /data/{factor}.sumstats.gz 0.1 0.02 5 1e-6 0.1 0.01 1.0 0.01 0.0 0.01"
                for factor in factors
            ]
            text = "Summary of Genetic Correlation Results\n" + header + "\n" + "\n".join(rows) + "\n\nAnalysis finished\n"
            path.write_text(text, encoding="utf-8")
            self.assertTrue(RUNNER.valid_rg_log(path, factors))
            path.write_text(text.replace(rows[-1] + "\n", ""), encoding="utf-8")
            self.assertFalse(RUNNER.valid_rg_log(path, factors))
