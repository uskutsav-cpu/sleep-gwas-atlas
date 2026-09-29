from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/53_build_sleep_power_screen.py"
SPEC = importlib.util.spec_from_file_location("sleep_power_screen", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SleepPowerScreenTests(unittest.TestCase):
    def test_case_control_effective_n_uses_conventional_four_factor(self):
        value = MODULE.case_control_neff(34_184, 305_742)
        self.assertEqual(round(value), 122_985)
        candidate = MODULE.case_control_neff(15_962, 382_950)
        self.assertEqual(round(candidate), 61_293)
        self.assertLess(candidate, value * 0.51)
        self.assertAlmostEqual(candidate / value, 0.49837779, places=7)

    def test_local_eligibility_promotion_gate_is_outcome_blind_and_quantified(self):
        self.assertTrue(MODULE.eligibility_change(1_204, 1_454, 2_495)["passes_promotion_gate"])
        self.assertFalse(MODULE.eligibility_change(1_204, 1_428, 2_495)["passes_promotion_gate"])
        self.assertTrue(MODULE.eligibility_change(1_204, 1_527, 2_495)["passes_promotion_gate"])

    def test_sleep_panel_scope_is_read_from_locked_configuration(self):
        rows = MODULE.read_tsv(MODULE.INPUTS["sleep_panel"])
        self.assertEqual(len([row for row in rows if row["domain"] == "sleep"]), 12)
