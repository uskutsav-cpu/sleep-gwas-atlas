from __future__ import annotations

import unittest
from pathlib import Path


import importlib.util

SPEC = importlib.util.spec_from_file_location(
    "latent_correction_sensitivity",
    Path(__file__).resolve().parents[1] / "scripts" / "40_build_latent_correction_sensitivity.py",
)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


def full_family() -> list[dict[str, str]]:
    return [
        {
            "sleep_trait": sleep,
            "disease_trait": factor,
            "rg": "0.2",
            "se": "0.03",
            "p_value": "0.0001" if sleep == "insomnia" and factor == "frailty_general" else "0.2",
            "q_value": "0.0084" if sleep == "insomnia" and factor == "frailty_general" else "0.5",
            "analysis_family": "secondary_sleep_x_latent_frailty",
            "family_denominator": "84",
            "cohort_overlap_status": "EXPECTED_OR_POSSIBLE; exact overlap UNKNOWN",
        }
        for sleep in module.EXPECTED_SLEEP
        for factor in module.EXPECTED_FACTORS
    ]


class LatentCorrectionSensitivityTests(unittest.TestCase):
    def test_same_family_bonferroni_preserves_bh_and_marks_unknown_overlap(self) -> None:
        rows = module.derive(full_family())
        significant = next(r for r in rows if r["sleep_trait"] == "insomnia" and r["frailty_factor"] == "frailty_general")
        self.assertEqual(significant["bh_q_secondary_84"], "0.0084")
        self.assertEqual(significant["bonferroni_p_secondary_84"], "0.0084")
        self.assertEqual(significant["bh_significant_at_0.05"], "TRUE")
        self.assertEqual(significant["bonferroni_significant_at_0.05"], "TRUE")
        self.assertEqual(significant["exact_participant_overlap"], "UNKNOWN")
        nonsignificant = next(r for r in rows if r["sleep_trait"] == "shortsleep" and r["frailty_factor"] == "frailty_factor_1")
        self.assertEqual(nonsignificant["bh_significant_at_0.05"], "FALSE")
        self.assertEqual(nonsignificant["bonferroni_significant_at_0.05"], "FALSE")

    def test_rejects_incomplete_or_changed_family(self) -> None:
        rows = full_family()
        with self.assertRaisesRegex(ValueError, "12 x 7"):
            module.derive(rows[:-1])
        rows = full_family()
        rows[0]["family_denominator"] = "396"
        with self.assertRaisesRegex(ValueError, "denominator"):
            module.derive(rows)

    def test_rejects_nonfinite_raw_p(self) -> None:
        rows = full_family()
        rows[0]["p_value"] = "nan"
        with self.assertRaisesRegex(ValueError, "p/q"):
            module.derive(rows)


if __name__ == "__main__":
    unittest.main()
