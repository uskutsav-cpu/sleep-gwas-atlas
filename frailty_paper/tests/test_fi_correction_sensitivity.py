from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "frailty_paper/scripts/35_build_fi_correction_sensitivity.py"
spec = importlib.util.spec_from_file_location("fi_correction_sensitivity", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


class FICorrectionSensitivityTests(unittest.TestCase):
    def setUp(self):
        self.rows = [
            {"pair_id": f"{trait}__frailty", "sleep_trait": trait,
             "non_sleep_trait": "frailty", "global_rg": "0.1",
             "global_rg_se": "0.02", "global_rg_p": "1e-5",
             "global_rg_fdr_all_396": "0.01" if i < 9 else "0.2"}
            for i, trait in enumerate(sorted(module.EXPECTED_TRAITS))
        ]

    def test_preserves_primary_bh_and_applies_same_family_bonferroni(self):
        result = module.build_rows(self.rows, 396)
        self.assertEqual(len(result), 12)
        self.assertEqual(sum(r["bh_all_396_significant_at_0.05"] == "TRUE" for r in result), 9)
        self.assertEqual(sum(r["bonferroni_all_396_significant_at_0.05"] == "TRUE" for r in result), 12)
        self.assertTrue(all(r["bonferroni_family_size"] == "396" for r in result))
        self.assertTrue(all("primary BH q-values unchanged" in r["comparison_interpretation"] for r in result))

    def test_rejects_changed_family_or_incomplete_trait_panel(self):
        with self.assertRaisesRegex(ValueError, "denominator"):
            module.build_rows(self.rows, 12)
        with self.assertRaisesRegex(ValueError, "exactly the 12 locked"):
            module.build_rows(self.rows[:-1], 396)


if __name__ == "__main__":
    unittest.main()
