from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "frailty_paper/scripts/34_build_physical_component_qc_summary.py"
spec = importlib.util.spec_from_file_location("physical_component_qc_summary", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


class PhysicalComponentQCSummaryTests(unittest.TestCase):
    def test_all_components_pass_structure_but_remain_analysis_blocked(self):
        rows, hashes = module.assemble(REPO)
        self.assertEqual(len(rows), 5)
        self.assertEqual({row["resource_id"] for row in rows}, module.EXPECTED)
        self.assertEqual(len(hashes), 3)
        for row in rows:
            self.assertEqual(row["structural_qc"].split()[0], "PASS")
            self.assertEqual(row["analysis_eligibility"], "BLOCKED")
            self.assertEqual(row["duplicate_key_groups"], "0")
            self.assertEqual(row["duplicate_excess_rows"], "0")
            self.assertIn("UNVERIFIED", row["effect_allele_mapping"])

    def test_fail_closed_if_expected_component_is_missing(self):
        original_read_rows = module.read_rows

        def read_rows_without_weight_loss(path):
            rows = original_read_rows(path)
            if path.name == "all_acquired_resources.tsv":
                return [row for row in rows if row.get("resource_id") != "zenodo_14011550_weight_loss"]
            return rows

        with patch.object(module, "read_rows", side_effect=read_rows_without_weight_loss):
            with self.assertRaisesRegex(ValueError, "resource IDs mismatch"):
                module.assemble(REPO)


if __name__ == "__main__":
    unittest.main()
