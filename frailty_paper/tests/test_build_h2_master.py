import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "frailty_paper/scripts/31_build_h2_master.py"
spec = importlib.util.spec_from_file_location("build_h2_master", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class BuildH2MasterTests(unittest.TestCase):
    def test_exact_frozen_trait_family_and_provenance(self):
        rows, hashes = mod.assemble(ROOT)
        self.assertEqual(len(rows), 20)
        self.assertEqual(len({row["trait_id"] for row in rows}), 20)
        self.assertEqual({r["analysis_family"] for r in rows}, {"frailty_index", "sleep_panel", "latent_factors"})
        self.assertEqual(sum(r["analysis_family"] == "sleep_panel" for r in rows), 12)
        self.assertEqual(sum(r["analysis_family"] == "latent_factors" for r in rows), 7)
        self.assertTrue(all(r["verdict"] == "PASS" and r["raw_sha256"] for r in rows))
        self.assertTrue(all(r["input_log"] in hashes for r in rows))
        self.assertTrue(any(r["scale"] == "liability" for r in rows if r["analysis_family"] == "sleep_panel"))
        self.assertEqual(next(r["scale"] for r in rows if r["trait_id"] == "frailty"), "observed")

    def test_source_row_deletion_fails_closed(self):
        original = mod.SOURCES["sleep_panel"]
        mod.SOURCES["sleep_panel"] = (original[0], original[1] | {"missing_trait"})
        try:
            with self.assertRaisesRegex(ValueError, "trait IDs mismatch"):
                mod.assemble(ROOT)
        finally:
            mod.SOURCES["sleep_panel"] = original


if __name__ == "__main__":
    unittest.main()
