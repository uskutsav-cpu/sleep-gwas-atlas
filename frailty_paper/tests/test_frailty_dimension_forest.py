import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "frailty_paper/scripts/33_plot_insomnia_sleep_apnea_frailty_forest.py"
spec = importlib.util.spec_from_file_location("frailty_dimension_forest", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class FrailtyDimensionForestTests(unittest.TestCase):
    def test_relative_output_directory_resolves_from_working_directory(self):
        relative = Path("work") / "pinned-figures"
        self.assertEqual(mod.resolve_output_dir(relative, ROOT), (Path.cwd() / relative).resolve())

    def test_available_endpoints_and_separate_fdr_families(self):
        data = mod.collect(ROOT)
        self.assertEqual(set(data), {"insomnia", "sleep_apnea"})
        for rows in data.values():
            self.assertEqual(len(rows), 8)
            self.assertEqual(rows[0]["q_family"], "all 396 frozen atlas pairs")
            self.assertTrue(all(r["q_family"] == "84 latent-factor pairs" for r in rows[1:]))
            self.assertEqual([r["endpoint"] for r in rows], mod.ENDPOINTS)


if __name__ == "__main__":
    unittest.main()
