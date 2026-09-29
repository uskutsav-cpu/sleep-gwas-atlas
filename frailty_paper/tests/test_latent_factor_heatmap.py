import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "frailty_paper/scripts/32_plot_latent_factor_heatmap.py"
spec = importlib.util.spec_from_file_location("latent_heatmap", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class LatentFactorHeatmapTests(unittest.TestCase):
    def test_output_labels_support_relative_and_external_paths(self):
        relative_outdir = ROOT / "work" / "pinned-figures"
        self.assertEqual(
            mod.output_label(relative_outdir / "plot.png", ROOT),
            "work/pinned-figures/plot.png",
        )
        self.assertEqual(
            mod.output_label(Path(tempfile.gettempdir()) / "plot.png", ROOT),
            str((Path(tempfile.gettempdir()) / "plot.png").resolve()),
        )

    def test_fixed_84_pair_family_and_matrix(self):
        path = ROOT / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
        rg, q, rows = mod.load_family(path)
        self.assertEqual(len(rows), 84)
        self.assertEqual(rg.shape, (12, 7))
        self.assertEqual(q.shape, (12, 7))
        self.assertTrue((q >= 0).all() and (q <= 1).all())

    def test_missing_pair_fails_closed(self):
        source = ROOT / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
        rows = source.read_text().splitlines()
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory) / "partial.tsv"
            partial.write_text("\n".join(rows[:-1]) + "\n")
            with self.assertRaisesRegex(ValueError, "exactly 84"):
                mod.load_family(partial)


if __name__ == "__main__":
    unittest.main()
