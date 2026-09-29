import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "36_compare_latent_hm3_alleles.py"
SPEC = importlib.util.spec_from_file_location("latent_hm3_alleles", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class LatentHm3AlleleTests(unittest.TestCase):
    def test_orientation_classes(self):
        self.assertEqual(MODULE.classify_alleles("A", "C", "A", "C"), "same_orientation")
        self.assertEqual(MODULE.classify_alleles("C", "A", "A", "C"), "swapped_orientation")
        self.assertEqual(MODULE.classify_alleles("T", "G", "A", "C"), "complemented_same")
        self.assertEqual(MODULE.classify_alleles("G", "T", "A", "C"), "complemented_swapped")
        self.assertEqual(MODULE.classify_alleles("A", "T", "A", "T"), "palindromic_ambiguous")
        self.assertEqual(MODULE.classify_alleles("A", "G", "A", "C"), "allele_mismatch")
        self.assertEqual(MODULE.classify_alleles("AT", "C", "A", "C"), "allele_mismatch")


if __name__ == "__main__":
    unittest.main()
