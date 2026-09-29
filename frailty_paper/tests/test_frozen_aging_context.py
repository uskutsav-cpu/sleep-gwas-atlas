import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/29_audit_frozen_aging_context_rg.py"
SPEC = importlib.util.spec_from_file_location("frozen_aging_context", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FrozenAgingContextTests(unittest.TestCase):
    def setUp(self):
        self.sleeps = [f"sleep_{i}" for i in range(12)]
        self.context = ("healthspan", "parental_lifespan")
        self.traits = [
            {"trait_id": trait, "label": trait.title(), "source_id": f"source_{trait}",
             "dataset_version": "v1", "ancestry": "EUR", "build": "hg19",
             "ncase": "NA", "ncontrol": "NA", "n_total": "1000", "h2_scale": "observed",
             "h2": "0.1", "h2_se": "0.01", "h2_z": "10", "ldsc_intercept": "1.0",
             "phase1_verdict": "PASS"}
            for trait in self.context
        ]
        self.pairs = [
            {"pair_id": f"{sleep}__{context}", "sleep_trait": sleep,
             "non_sleep_trait": context, "analysis_tier": "PRIMARY_PHASE1",
             "interpretation_status": "PRIMARY", "global_rg": "0.1", "global_rg_se": "0.02",
             "global_rg_z": "5", "global_rg_p": "0.1", "global_rg_fdr_all_396": "0.2",
             "global_rg_fdr_primary_372": "0.2", "global_rg_primary_significant": "FALSE",
             "effect_direction": "POSITIVE"}
            for sleep in self.sleeps for context in self.context
        ]

    def test_extract_is_exact_sleep_by_context_cartesian_family_and_preserves_q(self):
        rows = MODULE.select_context_rows(self.sleeps, self.pairs, self.traits, self.context)
        self.assertEqual(len(rows), 24)
        self.assertEqual(len({r["pair_id"] for r in rows}), 24)
        self.assertEqual({r["global_rg_fdr_all_396"] for r in rows}, {"0.2"})
        self.assertEqual({r["fdr_denominator"] for r in rows}, {"396 locked sleep-by-non-sleep pairs"})
        self.assertTrue(all(r["reuse_status"] == "READ_ONLY_FROZEN_ATLAS_RESULT" for r in rows))

    def test_missing_frozen_pair_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "family mismatch"):
            MODULE.select_context_rows(self.sleeps, self.pairs[:-1], self.traits, self.context)

    def test_missing_context_trait_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "context traits are absent"):
            MODULE.select_context_rows(self.sleeps, self.pairs, self.traits[:-1], self.context)


if __name__ == "__main__":
    unittest.main()
