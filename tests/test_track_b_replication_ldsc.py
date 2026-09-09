import importlib.util
import math
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/118_run_track_b_pair_b_ldsc.py"
SPEC = importlib.util.spec_from_file_location("track_b_pair_b_ldsc", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class TrackBReplicationLDSCTests(unittest.TestCase):
    def test_preflight_is_state_aware_and_fails_closed_without_ingest(self) -> None:
        MODULE.preflight()
        original = MODULE.REPLICATION_INPUT
        try:
            MODULE.REPLICATION_INPUT = ROOT / "work/fixtures/definitely_missing_pair_b.sumstats.gz"
            with self.assertRaisesRegex(SystemExit, "BLOCKED_BY_DATA"):
                MODULE.preflight()
        finally:
            MODULE.REPLICATION_INPUT = original

    def test_h2_parser_enforces_prespecified_gate(self) -> None:
        text = """Read summary statistics for 1000000 SNPs.
After merging with regression SNP LD, 120000 SNPs remain.
Total Observed scale h2: 0.12 (0.02)
Intercept: 1.04 (0.01)
"""
        parsed = MODULE.parse_h2_log(text)
        self.assertEqual(parsed["h2_gate_status"], "PASS")
        self.assertEqual(parsed["h2_z"], 6.0)
        failing = text.replace("0.12 (0.02)", "0.06 (0.02)")
        self.assertEqual(MODULE.parse_h2_log(failing)["h2_gate_status"], "FAIL")

    def test_rg_parser_recomputes_p_and_requires_frozen_order(self) -> None:
        text = f"""Computing rg for phenotype 2/2
Reading summary statistics from {MODULE.REPLICATION_INPUT} ...
Read summary statistics for 900000 SNPs.
After merging with summary statistics, 880000 SNPs remain.
870000 SNPs with valid alleles.
Summary of Genetic Correlation Results
p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se
{MODULE.SLEEP_INPUT} {MODULE.REPLICATION_INPUT} 0.2 0.05 4.0 0.0001 0.1 0.01 1.02 0.01 0.01 0.005

"""
        parsed = MODULE.parse_rg_log(text)
        self.assertAlmostEqual(parsed["rg_p"], math.erfc(4 / math.sqrt(2)))
        self.assertEqual(parsed["rg_valid_alleles"], 870000)

    def test_result_classification_is_direction_and_power_conservative(self) -> None:
        pair = {"discovery_rg": "0.30"}
        self.assertEqual(
            MODULE.classify(pair, {"rg": 0.2, "rg_p": 0.01}),
            ("CONCORDANT", "DIRECTIONAL_REPLICATION"),
        )
        self.assertEqual(
            MODULE.classify(pair, {"rg": -0.2, "rg_p": 0.01}),
            ("OPPOSITE", "OPPOSITE_DIRECTION"),
        )
        self.assertEqual(
            MODULE.classify(pair, {"rg": 0.1, "rg_p": 0.2})[1],
            "UNDERPOWERED",
        )

    def test_h2_failure_record_has_no_inferred_rg(self) -> None:
        h2 = MODULE.parse_h2_log("""Read summary statistics for 1000000 SNPs.
After merging with regression SNP LD, 120000 SNPs remain.
Total Observed scale h2: 0.06 (0.02)
Intercept: 1.04 (0.01)
""")
        record = {
            "pair_id": "B", "sleep_trait": "insomnia", "external_trait": "adhd",
            "replication_source_id": "finngen_r13_F5_ADHD", **h2,
            "rg": "NA", "rg_se": "NA", "rg_z": "NA", "rg_p": "NA",
            "cross_trait_intercept": "NA", "cross_trait_intercept_se": "NA",
            "rg_input_snps": "NA", "rg_overlap_after_merge": "NA", "rg_valid_alleles": "NA",
            "direction_vs_discovery": "NOT_APPLICABLE",
            "replication_class": "UNDERPOWERED", "claim_limit": "test",
        }
        rendered = MODULE.output_text(record)
        self.assertIn("UNDERPOWERED", rendered)
        self.assertIn("\tNA\tNA\tNA\tNA\t", rendered)


if __name__ == "__main__":
    unittest.main()
