from __future__ import annotations

import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/19_scan_acquired_gwas.py"
SPEC = importlib.util.spec_from_file_location("gwas_source_scan", SCRIPT)
SCAN = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(SCAN)


class GwasSourceScanTests(unittest.TestCase):
    def scan(self, text: str) -> dict[str, object]:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "source.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write(text)
            return SCAN.scan_file(path)

    def test_standard_summary_statistics_header_and_domains(self) -> None:
        result = self.scan(
            "variant_id\tchromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\tp_value\teaf\n"
            "rs1\t1\t100\tA\tG\t0.1\t0.02\t0.5\t0.3\n"
        )
        self.assertEqual(result["stream_status"], "GZIP_EOF_OK")
        self.assertEqual(result["header_status"], "REQUIRED_COLUMNS_PRESENT")
        self.assertEqual(result["data_rows"], 1)
        self.assertEqual(result["rows_invalid_eaf"], 0)

    def test_log10p_is_supported_for_deposit_style_files(self) -> None:
        result = self.scan(
            "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tBETA\tSE\tLOG10P\n"
            "1\t100\trs1\tA\tG\t0.3\t0.1\t0.02\t2.0\n"
        )
        self.assertEqual(result["header_status"], "REQUIRED_COLUMNS_PRESENT")
        self.assertEqual(result["rows_missing_p"], 0)
        self.assertEqual(result["rows_invalid_p_or_log10p"], 0)

    def test_a2_is_recognized_as_non_effect_allele(self) -> None:
        result = self.scan(
            "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\n"
            "rs1\t1\t100\tA\tG\t0.1\t0.02\t0.5\n"
        )
        self.assertEqual(result["header_status"], "REQUIRED_COLUMNS_PRESENT")
        self.assertEqual(result["rows_missing_alleles"], 0)

    def test_impossible_values_are_counted_but_not_filtered(self) -> None:
        result = self.scan(
            "rsid\tchr\tpos\ta1\ta0\tbeta1\tse\tp\tfreq1\n"
            "rs1\t1\t0\tA\tG\t0.1\t-0.02\t1.2\t1.2\n"
        )
        self.assertEqual(result["data_rows"], 1)
        self.assertEqual(result["rows_invalid_position"], 1)
        self.assertEqual(result["rows_invalid_beta_or_se"], 1)
        self.assertEqual(result["rows_invalid_p_or_log10p"], 1)
        self.assertEqual(result["rows_invalid_eaf"], 1)
        self.assertIn('"data_row":1', result["anomaly_examples"])

    def test_sleep_trait_specific_beta_headers_are_recognized(self) -> None:
        result = self.scan(
            "SNP\tCHR\tBP\tALLELE1\tALLELE0\tA1FREQ\tBETA_SLEEPDURATION\tSE_SLEEPDURATION\tP_SLEEPDURATION\n"
            "rs1\t1\t100\tA\tG\t0.3\t0.1\t0.02\t0.5\n"
        )
        self.assertEqual(result["header_status"], "REQUIRED_COLUMNS_PRESENT")
        self.assertEqual(result["rows_missing_beta_or_se"], 0)

    def test_odds_ratio_is_log_scaled_for_scan_without_rewriting_source(self) -> None:
        result = self.scan(
            "SNP\tCHR\tBP\tA1\tA2\tOR\tSE\tP\n"
            "rs1\t1\t100\tA\tG\t1.1\t0.02\t0.5\n"
        )
        self.assertEqual(result["header_status"], "REQUIRED_COLUMNS_PRESENT;OR_TO_LOG_EFFECT_FOR_SCAN")
        self.assertEqual(result["rows_missing_beta_or_se"], 0)

    def test_finngen_sleep_apnea_headers_are_recognized(self) -> None:
        result = self.scan(
            "#chrom\tpos\tref\talt\trsids\tsebeta\tpval\tbeta\taf_alt\n"
            "1\t100\tA\tG\trs1\t0.02\t0.5\t0.1\t0.3\n"
        )
        self.assertEqual(result["header_status"], "REQUIRED_COLUMNS_PRESENT")
        self.assertEqual(result["rows_missing_alleles"], 0)


if __name__ == "__main__":
    unittest.main()
