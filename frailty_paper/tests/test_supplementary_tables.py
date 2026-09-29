import importlib.util
from pathlib import Path
import re
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "frailty_paper/scripts/39_build_supplementary_tables.py"
spec = importlib.util.spec_from_file_location("supplementary_tables", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class SupplementaryTableBuilderTests(unittest.TestCase):
    def test_current_source_linked_table_counts(self):
        tables, counts = mod.build_tables(ROOT)
        self.assertEqual(counts["table_s1_gwas_metadata.tsv"], 33)
        ffs = [r for r in tables["table_s1_gwas_metadata.tsv"] if r["accession"] == "GCST90295968"]
        self.assertEqual(len(ffs), 1)
        self.assertEqual(ffs[0]["sample_size"], "386565")
        self.assertIn("GRCh37/hg19", ffs[0]["genome_build"])
        self.assertIn("ordinal FFS", ffs[0]["effect_type"])
        self.assertEqual(counts["table_s4_snp_heritability.tsv"], 20)
        self.assertEqual(counts["table_s5_global_rg.tsv"], 168)
        self.assertEqual(counts["table_s7_latent_factor_rg.tsv"], 84)
        self.assertEqual(counts["table_s15_correction_sensitivity.tsv"], 96)
        self.assertEqual(counts["table_s16_software_resources_versions.tsv"], 53)
        self.assertGreaterEqual(counts["table_s16_software_resources_versions.tsv"], 30)

    def test_global_rg_keeps_families_and_claim_limits_separate(self):
        tables, _ = mod.build_tables(ROOT)
        rows = tables["table_s5_global_rg.tsv"]
        self.assertEqual(
            {family: sum(r["analysis_family"] == family for r in rows) for family in {
                "frozen_atlas_sleep_x_FI",
                "secondary_sleep_x_latent_frailty",
                "frozen_atlas_sleep_x_aging_context",
            }},
            {
                "frozen_atlas_sleep_x_FI": 12,
                "secondary_sleep_x_latent_frailty": 84,
                "frozen_atlas_sleep_x_aging_context": 72,
            },
        )
        latent = [r for r in rows if r["analysis_family"] == "secondary_sleep_x_latent_frailty"]
        self.assertTrue(all(r["family_denominator"] == "84" for r in latent))
        self.assertTrue(all(r["interpretation_status"] == "SENSITIVITY_ONLY" for r in latent))
        aging = [r for r in rows if r["analysis_family"] == "frozen_atlas_sleep_x_aging_context"]
        self.assertTrue(all(r["replication_status"] == "READ_ONLY_ATLAS_REUSE_NOT_REPLICATION" for r in aging))

    def test_latent_rg_flags_indicator_dependence(self):
        tables, _ = mod.build_tables(ROOT)
        rows = tables["table_s7_latent_factor_rg.tsv"]
        direct = [r for r in rows if r["sleep_trait"] == "insomnia" and r["disease_trait"] == "frailty_general"]
        self.assertEqual(len(direct), 1)
        self.assertEqual(direct[0]["indicator_overlap_status"], "DIRECT_CONSTRUCT_INDICATOR_MATCH_PART_WHOLE_DEPENDENCE")
        self.assertIn("not independent validation", direct[0]["indicator_overlap_note"])
        related = [r for r in rows if r["sleep_trait"] == "sleepiness"]
        self.assertEqual(len(related), 7)
        self.assertTrue(all("RELATED_SLEEP_INDICATOR" in r["indicator_overlap_status"] for r in related))

    def test_index_marks_unavailable_results_explicitly(self):
        _, counts = mod.build_tables(ROOT)
        status = {r["table"]: r for r in mod.table_status(counts)}
        self.assertEqual(status["S2"]["status"], "NOT_AVAILABLE_UNSCREENED")
        self.assertEqual(status["S6"]["status"], "BLOCKED_SOURCE_PROVENANCE")
        self.assertEqual(status["S9"]["status"], "NOT_JUSTIFIED")
        self.assertEqual(status["S14"]["status"], "NOT_JUSTIFIED")
        self.assertEqual(status["S16"]["status"], "PARTIAL_VERSION_INVENTORY")

    def test_s15_keeps_correction_families_and_denominators_separate(self):
        tables, _ = mod.build_tables(ROOT)
        rows = tables["table_s15_correction_sensitivity.tsv"]
        self.assertEqual(len(rows), 96)
        self.assertEqual(
            {family: sum(r["analysis_family"] == family for r in rows) for family in {
                "primary_sleep_x_FI", "secondary_sleep_x_latent_frailty",
            }},
            {"primary_sleep_x_FI": 12, "secondary_sleep_x_latent_frailty": 84},
        )
        fi = [r for r in rows if r["analysis_family"] == "primary_sleep_x_FI"]
        latent = [r for r in rows if r["analysis_family"] == "secondary_sleep_x_latent_frailty"]
        self.assertTrue(all(r["bh_family_denominator"] == "396" and r["bonferroni_family_size"] == "396" for r in fi))
        self.assertTrue(all(r["bh_family_denominator"] == "84" and r["bonferroni_family_size"] == "84" for r in latent))
        self.assertTrue(all(r["exact_participant_overlap"] == "UNKNOWN" for r in latent))
        self.assertEqual(sum(r["bonferroni_significant_at_0.05"] == "TRUE" for r in latent), 37)

    def test_software_table_separates_pins_from_observed_runs(self):
        rows = mod.software_version_table(ROOT)
        by_name = {r["component"]: r for r in rows}
        self.assertNotIn("3.13.5", by_name["python_workflow"]["observed_version_or_commit"])
        self.assertIn("3.0.1", by_name["ldsc_CBIIT"]["observed_version_or_commit"])
        self.assertEqual(by_name["matplotlib"]["status"], "EXECUTED_MATCHES_PIN")
        self.assertIn("3.11.11", by_name["python_workflow"]["observed_version_or_commit"])
        self.assertEqual(by_name["python_workflow"]["status"], "EXECUTED_MATCHES_PIN")
        self.assertEqual(by_name["snakemake"]["observed_version_or_commit"], "8.30.0")
        self.assertEqual(by_name["lxml"]["status"], "EXECUTED_MATCHES_PIN")
        self.assertEqual(by_name["python_tests"]["observed_version_or_commit"], "3.11.11: 91/91 passed")
        self.assertEqual(by_name["python_clean_venv_tests"]["observed_version_or_commit"], "3.11.11: 91/91 passed; pip check clean")
        self.assertEqual(by_name["python_clean_venv_tests"]["status"], "FRESH_ISOLATED_VENV_SETUP_AND_REVALIDATION_PASS")
        self.assertEqual(by_name["python_integrated_preflight_tests"]["observed_version_or_commit"], "94/94 tests passed; 350/350 registered resource files verified; integrated preflight exit status 0")
        self.assertEqual(by_name["python_integrated_preflight_tests"]["status"], "INTEGRATED_PREFLIGHT_PASS_AT_CAPTURED_SNAPSHOT")
        self.assertEqual(by_name["python_integrated_preflight_tests_latest"]["observed_version_or_commit"], "97/97 tests passed; 351/351 registered resource files verified; 935 missing-abstract PMIDs source-checked; 28 manuscript claims passed; source HEAD 6b906ba")
        self.assertEqual(by_name["python_integrated_preflight_tests_latest"]["status"], "INTEGRATED_PREFLIGHT_PASS_AT_CAPTURED_SNAPSHOT")
        self.assertEqual(by_name["python_post_s16_update_tests"]["observed_version_or_commit"], "97/97 tests passed after S16 builder and regression update; source HEAD 6b906ba")
        self.assertEqual(by_name["python_post_s16_update_tests"]["status"], "FULL_PACKAGE_TEST_SUITE_PASS_AT_UPDATED_SOURCE")
        self.assertEqual(by_name["python_latest_package_suite_tests"]["observed_version_or_commit"], "Python 3.11.11: 169/169 tests passed at 2026-09-27 03:20 UTC")
        self.assertEqual(by_name["python_latest_package_suite_tests"]["status"], "FULL_PACKAGE_TEST_SUITE_PASS_AT_CAPTURED_SNAPSHOT")
        self.assertIn("test_suite_frailty_py311_2026-09-27_0320.log", by_name["python_latest_package_suite_tests"]["evidence"])
        current_suite = by_name["python_current_canonical_package_suite_tests"]
        self.assertIn("Python 3.13.", current_suite["observed_version_or_commit"])
        self.assertIn("source HEAD", current_suite["observed_version_or_commit"])
        self.assertIn("make_test_unittest_2026-09-28_2112.log", current_suite["evidence"])
        suite_log = (ROOT / current_suite["evidence"]).read_text(encoding="utf-8")
        suite_run = re.search(r"Ran (\d+) tests? in ([0-9.]+)s", suite_log)
        self.assertIsNotNone(suite_run)
        self.assertIn(f"{suite_run.group(1)}/{suite_run.group(1)} tests passed", current_suite["observed_version_or_commit"])
        builder_hash = re.search(r"^Supplementary builder SHA-256: ([0-9a-f]{64})$", suite_log, re.MULTILINE)
        tests_hash = re.search(r"^Supplementary builder tests SHA-256: ([0-9a-f]{64})$", suite_log, re.MULTILINE)
        exact_source = (
            builder_hash is not None and tests_hash is not None
            and builder_hash.group(1) == mod.sha256(SCRIPT)
            and tests_hash.group(1) == mod.sha256(Path(__file__).resolve())
        )
        expected_status = "FULL_PACKAGE_TEST_SUITE_PASS_AT_CAPTURED_SOURCE" if exact_source else "PRIOR_SOURCE_SNAPSHOT_NOT_BOUND_TO_CURRENT_WORKTREE"
        self.assertEqual(current_suite["status"], expected_status)
        self.assertEqual(by_name["current_manuscript_validation_audits"]["status"], "CURRENT_MANUSCRIPT_AUDITS_PASS")
        self.assertIn("29 citation uses/27 entries", by_name["current_manuscript_validation_audits"]["observed_version_or_commit"])
        self.assertIn("28/28", by_name["current_manuscript_validation_audits"]["observed_version_or_commit"])
        self.assertIn("92 rows/93 references", by_name["current_manuscript_validation_audits"]["observed_version_or_commit"])
        self.assertIn("91-test workflow validation", by_name["python_workflow"]["observed_version_or_commit"])
        self.assertEqual(by_name["workflow_environment_snapshot"]["status"], "PLATFORM_SCOPED_SNAPSHOT_TESTED")
        self.assertIn("88 Conda packages SHA256-pinned; 45 PyPI packages version-pinned", by_name["workflow_environment_snapshot"]["observed_version_or_commit"])
        self.assertIn("pinned_workflow_environment_validation_2026-09-23.json", by_name["python_workflow"]["evidence"])
        self.assertIn("pinned_workflow_environment_validation_2026-09-23.json", by_name["matplotlib"]["evidence"])
        self.assertEqual(by_name["R"]["status"], "PINNED_NOT_USED_OR_NOT_EVIDENCED_HERE")

    def test_s1_manifest_provenance_hash_covers_only_rows_used_by_s1(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "manifest.tsv"
            path.write_text(
                "resource_id\tresource_type\tnotes\n"
                "gwas1\tGWAS summary statistics\tunchanged\n"
                "metadata1\tResource metadata\tbefore\n",
                encoding="utf-8",
            )
            hash_before, rows_before = mod.gwas_manifest_scope(path)
            path.write_text(
                "resource_id\tresource_type\tnotes\n"
                "gwas1\tGWAS summary statistics\tunchanged\n"
                "metadata1\tResource metadata\tafter\n",
                encoding="utf-8",
            )
            hash_after, rows_after = mod.gwas_manifest_scope(path)
            self.assertEqual((hash_before, rows_before), (hash_after, rows_after))
            path.write_text(
                "resource_id\tresource_type\tnotes\n"
                "gwas1\tGWAS summary statistics\tchanged\n"
                "metadata1\tResource metadata\tafter\n",
                encoding="utf-8",
            )
            hash_changed, rows_changed = mod.gwas_manifest_scope(path)
            self.assertNotEqual(hash_after, hash_changed)
            self.assertEqual(rows_changed, 1)


if __name__ == "__main__":
    unittest.main()
