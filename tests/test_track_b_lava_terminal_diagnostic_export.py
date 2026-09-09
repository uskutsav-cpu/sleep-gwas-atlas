import dataclasses
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    specification = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


EXPORTER = load_module(
    "track_b_lava_terminal_diagnostic_export_test",
    "scripts/142_export_track_b_lava_terminal_diagnostics.py",
)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.write_bytes(EXPORTER.tsv_bytes(fields, rows))


class SyntheticTerminalFixture:
    def __init__(self, directory: Path):
        self.terminal = directory / "finalize_unit_all.fixture"
        self.aggregate = directory / "aggregate_unit_all.fixture"
        self.terminal.mkdir()
        self.aggregate.mkdir()
        self.policy = {
            "analysis_id": "track-b-v1.0-local",
            "expected_loci": 2,
            "trait_order": [
                "snoring", "parental_lifespan", "insomnia", "adhd", "frailty",
                "bmi", "sleep_apnea", "mdd",
            ],
            "pair_order": ["A", "B", "CONTROL"],
            "planned_bivariate_pair_locus_family_max": 6,
            "univariate_p_threshold": 0.01,
            "bivariate_multiple_testing": (
                "BH FDR across every locally h2-eligible pair-locus row in the three frozen pairs; "
                "failed eligible attempts contribute p=1"
            ),
            "conditional_multiple_testing": (
                "BH FDR across every model-locus row eligible for conditional testing; failed "
                "eligible attempts contribute p=1 and conditioner-ineligible rows are excluded"
            ),
        }
        self.status = self._status_rows()
        self.univ = self._univ_rows()
        self.bivar = self._bivar_rows()
        self.conditional = self._conditional_rows()
        self.attestation = self._attestation()
        self.checkpoints = (
            {
                "locus_index": 1, "locus": "1", "chromosome": "1", "pair": "B",
                "qc": "CONDITIONER_LOCAL_H2_INELIGIBLE", "receipt_sha256": "1" * 64,
                "result_sha256": "2" * 64,
            },
            {
                "locus_index": 2, "locus": "2", "chromosome": "2", "pair": "NONE",
                "qc": "NOT_BIVARIATE_FDR_ELIGIBLE", "receipt_sha256": "3" * 64,
                "result_sha256": "4" * 64,
            },
        )
        write_tsv(
            self.terminal / "lava_locus_status.tsv", EXPORTER.RESULT_VALIDATOR.LEGACY.STATUS_FIELDS,
            self.status,
        )
        write_tsv(
            self.terminal / "lava_univariate.tsv", EXPORTER.RESULT_VALIDATOR.LEGACY.UNIV_FIELDS,
            self.univ,
        )
        write_tsv(
            self.terminal / "lava_bivariate.tsv", EXPORTER.RESULT_VALIDATOR.LEGACY.BIVAR_FIELDS,
            self.bivar,
        )
        write_tsv(
            self.terminal / "lava_conditional.tsv",
            EXPORTER.RESULT_VALIDATOR.LEGACY.RAW_CONDITIONAL_FIELDS, self.conditional,
        )
        write_tsv(self.aggregate / "conditional_candidates.tsv", EXPORTER.CANDIDATE_FIELDS, [])

    @staticmethod
    def _status_rows() -> list[dict[str, str]]:
        rows = []
        for locus, chromosome, state, n_snps, components, tested, eligible in (
            ("1", "1", "PROCESSED", "100", "20", "7", "1"),
            ("2", "2", "PROCESS_FAILED", "NA", "NA", "0", "0"),
        ):
            rows.append({
                "LOC": locus, "CHR": chromosome, "START": str(int(locus) * 100),
                "STOP": str(int(locus) * 100 + 99), "status": state, "n_snps": n_snps,
                "K": components, "univariate_tested": tested,
                "eligible_bivariate_pairs": eligible, "elapsed_seconds": "1.5",
                "analysis_fingerprint": "a" * 64, "lava_version": "0.1.5",
                "reference_prefix": "fixture/reference",
            })
        return rows

    def _univ_rows(self) -> list[dict[str, str]]:
        rows = []
        for locus, chromosome in (("1", "1"), ("2", "2")):
            for trait in self.policy["trait_order"]:
                if locus == "2":
                    state, p_value, h2, error = "LOCUS_PROCESS_FAILED", "NA", "NA", "fixture process failure"
                elif trait == "mdd":
                    state, p_value, h2, error = "PHENOTYPE_DROPPED", "NA", "NA", "fixture negative variance"
                else:
                    state = "TESTED"
                    p_value = {"insomnia": "0.001", "adhd": "0.002"}.get(trait, "0.5")
                    h2, error = "0.1", "NA"
                rows.append({
                    "LOC": locus, "CHR": chromosome, "START": str(int(locus) * 100),
                    "STOP": str(int(locus) * 100 + 99), "phen": trait, "h2.obs": h2,
                    "h2.latent": "NA", "ascertained": "NA" if state != "TESTED" else "FALSE",
                    "p": p_value, "analysis_status": state, "error": error,
                    "n_snps": "100" if locus == "1" else "NA",
                    "K": "20" if locus == "1" else "NA", "univariate_test_family_n": "16",
                    "univariate_p_threshold": "0.01", "p_bonferroni": "NA",
                    "p_fdr": "NA",
                })
        family_p = [float(row["p"]) if row["analysis_status"] == "TESTED" else 1.0 for row in rows]
        adjusted = EXPORTER.bh(family_p)
        for row, value in zip(rows, adjusted, strict=True):
            if row["analysis_status"] == "TESTED":
                row["p_fdr"] = repr(value)
                row["p_bonferroni"] = repr(min(1.0, float(row["p"]) * len(rows)))
        return rows

    @staticmethod
    def _bivar_rows() -> list[dict[str, str]]:
        return [{
            "LOC": "1", "CHR": "1", "START": "100", "STOP": "199", "pair_id": "B",
            "trait1": "insomnia", "trait2": "adhd", "discovery_rg": "0.4",
            "discovery_SE": "0.03", "discovery_P": "1e-20", "discovery_FDR": "1e-19",
            "local_covariance": "0.02", "rho": "0.3", "rho.lower": "0.1",
            "rho.upper": "0.5", "r2": "0.09", "r2.lower": "0.01", "r2.upper": "0.25",
            "p": "0.01", "analysis_status": "TESTED", "error": "NA",
            "bivariate_test_family_n": "1", "p_fdr": "0.01", "fdr_significant": "TRUE",
        }]

    @staticmethod
    def _conditional_rows() -> list[dict[str, str]]:
        return [{
            "LOC": "1", "CHR": "1", "START": "100", "STOP": "199", "pair_id": "B",
            "trait1": "insomnia", "trait2": "adhd", "conditional_model_id": "B_MDD_ONLY",
            "covariates": "mdd", "pcor": "NA", "ci.lower": "NA", "ci.upper": "NA",
            "p": "NA", "r2.trait1_z": "NA", "r2.trait2_z": "NA",
            "analysis_status": "CONDITIONER_LOCAL_H2_INELIGIBLE",
            "error": "fixture conditioner did not pass local h2", "conditional_test_family_n": "0",
            "p_fdr": "NA", "fdr_significant": "FALSE",
        }]

    @staticmethod
    def _attestation() -> dict:
        reasons = [
            {
                "code": "LOCUS_FAILURE_FRACTION_EXCEEDED", "metric": "locus_failure_fraction",
                "numerator": 1, "denominator": 2, "observed_fraction": 0.5,
                "observed_fraction_decimal": "0.5", "comparison": ">",
                "maximum_allowed_fraction": 0.1, "maximum_allowed_fraction_decimal": "0.1",
            },
            {
                "code": "UNIVARIATE_UNTESTED_FRACTION_EXCEEDED",
                "metric": "univariate_untested_fraction", "numerator": 9, "denominator": 16,
                "observed_fraction": 0.5625, "observed_fraction_decimal": "0.5625",
                "comparison": ">", "maximum_allowed_fraction": 0.1,
                "maximum_allowed_fraction_decimal": "0.1",
            },
        ]
        return {
            "schema_version": "track-b-lava-terminal-qc.1", "state": "TERMINAL_FAILED_QC",
            "analysis_id": "track-b-v1.0-local", "source_discovery_fingerprint": "a" * 64,
            "continuation_execution_fingerprint": "b" * 64, "full_family_complete": True,
            "scientific_validation_passed": False, "canonical_publication_allowed": False,
            "terminal": True,
            "qc": {
                "counts": {
                    "locus": {"failed": 1, "family": 2},
                    "univariate": {"failed_or_untested": 9, "family": 16},
                    "bivariate": {"failed": 0, "eligible_family": 1},
                    "conditional": {
                        "failed": 0, "eligible_family": 0, "ineligible": 1,
                        "reported_family": 1,
                    },
                },
                "reasons": reasons,
            },
        }

    def sealed(self) -> object:
        return EXPORTER.SealedTerminal(
            source_fingerprint="a" * 64, continuation_fingerprint="b" * 64,
            terminal_target=self.terminal, terminal_receipt={},
            terminal_attestation=self.attestation, aggregate_target=self.aggregate,
            aggregate_receipt={}, conditional_checkpoints=self.checkpoints,
            conditional_family_sha256="c" * 64,
            candidate_path=self.aggregate / "conditional_candidates.tsv",
        )


class TrackBLavaTerminalDiagnosticExportTests(unittest.TestCase):
    def validated_fixture(self, temporary: str):
        fixture = SyntheticTerminalFixture(Path(temporary))
        with mock.patch.object(
            EXPORTER.CONTRACT, "legacy_contract",
            return_value=type("Legacy", (), {"load_policy": lambda unused: fixture.policy})(),
        ):
            family = EXPORTER.validate_family(fixture.sealed())
        return fixture, family

    def test_terminal_attestation_cannot_be_relabelled_as_scientific_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = SyntheticTerminalFixture(Path(temporary))
            EXPORTER._require_terminal_attestation(fixture.attestation, "a" * 64, "b" * 64)
            for field, value in (
                ("scientific_validation_passed", True),
                ("canonical_publication_allowed", True),
                ("full_family_complete", False),
                ("state", "PASS"),
            ):
                changed = json.loads(json.dumps(fixture.attestation))
                changed[field] = value
                with self.assertRaises(EXPORTER.ExportError):
                    EXPORTER._require_terminal_attestation(changed, "a" * 64, "b" * 64)

    def test_synthetic_family_derives_counts_and_preserves_exact_bh_denominators(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            _, family = self.validated_fixture(temporary)
            self.assertEqual(family.profile["official_loci"], 2)
            self.assertEqual(family.profile["validated_conditional_checkpoints"], 2)
            self.assertEqual(family.profile["process_failed_loci"], 1)
            self.assertEqual(family.profile["phenotype_dropped"], 1)
            self.assertEqual(family.profile["bivariate_bh_family_n"], 1)
            self.assertEqual(family.profile["bivariate_tested"], 1)
            self.assertEqual(family.profile["conditional_candidate_loci"], 0)
            self.assertEqual(family.profile["conditional_eligible_bh_family_n"], 0)

    def test_full_pair_locus_grid_distinguishes_tested_ineligible_and_failed_locus(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            _, family = self.validated_fixture(temporary)
            grid = EXPORTER.build_bivariate_grid(family)
            self.assertEqual(len(grid), 6)
            statuses = [row["analysis_status"] for row in grid]
            self.assertEqual(statuses.count("TESTED"), 1)
            self.assertEqual(statuses.count("NOT_ELIGIBLE_LOCAL_H2"), 2)
            self.assertEqual(statuses.count("LOCUS_PROCESS_FAILED"), 3)
            for row in grid:
                if row["analysis_status"] not in {"TESTED", "BIVARIATE_FAILED"}:
                    self.assertEqual(row["p"], "NA")
                    self.assertEqual(row["correction_only_family_p"], "NA")
                    self.assertEqual(row["fdr_significant"], "NA")

    def test_failed_eligible_test_keeps_diagnostic_p_separate_from_correction_only_one(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            _, family = self.validated_fixture(temporary)
            failed = dict(family.bivar[0])
            failed.update({
                "analysis_status": "BIVARIATE_FAILED", "p": "0.2", "p_fdr": "NA",
                "fdr_significant": "FALSE", "error": "fixture numerical failure",
            })
            changed = dataclasses.replace(family, bivar=(failed,))
            grid = EXPORTER.build_bivariate_grid(changed)
            tested_scope = next(row for row in grid if row["pair_id"] == "B" and row["LOC"] == "1")
            self.assertEqual(tested_scope["analysis_status"], "BIVARIATE_FAILED")
            self.assertEqual(tested_scope["p"], "0.2")
            self.assertEqual(tested_scope["correction_only_family_p"], "1")
            self.assertIn("NOT_AN_OBSERVED_NULL", tested_scope["correction_only_role"])
            fdr = EXPORTER.build_bivariate_fdr_rows(changed)
            self.assertEqual(len(fdr), 1)
            self.assertEqual(fdr[0]["p"], "0.2")
            self.assertEqual(fdr[0]["correction_only_family_p"], "1")

    def test_rendered_outputs_are_diagnostic_and_preserve_inapplicable_states(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture, family = self.validated_fixture(temporary)
            outputs, derived = EXPORTER.render_outputs(family, fixture.sealed())
            names = {path.name for path in outputs}
            self.assertEqual(names, {
                "local_univariate_results.tsv", "local_bivariate_results.tsv",
                "local_bivariate_fdr.tsv", "local_failures.tsv",
                "local_family_accounting.tsv", "LAVA_DISCOVERY_SUMMARY.md",
            })
            report = outputs[EXPORTER.SUMMARY_PATH].decode("utf-8")
            self.assertIn("TERMINAL_FAILED_QC", report)
            self.assertIn("not a canonical scientific PASS", report)
            self.assertIn("establish causality", report)
            self.assertIn("No skipped locus is called null", report)
            self.assertEqual(derived["bivariate_grid_status_counts"]["TESTED"], 1)
            failures = outputs[EXPORTER.OUTPUT_DIR / "local_failures.tsv"].decode("utf-8")
            self.assertIn("INELIGIBLE_NOT_NULL", failures)
            self.assertIn("INAPPLICABLE_NOT_NULL", failures)

    def test_preflight_fails_closed_before_terminal_ready_exists(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            missing = Path(temporary) / "terminal-qc" / "unit_all"
            forbidden = Path(temporary) / "lava_results.provenance.json"
            with (
                mock.patch.object(EXPORTER, "FORBIDDEN_CANONICAL_PROVENANCE", forbidden),
                mock.patch.object(EXPORTER.SUPERVISOR, "bundle_path", return_value=missing),
            ):
                with self.assertRaisesRegex(EXPORTER.ExportError, "does not exist"):
                    EXPORTER.deep_validate_terminal("a" * 64, "b" * 64)

    def test_no_replace_helper_rejects_differing_existing_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "artifact.tsv"
            EXPORTER._publish_or_verify(path, b"sealed\n")
            EXPORTER._publish_or_verify(path, b"sealed\n")
            with self.assertRaisesRegex(EXPORTER.ExportError, "refusing to replace"):
                EXPORTER._publish_or_verify(path, b"different\n")
            self.assertEqual(path.read_bytes(), b"sealed\n")


if __name__ == "__main__":
    unittest.main()
