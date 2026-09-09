import csv
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch


ROOT = Path(__file__).resolve().parents[1]


def load_validator():
    spec = importlib.util.spec_from_file_location("track_b_lava_validator_test", ROOT / "scripts/121_validate_track_b_lava.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load Track B LAVA validator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class FakeSupervisor:
    def __init__(self, checkpoint_dir: Path, fingerprint: str, loci: int) -> None:
        self.checkpoint_dir = checkpoint_dir
        self.fingerprint = fingerprint
        self.loci = loci

    def bundle_path(self, phase: str, index: int | None, fingerprint: str) -> Path:
        assert fingerprint == self.fingerprint
        name = f"locus_{index:04d}" if index is not None else "unit_all"
        return self.checkpoint_dir / fingerprint / phase / name

    @staticmethod
    def _safe_attempt_target(bundle: Path) -> Path:
        return bundle

    def validate_bundle(self, phase: str, index: int | None, fingerprint: str):
        bundle = self.bundle_path(phase, index, fingerprint)
        if not bundle.is_dir():
            return None
        return json.loads((bundle / "receipt.json").read_text(encoding="utf-8"))


def write_contract_checkpoint_family(checkpoint_dir: Path, loci: int, fingerprint: str) -> FakeSupervisor:
    supervisor = FakeSupervisor(checkpoint_dir, fingerprint, loci)
    for phase in ("discovery", "aggregate-discovery", "conditional", "finalize"):
        indices = range(1, loci + 1) if phase in {"discovery", "conditional"} else (None,)
        for index in indices:
            bundle = supervisor.bundle_path(phase, index, fingerprint)
            bundle.mkdir(parents=True, exist_ok=True)
            receipt = {
                "peak_rss_bytes": 1024,
                "runtime_sec": 0.1,
                "marker": {"qc": "SYNTHETIC_TEST_CHECKPOINT"},
            }
            (bundle / "receipt.json").write_text(json.dumps(receipt) + "\n", encoding="utf-8")
            (bundle / "READY").write_text("synthetic ready checkpoint\n", encoding="utf-8")
    return supervisor


def synthetic_policy() -> dict[str, object]:
    return {
        "analysis_id": "track-b-synthetic-test",
        "lava_version": "0.1.5",
        "expected_loci": 6,
        "planned_univariate_tests": 6 * 8,
        "trait_order": [
            "snoring", "parental_lifespan", "insomnia", "adhd",
            "frailty", "bmi", "sleep_apnea", "mdd",
        ],
        "pair_order": ["A", "B", "CONTROL"],
        "conditional_models": {
            "A": [
                {"conditional_model_id": "A_BMI_ONLY", "covariates": ["bmi"]},
                {"conditional_model_id": "A_SLEEP_APNEA_ONLY", "covariates": ["sleep_apnea"]},
            ],
            "B": [{"conditional_model_id": "B_MDD_ONLY", "covariates": ["mdd"]}],
            "CONTROL": [],
        },
        "reference_prefix": "ref/synthetic/lava",
        "univariate_p_threshold": 0.01,
        "bivariate_fdr_alpha": 0.05,
        "maximum_locus_failure_fraction": 1.0,
        "maximum_univariate_untested_fraction": 1.0,
        "maximum_bivariate_failure_fraction": 1.0,
        "maximum_conditional_failure_fraction": 1.0,
        "runtime": {"conditional_max_r2": 0.95},
        "interpretation_policy": {
            "conditional_attenuation_fraction_partial": 0.25,
            "conditional_attenuation_fraction_full": 0.75,
        },
    }


def build_synthetic_result_family(validator, root: Path, fingerprint: str) -> tuple[dict[str, object], list[Path]]:
    policy = synthetic_policy()
    loci = [
        {"LOC": "101", "CHR": "1", "START": "100", "STOP": "199"},
        {"LOC": "102", "CHR": "1", "START": "200", "STOP": "299"},
        {"LOC": "103", "CHR": "1", "START": "300", "STOP": "399"},
        {"LOC": "104", "CHR": "2", "START": "100", "STOP": "199"},
        {"LOC": "105", "CHR": "2", "START": "200", "STOP": "299"},
        {"LOC": "106", "CHR": "2", "START": "300", "STOP": "399"},
    ]
    locus_path = root / "ref/loci.locfile"
    locus_path.parent.mkdir(parents=True, exist_ok=True)
    locus_path.write_text(
        "LOC CHR START STOP\n"
        + "".join(f"{row['LOC']} {row['CHR']} {row['START']} {row['STOP']}\n" for row in loci),
        encoding="utf-8",
    )

    binary_traits = {"snoring", "insomnia", "adhd", "sleep_apnea", "mdd"}
    input_rows = [
        {"phenotype": trait, "cases": "100" if trait in binary_traits else "NA"}
        for trait in policy["trait_order"]
    ]
    write_tsv(root / "results/track_b/lava_input_info.tsv", ["phenotype", "cases"], input_rows)

    pairs = [
        {"pair_id": "A", "trait1": "snoring", "trait2": "parental_lifespan", "discovery_rg": "-0.16", "discovery_SE": "0.03", "discovery_P": "1e-8", "discovery_FDR": "3e-8"},
        {"pair_id": "B", "trait1": "insomnia", "trait2": "adhd", "discovery_rg": "0.4", "discovery_SE": "0.03", "discovery_P": "1e-30", "discovery_FDR": "3e-30"},
        {"pair_id": "CONTROL", "trait1": "insomnia", "trait2": "frailty", "discovery_rg": "0.64", "discovery_SE": "0.02", "discovery_P": "1e-100", "discovery_FDR": "3e-100"},
    ]
    pair_fields = [
        "pair_id", "trait1", "trait2", "discovery_rg", "discovery_SE", "discovery_P", "discovery_FDR",
    ]
    write_tsv(root / "results/track_b/lava_pair_manifest.tsv", pair_fields, pairs)
    models = [
        {"pair_id": "A", "conditional_model_id": "A_BMI_ONLY", "covariates": "bmi"},
        {"pair_id": "A", "conditional_model_id": "A_SLEEP_APNEA_ONLY", "covariates": "sleep_apnea"},
        {"pair_id": "B", "conditional_model_id": "B_MDD_ONLY", "covariates": "mdd"},
        {"pair_id": "CONTROL", "conditional_model_id": "NONE", "covariates": "NONE"},
    ]
    write_tsv(
        root / "results/track_b/local_conditional_manifest.tsv",
        ["pair_id", "conditional_model_id", "covariates"], models,
    )

    coordinate_index = {row["LOC"]: row for row in loci}
    status_specs = {
        "101": ("PROCESSED", "100", "10", "8", "1"),
        "102": ("PROCESS_FAILED", "NA", "NA", "0", "0"),
        "103": ("PROCESSED", "120", "12", "8", "0"),
        "104": ("PROCESSED", "130", "13", "8", "1"),
        "105": ("PROCESSED", "140", "14", "8", "1"),
        "106": ("PROCESSED", "150", "15", "8", "1"),
    }
    status_rows = []
    for coordinates in loci:
        status, n_snps, components, univariate_tested, eligible_pairs = status_specs[coordinates["LOC"]]
        status_rows.append({
            **coordinates, "status": status, "n_snps": n_snps, "K": components,
            "univariate_tested": univariate_tested, "eligible_bivariate_pairs": eligible_pairs,
            "elapsed_seconds": "0.01", "analysis_fingerprint": fingerprint,
            "lava_version": policy["lava_version"], "reference_prefix": policy["reference_prefix"],
        })

    significant_traits = {
        "101": {"snoring", "parental_lifespan", "bmi"},
        "102": set(),
        "103": set(),
        "104": {"insomnia", "adhd", "mdd"},
        "105": {"snoring", "parental_lifespan"},
        "106": {"insomnia", "adhd", "mdd"},
    }
    univ_rows = []
    univ_family_p = []
    for coordinates in loci:
        failed_locus = coordinates["LOC"] == "102"
        n_snps, components = status_specs[coordinates["LOC"]][1:3]
        for trait in policy["trait_order"]:
            if failed_locus:
                row = {
                    **coordinates, "phen": trait, "h2.obs": "NA", "h2.latent": "NA",
                    "ascertained": "NA", "p": "NA", "analysis_status": "LOCUS_PROCESS_FAILED",
                    "error": "synthetic process failure", "n_snps": "NA", "K": "NA",
                    "univariate_test_family_n": str(policy["planned_univariate_tests"]),
                    "univariate_p_threshold": "0.01", "p_bonferroni": "NA", "p_fdr": "NA",
                }
                univ_family_p.append(1.0)
            else:
                p_value = 0.001 if trait in significant_traits[coordinates["LOC"]] else 0.5
                row = {
                    **coordinates, "phen": trait, "h2.obs": "0.1",
                    "h2.latent": "0.12" if trait in binary_traits else "NA",
                    "ascertained": "TRUE" if trait in binary_traits else "FALSE", "p": f"{p_value:.12g}",
                    "analysis_status": "TESTED", "error": "NA", "n_snps": n_snps, "K": components,
                    "univariate_test_family_n": str(policy["planned_univariate_tests"]),
                    "univariate_p_threshold": "0.01",
                    "p_bonferroni": f"{min(1.0, p_value * policy['planned_univariate_tests']):.12g}",
                    "p_fdr": "PENDING",
                }
                univ_family_p.append(p_value)
            univ_rows.append(row)
    for row, adjusted in zip(univ_rows, validator.bh(univ_family_p), strict=True):
        if row["analysis_status"] == "TESTED":
            row["p_fdr"] = f"{adjusted:.12g}"

    pair_index = {row["pair_id"]: row for row in pairs}
    bivar_specs = [
        ("101", "A", "TESTED", 0.001, 0.4),
        ("104", "B", "TESTED", 0.002, -0.3),
        ("105", "A", "BIVARIATE_FAILED", None, None),
        ("106", "B", "TESTED", 0.003, 0.5),
    ]
    bivar_rows = []
    bivar_family_p = []
    for locus, pair_id, analysis_status, p_value, rho in bivar_specs:
        coordinates, pair = coordinate_index[locus], pair_index[pair_id]
        row = {
            **coordinates, "pair_id": pair_id, "trait1": pair["trait1"], "trait2": pair["trait2"],
            "discovery_rg": pair["discovery_rg"], "discovery_SE": pair["discovery_SE"],
            "discovery_P": pair["discovery_P"], "discovery_FDR": pair["discovery_FDR"],
            "local_covariance": "NA", "rho": "NA", "rho.lower": "NA", "rho.upper": "NA",
            "r2": "NA", "r2.lower": "NA", "r2.upper": "NA", "p": "NA",
            "analysis_status": analysis_status,
            "error": "synthetic bivariate failure" if analysis_status == "BIVARIATE_FAILED" else "NA",
            "bivariate_test_family_n": str(len(bivar_specs)), "p_fdr": "NA", "fdr_significant": "FALSE",
        }
        if analysis_status == "TESTED":
            row.update({
                "local_covariance": "0.02", "rho": f"{rho:.12g}",
                "rho.lower": f"{rho - 0.1:.12g}", "rho.upper": f"{rho + 0.1:.12g}",
                "r2": f"{rho * rho:.12g}", "r2.lower": "0.01", "r2.upper": "0.36",
                "p": f"{p_value:.12g}",
            })
            bivar_family_p.append(p_value)
        else:
            bivar_family_p.append(1.0)
        bivar_rows.append(row)
    for row, adjusted in zip(bivar_rows, validator.bh(bivar_family_p), strict=True):
        if row["analysis_status"] == "TESTED":
            row["p_fdr"] = f"{adjusted:.12g}"
            row["fdr_significant"] = str(adjusted <= policy["bivariate_fdr_alpha"]).upper()

    def conditional_row(locus: str, pair_id: str, model_id: str, covariates: str) -> dict[str, str]:
        coordinates, pair = coordinate_index[locus], pair_index[pair_id]
        return {
            **coordinates, "pair_id": pair_id, "trait1": pair["trait1"], "trait2": pair["trait2"],
            "conditional_model_id": model_id, "covariates": covariates,
            "pcor": "NA", "ci.lower": "NA", "ci.upper": "NA", "p": "NA",
            "r2.trait1_z": "NA", "r2.trait2_z": "NA", "analysis_status": "PENDING",
            "error": "NA", "conditional_test_family_n": "3", "p_fdr": "NA", "fdr_significant": "FALSE",
        }

    unstable = conditional_row("101", "A", "A_BMI_ONLY", "bmi")
    unstable.update({
        "pcor": "0.2", "ci.lower": "0.1", "ci.upper": "0.3", "r2.trait1_z": "0.95",
        "r2.trait2_z": "0.2", "analysis_status": "CONDITIONAL_UNSTABLE_MAX_R2",
        "error": "synthetic max-r2 instability",
    })
    ineligible = conditional_row("101", "A", "A_SLEEP_APNEA_ONLY", "sleep_apnea")
    ineligible.update({
        "analysis_status": "CONDITIONER_LOCAL_H2_INELIGIBLE",
        "error": "synthetic conditioner failed the local h2 gate",
    })
    tested = conditional_row("104", "B", "B_MDD_ONLY", "mdd")
    tested.update({
        "pcor": "-0.25", "ci.lower": "-0.4", "ci.upper": "-0.1", "p": "0.001",
        "r2.trait1_z": "0.2", "r2.trait2_z": "0.3", "analysis_status": "TESTED",
        "p_fdr": f"{validator.bh([1.0, 0.001, 1.0])[1]:.12g}", "fdr_significant": "TRUE",
    })
    failed = conditional_row("106", "B", "B_MDD_ONLY", "mdd")
    failed.update({
        "pcor": "0.2", "ci.lower": "0.1", "ci.upper": "0.3", "p": "0.2",
        "r2.trait1_z": "0.4", "r2.trait2_z": "0.5", "analysis_status": "CONDITIONAL_FAILED",
        "error": "synthetic conditional failure",
    })
    conditional_rows = [unstable, ineligible, tested, failed]

    staged_dir = root / "checkpoints" / fingerprint / "finalize" / "unit_all"
    staged = [
        staged_dir / "lava_locus_status.tsv", staged_dir / "lava_univariate.tsv",
        staged_dir / "lava_bivariate.tsv", staged_dir / "lava_conditional.tsv",
        staged_dir / "04_lava_local_results.tsv", staged_dir / "05_local_conditional_results.tsv",
    ]
    write_tsv(staged[0], validator.STATUS_FIELDS, status_rows)
    write_tsv(staged[1], validator.UNIV_FIELDS, univ_rows)
    write_tsv(staged[2], validator.BIVAR_FIELDS, bivar_rows)
    write_tsv(staged[3], validator.RAW_CONDITIONAL_FIELDS, conditional_rows)
    return policy, staged


def build_zero_hit_result_family(validator, root: Path, fingerprint: str) -> tuple[dict[str, object], list[Path]]:
    policy = synthetic_policy()
    policy["expected_loci"] = 1
    policy["planned_univariate_tests"] = 8
    coordinates = {"LOC": "201", "CHR": "3", "START": "100", "STOP": "199"}
    locus_path = root / "ref/loci.locfile"
    locus_path.parent.mkdir(parents=True, exist_ok=True)
    locus_path.write_text("LOC CHR START STOP\n201 3 100 199\n", encoding="utf-8")

    binary_traits = {"snoring", "insomnia", "adhd", "sleep_apnea", "mdd"}
    write_tsv(
        root / "results/track_b/lava_input_info.tsv", ["phenotype", "cases"],
        [
            {"phenotype": trait, "cases": "100" if trait in binary_traits else "NA"}
            for trait in policy["trait_order"]
        ],
    )
    pairs = [
        {"pair_id": "A", "trait1": "snoring", "trait2": "parental_lifespan", "discovery_rg": "-0.16", "discovery_SE": "0.03", "discovery_P": "1e-8", "discovery_FDR": "3e-8"},
        {"pair_id": "B", "trait1": "insomnia", "trait2": "adhd", "discovery_rg": "0.4", "discovery_SE": "0.03", "discovery_P": "1e-30", "discovery_FDR": "3e-30"},
        {"pair_id": "CONTROL", "trait1": "insomnia", "trait2": "frailty", "discovery_rg": "0.64", "discovery_SE": "0.02", "discovery_P": "1e-100", "discovery_FDR": "3e-100"},
    ]
    write_tsv(
        root / "results/track_b/lava_pair_manifest.tsv",
        ["pair_id", "trait1", "trait2", "discovery_rg", "discovery_SE", "discovery_P", "discovery_FDR"],
        pairs,
    )
    write_tsv(
        root / "results/track_b/local_conditional_manifest.tsv",
        ["pair_id", "conditional_model_id", "covariates"],
        [
            {"pair_id": "A", "conditional_model_id": "A_BMI_ONLY", "covariates": "bmi"},
            {"pair_id": "A", "conditional_model_id": "A_SLEEP_APNEA_ONLY", "covariates": "sleep_apnea"},
            {"pair_id": "B", "conditional_model_id": "B_MDD_ONLY", "covariates": "mdd"},
            {"pair_id": "CONTROL", "conditional_model_id": "NONE", "covariates": "NONE"},
        ],
    )
    status = [{
        **coordinates, "status": "PROCESSED", "n_snps": "100", "K": "10",
        "univariate_tested": "8", "eligible_bivariate_pairs": "0", "elapsed_seconds": "0.01",
        "analysis_fingerprint": fingerprint, "lava_version": policy["lava_version"],
        "reference_prefix": policy["reference_prefix"],
    }]
    univ = [
        {
            **coordinates, "phen": trait, "h2.obs": "0.1",
            "h2.latent": "0.12" if trait in binary_traits else "NA",
            "ascertained": "TRUE" if trait in binary_traits else "FALSE", "p": "0.5",
            "analysis_status": "TESTED", "error": "NA", "n_snps": "100", "K": "10",
            "univariate_test_family_n": "8", "univariate_p_threshold": "0.01",
            "p_bonferroni": "1", "p_fdr": "0.5",
        }
        for trait in policy["trait_order"]
    ]
    staged_dir = root / "checkpoints" / fingerprint / "finalize" / "unit_all"
    staged = [
        staged_dir / "lava_locus_status.tsv", staged_dir / "lava_univariate.tsv",
        staged_dir / "lava_bivariate.tsv", staged_dir / "lava_conditional.tsv",
        staged_dir / "04_lava_local_results.tsv", staged_dir / "05_local_conditional_results.tsv",
    ]
    write_tsv(staged[0], validator.STATUS_FIELDS, status)
    write_tsv(staged[1], validator.UNIV_FIELDS, univ)
    write_tsv(staged[2], validator.BIVAR_FIELDS, [])
    write_tsv(staged[3], validator.RAW_CONDITIONAL_FIELDS, [])
    return policy, staged


class TrackBLAVARunnerTests(unittest.TestCase):
    def test_runtime_policy_is_pre_result_and_exact_scope(self) -> None:
        policy = json.loads((ROOT / "config/track_b_local_analysis_policy.json").read_text(encoding="utf-8"))
        self.assertEqual(policy["planned_univariate_tests"], 2495 * 8)
        self.assertEqual(policy["planned_bivariate_pair_locus_family_max"], 2495 * 3)
        self.assertIn("FDR<=0.05", policy["conditional_execution_gate"])
        runtime = {r["key"]: r["value"] for r in read_tsv(ROOT / "results/track_b/lava_runtime_policy.tsv")}
        self.assertEqual(runtime["expected_traits"], "8")
        self.assertEqual(runtime["expected_pairs"], "3")
        self.assertEqual(runtime["conditional_max_r2"], "0.95")

    def test_r_runner_parses_and_uses_exact_track_b_inputs(self) -> None:
        result = subprocess.run(
            [str(ROOT / ".r-env/bin/Rscript"), "-e", "parse(file='scripts/120_run_track_b_lava.R')"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        source = (ROOT / "scripts/120_run_track_b_lava.R").read_text(encoding="utf-8")
        self.assertIn('!identical(pairs$pair_id, c("A", "B", "CONTROL"))', source)
        self.assertIn("run.pcor", source)
        self.assertIn("CONDITIONER_LOCAL_H2_INELIGIBLE", source)
        self.assertIn("run.pcor returned an incomplete partial estimate", source)
        self.assertIn("CONDITIONAL_UNSTABLE_MAX_R2", source)
        self.assertIn("A_SLEEP_APNEA_ONLY", source)
        self.assertIn('if ("h2.latent" %in% names(result))', source)
        self.assertIn("scripts/121_validate_track_b_lava.py", source)
        self.assertNotIn('"--input-info"', source)
        self.assertNotIn('"--checkpoint-dir"', source)

    def test_contract_uses_measured_per_locus_ram_not_an_assumed_gate(self) -> None:
        contract_source = (ROOT / "scripts/119_track_b_lava_contract.py").read_text(encoding="utf-8")
        policy = json.loads((ROOT / "config/track_b_local_analysis_policy.json").read_text(encoding="utf-8"))
        self.assertNotIn("validate_execution_compute", contract_source)
        self.assertNotIn("minimum_ram_bytes", contract_source)
        self.assertNotIn("minimum_ram_bytes", policy)
        self.assertEqual(policy["ram_aware_execution"]["maximum_loci_per_process"], 1)

    def test_validator_enforces_exact_families_and_conditional_gate(self) -> None:
        result = subprocess.run(
            ["python3", "-m", "py_compile", "scripts/119_track_b_lava_contract.py", "scripts/121_validate_track_b_lava.py"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        source = (ROOT / "scripts/121_validate_track_b_lava.py").read_text(encoding="utf-8")
        self.assertIn("univariate result is not the exact 2,495 x 8 family", source)
        self.assertIn("bivariate result differs from the locally eligible frozen three-pair family", source)
        self.assertIn("conditional result does not exactly cover FDR-supported A/B loci", source)

    def test_contract_does_not_expose_unvalidated_result_sealing(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/119_track_b_lava_contract.py", "--seal-results"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("TRACK_B_LAVA_RESULTS_SEALED", result.stdout + result.stderr)

    def test_conditional_interpretation_boundaries_and_unstable_status(self) -> None:
        validator = load_validator()
        self.assertEqual(
            validator.conditional_interpretation("TESTED", 0.4, 0.31, True, 0.25, 0.75)[0],
            "INDEPENDENT_LOCAL_COMPONENT",
        )
        self.assertEqual(
            validator.conditional_interpretation("TESTED", 0.4, 0.30, True, 0.25, 0.75)[0],
            "PARTIALLY_ATTENUATED",
        )
        self.assertEqual(
            validator.conditional_interpretation("TESTED", 0.4, 0.10, False, 0.25, 0.75)[0],
            "FULLY_ATTENUATED",
        )
        interpretation, _, qc = validator.conditional_interpretation(
            "CONDITIONAL_UNSTABLE_MAX_R2", 0.4, 0.2, False, 0.25, 0.75,
        )
        self.assertEqual((interpretation, qc), ("UNSTABLE", "CONDITIONAL_UNSTABLE_MAX_R2"))

    def test_failed_attempts_conservatively_expand_bh_family(self) -> None:
        validator = load_validator()
        self.assertEqual(validator.bh([0.01, 1.0]), [0.02, 1.0])

    def test_runner_is_two_pass_fresh_process_and_full_family_corrected(self) -> None:
        runner = (ROOT / "scripts/120_run_track_b_lava.R").read_text(encoding="utf-8")
        supervisor = (ROOT / "scripts/131_run_track_b_lava_sequential.py").read_text(encoding="utf-8")
        self.assertIn('"aggregate-discovery"', runner)
        self.assertIn('"conditional"', runner)
        self.assertIn("rep.int(1, nrow(univ))", runner)
        self.assertIn("rep.int(1, nrow(bivar))", runner)
        self.assertIn("p.adjust", runner)
        self.assertIn("gc(full = TRUE)", runner)
        self.assertIn("for index in range(1, 2496)", supervisor)
        self.assertIn("RUSAGE_CHILDREN", supervisor)
        self.assertIn("checkpoint ready link", supervisor)
        self.assertIn("RAM_BENCHMARK.provenance.json", supervisor)

    def test_staging_directory_requires_a_lowercase_sha256_identity(self) -> None:
        contract = load_validator().contract
        with tempfile.TemporaryDirectory() as directory:
            staging_root = Path(directory) / "staging"
            with patch.object(contract, "STAGING_ROOT", staging_root):
                fingerprint = "a" * 64
                self.assertEqual(contract.staging_dir(fingerprint), staging_root / f"lava_{fingerprint}")
                for invalid in ("a" * 63, "A" * 64, "g" * 64, "../" + "a" * 64):
                    with self.subTest(invalid=invalid), self.assertRaisesRegex(SystemExit, "invalid.*fingerprint"):
                        contract.staging_dir(invalid)

    def test_staged_artifact_identity_rejects_path_replacement_while_hashing(self) -> None:
        validator = load_validator()
        real_sha256 = validator.hashlib.sha256
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "staged.tsv"
            replacement = root / "replacement.tsv"
            source.write_bytes(b"original staged bytes\n")
            replacement.write_bytes(b"different and longer replacement bytes\n")

            class ReplacingDigest:
                def __init__(self) -> None:
                    self.inner = real_sha256()
                    self.replaced = False

                def update(self, block: bytes) -> None:
                    self.inner.update(block)
                    if not self.replaced:
                        replacement.replace(source)
                        self.replaced = True

                def hexdigest(self) -> str:
                    return self.inner.hexdigest()

            with (
                patch.object(validator, "ROOT", root),
                patch.object(validator.hashlib, "sha256", side_effect=ReplacingDigest),
                self.assertRaisesRegex(SystemExit, "changed while hashing"),
            ):
                validator.artifact_identity(source)

    def test_contract_hash_rejects_path_replacement_while_hashing(self) -> None:
        contract = load_validator().contract
        real_sha256 = contract.hashlib.sha256
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "canonical.tsv"
            replacement = root / "replacement.tsv"
            source.write_bytes(b"original canonical bytes\n")
            replacement.write_bytes(b"different and longer canonical replacement bytes\n")

            class ReplacingDigest:
                def __init__(self) -> None:
                    self.inner = real_sha256()
                    self.replaced = False

                def update(self, block: bytes) -> None:
                    self.inner.update(block)
                    if not self.replaced:
                        replacement.replace(source)
                        self.replaced = True

                def hexdigest(self) -> str:
                    return self.inner.hexdigest()

            with (
                patch.object(contract, "ROOT", root),
                patch.object(contract.hashlib, "sha256", side_effect=ReplacingDigest),
                self.assertRaisesRegex(SystemExit, "changed while hashing"),
            ):
                contract.sha256(source)

    def test_synthetic_result_family_validates_publishes_and_seals_without_replacement(self) -> None:
        validator = load_validator()
        contract = validator.contract
        fingerprint = "b" * 64
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            policy = synthetic_policy()
            policy_path = root / "config/policy.json"
            input_lock = root / "results/track_b/input.lock.json"
            reference = root / "ref/reference.json"
            chromosome_lock = root / "ref/chromosome-inputs.json"
            checkpoint_dir = root / "checkpoints"
            result_dir = root / "canonical/local"
            result_lock = result_dir / "lava_results.provenance.json"
            locus_file = root / "ref/loci.locfile"
            canonical = [
                result_dir / "lava_locus_status.tsv", result_dir / "lava_univariate.tsv",
                result_dir / "lava_bivariate.tsv", result_dir / "lava_conditional.tsv",
                root / "canonical/04_lava_local_results.tsv",
                root / "canonical/05_local_conditional_results.tsv",
            ]
            policy_path.parent.mkdir(parents=True, exist_ok=True)
            policy_path.write_text(json.dumps(policy) + "\n", encoding="utf-8")
            input_lock.parent.mkdir(parents=True, exist_ok=True)
            input_lock.write_text("synthetic input lock\n", encoding="utf-8")
            reference.parent.mkdir(parents=True, exist_ok=True)
            reference.write_text("synthetic reference provenance\n", encoding="utf-8")
            chromosome_lock.write_text("synthetic chromosome-input provenance\n", encoding="utf-8")
            fake_supervisor = write_contract_checkpoint_family(
                checkpoint_dir, policy["expected_loci"], fingerprint,
            )

            with (
                patch.object(validator, "ROOT", root),
                patch.multiple(
                    contract,
                    ROOT=root,
                    POLICY=policy_path,
                    INPUT_LOCK=input_lock,
                    REFERENCE=reference,
                    CHROMOSOME_INPUT_LOCK=chromosome_lock,
                    CHECKPOINT_DIR=checkpoint_dir,
                    RESULT_DIR=result_dir,
                    LOCUS_FILE=locus_file,
                    LOW_LEVEL_RESULTS=canonical[:4],
                    SUMMARY_RESULTS=canonical[4:],
                    RESULTS=canonical,
                    RESULT_LOCK=result_lock,
                    STAGING_ROOT=root / "staging",
                    load_policy=Mock(return_value=policy),
                    run_fingerprint=Mock(return_value=fingerprint),
                    load_supervisor=Mock(return_value=fake_supervisor),
                    validate_reference=Mock(return_value={}),
                    fully_verify_chromosome_inputs=Mock(return_value=None),
                ),
            ):
                built_policy, staged = build_synthetic_result_family(validator, root, fingerprint)
                self.assertEqual(built_policy, policy)
                self.assertEqual(staged, contract.staged_results(fingerprint))

                validated_identity = validator.validate(staged, build_summary_files=True)
                self.assertEqual(set(validated_identity), set(staged))
                local_summary = read_tsv(staged[4])
                conditional_summary = read_tsv(staged[5])
                self.assertEqual(len(local_summary), 3 * policy["expected_loci"])
                self.assertTrue(all(
                    row["QC"] == "INSUFFICIENT_LOCAL_H2" and row["interpretation_status"] == "UNDERPOWERED"
                    for row in local_summary if row["locus_id"] == "103"
                ))
                self.assertTrue(all(
                    row["QC"] == "LOCUS_PROCESS_FAILED"
                    for row in local_summary if row["locus_id"] == "102"
                ))
                self.assertEqual(
                    {row["QC"] for row in conditional_summary},
                    {"PASS", "CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2", "CONDITIONER_LOCAL_H2_INELIGIBLE"},
                )

                original_status = staged[0].read_bytes()
                staged[0].write_bytes(original_status + b"tampered after validation\n")
                with self.assertRaisesRegex(SystemExit, "changed after validation"):
                    validator.publish(staged, validated_identity)
                self.assertFalse(result_lock.exists())
                self.assertFalse(any(path.exists() for path in canonical))
                staged[0].write_bytes(original_status)

                canonical[2].parent.mkdir(parents=True, exist_ok=True)
                canonical[2].write_text("preexisting output must survive\n", encoding="utf-8")
                with self.assertRaisesRegex(SystemExit, "canonical output differs.*overwrite is forbidden"):
                    validator.publish(staged, validated_identity)
                self.assertEqual(canonical[2].read_text(encoding="utf-8"), "preexisting output must survive\n")
                self.assertFalse(result_lock.exists())
                self.assertEqual([path for path in canonical if path.exists()], canonical[:3])
                canonical[2].unlink()

                # Matching prefix files survive the interruption and are verified, not replaced.
                validator.publish(staged, validated_identity)
                self.assertTrue(result_lock.is_file())
                for source, destination in zip(staged, canonical, strict=True):
                    self.assertEqual(validator.artifact_identity(destination), validated_identity[source])
                result_payload = json.loads(result_lock.read_text(encoding="utf-8"))
                self.assertEqual(result_payload["run_fingerprint"], fingerprint)
                self.assertEqual(result_payload["checkpoint_locus_count"], policy["expected_loci"])
                self.assertEqual(result_payload["checkpoint_artifact_count"], policy["expected_loci"] * 2 + 2)
                self.assertEqual(
                    [(item["bytes"], item["sha256"]) for item in result_payload["results"]],
                    [validated_identity[source] for source in staged],
                )
                contract.validate_results()

    def test_zero_hit_family_publishes_header_only_bivariate_and_conditional_results(self) -> None:
        validator = load_validator()
        contract = validator.contract
        fingerprint = "c" * 64
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            policy = synthetic_policy()
            policy["expected_loci"] = 1
            policy["planned_univariate_tests"] = 8
            policy_path = root / "config/policy.json"
            input_lock = root / "results/track_b/input.lock.json"
            reference = root / "ref/reference.json"
            chromosome_lock = root / "ref/chromosome-inputs.json"
            checkpoint_dir = root / "checkpoints"
            result_dir = root / "canonical/local"
            result_lock = result_dir / "lava_results.provenance.json"
            locus_file = root / "ref/loci.locfile"
            canonical = [
                result_dir / "lava_locus_status.tsv", result_dir / "lava_univariate.tsv",
                result_dir / "lava_bivariate.tsv", result_dir / "lava_conditional.tsv",
                root / "canonical/04_lava_local_results.tsv",
                root / "canonical/05_local_conditional_results.tsv",
            ]
            policy_path.parent.mkdir(parents=True, exist_ok=True)
            policy_path.write_text(json.dumps(policy) + "\n", encoding="utf-8")
            input_lock.parent.mkdir(parents=True, exist_ok=True)
            input_lock.write_text("synthetic zero-hit input lock\n", encoding="utf-8")
            reference.parent.mkdir(parents=True, exist_ok=True)
            reference.write_text("synthetic zero-hit reference provenance\n", encoding="utf-8")
            chromosome_lock.write_text("synthetic zero-hit chromosome-input provenance\n", encoding="utf-8")
            fake_supervisor = write_contract_checkpoint_family(checkpoint_dir, 1, fingerprint)

            with (
                patch.object(validator, "ROOT", root),
                patch.multiple(
                    contract,
                    ROOT=root,
                    POLICY=policy_path,
                    INPUT_LOCK=input_lock,
                    REFERENCE=reference,
                    CHROMOSOME_INPUT_LOCK=chromosome_lock,
                    CHECKPOINT_DIR=checkpoint_dir,
                    RESULT_DIR=result_dir,
                    LOCUS_FILE=locus_file,
                    LOW_LEVEL_RESULTS=canonical[:4],
                    SUMMARY_RESULTS=canonical[4:],
                    RESULTS=canonical,
                    RESULT_LOCK=result_lock,
                    STAGING_ROOT=root / "staging",
                    load_policy=Mock(return_value=policy),
                    run_fingerprint=Mock(return_value=fingerprint),
                    load_supervisor=Mock(return_value=fake_supervisor),
                    validate_reference=Mock(return_value={}),
                    fully_verify_chromosome_inputs=Mock(return_value=None),
                ),
            ):
                built_policy, staged = build_zero_hit_result_family(validator, root, fingerprint)
                self.assertEqual(built_policy, policy)
                validated_identity = validator.validate(staged, build_summary_files=True)
                self.assertEqual(read_tsv(staged[2]), [])
                self.assertEqual(read_tsv(staged[3]), [])
                self.assertEqual(read_tsv(staged[5]), [])
                local_summary = read_tsv(staged[4])
                self.assertEqual(len(local_summary), 3)
                self.assertTrue(all(row["QC"] == "INSUFFICIENT_LOCAL_H2" for row in local_summary))

                validator.publish(staged, validated_identity)
                self.assertTrue(result_lock.is_file())
                self.assertEqual(read_tsv(canonical[2]), [])
                self.assertEqual(read_tsv(canonical[3]), [])
                self.assertEqual(read_tsv(canonical[5]), [])
                contract.validate_results()

    def test_summary_is_full_three_pair_locus_grid(self) -> None:
        validator = load_validator()
        policy = validator.contract.load_policy()
        loci, locus_index = validator.load_loci(policy)
        pairs = read_tsv(ROOT / "results/track_b/lava_pair_manifest.tsv")
        status = [{"LOC": locus, "status": "PROCESSED", "n_snps": "100"} for locus in loci]
        status_index = {row["LOC"]: row for row in status}
        univ_index = {
            (locus, trait): {"LOC": locus, "phen": trait, "h2.obs": "0.1", "analysis_status": "TESTED"}
            for locus in loci for trait in policy["trait_order"]
        }
        local = {
            "LOC": loci[0], "pair_id": "A", "analysis_status": "TESTED", "fdr_significant": "TRUE",
            "local_covariance": "0.02", "rho": "0.4", "rho.lower": "0.2", "rho.upper": "0.6",
            "p": "0.001", "p_fdr": "0.01",
        }
        context = {
            "policy": policy, "loci": loci, "locus_index": locus_index, "pairs": pairs,
            "bivar": [local], "conditional": [], "univ_index": univ_index,
            "bivar_index": {(loci[0], "A"): local}, "status_index": status_index,
        }
        local_rows, conditional_rows = validator.build_summaries(context)
        self.assertEqual(len(local_rows), 3 * 2495)
        self.assertEqual(conditional_rows, [])
        self.assertEqual(local_rows[0]["interpretation_status"], "STRONG_LOCAL_SHARING")
        self.assertEqual(local_rows[2495]["interpretation_status"], "UNDERPOWERED")
        self.assertEqual(local_rows[0]["local_h2_scale"], "OBSERVED")
        self.assertIn("CI_EQUIVALENT", local_rows[0]["SE_method"])

    def test_publisher_refuses_preexisting_unsealed_destination(self) -> None:
        validator = load_validator()
        with tempfile.TemporaryDirectory(dir=ROOT / "work") as directory:
            root = Path(directory)
            staged = [root / "staged.tsv"]
            staged[0].write_text("x\n", encoding="utf-8")
            destination = root / "canonical.tsv"
            destination.write_text("preserve\n", encoding="utf-8")
            result_lock = root / "seal.json"
            original_results = validator.contract.RESULTS
            original_lock = validator.contract.RESULT_LOCK
            try:
                validator.contract.RESULTS = [destination]
                validator.contract.RESULT_LOCK = result_lock
                with self.assertRaises(SystemExit):
                    validator.publish(staged)
                self.assertEqual(destination.read_text(encoding="utf-8"), "preserve\n")
                self.assertFalse(result_lock.exists())
            finally:
                validator.contract.RESULTS = original_results
                validator.contract.RESULT_LOCK = original_lock


if __name__ == "__main__":
    unittest.main()
