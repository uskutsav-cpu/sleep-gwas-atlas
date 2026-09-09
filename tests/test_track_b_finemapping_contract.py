import csv
import hashlib
import importlib.util
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "config/track_b_finemapping_policy.json"
TRACK_B = ROOT / "results/track_b"
FINEMAPPING = TRACK_B / "finemapping"


def load_continuation():
    path = ROOT / "scripts/145_build_track_b_finemapping_continuation_gate.py"
    specification = importlib.util.spec_from_file_location(
        "track_b_finemapping_continuation_for_historical_test", path,
    )
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


CONTINUATION = load_continuation()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBFineMappingContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.policy = json.loads(POLICY.read_text(encoding="utf-8"))
        cls.lock = json.loads(
            (FINEMAPPING / "contract.lock.json").read_text(encoding="utf-8")
        )

    def test_historical_lock_stays_immutable_and_continuation_discloses_drift(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/130_track_b_finemapping_contract.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("fine-mapping contract drifted", result.stdout + result.stderr)
        current_policy_hash = hashlib.sha256(POLICY.read_bytes()).hexdigest()
        current_script_hash = hashlib.sha256(
            (ROOT / "scripts/130_track_b_finemapping_contract.py").read_bytes()
        ).hexdigest()
        self.assertNotEqual(self.lock["policy_sha256"], current_policy_hash)
        self.assertNotEqual(self.lock["script_sha256"], current_script_hash)
        snapshot = CONTINUATION.validate_base_contract_snapshot()
        lineage = snapshot["lineage"]
        self.assertTrue(lineage["historical_lock_preserved"])
        self.assertFalse(lineage["policy_identity_matches_historical_lock"])
        self.assertFalse(lineage["contract_script_identity_matches_historical_lock"])
        self.assertIn("HISTORICAL_LOCK_NOT_REWRITTEN", lineage["mismatch_handling"])

    def test_locus_entry_is_the_nonduplicated_method_union(self) -> None:
        locus = self.policy["locus_entry_contract"]
        self.assertFalse(locus["intersection_required"])
        self.assertTrue(locus["consensus_only_selection_forbidden"])
        self.assertTrue(locus["result_ranked_maximum_forbidden"])
        self.assertTrue(locus["global_rg_filter_forbidden"])
        self.assertIn("nonduplicated union", locus["family_definition"])
        self.assertIn("same pair and block", locus["deduplication_rule"])
        self.assertIn("Multiple eligible signals", locus["multi_signal_rule"])

        upstream = self.policy["future_upstream_results"]
        self.assertIn("FDR<=0.05", upstream["local_sharing"]["eligible_row_rule"])
        self.assertIn("do not require pleiotropy", upstream["local_sharing"]["eligible_row_rule"])
        self.assertEqual(
            upstream["local_sharing"]["required_schema"][-3:],
            ["SNP_count", "QC", "interpretation_status"],
        )
        pleiotropy = upstream["pleiotropy"]
        self.assertEqual(pleiotropy["method_union_loci"], "results/track_b/08_shared_loci.tsv")
        self.assertEqual(
            pleiotropy["required_method_labels"],
            ["PLACO_AND_CONJFDR", "PLACO_ONLY", "CONJFDR_ONLY"],
        )
        self.assertIn("blocked or not tested", pleiotropy["blocked_method_label_rule"])
        self.assertEqual(
            pleiotropy["method_union_required_schema_prefix"],
            [
                "pair_id", "lead_variant", "chr", "position", "PLACO_P", "FDR",
                "trait1_P", "trait2_P", "locus_start", "locus_end",
                "independent_signal", "annotations",
            ],
        )

    def test_primary_and_control_roles_are_frozen(self) -> None:
        scope = self.policy["pair_scope"]
        self.assertEqual(scope["pair_order"], ["A", "B", "CONTROL"])
        self.assertEqual(scope["primary_pairs"], ["A", "B"])
        self.assertEqual(scope["control_pairs"], ["CONTROL"])
        by_pair = {row["pair_id"]: row for row in scope["pairs"]}
        self.assertEqual(by_pair["A"]["family_role"], "PRIMARY_DISCOVERY")
        self.assertEqual(by_pair["B"]["family_role"], "PRIMARY_DISCOVERY")
        self.assertEqual(
            by_pair["CONTROL"]["family_role"], "POSITIVE_CONTROL_NON_NOVELTY"
        )
        self.assertIn("excluded from novelty", scope["control_policy"])

    def test_full_dense_inputs_and_pinned_engine_are_bound(self) -> None:
        immutable = self.lock["immutable_input_sha256"]
        for path in (
            "results/track_b/pair_manifest.tsv",
            "results/track_b/pair_manifest.lock.json",
            "results/track_b/03_dense_input_qc.tsv",
            "results/track_b/03_dense_input_qc.lock.json",
            "config/fine_mapping_sources.tsv",
            "discovery_extension/config/fine_mapping_method_references.tsv",
            "discovery_extension/scripts/35_run_susie_coloc.R",
            "scripts/fine_mapping_contract.py",
            "results/track_b/pleiotropy/contract.lock.json",
            "results/track_b/pleiotropy/input_gate.lock.json",
            ".r-env/share/finemapping/source_archives/susieR_0.14.2.tar.gz",
            ".r-env/share/finemapping/source_archives/coloc_v5.2.3.tar.gz",
        ):
            self.assertIn(path, immutable)
        # These dependencies were added after the historical lock. They are
        # current pre-result policy inputs and are bound by the additive 145
        # continuation instead of rewriting the old evidence.
        for path in (
            "config/track_b_local_analysis_policy.json",
            "results/track_b/local_analysis_input.lock.json",
        ):
            self.assertNotIn(path, immutable)
            self.assertIn(path, self.policy["immutable_pre_result_inputs"].values())
        continuation = CONTINUATION.validate_base_contract_snapshot()
        self.assertIn("local_analysis_policy", continuation["identities"])
        self.assertIn("local_analysis_input_lock", continuation["identities"])
        self.assertNotIn(
            "results/track_b/pleiotropy/contract.lock.json",
            self.lock["future_required_paths_not_hashed"],
        )
        self.assertNotIn(
            "results/track_b/pleiotropy/input_gate.lock.json",
            self.lock["future_required_paths_not_hashed"],
        )
        self.assertEqual(len(self.lock["dense_input_sha256"]), 5)
        self.assertEqual(
            set(self.lock["dense_input_sha256"]),
            {
                "data/harmonized/snoring.harmonized.tsv.gz",
                "data/harmonized/parental_lifespan.harmonized.tsv.gz",
                "data/harmonized/insomnia.harmonized.tsv.gz",
                "data/harmonized/adhd.harmonized.tsv.gz",
                "data/harmonized/frailty.harmonized.tsv.gz",
            },
        )
        dense = self.policy["dense_input_contract"]
        self.assertTrue(dense["hapmap3_only_forbidden"])
        self.assertEqual(dense["minimum_variant_rows_per_trait"], 5_000_000)
        self.assertIsNone(dense["maximum_locus_variants"])
        self.assertIn("no count cap", dense["locus_size_rule"])
        software = self.policy["software_contract"]
        self.assertEqual(software["susieR"]["version"], "0.14.2")
        self.assertEqual(software["coloc"]["version"], "5.2.3")
        self.assertIn(
            "generic atlas consensus-only",
            software["shared_python_contract_use_rule"],
        )

    def test_signed_ld_prior_grid_and_diagnostics_are_exact(self) -> None:
        ld = self.policy["signed_ld_contract"]
        self.assertEqual(ld["ancestry"], "EUR")
        self.assertEqual(ld["build"], "GRCh37/hg19")
        self.assertTrue(ld["exact_variant_order_sha256_required"])
        self.assertIn("unsigned r2", ld["allele_rule"])

        coloc = self.policy["colocalization_contract"]
        self.assertEqual(coloc["p1"], 1e-4)
        self.assertEqual(coloc["p2"], 1e-4)
        self.assertEqual(coloc["p12_primary"], 1e-5)
        self.assertEqual(coloc["p12_sensitivity_grid"], [1e-6, 5e-6, 1e-5, 5e-5])
        self.assertIn("same signal1/signal2 pair", coloc["prior_robust_rule"])

        diagnostics = self.policy["diagnostic_contract"]
        self.assertIn("RSS-LD inconsistency s", diagnostics["required_checks"])
        self.assertIn("kriging allele-switch outliers", diagnostics["required_checks"])
        self.assertIn("LD positive semidefiniteness", diagnostics["required_checks"])
        self.assertIn("FAILED_LD_QC", diagnostics["allowed_statuses"])
        self.assertIn("FAILED_MODEL_CONVERGENCE", diagnostics["allowed_statuses"])

        resources = self.policy["ram_aware_execution_contract"]
        self.assertEqual(
            resources["execution_unit"], "ONE_VALIDATED_PAIR_LOCUS_PER_FRESH_OS_PROCESS"
        )
        self.assertEqual(resources["maximum_concurrent_loci_per_worker"], 1)
        self.assertTrue(resources["full_locus_variant_universe_required"])
        self.assertTrue(resources["full_ancestry_matched_signed_ld_required"])
        self.assertTrue(resources["locus_splitting_forbidden"])
        self.assertTrue(resources["variant_thinning_for_compute_forbidden"])
        estimate = resources["admission_estimate"]
        self.assertEqual(estimate["ld_scalar_bytes"], 8)
        self.assertEqual(estimate["ld_live_copy_multiplier"], 8)
        self.assertEqual(estimate["absolute_minimum_memory_bytes"], 4 * 1024**3)
        self.assertIn("larger host", resources["insufficient_memory_rule"])
        self.assertIn("exclusive no-replace", resources["atomic_resume"]["staging_rule"])
        self.assertIn("never overwritten", resources["atomic_resume"]["resume_rule"])
        self.assertTrue({
            "elapsed_monotonic_seconds", "peak_rss_bytes", "rss_measurement_backend",
            "task_sha256", "variant_order_sha256", "signed_ld_sha256",
            "output_manifest_sha256",
        }.issubset(resources["resource_metrics_schema"]))
        self.assertEqual(
            self.lock["execution_unit"],
            "ONE_VALIDATED_PAIR_LOCUS_PER_FRESH_OS_PROCESS",
        )
        self.assertEqual(self.lock["maximum_concurrent_loci_per_worker"], 1)
        self.assertEqual(self.lock["locus_splitting"], "FORBIDDEN")

    def test_exact_10_11_12_schemas_and_classifications(self) -> None:
        outputs = self.policy["required_science_outputs"]
        self.assertEqual(
            [outputs[key]["path"] for key in ("output_10", "output_11", "output_12")],
            [
                "results/track_b/10_finemap_trait1.tsv",
                "results/track_b/11_finemap_trait2.tsv",
                "results/track_b/12_trait_trait_coloc.tsv",
            ],
        )
        self.assertEqual(outputs["output_10"]["schema"], outputs["output_11"]["schema"])
        self.assertEqual(len(outputs["output_10"]["schema"]), len(set(outputs["output_10"]["schema"])))
        self.assertEqual(len(outputs["output_12"]["schema"]), len(set(outputs["output_12"]["schema"])))
        self.assertTrue({
            "PIP", "credible_set_ids", "credible_set_sizes", "credible_set_min_abs_corr",
            "rss_ld_s", "kriging_allele_switch_outlier", "fine_mapping_classifications",
        }.issubset(outputs["output_10"]["schema"]))
        self.assertTrue({
            "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4",
            "same_signal_pair_prior_robust", "classification",
        }.issubset(outputs["output_12"]["schema"]))

        classifications = self.policy["classification_contract"]
        self.assertEqual(
            classifications["fine_mapping_labels_in_order"],
            [
                "HIGH_PIP_VARIANT", "SMALL_CREDIBLE_SET", "DIFFUSE_SIGNAL",
                "LD_UNCERTAINTY", "MODEL_INSTABILITY",
            ],
        )
        self.assertEqual(
            classifications["trait_coloc_allowed"],
            [
                "STRONG_SHARED_SIGNAL", "MODERATE_SHARED_SIGNAL", "DISTINCT_SIGNALS",
                "INCONCLUSIVE", "INVALID_INPUT",
            ],
        )

    def test_current_readiness_is_blocked_without_science_placeholders(self) -> None:
        rows = read_tsv(FINEMAPPING / "readiness.tsv")
        self.assertEqual([row["gate_order"] for row in rows], [str(i) for i in range(1, 11)])
        self.assertEqual(rows[-1]["gate_id"], "OVERALL")
        self.assertEqual(rows[-1]["status"], "BLOCKED_UPSTREAM")
        self.assertTrue(all(row["result_accessed"] == "FALSE" for row in rows))
        self.assertTrue(all("science" not in row["status"].lower() for row in rows))
        self.assertEqual(
            self.lock["artifact_role"],
            "PRE_RESULT_FINE_MAPPING_AND_TRAIT_COLOC_CONTRACT_NOT_SCIENTIFIC_RESULT",
        )
        self.assertFalse(self.lock["results_accessed_before_contract_freeze"])
        self.assertFalse(self.lock["science_outputs_created_by_contract"])
        resource_gate = next(row for row in rows if row["gate_id"] == "FRESH_PROCESS_RESOURCE_EXECUTION")
        self.assertEqual(resource_gate["status"], "DEFERRED_UPSTREAM")
        self.assertIn("NO_SPLITTING", resource_gate["observed_state"])
        for key in ("output_10", "output_11", "output_12"):
            self.assertFalse((ROOT / self.policy["required_science_outputs"][key]["path"]).exists())

    def test_zero_family_and_publication_semantics_are_fail_closed(self) -> None:
        locus = self.policy["locus_entry_contract"]
        self.assertEqual(
            locus["zero_family_terminal_status"],
            "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY",
        )
        self.assertIn("both complete sealed upstream families", locus["zero_family_rule"])
        self.assertIn("do not create header-only 10/11/12", locus["zero_family_rule"])
        execution = self.policy["execution_and_publication"]
        self.assertTrue(execution["canonical_output_overwrite_forbidden"])
        self.assertTrue(execution["staged_validation_required"])
        self.assertTrue(execution["exclusive_publication_required"])
        self.assertTrue(execution["generic_atlas_science_output_substitution_forbidden"])
        self.assertIn("header-only publication is forbidden", execution["publication_rule"])
        self.assertFalse(self.lock["zero_family_creates_10_11_12"])
        statuses = self.policy["result_status_semantics"]
        self.assertEqual(statuses["fine_mapping_allowed"], ["COMPLETE", "FAILED_QC"])
        self.assertEqual(statuses["trait_coloc_allowed"], ["TESTED", "INVALID_INPUT"])
        self.assertIn("Every locked pair-locus", statuses["complete_family_rule"])


if __name__ == "__main__":
    unittest.main()
