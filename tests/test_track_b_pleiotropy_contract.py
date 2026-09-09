import csv
import hashlib
import importlib.util
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "config/track_b_pleiotropy_policy.json"
GATE_DIR = ROOT / "results/track_b/pleiotropy"


def load_script(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, ROOT / "scripts" / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONTRACT = load_script("123_track_b_pleiotropy_contract.py", "track_b_pleiotropy_contract_test")
INPUT_GATE = load_script(
    "124_build_track_b_pleiotropy_input_gate.py", "track_b_pleiotropy_input_gate_test",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBPleiotropyContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.policy = json.loads(POLICY.read_text(encoding="utf-8"))

    def test_contract_and_full_input_gate_verify(self) -> None:
        contract_result = subprocess.run(
            ["python3", "scripts/123_track_b_pleiotropy_contract.py", "--verify-contract"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(contract_result.returncode, 0, contract_result.stdout + contract_result.stderr)
        gate_result = subprocess.run(
            ["python3", "scripts/124_build_track_b_pleiotropy_input_gate.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        if gate_result.returncode == 0:
            self.assertIn("full_gzip_scans=5", gate_result.stdout)
            return

        # The no-replace V1 readiness snapshot legitimately retains the
        # transient "active LAVA download" blocker that was true at freeze
        # time.  Once the download finishes, live regeneration must differ;
        # validate that this is the *only* drift and that the frozen bytes are
        # still exactly hash-bound rather than weakening or rewriting history.
        self.assertIn(
            "readiness_gate.tsv", gate_result.stdout + gate_result.stderr,
        )
        expected = INPUT_GATE.build(ROOT)
        frozen_input = (ROOT / INPUT_GATE.INPUT_GATE).read_text(encoding="utf-8")
        frozen_readiness = (ROOT / INPUT_GATE.READINESS_GATE).read_text(encoding="utf-8")
        frozen_lock = json.loads((ROOT / INPUT_GATE.INPUT_LOCK).read_text(encoding="utf-8"))
        self.assertEqual(frozen_input, expected[INPUT_GATE.INPUT_GATE])
        transient = ";ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION"
        self.assertEqual(
            frozen_readiness.replace(transient, ""),
            expected[INPUT_GATE.READINESS_GATE],
        )
        self.assertEqual(
            frozen_lock["input_gate_sha256"],
            hashlib.sha256(frozen_input.encode("utf-8")).hexdigest(),
        )
        self.assertEqual(
            frozen_lock["readiness_gate_sha256"],
            hashlib.sha256(frozen_readiness.encode("utf-8")).hexdigest(),
        )

    def test_primary_family_and_control_are_exact_and_disjoint(self) -> None:
        self.assertEqual(self.policy["primary_pair_family"], ["A", "B"])
        self.assertEqual(self.policy["control_pair_family"], ["CONTROL"])
        self.assertTrue(
            set(self.policy["primary_pair_family"]).isdisjoint(self.policy["control_pair_family"])
        )
        observed = [
            (row["pair_id"], row["trait1"], row["trait2"], row["family_role"])
            for row in self.policy["pairs"]
        ]
        self.assertEqual(observed, [
            ("A", "snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
            ("B", "insomnia", "adhd", "PRIMARY_DISCOVERY"),
            ("CONTROL", "insomnia", "frailty", "POSITIVE_CONTROL"),
        ])

    def test_exact_placo_and_conjfdr_threshold_families(self) -> None:
        placo = self.policy["placo_plus"]
        self.assertEqual(placo["primary_pair_family_headline_threshold"], 2.5e-8)
        self.assertEqual(placo["per_pair_genome_wide_threshold"], 5e-8)
        self.assertEqual(placo["control_genome_wide_threshold"], 5e-8)
        self.assertEqual(placo["within_pair_bh_alpha"], 0.05)
        self.assertEqual(
            2 * placo["primary_pair_family_headline_threshold"],
            placo["per_pair_genome_wide_threshold"],
        )
        self.assertTrue(placo["full_p_ledger_required"])
        self.assertIn("every aligned eligible genome-wide variant", placo["within_pair_bh_family"])
        self.assertIn("P=1", placo["within_pair_bh_family"])
        self.assertEqual(self.policy["conjfdr"]["threshold"], 0.05)
        self.assertTrue(self.policy["conjfdr"]["independent_inclusion"])
        self.assertFalse(self.policy["evidence_union"]["intersection_required"])

    def test_pair_materialization_plausibility_floors_are_machine_checkable(self) -> None:
        plausibility = self.policy["dense_input_contract"]["pair_materialization_plausibility"]
        self.assertEqual(plausibility["minimum_aligned_eligible_variants"], 1_000_000)
        self.assertEqual(plausibility["minimum_fraction_of_smaller_dense_input"], 0.8)
        self.assertEqual(plausibility["required_autosomes"], list(range(1, 23)))
        self.assertTrue(plausibility["nonzero_eligible_variants_on_every_required_autosome"])
        self.assertTrue(plausibility["per_autosome_provenance_required"])
        self.assertEqual(plausibility["per_autosome_provenance_fields"], [
            "CHR", "trait1_rows", "trait2_rows", "exact_rsid_coordinate_matches",
            "allele_matches", "eligible_written",
        ])
        self.assertIn("FAILED_QC", plausibility["failure_rule"])

    def test_canonical_07_08_09_and_control_contract(self) -> None:
        paths = self.policy["future_result_paths"]
        self.assertEqual(paths["placo_plus_variants_07"], "results/track_b/07_placo_plus_variants.tsv")
        self.assertEqual(paths["shared_loci_08"], "results/track_b/08_shared_loci.tsv")
        self.assertEqual(paths["pleiotropy_comparison_09"], "results/track_b/09_pleiotropy_comparison.tsv")
        self.assertEqual(paths["control_directory"], "results/track_b/control_insomnia_frailty/")
        self.assertEqual(
            self.policy["required_future_result_schemas"]["shared_loci_08_required_prefix"],
            [
                "pair_id", "lead_variant", "chr", "position", "PLACO_P", "FDR",
                "trait1_P", "trait2_P", "locus_start", "locus_end",
                "independent_signal", "annotations",
            ],
        )

    def test_comparison_labels_do_not_misclassify_blocked_methods(self) -> None:
        evidence = self.policy["evidence_union"]
        self.assertEqual(
            evidence["comparison_labels_in_order"],
            ["PLACO_AND_CONJFDR", "PLACO_ONLY", "CONJFDR_ONLY"],
        )
        self.assertIn("only when both methods completed", evidence["comparison_only_label_rule"])
        self.assertIn("do not emit an ONLY label", evidence["comparison_only_label_rule"])
        states = self.policy["method_state_semantics"]
        self.assertIn("not P=1", states["blocked_rule"])
        self.assertIn("Never publish an empty scientific result table", states["partial_completion_rule"])
        self.assertIn("complete checksum-sealed genome-wide scan", states["zero_hit_rule"])

    def test_full_eur_hg19_ld_and_signal_rules_are_frozen(self) -> None:
        ld = self.policy["ld_reference_policy"]
        self.assertEqual((ld["ancestry"], ld["build"], ld["sample_count"]), ("EUR", "GRCh37/hg19", 503))
        self.assertEqual(ld["source_archive_expected_bytes"], 511_626_945)
        self.assertEqual(
            ld["source_archive_sha256"],
            "83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954",
        )
        self.assertEqual(ld["materialized_directory"], "ref/track_b/pleiotropy/g1000_eur")
        self.assertEqual(
            ld["materialized_prefix"], "ref/track_b/pleiotropy/g1000_eur/g1000_eur",
        )
        self.assertEqual(ld["reference_chromosomes"], list(range(1, 24)))
        self.assertEqual(ld["analysis_autosomes"], list(range(1, 23)))
        self.assertEqual(ld["nonanalysis_chromosomes_retained_for_source_fidelity"], [23])
        self.assertEqual(ld["variant_count"], 22_665_064)
        self.assertEqual(ld["expected_autosomal_variant_count"], 22_132_657)
        self.assertEqual(sum(ld["expected_reference_chromosome_variant_counts"].values()), 22_665_064)
        self.assertEqual(
            sum(
                ld["expected_reference_chromosome_variant_counts"][str(chromosome)]
                for chromosome in range(1, 23)
            ),
            22_132_657,
        )
        self.assertEqual(ld["expected_reference_chromosome_variant_counts"]["23"], 532_407)
        self.assertEqual(ld["source_members"]["bed"]["expected_bytes"], 2_855_798_067)
        self.assertEqual(ld["source_members"]["bim"]["expected_bytes"], 658_724_714)
        self.assertEqual(ld["source_members"]["fam"]["expected_bytes"], 12_575)
        self.assertEqual(ld["source_members"]["fam"]["records"], 503)
        self.assertEqual(ld["staging_safety_bytes"], 1_073_741_824)
        self.assertEqual(ld["minimum_free_storage_bytes"], 4_588_277_180)
        self.assertIn("do not derive or relabel an autosome-only BED/BIM/FAM", ld["raw_reference_fidelity_rule"])
        self.assertIn("chr23/X is always forbidden", ld["analysis_chromosome_rule"])
        self.assertIn("Byte-identical extraction", ld["materialization_transform_rule"])
        self.assertIn("no LD-block splitting", ld["materialization_transform_rule"])
        self.assertIn("SNP reduction", ld["materialization_transform_rule"])
        self.assertTrue(ld["full_genome_wide_reference_required"])
        self.assertTrue(ld["hapmap3_subset_forbidden"])
        self.assertEqual(ld["clump_r2_strict_upper_bound"], 0.1)
        self.assertEqual(ld["clump_window_kb"], 1000)
        self.assertIn("r2<0.1 within 1 Mb", ld["lead_independence_rule"])
        self.assertIn("r2>=0.6", ld["cross_method_locus_rule"])
        self.assertIn("distance-only", ld["cross_method_locus_rule"])

    def test_gate_records_full_live_hash_count_and_schema(self) -> None:
        rows = read_tsv(GATE_DIR / "input_gate.tsv")
        self.assertEqual([row["pair_id"] for row in rows], ["A", "B", "CONTROL"])
        self.assertTrue(all(row["input_gate_status"] == "READY_DENSE_INPUTS" for row in rows))
        self.assertTrue(all(row["gzip_integrity"] == "FULL_DECOMPRESSION_CRC_PASS_BOTH" for row in rows))
        expected_schema = "SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N"
        for row in rows:
            self.assertEqual(row["ancestry"], "EUR")
            self.assertEqual(row["analysis_build"], "hg19")
            self.assertEqual(row["trait1_schema"], expected_schema)
            self.assertEqual(row["trait2_schema"], expected_schema)
            self.assertGreaterEqual(int(row["trait1_rows"]), 5_000_000)
            self.assertGreaterEqual(int(row["trait2_rows"]), 5_000_000)
            self.assertRegex(row["trait1_sha256"], r"^[0-9a-f]{64}$")
            self.assertRegex(row["trait2_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(row["minimum_aligned_eligible_variants"], "1000000")
            self.assertEqual(row["minimum_fraction_of_smaller_dense_input"], "0.8")
            self.assertEqual(row["required_autosomes"], ",".join(str(value) for value in range(1, 23)))
            self.assertEqual(row["per_autosome_provenance_required"], "TRUE")
            self.assertEqual(row["pair_alignment_status"], "NOT_RUN_PRE_RESULT;EXACT_GLOBAL_AND_22_AUTOSOME_COUNTS_REQUIRED_AT_MATERIALIZATION")
        lock = json.loads((GATE_DIR / "input_gate.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(lock["full_gzip_integrity_scans_completed"], 5)
        self.assertEqual(lock["unique_dense_trait_count"], 5)
        self.assertTrue(lock["all_pair_dense_input_gates_pass"])
        self.assertEqual(lock["scientific_result_count"], 0)
        for identity in lock["live_dense_inputs"].values():
            self.assertGreaterEqual(identity["rows"], 5_000_000)
            self.assertEqual(identity["schema"], expected_schema)
            self.assertEqual(identity["gzip_integrity"], "FULL_DECOMPRESSION_CRC_PASS")

    def test_readiness_is_not_presented_as_science(self) -> None:
        rows = read_tsv(GATE_DIR / "readiness_gate.tsv")
        self.assertEqual([row["component_id"] for row in rows], [
            "PLACO_PLUS", "CONJFDR", "METHOD_UNION_08_09",
        ])
        self.assertTrue(all(row["overall_status"] == "BLOCKED_BY_MULTIPLE_GATES" for row in rows))
        self.assertTrue(all(row["claim_status"] == "NO_SCIENTIFIC_RESULT" for row in rows))
        self.assertTrue(all("NOT_TESTED_METHOD_BLOCKED" in row["blocked_result_semantics"] for row in rows))
        self.assertTrue(all("NO_ONLY_LABEL_IF_OTHER_METHOD_BLOCKED" in row["blocked_result_semantics"] for row in rows))
        for path in (
            ROOT / "results/track_b/07_placo_plus_variants.tsv",
            ROOT / "results/track_b/08_shared_loci.tsv",
            ROOT / "results/track_b/09_pleiotropy_comparison.tsv",
        ):
            self.assertFalse(path.exists(), f"pre-result contract must not manufacture {path}")

    def test_contract_lock_binds_existing_pair_and_dense_locks(self) -> None:
        lock = json.loads((GATE_DIR / "contract.lock.json").read_text(encoding="utf-8"))
        self.assertFalse(lock["pleiotropy_results_accessed_before_contract_freeze"])
        self.assertEqual(lock["primary_pair_family"], ["A", "B"])
        self.assertEqual(lock["control_pair_family"], ["CONTROL"])
        self.assertTrue(lock["method_union_not_intersection"])
        self.assertTrue(lock["generic_396_pair_substitution_forbidden"])
        self.assertEqual(
            lock["script_sha256"]["scripts/127_prepare_track_b_pleiotropy_ld.py"],
            CONTRACT.sha256(ROOT / "scripts/127_prepare_track_b_pleiotropy_ld.py"),
        )
        self.assertEqual(
            lock["upstream"]["pair_manifest_sha256"],
            CONTRACT.sha256(ROOT / "results/track_b/pair_manifest.tsv"),
        )
        self.assertEqual(
            lock["upstream"]["dense_qc_sha256"],
            CONTRACT.sha256(ROOT / "results/track_b/03_dense_input_qc.tsv"),
        )
        input_lock = json.loads((GATE_DIR / "input_gate.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(
            input_lock["ld_materializer_script_sha256"],
            lock["script_sha256"]["scripts/127_prepare_track_b_pleiotropy_ld.py"],
        )


if __name__ == "__main__":
    unittest.main()
