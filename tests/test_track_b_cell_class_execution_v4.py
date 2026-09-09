#!/usr/bin/env python3
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E = load("track_b_cell_class_enrichment_v4_tested", "scripts/164_run_track_b_cell_class_enrichment.py")
C = load("track_b_cell_class_replication_v4_tested", "scripts/165_classify_track_b_cell_replication.py")
B = load("track_b_cell_class_contract_v4_tested", "scripts/163_build_track_b_cell_class_contract.py")


def write(root: Path, relative: str, payload: bytes) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return path


def identity_block(root: Path, relative: str, rows: int) -> dict[str, object]:
    payload = (root / relative).read_bytes()
    return {
        "path": relative, "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(),
        "rows": rows, "fields": E.ELIGIBLE_FIELDS,
    }


def gate_fixture(
    root: Path,
    *,
    ancestry: str = "EUR",
    selected_indexes: list[int] | None = None,
) -> tuple[dict[str, object], str]:
    pair_rel = "results/track_b/pair_manifest.tsv"
    panel_rel = "config/analysis_panel.tsv"
    write(root, pair_rel, b"pair_id\tsleep_trait\texternal_trait\nA\tSleepTrait\tExternalTrait\n")
    write(
        root, panel_rel,
        f"trait_id\tsource_id\tancestry\tsource_status\nSleepTrait\tSRC_SLEEP\t{ancestry}\tSOURCE_VERIFIED\n".encode(),
    )
    candidates = [
        {"pair_id": "A", "target_trait_id": "SleepTrait", "locus_id": "chr1:1-2", "gene_id": "ENSG00000000001", "genetic_evidence_tier": "G4", "brain_tissue_gate": "PASS", "eligibility_status": "ELIGIBLE"},
        {"pair_id": "A", "target_trait_id": "SleepTrait", "locus_id": "chr2:1-2", "gene_id": "ENSG00000000002", "genetic_evidence_tier": "G5", "brain_tissue_gate": "PASS", "eligibility_status": "ELIGIBLE"},
        {"pair_id": "A", "target_trait_id": "SleepTrait", "locus_id": "chr3:1-2", "gene_id": "ENSG00000000003", "genetic_evidence_tier": "G3", "brain_tissue_gate": "PASS", "eligibility_status": "INELIGIBLE"},
        {"pair_id": "A", "target_trait_id": "SleepTrait", "locus_id": "chr4:1-2", "gene_id": "NA", "genetic_evidence_tier": "FAILED", "brain_tissue_gate": "FAILED", "eligibility_status": "FAILED_UPSTREAM"},
    ]
    selected = [candidates[index] for index in (selected_indexes if selected_indexes is not None else [0, 1])]
    candidate_rel = "results/track_b/upstream/complete_candidates.tsv"
    eligible_rel = "results/track_b/upstream/eligible.tsv"
    source_rel = "results/track_b/upstream/source_contract.json"
    gate_rel = "results/track_b/upstream/gate.json"
    write(root, candidate_rel, E.tsv_bytes(E.ELIGIBLE_FIELDS, candidates))
    write(root, eligible_rel, E.tsv_bytes(E.ELIGIBLE_FIELDS, selected))
    candidate_block = identity_block(root, candidate_rel, len(candidates))
    eligible_block = identity_block(root, eligible_rel, len(selected))
    source = {
        "schema_version": E.UPSTREAM_SOURCE_SCHEMA, "contract_kind": E.UPSTREAM_SOURCE_KIND,
        "terminal_status": E.UPSTREAM_SOURCE_STATUS, "pair_id": "A",
        "target_trait_id": "SleepTrait", "target_trait_role": "sleep",
        "selection_rule": "PREDECLARED_G4_G5_AND_BRAIN_TISSUE_PASS_ONLY",
        "complete_candidate_family": candidate_block, "eligible_family": eligible_block,
        "failed_units_preserved": True,
    }
    write(root, source_rel, E.canonical_json(source))
    source_payload = (root / source_rel).read_bytes()
    gate = {
        "schema_version": E.GATE_SCHEMA, "gate_kind": "REAL_UPSTREAM_ELIGIBLE_LOCUS_GENE_FAMILY",
        "terminal_status": E.GATE_STATUS, "pair_id": "A", "target_trait_id": "SleepTrait",
        "target_trait_role": "sleep", "brain_tissue_gate": "PASS",
        "eligible_genetic_tiers": ["G4", "G5"],
        "selection_rule": "PREDECLARED_G4_G5_AND_BRAIN_TISSUE_PASS_ONLY",
        "complete_candidate_family": candidate_block, "eligible_family": eligible_block,
        "upstream_contract": {
            "path": source_rel, "bytes": len(source_payload),
            "sha256": hashlib.sha256(source_payload).hexdigest(),
            "schema_version": E.UPSTREAM_SOURCE_SCHEMA, "contract_kind": E.UPSTREAM_SOURCE_KIND,
            "terminal_status": E.UPSTREAM_SOURCE_STATUS,
        },
        "failed_units_preserved": True,
    }
    write(root, gate_rel, E.canonical_json(gate))
    policy = {
        "upstream_gate": {
            "pair_manifest_path": pair_rel, "analysis_panel_path": panel_rel,
            "maximum_family_bytes": 67108864,
        }
    }
    return policy, gate_rel


def gsa_payload(classes: tuple[str, ...], *, negative: str | None = None, na_beta: str | None = None) -> bytes:
    lines = ["VARIABLE TYPE NGENES BETA SE P"]
    for index, class_id in enumerate(classes):
        beta = "NA" if class_id == na_beta else ("-1" if class_id == negative else "1")
        p_value = "0" if index == 0 else "0.001"
        lines.append(f"{class_id} COVAR 15000 {beta} 0.5 {p_value}")
    return ("\n".join(lines) + "\n").encode()


def normalized_family(
    root: Path,
    classes: tuple[str, ...],
    dataset_id: str,
    *,
    negative: str | None = None,
) -> list[dict[str, object]]:
    path = root / f"{dataset_id}.gsa.out"
    path.write_bytes(gsa_payload(classes, negative=negative))
    counts = E.EXPECTED_CLASS_COUNTS if dataset_id == "Allen_Human_MTG_level2" else {name: 1 for name in classes}
    aggregation = (
        "ARITHMETIC_MEAN_ALL_MAPPED_SUBTYPES_MATH_FSUM_SOURCE_ORDER"
        if dataset_id == "Allen_Human_MTG_level2" else "SOURCE_CLASS_COLUMN_NO_AGGREGATION"
    )
    return E.normalize_gsa(path, classes, dataset_id, "SleepTrait", aggregation, counts)


class FrozenLineageAndMappingTests(unittest.TestCase):
    def test_v1_v2_v3_are_still_byte_exact(self) -> None:
        expected = {
            B.V1_LOCK_REL: B.V1_LOCK_SHA256,
            B.V2_LOCK_REL: B.V2_LOCK_SHA256,
            E.PARENT_LOCK_REL: E.PARENT_LOCK_SHA256,
            "config/track_b_mechanism_followup_policy_v3.json": "30f6cd67b1227421664f55b3cb20f5c2139f9f0293365d687c42fa03309de0a4",
            "scripts/162_build_track_b_mechanism_readiness_v3.py": "28569632ecbeb842f07d0214a6de7912f18d089ce7f813f6aec16f867c7c1f1b",
        }
        for relative, digest in expected.items():
            self.assertEqual(E.stable_file(REPO, relative).sha256, digest)

    def test_policy_is_canonical_result_blind_and_executor_pinned(self) -> None:
        policy, record = E.load_policy(REPO)
        self.assertEqual(E.canonical_json(policy), (REPO / E.POLICY_REL).read_bytes())
        self.assertTrue(policy["result_blind"])
        self.assertFalse(policy["future_results_accessed"])
        self.assertEqual(policy["executors"]["enrichment"]["sha256"], E.stable_file(REPO, E.SELF_REL).sha256)
        self.assertEqual(policy["executors"]["classifier"]["sha256"], E.stable_file(REPO, E.CLASSIFIER_REL).sha256)
        self.assertEqual(record.sha256, hashlib.sha256(E.canonical_json(policy)).hexdigest())

    def test_real_source_builds_exact_complete_mapping_and_alignment(self) -> None:
        artifacts, _lock = B.construct_artifacts(REPO)
        mapping = E.parse_tsv_payload(artifacts[E.MAPPING_REL], E.MAPPING_FIELDS, "mapping")
        self.assertEqual(len(mapping), 75)
        self.assertEqual([int(row["source_position"]) for row in mapping], list(range(1, 76)))
        self.assertEqual(len({row["source_column"] for row in mapping}), 75)
        self.assertNotIn("Average", {row["source_column"] for row in mapping})
        self.assertEqual(
            {name: sum(row["broad_class"] == name for row in mapping) for name in E.CLASS_ORDER},
            E.EXPECTED_CLASS_COUNTS,
        )
        alignment = E.parse_tsv_payload(artifacts[E.ALIGNMENT_REL], E.ALIGNMENT_FIELDS, "alignment")
        hybrid = next(row for row in alignment if row["gse_class"] == "hybrid")
        neurons = next(row for row in alignment if row["gse_class"] == "neurons")
        self.assertEqual(hybrid["allen_broad_classes"], "NA")
        self.assertEqual(hybrid["gse_bh_denominator_inclusion"], "YES")
        self.assertEqual(neurons["allen_broad_classes"], "EXCITATORY_NEURONS;INHIBITORY_NEURONS")
        self.assertEqual(neurons["independent_replication_credit_rule"], "ONE_SHARED_TARGET_MAXIMUM")
        self.assertNotIn("PERICYTE", artifacts[E.MAPPING_REL].decode().upper())
        self.assertNotIn("MOUSE", artifacts[E.ALIGNMENT_REL].decode().upper())

    def test_mapping_refuses_missing_duplicate_and_unmapped_columns(self) -> None:
        columns: list[str] = []
        for prefix, broad in E.PREFIX_TO_CLASS:
            columns.extend(f"{prefix}{index}" for index in range(E.EXPECTED_CLASS_COUNTS[broad]))
        self.assertEqual(len(E.class_mapping_rows(columns)), 75)
        with self.assertRaisesRegex(E.ContractError, "exactly 75"):
            E.class_mapping_rows(columns[:-1])
        duplicated = list(columns)
        duplicated[-1] = duplicated[0]
        with self.assertRaisesRegex(E.ContractError, "unique"):
            E.class_mapping_rows(duplicated)
        unmapped = list(columns)
        unmapped[-1] = "Pericyte_RGS5"
        with self.assertRaisesRegex(E.ContractError, "map exactly once"):
            E.class_mapping_rows(unmapped)


class CompleteFamilyStatisticsTests(unittest.TestCase):
    def test_bh_requires_exactly_seven_and_retains_literal_zero(self) -> None:
        with self.assertRaisesRegex(E.ContractError, "seven"):
            E.bh_adjust(list("abcdef"), [0.1] * 6)
        adjusted = E.bh_adjust(list("abcdefg"), [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 1.0])
        self.assertEqual(adjusted["a"], 0.0)

    def test_na_and_missing_rows_are_retained_with_p_one_for_bh(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "test.gsa.out"
            path.write_bytes(gsa_payload(E.GSE_ORDER, na_beta="endothelial"))
            rows = E.normalize_gsa(
                path, E.GSE_ORDER, "GSE67835_Human_Cortex_woFetal", "SleepTrait",
                "SOURCE_CLASS_COLUMN_NO_AGGREGATION", {name: 1 for name in E.GSE_ORDER},
            )
            self.assertEqual(len(rows), 7)
            by_id = {row["class_id"]: row for row in rows}
            self.assertEqual(by_id["astrocytes"]["p_value"], "0")
            self.assertEqual(by_id["astrocytes"]["p_value_for_bh"], "0")
            self.assertEqual(by_id["endothelial"]["analysis_status"], "FAILED_NA_OR_NONFINITE")
            self.assertEqual(by_id["endothelial"]["p_value_for_bh"], "1")
            path.write_bytes(b"VARIABLE TYPE NGENES BETA SE P\nastrocytes COVAR 15000 1 1 0.01\n")
            rows = E.normalize_gsa(
                path, E.GSE_ORDER, "GSE67835_Human_Cortex_woFetal", "SleepTrait",
                "SOURCE_CLASS_COLUMN_NO_AGGREGATION", {name: 1 for name in E.GSE_ORDER},
            )
            self.assertEqual(len(rows), 7)
            self.assertEqual(rows[-1]["analysis_status"], "FAILED_MISSING_RESULT_ROW")
            self.assertEqual(rows[-1]["p_value_for_bh"], "1")

    def test_subprocess_failure_forces_all_seven_rows_out_of_inference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "failed.gsa.out"
            path.write_bytes(gsa_payload(E.GSE_ORDER))
            rows = E.normalize_gsa(
                path, E.GSE_ORDER, "GSE67835_Human_Cortex_woFetal", "SleepTrait",
                "SOURCE_CLASS_COLUMN_NO_AGGREGATION", {name: 1 for name in E.GSE_ORDER},
                forced_failure="FAILED_MAGMA_EXIT",
            )
            self.assertEqual({row["analysis_status"] for row in rows}, {"FAILED_MAGMA_EXIT"})
            self.assertEqual({row["p_value_for_bh"] for row in rows}, {"1"})
            self.assertEqual({row["bh_fdr"] for row in rows}, {"1"})

    def test_classifier_requires_two_corrected_positive_atlas_results(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            allen = normalized_family(root, E.CLASS_ORDER, "Allen_Human_MTG_level2")
            gse = normalized_family(root, E.GSE_ORDER, "GSE67835_Human_Cortex_woFetal")
            alignment = E.parse_tsv_payload(
                E.tsv_bytes(E.ALIGNMENT_FIELDS, E.alignment_rows()), E.ALIGNMENT_FIELDS, "alignment",
            )
            rows, summary = C.classify_replication("SleepTrait", allen, gse, alignment)
            self.assertEqual(len(rows), 8)
            self.assertEqual(summary["independent_one_to_one_replication_count"], 5)
            self.assertEqual(summary["partial_coarse_neuron_allen_row_count"], 2)
            self.assertEqual(summary["partial_coarse_neuron_shared_gse_target_count"], 1)
            neurons = [row for row in rows if row["gse_class"] == "neurons"]
            self.assertEqual({row["replication_classification"] for row in neurons}, {"PARTIAL_COARSE_NEURON"})
            self.assertEqual({row["independent_replication_credit"] for row in neurons}, {0})
            hybrid = next(row for row in rows if row["gse_class"] == "hybrid")
            self.assertEqual(hybrid["allen_class"], "NA")
            self.assertEqual(hybrid["replication_classification"], "UNALIGNED_HYBRID_DENOMINATOR_RETAINED")

    def test_negative_gse_neurons_cannot_replicate_either_neuron_class(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            allen = normalized_family(root, E.CLASS_ORDER, "Allen_Human_MTG_level2")
            gse = normalized_family(root, E.GSE_ORDER, "GSE67835_Human_Cortex_woFetal", negative="neurons")
            alignment = E.parse_tsv_payload(
                E.tsv_bytes(E.ALIGNMENT_FIELDS, E.alignment_rows()), E.ALIGNMENT_FIELDS, "alignment",
            )
            rows, summary = C.classify_replication("SleepTrait", allen, gse, alignment)
            neurons = [row for row in rows if row["gse_class"] == "neurons"]
            self.assertEqual({row["replication_classification"] for row in neurons}, {"FAILED_REPLICATION_DIRECTION"})
            self.assertEqual(summary["partial_coarse_neuron_shared_gse_target_count"], 0)


class UpstreamGateTests(unittest.TestCase):
    def test_complete_candidate_family_deterministically_reproduces_eligible_family(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            policy, gate_rel = gate_fixture(root)
            gate = E.validate_upstream_gate(root, gate_rel, policy)
            self.assertEqual(gate["target_source_id"], "SRC_SLEEP")
            self.assertEqual(gate["eligible_family"]["rows"], 2)

    def test_partial_result_selected_family_cannot_become_pseudo_ready(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            policy, gate_rel = gate_fixture(root, selected_indexes=[0])
            with self.assertRaisesRegex(E.ContractError, "exact deterministic subset"):
                E.validate_upstream_gate(root, gate_rel, policy)

    def test_wrong_ancestry_and_symlink_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            policy, gate_rel = gate_fixture(root, ancestry="AFR")
            with self.assertRaisesRegex(E.ContractError, "source-verified EUR"):
                E.validate_upstream_gate(root, gate_rel, policy)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            policy, gate_rel = gate_fixture(root)
            candidate = root / "results/track_b/upstream/complete_candidates.tsv"
            payload = candidate.read_bytes()
            target = write(root, "results/track_b/upstream/other.tsv", payload)
            candidate.unlink()
            candidate.symlink_to(target)
            with self.assertRaisesRegex(E.ContractError, "regular file"):
                E.validate_upstream_gate(root, gate_rel, policy)

    def test_execution_does_not_create_temp_or_output_before_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                mock.patch.object(E, "validate_runtime_contract", return_value={"placeholder": True}),
                mock.patch.object(E, "validate_upstream_gate", side_effect=E.ContractError("gate absent")),
                mock.patch.object(E.tempfile, "TemporaryDirectory") as temporary_directory,
                mock.patch.object(E, "publish_directory_no_replace") as publish,
            ):
                with self.assertRaisesRegex(E.ContractError, "gate absent"):
                    E.run_production(root, "results/gate.json", "results/out")
            temporary_directory.assert_not_called()
            publish.assert_not_called()
            self.assertEqual(list(root.iterdir()), [])


class NoReplaceAndVerifyTests(unittest.TestCase):
    def fake_artifacts(self) -> dict[str, bytes]:
        paths = [
            E.MAPPING_REL, E.ALIGNMENT_REL, E.AMENDMENT_REL, E.READINESS_REL,
            E.REPORT_REL, E.SUPERSESSION_REL, E.CONTRACT_LOCK_REL,
        ]
        return {path: f"fixture:{path}\n".encode() for path in paths}

    def test_freeze_resumes_exact_partial_family_and_publishes_lock_last(self) -> None:
        artifacts = self.fake_artifacts()
        expected_order = list(artifacts)
        for stop_after in range(1, len(expected_order) + 1):
            with self.subTest(stop_after=stop_after), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                seen: list[str] = []

                def stop(relative: str) -> None:
                    seen.append(relative)
                    if len(seen) == stop_after:
                        raise RuntimeError("simulated crash")

                with self.assertRaisesRegex(RuntimeError, "simulated crash"):
                    B.freeze_no_replace(root, artifacts, after_link=stop)
                self.assertEqual(seen, expected_order[:stop_after])
                B.freeze_no_replace(root, artifacts)
                for relative, payload in artifacts.items():
                    self.assertEqual((root / relative).read_bytes(), payload)
                self.assertEqual(expected_order[-1], E.CONTRACT_LOCK_REL)

    def test_freeze_rejects_nonidentical_and_unknown_targets(self) -> None:
        artifacts = self.fake_artifacts()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write(root, E.MAPPING_REL, b"wrong\n")
            with self.assertRaisesRegex(B.E.ContractError, "nonidentical"):
                B.freeze_no_replace(root, artifacts)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write(root, "results/track_b/mechanism_followup/v4/unknown.txt", b"unknown\n")
            with self.assertRaisesRegex(B.E.ContractError, "unknown"):
                B.freeze_no_replace(root, artifacts)

    def test_readiness_is_resource_ready_but_not_scientific_result_ready(self) -> None:
        rows = B.readiness_rows()
        self.assertEqual({row["resource_status"] for row in rows}, {"READY"})
        self.assertEqual({row["status"] for row in rows}, {"NOT_APPLICABLE_UNTIL_UPSTREAM"})
        self.assertTrue(all("complete candidate" in row["required_upstream_gate"] for row in rows))

    @unittest.skipUnless((REPO / E.CONTRACT_LOCK_REL).exists(), "V4 not frozen yet")
    def test_real_verify_is_side_effect_free(self) -> None:
        tracked = [
            REPO / E.POLICY_REL, REPO / B.SELF_REL, REPO / E.SELF_REL,
            REPO / E.CLASSIFIER_REL, REPO / "results/track_b/mechanism_followup/v4",
            REPO / E.SUPERSESSION_REL,
        ]

        def snapshot() -> dict[str, tuple[int, int, int, int, int, str]]:
            paths: list[Path] = []
            for item in tracked:
                if item.is_dir():
                    paths.extend(sorted(path for path in item.rglob("*") if path.is_file()))
                elif item.exists():
                    paths.append(item)
            return {
                str(path.relative_to(REPO)): (
                    path.stat().st_dev, path.stat().st_ino, path.stat().st_mode,
                    path.stat().st_size, path.stat().st_mtime_ns,
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                )
                for path in paths
            }

        before = snapshot()
        result = subprocess.run(
            [sys.executable, str(REPO / B.SELF_REL), "--verify", "--root", str(REPO)],
            cwd=REPO, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(before, snapshot())


if __name__ == "__main__":
    unittest.main()
