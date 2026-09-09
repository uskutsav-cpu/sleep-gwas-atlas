#!/usr/bin/env python3
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/160_build_track_b_mechanism_readiness.py"
SPEC = importlib.util.spec_from_file_location("track_b_mechanism_readiness_tested", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def write(root: Path, relative: str, payload: bytes = b"payload\n") -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return path


def production_policy() -> dict[str, object]:
    return json.loads((REPO / M.POLICY_REL).read_text(encoding="utf-8"))


def fuma_manifest() -> dict[str, object]:
    rows = [
        ("Allen_Human_MTG_level2", "brain", "human", "middle_temporal_gyrus", 75),
        ("GSE67835_Human_Cortex_woFetal", "brain", "human", "cerebral_cortex", 7),
        ("GSE89232_Human_Blood", "immune", "human", "blood_dendritic_cell_compartment", 4),
        ("PBMC_10x_68k", "immune", "human", "peripheral_blood_mononuclear_cells", 11),
        ("GSE81547_Human_Pancreas", "metabolic", "human", "pancreas", 7),
        ("GSE84133_Human_Pancreas", "metabolic", "human", "pancreas", 14),
        ("GSE98816_Mouse_Brain_Vascular", "vascular", "mouse_mapped_to_human_Ensembl_gene_ID", "brain_vasculature", 15),
        ("GSE99235_Mouse_Lung_Vascular", "vascular", "mouse_mapped_to_human_Ensembl_gene_ID", "lung_vasculature", 17),
    ]
    return {"datasets": [
        {"dataset_id": identity, "domain": domain, "species": species, "tissue": tissue, "cell_type_count": count}
        for identity, domain, species, tissue, count in rows
    ]}


def valid_build_manifests() -> dict[str, dict[str, object]]:
    return {
        "molecular": {"reference_build": {"analysis_build": "GRCh37", "molecular_source_build": "GRCh38"}},
        "interpretation_screen_registry_v4": {"genome_build": "GRCh38"},
        "interpretation_catlas_adult_v4": {"genome_build": "GRCh38", "life_stage": "Adult", "source_counts": {"adult_nuclei": 615998}},
        "interpretation_abc_2021": {"genome_build": "GRCh37"},
        "interpretation_pchic_2016": {"genome_build": "GRCh37"},
        "interpretation_causal_runtime": {"ld_reference": {"build": "GRCh37", "ancestry": "European", "sample_count": 503}},
    }


def fake_artifacts() -> tuple[dict[str, bytes], dict[str, object]]:
    outputs = {relative: f"content:{relative}\n".encode() for relative in M.OUTPUT_RELS}
    lock = {
        "schema_version": M.SCHEMA,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK",
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": M.POLICY_REL,
        "policy_sha256": "0" * 64,
        "script_sha256": "1" * 64,
        "analysis_family_count": len(M.EXPECTED_ANALYSIS_IDS),
        "evidence_record_count": 1,
        "evidence_bundle_sha256": "2" * 64,
        "outputs": {
            relative: {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
            for relative, payload in sorted(outputs.items())
        },
    }
    artifacts = dict(outputs)
    artifacts[M.LOCK_REL] = M.canonical_json(lock)
    return artifacts, lock


class SecureEvidenceTests(unittest.TestCase):
    def test_missing_required_file_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            registry = M.EvidenceRegistry(Path(temporary))
            with self.assertRaisesRegex(M.ContractError, "required evidence is missing"):
                registry.add("ref/missing.tsv")

    def test_missing_optional_file_is_explicit_not_null(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            registry = M.EvidenceRegistry(Path(temporary))
            identity = registry.add("ref/missing.tsv", required=False)
            row = registry.get(identity)
            self.assertEqual(row.kind, "MISSING")
            self.assertEqual(row.sha256, "NA")
            self.assertEqual(row.validation_status, "MISSING")

    def test_truncated_file_fails_size_pin(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write(root, "ref/value.tsv", b"short")
            with self.assertRaisesRegex(M.ContractError, "byte-size mismatch"):
                M.EvidenceRegistry(root).add("ref/value.tsv", expected_bytes=99)

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_symlink_file_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = write(root, "real/value.tsv")
            (root / "ref").mkdir()
            os.symlink(target, root / "ref/value.tsv")
            with self.assertRaisesRegex(M.ContractError, "symlink is forbidden"):
                M.stable_file(root, "ref/value.tsv")

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_symlink_inside_tree_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tree = root / "tree"
            tree.mkdir()
            target = write(root, "outside/value")
            os.symlink(target, tree / "link")
            with self.assertRaisesRegex(M.ContractError, "symlink is forbidden"):
                M.stable_tree(root, "tree")

    def test_in_read_drift_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = write(root, "large.bin", b"a" * (9 * 1024 * 1024))
            changed = False

            def mutate(block: int) -> None:
                nonlocal changed
                if block == 1 and not changed:
                    with path.open("ab") as handle:
                        handle.write(b"x")
                    changed = True

            with self.assertRaisesRegex(M.ContractError, "drifted while hashing"):
                M.stable_file(root, "large.bin", chunk_hook=mutate)

    def test_future_result_paths_are_forbidden(self) -> None:
        for relative in M.FORBIDDEN_RESULT_PREFIXES:
            with self.assertRaisesRegex(M.ContractError, "future science/result path"):
                M.safe_relative(relative + "/leak.tsv")

    def test_parent_traversal_and_absolute_paths_are_forbidden(self) -> None:
        for relative in ("../escape", "/absolute", "a/../../escape"):
            with self.assertRaises(M.ContractError):
                M.safe_relative(relative)


class ScientificSemanticsTests(unittest.TestCase):
    def test_wrong_species_fails_human_replication(self) -> None:
        manifest = fuma_manifest()
        for row in manifest["datasets"]:
            if row["dataset_id"] == "GSE67835_Human_Cortex_woFetal":
                row["species"] = "mouse_mapped_to_human_Ensembl_gene_ID"
        with self.assertRaisesRegex(M.ContractError, "species/context/count mismatch"):
            M.validate_fuma_species(manifest)

    def test_mouse_vascular_is_not_relabelled_human(self) -> None:
        manifest = fuma_manifest()
        observed = M.validate_fuma_species(manifest)
        self.assertNotEqual(observed["GSE98816_Mouse_Brain_Vascular"]["species"], "human")
        self.assertNotEqual(observed["GSE99235_Mouse_Lung_Vascular"]["species"], "human")

    def test_wrong_build_fails(self) -> None:
        manifests = valid_build_manifests()
        manifests["interpretation_catlas_adult_v4"]["genome_build"] = "GRCh37"
        with self.assertRaisesRegex(M.ContractError, "wrong genome build"):
            M.validate_build_semantics(manifests)

    def test_complete_analysis_family_required(self) -> None:
        policy = production_policy()
        policy["analyses"] = policy["analyses"][:-1]
        with self.assertRaisesRegex(M.ContractError, "incomplete, duplicated, or out of order"):
            M.validate_policy_shape(policy)

    def test_duplicate_analysis_family_rejected(self) -> None:
        policy = production_policy()
        policy["analyses"][-1] = copy.deepcopy(policy["analyses"][-2])
        with self.assertRaisesRegex(M.ContractError, "incomplete, duplicated, or out of order"):
            M.validate_policy_shape(policy)

    def test_policy_cannot_assert_pseudo_ready(self) -> None:
        policy = production_policy()
        policy["analyses"][0]["status"] = "READY"
        with self.assertRaisesRegex(M.ContractError, "unexpected schema"):
            M.validate_policy_shape(policy)

    def test_status_is_derived_and_data_blocker_precedes_software(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write(root, "evidence")
            registry = M.EvidenceRegistry(root)
            eid = registry.add("evidence")
            facts = {
                name: M.Fact(name, True, blocker_type, "ready", (eid,))
                for name, blocker_type in M.FACT_TYPES.items()
            }
            facts["qtl_payloads"] = M.Fact("qtl_payloads", False, "DATA", "payload absent, never null", (eid,))
            facts["tabix"] = M.Fact("tabix", False, "SOFTWARE", "tabix absent", (eid,))
            rows = M.evaluate_analyses(production_policy(), facts, registry)
            by_id = {row["analysis_id"]: row for row in rows}
            self.assertEqual(by_id["P17_CELL_EQTL"]["status"], "BLOCKED_BY_DATA")
            self.assertEqual(by_id["P17_CELL_EQTL"]["blocker_types"], "DATA;SOFTWARE")
            self.assertIn("never null", by_id["P17_CELL_EQTL"]["status_reason"])
            self.assertEqual(by_id["P14_BRAIN_CELL_CLASS_DISCOVERY"]["resource_status"], "READY")
            self.assertEqual(by_id["P14_BRAIN_CELL_CLASS_DISCOVERY"]["status"], "NOT_APPLICABLE_UNTIL_UPSTREAM")
            self.assertNotIn("READY", {row["status"] for row in rows})

    def test_all_records_have_ram_output_and_claim_limits(self) -> None:
        policy = production_policy()
        M.validate_policy_shape(policy)
        for row in policy["analyses"]:
            self.assertTrue(row["ram_decomposition"].strip())
            self.assertTrue(row["output"].startswith("results/track_b/mechanism_followup/science/"))
            self.assertTrue(row["claim_limit"].strip())


class FreezeVerifyTests(unittest.TestCase):
    def test_no_replace_preserves_existing_target(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, _lock = fake_artifacts()
            existing = write(root, M.OUTPUT_RELS[0], b"user-owned\n")
            with self.assertRaisesRegex(M.ContractError, "no-replace freeze refused"):
                M.freeze_no_replace(root, artifacts)
            self.assertEqual(existing.read_bytes(), b"user-owned\n")
            self.assertFalse((root / M.LOCK_REL).exists())

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks unavailable")
    def test_symlink_output_directory_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "results/track_b").mkdir(parents=True)
            elsewhere = root / "elsewhere"
            elsewhere.mkdir()
            os.symlink(elsewhere, root / M.OUTPUT_DIR_REL)
            with self.assertRaisesRegex(M.ContractError, "unsafe output directory"):
                M.freeze_no_replace(root, fake_artifacts()[0])

    def test_verify_is_side_effect_free(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, lock = fake_artifacts()
            M.freeze_no_replace(root, artifacts)

            def snapshot() -> dict[str, tuple[int, str]]:
                return {
                    path.relative_to(root).as_posix(): (path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest())
                    for path in sorted(root.rglob("*")) if path.is_file()
                }

            before = snapshot()
            with mock.patch.object(M, "construct_artifacts", return_value=(artifacts, lock)):
                M.verify_frozen(root)
            self.assertEqual(before, snapshot())

    def test_tampered_output_fails_before_recompute(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, _lock = fake_artifacts()
            M.freeze_no_replace(root, artifacts)
            (root / M.OUTPUT_RELS[1]).write_bytes(b"tampered\n")
            with self.assertRaisesRegex(M.ContractError, "differs from lock"):
                M.verify_frozen(root)

    def test_canonical_serializers_are_deterministic(self) -> None:
        value = {"z": 1, "a": [2, 3]}
        self.assertEqual(M.canonical_json(value), M.canonical_json(value))
        rows = [{"a": "x", "b": "y\nz"}]
        self.assertEqual(M.tsv_bytes(("a", "b"), rows), b"a\tb\nx\ty z\n")


if __name__ == "__main__":
    unittest.main()
