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
SCRIPT = REPO / "scripts/162_build_track_b_mechanism_readiness_v3.py"
SPEC = importlib.util.spec_from_file_location("track_b_mechanism_readiness_v3_tested", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def write(root: Path, relative: str, payload: bytes) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return path


def fake_artifacts() -> tuple[dict[str, bytes], dict[str, object]]:
    output_paths = list(M.OUTPUT_RELS) + [M.SUPERSESSION_REL]
    outputs = {relative: f"v3:{relative}\n".encode() for relative in output_paths}
    lock = {
        "schema_version": M.SCHEMA,
        "contract_revision": 3,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK",
        "authoritative_for_execution": True,
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": M.POLICY_REL,
        "policy_sha256": "0" * 64,
        "effective_policy_sha256": "1" * 64,
        "script_path": M.SCRIPT_REL,
        "script_sha256": "2" * 64,
        "superseded_lock_path": M.BASE_LOCK_REL,
        "superseded_lock_sha256": M.BASE_LOCK_SHA256,
        "retained_revision_1_lock_path": M.M2.OLD_LOCK_REL,
        "retained_revision_1_lock_sha256": "3" * 64,
        "analysis_family_count": len(M.EXPECTED_ANALYSIS_IDS),
        "evidence_record_count": 1,
        "evidence_bundle_sha256": "4" * 64,
        "outputs": {
            relative: {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
            for relative, payload in sorted(outputs.items())
        },
    }
    artifacts = dict(outputs)
    artifacts[M.LOCK_REL] = M.B.canonical_json(lock)
    return artifacts, lock


def gtex_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for index in range(49):
        group = f"Context_{index:02d}"
        rows.append({
            "study": "GTEx_V8",
            "qtl_group": group,
            "tissue_ontology_id": f"UBER_{index:07d}",
            "tissue_ontology_term": f"term {index}",
            "tissue_label": f"Tissue {index}",
            "condition_label": "naive",
            "quant_method": "ge",
            "ftp_path": f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/{group}.tsv.gz",
        })
    return rows


def audit_fixture() -> tuple[list[dict[str, object]], list[str]]:
    _overlay, effective, _hash, _lock = M.load_effective_policy(REPO)
    gtex_groups = [f"Context_{index:02d}" for index in range(49)]
    ready_ids = {
        "P15_BRAIN_CELL_SUBTYPE_DISCOVERY", "P20_CATLAS_SCATAC",
        "P20_SCREEN_CHROMATIN", "P21_ABC_ENHANCER_GENE",
        "P21_PCHIC_ENHANCER_GENE", "P26_POSITIVE_CONTROL", "P27_PATHWAYS",
    }
    software_ids = {
        "P07_CONJFDR", "P14_BRAIN_CELL_CLASS_DISCOVERY",
        "P16_HUMAN_CELL_CLASS_REPLICATION",
        "P24_LOCKED_ATLAS_CROSS_SLEEP", "P25_PAIR_A_VS_B",
        "P29_BIDIRECTIONAL_MR",
    }
    rows: list[dict[str, object]] = []
    for declared in effective["analyses"]:
        identity = declared["analysis_id"]
        resource, status = "BLOCKED_BY_DATA", "BLOCKED_BY_DATA"
        if identity in ready_ids:
            resource, status = "READY", "NOT_APPLICABLE_UNTIL_UPSTREAM"
        elif identity in software_ids:
            resource, status = "BLOCKED_BY_SOFTWARE", "BLOCKED_BY_SOFTWARE"
        paths = {f"evidence/{identity}"}
        if identity in M.INTERPRETATION_CODE_ANALYSIS_IDS:
            paths.update(M.INTERPRETATION_CODE_PATHS)
        if identity in {"P17_CELL_EQTL", "P19_THREE_WAY_EQTL", "P22_REGULATORY_CHAIN_EQTL"}:
            paths.update(f"{M.QTL_ROOT_REL}/{value}.tsv.gz" for value in M.CELL_EQTL_IDS)
        if identity in {"P18_CELL_SQTL", "P19_THREE_WAY_SQTL", "P22_REGULATORY_CHAIN_SQTL"}:
            paths.update(f"{M.QTL_ROOT_REL}/{value}.tsv.gz" for value in M.CELL_SQTL_IDS)
        if identity == "P13_GTEX_TISSUE_EQTL":
            paths.add(M.GTEX_METADATA_REL)
            for group in gtex_groups:
                paths.add(f"ref/molecular/eqtl_catalogue_r7/imported/GTEx_V8/ge/{group}.tsv.gz")
                paths.add(f"ref/molecular/eqtl_catalogue_r7/imported/GTEx_V8/ge/{group}.tsv.gz.tbi")
        rows.append({
            "analysis_id": identity,
            "requirements": declared["requirements"],
            "resource_status": resource,
            "status": status,
            "status_reason": "blocked or awaiting upstream; never interpreted as negative evidence",
            "exact_unblock_condition": declared["unblock"],
            "evidence_paths": ";".join(sorted(paths)),
            "planned_output": declared["output"],
        })
    return rows, gtex_groups


class RevisionThreeLineageTests(unittest.TestCase):
    def test_v1_and_v2_are_preserved_byte_exact(self) -> None:
        expected = {
            M.M2.B.POLICY_REL: "503d34e1cdebac5c1e427c41a71b9c320bde948828ddf18df9767d7ac4e06c7e",
            M.M2.BASE_REL: M.M2.BASE_SHA256,
            M.M2.OLD_LOCK_REL: "4ab3c004d109fa856936a5411f4495975b176333916bbe3f70b7974aac9de40f",
            M.BASE_POLICY_REL: M.BASE_POLICY_SHA256,
            M.BASE_REL: M.BASE_SHA256,
            M.BASE_LOCK_REL: M.BASE_LOCK_SHA256,
        }
        for relative, digest in expected.items():
            self.assertEqual(M.B.stable_file(REPO, relative).sha256, digest)

    def test_both_superseded_families_are_complete_and_byte_exact(self) -> None:
        overlay = M.B.read_json(REPO, M.POLICY_REL)
        lock = M.validate_superseded_families(REPO, overlay)
        self.assertEqual(set(lock["outputs"]), set(M.M2.OUTPUT_RELS) | {M.M2.SUPERSESSION_REL})

    def test_effective_family_is_exact_and_modality_separated(self) -> None:
        _overlay, effective, _hash, _lock = M.load_effective_policy(REPO)
        self.assertEqual([row["analysis_id"] for row in effective["analyses"]], list(M.EXPECTED_ANALYSIS_IDS))
        by_id = {row["analysis_id"]: row for row in effective["analyses"]}
        self.assertEqual(by_id["P17_CELL_EQTL"]["requirements"], ["cell_eqtl_metadata", "cell_eqtl_payloads", "tabix"])
        self.assertEqual(by_id["P18_CELL_SQTL"]["requirements"], ["cell_sqtl_metadata", "cell_sqtl_payloads", "tabix"])
        self.assertIn("brain_class_mapping_executor", by_id["P14_BRAIN_CELL_CLASS_DISCOVERY"]["requirements"])
        self.assertEqual(by_id["P14_HUMAN_BRAIN_PERICYTE_DISCOVERY"]["requirements"], ["human_brain_pericyte_reference"])
        self.assertEqual(
            by_id["P13_SLDSC_GTEX"]["requirements"],
            ["ldsc_source_bundle", "sldsc_broad_tissue_coverage", "ldsc_runtime", "ldsc_static_reference", "interpretation_code_family"],
        )
        self.assertIn("Lung index 36", by_id["P13_SLDSC_GTEX"]["required_inputs"])
        self.assertIn("Muscle_Skeletal index 38", by_id["P13_SLDSC_GTEX"]["required_inputs"])
        self.assertIn("pair_comparison_executor", by_id["P25_PAIR_A_VS_B"]["requirements"])
        self.assertEqual(len({row["output"] for row in effective["analyses"]}), 31)

    def test_overlay_is_explicitly_result_blind(self) -> None:
        overlay = M.B.read_json(REPO, M.POLICY_REL)
        self.assertIs(overlay["result_blind"], True)
        self.assertIs(overlay["future_results_accessed"], False)
        self.assertEqual(
            overlay["supersedes"]["retention"],
            "PRESERVE_V1_AND_V2_BYTE_EXACT_AS_AUDIT_TRAILS_DO_NOT_USE_FOR_EXECUTION",
        )


class RevisionThreeScientificBoundaryTests(unittest.TestCase):
    def test_every_requirement_blocks_when_any_exact_unit_is_unavailable(self) -> None:
        _overlay, effective, _hash, _lock = M.load_effective_policy(REPO)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write(root, "evidence/pin", b"immutable fixture\n")
            registry = M.B.EvidenceRegistry(root)
            evidence_id = registry.add("evidence/pin")
            all_ready = {
                name: M.B.Fact(name, True, blocker, "fixture ready", (evidence_id,))
                for name, blocker in M.FACT_TYPES.items()
            }
            fully_ready = {
                row["analysis_id"]: row
                for row in M.evaluate_analyses(effective, all_ready, registry)
            }
            for identity in M.EXPECTED_ANALYSIS_IDS:
                self.assertEqual(fully_ready[identity]["resource_status"], "READY")
            for declared in effective["analyses"]:
                identity = declared["analysis_id"]
                for missing in declared["requirements"]:
                    with self.subTest(analysis_id=identity, missing=missing):
                        facts = dict(all_ready)
                        facts[missing] = M.B.Fact(
                            missing, False, M.FACT_TYPES[missing],
                            "fixture unavailable", (evidence_id,),
                        )
                        observed = {
                            row["analysis_id"]: row
                            for row in M.evaluate_analyses(effective, facts, registry)
                        }
                        self.assertNotEqual(observed[identity]["resource_status"], "READY")

    def test_every_unblock_predicate_is_bound_to_its_exact_units(self) -> None:
        _overlay, effective, _hash, _lock = M.load_effective_policy(REPO)
        for declared in effective["analyses"]:
            identity = declared["analysis_id"]
            with self.subTest(analysis_id=identity):
                mutated = copy.deepcopy(effective)
                target = next(row for row in mutated["analyses"] if row["analysis_id"] == identity)
                target["unblock"] = "unrelated partial-family availability"
                with self.assertRaisesRegex(M.B.ContractError, "unblock predicate"):
                    M.validate_effective_policy(mutated)

    def test_exact_gtex_context_metadata_accepts_only_49_ge_rows(self) -> None:
        rows = gtex_rows()
        with mock.patch.object(M.B, "read_tsv", return_value=(M.GTEX_FIELDS, rows)):
            self.assertEqual(M.validate_gtex_contexts(Path("/unused")), [row["qtl_group"] for row in rows])
        wrong = copy.deepcopy(rows)
        wrong[0]["quant_method"] = "leafcutter"
        with mock.patch.object(M.B, "read_tsv", return_value=(M.GTEX_FIELDS, wrong)):
            with self.assertRaisesRegex(M.B.ContractError, "wrong release/method/path/schema"):
                M.validate_gtex_contexts(Path("/unused"))

    def test_present_unpinned_payload_never_becomes_pseudo_ready(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            relative = f"{M.QTL_ROOT_REL}/QTD000559.tsv.gz"
            write(root, relative, b"appeared but unvalidated\n")
            registry = M.B.EvidenceRegistry(root)
            evidence_id = registry.add(relative, required=False)
            fact = M._always_blocked_fact(
                "cell_eqtl_payloads", "DATA", (evidence_id,), registry,
                "absent", "present but unvalidated",
            )
            self.assertFalse(fact.ready)
            self.assertEqual(fact.reason, "present but unvalidated")

    def test_sldsc_missing_tissues_are_an_explicit_fail_closed_unit(self) -> None:
        overlay, effective, _hash, _lock = M.load_effective_policy(REPO)
        probe = str(overlay["new_absence_probes"]["sldsc_broad_tissue_coverage"])
        self.assertFalse((REPO / probe).exists())
        row = next(value for value in effective["analyses"] if value["analysis_id"] == "P13_SLDSC_GTEX")
        self.assertIn("sldsc_broad_tissue_coverage", row["requirements"])
        self.assertIn("all 18 tissues", row["ram_decomposition"])
        self.assertIn("airway-specific", row["claim_limit"])

    def test_qtl_evidence_cross_contamination_fails_closed(self) -> None:
        rows, groups = audit_fixture()
        self.assertEqual(M.audit_analysis_family(rows, groups)["analysis_ids_audited"], 31)
        p17 = next(row for row in rows if row["analysis_id"] == "P17_CELL_EQTL")
        p17["evidence_paths"] += f";{M.QTL_ROOT_REL}/QTD000563.tsv.gz"
        with self.assertRaisesRegex(M.B.ContractError, "eQTL/sQTL evidence contamination"):
            M.audit_analysis_family(rows, groups)

    def test_wrong_build_and_wrong_species_fail_closed(self) -> None:
        manifests = {
            "molecular": M.B.read_json(REPO, "config/molecular_analysis_policy.json"),
            "interpretation_screen_registry_v4": M.B.read_json(REPO, "config/interpretation_screen_registry_v4.json"),
            "interpretation_catlas_adult_v4": M.B.read_json(REPO, "config/interpretation_catlas_adult_v4.json"),
            "interpretation_abc_2021": M.B.read_json(REPO, "config/interpretation_abc_2021.json"),
            "interpretation_pchic_2016": M.B.read_json(REPO, "config/interpretation_pchic_2016.json"),
            "interpretation_causal_runtime": M.B.read_json(REPO, "config/interpretation_causal_runtime.json"),
        }
        wrong_build = copy.deepcopy(manifests)
        wrong_build["interpretation_screen_registry_v4"]["genome_build"] = "GRCh37"
        with self.assertRaisesRegex(M.B.ContractError, "wrong genome build"):
            M.B.validate_build_semantics(wrong_build)
        fuma = M.B.read_json(REPO, "config/interpretation_fuma_scrna.json")
        wrong_species = copy.deepcopy(fuma)
        next(row for row in wrong_species["datasets"] if row["dataset_id"] == "GSE67835_Human_Cortex_woFetal")["species"] = "mouse"
        with self.assertRaisesRegex(M.B.ContractError, "species/context/count mismatch"):
            M.B.validate_fuma_species(wrong_species)

    def test_symlink_and_drift_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = write(root, "real/data.tsv", b"one\n")
            link = root / "probe.tsv"
            link.symlink_to(target)
            registry = M.B.EvidenceRegistry(root)
            with self.assertRaisesRegex(M.B.ContractError, "symlink"):
                registry.add("probe.tsv", required=False)
            registry.add("real/data.tsv")
            target.write_bytes(b"two\n")
            with self.assertRaisesRegex(M.B.ContractError, "drifted"):
                registry.recheck()


class RevisionThreeFreezeTests(unittest.TestCase):
    def test_nonidentical_existing_target_is_never_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, _lock = fake_artifacts()
            target = write(root, M.SUPERSESSION_REL, b"existing\n")
            with self.assertRaisesRegex(M.B.ContractError, "refused nonidentical"):
                M.freeze_no_replace(root, artifacts)
            self.assertEqual(target.read_bytes(), b"existing\n")
            self.assertFalse((root / M.LOCK_REL).exists())

    def test_exact_partial_publication_resumes_after_every_link_boundary(self) -> None:
        ordered = list(M.OUTPUT_RELS) + [M.SUPERSESSION_REL, M.LOCK_REL]
        for stop_after in ordered:
            with self.subTest(stop_after=stop_after), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                artifacts, _lock = fake_artifacts()

                def interrupt(relative: str) -> None:
                    if relative == stop_after:
                        raise RuntimeError("simulated process loss")

                with self.assertRaisesRegex(RuntimeError, "simulated process loss"):
                    M.freeze_no_replace(root, artifacts, after_link=interrupt)
                M.freeze_no_replace(root, artifacts)
                for relative, payload in artifacts.items():
                    self.assertEqual((root / relative).read_bytes(), payload)

    def test_verify_is_side_effect_free_and_tamper_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, lock = fake_artifacts()
            M.freeze_no_replace(root, artifacts)

            def snapshot() -> dict[str, tuple[int, int, str]]:
                return {
                    path.relative_to(root).as_posix(): (
                        path.stat().st_ino, path.stat().st_size,
                        hashlib.sha256(path.read_bytes()).hexdigest(),
                    )
                    for path in sorted(root.rglob("*")) if path.is_file()
                }

            before = snapshot()
            with mock.patch.object(M, "construct_artifacts", return_value=(artifacts, lock)):
                M.verify_frozen(root)
            self.assertEqual(before, snapshot())
            (root / M.OUTPUT_RELS[0]).write_bytes(b"tampered\n")
            with self.assertRaisesRegex(M.B.ContractError, "differs from lock"):
                M.verify_frozen(root)


if __name__ == "__main__":
    unittest.main()
