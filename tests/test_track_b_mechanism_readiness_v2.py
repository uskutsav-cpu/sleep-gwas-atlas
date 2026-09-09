#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/161_build_track_b_mechanism_readiness_v2.py"
SPEC = importlib.util.spec_from_file_location("track_b_mechanism_readiness_v2_tested", SCRIPT)
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
    paths = list(M.OUTPUT_RELS) + [M.SUPERSESSION_REL]
    outputs = {relative: f"v2:{relative}\n".encode() for relative in paths}
    lock = {
        "schema_version": M.SCHEMA,
        "contract_revision": 2,
        "contract_kind": "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK",
        "result_blind": True,
        "future_results_accessed": False,
        "policy_path": M.POLICY_REL,
        "policy_sha256": "0" * 64,
        "effective_policy_sha256": "1" * 64,
        "script_sha256": "2" * 64,
        "superseded_lock_path": M.OLD_LOCK_REL,
        "superseded_lock_sha256": "3" * 64,
        "analysis_family_count": len(M.B.EXPECTED_ANALYSIS_IDS),
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


class RevisionTwoLineageTests(unittest.TestCase):
    def test_revision_one_is_preserved_byte_exact(self) -> None:
        expected = {
            M.B.POLICY_REL: "503d34e1cdebac5c1e427c41a71b9c320bde948828ddf18df9767d7ac4e06c7e",
            M.BASE_REL: M.BASE_SHA256,
            M.OLD_LOCK_REL: "4ab3c004d109fa856936a5411f4495975b176333916bbe3f70b7974aac9de40f",
        }
        for relative, digest in expected.items():
            self.assertEqual(M.B.stable_file(REPO, relative).sha256, digest)

    def test_superseded_family_is_complete_and_byte_exact(self) -> None:
        overlay = M.B.read_json(REPO, M.POLICY_REL)
        lock = M.validate_superseded_family(REPO, overlay)
        self.assertEqual(set(lock["outputs"]), set(M.OLD_OUTPUT_RELS))

    def test_exact_effective_policy_diff_is_leafcutter_only(self) -> None:
        _overlay, effective, _hash = M.load_effective_policy(REPO)
        base = M.B.read_json(REPO, M.B.POLICY_REL)
        self.assertEqual(effective["required_qtl_datasets"]["eQTL"], base["required_qtl_datasets"]["eQTL"])
        self.assertEqual(effective["required_qtl_datasets"]["sQTL"], ["QTD000563", "QTD000573"])
        p18 = next(row for row in effective["analyses"] if row["analysis_id"] == "P18_CELL_SQTL")
        self.assertIn("leafcutter", p18["scope"])
        self.assertIn("QTD000563 and QTD000573", p18["required_inputs"])
        self.assertNotIn("QTD000560", p18["required_inputs"])
        self.assertNotIn("QTD000570", p18["required_inputs"])
        base_other = [row for row in base["analyses"] if row["analysis_id"] != "P18_CELL_SQTL"]
        effective_other = [row for row in effective["analyses"] if row["analysis_id"] != "P18_CELL_SQTL"]
        self.assertEqual(base_other, effective_other)

    def test_overlay_is_explicitly_result_blind(self) -> None:
        overlay = M.B.read_json(REPO, M.POLICY_REL)
        self.assertIs(overlay["result_blind"], True)
        self.assertIs(overlay["future_results_accessed"], False)
        self.assertEqual(overlay["supersedes"]["retention"], "PRESERVE_BYTE_EXACT_AS_AUDIT_TRAIL_DO_NOT_USE_FOR_EXECUTION")


class RevisionTwoFreezeTests(unittest.TestCase):
    def test_no_replace_refuses_existing_supersession_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, _lock = fake_artifacts()
            receipt = write(root, M.SUPERSESSION_REL, b"existing receipt\n")
            with self.assertRaisesRegex(M.B.ContractError, "no-replace revision-2 freeze refused"):
                M.freeze_no_replace(root, artifacts)
            self.assertEqual(receipt.read_bytes(), b"existing receipt\n")
            self.assertFalse((root / M.LOCK_REL).exists())

    def test_freeze_and_verify_are_side_effect_free_after_freeze(self) -> None:
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

    def test_tamper_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts, _lock = fake_artifacts()
            M.freeze_no_replace(root, artifacts)
            (root / M.SUPERSESSION_REL).write_bytes(b"tampered\n")
            with self.assertRaisesRegex(M.B.ContractError, "differs from lock"):
                M.verify_frozen(root)


if __name__ == "__main__":
    unittest.main()
