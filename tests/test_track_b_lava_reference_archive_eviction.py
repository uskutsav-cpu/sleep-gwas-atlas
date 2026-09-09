import csv
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]


def load_module():
    path = REPO / "scripts/146_manage_lava_reference_archives.py"
    specification = importlib.util.spec_from_file_location("lava_archive_eviction_test", path)
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


GUARD = load_module()
SOURCE = GUARD.PINNED_SOURCE_FINGERPRINT
CONTINUATION = GUARD.PINNED_CONTINUATION_FINGERPRINT


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class Fixture:
    def __init__(self, root: Path, *, terminal: bool = True, integrated: bool = True):
        self.root = root
        integration = (
            'GUARD_PATH = "scripts/146_manage_lava_reference_archives.py"\n'
            f'ACCEPTED_ARCHIVE_STATES = {{"{GUARD.ARCHIVES_PRESENT}", '
            f'"{GUARD.ARCHIVES_EVICTED}"}}\n'
            "def verify_archive_state(guard):\n"
            "    return guard.verify_reference_state(None)\n"
            if integrated else ""
        )
        scripts = {
            "146_manage_lava_reference_archives.py": "# fixture archive guard implementation\n",
            "lava_contract.py": "require_live=True\narchives + payloads\n",
            "119_track_b_lava_contract.py": "lava_contract.validate_reference\n",
            "134_track_b_lava_continuation_contract.py": (
                "_require_live_source_fingerprint\nlegacy_contract().run_fingerprint()\n"
            ),
            "139_build_track_b_placo_terminal_gate_v2.py": (
                "def deep_lava_terminal(): pass\nvalidate_source_bundles=True\n" + integration
            ),
            "140_materialize_track_b_placo_pair_v2.py": "GATE.verify_gate()\n",
            "141_run_track_b_placo_pair_v2.R": "# worker has no archive dependency\n",
            "143_run_track_b_placo_sequential_v2.py": "GATE.verify_gate()\n" + integration,
            "145_build_track_b_finemapping_continuation_gate.py": "# additive gate\n" + integration,
            "fine_mapping_contract.py": "lava_contract.validate_reference\n",
            "58_materialize_finemapping_locus.py": 'suffixes = [".info", ".bcor"]\n',
        }
        for name, content in scripts.items():
            path = root / "scripts" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")

        archive_specs = [
            ("01_02", "1-2", "chr1-2", b"archive-1"),
            ("03_04", "3-4", "chr3-4", b"archive-2"),
            ("05_06", "5-6", "chr5-6", b"archive-3"),
            ("07_09", "7-9", "chr7-9", b"archive-4"),
            ("10_12", "10-12", "chr10-12", b"archive-5"),
            ("13_16", "13-16", "chr13-16", b"archive-6"),
            ("17_23", "17-22,X", "chr17-23", b"archive-7"),
        ]
        source_rows = []
        download_rows = []
        archive_records = []
        self.archive_paths = []
        for suffix, chromosomes, label, content in archive_specs:
            filename = f"lava-ukb-v1.1_{label}.zip"
            relative = Path("ref/lava/ukb_v1.1/archives") / filename
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            self.archive_paths.append(path)
            source = {
                "archive_id": f"lava_ukb_v1.1_{suffix}",
                "chromosomes": chromosomes,
                "url": f"https://example.test/{label}",
                "archive_bytes": str(len(content)),
                "archive_filename": filename,
            }
            source_rows.append(source)
            download_rows.append(source | {"sha256": sha(path)})
            archive_records.append({
                "path": str(relative), "bytes": len(content), "sha256": sha(path),
            })
        self.archive_bytes = sum(path.stat().st_size for path in self.archive_paths)
        write_tsv(root / GUARD.SOURCE_REGISTRY, GUARD.SOURCE_FIELDS, source_rows)
        write_tsv(root / GUARD.DOWNLOAD_MANIFEST, GUARD.DOWNLOAD_FIELDS, download_rows)

        extracted_rows = []
        extracted_records = []
        self.payload_paths = []
        for chromosome in range(1, 23):
            for file_type in ("info", "bcor"):
                filename = f"lava-ukb-v1.1_chr{chromosome}.{file_type}"
                relative = Path("ref/lava/ukb_v1.1") / filename
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(f"payload-{chromosome}-{file_type}".encode())
                self.payload_paths.append(path)
                extracted_rows.append({
                    "chromosome": chromosome, "file_type": file_type, "path": filename,
                    "bytes": path.stat().st_size, "sha256": sha(path),
                })
                extracted_records.append({
                    "path": str(relative), "bytes": path.stat().st_size, "sha256": sha(path),
                })
        write_tsv(root / GUARD.EXTRACTED_MANIFEST, GUARD.EXTRACTED_FIELDS, extracted_rows)

        source_identity = {
            "path": str(GUARD.SOURCE_REGISTRY),
            "bytes": (root / GUARD.SOURCE_REGISTRY).stat().st_size,
            "sha256": sha(root / GUARD.SOURCE_REGISTRY),
        }
        policy = {
            "analysis_id": "atlas-v1.0-lava",
            "lava_version": "0.1.5",
            "reference": "LAVA UK Biobank LD v1.1",
            "reference_population": "fixture EUR",
            "reference_expected_archives": 7,
            "reference_expected_archive_bytes": self.archive_bytes,
            "reference_expected_extracted_files": 44,
            "reference_source_registry": str(GUARD.SOURCE_REGISTRY),
            "reference_source_registry_sha256": source_identity["sha256"],
            "reference_download_manifest": str(GUARD.DOWNLOAD_MANIFEST),
            "reference_extracted_manifest": str(GUARD.EXTRACTED_MANIFEST),
            "reference_payload_root": str(GUARD.REFERENCE_ROOT),
            "reference_provenance": str(GUARD.REFERENCE_PROVENANCE),
        }
        write_json(root / GUARD.POLICY, policy)
        provenance = {
            "schema_version": "sleep-atlas-lava-reference.1",
            "analysis_id": policy["analysis_id"],
            "reference": policy["reference"],
            "reference_population": policy["reference_population"],
            "source_registry": str(GUARD.SOURCE_REGISTRY),
            "source_registry_sha256": source_identity["sha256"],
            "download_manifest": str(GUARD.DOWNLOAD_MANIFEST),
            "download_manifest_sha256": sha(root / GUARD.DOWNLOAD_MANIFEST),
            "extracted_manifest": str(GUARD.EXTRACTED_MANIFEST),
            "extracted_manifest_sha256": sha(root / GUARD.EXTRACTED_MANIFEST),
            "archive_count": 7,
            "archive_bytes": self.archive_bytes,
            "extracted_file_count": 44,
            "extracted_bytes": sum(item["bytes"] for item in extracted_records),
            "archives": archive_records,
            "extracted_files": extracted_records,
            "verification": "SHA-256 verified after official HTTPS acquisition",
            "contract_script_sha256": sha(root / GUARD.LEGACY_CONTRACT),
        }
        write_json(root / GUARD.REFERENCE_PROVENANCE, provenance)

        reference_identity = {
            "path": str(GUARD.REFERENCE_PROVENANCE),
            "bytes": (root / GUARD.REFERENCE_PROVENANCE).stat().st_size,
            "sha256": sha(root / GUARD.REFERENCE_PROVENANCE),
        }
        entries = [
            {
                "path": item["path"], "bytes": item["bytes"], "sha256": item["sha256"],
                "role": "REFERENCE_CHROMOSOME", "chromosome": index // 2 + 1,
            }
            for index, item in enumerate(extracted_records)
        ]
        entries.append(reference_identity | {"role": "COMMON_INPUT", "chromosome": None})
        entries.extend(
            {"path": f"common/{index}", "bytes": 1, "sha256": "1" * 64,
             "role": "COMMON_INPUT", "chromosome": None}
            for index in range(11)
        )
        entries.extend(
            {"path": f"shard/{index}", "bytes": 1, "sha256": "2" * 64,
             "role": "CHROMOSOME_SHARD", "chromosome": index // 8 + 1}
            for index in range(176)
        )
        entries.extend(
            {"path": f"input/{index}", "bytes": 1, "sha256": "3" * 64,
             "role": "CHROMOSOME_INPUT_INFO", "chromosome": index + 1}
            for index in range(22)
        )
        baseline_path = (
            root / "results/track_b/checkpoints/lava" / SOURCE
            / "active_input_identity_baseline.json"
        )
        write_json(baseline_path, {
            "schema_version": 1,
            "analysis_id": "track-b-v1.0-local",
            "execution_fingerprint": SOURCE,
            "content_verification": "FULL_SHA256_PREFLIGHT_BEFORE_BASELINE_PUBLICATION",
            "entries": entries,
        })
        baseline_record = {
            "path": str(baseline_path.relative_to(root)),
            "bytes": baseline_path.stat().st_size,
            "sha256": sha(baseline_path),
        }
        checkpoints = [{"locus_index": index} for index in range(1, 2496)]
        source_lock_path = root / GUARD.SOURCE_LOCK_ROOT / f"{SOURCE}.json"
        write_json(source_lock_path, {
            "schema_version": GUARD.SOURCE_LOCK_SCHEMA,
            "analysis_id": "track-b-v1.0-local",
            "source_discovery_fingerprint": SOURCE,
            "checkpoint_count": 2495,
            "checkpoint_family_sha256": GUARD.digest_json(checkpoints),
            "checkpoints": checkpoints,
            "active_input_baseline": baseline_record,
        })
        source_lock_record = {
            "path": str(source_lock_path.relative_to(root)),
            "bytes": source_lock_path.stat().st_size,
            "sha256": sha(source_lock_path),
        }
        lineage_path = root / GUARD.CONTINUATION_RUN_ROOT / CONTINUATION / "lineage.json"
        write_json(lineage_path, {
            "analysis_id": "track-b-v1.0-local",
            "source_discovery_fingerprint": SOURCE,
            "continuation_execution_fingerprint": CONTINUATION,
            "contract_payload_sha256": CONTINUATION,
            "source_family_lock": source_lock_record,
        })
        self.run = lineage_path.parent
        self.add_receipt("aggregate-discovery", None)
        if terminal:
            self.add_receipt("terminal-qc", None)
            terminal_receipt = self.receipt_identity("terminal-qc")
            terminal_ready = self.ready_identity("terminal-qc")
            write_json(root / GUARD.PLACO_TERMINAL_GATE, {
                "schema_version": "sleep-atlas-track-b-pleiotropy-terminal-gate.2",
                "lava_scientific_evidence": {"status": "FAILED_QC_NOT_CONSUMED"},
                "lava_terminal": {
                    "source_discovery_fingerprint": SOURCE,
                    "continuation_execution_fingerprint": CONTINUATION,
                    "receipt": terminal_receipt,
                    "ready": terminal_ready,
                },
            })

    def _attempt(self, phase: str) -> Path:
        return self.run / ".attempts" / f"{phase}_unit"

    def add_receipt(self, phase: str, locus_index: int | None) -> None:
        attempt = self._attempt(phase)
        write_json(attempt / "receipt.json", {
            "phase": phase, "locus_index": locus_index,
            "active_inputs": [{"path": "config/safe", "bytes": 1, "sha256": "4" * 64}],
        })
        write_json(attempt / "READY", {
            "state": "READY", "receipt_sha256": sha(attempt / "receipt.json"),
        })

    def receipt_identity(self, phase: str) -> dict[str, object]:
        path = self._attempt(phase) / "receipt.json"
        return {"path": str(path.relative_to(self.root)), "bytes": path.stat().st_size, "sha256": sha(path)}

    def ready_identity(self, phase: str) -> dict[str, object]:
        path = self._attempt(phase) / "READY"
        return {"path": str(path.relative_to(self.root)), "bytes": path.stat().st_size, "sha256": sha(path)}


class TrackBLavaReferenceArchiveEvictionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        self.byte_patch = patch.object(GUARD, "EXPECTED_ARCHIVE_BYTES", self.fixture.archive_bytes)
        self.byte_patch.start()

    def tearDown(self):
        self.byte_patch.stop()
        self.temporary.cleanup()

    def test_present_state_rehashes_exact_archive_and_payload_families(self):
        state = GUARD.verify_reference_state(self.root)
        self.assertEqual(state["state"], GUARD.ARCHIVES_PRESENT)
        self.assertEqual(state["archive_count_present"], 7)
        self.assertEqual(state["extracted_payload_count"], 44)

    def test_archive_or_payload_tampering_fails_closed(self):
        self.fixture.archive_paths[2].write_bytes(b"wrong")
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "archive differs"):
            GUARD.validate_reference_metadata(self.root, require_archives=True)
        self.fixture.archive_paths[2].write_bytes(b"archive-3")
        self.fixture.payload_paths[4].write_bytes(b"tampered")
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "payload differs"):
            GUARD.validate_reference_metadata(self.root, require_archives=True)

    def test_source_baseline_archive_dependency_is_rejected(self):
        lock = json.loads((self.root / GUARD.SOURCE_LOCK_ROOT / f"{SOURCE}.json").read_text())
        baseline_path = self.root / lock["active_input_baseline"]["path"]
        baseline = json.loads(baseline_path.read_text())
        baseline["entries"][0]["path"] = str(GUARD.REFERENCE_ROOT / "archives/bad.zip")
        write_json(baseline_path, baseline)
        lock["active_input_baseline"] = {
            "path": str(baseline_path.relative_to(self.root)),
            "bytes": baseline_path.stat().st_size,
            "sha256": sha(baseline_path),
        }
        write_json(self.root / GUARD.SOURCE_LOCK_ROOT / f"{SOURCE}.json", lock)
        reference = GUARD.validate_reference_metadata(self.root, require_archives=True)
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "depends on acquisition archives"):
            GUARD.validate_lava_dependency_baseline(self.root, reference)

    def test_plan_is_no_replace_and_records_operational_dependency(self):
        plan = GUARD.create_or_verify_plan(self.root)
        plan_id = GUARD.digest_json(plan)
        relative = GUARD.plan_path(plan_id)
        before = (self.root / relative).read_bytes()
        self.assertEqual(plan["state"], "AUDITED_NO_FILES_REMOVED")
        self.assertFalse(
            plan["downstream_dependency_audit"]["placo_current_operational_verifier"][
                "archive_bytes_required"
            ]
        )
        self.assertEqual(GUARD.create_or_verify_plan(self.root), plan)
        self.assertEqual((self.root / relative).read_bytes(), before)
        changed = json.loads(before)
        changed["state"] = "changed"
        write_json(self.root / relative, changed)
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "plan differs"):
            GUARD.create_or_verify_plan(self.root)

    def test_explicit_eviction_removes_only_allowlist_and_verifies_receipt(self):
        unrelated = self.root / GUARD.REFERENCE_ROOT / "archives/keep.txt"
        unrelated.write_text("keep", encoding="utf-8")
        plan = GUARD.create_or_verify_plan(self.root)
        plan_id = GUARD.digest_json(plan)
        receipt = GUARD.execute_eviction(
            root=self.root,
            plan_id=plan_id,
            confirmation=GUARD.CONFIRMATION,
            acknowledge_operational_breakage=True,
        )
        self.assertEqual(receipt["state"], GUARD.ARCHIVES_EVICTED)
        self.assertEqual(receipt["plan_id"], plan_id)
        self.assertTrue((self.root / GUARD.intent_path(plan_id)).is_file())
        self.assertTrue((self.root / GUARD.receipt_path(plan_id)).is_file())
        self.assertTrue(unrelated.is_file())
        self.assertTrue(all(not path.exists() for path in self.fixture.archive_paths))
        self.assertTrue(all(path.is_file() for path in self.fixture.payload_paths))
        state = GUARD.verify_reference_state(self.root)
        self.assertEqual(state["state"], GUARD.ARCHIVES_EVICTED)

    def test_interrupted_eviction_resumes_only_after_durable_intent(self):
        plan = GUARD.create_or_verify_plan(self.root)
        plan_id = GUARD.digest_json(plan)
        reference = GUARD.validate_reference_metadata(self.root, require_archives=True)
        baseline = GUARD.validate_lava_dependency_baseline(self.root, reference)
        intent = GUARD._intent_payload(self.root, plan_id, plan, reference, baseline)
        GUARD.publish_no_replace(self.root, GUARD.intent_path(plan_id), intent)
        self.fixture.archive_paths[0].unlink()
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, GUARD.PARTIAL_STATE):
            GUARD.verify_reference_state(self.root)
        receipt = GUARD.execute_eviction(
            root=self.root,
            plan_id=plan_id,
            confirmation=GUARD.CONFIRMATION,
            acknowledge_operational_breakage=True,
        )
        self.assertEqual(receipt["state"], GUARD.ARCHIVES_EVICTED)

    def test_eviction_is_forbidden_without_terminal_lineage(self):
        other_root = Path(tempfile.mkdtemp(dir=self.root.parent))
        try:
            fixture = Fixture(other_root, terminal=False)
            with patch.object(GUARD, "EXPECTED_ARCHIVE_BYTES", fixture.archive_bytes):
                plan = GUARD.create_or_verify_plan(other_root)
                with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "not terminal"):
                    GUARD.execute_eviction(
                        root=other_root,
                        plan_id=GUARD.digest_json(plan),
                        confirmation=GUARD.CONFIRMATION,
                        acknowledge_operational_breakage=True,
                    )
                self.assertTrue(all(path.is_file() for path in fixture.archive_paths))
        finally:
            import shutil
            shutil.rmtree(other_root)

    def test_eviction_is_forbidden_until_consumer_gates_accept_both_states(self):
        other_root = Path(tempfile.mkdtemp(dir=self.root.parent))
        try:
            fixture = Fixture(other_root, integrated=False)
            with patch.object(GUARD, "EXPECTED_ARCHIVE_BYTES", fixture.archive_bytes):
                plan = GUARD.create_or_verify_plan(other_root)
                self.assertEqual(
                    plan["downstream_dependency_audit"]["eviction_compatibility"]["state"],
                    "BLOCKED_ACTIVE_CONSUMER_GATES_DO_NOT_ACCEPT_FINAL_RECEIPT",
                )
                with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "consumer gates"):
                    GUARD.execute_eviction(
                        root=other_root,
                        plan_id=GUARD.digest_json(plan),
                        confirmation=GUARD.CONFIRMATION,
                        acknowledge_operational_breakage=True,
                )
                self.assertTrue(all(path.is_file() for path in fixture.archive_paths))
                self.assertFalse((other_root / GUARD.intent_path(GUARD.digest_json(plan))).exists())
        finally:
            import shutil
            shutil.rmtree(other_root)

    def test_final_receipt_detects_extracted_payload_tampering(self):
        plan = GUARD.create_or_verify_plan(self.root)
        plan_id = GUARD.digest_json(plan)
        GUARD.execute_eviction(
            root=self.root,
            plan_id=plan_id,
            confirmation=GUARD.CONFIRMATION,
            acknowledge_operational_breakage=True,
        )
        self.fixture.payload_paths[-1].write_bytes(b"changed")
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "payload differs"):
            GUARD.verify_reference_state(self.root)

    def test_family_wide_hardlink_preflight_happens_before_intent_or_unlink(self):
        plan = GUARD.create_or_verify_plan(self.root)
        plan_id = GUARD.digest_json(plan)
        alias = self.root / GUARD.REFERENCE_ROOT / "archives/archive-hardlink-alias"
        os.link(self.fixture.archive_paths[-1], alias)
        with self.assertRaisesRegex(GUARD.ArchiveEvictionError, "private regular file"):
            GUARD.execute_eviction(
                root=self.root,
                plan_id=plan_id,
                confirmation=GUARD.CONFIRMATION,
                acknowledge_operational_breakage=True,
            )
        self.assertTrue(all(path.is_file() for path in self.fixture.archive_paths))
        self.assertFalse((self.root / GUARD.intent_path(plan_id)).exists())

    def test_public_content_addressed_lifecycle_is_restart_safe_and_self_bound(self):
        plan = GUARD.create_or_verify_plan(self.root)
        plan_id = GUARD.digest_json(plan)
        self.assertRegex(plan_id, r"^[0-9a-f]{64}$")
        self.assertEqual(
            GUARD._load_plan(self.root, plan_id, require_live_environment=True), plan,
        )
        receipt = GUARD.execute_eviction(
            root=self.root,
            plan_id=plan_id,
            confirmation=GUARD.CONFIRMATION,
            acknowledge_operational_breakage=True,
        )
        intent = json.loads((self.root / GUARD.intent_path(plan_id)).read_text())
        self.assertEqual(intent["plan_id"], plan_id)
        self.assertEqual(receipt["plan_id"], plan_id)
        self.assertEqual(GUARD.verify_final_receipt(self.root, plan_id), receipt)
        # An idempotent retry verifies the existing final receipt; it cannot
        # publish or unlink a second lineage.
        self.assertEqual(
            GUARD.execute_eviction(
                root=self.root,
                plan_id=plan_id,
                confirmation=GUARD.CONFIRMATION,
                acknowledge_operational_breakage=True,
            ),
            receipt,
        )

    def test_plan_can_bind_future_consumers_without_changing_core_script_set(self):
        relative = Path("scripts/147_post_lava_consumer.py")
        (self.root / relative).write_text(
            'GUARD_PATH = "scripts/146_manage_lava_reference_archives.py"\n'
            f'ACCEPTED_ARCHIVE_STATES = {{"{GUARD.ARCHIVES_PRESENT}", '
            f'"{GUARD.ARCHIVES_EVICTED}"}}\n'
            "def verify_archive_state(guard):\n"
            "    return guard.verify_reference_state(None)\n",
            encoding="utf-8",
        )
        plan = GUARD.create_or_verify_plan(self.root, [relative])
        audit = plan["downstream_dependency_audit"]["additional_consumer_audit"]
        self.assertEqual(audit["paths"], [str(relative)])
        self.assertTrue(audit["integrations"][str(relative)])
        self.assertEqual(
            audit["script_identities"][str(relative)]["sha256"], sha(self.root / relative),
        )


class TrackBLavaReferenceArchiveRealStateSmokeTest(unittest.TestCase):
    def test_importable_two_state_api_on_real_repository(self):
        # This catches stale module-level lifecycle names and exercises the real
        # manifests, immutable provenance, seven archives and all 44 payloads.
        state = GUARD.verify_reference_state(REPO)
        self.assertIn(state["state"], {GUARD.ARCHIVES_PRESENT, GUARD.ARCHIVES_EVICTED})
        self.assertEqual(state["extracted_payload_count"], 44)


if __name__ == "__main__":
    unittest.main()
