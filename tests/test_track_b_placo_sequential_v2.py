#!/usr/bin/env python3

from __future__ import annotations

from contextlib import nullcontext
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load_coordinator():
    path = ROOT / "scripts/143_run_track_b_placo_sequential_v2.py"
    specification = importlib.util.spec_from_file_location("track_b_placo_sequential_v2_test", path)
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load sequential PLACO V2 coordinator")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


COORDINATOR = load_coordinator()


class TrackBPLACOSequentialV2Tests(unittest.TestCase):
    def test_frozen_and_additive_predecessors_remain_byte_identical(self) -> None:
        expected = {
            "scripts/123_track_b_pleiotropy_contract.py": "7001146c8548cbb2a6d117f3b5bc7b247d22b3b307b2b67d4cb020f7d64830a2",
            "scripts/124_build_track_b_pleiotropy_input_gate.py": "72dbf1f51a54064247bfe1f2ed06b5ec12ace18fc19445514357a093855c3e56",
            "scripts/125_materialize_track_b_placo_pair.py": "0ec584da9823ba3af8e36e663bb7e10cc6da694e7b95162c27fb1b64e95e9250",
            "scripts/126_run_track_b_placo_pair.R": "d418e36acee10cea4d552beb8feb0c86057a0dc2b2be8bfb436eee08715bf12c",
            "scripts/127_prepare_track_b_pleiotropy_ld.py": "f976ffeb59d4de64291e7909b2f12c2fc2e158a5037cee30b63d136e5cd96b9d",
            "scripts/139_build_track_b_placo_terminal_gate_v2.py": "ebf0ea3bbf0909fe68eeddecd5b05e6b5f02ccc27e53fa432206afd6fa43dca8",
            "scripts/140_materialize_track_b_placo_pair_v2.py": "f47044f53462acd1a2e421a9346977b067d39e074aeaf58f8d3065e9042a0094",
            "scripts/141_run_track_b_placo_pair_v2.R": "c16e444a201139495261dece4656ded10a06b52659478d5e6d335afb1eb301a7",
        }
        for relative, wanted in expected.items():
            with self.subTest(relative=relative):
                observed = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
                self.assertEqual(observed, wanted)

    def test_contract_freezes_sequential_science_and_no_ld(self) -> None:
        gate = {
            "analysis_id": "track-b-v1.0-pleiotropy",
            "resource_envelope": {"ram_safety_reserve_bytes": 1024**3},
        }

        def fake_identity(path: Path) -> dict[str, object]:
            return {"path": str(path), "bytes": 1, "sha256": "a" * 64,
                    "device": 1, "inode": 1}

        reference_state = {
            "state": COORDINATOR.ARCHIVES_PRESENT,
            "reference_provenance_identity": {
                "path": "ref/lava/ukb_v1.1/reference.provenance.json",
                "bytes": 1,
                "sha256": "b" * 64,
            },
            "archive_family_sha256": "c" * 64,
            "extracted_payload_count": 44,
            "extracted_family_sha256": "d" * 64,
        }
        with mock.patch.object(COORDINATOR, "stable_identity", side_effect=fake_identity):
            payload = COORDINATOR.execution_contract_payload(gate, reference_state)
        self.assertEqual(payload["pair_order"], ["A", "B", "CONTROL"])
        self.assertEqual(payload["workers"], 1)
        self.assertEqual(payload["benchmark_variants"], 2000)
        self.assertEqual(payload["ram_limit_bytes"], 8 * 1024**3)
        self.assertIn("ALL_VALID_ALIGNED_GENOME_WIDE_VARIANTS", payload["nuisance_scope"])
        self.assertIn("COMPLETE_FIXED_ALIGNED_PAIR_FAMILY", payload["within_pair_bh"])
        self.assertFalse(payload["ld"]["extraction_during_full_p_scans"])
        self.assertIn("BLOCKED", payload["ld"]["locus_publication"])
        self.assertIn("600_BYTES", payload["storage"]["materialization_admission"])
        self.assertIn("8_GIB", payload["storage"]["post_materialization_runner_floor"])
        self.assertEqual(payload["ram_benchmark_fields"], COORDINATOR.RAM_FIELDS)
        self.assertEqual(
            payload["lava_reference_archive_state_guard"]["extracted_payload_count"], 44,
        )

    def test_coordinator_revalidates_both_archive_states_against_terminal_gate(self) -> None:
        base = {
            "reference_provenance_identity": {
                "path": "ref/lava/ukb_v1.1/reference.provenance.json",
                "bytes": 1,
                "sha256": "a" * 64,
            },
            "archive_family_sha256": "b" * 64,
            "extracted_payload_count": 44,
            "extracted_family_sha256": "c" * 64,
        }

        def fake_identity(_path: Path) -> dict[str, object]:
            return {
                "path": "scripts/146_manage_lava_reference_archives.py",
                "bytes": 1,
                "sha256": "d" * 64,
                "device": 1,
                "inode": 2,
            }

        with mock.patch.object(COORDINATOR, "stable_identity", side_effect=fake_identity):
            frozen = COORDINATOR.reference_archive_state_invariant({
                **base, "state": COORDINATOR.ARCHIVES_PRESENT,
            })
        gate = {"lava_terminal": {"reference_archive_state_guard": frozen}}
        for state in (COORDINATOR.ARCHIVES_PRESENT, COORDINATOR.ARCHIVES_EVICTED):
            with (
                mock.patch.object(
                    COORDINATOR.ARCHIVE_STATE_GUARD, "verify_reference_state",
                    return_value={**base, "state": state},
                ),
                mock.patch.object(COORDINATOR, "stable_identity", side_effect=fake_identity),
            ):
                self.assertEqual(
                    COORDINATOR.verify_reference_archive_state(gate)["state"], state,
                )
        with (
            mock.patch.object(
                COORDINATOR.ARCHIVE_STATE_GUARD, "verify_reference_state",
                return_value={**base, "state": "PARTIAL_ARCHIVE_EVICTION_RESTART_REQUIRED"},
            ),
            self.assertRaisesRegex(COORDINATOR.CoordinatorError, "non-permitted state"),
        ):
            COORDINATOR.verify_reference_archive_state(gate)

    def test_ram_admission_is_capped_at_eight_gib_even_on_larger_host(self) -> None:
        gate = {"resource_envelope": {"ram_safety_reserve_bytes": 1024**3}}
        safe = {"workers": "1", "peak_process_tree_rss_bytes": str(7 * 1024**3)}
        admitted = COORDINATOR.admit_ram(safe, gate, physical_memory_bytes=64 * 1024**3)
        self.assertEqual(admitted["limit_bytes"], 8 * 1024**3)
        unsafe = {"workers": "1", "peak_process_tree_rss_bytes": str(7 * 1024**3 + 1)}
        with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "additional_required=1"):
            COORDINATOR.admit_ram(unsafe, gate, physical_memory_bytes=64 * 1024**3)

    def test_ram_admission_rejects_worker_mismatch_and_uses_smaller_host(self) -> None:
        gate = {"resource_envelope": {"ram_safety_reserve_bytes": 1024**3}}
        with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "one-worker"):
            COORDINATOR.admit_ram(
                {"workers": "2", "peak_process_tree_rss_bytes": str(1024**3)},
                gate, physical_memory_bytes=8 * 1024**3,
            )
        with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "BLOCKED_BY_COMPUTE"):
            COORDINATOR.admit_ram(
                {"workers": "1", "peak_process_tree_rss_bytes": str(4 * 1024**3)},
                gate, physical_memory_bytes=4 * 1024**3,
            )

    def test_storage_forecast_quantifies_exact_shortfall(self) -> None:
        context = {"sources": [{"rows": 10}, {"rows": 20}]}
        required = COORDINATOR.ENGINE.MINIMUM_FREE_BYTES + 600 * 30
        row = COORDINATOR.storage_forecast(context, required - 123)
        self.assertEqual(row["required_bytes"], required)
        self.assertEqual(row["additional_bytes_required"], 123)
        self.assertEqual(row["headroom_bytes"], 0)

    def test_post_materialization_runner_uses_frozen_eight_gib_floor(self) -> None:
        context = {"pair_id": "A"}
        paths = {"task": Path("task.tsv")}
        with (
            mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
            mock.patch.object(
                COORDINATOR.ENGINE, "one_tsv_row",
                return_value={"minimum_free_bytes": str(8 * 1024**3)},
            ),
            mock.patch.object(COORDINATOR.shutil, "disk_usage") as usage,
        ):
            usage.return_value.free = 8 * 1024**3 - 9
            with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "additional_required=9"):
                COORDINATOR.require_runner_storage(context)
            usage.return_value.free = 8 * 1024**3
            observed = COORDINATOR.require_runner_storage(context)
        self.assertEqual(observed["required_bytes"], 8 * 1024**3)

    def test_cleanup_allowlist_is_narrow(self) -> None:
        allowed = (
            "aligned.tsv.gz",
            "materialization.provenance.json",
            "checkpoints/nuisance.rds",
            "checkpoints/shard_000001.rds",
            "checkpoints/shard_000001.rds.sha256.tsv",
            ".aligned.123.tmp.gz",
        )
        for relative in allowed:
            with self.subTest(relative=relative):
                self.assertTrue(COORDINATOR.allowed_cleanup_relative(relative))
        forbidden = (
            "../../data/harmonized/insomnia.harmonized.tsv.gz",
            "canonical.full.tsv.gz",
            "results/track_b/07_placo_plus_variants.tsv",
            "checkpoints/arbitrary.rds",
            "g1000_eur.bed",
        )
        for relative in forbidden:
            with self.subTest(relative=relative):
                self.assertFalse(COORDINATOR.allowed_cleanup_relative(relative))

    def test_cleanup_inventory_fails_closed_on_unknown_file_and_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / "work/track_b_pleiotropy/placo/fingerprint/A"
            base.mkdir(parents=True)
            (base / "aligned.tsv.gz").write_bytes(b"aligned")
            (base / "mystery.bin").write_bytes(b"unknown")
            with mock.patch.object(COORDINATOR, "ROOT", root):
                with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "unknown"):
                    COORDINATOR.collect_cleanup_candidates(base)
            (base / "mystery.bin").unlink()
            outside = root / "outside"
            outside.write_bytes(b"outside")
            (base / "bad-link").symlink_to(outside)
            with mock.patch.object(COORDINATOR, "ROOT", root):
                with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "non-regular"):
                    COORDINATOR.collect_cleanup_candidates(base)

    def test_no_replace_publication_rolls_back_its_exact_inode_on_fsync_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / "results/receipt.json"
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR, "fsync_directory", side_effect=OSError("fsync failed")),
            ):
                with self.assertRaisesRegex(OSError, "fsync failed"):
                    COORDINATOR.publish_bytes_no_replace(destination, b"immutable\n")
            self.assertFalse(destination.exists())

    def test_ram_provenance_can_restore_missing_tsv_after_process_death(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gate = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            gate.parent.mkdir(parents=True)
            gate.write_text("{\"gate\":true}\n", encoding="utf-8")
            task = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / ("f" * 64) / "A/task.tsv"
            task.parent.mkdir(parents=True)
            task.write_text("task\nrow\n", encoding="utf-8")
            context = {
                "pair_id": "A", "fingerprint": "f" * 64,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value={"task": task}),
                mock.patch.object(COORDINATOR.ENGINE, "one_tsv_row", return_value={"aligned_input_rows": "123"}),
                mock.patch.object(COORDINATOR, "rebuild_ram_aggregate"),
            ):
                COORDINATOR.record_ram(
                    context, "full_scan", analysis="scan", peak_bytes=1024**3,
                    runtime_seconds=12.5, exit_status=0, output_hash="a" * 64,
                    rss_measurement_method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
                )
                receipt = COORDINATOR.ram_receipt_path(context, "full_scan")
                receipt.unlink()
                self.assertFalse(receipt.exists())
                loaded = COORDINATOR.load_ram_measurement(context, "full_scan")
                self.assertIsNotNone(loaded)
                self.assertTrue(receipt.is_file())

    def test_ram_provenance_rejects_exact_peak_row_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gate = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            gate.parent.mkdir(parents=True)
            gate.write_text("{\"gate\":true}\n", encoding="utf-8")
            fingerprint = "d" * 64
            task = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / "A/task.tsv"
            task.parent.mkdir(parents=True)
            task.write_text("task\nrow\n", encoding="utf-8")
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value={"task": task}),
                mock.patch.object(
                    COORDINATOR.ENGINE, "one_tsv_row",
                    return_value={"aligned_input_rows": "123"},
                ),
                mock.patch.object(COORDINATOR, "rebuild_ram_aggregate"),
            ):
                COORDINATOR.record_ram(
                    context, "full_scan", analysis="scan", peak_bytes=1024**3,
                    runtime_seconds=2.0, exit_status=0, output_hash="a" * 64,
                    rss_measurement_method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
                )
                provenance_path = COORDINATOR.ram_receipt_path(
                    context, "full_scan",
                ).with_suffix(".provenance.json")
                provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
                provenance["peak_process_tree_rss_bytes"] = 2 * 1024**3
                provenance_path.write_text(
                    json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8",
                )
                with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "TSV row disagrees"):
                    COORDINATOR.load_ram_measurement(context, "full_scan")

    def test_ram_aggregate_federates_three_pair_specific_fingerprints(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gate = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            contract = root / COORDINATOR.EXECUTION_CONTRACT
            gate.parent.mkdir(parents=True)
            gate.write_text("{\"gate\":true}\n", encoding="utf-8")
            contract.write_text("{\"contract\":true}\n", encoding="utf-8")
            contexts = []
            tasks: dict[str, Path] = {}
            for pair, digit in zip(COORDINATOR.PAIR_ORDER, "abc", strict=True):
                fingerprint = digit * 64
                task = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / pair / "task.tsv"
                task.parent.mkdir(parents=True)
                task.write_text("task\nrow\n", encoding="utf-8")
                tasks[pair] = task
                contexts.append({
                    "pair_id": pair, "fingerprint": fingerprint,
                    "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
                })

            def paths(context):
                return {"task": tasks[context["pair_id"]]}

            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", side_effect=paths),
                mock.patch.object(
                    COORDINATOR.ENGINE, "one_tsv_row",
                    return_value={"aligned_input_rows": "99"},
                ),
            ):
                with mock.patch.object(COORDINATOR, "rebuild_ram_aggregate"):
                    for index, context in enumerate(contexts, start=1):
                        COORDINATOR.record_ram(
                            context, "matched_worker_benchmark", analysis=f"scan-{index}",
                            peak_bytes=index * 1024**3, runtime_seconds=float(index),
                            exit_status=0, output_hash=f"{index:x}" * 64,
                            rss_measurement_method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
                        )
                COORDINATOR.rebuild_ram_aggregate()
                fields, rows = COORDINATOR.ENGINE.read_tsv(root / COORDINATOR.RAM_AGGREGATE)
                provenance = json.loads(
                    (root / COORDINATOR.RAM_AGGREGATE_PROVENANCE).read_text(encoding="utf-8")
                )
            self.assertEqual(fields, COORDINATOR.RAM_FIELDS)
            self.assertEqual([row["pair"] for row in rows], list(COORDINATOR.PAIR_ORDER))
            self.assertEqual(
                provenance["execution_fingerprints_by_pair"],
                {context["pair_id"]: context["fingerprint"] for context in contexts},
            )

    def test_every_documented_crash_temporary_including_empty_is_audited(self) -> None:
        relatives = [
            ".aligned.1.tmp.gz", ".provenance.2.tmp.json", ".task.3.tmp.tsv",
            ".benchmark.tsv.4.tmp", ".staged.full.tsv.gz.5.tmp",
            ".staged.provenance.json.6.tmp", "benchmark.raw.tsv.7.tmp",
            "staged.full.tsv.gz.8.tmp", "run.summary.tsv.9.tmp",
            "checkpoints/nuisance.rds.10.tmp",
            "checkpoints/nuisance.sha256.tsv.11.tmp",
            "checkpoints/shard_000001.rds.12.tmp",
            "checkpoints/shard_000001.rds.sha256.tsv.13.tmp",
            ".staged.full.tsv.gz.14.validate.sqlite",
            ".staged.full.tsv.gz.14.validate.sqlite-journal",
            ".staged.full.tsv.gz.14.validate.sqlite-wal",
            ".staged.full.tsv.gz.14.validate.sqlite-shm",
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "e" * 64
            base = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / "A"
            (base / "checkpoints").mkdir(parents=True)
            scientific = base / "checkpoints/shard_000002.rds"
            scientific.write_bytes(b"scientific checkpoint")
            for index, relative in enumerate(relatives):
                path = base / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"" if index % 2 == 0 else b"crash bytes")
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value={"base": base}),
            ):
                COORDINATOR.repair_crash_temporaries(context)
            self.assertTrue(scientific.is_file())
            self.assertTrue(all(not (base / relative).exists() for relative in relatives))
            receipts = list(
                (root / COORDINATOR.RECORD_ROOT / fingerprint / "A/repairs").glob(
                    "documented_crash_temporaries.*.json"
                )
            )
            self.assertEqual(len(receipts), 1)
            payload = json.loads(receipts[0].read_text(encoding="utf-8"))
            self.assertEqual(len(payload["removed"]), len(relatives))

    def test_live_rss_guard_terminates_owned_process_group(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR, "process_group_rss_bytes", return_value=2),
            ):
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    measurement = COORDINATOR.monitored_command(
                        [sys.executable, "-c", "import time; time.sleep(30)"],
                        maximum_process_tree_rss_bytes=1,
                    )
            self.assertNotEqual(measurement["exit_status"], 0)
            self.assertEqual(measurement["termination_reason"], "RAM_LIMIT_EXCEEDED")
            self.assertEqual(measurement["peak_bytes"], 2)
            self.assertIn(b"TRACK_B_PLACO_LIVE_RSS_MONITOR_TERMINATED", measurement["log"])

    def test_worker_lease_payload_binds_pair_run_and_inode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = {"pair_id": "A", "fingerprint": "a" * 64}
            with mock.patch.object(COORDINATOR, "ROOT", root):
                with COORDINATOR.worker_lease(context):
                    active = COORDINATOR.require_active_worker_lease(context)
                    payload = json.loads(
                        COORDINATOR.base64.b64decode(active["payload_b64"]).decode("utf-8")
                    )
                    self.assertEqual(payload["pair_id"], "A")
                    self.assertEqual(payload["run_fingerprint"], "a" * 64)
                    self.assertEqual(payload["device"], active["device"])
                    self.assertEqual(payload["inode"], active["inode"])
                    environment = COORDINATOR.inherited_worker_environment({
                        "PYTHONWARNINGS": "error::ResourceWarning",
                        "PYTHONDEVMODE": "1",
                        "TRACK_B_UNRELATED": "preserved",
                    })
                    self.assertNotIn("PYTHONWARNINGS", environment)
                    self.assertNotIn("PYTHONDEVMODE", environment)
                    self.assertEqual(environment["TRACK_B_UNRELATED"], "preserved")
                    self.assertEqual(
                        environment[COORDINATOR.WORKER_LEASE_DIGEST_ENV],
                        active["payload_sha256"],
                    )
                    with self.assertRaisesRegex(
                        COORDINATOR.CoordinatorError, "different pair or run fingerprint",
                    ):
                        COORDINATOR.require_active_worker_lease(
                            {"pair_id": "A", "fingerprint": "b" * 64},
                        )

    def test_record_temp_repair_requires_the_affected_pair_lease(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context_a = {"pair_id": "A", "fingerprint": "a" * 64}
            context_b = {"pair_id": "B", "fingerprint": "b" * 64}
            with mock.patch.object(COORDINATOR, "ROOT", root):
                with COORDINATOR.worker_lease(context_a):
                    with self.assertRaisesRegex(
                        COORDINATOR.CoordinatorError, "different pair or run fingerprint",
                    ):
                        COORDINATOR.repair_record_publication_temporaries(context_b)

    @unittest.skipUnless(hasattr(os, "fork"), "requires POSIX fork semantics")
    def test_worker_lease_survives_coordinator_sigkill_until_child_exits(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            child_pid_file = root / "child.pid"
            context = {"pair_id": "A", "fingerprint": "a" * 64}
            with mock.patch.object(COORDINATOR, "ROOT", root):
                supervisor = os.fork()
                if supervisor == 0:  # pragma: no cover - asserted by parent process
                    try:
                        with COORDINATOR.worker_lease(context):
                            child = subprocess.Popen(
                                [sys.executable, "-c", "import time; time.sleep(30)"],
                                pass_fds=COORDINATOR.inherited_worker_lease_fds(),
                                env=COORDINATOR.inherited_worker_environment(),
                            )
                            child_pid_file.write_text(str(child.pid), encoding="ascii")
                            with child_pid_file.open("rb") as handle:
                                os.fsync(handle.fileno())
                            os.kill(os.getpid(), 9)
                    finally:
                        os._exit(91)
                _, status = os.waitpid(supervisor, 0)
                self.assertTrue(os.WIFSIGNALED(status))
                deadline = time.monotonic() + 3
                while not child_pid_file.exists() and time.monotonic() < deadline:
                    time.sleep(0.02)
                child_pid = int(child_pid_file.read_text(encoding="ascii"))
                try:
                    with self.assertRaisesRegex(
                        COORDINATOR.CoordinatorError, "active or orphaned",
                    ):
                        with COORDINATOR.worker_lease(context):
                            pass
                finally:
                    try:
                        os.kill(child_pid, 9)
                    except ProcessLookupError:
                        pass
                acquired = False
                deadline = time.monotonic() + 3
                while not acquired and time.monotonic() < deadline:
                    try:
                        with COORDINATOR.worker_lease(context):
                            acquired = True
                    except COORDINATOR.CoordinatorError:
                        time.sleep(0.02)
                self.assertTrue(acquired, "orphan lease did not release after exact child exit")

    def test_r_runner_rejects_wrong_run_lease_before_scientific_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "scripts").mkdir(parents=True)
            (root / "scripts/143_run_track_b_placo_sequential_v2.py").write_bytes(
                (ROOT / "scripts/143_run_track_b_placo_sequential_v2.py").read_bytes()
            )
            task_path = root / "task.tsv"
            row = {field: "x" for field in COORDINATOR.ENGINE.TASK_FIELDS}
            row.update({
                "analysis_id": "track-b-v1.0-pleiotropy",
                "pair_id": "A", "trait1": "snoring", "trait2": "parental_lifespan",
                "family_role": "PRIMARY_DISCOVERY", "run_fingerprint": "b" * 64,
            })
            task_path.write_bytes(COORDINATOR.tsv_bytes(COORDINATOR.ENGINE.TASK_FIELDS, [row]))
            rscript = ROOT / ".r-env/bin/Rscript"
            runner = ROOT / "scripts/141_run_track_b_placo_pair_v2.R"
            with mock.patch.object(COORDINATOR, "ROOT", root):
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    result = subprocess.run(
                        [str(rscript), str(runner), str(task_path), "--root", str(root), "--execute"],
                        check=False, capture_output=True, text=True,
                        pass_fds=COORDINATOR.inherited_worker_lease_fds(),
                        env=COORDINATOR.inherited_worker_environment(),
                    )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not authorized for this exact pair/run and inode", result.stderr)
            row["run_fingerprint"] = "a" * 64
            task_path.write_bytes(COORDINATOR.tsv_bytes(COORDINATOR.ENGINE.TASK_FIELDS, [row]))
            with mock.patch.object(COORDINATOR, "ROOT", root):
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    matched = subprocess.run(
                        [str(rscript), str(runner), str(task_path), "--root", str(root), "--execute"],
                        check=False, capture_output=True, text=True,
                        pass_fds=COORDINATOR.inherited_worker_lease_fds(),
                        env=COORDINATOR.inherited_worker_environment(),
                    )
            self.assertNotEqual(matched.returncode, 0)
            self.assertNotIn("not authorized for this exact pair/run and inode", matched.stderr)

    def test_process_group_rss_includes_reparented_group_members(self) -> None:
        listing = "10 77 100\n11 77 250\n12 88 999\n"
        completed = subprocess.CompletedProcess([], 0, stdout=listing, stderr="")
        with mock.patch.object(COORDINATOR.subprocess, "run", return_value=completed):
            self.assertEqual(COORDINATOR.process_group_rss_bytes(77), 350 * 1024)

    def test_unmeasurable_live_rss_terminates_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR, "process_group_rss_bytes", return_value=0),
            ):
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    measurement = COORDINATOR.monitored_command(
                        [sys.executable, "-c", "import time; time.sleep(30)"],
                        maximum_process_tree_rss_bytes=7 * 1024**3,
                    )
            self.assertNotEqual(measurement["exit_status"], 0)
            self.assertEqual(
                measurement["termination_reason"], "RSS_SAMPLER_RETURNED_NONPOSITIVE",
            )

    def test_matched_benchmark_guard_uses_owned_process_group_and_records_failure(self) -> None:
        limit = 7 * 1024**3
        captured: dict[str, object] = {}

        def fake_benchmark(_context, _variants, _workers):
            process = COORDINATOR.ENGINE.subprocess.Popen(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                stdout=COORDINATOR.subprocess.PIPE,
                stderr=COORDINATOR.subprocess.PIPE,
                text=True,
            )
            COORDINATOR.ENGINE.process_tree_rss_bytes(process.pid)
            process.communicate()
            raise SystemExit("benchmark killed")

        with tempfile.TemporaryDirectory() as directory:
            with (
                mock.patch.object(COORDINATOR, "ROOT", Path(directory)),
                mock.patch.object(COORDINATOR.ENGINE, "physical_memory_bytes", return_value=8 * 1024**3),
                mock.patch.object(COORDINATOR, "process_group_rss_bytes", return_value=limit + 1),
                mock.patch.object(COORDINATOR.ENGINE, "benchmark_pair", side_effect=fake_benchmark),
                mock.patch.object(COORDINATOR, "benchmark_failure_measurement") as failure,
            ):
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    with self.assertRaisesRegex(SystemExit, "benchmark killed"):
                        COORDINATOR.run_benchmark_with_live_guard(
                            {}, {"resource_envelope": {"ram_safety_reserve_bytes": 1024**3}},
                        )
        captured.update(failure.call_args.kwargs)
        self.assertEqual(captured["termination_reason"], "RAM_LIMIT_EXCEEDED")
        self.assertEqual(captured["guard_limit_bytes"], limit)
        self.assertTrue(captured["process_started"])

    def test_benchmark_leader_exit_cannot_block_on_descendant_held_output(self) -> None:
        captured: dict[str, object] = {}

        def fake_benchmark(_context, _variants, _workers):
            leader = (
                "import subprocess,sys; "
                "child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)']); "
                "print(f'descendant_pid={child.pid}', flush=True)"
            )
            process = COORDINATOR.ENGINE.subprocess.Popen(
                [sys.executable, "-c", leader],
                stdout=COORDINATOR.ENGINE.subprocess.PIPE,
                stderr=COORDINATOR.ENGINE.subprocess.PIPE,
                text=True,
            )
            while process.poll() is None:
                COORDINATOR.ENGINE.process_tree_rss_bytes(process.pid)
                time.sleep(0.01)
            # The production V1 engine has no timeout here.  A pipe would wait
            # for the descendant's inherited writer; the coordinator facade
            # must instead terminate the group and fail before any publication.
            process.communicate()
            return {"unexpected": "success"}

        with tempfile.TemporaryDirectory() as directory:
            with (
                mock.patch.object(COORDINATOR, "ROOT", Path(directory)),
                mock.patch.object(COORDINATOR.ENGINE, "physical_memory_bytes", return_value=8 * 1024**3),
                mock.patch.object(COORDINATOR, "process_group_rss_bytes", return_value=1),
                mock.patch.object(COORDINATOR.ENGINE, "benchmark_pair", side_effect=fake_benchmark),
                mock.patch.object(COORDINATOR, "benchmark_failure_measurement") as failure,
            ):
                started = time.monotonic()
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    with self.assertRaisesRegex(
                        COORDINATOR.CoordinatorError, "leader exited while its owned process group remained",
                    ):
                        COORDINATOR.run_benchmark_with_live_guard(
                            {}, {"resource_envelope": {"ram_safety_reserve_bytes": 1024**3}},
                        )
                elapsed = time.monotonic() - started
        captured.update(failure.call_args.kwargs)
        self.assertLess(elapsed, 3.0)
        self.assertEqual(
            captured["termination_reason"],
            "BENCHMARK_LEADER_EXITED_WITH_OWNED_DESCENDANTS",
        )
        self.assertRegex(captured["process_log"].decode(), r"descendant_pid=[0-9]+")

    def test_unattested_complete_benchmark_is_audited_before_restart(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "d" * 64
            base = root / "work/track_b_pleiotropy/placo" / fingerprint / "A"
            checkpoint = base / "checkpoints"
            checkpoint.mkdir(parents=True)
            paths = {
                "base": base,
                "benchmark_raw": base / "benchmark.raw.tsv",
                "benchmark": base / "benchmark.tsv",
                "staged_ledger": base / "staged.full.tsv.gz",
                "run_summary": base / "run.summary.tsv",
                "staged_provenance": base / "staged.provenance.json",
                "checkpoint_dir": checkpoint,
            }
            for path in (
                paths["benchmark_raw"], paths["benchmark"],
                checkpoint / "nuisance.rds", checkpoint / "nuisance.sha256.tsv",
            ):
                path.write_bytes(b"unattested benchmark lineage\n")
            context = {"pair_id": "A", "fingerprint": fingerprint}
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
                mock.patch.object(COORDINATOR, "require_runner_storage"),
                mock.patch.object(
                    COORDINATOR, "run_benchmark_with_live_guard",
                    side_effect=COORDINATOR.CoordinatorError("fresh monitored rerun reached"),
                ) as rerun,
            ):
                with self.assertRaisesRegex(
                    COORDINATOR.CoordinatorError, "fresh monitored rerun reached",
                ):
                    COORDINATOR.ensure_benchmark(context, {})
            rerun.assert_called_once()
            self.assertFalse(any(path.exists() for path in (
                paths["benchmark_raw"], paths["benchmark"],
                checkpoint / "nuisance.rds", checkpoint / "nuisance.sha256.tsv",
            )))
            receipts = list(
                (root / COORDINATOR.RECORD_ROOT / fingerprint / "A/repairs").glob(
                    "unadmitted_benchmark_lineage.*.json"
                )
            )
            self.assertEqual(len(receipts), 1)

    def test_benchmark_ps_sampler_cannot_replace_owned_r_process(self) -> None:
        reserve = 8 * 1024**3 - 1

        def nested_ps_sampler(pid):
            COORDINATOR.ENGINE.subprocess.run(
                [sys.executable, "-c", "pass"], check=False, capture_output=True, text=True,
            )
            return 2

        def fake_benchmark(_context, _variants, _workers):
            process = COORDINATOR.ENGINE.subprocess.Popen(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                stdout=COORDINATOR.ENGINE.subprocess.PIPE,
                stderr=COORDINATOR.ENGINE.subprocess.PIPE,
                text=True,
            )
            try:
                # This is the real engine sampler, which invokes
                # ENGINE.subprocess.run(["ps", ...]).  That nested ps launch
                # must not replace ownership of the long-lived R stand-in.
                COORDINATOR.ENGINE.process_tree_rss_bytes(process.pid)
                process.communicate(timeout=3)
            except COORDINATOR.subprocess.TimeoutExpired:
                process.kill()
                process.communicate()
                self.fail("live guard terminated ps instead of its owned R process")
            # Adversarially return success even though the sampler crossed the
            # one-byte limit; the outer wrapper must still reject.
            return {"unexpected": "success"}

        with tempfile.TemporaryDirectory() as directory:
            with (
                mock.patch.object(COORDINATOR, "ROOT", Path(directory)),
                mock.patch.object(COORDINATOR.ENGINE, "physical_memory_bytes", return_value=8 * 1024**3),
                mock.patch.object(COORDINATOR, "process_group_rss_bytes", side_effect=nested_ps_sampler),
                mock.patch.object(COORDINATOR.ENGINE, "benchmark_pair", side_effect=fake_benchmark),
                mock.patch.object(COORDINATOR, "benchmark_failure_measurement") as failure,
            ):
                with COORDINATOR.worker_lease({"pair_id": "A", "fingerprint": "a" * 64}):
                    with self.assertRaisesRegex(
                        COORDINATOR.CoordinatorError, "returned after fail-closed termination",
                    ):
                        COORDINATOR.run_benchmark_with_live_guard(
                            {}, {"resource_envelope": {"ram_safety_reserve_bytes": reserve}},
                        )
        self.assertEqual(failure.call_args.kwargs["termination_reason"], "RAM_LIMIT_EXCEEDED")

    def test_over_envelope_benchmark_routes_to_durable_rejection_event(self) -> None:
        paths = {
            "benchmark": mock.MagicMock(), "benchmark_raw": mock.MagicMock(),
            "checkpoint_dir": mock.MagicMock(), "task": mock.MagicMock(),
        }
        paths["benchmark"].exists.return_value = True
        row = {
            "workers": "1", "peak_process_tree_rss_bytes": str(7 * 1024**3 + 1),
            "wrapper_wall_seconds": "1", "rss_measurement_method": "AGGREGATE_PROCESS_TREE_PS_SAMPLED",
            "global_nuisance_checkpoint_reused": "FALSE",
        }
        fields = COORDINATOR.ENGINE.BENCHMARK_RAW_FIELDS + [
            "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
            "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
            "peak_rss_fraction_of_physical_memory",
        ]
        calls: list[str] = []
        with (
            mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
            mock.patch.object(COORDINATOR, "repair_partial_benchmark"),
            mock.patch.object(COORDINATOR.BRIDGE, "validate_benchmark", return_value=row),
            mock.patch.object(
                COORDINATOR, "admit_ram",
                side_effect=COORDINATOR.CoordinatorError("BLOCKED_BY_COMPUTE"),
            ),
            mock.patch.object(COORDINATOR.ENGINE, "read_tsv", return_value=(fields, [row])),
            mock.patch.object(
                COORDINATOR, "record_benchmark_rejection",
                side_effect=lambda *_: calls.append("ram"),
            ),
        ):
            with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "BLOCKED_BY_COMPUTE"):
                COORDINATOR.ensure_benchmark(
                    {"pair_id": "A", "fingerprint": "a" * 64}, {},
                )
        self.assertEqual(calls, ["ram"])

    def test_sigkill_window_reconciles_terminal_failure_event_without_retry(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "e" * 64
            gate = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            contract = root / COORDINATOR.EXECUTION_CONTRACT
            gate.parent.mkdir(parents=True)
            gate.write_text("{\"gate\":true}\n", encoding="utf-8")
            contract.write_text("{\"contract\":true}\n", encoding="utf-8")
            task = (
                root / COORDINATOR.ENGINE.PLACO_WORK_ROOT
                / fingerprint / "A" / "task.tsv"
            )
            task.parent.mkdir(parents=True)
            task.write_text(
                "aligned_input_sha256\taligned_input_rows\n"
                + "a" * 64 + "\t10\n",
                encoding="utf-8",
            )
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(
                    COORDINATOR.ENGINE, "materialized_paths", return_value={"task": task},
                ),
            ):
                # Model a coordinator SIGKILL in the only dangerous window:
                # the self-contained event is durable, but neither its attempt
                # log nor RAM provenance has yet been derived.
                with COORDINATOR.worker_lease(context):
                    with (
                        mock.patch.object(
                            COORDINATOR, "reconcile_failure_event",
                            side_effect=COORDINATOR.CoordinatorError(
                                "simulated coordinator SIGKILL after event commit"
                            ),
                        ),
                        self.assertRaisesRegex(
                            COORDINATOR.CoordinatorError, "simulated coordinator SIGKILL",
                        ),
                    ):
                        COORDINATOR.publish_failure_event(
                            context,
                            failure_class="FULL_SCAN_RAM_ENVELOPE_REJECTION",
                            phase_prefix=COORDINATOR.FULL_SCAN_REJECTION_PREFIX,
                            analysis="PLACO_PLUS_FULL_P_SCAN_RAM_ENVELOPE_REJECTED",
                            attempt_log_label="full_scan_ram_rejected",
                            process_log=b"worker completed but exceeded the frozen cap\n",
                            peak_bytes=7 * 1024**3 + 1,
                            runtime_seconds=2.0,
                            public_exit_status=78,
                            owned_child_exit_status=0,
                            rss_measurement_method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
                            output_hash="b" * 64,
                            termination_reason="RAM_LIMIT_EXCEEDED_AFTER_FULL_SCAN",
                            process_group_termination_status=(
                                "NOT_REQUIRED_PROCESS_ALREADY_EXITED"
                            ),
                            terminal_compute_rejection=True,
                            details={"ram_guard_limit_bytes": 7 * 1024**3},
                        )
                attempts = (
                    root / COORDINATOR.RECORD_ROOT / fingerprint / "A" / "attempts"
                )
                events = list(attempts.glob("failure_event.*.json"))
                self.assertEqual(len(events), 1)
                self.assertFalse(
                    (root / COORDINATOR.RECORD_ROOT / fingerprint / "A" / "ram").exists()
                )

                # A fresh coordinator deterministically completes the exact
                # log/RAM derivation and preserves public rejection status 78
                # separately from the successful owned child status 0.
                with COORDINATOR.worker_lease(context):
                    COORDINATOR.reconcile_failure_events(context)
                ram_root = root / COORDINATOR.RECORD_ROOT / fingerprint / "A" / "ram"
                provenance_paths = list(
                    ram_root.glob("full_scan_ram_rejected_*.provenance.json")
                )
                self.assertEqual(len(provenance_paths), 1)
                provenance = json.loads(provenance_paths[0].read_text(encoding="utf-8"))
                self.assertEqual(provenance["exit_status"], 78)
                self.assertEqual(provenance["output_hash"], "b" * 64)
                event = json.loads(events[0].read_text(encoding="utf-8"))["event"]
                self.assertEqual(event["owned_child_exit_status"], 0)
                self.assertEqual(event["process_group_termination_status"], (
                    "NOT_REQUIRED_PROCESS_ALREADY_EXITED"
                ))
                _, aggregate_rows = COORDINATOR.ENGINE.read_tsv(
                    root / COORDINATOR.RAM_AGGREGATE
                )
                self.assertEqual(aggregate_rows[0]["exit_status"], "78")

                # The reconciled over-cap event is terminal: the full scan is
                # rejected before command construction or worker launch.
                with COORDINATOR.worker_lease(context):
                    with self.assertRaisesRegex(
                        COORDINATOR.CoordinatorError, "BLOCKED_BY_COMPUTE",
                    ):
                        COORDINATOR.run_full_scan_monitored(context, {}, {})

    def test_zero_rss_start_failure_is_valid_immutable_failure_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "7" * 64
            gate = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            gate.parent.mkdir(parents=True)
            gate.write_text("{\"gate\":true}\n", encoding="utf-8")
            task = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / "A/task.tsv"
            task.parent.mkdir(parents=True)
            task.write_text("task\nrow\n", encoding="utf-8")
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value={"task": task}),
                mock.patch.object(
                    COORDINATOR.ENGINE, "one_tsv_row", return_value={"aligned_input_rows": "1"},
                ),
                mock.patch.object(COORDINATOR, "rebuild_ram_aggregate"),
            ):
                row = COORDINATOR.record_ram(
                    context, "benchmark_start_failed", analysis="failed", peak_bytes=0,
                    runtime_seconds=0.01, exit_status=70, output_hash="f" * 64,
                    rss_measurement_method="NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE",
                )
                loaded = COORDINATOR.load_ram_measurement(context, "benchmark_start_failed")
            self.assertEqual(row["peak_ram_gb"], "0")
            self.assertEqual(loaded["peak_process_tree_rss_bytes"], 0)

    def test_complete_stage_without_monitor_receipt_is_audited_before_rebuild(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "6" * 64
            base = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / "A"
            base.mkdir(parents=True)
            staged = base / "staged.full.tsv.gz"
            summary = base / "run.summary.tsv"
            staged.write_bytes(b"complete stage")
            summary.write_bytes(b"complete summary")
            paths = {
                "staged_ledger": staged, "run_summary": summary,
                "staged_provenance": base / "staged.provenance.json",
                "canonical_ledger": root / "results/A.full.tsv.gz",
                "canonical_provenance": root / "results/A.provenance.json",
                "task": base / "task.tsv",
                "checkpoint_dir": base / "checkpoints",
            }
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }

            def audit(_context, _label, targets, _reason):
                for target in targets:
                    target.unlink()

            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
                mock.patch.object(COORDINATOR, "load_ram_measurement", return_value=None),
                mock.patch.object(COORDINATOR, "audit_and_unlink_regenerable", side_effect=audit) as repair,
                mock.patch.object(
                    COORDINATOR, "require_runner_storage",
                    side_effect=COORDINATOR.CoordinatorError("STOP_AFTER_REPAIR"),
                ),
            ):
                with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "STOP_AFTER_REPAIR"):
                    COORDINATOR.run_full_scan_monitored(context, {}, {})
            repair.assert_called_once()
            self.assertFalse(staged.exists())
            self.assertFalse(summary.exists())

    def test_run_order_is_fixed_and_one_repository_lock_spans_family(self) -> None:
        calls: list[tuple[str, str]] = []

        def context(pair: str):
            calls.append(("configure", pair))
            return {"pair_id": pair}, {"gate": True}

        def process(value, _gate, _contract):
            calls.append(("process", value["pair_id"]))
            return "DONE"

        lock = mock.MagicMock(return_value=nullcontext())
        with (
            mock.patch.object(COORDINATOR, "verify_execution_contract", return_value={}),
            mock.patch.object(COORDINATOR.BRIDGE, "execution_lock", lock),
            mock.patch.object(COORDINATOR.BRIDGE, "configure_context", side_effect=context),
            mock.patch.object(COORDINATOR, "worker_lease", return_value=nullcontext()),
            mock.patch.object(COORDINATOR, "repair_record_publication_temporaries"),
            mock.patch.object(COORDINATOR, "reconcile_failure_events"),
            mock.patch.object(COORDINATOR, "process_pair", side_effect=process),
        ):
            observed = COORDINATOR.run_sequential()
        self.assertEqual([row["pair"] for row in observed], ["A", "B", "CONTROL"])
        self.assertEqual(
            calls,
            [
                ("configure", "A"), ("configure", "B"), ("configure", "CONTROL"),
                ("process", "A"), ("process", "B"), ("process", "CONTROL"),
            ],
        )
        lock.assert_called_once_with("sequential-full-p-family", "A,B,CONTROL")

    def test_existing_cleanup_plan_resumes_without_reanalysis(self) -> None:
        context = {"pair_id": "A"}
        plan = mock.MagicMock()
        plan.exists.return_value = True
        plan.is_symlink.return_value = False
        done = mock.MagicMock()
        done.exists.return_value = False
        done.is_symlink.return_value = False
        with (
            mock.patch.object(COORDINATOR, "worker_lease", return_value=nullcontext()),
            mock.patch.object(COORDINATOR, "require_active_worker_lease", return_value={}),
            mock.patch.object(COORDINATOR, "repair_record_publication_temporaries"),
            mock.patch.object(COORDINATOR, "reconcile_failure_events"),
            mock.patch.object(COORDINATOR, "cleanup_bundle_path", return_value=plan),
            mock.patch.object(COORDINATOR, "cleanup_complete_path", return_value=done),
            mock.patch.object(COORDINATOR, "complete_cleanup", return_value={}) as cleanup,
            mock.patch.object(COORDINATOR, "ensure_materialized") as materialize,
            mock.patch.object(COORDINATOR, "run_full_scan_monitored") as scan,
        ):
            outcome = COORDINATOR.process_pair(context, {}, {})
        self.assertEqual(outcome, "RESUMED_CLEANUP_COMPLETE")
        cleanup.assert_called_once()
        materialize.assert_not_called()
        scan.assert_not_called()

    def test_partial_materialization_is_audited_then_regenerated(self) -> None:
        paths = {
            "aligned": mock.MagicMock(), "provenance": mock.MagicMock(),
            "task": mock.MagicMock(), "checkpoint_dir": mock.MagicMock(),
        }
        paths["aligned"].exists.return_value = True
        paths["aligned"].is_symlink.return_value = False
        for key in ("provenance", "task"):
            paths[key].exists.return_value = False
            paths[key].is_symlink.return_value = False
        paths.update({
            key: mock.MagicMock() for key in (
                "benchmark_raw", "benchmark", "staged_ledger", "run_summary",
                "staged_provenance", "canonical_ledger", "canonical_provenance",
            )
        })
        for key in (
            "benchmark_raw", "benchmark", "staged_ledger", "run_summary",
            "staged_provenance", "canonical_ledger", "canonical_provenance",
        ):
            paths[key].exists.return_value = False
            paths[key].is_symlink.return_value = False
        paths["checkpoint_dir"].exists.return_value = False
        context = {"pair_id": "A"}
        with (
            mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
            mock.patch.object(COORDINATOR, "audit_and_unlink_regenerable") as audit,
            mock.patch.object(COORDINATOR, "require_storage", return_value={}),
            mock.patch.object(COORDINATOR.ENGINE, "resource_preflight", return_value={}),
            mock.patch.object(COORDINATOR.ENGINE, "materialize_pair", return_value={"ok": True}) as materialize,
        ):
            observed = COORDINATOR.ensure_materialized(context)
        self.assertEqual(observed, {"ok": True})
        audit.assert_called_once()
        materialize.assert_called_once_with(context)

    def test_partial_canonical_family_never_triggers_replacement(self) -> None:
        context = {"pair_id": "B"}
        plan = mock.MagicMock()
        plan.exists.return_value = False
        plan.is_symlink.return_value = False
        done = mock.MagicMock()
        done.exists.return_value = False
        done.is_symlink.return_value = False
        with (
            mock.patch.object(COORDINATOR, "worker_lease", return_value=nullcontext()),
            mock.patch.object(COORDINATOR, "require_active_worker_lease", return_value={}),
            mock.patch.object(COORDINATOR, "repair_record_publication_temporaries"),
            mock.patch.object(COORDINATOR, "reconcile_failure_events"),
            mock.patch.object(COORDINATOR, "cleanup_bundle_path", return_value=plan),
            mock.patch.object(COORDINATOR, "cleanup_complete_path", return_value=done),
            mock.patch.object(COORDINATOR, "repair_crash_temporaries"),
            mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value={}),
            mock.patch.object(
                COORDINATOR, "reconcile_partial_canonical_publication",
                side_effect=COORDINATOR.CoordinatorError("not an exact crash prefix"),
            ),
            mock.patch.object(COORDINATOR, "ensure_materialized") as materialize,
        ):
            with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "not an exact crash prefix"):
                COORDINATOR.process_pair(context, {}, {})
        materialize.assert_not_called()

    def test_exact_canonical_publication_prefixes_reconcile_after_each_link(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "9" * 64
            base = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / "A"
            canonical_root = root / "results/track_b/pleiotropy/results/placo"
            base.mkdir(parents=True)
            canonical_root.mkdir(parents=True)
            staged_ledger = base / "staged.full.tsv.gz"
            staged_ledger.write_bytes(b"deeply validated complete ledger")
            paths = {
                "staged_ledger": staged_ledger,
                "canonical_ledger": canonical_root / "A.full.tsv.gz",
                "canonical_provenance": canonical_root / "A.provenance.json",
                "staged_provenance": base / "staged.provenance.json",
            }
            gate_path = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            contract_path = root / COORDINATOR.EXECUTION_CONTRACT
            gate_path.parent.mkdir(parents=True, exist_ok=True)
            gate_lock = {"gate": True}
            contract = {"schema_version": COORDINATOR.CONTRACT_SCHEMA}
            gate_path.write_text(json.dumps(gate_lock) + "\n", encoding="utf-8")
            contract_path.write_text(json.dumps(contract) + "\n", encoding="utf-8")
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            benchmark = {"workers": "1", "peak_process_tree_rss_bytes": "1"}
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
                mock.patch.object(
                    COORDINATOR, "validate_admit_fresh_benchmark", return_value=benchmark,
                ),
                mock.patch.object(COORDINATOR, "require_full_scan_ram_evidence", return_value={}),
                mock.patch.object(COORDINATOR.BRIDGE, "revalidate_publication_state"),
                mock.patch.object(COORDINATOR, "validate_completed_pair", return_value={}),
                mock.patch.object(COORDINATOR.GATE, "verify_gate", return_value=gate_lock),
            ):
                # SIGKILL after link 1: only the exact staged-ledger hardlink is legal.
                paths["canonical_ledger"].hardlink_to(staged_ledger)
                self.assertEqual(
                    COORDINATOR.reconcile_partial_canonical_publication(
                        context, gate_lock, contract,
                    ),
                    "ABSENT",
                )
                self.assertFalse(paths["canonical_ledger"].exists())
                self.assertTrue(staged_ledger.exists())

                # SIGKILL after link 2: recover link 3 from the exact canonical
                # provenance inode, then require deep validation.
                paths["canonical_ledger"].hardlink_to(staged_ledger)
                paths["canonical_provenance"].write_bytes(b"exact result provenance")
                self.assertEqual(
                    COORDINATOR.reconcile_partial_canonical_publication(
                        context, gate_lock, contract,
                    ),
                    "COMPLETE",
                )
                self.assertTrue(paths["canonical_provenance"].samefile(paths["staged_provenance"]))

                # SIGKILL after link 3: the complete prefix is left for normal
                # deep validation and is never republished or replaced.
                self.assertEqual(
                    COORDINATOR.reconcile_partial_canonical_publication(
                        context, gate_lock, contract,
                    ),
                    "COMPLETE",
                )

    def test_non_hardlinked_canonical_prefix_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "8" * 64
            base = root / COORDINATOR.ENGINE.PLACO_WORK_ROOT / fingerprint / "A"
            canonical_root = root / "results/track_b/pleiotropy/results/placo"
            base.mkdir(parents=True)
            canonical_root.mkdir(parents=True)
            staged = base / "staged.full.tsv.gz"
            staged.write_bytes(b"same bytes, different inode")
            canonical = canonical_root / "A.full.tsv.gz"
            canonical.write_bytes(staged.read_bytes())
            paths = {
                "staged_ledger": staged, "canonical_ledger": canonical,
                "canonical_provenance": canonical_root / "A.provenance.json",
                "staged_provenance": base / "staged.provenance.json",
            }
            gate_path = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            contract_path = root / COORDINATOR.EXECUTION_CONTRACT
            gate_path.parent.mkdir(parents=True, exist_ok=True)
            gate_path.write_text("{\"gate\":true}\n", encoding="utf-8")
            contract = {"schema_version": COORDINATOR.CONTRACT_SCHEMA}
            contract_path.write_text(json.dumps(contract) + "\n", encoding="utf-8")
            context = {
                "pair_id": "A", "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
                mock.patch.object(
                    COORDINATOR, "validate_admit_fresh_benchmark", return_value={},
                ),
                mock.patch.object(COORDINATOR, "require_full_scan_ram_evidence", return_value={}),
                mock.patch.object(COORDINATOR.BRIDGE, "revalidate_publication_state"),
            ):
                with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "exact staged hardlink"):
                    COORDINATOR.reconcile_partial_canonical_publication(
                        context, {"gate": True}, contract,
                    )
            self.assertTrue(canonical.exists())
            self.assertTrue(staged.exists())

    def test_terminal_gate_exact_prefix_is_restart_safe_and_out_of_order_fails(self) -> None:
        gate = COORDINATOR.GATE
        readiness = "schema\tstatus\nrow\tPASS\n"
        payload = {"schema_version": gate.GATE_SCHEMA, "status": "READY"}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            readiness_path = root / gate.READINESS_GATE
            lock_path = root / gate.TERMINAL_GATE_LOCK
            readiness_path.parent.mkdir(parents=True)
            with (
                mock.patch.object(gate, "ROOT", root),
                mock.patch.object(gate, "build_expected", return_value=(readiness, payload)),
                mock.patch.object(gate.PLEIOTROPY, "validate_policy", return_value={}),
                mock.patch.object(gate, "placo_result_evidence", return_value=[]),
                mock.patch.object(gate, "assert_no_canonical_lava_results"),
            ):
                # Resume after link 1.
                readiness_path.write_text(readiness, encoding="utf-8")
                gate.publish_family(readiness, payload)
                self.assertEqual(gate.verify_gate(), payload)
                # Resume after link 2 is an exact idempotent verification.
                gate.publish_family(readiness, payload)

            readiness_path.unlink()
            with (
                mock.patch.object(gate, "ROOT", root),
                mock.patch.object(gate, "build_expected", return_value=(readiness, payload)),
            ):
                with self.assertRaisesRegex(gate.GateError, "exact publication prefix"):
                    gate.publish_family(readiness, payload)

    def test_hardlink_reclaim_accounting_counts_physical_bytes_once(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / "work/pair"
            outside = root / "results/canonical"
            base.mkdir(parents=True)
            outside.mkdir(parents=True)
            aligned = base / "aligned.tsv.gz"
            aligned.write_bytes(b"abc")
            (outside / "A.full.tsv.gz").hardlink_to(aligned)
            task = base / "task.tsv"
            task.write_bytes(b"12345")
            (base / "materialization.provenance.json").hardlink_to(task)
            with mock.patch.object(COORDINATOR, "ROOT", root):
                candidates = COORDINATOR.collect_cleanup_candidates(base)
                accounting = COORDINATOR.cleanup_reclaim_accounting(candidates)
                COORDINATOR.validate_cleanup_reclaim_accounting(candidates, accounting)
            self.assertEqual(accounting["authorized_logical_bytes"], 13)
            self.assertEqual(accounting["unique_work_inode_bytes"], 8)
            self.assertEqual(accounting["logical_hardlink_alias_bytes"], 5)
            self.assertEqual(accounting["physically_reclaimable_bytes_at_authorization"], 5)
            self.assertEqual(accounting["externally_retained_hardlink_bytes_at_authorization"], 3)

    def test_preflight_forecasts_pending_prior_canonical_reserve(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contexts = []
            for pair, rows in zip(COORDINATOR.PAIR_ORDER, (10, 20, 30), strict=True):
                contexts.append(({
                    "pair_id": pair, "fingerprint": pair.lower().replace("control", "c") * 64,
                    "sources": [{"rows": rows}, {"rows": rows}],
                }, {"reference_and_publication": {"ld_locus_publication_allowed": False}}))

            def paths(context):
                base = root / "work" / context["pair_id"]
                return {
                    "canonical_ledger": base / "canonical.tsv.gz",
                    "canonical_provenance": base / "canonical.json",
                    "staged_ledger": base / "staged.full.tsv.gz",
                    "aligned": base / "aligned.tsv.gz",
                    "provenance": base / "materialization.provenance.json",
                    "task": base / "task.tsv",
                }

            usage = mock.MagicMock()
            usage.free = 20 * 1024**3
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(
                    COORDINATOR, "verify_execution_contract",
                    return_value={"ram_limit_bytes": 8 * 1024**3, "ram_safety_reserve_bytes": 1024**3},
                ),
                mock.patch.object(COORDINATOR.BRIDGE, "configure_context", side_effect=contexts),
                mock.patch.object(COORDINATOR.BRIDGE, "execution_lock", return_value=nullcontext()),
                mock.patch.object(COORDINATOR, "worker_lease", return_value=nullcontext()),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", side_effect=paths),
                mock.patch.object(COORDINATOR.shutil, "disk_usage", return_value=usage),
            ):
                rows = COORDINATOR.preflight_rows()
            self.assertEqual(rows[0]["pending_prior_canonical_planning_reserve_bytes"], 0)
            self.assertEqual(rows[1]["pending_prior_canonical_planning_reserve_bytes"], 600 * 20)
            self.assertEqual(rows[2]["pending_prior_canonical_planning_reserve_bytes"], 600 * 60)

    def test_reproducibility_bundle_requires_ram_and_all_small_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / "work/pair"
            base.mkdir(parents=True)
            paths = {"base": base}
            context = {"pair_id": "A", "fingerprint": "f" * 64,
                       "policy": {"analysis_id": "track-b-v1.0-pleiotropy"}}
            validated = {
                "counts": {}, "canonical_ledger": {"sha256": "a" * 64},
                "canonical_provenance": {},
                "full_scan_ram": {
                    "phase": "full_scan", "row": {"analysis": COORDINATOR.FULL_SCAN_ANALYSIS},
                    "exit_status": 0, "output_hash": "a" * 64,
                },
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
                mock.patch.object(COORDINATOR, "rebuild_ram_aggregate"),
            ):
                with self.assertRaisesRegex(COORDINATOR.CoordinatorError, "required reproducibility"):
                    COORDINATOR.build_cleanup_bundle(context, validated, {})

    def test_validated_cleanup_is_restartable_and_retains_canonical_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = "f" * 64
            pair = "A"
            base = root / "work/track_b_pleiotropy/placo" / fingerprint / pair
            checkpoint = base / "checkpoints"
            checkpoint.mkdir(parents=True)
            benchmark_fields = COORDINATOR.ENGINE.BENCHMARK_RAW_FIELDS + [
                "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
                "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
                "peak_rss_fraction_of_physical_memory",
            ]
            benchmark_row = {field: "fixture" for field in benchmark_fields}
            benchmark_row.update({
                "global_nuisance_checkpoint_reused": "FALSE", "workers": "1",
                "peak_process_tree_rss_bytes": str(1024**3), "wrapper_wall_seconds": "1",
                "input_sha256": "b" * 64,
            })
            contents = {
                "aligned.tsv.gz": b"regenerable aligned family\n",
                "materialization.provenance.json": b"{\"materialized\":true}\n",
                "task.tsv": b"task\nrow\n",
                "benchmark.raw.tsv": b"raw\nrow\n",
                "benchmark.tsv": COORDINATOR.tsv_bytes(benchmark_fields, [benchmark_row]),
                "staged.full.tsv.gz": b"complete full ledger\n",
                "run.summary.tsv": b"summary\nrow\n",
                "staged.provenance.json": b'{"complete_family_counts":{"rows":1}}\n',
                "checkpoints/nuisance.rds": b"global nuisance all variants\n",
                "checkpoints/nuisance.sha256.tsv": b"hash\nvalue\n",
                "checkpoints/shard_000001.rds": b"complete shard\n",
                "checkpoints/shard_000001.rds.sha256.tsv": b"hash\nvalue\n",
            }
            for relative, content in contents.items():
                path = base / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
            canonical_root = root / "results/track_b/pleiotropy/results/placo"
            canonical_root.mkdir(parents=True)
            canonical_ledger = canonical_root / "A.full.tsv.gz"
            canonical_provenance = canonical_root / "A.provenance.json"
            # Production publishes these as the exact staged inodes.
            canonical_ledger.hardlink_to(base / "staged.full.tsv.gz")
            canonical_provenance.hardlink_to(base / "staged.provenance.json")
            contract_path = root / COORDINATOR.EXECUTION_CONTRACT
            gate_path = root / COORDINATOR.GATE.TERMINAL_GATE_LOCK
            contract_path.parent.mkdir(parents=True, exist_ok=True)
            contract_path.write_text("{\"contract\":true}\n", encoding="utf-8")
            gate_path.write_text("{\"gate\":true}\n", encoding="utf-8")
            context = {
                "pair_id": pair, "fingerprint": fingerprint,
                "policy": {"analysis_id": "track-b-v1.0-pleiotropy"},
            }
            paths = {
                "base": base,
                "task": base / "task.tsv",
                "benchmark": base / "benchmark.tsv",
                "checkpoint_dir": checkpoint,
                "canonical_ledger": canonical_ledger,
                "canonical_provenance": canonical_provenance,
            }
            validated = {
                "counts": {"rows": 1},
                "canonical_ledger": None,
                "canonical_provenance": None,
            }
            contract = {
                "schema_version": COORDINATOR.CONTRACT_SCHEMA,
                "ram_limit_bytes": 8 * 1024**3,
                "ram_safety_reserve_bytes": 1024**3,
            }
            with (
                mock.patch.object(COORDINATOR, "ROOT", root),
                mock.patch.object(COORDINATOR.ENGINE, "materialized_paths", return_value=paths),
                mock.patch.object(
                    COORDINATOR.ENGINE, "one_tsv_row",
                    return_value={"aligned_input_rows": "1"},
                ),
                mock.patch.object(COORDINATOR, "rebuild_ram_aggregate"),
            ):
                benchmark_hash = COORDINATOR.benchmark_bundle_hash(paths, benchmark_row)
                COORDINATOR.record_ram(
                    context, "matched_worker_benchmark",
                    analysis=COORDINATOR.MATCHED_BENCHMARK_ANALYSIS,
                    peak_bytes=1024**3, runtime_seconds=1.0, exit_status=0,
                    output_hash=benchmark_hash,
                    rss_measurement_method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
                )
                ledger_hash = COORDINATOR.stable_identity(canonical_ledger)["sha256"]
                COORDINATOR.record_ram(
                    context, "full_scan", analysis=COORDINATOR.FULL_SCAN_ANALYSIS,
                    peak_bytes=1024**3, runtime_seconds=1.0, exit_status=0,
                    output_hash=str(ledger_hash),
                    rss_measurement_method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
                )
                validated["canonical_ledger"] = COORDINATOR.stable_identity(canonical_ledger)
                validated["canonical_provenance"] = COORDINATOR.stable_identity(canonical_provenance)
                validated["full_scan_ram"] = COORDINATOR.load_ram_measurement(
                    context, "full_scan",
                )
                bundle = COORDINATOR.build_cleanup_bundle(context, validated, contract)
                plan = COORDINATOR.cleanup_bundle_path(context)
                COORDINATOR.publish_bytes_no_replace(
                    plan, json.dumps(bundle, indent=2, sort_keys=True).encode() + b"\n",
                )
                unknown = base / "unknown-after-plan.bin"
                unknown.write_bytes(b"must never be removed")
                with self.assertRaisesRegex(
                    COORDINATOR.CoordinatorError, "unlisted work artifact appeared",
                ):
                    COORDINATOR.complete_cleanup(context, contract)
                self.assertTrue(unknown.is_file())
                self.assertTrue((base / "aligned.tsv.gz").is_file())
                unknown.unlink()
                unknown_directory = base / "unknown-empty-directory"
                unknown_directory.mkdir()
                with self.assertRaisesRegex(
                    COORDINATOR.CoordinatorError, "unlisted work directory appeared",
                ):
                    COORDINATOR.complete_cleanup(context, contract)
                self.assertTrue(unknown_directory.is_dir())
                self.assertTrue((base / "aligned.tsv.gz").is_file())
                unknown_directory.rmdir()
                completion = COORDINATOR.complete_cleanup(context, contract)
                self.assertEqual(completion["removed_artifact_count"], len(contents))
                self.assertTrue(canonical_ledger.is_file())
                self.assertTrue(canonical_provenance.is_file())
                self.assertFalse(base.exists())
                COORDINATOR.validate_clean_completed_pair(context, contract)

    def test_source_contains_no_ld_materializer_execution(self) -> None:
        source = (ROOT / "scripts/143_run_track_b_placo_sequential_v2.py").read_text(encoding="utf-8")
        self.assertNotIn("materialize_reference(", source)
        self.assertNotIn("127_prepare_track_b_pleiotropy_ld.py\"),", source)
        self.assertIn('"extraction_during_full_p_scans": False', source)
        self.assertIn('"canonical_locus_claim_allowed": False', source)


if __name__ == "__main__":
    unittest.main()
