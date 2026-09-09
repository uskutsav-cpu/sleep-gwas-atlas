from contextlib import nullcontext
import hashlib
import importlib.util
import json
import math
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load(name: str, relative: str):
    specification = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


GATE = load("track_b_placo_terminal_gate_v2_test", "scripts/139_build_track_b_placo_terminal_gate_v2.py")
BRIDGE = load("track_b_placo_materializer_v2_test", "scripts/140_materialize_track_b_placo_pair_v2.py")


def terminal_payload() -> dict:
    return {
        "schema_version": "track-b-lava-terminal-qc.1",
        "state": "TERMINAL_FAILED_QC",
        "analysis_id": "track-b-v1.0-local",
        "full_family_complete": True,
        "scientific_validation_passed": False,
        "canonical_publication_allowed": False,
        "terminal": True,
        "qc": {
            "counts": {
                "locus": {"failed": 63, "family": 2495},
                "univariate": {"failed_or_untested": 4510, "family": 19960},
            },
            "thresholds": {
                "maximum_locus_failure_fraction": 0.01,
                "maximum_univariate_untested_fraction": 0.05,
            },
            "reasons": [
                {"code": "LOCUS_FAILURE_FRACTION_EXCEEDED"},
                {"code": "UNIVARIATE_UNTESTED_FRACTION_EXCEEDED"},
            ],
        },
    }


class TrackBPLACOTerminalGateV2Tests(unittest.TestCase):
    @staticmethod
    def archive_state(state: str) -> dict:
        return {
            "state": state,
            "reference_provenance_identity": {
                "path": "ref/lava/ukb_v1.1/reference.provenance.json",
                "bytes": 1,
                "sha256": "a" * 64,
            },
            "archive_family_sha256": "b" * 64,
            "extracted_payload_count": 44,
            "extracted_family_sha256": "c" * 64,
            "eviction_receipt": None,
        }

    def test_reference_archive_guard_accepts_exactly_two_states(self) -> None:
        for state in (GATE.ARCHIVES_PRESENT, GATE.ARCHIVES_EVICTED):
            with mock.patch.object(
                GATE.ARCHIVE_STATE_GUARD, "verify_reference_state",
                return_value=self.archive_state(state),
            ):
                self.assertEqual(GATE.verified_archive_state()["state"], state)
        bad = self.archive_state("PARTIAL_ARCHIVE_EVICTION_RESTART_REQUIRED")
        with (
            mock.patch.object(
                GATE.ARCHIVE_STATE_GUARD, "verify_reference_state", return_value=bad,
            ),
            self.assertRaisesRegex(GATE.GateError, "non-permitted state"),
        ):
            GATE.verified_archive_state()

    def test_evicted_lineage_skips_only_legacy_live_archive_fingerprint(self) -> None:
        source = GATE.LAVA_CONTRACT.PINNED_SOURCE_FINGERPRINT
        continuation = "d" * 64
        state = self.archive_state(GATE.ARCHIVES_EVICTED)
        with (
            mock.patch.object(
                GATE.LAVA_CONTRACT, "validate_source_lock",
            ) as validate_source_lock,
            mock.patch.object(
                GATE.LAVA_CONTRACT, "continuation_fingerprint", return_value=continuation,
            ) as continuation_fingerprint,
            mock.patch.object(
                GATE.LAVA_CONTRACT, "validate_lineage_record",
            ) as validate_lineage_record,
            mock.patch.object(GATE.LAVA_CONTRACT, "validate_lineage") as validate_lineage,
        ):
            self.assertEqual(
                GATE.validated_lava_continuation(state), (source, continuation),
            )
        validate_source_lock.assert_called_once_with(
            source, validate_bundles=True, revalidate_semantics=True,
            require_live_fingerprint=False,
        )
        continuation_fingerprint.assert_called_once_with(
            source, validate_source_bundles=False,
        )
        validate_lineage_record.assert_called_once_with(source, continuation)
        validate_lineage.assert_not_called()

    def test_present_lineage_retains_deep_legacy_validation(self) -> None:
        source = GATE.LAVA_CONTRACT.PINNED_SOURCE_FINGERPRINT
        continuation = "e" * 64
        state = self.archive_state(GATE.ARCHIVES_PRESENT)
        with (
            mock.patch.object(
                GATE.LAVA_CONTRACT, "continuation_fingerprint", return_value=continuation,
            ) as continuation_fingerprint,
            mock.patch.object(GATE.LAVA_CONTRACT, "validate_lineage") as validate_lineage,
            mock.patch.object(GATE.LAVA_CONTRACT, "validate_source_lock") as validate_source_lock,
            mock.patch.object(GATE.LAVA_CONTRACT, "validate_lineage_record") as validate_lineage_record,
        ):
            self.assertEqual(
                GATE.validated_lava_continuation(state), (source, continuation),
            )
        continuation_fingerprint.assert_called_once_with(
            source, validate_source_bundles=True,
        )
        validate_lineage.assert_called_once_with(source, continuation, deep=True)
        validate_source_lock.assert_not_called()
        validate_lineage_record.assert_not_called()

    def test_archive_gate_binding_is_state_invariant(self) -> None:
        present = self.archive_state(GATE.ARCHIVES_PRESENT)
        evicted = self.archive_state(GATE.ARCHIVES_EVICTED)
        evicted["eviction_receipt"] = {"path": "receipt", "sha256": "f" * 64}
        with mock.patch.object(
            GATE, "identity", return_value={"path": "scripts/146", "bytes": 1, "sha256": "9" * 64},
        ):
            self.assertEqual(
                GATE.archive_state_invariant(present),
                GATE.archive_state_invariant(evicted),
            )

    def test_zero_and_truncated_gate_temporaries_are_audited_then_removed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gate_root = root / GATE.GATE_ROOT
            gate_root.mkdir(parents=True)
            zero = gate_root / f".{GATE.READINESS_GATE.name}.111.tmp"
            truncated = gate_root / f".{GATE.TERMINAL_GATE_LOCK.name}.222.tmp"
            zero.write_bytes(b"")
            truncated.write_bytes(b"{truncated")
            with mock.patch.object(GATE, "ROOT", root):
                GATE.reconcile_gate_temporaries([b"readiness\n", b"{}\n"])
            self.assertFalse(zero.exists())
            self.assertFalse(truncated.exists())
            receipts = list((gate_root / "publication_reconciliation").glob("*.json"))
            self.assertEqual(len(receipts), 1)
            payload = json.loads(receipts[0].read_text(encoding="utf-8"))
            self.assertEqual({item["bytes"] for item in payload["artifacts"]}, {0, 10})

    def test_terminal_failed_qc_is_terminal_but_never_passed_evidence(self) -> None:
        GATE.validate_terminal_payload(terminal_payload())
        for key, value in (
            ("scientific_validation_passed", True),
            ("canonical_publication_allowed", True),
            ("terminal", False),
        ):
            payload = terminal_payload()
            payload[key] = value
            with self.assertRaisesRegex(GATE.GateError, "complete TERMINAL_FAILED_QC"):
                GATE.validate_terminal_payload(payload)

    def test_readiness_is_full_p_only_and_keeps_locus_publication_blocked(self) -> None:
        text = GATE.readiness_text()
        self.assertIn("PASS_TERMINAL_FAILED_QC_ATTESTED", text)
        self.assertIn("FAILED_QC_NOT_CONSUMED", text)
        self.assertIn("READY_CONDITIONAL_FULL_P_ONLY", text)
        self.assertIn("BLOCKED_FULL_EUR_HG19_LD_AND_FAMILY_COLLATOR", text)
        self.assertNotIn("LAVA_RESULTS_PASSED", text)

    def test_bridge_accepts_only_the_exact_inherited_v1_snapshot(self) -> None:
        row = {
            "dense_input_gate": "PASS_ALL_THREE_PAIRS",
            "replication_gate": "PASS_TERMINAL_PRIMARY_REPLICATION_FAMILY",
            "local_analysis_gate": "BLOCKED_LOCAL_ANALYSIS_NOT_TERMINAL",
            "software_gate": "READY_PINNED_PLACO_SOURCE",
            "claim_status": "NO_SCIENTIFIC_RESULT",
        }
        BRIDGE.inherited_v1_upstream(row)
        changed = dict(row, local_analysis_gate="PASS_TERMINAL")
        with self.assertRaisesRegex(SystemExit, "UNEXPECTED_V1_LOCAL_GATE_STATE"):
            BRIDGE.inherited_v1_upstream(changed)

    def test_bridge_direct_mutation_cli_is_retired_before_context_access(self) -> None:
        for action in ("--materialize", "--benchmark", "--run"):
            with (
                self.subTest(action=action),
                mock.patch("sys.argv", ["140", "A", action, "--execute"]),
                mock.patch.object(BRIDGE, "configure_context") as configure,
                self.assertRaisesRegex(BRIDGE.BridgeError, "direct PLACO V2 mutation is retired"),
            ):
                BRIDGE.main()
            configure.assert_not_called()

    def test_v2_fingerprint_binds_terminal_gate_bridge_and_failed_qc_semantics(self) -> None:
        fake_gate = {
            "pair_scope": ["A", "B", "CONTROL"],
            "lava_scientific_evidence": {
                "status": "FAILED_QC_NOT_CONSUMED",
                "scientific_validation_passed": False,
                "passed_lava_results_used_by_placo": False,
                "consumption": "ORDERING_TERMINALITY_ATTESTATION_ONLY",
            },
            "resource_envelope": {
                "maximum_workers": 4,
                "ram_safety_reserve_bytes": 1024**3,
            },
        }
        fake_base = {
            "root": ROOT,
            "pair_id": "A",
            "trait1": "snoring",
            "trait2": "parental_lifespan",
            "family_role": "PRIMARY_DISCOVERY",
            "policy": {
                "analysis_id": "track-b-v1.0-pleiotropy",
                "placo_plus": {"source_sha256": "a" * 64},
            },
            "sources": [{"trait_id": "snoring", "absolute_path": ROOT / "unused", "sha256": "b" * 64}],
            "placo_source": ROOT / "unused",
            "fingerprint": "c" * 64,
            "fingerprint_payload": {},
        }
        original = BRIDGE.ENGINE.INPUT_LOCK
        try:
            with (
                mock.patch.object(BRIDGE.GATE, "verify_gate", return_value=fake_gate),
                mock.patch.object(BRIDGE.ENGINE, "production_context", return_value=fake_base),
                mock.patch.object(BRIDGE.ENGINE, "sha256", side_effect=lambda path: "d" * 64),
            ):
                context, observed_gate = BRIDGE.configure_context("A")
        finally:
            BRIDGE.ENGINE.INPUT_LOCK = original
        self.assertIs(observed_gate, fake_gate)
        self.assertRegex(context["fingerprint"], r"^[0-9a-f]{64}$")
        payload = context["fingerprint_payload"]
        self.assertEqual(payload["lava_scientific_evidence"], "FAILED_QC_NOT_CONSUMED")
        self.assertEqual(payload["v2_terminal_gate_lock_sha256"], "d" * 64)
        self.assertEqual(payload["v2_materializer_bridge_sha256"], "d" * 64)
        self.assertEqual(payload["v1_runner_sha256"], "d" * 64)
        self.assertEqual(payload["v2_worker_bound_runner_sha256"], "d" * 64)
        self.assertFalse(payload["ld_locus_publication_allowed"])

    def benchmark_fixture(self, *, workers: int, method: str) -> tuple[dict, dict, dict, dict]:
        task = {
            "analysis_id": "track-b-v1.0-pleiotropy", "pair_id": "A",
            "run_fingerprint": "e" * 64, "aligned_input_sha256": "f" * 64,
            "aligned_input_rows": "4320000", "marginal_p_threshold": "0.0001",
            "policy_sha256": "1" * 64, "contract_lock_sha256": "2" * 64,
            "placo_source_sha256": "3" * 64, "minimum_free_bytes": str(8 * 1024**3),
        }
        row = {field: "x" for field in BRIDGE.ENGINE.BENCHMARK_RAW_FIELDS + [
            "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
            "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
            "peak_rss_fraction_of_physical_memory",
        ]}
        row.update({
            "analysis_id": task["analysis_id"], "pair_id": task["pair_id"],
            "run_fingerprint": task["run_fingerprint"], "input_sha256": task["aligned_input_sha256"],
            "task_sha256": "4" * 64, "policy_sha256": task["policy_sha256"],
            "contract_lock_sha256": task["contract_lock_sha256"], "input_gate_lock_sha256": "5" * 64,
            "materializer_sha256": "6" * 64, "runner_sha256": "7" * 64,
            "placo_source_sha256": task["placo_source_sha256"],
            "input_rows": task["aligned_input_rows"],
            "nuisance_checkpoint_sha256": "8" * 64,
            "nuisance_null_rows_variance": "4200000",
            "nuisance_null_rows_correlation": "4100000",
            "VarZ1": "1.05", "VarZ2": "0.95", "CorZ": "0.1",
            "global_nuisance_input_rows": task["aligned_input_rows"],
            "global_nuisance_elapsed_seconds": "0.5",
            "global_nuisance_checkpoint_reused": "FALSE",
            "benchmark_variants": "2000", "numerical_failures": "0",
            "single_variant_testing_elapsed_seconds": "3.3333333333333335",
            "measured_runner_elapsed_seconds": "4.5",
            "scientific_equivalence": (
                "GLOBAL_OFFICIAL_NUISANCE_ESTIMATED_ON_ALL_VALID_VARIANTS_"
                "THEN_IDENTICAL_SINGLE_VARIANT_PLACO_PLUS_CALLS_SAMPLED"
            ),
            "workers": str(workers), "peak_process_tree_rss_bytes": str(2 * 1024**3),
            "rss_measurement_method": method, "projected_full_family_seconds": "7200",
            "projected_full_family_hours": "2", "variants_per_second": "600",
            "wrapper_wall_seconds": "5", "rss_sample_interval_seconds": "0.05",
            "host_physical_memory_bytes": str(8 * 1024**3),
            "host_free_bytes_at_benchmark": str(20 * 1024**3),
            "peak_rss_fraction_of_physical_memory": "0.25",
            "r_version": "R test", "data_table_version": "1.0", "runtime_platform": "test",
        })
        nuisance_lock = {
            "run_fingerprint": task["run_fingerprint"],
            "input_sha256": task["aligned_input_sha256"], "task_sha256": "4" * 64,
            "placo_source_sha256": task["placo_source_sha256"], "runner_sha256": "7" * 64,
            "workers": str(workers),
            "marginal_p_threshold": task["marginal_p_threshold"],
            "variance_null_rows": row["nuisance_null_rows_variance"],
            "correlation_null_rows": row["nuisance_null_rows_correlation"],
            "nuisance_rds_sha256": "8" * 64,
        }
        gate = {"resource_envelope": {"maximum_workers": 4, "ram_safety_reserve_bytes": 1024**3}}
        return task, row, nuisance_lock, gate

    def test_multiworker_scan_requires_aggregate_rss_but_one_worker_fallback_is_legal(self) -> None:
        context = {"pair_id": "A"}
        checkpoint_dir = ROOT / "checkpoints"
        paths = {
            "benchmark": ROOT / "benchmark.tsv", "task": ROOT / "task.tsv",
            "checkpoint_dir": checkpoint_dir,
        }
        for workers, method, accepted in (
            (1, "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE", True),
            (4, "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE", False),
            (4, "AGGREGATE_PROCESS_TREE_PS_SAMPLED", True),
        ):
            with self.subTest(workers=workers, method=method):
                task, row, nuisance_lock, gate = self.benchmark_fixture(workers=workers, method=method)
                fields = list(row)

                def digest(path):
                    if path == paths["task"]:
                        return "4" * 64
                    if path == ROOT / GATE.TERMINAL_GATE_LOCK:
                        return "5" * 64
                    if path == BRIDGE.ENGINE_PATH:
                        return "6" * 64
                    if path == ROOT / BRIDGE.ENGINE.RUNNER:
                        return "7" * 64
                    if path == checkpoint_dir / "nuisance.rds":
                        return "8" * 64
                    return "0" * 64

                def read_tsv(path):
                    if path == paths["benchmark"]:
                        return fields, [row]
                    if path == checkpoint_dir / "nuisance.sha256.tsv":
                        return BRIDGE.NUISANCE_LOCK_FIELDS, [nuisance_lock]
                    raise AssertionError(path)

                with (
                    mock.patch.object(BRIDGE.ENGINE, "materialized_paths", return_value=paths),
                    mock.patch.object(BRIDGE.ENGINE, "verify_materialized"),
                    mock.patch.object(BRIDGE.ENGINE, "read_tsv", side_effect=read_tsv),
                    mock.patch.object(BRIDGE.ENGINE, "one_tsv_row", return_value=task),
                    mock.patch.object(BRIDGE.ENGINE, "sha256", side_effect=digest),
                    mock.patch.object(BRIDGE.ENGINE, "physical_memory_bytes", return_value=8 * 1024**3),
                    mock.patch.object(BRIDGE.shutil, "disk_usage", return_value=SimpleNamespace(free=20 * 1024**3)),
                ):
                    if accepted:
                        result = BRIDGE.validate_benchmark(context, gate, required_workers=workers)
                        self.assertEqual(int(result["workers"]), workers)
                        self.assertTrue(math.isfinite(float(result["projected_full_family_seconds"])))
                    else:
                        with self.assertRaisesRegex(BRIDGE.BridgeError, "aggregate process-tree RSS"):
                            BRIDGE.validate_benchmark(context, gate, required_workers=workers)

    def test_benchmark_rejects_partial_global_nuisance_scope(self) -> None:
        context = {"pair_id": "A"}
        checkpoint_dir = ROOT / "checkpoints"
        paths = {
            "benchmark": ROOT / "benchmark.tsv", "task": ROOT / "task.tsv",
            "checkpoint_dir": checkpoint_dir,
        }
        task, row, nuisance_lock, gate = self.benchmark_fixture(
            workers=1, method="AGGREGATE_PROCESS_TREE_PS_SAMPLED",
        )
        row["global_nuisance_input_rows"] = "2000"

        def digest(path):
            return {
                paths["task"]: "4" * 64,
                ROOT / GATE.TERMINAL_GATE_LOCK: "5" * 64,
                BRIDGE.ENGINE_PATH: "6" * 64,
                ROOT / BRIDGE.ENGINE.RUNNER: "7" * 64,
                checkpoint_dir / "nuisance.rds": "8" * 64,
            }.get(path, "0" * 64)

        def read_tsv(path):
            if path == paths["benchmark"]:
                return list(row), [row]
            return BRIDGE.NUISANCE_LOCK_FIELDS, [nuisance_lock]

        with (
            mock.patch.object(BRIDGE.ENGINE, "materialized_paths", return_value=paths),
            mock.patch.object(BRIDGE.ENGINE, "verify_materialized"),
            mock.patch.object(BRIDGE.ENGINE, "read_tsv", side_effect=read_tsv),
            mock.patch.object(BRIDGE.ENGINE, "one_tsv_row", return_value=task),
            mock.patch.object(BRIDGE.ENGINE, "sha256", side_effect=digest),
        ):
            with self.assertRaisesRegex(BRIDGE.BridgeError, "global nuisance fit"):
                BRIDGE.validate_benchmark(context, gate, required_workers=1)

    def test_full_scan_is_forced_through_stage_only_before_v2_publication(self) -> None:
        context = {"pair_id": "A"}
        gate = {"resource_envelope": {"maximum_workers": 4}}
        paths = {"task": ROOT / "task.tsv", "run_summary": ROOT / "summary.tsv"}
        benchmark = {
            "nuisance_checkpoint_sha256": "8" * 64,
            "nuisance_null_rows_variance": "25",
            "nuisance_null_rows_correlation": "24",
            "VarZ1": "1.1", "VarZ2": "0.9", "CorZ": "0.1",
            "r_version": "R test", "data_table_version": "1.0", "runtime_platform": "test",
        }
        summary = {**benchmark, "workers": "1"}
        published = {
            "reference_sha256": "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING",
            "output_rows": 5_000_000, "terminal_result_state": "COMPLETE_WITH_HITS",
            "execution_resource_envelope": {"observed_workers": 1},
        }

        def one_tsv_row(path):
            return {} if path == paths["task"] else summary

        def publish_run(*_args):
            # The real V1 publisher invokes this only after its full semantic
            # ledger scan.  Exercise the injected V2 boundary callbacks.
            BRIDGE.ENGINE.exclusive_family([ROOT / "source"], [ROOT / "destination"])
            return published

        def publication_boundary_probe(_sources, _destinations, **callbacks):
            callbacks["before_publication"]()
            callbacks["after_publication"]()

        with (
            mock.patch.object(BRIDGE, "execution_lock", return_value=nullcontext()),
            mock.patch.object(BRIDGE, "validate_v2_paths"),
            mock.patch.object(BRIDGE, "validate_benchmark", return_value=benchmark) as validate_benchmark,
            mock.patch.object(BRIDGE.ENGINE, "resource_preflight"),
            mock.patch.object(BRIDGE.ENGINE, "materialized_paths", return_value=paths),
            mock.patch.object(BRIDGE.ENGINE, "verify_materialized"),
            mock.patch.object(BRIDGE.ENGINE, "relative", return_value="work/task.tsv"),
            mock.patch.object(BRIDGE.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout="", stderr="")) as run,
            mock.patch.object(BRIDGE.ENGINE, "one_tsv_row", side_effect=one_tsv_row),
            mock.patch.object(BRIDGE.ENGINE, "validate_task_against_context"),
            mock.patch.object(BRIDGE.ENGINE, "publish_run", side_effect=publish_run),
            mock.patch.object(BRIDGE, "exclusive_family_v2", side_effect=publication_boundary_probe),
            mock.patch.object(BRIDGE.GATE, "verify_gate", return_value=gate) as verify_gate,
            mock.patch.object(BRIDGE.GATE, "assert_no_canonical_lava_results") as no_lava,
        ):
            result = BRIDGE.run_full_scan(context, gate, 1)
        command = run.call_args.args[0]
        self.assertIn("--stage-only", command)
        self.assertIn("--execute", command)
        self.assertIn(str(ROOT / BRIDGE.V2_RUNNER), command)
        self.assertEqual(validate_benchmark.call_count, 5)
        self.assertEqual(verify_gate.call_count, 4)
        self.assertEqual(no_lava.call_count, 4)
        self.assertEqual(result, published)

    def test_full_scan_rejects_staged_worker_resume_that_differs_from_admission(self) -> None:
        paths = {"run_summary": ROOT / "summary.tsv"}
        benchmark = {
            "nuisance_checkpoint_sha256": "8" * 64,
            "nuisance_null_rows_variance": "25", "nuisance_null_rows_correlation": "24",
            "VarZ1": "1.1", "VarZ2": "0.9", "CorZ": "0.1",
            "r_version": "R test", "data_table_version": "1.0", "runtime_platform": "test",
        }
        summary = {**benchmark, "workers": "4"}
        with mock.patch.object(BRIDGE.ENGINE, "one_tsv_row", return_value=summary):
            with self.assertRaisesRegex(BRIDGE.BridgeError, "admitted worker count"):
                BRIDGE.validate_stage_against_benchmark(paths, benchmark, workers=1)

    def test_gate_publication_rejects_symlinked_namespace(self) -> None:
        with tempfile.TemporaryDirectory() as directory, tempfile.TemporaryDirectory() as outside:
            root = Path(directory)
            (root / "results").symlink_to(Path(outside), target_is_directory=True)
            with mock.patch.object(GATE, "ROOT", root):
                with self.assertRaisesRegex(GATE.GateError, "symbolic link"):
                    GATE.publish_family("ready\n", {"probe": True})
            self.assertEqual(list(Path(outside).iterdir()), [])

    def test_inode_safe_rollback_never_unlinks_a_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "published.tsv"
            path.write_text("ours\n", encoding="utf-8")
            ours = os.lstat(path)
            replacement = root / "replacement.tsv"
            replacement.write_text("replacement\n", encoding="utf-8")
            os.replace(replacement, path)
            with mock.patch.object(GATE, "ROOT", root):
                GATE.unlink_if_identity(path, (ours.st_dev, ours.st_ino))
            self.assertEqual(path.read_text(encoding="utf-8"), "replacement\n")

    def test_canonical_publication_boundary_preserves_replacement_on_rollback(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sources = [root / "stage.tsv", root / "provenance.json"]
            destinations = [root / "canonical/result.tsv", root / "canonical/result.json"]
            for index, path in enumerate(sources):
                path.write_text(f"source-{index}\n", encoding="utf-8")

            def replace_then_fail():
                replacement = root / "replacement.tsv"
                replacement.write_text("replacement\n", encoding="utf-8")
                os.replace(replacement, destinations[0])
                raise BRIDGE.BridgeError("post-publication revalidation failed")

            committed = []
            with mock.patch.object(BRIDGE.GATE, "ROOT", root):
                with self.assertRaisesRegex(BRIDGE.BridgeError, "post-publication"):
                    BRIDGE.exclusive_family_v2(
                        sources,
                        destinations,
                        before_publication=lambda: None,
                        after_publication=replace_then_fail,
                        committed=committed,
                    )
            self.assertEqual(destinations[0].read_text(encoding="utf-8"), "replacement\n")
            self.assertFalse(destinations[1].exists())
            self.assertEqual(committed, [])

    def test_canonical_publication_supports_one_provenance_source_twice(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ledger = root / "stage.tsv"
            provenance = root / "provenance.json"
            ledger.write_text("ledger\n", encoding="utf-8")
            provenance.write_text("{}\n", encoding="utf-8")
            destinations = [
                root / "canonical/result.tsv",
                root / "canonical/result.json",
                root / "stage.provenance.json",
            ]
            committed = []
            with mock.patch.object(BRIDGE.GATE, "ROOT", root):
                BRIDGE.exclusive_family_v2(
                    [ledger, provenance, provenance],
                    destinations,
                    before_publication=lambda: None,
                    after_publication=lambda: None,
                    committed=committed,
                )
                self.assertEqual(len(committed), 3)
                self.assertEqual(os.lstat(destinations[1]).st_ino, os.lstat(destinations[2]).st_ino)
                BRIDGE.rollback_v2_publication(committed)
            self.assertFalse(any(path.exists() for path in destinations))

    def test_pre_result_scan_includes_work_stages_benchmarks_and_checkpoints(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / GATE.PLACO_WORK_ROOT / ("a" * 64) / "A"
            checkpoint = base / "checkpoints/shard_000001.rds"
            checkpoint.parent.mkdir(parents=True)
            checkpoint.write_bytes(b"checkpoint")
            (base / "staged.full.tsv.gz").write_bytes(b"ledger")
            (base / "benchmark.raw.tsv").write_text("benchmark\n", encoding="utf-8")
            (base / "staged.full.tsv.gz.123.tmp").write_bytes(b"crash ledger")
            (base / ".benchmark.tsv.123.tmp").write_text("crash benchmark\n", encoding="utf-8")
            (base / ".staged.provenance.json.123.tmp").write_text("{}\n", encoding="utf-8")
            with (
                mock.patch.object(GATE, "ROOT", root),
                mock.patch.object(GATE.PLEIOTROPY, "result_evidence", return_value=[]),
            ):
                evidence = GATE.placo_result_evidence({})
            self.assertEqual(
                {path.name for path in evidence},
                {
                    "shard_000001.rds", "staged.full.tsv.gz", "benchmark.raw.tsv",
                    "staged.full.tsv.gz.123.tmp", ".benchmark.tsv.123.tmp",
                    ".staged.provenance.json.123.tmp",
                },
            )

    def test_verify_gate_rejects_symlinked_readiness_even_with_matching_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory, tempfile.TemporaryDirectory() as outside:
            root = Path(directory)
            target = Path(outside) / "readiness.tsv"
            target.write_text("ready\n", encoding="utf-8")
            readiness = root / GATE.READINESS_GATE
            readiness.parent.mkdir(parents=True)
            readiness.symlink_to(target)
            with (
                mock.patch.object(GATE, "ROOT", root),
                mock.patch.object(GATE, "build_expected", return_value=("ready\n", {})),
            ):
                with self.assertRaisesRegex(GATE.GateError, "missing or drifted"):
                    GATE.verify_gate()

    def test_terminal_gate_rejects_any_partial_canonical_lava_family(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            partial = root / GATE.LAVA_CANONICAL_RESULTS[0]
            partial.parent.mkdir(parents=True)
            partial.write_text("partial scientific output\n", encoding="utf-8")
            with mock.patch.object(GATE, "ROOT", root):
                with self.assertRaisesRegex(GATE.GateError, "partial or sealed canonical LAVA"):
                    GATE.assert_no_canonical_lava_results()

    def test_repository_lock_is_cross_pair_nonblocking_and_stale_safe(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = root / "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/.execution.lock"
            with (
                mock.patch.object(BRIDGE.GATE, "ROOT", root),
                mock.patch.object(BRIDGE, "EXECUTION_LOCK", lock),
            ):
                with BRIDGE.execution_lock("materialize", "A"):
                    with self.assertRaisesRegex(BRIDGE.BridgeError, "another A/B/CONTROL"):
                        with BRIDGE.execution_lock("benchmark", "B"):
                            self.fail("a second pair acquired the repository lock")
                self.assertTrue(lock.is_file())
                with BRIDGE.execution_lock("full-scan", "CONTROL"):
                    self.assertIn("pair=CONTROL", lock.read_text(encoding="utf-8"))

    def test_repository_lock_rejects_hardlink_without_truncating_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "frozen-v1-artifact"
            target.write_bytes(b"must remain immutable\n")
            lock = root / "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/.execution.lock"
            lock.parent.mkdir(parents=True)
            os.link(target, lock)
            with (
                mock.patch.object(BRIDGE.GATE, "ROOT", root),
                mock.patch.object(BRIDGE, "EXECUTION_LOCK", lock),
            ):
                with self.assertRaisesRegex(BRIDGE.BridgeError, "private stable regular file"):
                    with BRIDGE.execution_lock("materialize", "A"):
                        self.fail("hard-linked lock was accepted")
            self.assertEqual(target.read_bytes(), b"must remain immutable\n")

    def test_repository_lock_rechecks_path_after_flock_before_writing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = root / "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/.execution.lock"
            lock.parent.mkdir(parents=True)
            lock.write_bytes(b"old lock\n")
            replacement = root / "replacement.lock"
            replacement.write_bytes(b"replacement must not be truncated\n")

            def replace_lock(_descriptor, _operation):
                os.replace(replacement, lock)

            with (
                mock.patch.object(BRIDGE.GATE, "ROOT", root),
                mock.patch.object(BRIDGE, "EXECUTION_LOCK", lock),
                mock.patch.object(BRIDGE.fcntl, "flock", side_effect=replace_lock),
            ):
                with self.assertRaisesRegex(BRIDGE.BridgeError, "changed while acquiring flock"):
                    with BRIDGE.execution_lock("benchmark", "B"):
                        self.fail("replaced lock path was accepted")
            self.assertEqual(lock.read_bytes(), b"replacement must not be truncated\n")

    def test_v2_path_validation_rejects_out_of_tree_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as directory, tempfile.TemporaryDirectory() as outside:
            root = Path(directory)
            redirect = root / "work"
            redirect.symlink_to(Path(outside), target_is_directory=True)
            paths = {"base": redirect / "fingerprint/A"}
            with (
                mock.patch.object(BRIDGE.GATE, "ROOT", root),
                mock.patch.object(BRIDGE.ENGINE, "materialized_paths", return_value=paths),
            ):
                with self.assertRaisesRegex(BRIDGE.BridgeError, "symbolic link"):
                    BRIDGE.validate_v2_paths({})

    def test_v1_sources_remain_byte_frozen_and_v2_runner_is_additive(self) -> None:
        expected = {
            "scripts/123_track_b_pleiotropy_contract.py": "7001146c8548cbb2a6d117f3b5bc7b247d22b3b307b2b67d4cb020f7d64830a2",
            "scripts/124_build_track_b_pleiotropy_input_gate.py": "72dbf1f51a54064247bfe1f2ed06b5ec12ace18fc19445514357a093855c3e56",
            "scripts/125_materialize_track_b_placo_pair.py": "0ec584da9823ba3af8e36e663bb7e10cc6da694e7b95162c27fb1b64e95e9250",
            "scripts/126_run_track_b_placo_pair.R": "d418e36acee10cea4d552beb8feb0c86057a0dc2b2be8bfb436eee08715bf12c",
            "scripts/127_prepare_track_b_pleiotropy_ld.py": "f976ffeb59d4de64291e7909b2f12c2fc2e158a5037cee30b63d136e5cd96b9d",
            ".r-env/share/placo/PLACO_v0.2.0.R": "fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124",
        }
        for relative, digest in expected.items():
            with self.subTest(relative=relative):
                self.assertEqual(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(), digest)
        v1 = (ROOT / BRIDGE.V1_RUNNER).read_text(encoding="utf-8")
        v2 = (ROOT / BRIDGE.V2_RUNNER).read_text(encoding="utf-8")
        for scientific_line in (
            "var_z <- var.placo(z_matrix, p_matrix, p.threshold = marginal_threshold)",
            "cor_z <- cor.pearson(z_matrix, p_matrix, p.threshold = marginal_threshold, returnMatrix = FALSE)",
            "placo.plus(c(dat$Z1[[index]], dat$Z2[[index]]), VarZ = var_z, CorZ = cor_z, AbsTol = absolute_tolerance)",
            'bh_q <- p.adjust(p_values, method = "BH", n = nrow(dat))',
        ):
            self.assertIn(scientific_line, v1)
            self.assertIn(scientific_line, v2)
        self.assertIn('"workers", "marginal_p_threshold"', v2)
        self.assertIn("nThread = workers", v2)
        self.assertIn("unlink_if_identity", v2)
        self.assertIn("V2 PLACO+ runner is stage-only", v2)

    def test_terminal_lock_binds_both_runners_and_runtime_launcher(self) -> None:
        policy = {"analysis_id": "track-b-v1.0-pleiotropy"}
        v1 = {
            "policy": policy, "upstream_identities": {}, "dense_identities": {},
            "placo_source": GATE.identity(ROOT / ".r-env/share/placo/PLACO_v0.2.0.R"),
        }
        payload = GATE.lock_payload(v1, {}, GATE.readiness_text())
        execution = payload["execution_code"]
        self.assertEqual(execution["v1_full_p_runner"]["path"], str(GATE.V1_RUNNER))
        self.assertEqual(execution["v2_worker_bound_full_p_runner"]["path"], str(GATE.V2_RUNNER))
        self.assertEqual(execution["r_runtime_launcher"]["path"], str(GATE.R_RUNTIME))


if __name__ == "__main__":
    unittest.main()
