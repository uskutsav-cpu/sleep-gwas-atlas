from __future__ import annotations

import importlib.util
import multiprocessing
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/51_run_lava_sensitivity_parallel.py"
SPEC = importlib.util.spec_from_file_location("frailty_lava_parallel", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _race_claim(root: str, barrier, results) -> None:
    barrier.wait()
    result = MODULE.claim_job(Path(root), 1, "accel_sleep_duration", 1)
    results.put(bool(result))


class LavaParallelSchedulerTests(unittest.TestCase):
    def test_coordinator_polls_launched_workers_to_reap_exited_children(self) -> None:
        class ExitedWorker:
            def __init__(self) -> None:
                self.polled = False

            def poll(self) -> int:
                self.polled = True
                return 0

        worker = ExitedWorker()
        MODULE.reap_worker_processes({123: worker})
        self.assertTrue(worker.polled)

    def test_queue_is_deterministic_and_covers_each_pair_locus_once(self) -> None:
        first = MODULE.jobs()
        self.assertEqual(first, MODULE.jobs())
        self.assertEqual(len(first), 29940)
        self.assertEqual(len(set(first)), 29940)
        self.assertEqual(first[0], (1, "accel_sleep_duration", 1))
        self.assertEqual(first[1], (2, "chronotype", 1))

    def test_atomic_claim_allows_only_one_concurrent_winner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            context = multiprocessing.get_context("fork")
            barrier = context.Barrier(2)
            results = context.Queue()
            processes = [context.Process(target=_race_claim, args=(tmp, barrier, results)) for _ in range(2)]
            for process in processes:
                process.start()
            for process in processes:
                process.join(5)
                self.assertEqual(process.exitcode, 0)
            self.assertEqual(sorted([results.get(timeout=1), results.get(timeout=1)]), [False, True])
            claim = Path(tmp) / "accel_sleep_duration/locus_0001.claim"
            self.assertTrue(claim.is_file())
            self.assertEqual(MODULE.recover_stale_claims(Path(tmp)), 0)

    def test_incomplete_claim_is_retained_while_workers_may_be_alive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            claim = Path(tmp) / "accel_sleep_duration/locus_0001.claim"
            claim.parent.mkdir()
            claim.write_text('{"process_group":')
            self.assertEqual(MODULE.recover_stale_claims(Path(tmp), {123}), 0)
            self.assertTrue(claim.exists())
            self.assertEqual(MODULE.recover_stale_claims(Path(tmp), set()), 1)
            self.assertFalse(claim.exists())

    def test_recovery_removes_a_claim_only_after_its_process_group_exits(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            child = subprocess.Popen(["/bin/sleep", "30"], start_new_session=True)
            claim_root = Path(tmp)
            claim = claim_root / "chronotype/locus_0001.claim"
            claim.parent.mkdir()
            claim.write_text(json.dumps({"process_group": child.pid}) + "\n")
            self.assertEqual(MODULE.recover_stale_claims(claim_root), 0)
            self.assertTrue(claim.exists())
            child.terminate()
            child.wait(timeout=5)
            self.assertEqual(MODULE.recover_stale_claims(claim_root), 1)
            self.assertFalse(claim.exists())

    def test_swap_monitor_parses_macos_vm_swapusage(self) -> None:
        result = subprocess.CompletedProcess(
            args=["sysctl", "vm.swapusage"], returncode=0,
            stdout="vm.swapusage: total = 11264.00M  used = 10157.94M  free = 1106.06M (encrypted)\n",
            stderr="")
        with patch.object(MODULE.subprocess, "run", return_value=result):
            self.assertEqual(MODULE.swap_usage_mib(), (10157.94, 11264.0))

    def test_six_workers_fall_back_atomically_when_swap_is_high(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "concurrency.json"
            control = {"desired_workers": 6}
            workers, reason = MODULE.maybe_reduce_workers(path, control, 6, 0, (8600.0, 10000.0))
            self.assertEqual(workers, 4)
            self.assertIn("swap use", reason)
            self.assertEqual(json.loads(path.read_text()), control)
            self.assertEqual(control["desired_workers"], 4)

    def test_six_workers_continue_below_swap_fallback_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 6}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 6, 0, (800.0, 10000.0))
            self.assertEqual((workers, reason), (6, None))
            self.assertEqual(control, {"desired_workers": 6})

    def test_six_workers_continue_at_current_swap_baseline_below_fallback_rung(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 6}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 6, 0, (4814.0, 6144.0))
            self.assertEqual((workers, reason), (6, None))

    def test_six_workers_fall_back_to_four_at_absolute_cap_with_expanded_swap_pool(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 6}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 6, 0, (5100.0, 8192.0))
            self.assertEqual(workers, 4)
            self.assertIn("reducing six workers to four", reason)

    def test_six_workers_fall_back_when_sustained_throughput_is_materially_worse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 6}
            workers, reason = MODULE.maybe_reduce_workers(
                Path(tmp) / "concurrency.json", control, 6, 0, None,
                observed_receipts_per_hour=400.0)
            self.assertEqual(workers, 4)
            self.assertIn("reducing six workers to four", reason)
            self.assertEqual(control["resource_fallback_events"][0]["from_workers"], 6)
            self.assertEqual(control["resource_fallback_events"][0]["to_workers"], 4)

    def test_six_workers_do_not_downshift_at_or_above_throughput_floor(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 6}
            floor = MODULE.FOUR_WORKER_BASELINE_RECEIPTS_PER_HOUR * MODULE.THROUGHPUT_FALLBACK_FRACTION
            workers, reason = MODULE.maybe_reduce_workers(
                Path(tmp) / "concurrency.json", control, 6, 0, None,
                observed_receipts_per_hour=floor)
            self.assertEqual((workers, reason), (6, None))
            self.assertEqual(control, {"desired_workers": 6})

    def test_four_worker_fallback_does_not_cascade_at_six_worker_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 4}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 4, 0, (5100.0, 6144.0))
            self.assertEqual((workers, reason), (4, None))

    def test_severe_preexisting_swap_blocks_launch_without_releasing_holds(self) -> None:
        reason = MODULE.launch_resource_block((11698.5, 12288.0))
        self.assertIn("pause holds retained", reason)

    def test_launch_resource_gate_allows_below_severe_swap_limit(self) -> None:
        self.assertIsNone(MODULE.launch_resource_block((10000.0, 12288.0)))

    def test_six_workers_fall_back_after_repeated_worker_failures(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 6}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 6, 2, None)
            self.assertEqual(workers, 4)
            self.assertIn("repeated worker failures", reason)

    def test_four_workers_step_down_only_at_deep_pressure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 4}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 4, 0, (9600.0, 10000.0))
            self.assertEqual(workers, 3)
            self.assertIn("reducing four workers to three", reason)
            self.assertEqual(control["resource_fallback_events"][0]["from_workers"], 4)
            self.assertEqual(control["resource_fallback_events"][0]["to_workers"], 3)

    def test_four_workers_step_down_after_swap_capacity_expands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 4}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 4, 0, (5900.0, 6144.0))
            self.assertEqual(workers, 3)
            self.assertIn("reducing four workers to three", reason)
            self.assertIn("5900/6144 MiB", reason)

    def test_fallback_preserves_a_newer_pause_control(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "concurrency.json"
            path.write_text(json.dumps({"desired_workers": 4, "pause_requested": True,
                                        "pause_id": "new-pause"}))
            stale_control = {"desired_workers": 4, "pause_requested": False}
            workers, _ = MODULE.maybe_reduce_workers(path, stale_control, 4, 0, (9600.0, 10000.0))
            persisted = json.loads(path.read_text())
            self.assertEqual(workers, 3)
            self.assertTrue(persisted["pause_requested"])
            self.assertEqual(persisted["pause_id"], "new-pause")
            self.assertTrue(stale_control["pause_requested"])

    def test_three_workers_step_down_to_two_only_under_deep_pressure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "concurrency.json"
            control = {"desired_workers": 3}
            workers, reason = MODULE.maybe_reduce_workers(path, control, 3, 0, (9600.0, 10000.0))
            self.assertEqual(workers, 2)
            self.assertIn("reducing three workers to two", reason)
            self.assertEqual(json.loads(path.read_text()), control)

    def test_three_workers_step_down_after_swap_capacity_expands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 3}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 3, 0, (5900.0, 6144.0))
            self.assertEqual(workers, 2)
            self.assertIn("reducing three workers to two", reason)
            self.assertIn("5900/6144 MiB", reason)

    def test_two_workers_do_not_step_down_further(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            control = {"desired_workers": 2}
            workers, reason = MODULE.maybe_reduce_workers(Path(tmp) / "concurrency.json",
                                                          control, 2, 10, (9900.0, 10000.0))
            self.assertEqual((workers, reason), (2, None))
            self.assertEqual(control, {"desired_workers": 2})


if __name__ == "__main__":
    unittest.main()
