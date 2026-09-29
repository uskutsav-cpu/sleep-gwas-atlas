"""Fail-closed tests for power-screen runtime process visibility."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from unittest import TestCase, mock


SCRIPT = Path(__file__).resolve().parents[1] / "audit_power_screen_runtime.py"
SPEC = importlib.util.spec_from_file_location("audit_power_screen_runtime", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class ProcessVisibilityTests(TestCase):
    def test_denied_process_table_is_unavailable_not_stopped(self) -> None:
        with mock.patch.object(
            AUDIT.subprocess,
            "run",
            side_effect=PermissionError("Operation not permitted: ps"),
        ):
            result = AUDIT.snapshot()

        self.assertEqual(result["status"], "PROCESS_INVENTORY_UNAVAILABLE")
        self.assertIsNone(result["active_workers"])
        self.assertEqual(
            result["process_inventory_error"],
            {
                "type": "PermissionError",
                "message": "Operation not permitted: ps",
            },
        )
        self.assertEqual(len(result["workers"]), 4)
        self.assertTrue(all(worker["pid"] is None for worker in result["workers"]))

    def test_readable_process_table_still_reports_active_workers(self) -> None:
        lines = [
            f"{1000 + worker} 1.0 1000 00:01 "
            f"R --file=run_power_optimized_sleep_univariate_v1.R "
            f"--args worker_{worker}.config.json"
            for worker in range(1, 5)
        ]
        completed = mock.Mock(stdout="\n".join(lines))
        with mock.patch.object(AUDIT.subprocess, "run", return_value=completed):
            result = AUDIT.snapshot()

        self.assertEqual(result["status"], "RUNNING")
        self.assertEqual(result["active_workers"], 4)
        self.assertIsNone(result["process_inventory_error"])

    def test_process_table_command_failure_is_unavailable_not_stopped(self) -> None:
        error = AUDIT.subprocess.CalledProcessError(
            1, ["ps", "-axo", "pid=,pcpu=,rss=,etime=,command="]
        )
        with mock.patch.object(AUDIT.subprocess, "run", side_effect=error):
            result = AUDIT.snapshot()

        self.assertEqual(result["status"], "PROCESS_INVENTORY_UNAVAILABLE")
        self.assertIsNone(result["active_workers"])
        self.assertEqual(result["process_inventory_error"]["type"], "CalledProcessError")
