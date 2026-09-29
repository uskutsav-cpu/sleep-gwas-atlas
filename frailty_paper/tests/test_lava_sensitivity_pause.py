from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/51_run_lava_sensitivity_parallel.py"
SPEC = importlib.util.spec_from_file_location("lava_sensitivity_parallel", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class LavaSensitivityPauseTests(unittest.TestCase):
    def test_pause_flag_prevents_worker_launch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "concurrency.json").write_text(json.dumps({
                "desired_workers": 4, "pause_requested": True,
            }))
            old_root = MODULE.OUTPUT_ROOT
            MODULE.OUTPUT_ROOT = root
            try:
                with patch.object(sys, "argv", [str(SCRIPT)]):
                    self.assertEqual(MODULE.main(), 0)
            finally:
                MODULE.OUTPUT_ROOT = old_root
            self.assertFalse((root / "runner.lock").exists())
            self.assertFalse((root / "workers").exists())

    def test_resume_releases_only_pause_holds_and_restores_four_workers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            claim = root / "claims" / "longsleep" / "locus_0001.claim"
            claim.parent.mkdir(parents=True)
            claim.write_text(json.dumps({"pause_hold": True}) + "\n")
            control_path = root / "concurrency.json"
            control = {"desired_workers": 4, "pause_requested": True}
            control_path.write_text(json.dumps(control))

            released = MODULE.release_pause_claims(root, control_path, control)

            self.assertEqual(released, 1)
            self.assertFalse(claim.exists())
            updated = json.loads(control_path.read_text())
            self.assertFalse(updated["pause_requested"])
            self.assertEqual(updated["desired_workers"], 4)

    def test_resume_refuses_non_pause_claims(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            claim = root / "claims" / "longsleep" / "locus_0001.claim"
            claim.parent.mkdir(parents=True)
            claim.write_text(json.dumps({"worker_pid": 123}) + "\n")
            control_path = root / "concurrency.json"
            control = {"desired_workers": 4, "pause_requested": True}

            with self.assertRaisesRegex(RuntimeError, "non-pause claims"):
                MODULE.release_pause_claims(root, control_path, control)

            self.assertTrue(claim.exists())

    def test_stale_recovery_preserves_pause_holds(self):
        with tempfile.TemporaryDirectory() as temporary:
            claim_root = Path(temporary) / "claims"
            claim = claim_root / "chronotype" / "locus_2377.claim"
            claim.parent.mkdir(parents=True)
            claim.write_text(json.dumps({
                "pause_hold": True, "process_group": 123,
            }) + "\n")

            with patch.object(MODULE, "process_group_alive", return_value=False):
                recovered = MODULE.recover_stale_claims(claim_root)

            self.assertEqual(recovered, 0)
            self.assertTrue(claim.exists())


if __name__ == "__main__":
    unittest.main()
