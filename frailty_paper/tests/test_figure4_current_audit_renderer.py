"""Regression checks for the dependency-light Figure 4 renderer."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "frailty_paper/scripts/52_refresh_figure4_current_audit.py"
AUDIT = ROOT / "frailty_paper/analysis/lava_full_receipt_integrity_2026-09-26_2307_checkpoint.json"
SPEC = importlib.util.spec_from_file_location("figure4_current_audit_renderer", SCRIPT)
assert SPEC and SPEC.loader
RENDERER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RENDERER)


class Figure4RendererTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.audit = json.loads(AUDIT.read_text(encoding="utf-8"))
        cls.svg = RENDERER.build_svg(cls.audit)

    def test_panel_tracks_the_independent_receipt_audit(self) -> None:
        self.assertIn("18,234/29,940 receipts", self.svg)
        self.assertIn("23:07:06–23:07:58 UTC", self.svg)
        self.assertIn("18,237", self.svg)

    def test_panel_preserves_locked_gate_and_failure_counts(self) -> None:
        self.assertIn("FAILED · 12/12", self.svg)
        self.assertIn("3,255 negative-variance", self.svg)
        self.assertIn("252 no-reference-SNP", self.svg)
        self.assertIn("1 other", self.svg)
        self.assertIn("frozen 1%", self.svg)
        self.assertIn("allowance; thresholds unchanged.", self.svg)

    def test_panel_disallows_local_inference_without_claiming_absence(self) -> None:
        self.assertIn("NO LOCAL INFERENCE", self.svg)
        self.assertIn("do not show sharing is absent.", self.svg)

    def test_inventory_note_reports_stable_inventory_and_runner_count(self) -> None:
        audit = {
            "receipt_file_list_changed_during_scan": False,
            "runner_verified_receipts_at_start": 27_828,
            "runner_verified_receipts_at_end": 27_828,
        }
        self.assertEqual(
            RENDERER.inventory_note(audit),
            "Receipt inventory and runner count remained stable during scanning.",
        )

    def test_inventory_note_reports_advanced_inventory_or_runner_count(self) -> None:
        audit = {
            "receipt_file_list_changed_during_scan": True,
            "runner_verified_receipts_at_start": 27_828,
            "runner_verified_receipts_at_end": 27_830,
        }
        self.assertEqual(
            RENDERER.inventory_note(audit),
            "Receipt inventory or runner count advanced during scanning.",
        )


if __name__ == "__main__":
    unittest.main()
