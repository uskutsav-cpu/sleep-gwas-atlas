from __future__ import annotations

import importlib.util
import tempfile
import unittest
from datetime import date
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/54_collect_pubmed_window.py"
SPEC = importlib.util.spec_from_file_location("collect_pubmed_window", SCRIPT)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class PubMedWindowTests(unittest.TestCase):
    def test_accepts_inclusive_ordered_date_window(self):
        start, end = module.date_window(date(2026, 9, 23), date(2026, 9, 25))
        self.assertEqual((start.isoformat(), end.isoformat()), ("2026-09-23", "2026-09-25"))

    def test_rejects_reversed_date_window(self):
        with self.assertRaisesRegex(ValueError, "on or before"):
            module.date_window(date(2026, 9, 25), date(2026, 9, 23))

    def test_refuses_to_write_into_active_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / "config.json"
            config.write_text('{"queries":{"q":"x"}}')
            with patch.object(module, "ACTIVE_CACHE", Path(directory) / "cache"):
                with self.assertRaisesRegex(ValueError, "separate cache"):
                    module.collect(config, Path(directory) / "cache" / "increment", date(2026, 9, 23), date(2026, 9, 25))


if __name__ == "__main__":
    unittest.main()
