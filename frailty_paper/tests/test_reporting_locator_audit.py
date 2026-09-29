from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/44_audit_reporting_locators.py"
SPEC = importlib.util.spec_from_file_location("audit_reporting_locators", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(AUDIT)


class ReportingLocatorAuditTests(unittest.TestCase):
    def make_repo(self, locator: str) -> Path:
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        repo = Path(temp.name)
        (repo / "frailty_paper/paper").mkdir(parents=True)
        (repo / "frailty_paper/paper/manuscript.md").write_text("# Title\n\n## Methods\n", encoding="utf-8")
        table_dir = repo / "frailty_paper/paper/reporting_checklists"
        table_dir.mkdir()
        for name in AUDIT.TABLES:
            value = locator if name == "strega.tsv" else "NOT PRESENT IN DRAFT (section absent)"
            (table_dir / name).write_text(
                "framework\titem_id\tmanuscript_page_line\nTEST\t1\t" + value + "\n",
                encoding="utf-8",
            )
        return repo

    def test_in_bounds_locator_passes_and_sources_are_hashed(self) -> None:
        result = AUDIT.audit(self.make_repo("manuscript.md:L1-L3"))
        self.assertEqual(result["result"], "PASS")
        self.assertEqual(result["checklist_rows"], 4)
        self.assertEqual(result["line_locators_checked"], 1)
        self.assertIn("frailty_paper/paper/manuscript.md", result["inputs_sha256"])

    def test_out_of_bounds_locator_fails(self) -> None:
        result = AUDIT.audit(self.make_repo("manuscript.md:L3-L4"))
        self.assertEqual(result["result"], "FAIL")
        self.assertTrue(any("outside current 3-line source" in error for error in result["errors"]))

    def test_locator_must_have_source_reference_or_explicit_absence(self) -> None:
        result = AUDIT.audit(self.make_repo("Methods section"))
        self.assertEqual(result["result"], "FAIL")
        self.assertTrue(any("no line reference or explicit absence status" in error for error in result["errors"]))


if __name__ == "__main__":
    unittest.main()
