from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/14_build_reporting_checklists.py"
SPEC = importlib.util.spec_from_file_location("build_reporting_checklists", SCRIPT)
CHECKLISTS = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(CHECKLISTS)


class ReportingChecklistTests(unittest.TestCase):
    def test_checklist_item_counts_and_required_subitems(self) -> None:
        prisma = CHECKLISTS.prisma_rows()
        self.assertEqual(len(prisma), 42)
        self.assertIn("10a", {row["item_id"] for row in prisma})
        self.assertIn("13f", {row["item_id"] for row in prisma})
        self.assertIn("16b", {row["item_id"] for row in prisma})
        self.assertIn("20d", {row["item_id"] for row in prisma})
        self.assertIn("24c", {row["item_id"] for row in prisma})
        self.assertEqual(len(CHECKLISTS.abstract_rows()), 12)
        self.assertEqual(len(CHECKLISTS.prismas_rows()), 16)
        self.assertEqual(len(CHECKLISTS.strega_rows()), 22)

    def test_every_item_has_location_status_and_gap(self) -> None:
        rows = (
            CHECKLISTS.prisma_rows() + CHECKLISTS.abstract_rows()
            + CHECKLISTS.prismas_rows() + CHECKLISTS.strega_rows()
        )
        for row in rows:
            with self.subTest(framework=row["framework"], item=row["item_id"]):
                self.assertTrue(row["planned_manuscript_section"])
                self.assertEqual(row["manuscript_page_line"], "NOT_DRAFTED")
                self.assertTrue(row["status"])
                self.assertTrue(row["remaining_action"])
                self.assertTrue(row["guideline_source"].startswith("https://"))

    def test_rebuild_preserves_existing_reviewed_locator(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "prisma.tsv"
            path.write_text(
                "framework\titem_id\tmanuscript_page_line\n"
                "PRISMA-2020\t1\tmanuscript.md:L13-L15 (current title and abstract)\n",
                encoding="utf-8",
            )
            rows = [CHECKLISTS.prisma_rows()[0]]
            result = CHECKLISTS.preserve_existing_locators(rows, path)
            self.assertEqual(result[0]["manuscript_page_line"], "manuscript.md:L13-L15 (current title and abstract)")


if __name__ == "__main__":
    unittest.main()
