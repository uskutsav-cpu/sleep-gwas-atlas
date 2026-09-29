from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/audit_europe_pmc_missing_abstracts.py"
SPEC = importlib.util.spec_from_file_location("europe_pmc_missing_abstract_audit", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class EuropePmcMissingAbstractAuditTests(unittest.TestCase):
    def test_heading_only_record_has_no_abstract_body(self) -> None:
        self.assertEqual(MODULE.abstract_body("<h4>Introduction</h4>"), "")

    def test_structured_abstract_keeps_body_and_excludes_headings(self) -> None:
        raw = "<h4>Introduction</h4><p>First &amp; important sentence.</p><h4>Methods</h4><p>Study details.</p>"
        self.assertEqual(MODULE.abstract_body(raw), "First & important sentence. Study details.")

    def test_plain_text_abstract_is_retained(self) -> None:
        self.assertEqual(MODULE.abstract_body("A plain abstract."), "A plain abstract.")

    def test_queue_selects_only_blank_abstract_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            queue = Path(temp) / "queue.tsv"
            queue.write_text("pmid\tabstract\n123\t\n456\tHas text\n", encoding="utf-8")
            self.assertEqual(MODULE.get_queue_pmids(queue), ["123"])

    def test_queue_rejects_duplicate_blank_abstract_pmids(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            queue = Path(temp) / "queue.tsv"
            queue.write_text("pmid\tabstract\n123\t\n123\t\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate PMIDs"):
                MODULE.get_queue_pmids(queue)

    def test_queue_requires_abstract_and_pmid_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            queue = Path(temp) / "queue.tsv"
            queue.write_text("pmid\ttitle\n123\tTitle\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "required pmid and abstract"):
                MODULE.get_queue_pmids(queue)


if __name__ == "__main__":
    unittest.main()
