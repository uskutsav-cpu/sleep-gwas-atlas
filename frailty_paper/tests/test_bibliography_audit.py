from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/42_audit_bibliography.py"
SPEC = importlib.util.spec_from_file_location("bibliography_audit", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class BibliographyAuditTests(unittest.TestCase):
    def test_parses_regular_entries_and_ignores_bibtex_directives(self):
        text = "@string{journal = {Example}}\n@article{alpha, title={A}}\n@book(beta, title={B})"
        self.assertEqual(MODULE.parse_bib_keys(text), ["alpha", "beta"])

    def test_reads_only_pandoc_citation_clusters(self):
        text = "Text [@alpha; @beta, p. 4; -@gamma]. Contact name@example.org."
        self.assertEqual(MODULE.parse_citation_keys(text), ["alpha", "beta", "gamma"])

    def test_audit_reports_missing_uncited_duplicate_and_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manuscript = root / "manuscript.md"
            references = root / "references.bib"
            manuscript.write_text("A claim [@present; @missing].\n", encoding="utf-8")
            references.write_text(
                "@article{present, title={A}}\n@article{present, title={B}}\n"
                "@article{unused, title={C}}\n",
                encoding="utf-8",
            )
            result = MODULE.audit(manuscript, references, SCRIPT)
            self.assertEqual(result["status"], "FAIL")
            self.assertEqual(result["duplicate_keys"], ["present"])
            self.assertEqual(result["missing_keys"], ["missing"])
            self.assertEqual(result["uncited_keys"], ["unused"])
            self.assertEqual(len(result["inputs_sha256"]), 2)
            self.assertEqual(len(result["script_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
