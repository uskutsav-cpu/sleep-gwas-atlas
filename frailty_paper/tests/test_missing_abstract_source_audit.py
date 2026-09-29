from __future__ import annotations

import csv
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/47_audit_missing_abstract_sources.py"
SPEC = importlib.util.spec_from_file_location("missing_abstract_source_audit", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MissingAbstractSourceAuditTests(unittest.TestCase):
    def make_fixture(self, second_record: str = "") -> tuple[Path, Path]:
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        review = Path(temp.name) / "review"
        xml_dir = review / "pubmed" / "query"
        xml_dir.mkdir(parents=True)
        (review / "screening").mkdir()
        xml = (
            "<PubmedArticleSet>"
            "<PubmedArticle><MedlineCitation><PMID>123</PMID><Article>"
            "<ArticleTitle>Title one</ArticleTitle></Article></MedlineCitation></PubmedArticle>"
            f"{second_record}</PubmedArticleSet>"
        ).encode()
        xml_path = xml_dir / "batch.xml"
        xml_path.write_bytes(xml)
        with (review / "pubmed_source_files.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["source_file", "bytes", "sha256"], delimiter="\t")
            writer.writeheader()
            writer.writerow({
                "source_file": "pubmed/query/batch.xml",
                "bytes": len(xml),
                "sha256": hashlib.sha256(xml).hexdigest(),
            })
        with (review / "screening/title_abstract_queue.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["pmid", "abstract"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"pmid": "123", "abstract": ""})
        return review, SCRIPT

    def test_passes_when_empty_abstract_is_absent_from_all_source_occurrences(self):
        second = (
            "<PubmedArticle><MedlineCitation><PMID>123</PMID><Article>"
            "<ArticleTitle>Duplicate title</ArticleTitle></Article></MedlineCitation></PubmedArticle>"
        )
        review, script = self.make_fixture(second)
        result = MODULE.audit(review, script)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["empty_abstract_pmids_found_in_source_xml"], 1)
        self.assertEqual(result["pmids_with_duplicate_source_occurrences"], 1)
        self.assertEqual(result["pmids_with_abstract_in_any_source_occurrence"], [])

    def test_fails_if_any_duplicate_source_has_an_abstract(self):
        second = (
            "<PubmedArticle><MedlineCitation><PMID>123</PMID><Article><Abstract>"
            "<AbstractText>Recovered abstract</AbstractText></Abstract>"
            "</Article></MedlineCitation></PubmedArticle>"
        )
        review, script = self.make_fixture(second)
        result = MODULE.audit(review, script)
        self.assertEqual(result["status"], "FAIL")
        self.assertEqual(result["pmids_with_abstract_in_any_source_occurrence"], ["123"])

    def test_explicit_external_root_allows_a_symlinked_pubmed_cache(self):
        review, script = self.make_fixture()
        local_cache = review / "pubmed"
        external_cache = review.parent / "external_cache"
        local_cache.rename(external_cache)
        local_cache.symlink_to(external_cache, target_is_directory=True)

        result = MODULE.audit(review, script, external_cache)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["allowed_external_storage_root"], str(external_cache.resolve()))


if __name__ == "__main__":
    unittest.main()
