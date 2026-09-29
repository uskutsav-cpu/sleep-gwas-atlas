from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/55_build_pubmed_window_review_records.py"
SPEC = importlib.util.spec_from_file_location("pubmed_window_review_records", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def article(pmid: str, title: str, doi: str = "") -> str:
    doi_node = f'<ArticleId IdType="doi">{doi}</ArticleId>' if doi else ""
    return (
        "<PubmedArticle><MedlineCitation><PMID>" + pmid + "</PMID><Article>"
        "<ArticleTitle>" + title + "</ArticleTitle><Abstract><AbstractText>Abstract "
        + pmid + "</AbstractText></Abstract><Journal><Title>Test Journal</Title>"
        "<JournalIssue><PubDate><Year>2026</Year></PubDate></JournalIssue>"
        "</Journal></Article></MedlineCitation><PubmedData><ArticleIdList>"
        + doi_node + "</ArticleIdList></PubmedData></PubmedArticle>"
    )


class PubmedWindowReviewRecordsTests(unittest.TestCase):
    def test_build_keeps_increment_separate_and_preserves_all_occurrences(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            increment = root / "increment"
            primary = increment / "primary_sleep_frailty/segment_2026-09-23_2026-09-25/batch_000000.xml"
            neighborhood = increment / "frailty_neighborhood/segment_2026-09-23_2026-09-25/batch_000000.xml"
            primary.parent.mkdir(parents=True)
            neighborhood.parent.mkdir(parents=True)
            primary.write_text(
                "<PubmedArticleSet>" + article("100", "Shared article")
                + article("200", "In prior snapshot") + "</PubmedArticleSet>"
            )
            neighborhood.write_text(
                "<PubmedArticleSet>" + article("100", "Shared article")
                + article("300", "New neighborhood article") + "</PubmedArticleSet>"
            )
            (increment / "window_manifest.json").write_text(json.dumps({
                "start_date": "2026-09-23", "end_date": "2026-09-25",
                "records_returned_across_queries": 4,
                "integration_status": "SEPARATE_CACHE_NOT_MERGED",
                "queries": [
                    {"query_id": "primary_sleep_frailty", "fetched": 2},
                    {"query_id": "genetic_sleep_frailty", "fetched": 0},
                    {"query_id": "frailty_neighborhood", "fetched": 2},
                ],
            }))
            reconciliation = root / "reconciliation.json"
            reconciliation.write_text(json.dumps({
                "checked_at_utc": "2026-09-25T00:00:00Z", "query_errors": 0,
                "missing_segment_id_occurrences": 0, "extra_segment_id_occurrences": 0,
                "pubmed_article_nodes": 4, "pubmed_book_article_nodes": 0,
            }))

            review = root / "review"
            review.mkdir()
            all_records = review / "all_records.tsv"
            with all_records.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=["pmid", "doi", "normalized_title"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"pmid": "200", "doi": "", "normalized_title": "in prior snapshot"})
            output_manifest = review / "record_build_outputs.tsv"
            with output_manifest.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=["file", "bytes", "sha256"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"file": "all_records.tsv", "bytes": all_records.stat().st_size,
                                 "sha256": hashlib.sha256(all_records.read_bytes()).hexdigest()})
            (review / "record_build_manifest.json").write_text(json.dumps({
                "xml_records": 1, "manual_import_records": 0,
                "total_source_records": 1, "unique_records": 1,
            }))
            frozen_hash = hashlib.sha256(all_records.read_bytes()).hexdigest()
            outdir = review / "reconciliation_2026-09-25_window"

            result = MODULE.build(increment, review, outdir, reconciliation,
                                  "pubmed_cache_increments/2026-09-23_2026-09-25")

            with (outdir / "pubmed_window_occurrences.tsv").open(encoding="utf-8", newline="") as stream:
                occurrences = list(csv.DictReader(stream, delimiter="\t"))
            with (outdir / "pubmed_window_unique_review_queue.tsv").open(encoding="utf-8", newline="") as stream:
                unique = list(csv.DictReader(stream, delimiter="\t"))
            raw_queue = (outdir / "pubmed_window_unique_review_queue.tsv").read_bytes()
            self.assertTrue(all(not line.endswith(b"\t") for line in raw_queue.splitlines()))
            self.assertEqual(len(occurrences), 4)
            self.assertEqual(len(unique), 3)
            self.assertEqual(result["duplicate_occurrences_within_window"], 1)
            self.assertEqual(result["prior_snapshot_matches"], 1)
            self.assertEqual(sum(row["window_duplicate_count"] == "2" for row in occurrences), 2)
            self.assertEqual(sum(row["frozen_snapshot_match"] == "YES" for row in unique), 1)
            self.assertTrue(all(row["screening_status"].endswith("UNSCREENED") for row in unique))
            self.assertTrue(all(not row["title_abstract_decision_r1"] and not row["title_abstract_decision_r2"]
                                and not row["adjudication"] for row in unique))
            self.assertEqual(hashlib.sha256(all_records.read_bytes()).hexdigest(), frozen_hash)


if __name__ == "__main__":
    unittest.main()
