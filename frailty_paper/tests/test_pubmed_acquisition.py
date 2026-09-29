from __future__ import annotations

import importlib.util
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import xml.etree.ElementTree as ET

import requests


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/08_build_review_records.py'
SPEC = importlib.util.spec_from_file_location('build_review_records', SCRIPT)
REVIEW = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(REVIEW)

SEARCH_SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/01_pubmed_search.py'
SEARCH_SPEC = importlib.util.spec_from_file_location('pubmed_search', SEARCH_SCRIPT)
SEARCH = importlib.util.module_from_spec(SEARCH_SPEC)
assert SEARCH_SPEC and SEARCH_SPEC.loader
SEARCH_SPEC.loader.exec_module(SEARCH)

ZENODO_SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/03_download_zenodo_record.py'
ZENODO_SPEC = importlib.util.spec_from_file_location('download_zenodo_record', ZENODO_SCRIPT)
ZENODO = importlib.util.module_from_spec(ZENODO_SPEC)
assert ZENODO_SPEC and ZENODO_SPEC.loader
ZENODO_SPEC.loader.exec_module(ZENODO)

IMPORT_SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/15_import_manual_review_exports.py'
IMPORT_SPEC = importlib.util.spec_from_file_location('import_manual_review_exports', IMPORT_SCRIPT)
IMPORTER = importlib.util.module_from_spec(IMPORT_SPEC)
assert IMPORT_SPEC and IMPORT_SPEC.loader
IMPORT_SPEC.loader.exec_module(IMPORTER)

VERIFY_SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/16_verify_acquired_resource_manifest.py'
VERIFY_SPEC = importlib.util.spec_from_file_location('verify_acquired_resource_manifest', VERIFY_SCRIPT)
VERIFY = importlib.util.module_from_spec(VERIFY_SPEC)
assert VERIFY_SPEC and VERIFY_SPEC.loader
VERIFY_SPEC.loader.exec_module(VERIFY)

SCREEN_SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/17_validate_review_screening.py'
SCREEN_SPEC = importlib.util.spec_from_file_location('validate_review_screening', SCREEN_SCRIPT)
SCREEN = importlib.util.module_from_spec(SCREEN_SPEC)
assert SCREEN_SPEC and SCREEN_SPEC.loader
SCREEN_SPEC.loader.exec_module(SCREEN)


class PubmedBookTitleTests(unittest.TestCase):
    def parse(self, xml: str) -> dict[str, str]:
        base = Path('/tmp/pubmed')
        source = base / 'frailty_neighborhood' / 'batch_000000.xml'
        return REVIEW.parse_article(ET.fromstring(xml), source, base, 1)

    def test_book_article_uses_chapter_title_and_collection(self) -> None:
        row = self.parse('''<PubmedBookArticle><BookDocument><PMID>1</PMID>
          <ArticleTitle>Chapter title</ArticleTitle><Book><BookTitle>Book title</BookTitle>
          <CollectionTitle>Series title</CollectionTitle></Book></BookDocument></PubmedBookArticle>''')
        self.assertEqual(row['title'], 'Chapter title')
        self.assertEqual(row['journal'], 'Series title')

    def test_book_title_fills_missing_chapter_title(self) -> None:
        row = self.parse('''<PubmedBookArticle><BookDocument><PMID>2</PMID><Book>
          <BookTitle>CADTH reimbursement report title</BookTitle>
          <CollectionTitle>CADTH Reviews</CollectionTitle></Book></BookDocument></PubmedBookArticle>''')
        self.assertEqual(row['title'], 'CADTH reimbursement report title')
        self.assertEqual(row['journal'], 'CADTH Reviews')

    def test_journal_article_title_is_unchanged(self) -> None:
        row = self.parse('''<PubmedArticle><MedlineCitation><PMID>3</PMID><Article>
          <ArticleTitle>Journal article title</ArticleTitle><Journal><Title>Journal name</Title></Journal>
          </Article></MedlineCitation></PubmedArticle>''')
        self.assertEqual(row['title'], 'Journal article title')
        self.assertEqual(row['journal'], 'Journal name')


class PubmedXmlDiscoveryTests(unittest.TestCase):
    def test_appledouble_sidecars_are_not_treated_as_xml_records(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            real = root / 'batch_000000.xml'
            sidecar = root / '._batch_000000.xml'
            real.write_text('<PubmedArticleSet/>', encoding='utf-8')
            sidecar.write_bytes(b'AppleDouble metadata')
            self.assertEqual(REVIEW.pubmed_xml_files(root), [real])

    def test_source_locator_uses_logical_cache_path(self) -> None:
        physical_root = Path('/Volumes/example/cache/pubmed_cache')
        xml = physical_root / 'frailty_neighborhood/segment_x/batch_000000.xml'
        self.assertEqual(
            REVIEW.pubmed_source_path(xml, physical_root, Path('pubmed')),
            'pubmed/frailty_neighborhood/segment_x/batch_000000.xml',
        )


class PubmedRetryTests(unittest.TestCase):
    def test_truncated_transport_response_is_retried(self) -> None:
        class Response:
            status_code = 200

        class Session:
            def __init__(self) -> None:
                self.calls = 0

            def get(self, *args, **kwargs):
                self.calls += 1
                if self.calls == 1:
                    raise requests.exceptions.ChunkedEncodingError('response ended prematurely')
                return Response()

        session = Session()
        with patch.object(SEARCH.time, 'sleep'):
            result = SEARCH.req(session, 'efetch.fcgi', {'db': 'pubmed'})
        self.assertIsInstance(result, Response)
        self.assertEqual(session.calls, 2)

    def test_esearch_retries_http_200_search_backend_outage(self) -> None:
        class Response:
            def __init__(self, payload):
                self.payload = payload

            def json(self):
                return self.payload

        session = object()
        with patch.object(SEARCH, 'req', side_effect=[
            Response({'esearchresult': {'ERROR': 'Search Backend failed: temporarily unavailable'}}),
            Response({'esearchresult': {'count': '2', 'idlist': ['1', '2']}}),
        ]) as request, patch.object(SEARCH.time, 'sleep') as sleep:
            result = SEARCH.esearch(session, 'sleep', '', '')
        self.assertEqual(result['count'], '2')
        self.assertEqual(request.call_count, 2)
        sleep.assert_called_once_with(1)

    def test_esearch_raises_clear_error_after_retryable_outages(self) -> None:
        class Response:
            def json(self):
                return {'esearchresult': {'ERROR': 'Search Backend failed: temporarily unavailable'}}

        with patch.object(SEARCH, 'req', return_value=Response()) as request, \
                patch.object(SEARCH.time, 'sleep') as sleep:
            with self.assertRaisesRegex(RuntimeError, 'PubMed ESearch failed: Search Backend failed'):
                SEARCH.esearch(object(), 'sleep', '', '')
        self.assertEqual(request.call_count, 6)
        self.assertEqual(sleep.call_count, 5)


class PubmedBatchCheckpointTests(unittest.TestCase):
    XML = b'''<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>123</PMID>
      <Article><ArticleTitle>Cached result</ArticleTitle><Journal><Title>Test</Title>
      </Journal></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>'''
    XML_124 = b'''<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>124</PMID>
      <Article><ArticleTitle>Refreshed result</ArticleTitle><Journal><Title>Test</Title>
      </Journal></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>'''

    def test_complete_cached_batch_is_reused_without_fetch(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            segment = Path(temp) / 'segment_2020-01-01_2020-12-31'
            segment.mkdir()
            (segment / 'batch_000000.xml').write_bytes(self.XML)
            rows = []
            with patch.object(SEARCH, 'esearch', side_effect=[
                {'count': '1', 'webenv': 'cached', 'querykey': '1'},
                {'count': '1', 'idlist': ['123']},
            ]), patch.object(SEARCH, 'req', side_effect=AssertionError('unexpected fetch')):
                count = SEARCH.fetch_segment(
                    object(), 'primary', 'sleep', SEARCH.dt.date(2020, 1, 1),
                    SEARCH.dt.date(2020, 12, 31), Path(temp), '', '', rows
                )
            self.assertEqual(count, 1)
            self.assertEqual([row['pmid'] for row in rows], ['123'])
            self.assertEqual((segment / 'esearch_ids.txt').read_text(encoding='utf-8'), '123\n')

    def test_same_count_with_changed_ids_refetches_cached_batch(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            segment = Path(temp) / 'segment_2020-01-01_2020-12-31'
            segment.mkdir()
            target = segment / 'batch_000000.xml'
            target.write_bytes(self.XML)

            class Response:
                content = PubmedBatchCheckpointTests.XML_124

            rows = []
            with patch.object(SEARCH, 'esearch', side_effect=[
                {'count': '1', 'webenv': 'fresh', 'querykey': '1'},
                {'count': '1', 'idlist': ['124']},
            ]), patch.object(SEARCH, 'req', return_value=Response()) as request:
                SEARCH.fetch_segment(
                    object(), 'primary', 'sleep', SEARCH.dt.date(2020, 1, 1),
                    SEARCH.dt.date(2020, 12, 31), Path(temp), '', '', rows
                )
            self.assertEqual(request.call_count, 1)
            self.assertEqual([row['pmid'] for row in rows], ['124'])
            self.assertEqual(SEARCH.parse_batch(target.read_bytes(), 'primary', 1)[0]['pmid'], '124')

    def test_batch_parser_retains_pubmed_book_articles(self) -> None:
        xml = b'''<PubmedArticleSet><PubmedBookArticle><BookDocument><PMID>456</PMID>
          <ArticleTitle>Book chapter</ArticleTitle><Book><BookTitle>Reference book</BookTitle>
          <PubDate><Year>2024</Year></PubDate></Book></BookDocument></PubmedBookArticle>
          </PubmedArticleSet>'''
        rows = SEARCH.parse_batch(xml, 'neighborhood', 1)
        self.assertEqual(rows[0]['pmid'], '456')
        self.assertEqual(rows[0]['title'], 'Book chapter')
        self.assertEqual(rows[0]['year'], '2024')

    def test_truncated_cached_batch_is_replaced_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            segment = Path(temp) / 'segment_2020-01-01_2020-12-31'
            segment.mkdir()
            target = segment / 'batch_000000.xml'
            target.write_bytes(b'<PubmedArticleSet><broken>')

            class Response:
                content = PubmedBatchCheckpointTests.XML

            rows = []
            with patch.object(SEARCH, 'esearch', side_effect=[
                {'count': '1', 'webenv': 'fresh', 'querykey': '1'},
                {'count': '1', 'idlist': ['123']},
            ]), patch.object(SEARCH, 'req', return_value=Response()):
                count = SEARCH.fetch_segment(
                    object(), 'primary', 'sleep', SEARCH.dt.date(2020, 1, 1),
                    SEARCH.dt.date(2020, 12, 31), Path(temp), '', '', rows
                )
            self.assertEqual(count, 1)
            self.assertEqual([row['pmid'] for row in rows], ['123'])
            self.assertEqual(len(SEARCH.parse_batch(target.read_bytes(), 'primary', 1)), 1)
            self.assertFalse(target.with_suffix('.xml.tmp').exists())


class ReviewOutputProtectionTests(unittest.TestCase):
    def test_blank_generated_queues_allow_rebuild(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            (out / 'screening').mkdir()
            (out / 'screening/title_abstract_queue.tsv').write_text(
                'screening_id\ttitle\ttitle_abstract_decision_r1\tnotes\n', encoding='utf-8'
            )
            (out / 'screening/fulltext_queue.tsv').write_text('screening_id\ttitle\n', encoding='utf-8')
            REVIEW.guard_review_outputs(out)

    def test_reviewer_decision_prevents_destructive_rebuild(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            (out / 'screening').mkdir()
            (out / 'screening/title_abstract_queue.tsv').write_text(
                'screening_id\ttitle\ttitle_abstract_decision_r1\tnotes\n'
                'record-1\tExample\tinclude\tPossible eligible cohort\n', encoding='utf-8'
            )
            with self.assertRaisesRegex(SystemExit, 'Refusing to overwrite reviewer/extraction data'):
                REVIEW.guard_review_outputs(out)


class ManualExportImportTests(unittest.TestCase):
    def test_ris_import_retains_accession_database_and_metadata(self) -> None:
        payload = '''TY  - JOUR
TI  - Sleep and frailty in older adults
AU  - Example, A
PY  - 2020
JO  - Journal of Sleep
DO  - 10.1234/EXAMPLE.1
AN  - WOS:000123456
AB  - A study abstract.
ER  -
'''
        row = IMPORTER.parse_ris(payload, 'Web of Science', Path('wos.ris'), 'review/manual_exports/wos.ris')[0]
        self.assertEqual(row['source_record_id'], 'WOS:000123456')
        self.assertEqual(row['source_database'], 'Web of Science')
        self.assertEqual(row['year'], '2020')
        self.assertEqual(row['normalized_title'], 'sleep and frailty in older adults')
        self.assertEqual(row['import_status'], 'READY')

    def test_scopus_csv_import_and_missing_title_are_explicit(self) -> None:
        payload = b'Title,Authors,Year,Source title,Abstract,DOI,EID\n"Sleep, frailty",Example A,2021,Journal,Abstract,10.1234/a,2-s2.0-abc\n,Example B,2022,Journal,,10.1234/b,2-s2.0-def\n'
        rows = IMPORTER.parse_csv(payload, 'Scopus', Path('scopus.csv'), 'review/manual_exports/scopus.csv')
        self.assertEqual(rows[0]['source_record_id'], '2-s2.0-abc')
        self.assertEqual(rows[0]['title'], 'Sleep, frailty')
        self.assertEqual(rows[1]['import_status'], 'MISSING_TITLE')
        self.assertEqual(len(rows), 2)

    def test_export_import_does_not_change_original_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder = root / 'frailty_paper/review/manual_exports'
            folder.mkdir(parents=True)
            source = folder / 'embase_20260922.ris'
            original = b'TY  - JOUR\nTI  - Example\nAN  - E123\nER  -\n'
            source.write_bytes(original)
            rows, manifest = IMPORTER.import_exports(folder, root)
            self.assertEqual(len(rows), 1)
            self.assertEqual(manifest[0]['sha256'], hashlib.sha256(original).hexdigest())
            self.assertEqual(source.read_bytes(), original)

    def test_stale_normalized_import_is_rejected_without_raw_exports(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder = root / 'frailty_paper/review/manual_exports'
            folder.mkdir(parents=True)
            (folder / 'manual_normalized_records.tsv').write_text('record_id\n', encoding='utf-8')
            with patch.object(sys, 'argv', ['importer', '--repo', str(root)]):
                with self.assertRaisesRegex(SystemExit, 'No raw RIS/CSV exports found'):
                    IMPORTER.main()

    def test_manual_and_pubmed_records_deduplicate_by_doi(self) -> None:
        pubmed = {'record_id':'p1','pmid':'1','doi':'10.1234/example','normalized_title':'example'}
        manual = {'record_id':'m1','pmid':'','doi':'https://doi.org/10.1234/EXAMPLE','normalized_title':'example'}
        groups = REVIEW.group_duplicate_records([pubmed, manual])
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0][0], 'DOI')
        self.assertEqual({row['record_id'] for row in groups[0][1]}, {'p1', 'm1'})


class AcquiredResourceManifestTests(unittest.TestCase):
    def write_manifest(self, root: Path, file_value: str, payload: bytes) -> Path:
        path = root / file_value
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        manifest = root / 'manifest.tsv'
        manifest.write_text(
            'resource_id\tfile\tbytes\tsha256\n'
            f'resource-1\t{file_value}\t{len(payload)}\t{hashlib.sha256(payload).hexdigest()}\n',
            encoding='utf-8',
        )
        return manifest

    def test_manifest_verifies_file_size_and_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self.write_manifest(root, 'data/source.tsv.gz', b'verified source')
            count, errors = VERIFY.verify_manifest(root, manifest)
            self.assertEqual(count, 1)
            self.assertEqual(errors, [])

    def test_manifest_reports_hash_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self.write_manifest(root, 'data/source.tsv.gz', b'original')
            (root / 'data/source.tsv.gz').write_bytes(b'changed!')
            count, errors = VERIFY.verify_manifest(root, manifest)
            self.assertEqual(count, 1)
            self.assertTrue(any('SHA-256 mismatch' in error for error in errors))

    def test_manifest_rejects_paths_outside_repository(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = root / 'manifest.tsv'
            manifest.write_text('resource_id\tfile\tbytes\tsha256\nr1\t../outside\t0\t' + '0' * 64 + '\n', encoding='utf-8')
            count, errors = VERIFY.verify_manifest(root, manifest)
            self.assertEqual(count, 1)
            self.assertTrue(any('repo-relative' in error for error in errors))


class ReviewScreeningValidationTests(unittest.TestCase):
    def make_review(self, root: Path) -> Path:
        review = root / 'review'
        (review / 'screening').mkdir(parents=True)
        (review / 'screening/title_abstract_queue.tsv').write_text(
            'screening_id\ttitle_abstract_decision_r1\treviewer_1\ttitle_abstract_decision_r2\treviewer_2\tadjudication\tadjudicator\n'
            's1\tinclude\tA\tinclude\tB\t\t\n', encoding='utf-8')
        (review / 'screening/fulltext_queue.tsv').write_text(
            'screening_id\tfulltext_location\tfulltext_decision_r1\treviewer_1\tfulltext_decision_r2\treviewer_2\tadjudication\tadjudicator\n'
            's1\tfulltext/s1.pdf\texclude\tA\texclude\tB\t\t\n', encoding='utf-8')
        (review / 'screening/fulltext_exclusions.tsv').write_text(
            'screening_id\tprimary_exclusion_reason\treviewer\ns1\tNo eligible sleep phenotype\tA\n', encoding='utf-8')
        return review

    def test_consistent_dual_screening_and_reasoned_exclusion_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            counts, errors = SCREEN.validate_review(self.make_review(Path(temp)))
            self.assertEqual(counts['title_abstract_decided'], 1)
            self.assertEqual(counts['fulltext_exclusions'], 1)
            self.assertEqual(errors, [])

    def test_disagreement_requires_adjudication(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            review = self.make_review(Path(temp))
            path = review / 'screening/fulltext_queue.tsv'
            path.write_text(path.read_text().replace('exclude\tB\t\t\n', 'include\tB\t\t\n'), encoding='utf-8')
            _, errors = SCREEN.validate_review(review)
            self.assertTrue(any('requires adjudication' in error for error in errors))

    def test_unresolved_or_missing_reason_exclusion_is_flagged(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            review = self.make_review(Path(temp))
            path = review / 'screening/fulltext_exclusions.tsv'
            path.write_text('screening_id\tprimary_exclusion_reason\treviewer\ns1\t\tA\n', encoding='utf-8')
            _, errors = SCREEN.validate_review(review)
            self.assertTrue(any('primary exclusion reason required' in error for error in errors))

    def test_title_abstract_include_must_route_to_fulltext(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            review = self.make_review(Path(temp))
            path = review / 'screening/fulltext_queue.tsv'
            path.write_text('screening_id\tfulltext_location\tfulltext_decision_r1\treviewer_1\tfulltext_decision_r2\treviewer_2\tadjudication\tadjudicator\n', encoding='utf-8')
            (review / 'screening/fulltext_exclusions.tsv').write_text(
                'screening_id\tprimary_exclusion_reason\treviewer\n', encoding='utf-8')
            _, errors = SCREEN.validate_review(review)
            self.assertTrue(any('missing from the full-text queue' in error for error in errors))

    def test_conflicting_titles_sharing_a_doi_are_not_collapsed(self) -> None:
        records = [
            {'record_id':'p1','pmid':'1','doi':'10.1234/reused','normalized_title':'sleep duration review'},
            {'record_id':'p2','pmid':'2','doi':'10.1234/reused','normalized_title':'healthy aging program'},
        ]
        self.assertEqual(len(REVIEW.group_duplicate_records(records)), 2)

    def test_generic_matching_title_does_not_merge_distinct_identifiers(self) -> None:
        records = [
            {'record_id':'p1','pmid':'1','doi':'10.1234/one','normalized_title':'sleep apnea'},
            {'record_id':'p2','pmid':'2','doi':'10.1234/two','normalized_title':'sleep apnea'},
        ]
        self.assertEqual(len(REVIEW.group_duplicate_records(records)), 2)

    def test_fulltext_location_is_protected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            (out / 'screening').mkdir()
            (out / 'screening/fulltext_queue.tsv').write_text(
                'screening_id\tfulltext_location\nrecord-1\tfulltext/record-1.pdf\n', encoding='utf-8'
            )
            with self.assertRaisesRegex(SystemExit, 'Refusing to overwrite reviewer/extraction data'):
                REVIEW.guard_review_outputs(out)


class ZenodoDownloadTests(unittest.TestCase):
    def test_unrecognized_partial_is_quarantined_before_verified_promotion(self) -> None:
        payload = b'checksum-pinned source bytes'
        expected_md5 = hashlib.md5(payload).hexdigest()
        with tempfile.TemporaryDirectory() as temp:
            dest = Path(temp) / 'summary.tsv.gz'
            part = dest.with_name(dest.name + '.part')
            part.write_bytes(b'unrecognized old partial')

            def fake_curl(command, check):
                output = Path(command[command.index('--output') + 1])
                output.write_bytes(payload)

            with patch.object(ZENODO.subprocess, 'run', side_effect=fake_curl):
                ZENODO.download('https://zenodo.example/file', dest, len(payload), expected_md5)

            self.assertEqual(dest.read_bytes(), payload)
            self.assertFalse(part.exists())
            self.assertEqual(part.with_name(part.name + '.unverified').read_bytes(), b'unrecognized old partial')


if __name__ == '__main__':
    unittest.main()
