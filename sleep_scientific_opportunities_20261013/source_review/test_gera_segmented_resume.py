"""Small synthetic tests for retained segments and aggregate resume storage."""
import hashlib
from pathlib import Path
import tempfile
import unittest

from acquire_gera_raw import (digest_prefix_segments, preserved_staging_bytes,
                              resume_staging_bytes)


class SegmentedResumeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(
            prefix='segmented-resume-',
            dir='/Volumes/Extreme SSD/sleep_scientific_opportunities_20261013/source_review')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def _segment(self, name, content, start):
        path = self.root / name
        path.write_bytes(content)
        return {'path': str(path), 'bytes': len(content),
                'sha256': hashlib.sha256(content).hexdigest(),
                'range_start': start, 'range_end': start + len(content) - 1,
                'mtime_ns': path.stat().st_mtime_ns}

    def _sidecar_bytes(self, paths):
        return sum((path.parent / ('._' + path.name)).stat().st_size
                   for path in paths
                   if (path.parent / ('._' + path.name)).exists())

    def test_prefix_hashes_segments_in_order_without_concatenating_them(self):
        first = self._segment('first.part', b'header\nrow1\n', 0)
        second = self._segment('second.part', b'row2\nrow3\n', first['bytes'])
        pieces = [first, second]
        md5, sha, identities = digest_prefix_segments(pieces, float('inf'), 0)
        expected = b'header\nrow1\nrow2\nrow3\n'
        self.assertEqual(md5, hashlib.md5(expected).hexdigest())
        self.assertEqual(sha, hashlib.sha256(expected).hexdigest())
        self.assertEqual([item['sha256'] for item in identities],
                         [first['sha256'], second['sha256']])
        self.assertFalse((self.root / 'combined.part').exists())

    def test_prefix_length_or_segment_digest_change_fails_closed(self):
        segment = self._segment('segment.part', b'original', 0)
        (self.root / 'segment.part').write_bytes(b'changed!')
        with self.assertRaisesRegex(RuntimeError, 'IDENTITY_CHANGED'):
            digest_prefix_segments([segment], float('inf'), 0)

    def test_staging_budget_counts_every_preserved_piece_and_reserve(self):
        first = self._segment('first.part', b'123456', 0)
        second = self._segment('second.part', b'1234', 6)
        paths = [Path(first['path']), Path(second['path'])]
        total = 10 + self._sidecar_bytes(paths)
        self.assertEqual(preserved_staging_bytes([first, second], total + 2, 2), total)
        with self.assertRaisesRegex(RuntimeError, 'GLOBAL_CAP'):
            preserved_staging_bytes([first, second], total + 1, 2)

    def test_staging_inventory_rejects_changed_file_and_duplicate(self):
        segment = self._segment('segment.part', b'123456', 0)
        segment['mtime_ns'] += 1
        with self.assertRaisesRegex(RuntimeError, 'INVENTORY_CHANGED'):
            preserved_staging_bytes([segment], 100, 0)
        segment['mtime_ns'] = (self.root / 'segment.part').stat().st_mtime_ns
        with self.assertRaisesRegex(RuntimeError, 'DUPLICATE'):
            preserved_staging_bytes([segment, segment], 100, 0)

    def test_live_budget_counts_segments_download_candidate_and_metadata(self):
        retained = self._segment('retained.part', b'123456', 0)
        stage = self.root / 'stage'; stage.mkdir()
        suffix = stage / 'suffix.part'; suffix.write_bytes(b'1234')
        headers = stage / 'HTTP_headers_private.txt'; headers.write_bytes(b'h')
        retained_path = Path(retained['path'])
        sidecars = self._sidecar_bytes([retained_path, suffix, headers])
        metadata = sidecars + headers.stat().st_size
        expected_total = 6 + 4 + headers.stat().st_size + sidecars
        total = resume_staging_bytes([retained], [suffix], stage,
                                     expected_total + 1, metadata, 'receipt.json')
        self.assertEqual(total, expected_total)
        with self.assertRaisesRegex(RuntimeError, 'SIZE_LIMIT'):
            resume_staging_bytes([retained], [suffix], stage,
                                 expected_total - 1, metadata, 'receipt.json')
        if metadata:
            with self.assertRaisesRegex(RuntimeError, 'METADATA_RESERVE'):
                resume_staging_bytes([retained], [suffix], stage,
                                     total + 1, metadata - 1, 'receipt.json')

    def test_live_budget_rejects_unexpected_files_and_metadata_overrun(self):
        stage = self.root / 'stage'; stage.mkdir()
        suffix = stage / 'suffix.part'; suffix.write_bytes(b'1234')
        surprise = stage / 'unregistered.part'; surprise.write_bytes(b'x')
        with self.assertRaisesRegex(RuntimeError, 'UNEXPECTED'):
            resume_staging_bytes([], [suffix], stage, 100, 10, 'receipt.json')
        surprise.unlink()
        large_header = stage / 'HTTP_headers_private.txt'
        large_header.write_bytes(b'x' * 11)
        with self.assertRaisesRegex(RuntimeError, 'METADATA_RESERVE'):
            resume_staging_bytes([], [suffix], stage, 100, 10, 'receipt.json')


if __name__ == '__main__':
    unittest.main()
