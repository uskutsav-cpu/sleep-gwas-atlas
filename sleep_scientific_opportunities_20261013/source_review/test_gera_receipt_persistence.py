"""Focused regression for the observed internal-ENOSPC receipt loss; no GWAS IO."""
import errno
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from acquire_gera_raw import persist_receipt


class ReceiptPersistenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(
            prefix='receipt-regression-',
            dir='/Volumes/Extreme SSD/sleep_scientific_opportunities_20261013/source_review')
        self.addCleanup(self.tmp.cleanup)
        self.ssd = Path(self.tmp.name) / 'ssd.json'
        self.internal = Path(self.tmp.name) / 'internal.json'
        self.receipt = {'status': 'FAILED_RESUME_CHECKPOINTS_PRESERVED',
                        'stop_reason': 'RESUME_INTERNAL_FLOOR',
                        'actual_md5': None, 'h2_or_rg_runs': 0}

    def test_normal_mirror_is_byte_identical(self):
        self.assertTrue(persist_receipt(self.receipt, self.internal, self.ssd))
        self.assertEqual(self.ssd.read_bytes(), self.internal.read_bytes())
        self.assertEqual(json.loads(self.ssd.read_text()), self.receipt)

    def test_internal_enospc_at_open_and_write_preserves_ssd(self):
        real_open = Path.open
        for failure_stage in ['open', 'write']:
            with self.subTest(failure_stage=failure_stage):
                ssd = self.ssd.with_name(failure_stage + '-ssd.json')
                internal = self.internal.with_name(failure_stage + '-internal.json')

                class FullInternalStream:
                    def __enter__(self):
                        real_open(internal, 'x').close()
                        return self

                    def __exit__(self, *args):
                        return False

                    def write(self, payload):
                        raise OSError(errno.ENOSPC, 'simulated full internal disk')

                def affected_open(path, *args, **kwargs):
                    if path == internal:
                        if failure_stage == 'open':
                            raise OSError(errno.ENOSPC, 'simulated full internal disk')
                        return FullInternalStream()
                    return real_open(path, *args, **kwargs)

                with patch.object(Path, 'open', affected_open):
                    self.assertFalse(persist_receipt(self.receipt, internal, ssd))
                self.assertEqual(json.loads(ssd.read_text()), self.receipt)
                self.assertEqual(self.receipt['status'], 'FAILED_RESUME_CHECKPOINTS_PRESERVED')

    def test_existing_receipts_are_never_overwritten(self):
        self.ssd.write_text('preserved SSD evidence')
        with self.assertRaises(FileExistsError):
            persist_receipt(self.receipt, self.internal, self.ssd)
        self.assertEqual(self.ssd.read_text(), 'preserved SSD evidence')
        self.assertFalse(self.internal.exists())
        self.ssd.unlink()
        self.internal.write_text('preserved internal evidence')
        with self.assertRaises(FileExistsError):
            persist_receipt(self.receipt, self.internal, self.ssd)
        self.assertEqual(self.internal.read_text(), 'preserved internal evidence')
        self.assertEqual(json.loads(self.ssd.read_text()), self.receipt)

    def test_other_internal_io_failures_remain_errors(self):
        real_open = Path.open

        def denied_open(path, *args, **kwargs):
            if path == self.internal:
                raise OSError(errno.EACCES, 'simulated permission failure')
            return real_open(path, *args, **kwargs)

        with patch.object(Path, 'open', denied_open):
            with self.assertRaises(PermissionError):
                persist_receipt(self.receipt, self.internal, self.ssd)
        self.assertEqual(json.loads(self.ssd.read_text()), self.receipt)


if __name__ == '__main__':
    unittest.main()
