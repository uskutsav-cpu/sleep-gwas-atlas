from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/46_build_hlma_portal_archive_crosswalk.py"
SPEC = importlib.util.spec_from_file_location("hlma_portal_archive_crosswalk", SCRIPT)
BUILDER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(BUILDER)


class HLMAPortalArchiveCrosswalkTests(unittest.TestCase):
    def test_archive_id_extracts_official_file_path_id(self):
        url = "https://download.cncb.ac.cn/OMIX/OMIX006021/OMIX006021-05.fragments.tar.gz"
        self.assertEqual(BUILDER.archive_id({"file_path": url}), ("OMIX006021-05", url))

    def test_archive_id_rejects_unexpected_host(self):
        with self.assertRaisesRegex(ValueError, "unexpected official archive URL"):
            BUILDER.archive_id({"file_path": "https://example.org/OMIX006021-05.tar.gz"})

    def test_canonical_note_removes_old_generated_suffix(self):
        prefix = "Observed on public listing; no archive checksum acquired."
        old = " Official portal API snapshot maps this listing row to OMIX006021-05; bytes unverified."
        result = BUILDER.canonical_note(prefix + old, "OMIX006021-05")
        self.assertTrue(result.startswith(prefix + " Official portal API snapshot reports archive path OMIX006021-05;"))
        self.assertNotIn("maps this listing row", result)
        self.assertEqual(BUILDER.canonical_note(result, "OMIX006021-05"), result)


if __name__ == "__main__":
    unittest.main()
