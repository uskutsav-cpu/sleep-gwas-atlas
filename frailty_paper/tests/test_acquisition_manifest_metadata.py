from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/45_audit_acquisition_manifest_metadata.py"
SPEC = importlib.util.spec_from_file_location("acquisition_manifest_metadata", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(AUDIT)
REQUIRED_COLUMNS = AUDIT.REQUIRED_COLUMNS


class AcquisitionManifestMetadataTests(unittest.TestCase):
    def write_manifest(self, rows, columns=REQUIRED_COLUMNS):
        directory = tempfile.TemporaryDirectory()
        path = Path(directory.name) / "manifest.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return directory, path

    def valid_row(self, **updates):
        row = {column: "UNKNOWN" for column in REQUIRED_COLUMNS}
        row.update({"resource_id": "resource-1", "file": "data.gz", "bytes": "1", "sha256": "a" * 64})
        row.update(updates)
        return row

    def test_accepts_complete_metadata_with_explicit_unknowns(self):
        directory, path = self.write_manifest([self.valid_row()])
        with directory:
            count, errors = AUDIT.audit_manifest(path)
        self.assertEqual(count, 1)
        self.assertEqual(errors, [])

    def test_rejects_missing_required_metadata(self):
        row = self.valid_row()
        row["genome_build"] = ""
        directory, path = self.write_manifest([row])
        with directory:
            _, errors = AUDIT.audit_manifest(path)
        self.assertIn("line 2: empty required field genome_build", errors)

    def test_rejects_duplicate_resource_ids(self):
        directory, path = self.write_manifest([self.valid_row(), self.valid_row()])
        with directory:
            _, errors = AUDIT.audit_manifest(path)
        self.assertIn("line 3: duplicate resource_id resource-1", errors)

    def test_rejects_zero_bytes_and_invalid_checksum(self):
        directory, path = self.write_manifest([self.valid_row(bytes="0", sha256="bad")])
        with directory:
            _, errors = AUDIT.audit_manifest(path)
        self.assertIn("line 2: byte count must be positive", errors)
        self.assertIn("line 2: invalid SHA-256", errors)

    def test_rejects_missing_columns(self):
        columns = [name for name in REQUIRED_COLUMNS if name != "license"]
        directory, path = self.write_manifest([], columns=columns)
        with directory:
            count, errors = AUDIT.audit_manifest(path)
        self.assertEqual(count, 0)
        self.assertEqual(errors, ["missing required columns: license"])

    def test_rejects_empty_manifest(self):
        directory, path = self.write_manifest([])
        with directory:
            count, errors = AUDIT.audit_manifest(path)
        self.assertEqual(count, 0)
        self.assertEqual(errors, ["manifest has no resource rows"])


if __name__ == "__main__":
    unittest.main()
