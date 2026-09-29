from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from build_independent_screening_packets import build_packets
from import_independent_screening_packets import import_packets


FIELDS = [
    "screening_id", "pmid", "doi", "year", "title", "abstract", "query_id",
    "duplicate_count", "duplicate_sources", "title_abstract_decision_r1",
    "title_abstract_reason_r1", "reviewer_1", "title_abstract_decision_r2",
    "title_abstract_reason_r2", "reviewer_2", "adjudication", "adjudicator", "notes",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


class IndependentScreeningPacketTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.queue = self.root / "queue.tsv"
        self.protocol = self.root / "protocol.md"
        self.packet_dir = self.root / "packets"
        self.protocol.write_text("Frozen test protocol\n", encoding="utf-8")
        rows = []
        for index in range(3):
            rows.append({field: "" for field in FIELDS} | {
                "screening_id": f"test:{index}", "pmid": str(100 + index),
                "year": "2020", "title": f"Sleep and frailty {index}",
                "abstract": "Example abstract", "query_id": "test-query",
                "duplicate_count": "1", "duplicate_sources": "test",
            })
        write_tsv(self.queue, FIELDS, rows)
        self.manifest = build_packets(self.queue, self.protocol, self.packet_dir, batch_size=2)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_packets_are_separate_and_complete(self) -> None:
        one, fields_one = read_tsv(self.packet_dir / "reviewer_1/batch_0001.tsv")
        two, fields_two = read_tsv(self.packet_dir / "reviewer_2/batch_0001.tsv")
        self.assertEqual(len(one), 2)
        self.assertEqual(len(two), 2)
        self.assertEqual(fields_one, fields_two)
        self.assertNotIn("title_abstract_decision_r2", fields_one)
        self.assertTrue(all(row["decision"] == "" for row in one + two))
        self.assertEqual(self.manifest["record_count"], 3)
        self.assertEqual(self.manifest["reviewers"]["reviewer_1"]["packet_count"], 2)

    def test_dry_run_and_atomic_import(self) -> None:
        packet = self.packet_dir / "reviewer_1/batch_0001.tsv"
        rows, fields = read_tsv(packet)
        rows[0].update(decision="include", reason="eligible", reviewer="Reviewer A")
        write_tsv(packet, fields, rows)

        result = import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")
        self.assertEqual(result["imported_decisions"], 1)
        self.assertFalse(result["applied"])
        queue_rows, _ = read_tsv(self.queue)
        self.assertEqual(queue_rows[0]["title_abstract_decision_r1"], "")

        result = import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1", apply=True)
        self.assertTrue(result["applied"])
        self.assertTrue(Path(result["backup_path"]).is_file())
        queue_rows, _ = read_tsv(self.queue)
        self.assertEqual(queue_rows[0]["title_abstract_decision_r1"], "include")
        self.assertEqual(queue_rows[0]["reviewer_1"], "Reviewer A")

        result = import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1", apply=True)
        self.assertEqual(result["already_present"], 1)
        self.assertFalse(result["applied"])

    def test_export_refuses_to_discard_existing_decisions(self) -> None:
        rows, fields = read_tsv(self.queue)
        rows[0]["title_abstract_decision_r1"] = "include"
        write_tsv(self.queue, fields, rows)
        with self.assertRaisesRegex(ValueError, "existing review data"):
            build_packets(self.queue, self.protocol, self.root / "new_packets", batch_size=2)

    def test_reviewer_identity_must_be_consistent_across_batches(self) -> None:
        first = self.packet_dir / "reviewer_1/batch_0001.tsv"
        rows, fields = read_tsv(first)
        rows[0].update(decision="include", reviewer="Reviewer A")
        write_tsv(first, fields, rows)
        second = self.packet_dir / "reviewer_1/batch_0002.tsv"
        rows, fields = read_tsv(second)
        rows[0].update(decision="exclude", reviewer="Reviewer B")
        write_tsv(second, fields, rows)
        with self.assertRaisesRegex(ValueError, "reviewer identity changed within reviewer_1"):
            import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")

    def test_reviewer_slots_must_have_distinct_global_identities(self) -> None:
        queue_rows, queue_fields = read_tsv(self.queue)
        queue_rows[0].update(title_abstract_decision_r1="include", reviewer_1="Reviewer A")
        write_tsv(self.queue, queue_fields, queue_rows)
        packet = self.packet_dir / "reviewer_2/batch_0002.tsv"
        rows, fields = read_tsv(packet)
        rows[0].update(decision="exclude", reviewer="Reviewer A")
        write_tsv(packet, fields, rows)
        with self.assertRaisesRegex(ValueError, "two independent decisions name the same reviewer"):
            import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_2")

    def test_protocol_change_invalidates_import(self) -> None:
        self.protocol.write_text("Amended protocol\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "protocol changed"):
            import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")


if __name__ == "__main__":
    unittest.main()
