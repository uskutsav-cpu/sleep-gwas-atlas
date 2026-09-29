from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/56_import_pubmed_window_screening.py"
SPEC = importlib.util.spec_from_file_location("pubmed_window_screening_import", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class PubmedWindowScreeningImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.expected_hashes = (MODULE.EXPECTED_PACKET_MANIFEST_SHA256,
                                MODULE.EXPECTED_QUEUE_SHA256,
                                MODULE.EXPECTED_PROTOCOL_SHA256)
        self.queue = self.root / "queue.tsv"
        self.protocol = self.root / "protocol.md"
        self.packet_dir = self.root / "packets"
        self.protocol.write_text("Frozen protocol\n")
        self.source_fields = ["screening_id", "pmid", "title"]
        self.response_fields = ["decision", "reason", "reviewer", "notes"]
        self.packet_fields = [*self.source_fields, "packet_slot", *self.response_fields]
        self.queue_fields = [*self.source_fields, "screening_status",
                             "title_abstract_decision_r1", "title_abstract_reason_r1", "reviewer_1",
                             "title_abstract_decision_r2", "title_abstract_reason_r2", "reviewer_2",
                             "adjudication", "adjudicator", "notes"]
        self.base_rows = [
            {"screening_id": "sid-1", "pmid": "1", "title": "Title one",
             "screening_status": "SUPPLEMENTAL_OUT_OF_WINDOW_UNSCREENED"},
            {"screening_id": "sid-2", "pmid": "2", "title": "Title two",
             "screening_status": "SUPPLEMENTAL_OUT_OF_WINDOW_UNSCREENED"},
        ]
        for row in self.base_rows:
            for field in self.queue_fields:
                row.setdefault(field, "")
        write_tsv(self.queue, self.queue_fields, self.base_rows)
        self.blank_sha = MODULE.sha256_file(self.queue)
        reviewers = {}
        for slot in ("reviewer_1", "reviewer_2"):
            rows = []
            for source in self.base_rows:
                rows.append({**{field: source[field] for field in self.source_fields},
                             "packet_slot": slot, "decision": "", "reason": "",
                             "reviewer": "", "notes": ""})
            rel = Path(slot) / "batch_0001.tsv"
            write_tsv(self.packet_dir / rel, self.packet_fields, rows)
            reviewers[slot] = {"packet_count": 1, "packets": [{
                "path": rel.as_posix(), "row_count": len(rows),
                "screening_ids_sha256": MODULE.stable_rows_sha256(rows, ("screening_id",)),
                "source_rows_sha256": MODULE.stable_rows_sha256(rows, tuple(self.source_fields)),
                "template_sha256": MODULE.sha256_file(self.packet_dir / rel),
            }]}
        manifest = {
            "schema_version": 1, "queue_sha256": self.blank_sha,
            "protocol_sha256": MODULE.sha256_file(self.protocol), "record_count": 2,
            "source_fields": self.source_fields, "response_fields": self.response_fields,
            "reviewers": reviewers,
        }
        (self.packet_dir / "packet_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        MODULE.EXPECTED_PACKET_MANIFEST_SHA256 = MODULE.sha256_file(self.packet_dir / "packet_manifest.json")
        MODULE.EXPECTED_QUEUE_SHA256 = self.blank_sha
        MODULE.EXPECTED_PROTOCOL_SHA256 = MODULE.sha256_file(self.protocol)

    def tearDown(self) -> None:
        (MODULE.EXPECTED_PACKET_MANIFEST_SHA256,
         MODULE.EXPECTED_QUEUE_SHA256,
         MODULE.EXPECTED_PROTOCOL_SHA256) = self.expected_hashes
        self.temp.cleanup()

    def set_decisions(self, slot: str, reviewer: str, decisions: list[str]) -> None:
        path = self.packet_dir / slot / "batch_0001.tsv"
        rows, fields = MODULE.read_tsv(path)
        for row, decision in zip(rows, decisions):
            row.update(decision=decision, reason="protocol rationale", reviewer=reviewer)
        write_tsv(path, fields, rows)

    def test_dry_run_then_dual_import_preserves_baseline_and_flags_disagreement(self) -> None:
        self.set_decisions("reviewer_1", "Reviewer A", ["include", "exclude"])
        result = MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")
        self.assertEqual(result["imported_decisions"], 2)
        self.assertFalse(result["applied"])
        self.assertEqual(MODULE.sha256_file(self.queue), self.blank_sha)

        result = MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1", apply=True)
        self.assertEqual(result["imported_decisions"], 2)
        self.assertEqual(MODULE.sha256_file(self.queue.with_name("queue.before_reviewers_" + self.blank_sha[:12] + ".tsv")), self.blank_sha)

        self.set_decisions("reviewer_2", "Reviewer B", ["include", "unclear"])
        result = MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_2", apply=True)
        self.assertEqual(result["complete_independent_pairs"], 2)
        self.assertEqual(result["unresolved_disagreement_ids"], ["sid-2"])
        rows, _ = MODULE.read_tsv(self.queue)
        self.assertEqual(rows[1]["adjudication"], "")
        self.assertEqual(rows[1]["reviewer_1"], "Reviewer A")
        self.assertEqual(rows[1]["reviewer_2"], "Reviewer B")

    def test_apply_refuses_incomplete_packet(self) -> None:
        path = self.packet_dir / "reviewer_1/batch_0001.tsv"
        rows, fields = MODULE.read_tsv(path)
        rows[0].update(decision="include", reason="protocol rationale", reviewer="Reviewer A")
        write_tsv(path, fields, rows)
        with self.assertRaisesRegex(ValueError, "incomplete reviewer packet"):
            MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1", apply=True)
        self.assertEqual(MODULE.sha256_file(self.queue), self.blank_sha)

    def test_packet_source_mutation_is_rejected(self) -> None:
        path = self.packet_dir / "reviewer_1/batch_0001.tsv"
        rows, fields = MODULE.read_tsv(path)
        rows[0]["title"] = "Edited title"
        write_tsv(path, fields, rows)
        with self.assertRaisesRegex(ValueError, "source rows changed"):
            MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")

    def test_packet_manifest_tampering_is_rejected(self) -> None:
        path = self.packet_dir / "packet_manifest.json"
        manifest = json.loads(path.read_text())
        manifest["record_count"] = 3
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "manifest hash differs"):
            MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")

    def test_same_human_cannot_fill_both_reviewer_slots(self) -> None:
        self.set_decisions("reviewer_1", "Reviewer A", ["include", "exclude"])
        MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1", apply=True)
        self.set_decisions("reviewer_2", "Reviewer A", ["include", "exclude"])
        with self.assertRaisesRegex(ValueError, "same reviewer"):
            MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_2", apply=True)

    def test_reviewer_identity_must_be_consistent_within_slot(self) -> None:
        self.set_decisions("reviewer_1", "Reviewer A", ["include", "exclude"])
        path = self.packet_dir / "reviewer_1/batch_0001.tsv"
        rows, fields = MODULE.read_tsv(path)
        rows[1]["reviewer"] = "Reviewer C"
        write_tsv(path, fields, rows)
        with self.assertRaisesRegex(ValueError, "identity changes"):
            MODULE.import_packets(self.packet_dir, self.queue, self.protocol, "reviewer_1")


if __name__ == "__main__":
    unittest.main()
