import csv
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRACK_B = ROOT / "results/track_b"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBLiteratureTests(unittest.TestCase):
    def test_audit_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/112_build_track_b_literature_audit.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_evidence_map_has_frozen_pairs_and_allowed_classes(self) -> None:
        rows = read_tsv(TRACK_B / "02_prior_evidence_map.tsv")
        self.assertGreaterEqual(len(rows), 20)
        observed = {(r["pair_id"], r["sleep_trait"], r["external_trait"]) for r in rows}
        self.assertEqual(observed, {
            ("A", "snoring", "parental_lifespan"),
            ("B", "insomnia", "adhd"),
            ("CONTROL", "insomnia", "frailty"),
        })
        allowed = {
            "DIRECT_PRIOR_GLOBAL", "DIRECT_PRIOR_LOCAL", "RELATED_ONLY",
            "PRIOR_SHARED_LOCUS", "PRIOR_CELL_MECHANISM", "UNDERREPORTED", "UNCERTAIN",
        }
        self.assertTrue({r["classification"] for r in rows}.issubset(allowed))
        self.assertEqual({r["classification"] for r in rows}, allowed)

    def test_all_requested_method_families_were_searched_for_each_pair(self) -> None:
        rows = read_tsv(TRACK_B / "02_literature_search_queries.tsv")
        by_pair: dict[str, set[str]] = {}
        for row in rows:
            by_pair.setdefault(row["pair_id"], set()).add(row["method_family"])
            self.assertIn("NOT_PROOF_OF_ABSENCE", row["interpretation"].upper().replace(" ", "_")) if row["audit_result"].startswith("NO_") else None
        self.assertEqual(set(by_pair), {"A", "B", "CONTROL"})
        self.assertTrue(all(len(methods) == 18 for methods in by_pair.values()))
        self.assertEqual(len(rows), 54)

    def test_primary_source_rows_have_identifiers_and_limits(self) -> None:
        rows = read_tsv(TRACK_B / "02_prior_evidence_map.tsv")
        primary = [r for r in rows if r["PMID"] != "NA"]
        self.assertTrue(primary)
        for row in primary:
            self.assertTrue(row["source_url"].startswith("https://"))
            self.assertNotEqual(row["limitations"], "NA")
        self.assertTrue(any("OVERLAP" in r["cohort_relation_to_atlas"] for r in rows))

    def test_pair_b_prior_art_blocks_generic_novelty_claim(self) -> None:
        rows = read_tsv(TRACK_B / "02_prior_evidence_map.tsv")
        pair_b = [r for r in rows if r["pair_id"] == "B"]
        self.assertTrue(any(r["method"].startswith("PLACO") for r in pair_b))
        self.assertTrue(any("PP4" in r["shared_locus_evidence"] for r in pair_b))
        self.assertTrue(any(r["classification"] == "PRIOR_CELL_MECHANISM" for r in pair_b))


if __name__ == "__main__":
    unittest.main()
