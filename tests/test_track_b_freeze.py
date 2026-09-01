import csv
import hashlib
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRACK_B = ROOT / "results/track_b"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


class TrackBFreezeTests(unittest.TestCase):
    def test_repository_checkpoint_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/110_track_b_checkpoint.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_pair_lock_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/111_freeze_track_b_pairs.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_candidate_family_is_complete_12_by_6(self) -> None:
        panel = read_tsv(ROOT / "config/analysis_panel.tsv")
        sleep = {row["trait_id"] for row in panel if row["domain"] == "sleep"}
        brain = {
            row["trait_id"] for row in panel
            if row["domain"] in {"neuro", "psychiatric"}
        }
        rows = read_tsv(TRACK_B / "01_brain_pair_selection.tsv")
        self.assertEqual(len(rows), 72)
        self.assertEqual(
            {(row["sleep_trait"], row["brain_trait"]) for row in rows},
            {(a, b) for a in sleep for b in brain},
        )
        self.assertEqual(
            brain, {"alz", "parkinson", "mdd", "scz", "bipolar", "adhd"}
        )

    def test_non_discoveries_are_not_given_invented_literature_status(self) -> None:
        rows = read_tsv(TRACK_B / "01_brain_pair_selection.tsv")
        non_discoveries = [row for row in rows if float(row["FDR"]) >= 0.05]
        self.assertTrue(non_discoveries)
        self.assertTrue(all(
            row["novelty_status"] == "NOT_APPLICABLE_GLOBAL_FDR_NOT_SIGNIFICANT"
            and row["direct_prior_rg_status"] == "NOT_AUDITED_NONDISCOVERY"
            for row in non_discoveries
        ))

    def test_exact_three_roles_and_no_post_hoc_pair_b_replacement(self) -> None:
        rows = read_tsv(TRACK_B / "pair_manifest.tsv")
        self.assertEqual([row["pair_id"] for row in rows], ["A", "B", "CONTROL"])
        observed = {
            row["pair_id"]: (row["sleep_trait"], row["external_trait"])
            for row in rows
        }
        self.assertEqual(observed["A"], ("snoring", "parental_lifespan"))
        self.assertEqual(observed["B"], ("insomnia", "adhd"))
        self.assertEqual(observed["CONTROL"], ("insomnia", "frailty"))
        lock = json.loads((TRACK_B / "pair_manifest.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(
            lock["pair_b_identity_sha256"], sha256_text("B\tinsomnia\tadhd\n")
        )
        self.assertFalse(lock["downstream_results_accessed_before_pair_freeze"])
        self.assertEqual(lock["pair_replacement_policy"].split(";")[0], "FORBIDDEN")

    def test_all_three_pairs_are_dense_ready_but_not_declared_mechanistic(self) -> None:
        rows = read_tsv(TRACK_B / "pair_manifest.tsv")
        self.assertTrue(all(
            row["dense_data_readiness"] == "READY_FULL_SUMSTATS_BOTH" for row in rows
        ))
        forbidden = ("causal gene", "cell mechanism", "shared causal signal")
        for row in rows:
            joined = " ".join(row.values()).lower()
            self.assertTrue(all(term not in joined for term in forbidden))


if __name__ == "__main__":
    unittest.main()
