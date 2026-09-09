import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRACK_B = ROOT / "results/track_b"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TrackBDenseQCTests(unittest.TestCase):
    def test_dense_qc_verifies(self) -> None:
        result = subprocess.run(
            ["python3", "scripts/115_build_track_b_dense_qc.py", "--verify"],
            cwd=ROOT, check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_exact_five_trait_family_and_schema(self) -> None:
        rows = read_tsv(TRACK_B / "03_dense_input_qc.tsv")
        self.assertEqual({row["trait_id"] for row in rows}, {
            "snoring", "parental_lifespan", "insomnia", "adhd", "frailty",
        })
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(
            row["observed_schema"] == "SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N"
            for row in rows
        ))
        self.assertTrue(all(int(row["rows_out"]) > 5_000_000 for row in rows))

    def test_all_dense_but_info_limitations_are_not_hidden(self) -> None:
        rows = read_tsv(TRACK_B / "03_dense_input_qc.tsv")
        by_trait = {row["trait_id"]: row for row in rows}
        self.assertTrue(by_trait["snoring"]["INFO_status"].startswith("ABSENT_IN_RELEASE"))
        self.assertTrue(by_trait["frailty"]["INFO_status"].startswith("ABSENT_IN_RELEASE"))
        self.assertTrue(by_trait["insomnia"]["INFO_status"].startswith("ROW_WISE_FILTERED"))
        self.assertEqual(
            set(json.loads((TRACK_B / "03_dense_input_qc.lock.json").read_text(encoding="utf-8"))["rowwise_info_limitations"]),
            {"snoring", "frailty"},
        )

    def test_fine_mapping_gate_does_not_claim_mechanism(self) -> None:
        rows = read_tsv(TRACK_B / "03_dense_input_qc.tsv")
        forbidden = ("causal gene", "shared causal signal", "cell mechanism")
        for row in rows:
            joined = " ".join(row.values()).lower()
            self.assertTrue(all(value not in joined for value in forbidden))


if __name__ == "__main__":
    unittest.main()
