from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/25_collate_fi_sleep_rg.py"
SPEC = importlib.util.spec_from_file_location("collate_fi_sleep_rg", SCRIPT)
COLLATOR = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(COLLATOR)


class FISleepRgCollatorTests(unittest.TestCase):
    def build_inputs(self, root: Path, changed_frozen_p: bool = False):
        pair_dir = root / "pairs"
        pair_dir.mkdir()
        frozen_rows = ["sleep_trait\tdisease_trait\tglobal_rg_p\tglobal_rg_fdr_all_396"]
        overlap_rows = ["frailty_endpoint\tsleep_trait\tcohort_overlap_status\texact_participant_overlap"]
        for index, trait in enumerate(sorted(COLLATOR.EXPECTED_SLEEP), start=1):
            p = 1e-8 * index
            frozen_p = p * (2 if changed_frozen_p and index == 1 else 1)
            (pair_dir / f"{trait}.tsv").write_text(
                "sleep_trait\tdisease_trait\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se\tinput_log\tfdr\n"
                f"{trait}\tfrailty\t0.1\t0.02\t5\t{p}\t0.1\t0.01\t1.02\t0.01\t0.01\t0.002\tlog\t{p}\n",
                encoding="utf-8",
            )
            frozen_rows.append(f"{trait}\tfrailty\t{frozen_p}\t0.01")
            overlap_rows.append(f"frailty\t{trait}\tCOHORT_OVERLAP_EXPECTED\tUNKNOWN")
        frozen = root / "frozen.tsv"
        frozen.write_text("\n".join(frozen_rows) + "\n", encoding="utf-8")
        overlap = root / "overlap.tsv"
        overlap.write_text("\n".join(overlap_rows) + "\n", encoding="utf-8")
        return pair_dir, frozen, overlap

    def test_uses_locked_full_family_q_after_exact_raw_p_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            pair_dir, frozen, overlap = self.build_inputs(Path(temp))
            result = COLLATOR.collate(pair_dir, frozen, overlap)
        self.assertEqual(len(result), 12)
        self.assertEqual(set(result["fdr_family_denominator"]), {396})
        self.assertTrue(result["global_rg_fdr_all_396"].notna().all())
        self.assertTrue(result["fdr_source"].str.contains("verified equal").all())

    def test_refuses_stale_q_values_when_raw_p_differs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            pair_dir, frozen, overlap = self.build_inputs(Path(temp), changed_frozen_p=True)
            with self.assertRaisesRegex(ValueError, "raw p-values differ"):
                COLLATOR.collate(pair_dir, frozen, overlap)
