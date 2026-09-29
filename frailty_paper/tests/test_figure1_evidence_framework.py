from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/37_plot_evidence_framework.py"
SPEC = importlib.util.spec_from_file_location("frailty_figure1", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def review_rows(*, pubmed: int = 61009, manual: int = 0, duplicates: int = 4892,
                queue: int = 56117, screened: int = 0, full_texts: int = 0,
                included: int = 0):
    return [
        {"stage": "Records identified from PubMed searches", "count": str(pubmed)},
        {"stage": "Records imported from other databases", "count": str(manual)},
        {"stage": "Duplicate records removed across all sources", "count": str(duplicates)},
        {"stage": "Records after deduplication", "count": str(queue)},
        {"stage": "Records screened", "count": str(screened)},
        {"stage": "Records excluded at title/abstract", "count": "0"},
        {"stage": "Reports sought for retrieval", "count": "0"},
        {"stage": "Reports not retrieved", "count": "0"},
        {"stage": "Reports assessed for eligibility", "count": str(full_texts)},
        {"stage": "Reports excluded after full text", "count": "0"},
        {"stage": "Studies included", "count": str(included)},
    ]


def source_rows(*, pubmed: int = 61009, queue: int = 56117):
    return [{"database": "PubMed", "records_imported_or_retrieved": str(pubmed),
             "unique_records_retained_by_priority": str(queue)}]


class Figure1ReviewCountTests(unittest.TestCase):
    def test_current_pubmed_only_counts_reconcile(self) -> None:
        counts = MODULE.reconcile_review_counts(review_rows(), source_rows())
        self.assertEqual(counts, {"pubmed": 61009, "manual": 0, "duplicates": 4892, "queue": 56117})

    def test_source_ledger_mismatch_is_rejected(self) -> None:
        with self.assertRaisesRegex(SystemExit, "disagree"):
            MODULE.reconcile_review_counts(review_rows(), source_rows(pubmed=60989, queue=56092))

    def test_nonconserving_duplicate_count_is_rejected(self) -> None:
        with self.assertRaisesRegex(SystemExit, "do not reconcile"):
            MODULE.reconcile_review_counts(review_rows(duplicates=4897), source_rows())

    def test_manual_records_require_revising_pubmed_only_copy(self) -> None:
        with self.assertRaisesRegex(SystemExit, "PubMed-only"):
            MODULE.reconcile_review_counts(review_rows(manual=1), source_rows())

    def test_screening_progress_requires_revising_figure_copy(self) -> None:
        with self.assertRaisesRegex(SystemExit, "Review status changed"):
            MODULE.reconcile_review_counts(review_rows(screened=1), source_rows())

    def test_full_text_progress_requires_revising_figure_copy(self) -> None:
        with self.assertRaisesRegex(SystemExit, "Review status changed"):
            MODULE.reconcile_review_counts(review_rows(full_texts=1), source_rows())

    def test_inclusion_progress_requires_revising_figure_copy(self) -> None:
        with self.assertRaisesRegex(SystemExit, "Review status changed"):
            MODULE.reconcile_review_counts(review_rows(included=1), source_rows())


if __name__ == "__main__":
    unittest.main()
