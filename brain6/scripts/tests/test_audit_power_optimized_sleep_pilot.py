import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "audit_power_optimized_sleep_pilot.py"
SPEC = importlib.util.spec_from_file_location("audit_power_optimized_sleep_pilot", SCRIPT)
AUDITOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(AUDITOR)


def _write_trait(root: Path, trait: str, tested: int) -> None:
    pair = root / "pairs" / trait
    (pair / "loci").mkdir(parents=True)
    (pair / "logs").mkdir()
    for index in range(1, 101):
        status = "TESTED" if index <= tested else "PHENOTYPE_DROPPED"
        receipt = {
            "trait": trait,
            "locus_index": index,
            "process_status": "PROCESSED",
            "sleep_status": status,
            "K": 10,
        }
        (pair / "loci" / f"locus_{index:04d}.json").write_text(json.dumps(receipt))
        (pair / "logs" / f"locus_{index:04d}.log").write_text("completed\n")


def _write_locus_file(path: Path) -> None:
    rows = ["LOC CHR START STOP"]
    rows.extend(f"{index} {(index - 1) % 22 + 1} 1 2" for index in range(1, 2496))
    path.write_text("\n".join(rows) + "\n")


def test_audit_computes_paired_eligibility_improvement(tmp_path):
    _write_trait(tmp_path, "longsleep", 60)
    _write_trait(tmp_path, "sleepdur", 70)
    locus_file = tmp_path / "loci.locfile"
    _write_locus_file(locus_file)

    result = AUDITOR.audit(tmp_path, locus_file)

    assert result["matched_locus_count"] == 100
    assert result["chromosomes_covered"] == 22
    assert result["eligibility_rate_difference_percentage_points"] == pytest.approx(10)
    assert result["relative_reduction_in_not_eligible"] == pytest.approx(0.25)
    assert result["traits"]["sleepdur"]["shared_reference_K"]["median"] == 10
    assert result["mutates_run_or_receipts"] is False


def test_audit_rejects_missing_receipt_log(tmp_path):
    _write_trait(tmp_path, "longsleep", 60)
    _write_trait(tmp_path, "sleepdur", 70)
    locus_file = tmp_path / "loci.locfile"
    _write_locus_file(locus_file)
    (tmp_path / "pairs" / "sleepdur" / "logs" / "locus_0001.log").unlink()

    with pytest.raises(ValueError, match="missing/empty paired log"):
        AUDITOR.audit(tmp_path, locus_file)
