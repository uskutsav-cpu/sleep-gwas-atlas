import csv
import gzip
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from publish_placo_results import FIELDS, atomic_copy, validate_variant_table


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(snp: str, p: str, q: str, across: str, status: str = "TESTED",
        headline: str = "False") -> dict[str, str]:
    return {"SNP": snp, "CHR": "1", "BP": "100", "Z1": "1.0", "Z2": "2.0",
            "P_PLACO": p, "Q_WITHIN_PAIR": q, "Q_BONFERRONI_ACROSS_PAIRS": across,
            "status": status, "headline": headline}


def test_validates_counts_and_frozen_headline_threshold(tmp_path):
    path = tmp_path / "variants.tsv.gz"
    write_table(path, [row("rs1", "1e-9", "1e-8", "5e-8", headline="True"),
                       row("rs2", "0.2", "0.3", "1.0")])
    assert validate_variant_table(path, 2, 1e-8) == {
        "total_rows": 2, "tested_rows": 2, "excluded_extreme_z": 0,
        "numerical_failures": 0, "headline_variants": 1,
    }


def test_rejects_headline_flag_that_disagrees_with_p(tmp_path):
    path = tmp_path / "variants.tsv.gz"
    write_table(path, [row("rs1", "1e-9", "1e-8", "5e-8")])
    with pytest.raises(ValueError, match="Headline flag disagrees"):
        validate_variant_table(path, 1, 1e-8)


def test_rejects_duplicate_snp_ids(tmp_path):
    path = tmp_path / "variants.tsv.gz"
    write_table(path, [row("rs1", "0.2", "0.3", "1.0"), row("rs1", "0.3", "0.4", "1.0")])
    with pytest.raises(ValueError, match="duplicate SNP"):
        validate_variant_table(path, 2, 1e-8)


def test_rejects_nonfinite_or_out_of_range_p_q(tmp_path):
    path = tmp_path / "variants.tsv.gz"
    write_table(path, [row("rs1", "nan", "0.3", "1.0")])
    with pytest.raises(ValueError, match="Invalid P/q"):
        validate_variant_table(path, 1, 1e-8)


def test_atomic_copy_is_idempotent_and_refuses_to_clobber(tmp_path):
    source = tmp_path / "source"
    target = tmp_path / "published" / "result.tsv.gz"
    source.write_bytes(b"complete-result")
    atomic_copy(source, target)
    atomic_copy(source, target)
    assert target.read_bytes() == source.read_bytes()
    source.write_bytes(b"different-result")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        atomic_copy(source, target)
