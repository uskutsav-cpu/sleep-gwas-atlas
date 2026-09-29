import csv
import gzip
import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_placo_family_outputs as placo_audit
from audit_placo_family_outputs import audit_pair, verify_bh


def test_independent_bh_and_five_track_bonferroni_are_recomputed(tmp_path):
    path = tmp_path / "variants.tsv.gz"
    fields = ["SNP", "P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"]
    values = [("a", .01, .025), ("b", .01, .025), ("c", .04, 1/15), ("d", .8, 1.0)]
    with gzip.open(path, "wt", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for snp, p, q in values:
            writer.writerow({"SNP": snp, "P_PLACO": p, "Q_WITHIN_PAIR": q,
                             "Q_BONFERRONI_ACROSS_PAIRS": min(1, q*5), "status": "TESTED", "headline": "False"})
        for snp, state in (("excluded", "EXCLUDED_EXTREME_Z"), ("failed", "NUMERICAL_FAILURE")):
            writer.writerow({"SNP": snp, "P_PLACO": "NA", "Q_WITHIN_PAIR": "NA",
                             "Q_BONFERRONI_ACROSS_PAIRS": "NA", "status": state, "headline": "False"})
    con = sqlite3.connect(":memory:")
    con.execute("CREATE TABLE result (pair TEXT, p REAL, q REAL, across REAL, status TEXT)")
    audit_pair(path, "pair", 6, {"total_rows": 6, "tested": 4, "excluded_extreme_z": 1,
        "numerical_failures": 1, "failure_rate": .2, "denominator": 5,
        "headline_threshold": 1e-8, "headline_variants": 0}, 1e-8, 5, con)
    assert verify_bh(con, ["pair"], 5, {"pair": 5}) == 4


def test_audit_rejects_wrong_denominator(tmp_path):
    path = tmp_path / "variants.tsv.gz"
    fields = ["SNP", "P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"]
    with gzip.open(path, "wt", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({"SNP": "a", "P_PLACO": .1, "Q_WITHIN_PAIR": .1,
                         "Q_BONFERRONI_ACROSS_PAIRS": .5, "status": "TESTED", "headline": "False"})
    con = sqlite3.connect(":memory:")
    con.execute("CREATE TABLE result (pair TEXT, p REAL, q REAL, across REAL, status TEXT)")
    bad_status = {"total_rows": 2, "tested": 1, "excluded_extreme_z": 0, "numerical_failures": 0,
                  "failure_rate": 0, "denominator": 1, "headline_threshold": 1e-8, "headline_variants": 0}
    with pytest.raises(ValueError, match="denominator mismatch"):
        audit_pair(path, "pair", 2, bad_status, 1e-8, 5, con)


@pytest.mark.parametrize("column,value,error", [
    ("q", .016, "Within-pair BH mismatch"),
    ("across", .076, "Across-track Bonferroni mismatch"),
])
def test_vectorized_adjustment_audit_rejects_tampered_values(column, value, error):
    con = sqlite3.connect(":memory:")
    con.execute("CREATE TABLE result (pair TEXT, p REAL, q REAL, across REAL, status TEXT)")
    # Expected q values are .015 for the tied p values and .7 for the last row.
    rows = [("pair", .01, .015, .075, "TESTED"),
            ("pair", .01, .015, .075, "TESTED"),
            ("pair", .7, .7, 1., "TESTED")]
    rows[0] = tuple(value if field == column else old
                    for field, old in zip(("pair", "p", "q", "across", "status"), rows[0]))
    con.executemany("INSERT INTO result VALUES(?,?,?,?,?)", rows)
    with pytest.raises(ValueError, match=error):
        verify_bh(con, ["pair"], 5, {"pair": 3})


def test_published_and_collated_output_hashes_are_each_read_once(tmp_path, monkeypatch):
    published = tmp_path / "published.tsv.gz"
    collated = tmp_path / "collated.tsv.gz"
    published.write_bytes(b"same immutable payload")
    collated.write_bytes(b"same immutable payload")
    real_sha256 = placo_audit.sha256
    calls = []

    def count_hashes(path):
        calls.append(Path(path))
        return real_sha256(Path(path))

    monkeypatch.setattr(placo_audit, "sha256", count_hashes)
    expected = real_sha256(published)
    assert placo_audit.verify_published_output(published, collated, expected) == expected
    assert calls == [published, collated]

    collated.write_bytes(b"different payload")
    with pytest.raises(ValueError, match="Published table differs"):
        placo_audit.verify_published_output(published, collated, expected)
