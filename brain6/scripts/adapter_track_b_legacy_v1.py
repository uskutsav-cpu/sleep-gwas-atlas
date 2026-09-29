#!/usr/bin/env python3
"""Stage and independently validate the authorized, exact historical B crosswalk.

This program never writes to the repository's protected Track B slot.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/track_b_admission_v1"
STAGE = Path("/Volumes/Extreme SSD/brain6-work/track-b-legacy-admission-v1")
SOURCE = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/results/track_b/pleiotropy/results/placo/B.full.tsv.gz")
FIELDS17 = ["analysis_id", "pair_id", "SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2", "T_PLACO_PLUS", "P_PLACO_PLUS", "PLACO_BH_Q", "within_pair_family_n", "analysis_status", "numerical_error"]
FIELDS10 = ["SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite_probability(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f"invalid probability: {value}")
    return number


def run() -> dict:
    criteria_path = OUT / "criteria.lock.json"
    criteria = json.loads(criteria_path.read_text())
    expected = criteria["expected"]
    preflight_path = OUT / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    validator_path = ROOT / "brain6/scripts/validate_track_b_legacy_admission_v1.py"
    if preflight["status"] != "PASS" or preflight["failed_checks"]:
        raise ValueError("admission preflight failed")
    if preflight["criteria_sha256"] != sha(criteria_path) or preflight["validator_sha256"] != sha(validator_path):
        raise ValueError("preflight source identity changed")
    if sha(SOURCE) != expected["ledger_sha256"]:
        raise ValueError("historical B source changed since preflight")
    stage_view = STAGE / "B.legacy_v3_adapter.v1.attempt2.tsv.gz"
    db_path = STAGE / "B.legacy_v3_adapter.v1.attempt2.audit.sqlite"
    if stage_view.exists() or db_path.exists():
        raise FileExistsError("adapter staging path already exists; refusing to overwrite")
    STAGE.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path)
    # This external volume permits ordinary writes but not SQLite journal sidecars.
    # The database is disposable audit scratch; the source and receipts stay immutable.
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA temp_store=FILE")
    con.execute("PRAGMA cache_size=-131072")
    con.execute("CREATE TABLE result (snp TEXT PRIMARY KEY, pair TEXT NOT NULL, p REAL NOT NULL, q REAL NOT NULL, across REAL NOT NULL, status TEXT NOT NULL)")
    rows = headlines = q_significant = 0
    last_coordinate = (0, 0, "")
    batch = []
    with gzip.open(SOURCE, "rt", newline="", encoding="utf-8") as source, stage_view.open("xb") as binary:
        with gzip.GzipFile(fileobj=binary, mode="wb", mtime=0, filename="") as compressed:
            import io
            text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
            reader = csv.DictReader(source, delimiter="\t", strict=True)
            if reader.fieldnames != FIELDS17:
                raise ValueError(f"unexpected historical schema: {reader.fieldnames}")
            writer = csv.writer(text, delimiter="\t", lineterminator="\n")
            writer.writerow(FIELDS10)
            for row in reader:
                rows += 1
                if None in row or len(row) != len(FIELDS17):
                    raise ValueError(f"historical field count differs at row {rows}")
                if row["analysis_id"] != "track-b-v1.0-pleiotropy" or row["pair_id"] != "B":
                    raise ValueError(f"cross-run row at {rows}")
                if row["analysis_status"] != "TESTED" or row["numerical_error"] != "NA":
                    raise ValueError(f"unexpected terminal status at {rows}")
                if row["within_pair_family_n"] != str(expected["rows"]):
                    raise ValueError(f"family denominator changed at {rows}")
                if not row["SNP"] or any(c in row["SNP"] for c in "\t\n\r"):
                    raise ValueError(f"invalid SNP at {rows}")
                chromosome, position = int(row["CHR"]), int(row["BP"])
                if not 1 <= chromosome <= 22 or position <= 0:
                    raise ValueError(f"invalid genomic coordinate at {rows}")
                coordinate = (chromosome, position, row["SNP"])
                if coordinate <= last_coordinate:
                    raise ValueError(f"unordered or repeated coordinate at {rows}")
                last_coordinate = coordinate
                if row["A1"] not in {"A", "C", "G", "T"} or row["A2"] not in {"A", "C", "G", "T"} or row["A1"] == row["A2"]:
                    raise ValueError(f"invalid allele at {rows}")
                for field in ("Z1", "Z2", "T_PLACO_PLUS"):
                    if not math.isfinite(float(row[field])):
                        raise ValueError(f"nonfinite {field} at {rows}")
                finite_probability(row["P1"])
                finite_probability(row["P2"])
                p = finite_probability(row["P_PLACO_PLUS"])
                q = finite_probability(row["PLACO_BH_Q"])
                if q + 1e-13 < p:
                    raise ValueError(f"q below raw P at {rows}")
                across = min(1.0, 5.0 * q)
                headline = p < criteria["frozen_statistics"]["headline_p_strictly_less_than"]
                headlines += headline
                q_significant += q < criteria["frozen_statistics"]["within_pair_bh_alpha"]
                writer.writerow((row["SNP"], row["CHR"], row["BP"], row["Z1"], row["Z2"],
                                 row["P_PLACO_PLUS"], row["PLACO_BH_Q"],
                                 format(across, ".17g"), "TESTED", str(headline)))
                batch.append((row["SNP"], "insomnia__adhd", p, q, across, "TESTED"))
                if len(batch) >= 20000:
                    con.executemany("INSERT INTO result VALUES(?,?,?,?,?,?)", batch)
                    con.commit()
                    batch.clear()
            if batch:
                con.executemany("INSERT INTO result VALUES(?,?,?,?,?,?)", batch)
                con.commit()
            text.flush()
            text.detach()
    if rows != expected["rows"]:
        raise ValueError(f"row count {rows} differs from {expected['rows']}")
    if con.execute("SELECT COUNT(*) FROM result").fetchone()[0] != rows:
        raise ValueError("SQLite uniqueness/row count mismatch")
    # Recompute BH from P alone, using upper ranks for ties and a reverse cumulative minimum.
    query = """
        WITH ranked AS (
            SELECT p, q, across,
                   RANK() OVER (ORDER BY p) + COUNT(*) OVER (PARTITION BY p) - 1 AS upper_rank
            FROM result
        ), scaled AS (
            SELECT p,q,across,MIN(1.0,p * ? / upper_rank) AS adjusted FROM ranked
        ), bh AS (
            SELECT p,q,across,
                   MIN(adjusted) OVER (ORDER BY p DESC ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS expected_q
            FROM scaled
        )
        SELECT COUNT(*),
               SUM(ABS(q-expected_q) > MAX(1e-13,1e-11*MAX(ABS(q),ABS(expected_q)))),
               SUM(ABS(across-MIN(1.0,5.0*expected_q)) > MAX(1e-13,1e-11*MAX(ABS(across),ABS(MIN(1.0,5.0*expected_q)))))
        FROM bh
    """
    checked, q_bad, across_bad = con.execute(query, (rows,)).fetchone()
    con.close()
    if (checked, q_bad, across_bad) != (rows, 0, 0):
        raise ValueError(f"independent BH failed: {(checked,q_bad,across_bad)}")
    # A second pass catches serialization/field mapping errors in the staged view.
    adapted_rows = adapted_headlines = 0
    with gzip.open(stage_view, "rt", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames != FIELDS10:
            raise ValueError("adapted schema differs")
        for row in reader:
            adapted_rows += 1
            if None in row or len(row) != len(FIELDS10) or row["status"] != "TESTED":
                raise ValueError(f"adapted row invalid at {adapted_rows}")
            p,q,across = (finite_probability(row[k]) for k in ("P_PLACO","Q_WITHIN_PAIR","Q_BONFERRONI_ACROSS_PAIRS"))
            if not math.isclose(across,min(1.0,5*q),rel_tol=1e-14,abs_tol=1e-15):
                raise ValueError(f"adapted Bonferroni mismatch at {adapted_rows}")
            if row["headline"] != str(p < criteria["frozen_statistics"]["headline_p_strictly_less_than"]):
                raise ValueError(f"adapted headline mismatch at {adapted_rows}")
            adapted_headlines += row["headline"] == "True"
    if (adapted_rows,adapted_headlines) != (rows,headlines):
        raise ValueError("adapted row/headline count differs")
    return {
        "status":"PASS", "criteria_sha256":sha(criteria_path), "preflight_sha256":sha(preflight_path),
        "adapter_sha256":sha(Path(__file__)), "source_ledger_sha256":sha(SOURCE),
        "staged_view_path":str(stage_view), "staged_view_sha256":sha(stage_view),
        "source_fields":FIELDS17,"adapted_fields":FIELDS10,"pair_id":"insomnia__adhd",
        "rows":rows,"unique_snps":rows,"independent_bh_checked_rows":checked,
        "independent_bh_q_mismatches":q_bad,"five_track_adjustment_mismatches":across_bad,
        "headline_variants":headlines,"within_pair_q_below_0_05":q_significant,
        "historical_abs_tol":criteria["frozen_statistics"]["historical_absolute_tolerance"],
        "v3_abs_tol":criteria["frozen_statistics"]["v3_absolute_tolerance"],
        "historical_method_parameter_policy":criteria["frozen_statistics"]["historical_abs_tol_difference_policy"],
    }


if __name__ == "__main__":
    result = run()
    receipt = OUT / "adapter_validation.json"
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if receipt.exists() and receipt.read_text() != payload:
        raise FileExistsError("refusing to overwrite an existing adapter receipt")
    if not receipt.exists():
        receipt.write_text(payload)
    print(json.dumps({k:result[k] for k in ("status","rows","headline_variants","within_pair_q_below_0_05","staged_view_sha256")}))
