#!/usr/bin/env python3
"""Independently audit published PLACO denominators and multiplicity values."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sqlite3
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCRATCH = Path("/private/tmp")
DEFAULT_OUTPUT_ROOT = Path("/Volumes/Extreme SSD/brain6-work/placo-family-v3-results")
FIELDS = {"SNP", "P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_published_output(published: Path, collated: Path, expected_sha256: str) -> str:
    """Hash each copy once, requiring both to match the frozen published digest."""
    published_sha256 = sha256(published)
    if published_sha256 != expected_sha256:
        raise ValueError(f"Published output hash mismatch: {published}")
    collated_sha256 = sha256(collated)
    if collated_sha256 != published_sha256:
        raise ValueError(f"Published table differs from QC-passed collation: {published}")
    return published_sha256


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-13)


def audit_pair(path: Path, pair: str, expected_rows: int, status: dict, threshold: float,
               n_pairs: int, con: sqlite3.Connection) -> dict:
    counts = {"TESTED": 0, "EXCLUDED_EXTREME_Z": 0, "NUMERICAL_FAILURE": 0}
    batch = []
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not FIELDS.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected PLACO schema in {path}: {reader.fieldnames}")
        for row in reader:
            state = row["status"]
            if state not in counts:
                raise ValueError(f"Unknown status {state!r} in {path}")
            counts[state] += 1
            raw = [row["P_PLACO"], row["Q_WITHIN_PAIR"], row["Q_BONFERRONI_ACROSS_PAIRS"]]
            if state == "TESTED":
                p, q, across = map(float, raw)
                if not all(math.isfinite(x) and 0 <= x <= 1 for x in (p, q, across)):
                    raise ValueError(f"Invalid P/q value in {path}: {row['SNP']}")
                if row["headline"] != str(p < threshold):
                    raise ValueError(f"Headline threshold mismatch in {path}: {row['SNP']}")
            else:
                if raw != ["NA", "NA", "NA"] or row["headline"] != "False":
                    raise ValueError(f"Non-tested row carries inferential values in {path}: {row['SNP']}")
                p = q = across = None
            batch.append((pair, p, q, across, state))
            if len(batch) == 100000:
                con.executemany("INSERT INTO result VALUES(?,?,?,?,?)", batch)
                batch.clear()
    if batch:
        con.executemany("INSERT INTO result VALUES(?,?,?,?,?)", batch)
    con.commit()
    total = sum(counts.values())
    denom = counts["TESTED"] + counts["NUMERICAL_FAILURE"]
    failure_rate = counts["NUMERICAL_FAILURE"] / denom if denom else math.inf
    expected_counts = {
        "total_rows": expected_rows, "tested": counts["TESTED"],
        "excluded_extreme_z": counts["EXCLUDED_EXTREME_Z"],
        "numerical_failures": counts["NUMERICAL_FAILURE"],
    }
    if total != expected_rows or denom != expected_rows - counts["EXCLUDED_EXTREME_Z"]:
        raise ValueError(f"PLACO family denominator mismatch for {pair}: {expected_counts}")
    for key, value in expected_counts.items():
        if int(status[key]) != value:
            raise ValueError(f"Pair status {key} mismatch for {pair}: {status[key]} != {value}")
    if int(status["denominator"]) != denom or not close(float(status["failure_rate"]), failure_rate):
        raise ValueError(f"Numerical failure rate mismatch for {pair}")
    if int(status["denominator"]) != denom or float(status["headline_threshold"]) != threshold:
        raise ValueError(f"Frozen denominator/threshold mismatch for {pair}")
    if int(status["headline_variants"]) != con.execute(
            "SELECT COUNT(*) FROM result WHERE pair=? AND status='TESTED' AND p<?", (pair, threshold)).fetchone()[0]:
        raise ValueError(f"Headline variant count mismatch for {pair}")
    return {"pair_id": pair, **expected_counts, "denominator": denom, "failure_rate": failure_rate,
            "headline_variants": int(status["headline_variants"]), "adjustment_tracks": n_pairs}


def verify_bh(con: sqlite3.Connection, pair_ids: list[str], n_pairs: int,
              denominators: dict[str, int]) -> int:
    """Recompute BH and five-track adjustment in SQLite without Python row iteration.

    Production outputs contain nearly 15 million tested cells. Keeping the rank,
    reverse cumulative minimum, and tolerance comparison inside SQLite avoids
    transferring each cell to Python while preserving the original rank and tie
    semantics.
    """
    con.execute("DROP TABLE IF EXISTS temp.bh_denominators")
    con.execute("CREATE TEMP TABLE bh_denominators (pair TEXT PRIMARY KEY, denominator INTEGER NOT NULL)")
    con.executemany("INSERT INTO bh_denominators VALUES(?,?)", sorted(denominators.items()))
    ctes = """
        WITH ranked AS (
            SELECT result.pair, result.p, result.q, result.across,
                   RANK() OVER(PARTITION BY result.pair ORDER BY result.p)
                   + COUNT(*) OVER(PARTITION BY result.pair, result.p) - 1 AS rank
            FROM result
            WHERE result.status='TESTED'
        ), scaled AS (
            SELECT ranked.pair, ranked.p, ranked.q, ranked.across,
                   CASE WHEN ranked.p * bh_denominators.denominator / ranked.rank < 1.0
                        THEN ranked.p * bh_denominators.denominator / ranked.rank
                        ELSE 1.0 END AS adjusted
            FROM ranked JOIN bh_denominators USING(pair)
        ), bh AS (
            SELECT pair, p, q, across,
                   MIN(adjusted) OVER(
                       PARTITION BY pair ORDER BY p DESC
                       ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
                   ) AS expected_q
            FROM scaled
        ), compared AS (
            SELECT pair, p,
                   ABS(q - expected_q) > MAX(1e-13, 1e-11 * MAX(ABS(q), ABS(expected_q))) AS bad_q,
                   ABS(across - MIN(1.0, expected_q * ?))
                       > MAX(1e-13, 1e-11 * MAX(ABS(across), ABS(MIN(1.0, expected_q * ?)))) AS bad_across
            FROM bh
        )
    """
    checked, pair_count, q_errors, across_errors, first_q_pair, first_q_p, first_across_pair, first_across_p = con.execute(
        ctes + """SELECT COUNT(*), COUNT(DISTINCT pair),
                    SUM(bad_q), SUM(bad_across),
                    MIN(CASE WHEN bad_q THEN pair END),
                    MIN(CASE WHEN bad_q THEN p END),
                    MIN(CASE WHEN bad_across THEN pair END),
                    MIN(CASE WHEN bad_across THEN p END)
                 FROM compared""",
        (n_pairs, n_pairs),
    ).fetchone()
    if pair_count != len(set(pair_ids)):
        observed_pairs = [row[0] for row in con.execute(
            "SELECT DISTINCT pair FROM result WHERE status='TESTED' ORDER BY pair")]
        raise ValueError(f"Expected tested pairs {sorted(pair_ids)}, observed {observed_pairs}")
    if q_errors:
        raise ValueError(f"Within-pair BH mismatch for {first_q_pair}/{first_q_p} ({q_errors} rows)")
    if across_errors:
        raise ValueError(f"Across-track Bonferroni mismatch for {first_across_pair}/{first_across_p} ({across_errors} rows)")
    return int(checked)


def run(scratch: Path = DEFAULT_SCRATCH, output_root: Path = DEFAULT_OUTPUT_ROOT) -> dict:
    family_path = ROOT / "brain6/config/placo_family_v3/family_lock.json"
    family = json.loads(family_path.read_text(encoding="utf-8"))
    required_tracks = 5
    if "5" not in family["across_pair_adjustment"]:
        raise ValueError("Frozen across-pair adjustment no longer specifies five tracks")
    master_path = ROOT / "brain6/results/placo/placo_master.tsv"
    with master_path.open(newline="", encoding="utf-8") as stream:
        master = {row["pair_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
    expected_pairs = set(family["pairs"])
    if set(master) != expected_pairs:
        raise ValueError(f"Published PLACO pair set differs from frozen v3 family: {sorted(master)}")
    reports = []
    output_hashes = {}
    source_receipt_hashes = {}
    with tempfile.TemporaryDirectory(prefix="brain6-placo-audit-", dir=scratch) as tmp:
        con = sqlite3.connect(Path(tmp) / "audit.sqlite")
        con.execute("PRAGMA journal_mode=OFF")
        con.execute("PRAGMA synchronous=OFF")
        con.execute("PRAGMA cache_size=-131072")
        con.execute("PRAGMA temp_store=FILE")
        con.execute("CREATE TABLE result (pair TEXT, p REAL, q REAL, across REAL, status TEXT)")
        for pair in sorted(expected_pairs):
            record = master[pair]
            if not record["stage_status"].startswith("COMPLETE_QC_PASS"):
                raise ValueError(f"Pair is not QC-passed: {pair}")
            out = ROOT / record["output_path"]
            receipt_path = output_root / pair / "collated"
            receipt = json.loads((receipt_path / "receipt.json").read_text(encoding="utf-8"))
            status = json.loads((receipt_path / "status.json").read_text(encoding="utf-8"))
            if receipt.get("scientific_status") != "PASS" or status.get("status") != "PASS":
                raise ValueError(f"Collation QC did not pass: {pair}")
            params = receipt["parameters"]
            if int(params["n_pairs"]) != required_tracks or int(params["expected_rows"]) != int(family["pairs"][pair]["expected_rows"]):
                raise ValueError(f"Frozen PLACO denominator differs for {pair}")
            if float(params["max_failure_rate"]) != float(family["max_pair_numerical_failure_rate"]):
                raise ValueError(f"Frozen failure ceiling differs for {pair}")
            output_hashes[pair] = verify_published_output(
                out, receipt_path / "results.tsv.gz", record["output_sha256"]
            )
            source_receipt_hashes[pair] = sha256(receipt_path / "receipt.json")
            reports.append(audit_pair(out, pair, int(params["expected_rows"]), status,
                                      float(family["headline_threshold"]), required_tracks, con))
            if reports[-1]["failure_rate"] > float(params["max_failure_rate"]):
                raise ValueError(f"Numerical failure ceiling exceeded: {pair}")
        denominators = {row["pair_id"]: row["denominator"] for row in reports}
        q_values_checked = verify_bh(con, sorted(expected_pairs), required_tracks, denominators)
        con.close()
    result = {
        "audited_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "PASS_FOUR_PAIR_QC; FIVE_TRACK_FAMILY_INCOMPLETE",
        "family_lock_sha256": sha256(family_path), "n_required_tracks": required_tracks,
        "published_v3_pairs": sorted(expected_pairs), "protected_track_b_published": "insomnia__adhd" in master,
        "published_output_sha256": output_hashes, "collation_receipt_sha256": source_receipt_hashes,
        "full_family_complete": False, "q_values_independently_verified": q_values_checked,
        "pairs": reports,
        "interpretation": "Four PLACO v3 pairs pass pair-level QC; five-track correction is verified, but the protected Track B output is not present and no full-family conclusion is made.",
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "brain6/results/placo/placo_v3_pair_qc_validation.json")
    args = parser.parse_args()
    report = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{args.output.name}.", dir=args.output.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    finally:
        temporary.unlink(missing_ok=True)
    print(json.dumps(report, indent=2))
