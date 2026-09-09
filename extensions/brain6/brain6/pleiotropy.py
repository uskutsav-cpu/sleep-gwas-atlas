"""Genome-wide PLACO chunk preparation and disk-backed full-family adjustment."""
from __future__ import annotations
import csv
import math
import sqlite3
from pathlib import Path
from .artifacts import transaction, verify_artifact
from .gwas import PAIR_FIELDS
from .io import (atomic_text, file_record, read_json, read_tsv, require, write_json, write_tsv)


def split_pair(pair_dir, root, relative, *, chunk_rows=20000):
    require(chunk_rows > 0, "Invalid chunk size")
    pair_dir = Path(pair_dir)
    receipt = verify_artifact(pair_dir)
    require(receipt.get("analysis_scope") != "GWAS_QTL_LOCUS", "Cis-only QTL data cannot estimate genome-wide PLACO parameters")
    with transaction(root, relative, stage="split_placo", inputs=[pair_dir / "pair.tsv.gz"],
                     parameters={"chunk_rows": chunk_rows}, synthetic=receipt["synthetic"]) as (work, meta):
        chunks = []; current = []
        def flush():
            name = f"chunk_{len(chunks):06d}.tsv"
            write_tsv(work / name, PAIR_FIELDS, current)
            chunks.append({"path": name, "rows": len(current)})
            current.clear()
        for row in read_tsv(pair_dir / "pair.tsv.gz", PAIR_FIELDS):
            current.append(row)
            if len(current) == chunk_rows:
                flush()
        if current:
            flush()
        require(chunks, "No chunks")
        write_json(work / "chunks.json", {"chunks": chunks, "total_rows": sum(x["rows"] for x in chunks)})
    return Path(root) / relative


def collate_placo(manifest, root, relative, *, expected_rows, n_pairs, max_failure_rate=.001,
                  synthetic=False):
    """Global per-pair BH over all eligible variants, not chunk-wise BH.

    Numerical failures are assigned p=1 ONLY for conservative denominator
    accounting. Their reported p/q remain NA and their status remains failed.
    """
    require(expected_rows > 0 and n_pairs >= 1, "Missing full denominator")
    entries = list(read_tsv(manifest, ["path", "sha256"]))
    from .io import check_hash
    for e in entries:
        check_hash(e["path"], e["sha256"])
    inputs = [manifest] + [e["path"] for e in entries]
    params = dict(expected_rows=expected_rows, n_pairs=n_pairs, max_failure_rate=max_failure_rate)
    with transaction(root, relative, stage="collate_placo", inputs=inputs, parameters=params,
                     synthetic=synthetic) as (work, meta):
        con = sqlite3.connect(work / "statistics.sqlite")
        con.execute("PRAGMA cache_size=-65536")
        con.execute("PRAGMA temp_store=FILE")
        con.execute("CREATE TABLE results (snp TEXT PRIMARY KEY, chr INTEGER,bp INTEGER,z1 REAL,z2 REAL,p REAL,status TEXT,q REAL)")
        batch = []; count = failures = excluded = 0
        for e in entries:
            for r in read_tsv(e["path"], ["SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO", "status"]):
                count += 1
                status = r["status"]
                require(status in {"TESTED", "NUMERICAL_FAILURE", "EXCLUDED_EXTREME_Z"}, "Unknown PLACO row status")
                p = float(r["P_PLACO"]) if status == "TESTED" else None
                if p is not None:
                    require(math.isfinite(p) and 0 <= p <= 1, "Invalid PLACO P")
                failures += status == "NUMERICAL_FAILURE"
                excluded += status == "EXCLUDED_EXTREME_Z"
                batch.append((r["SNP"], int(r["CHR"]), int(r["BP"]), float(r["Z1"]), float(r["Z2"]), p, status))
                if len(batch) >= 20000:
                    con.executemany("INSERT INTO results(snp,chr,bp,z1,z2,p,status) VALUES(?,?,?,?,?,?,?)", batch)
                    con.commit(); batch.clear()
        con.executemany("INSERT INTO results(snp,chr,bp,z1,z2,p,status) VALUES(?,?,?,?,?,?,?)", batch)
        con.commit()
        require(count == expected_rows, f"Missing/extra PLACO rows: expected {expected_rows}, got {count}")
        denom = count - excluded
        require(denom > 0, "No eligible PLACO tests")
        failure_rate = failures / denom
        con.execute("CREATE INDEX p_order ON results(p)")
        # Window rank handles ties using upper rank (same BH value for tied P).
        con.execute("""CREATE TEMP TABLE ranks AS
          SELECT snp,p,RANK() OVER(ORDER BY p)+COUNT(*) OVER(PARTITION BY p)-1 AS rank
          FROM results WHERE status='TESTED'""")
        q = 1.; batch = []
        for snp, p, rank in con.execute("SELECT snp,p,rank FROM ranks ORDER BY p DESC,snp DESC"):
            q = min(q, p * denom / rank)
            batch.append((q, snp))
            if len(batch) >= 20000:
                con.executemany("UPDATE results SET q=? WHERE snp=?", batch); batch.clear()
        con.executemany("UPDATE results SET q=? WHERE snp=?", batch); con.commit()
        threshold = 5e-8 / n_pairs
        fields = ["SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"]
        def export():
            for snp,ch,bp,z1,z2,p,status,q in con.execute("SELECT snp,chr,bp,z1,z2,p,status,q FROM results ORDER BY chr,bp,snp"):
                yield dict(zip(fields, [snp,ch,bp,z1,z2,"NA" if p is None else p,
                           "NA" if q is None else q,"NA" if q is None else min(1,q*n_pairs),
                           status, bool(p is not None and p < threshold)]))
        write_tsv(work / "results.tsv.gz", fields, export())
        hits = con.execute("SELECT COUNT(*) FROM results WHERE p<?", (threshold,)).fetchone()[0]
        status = "PASS" if failure_rate <= max_failure_rate else "FAILED_QC_NOT_CONSUMED"
        write_json(work / "status.json", dict(status=status, total_rows=count, tested=denom-failures,
                   numerical_failures=failures, excluded_extreme_z=excluded, denominator=denom,
                   failure_rate=failure_rate, headline_threshold=threshold, headline_variants=hits,
                   independent_loci="NOT_YET_ESTIMATED"))
        meta["scientific_status"] = status
        con.close()
    return Path(root) / relative
