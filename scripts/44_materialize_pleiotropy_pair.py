#!/usr/bin/env python3
"""Build one checksum-locked, allele-aligned full-resolution trait-pair file."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import sqlite3
from pathlib import Path


REQUIRED = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"}
OUTPUT_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2", "N1", "N2"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rows_out(qc: Path) -> int:
    for line in qc.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0] == "rows_out":
            return int(fields[1])
    raise SystemExit(f"ERROR: {qc} lacks rows_out")


def variant_rows(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"ERROR: {path} lacks fields: {sorted(missing)}")
        for row in reader:
            try:
                chromosome = int(row["CHR"])
                position = int(row["BP"])
                beta = float(row["BETA"])
                se = float(row["SE"])
                p_value = float(row["P"])
                sample_size = float(row["N"])
            except ValueError as exc:
                raise SystemExit(f"ERROR: invalid numeric value in {path}") from exc
            if not (
                row["SNP"].startswith("rs")
                and 1 <= chromosome <= 22
                and position > 0
                and len(row["A1"]) == len(row["A2"]) == 1
                and row["A1"].upper() in "ACGT"
                and row["A2"].upper() in "ACGT"
                and row["A1"].upper() != row["A2"].upper()
                and math.isfinite(beta)
                and math.isfinite(se)
                and se > 0
                and math.isfinite(p_value)
                and 0 < p_value <= 1
                and math.isfinite(sample_size)
                and sample_size > 0
            ):
                raise SystemExit(f"ERROR: invalid variant row in {path}: {row.get('SNP', '')}")
            yield (
                row["SNP"].lower(), chromosome, position, row["A1"].upper(),
                row["A2"].upper(), beta / se, p_value, sample_size,
            )


def load_table(connection: sqlite3.Connection, name: str, path: Path) -> int:
    connection.execute(
        f"CREATE TABLE {name} ("
        "snp TEXT PRIMARY KEY, chr INTEGER NOT NULL, bp INTEGER NOT NULL, "
        "a1 TEXT NOT NULL, a2 TEXT NOT NULL, z REAL NOT NULL, p REAL NOT NULL, n REAL NOT NULL"
        ") WITHOUT ROWID"
    )
    statement = f"INSERT INTO {name} VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
    count = 0
    batch = []
    try:
        for values in variant_rows(path):
            batch.append(values)
            if len(batch) == 50000:
                connection.executemany(statement, batch)
                count += len(batch)
                batch.clear()
        if batch:
            connection.executemany(statement, batch)
            count += len(batch)
        connection.commit()
    except sqlite3.IntegrityError as exc:
        raise SystemExit(f"ERROR: duplicate SNP in {path}") from exc
    return count


def align(a1: str, a2: str, b1: str, b2: str, z: float) -> tuple[float, str] | None:
    if (b1, b2) == (a1, a2):
        return z, "DIRECT"
    if (b1, b2) == (a2, a1):
        return -z, "SWAPPED"
    c1, c2 = b1.translate(COMPLEMENT), b2.translate(COMPLEMENT)
    if (c1, c2) == (a1, a2):
        return z, "COMPLEMENT"
    if (c1, c2) == (a2, a1):
        return -z, "COMPLEMENT_SWAPPED"
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/pleiotropy_pair_manifest.tsv")
    parser.add_argument("--lock", default="results/tables/pleiotropy_pair_manifest.lock.json")
    parser.add_argument("--policy", default="config/pleiotropy_analysis_policy.json")
    parser.add_argument("--out-dir", default="results/pleiotropy/inputs")
    parser.add_argument("--work-dir", default="work/pleiotropy")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest = root / args.manifest
    lock = json.loads((root / args.lock).read_text(encoding="utf-8"))
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if (
        lock.get("schema_version") != "sleep-atlas-pleiotropy-pairs.1"
        or lock.get("ready_pair_count") != policy["expected_sleep_non_sleep_pairs"]
        or lock.get("blocked_pair_count") != 0
        or lock.get("script_sha256") != sha256(root / "scripts/43_prepare_pleiotropy_pairs.py")
    ):
        raise SystemExit("ERROR: complete immutable pleiotropy pair family is required")
    if lock.get("analysis_id") != policy.get("analysis_id"):
        raise SystemExit("ERROR: pleiotropy policy and pair lock analysis IDs differ")
    if lock.get("policy_sha256") != sha256(policy_path):
        raise SystemExit("ERROR: pleiotropy analysis policy differs from its pair lock")
    if sha256(manifest) != lock["manifest_sha256"]:
        raise SystemExit("ERROR: pleiotropy pair manifest differs from its lock")
    matches = [row for row in read_tsv(manifest) if row["pair_id"] == args.pair_id]
    if len(matches) != 1:
        raise SystemExit("ERROR: pair_id must identify one locked pair")
    row = matches[0]
    if row["input_status"] != "READY_FULL_SUMSTATS":
        raise SystemExit(f"ERROR: pair is not full-input ready: {row['blocker']}")
    sources = [root / row["sleep_sumstats"], root / row["non_sleep_sumstats"]]
    expected_hashes = [row["sleep_sumstats_sha256"], row["non_sleep_sumstats_sha256"]]
    for path, expected in zip(sources, expected_hashes):
        if not path.is_file() or sha256(path) != expected:
            raise SystemExit(f"ERROR: locked source drifted: {path}")
    qcs = [path.with_name(path.name.replace(".harmonized.tsv.gz", ".qc.txt")) for path in sources]
    z_squared_maximum = float(policy["placo_z_squared_maximum"])
    if not math.isfinite(z_squared_maximum) or z_squared_maximum <= 0:
        raise SystemExit("ERROR: invalid locked PLACO squared-Z limit")
    estimated_work_bytes = 4 * 1024**3 + 180 * sum(rows_out(path) for path in qcs)
    observed_free = shutil.disk_usage(root).free
    print(f"Estimated pair materialization work space: {estimated_work_bytes} bytes")
    print(f"Observed free space: {observed_free} bytes")
    if not args.materialize:
        print("No pair materialization requested. Re-run with --materialize on the production host.")
        return 0
    if observed_free < estimated_work_bytes:
        raise SystemExit("ERROR: insufficient free space for atomic full-resolution pair materialization")

    work = root / args.work_dir
    work.mkdir(parents=True, exist_ok=True)
    database = work / f"{args.pair_id}.sqlite"
    destination = root / args.out_dir / f"{args.pair_id}.tsv.gz"
    provenance_path = root / args.out_dir / f"{args.pair_id}.provenance.json"
    if destination.exists() or provenance_path.exists():
        raise SystemExit("ERROR: immutable materialized pleiotropy pair already exists")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    database.unlink(missing_ok=True)
    temporary.unlink(missing_ok=True)
    counts = {"direct": 0, "swapped": 0, "complement": 0, "complement_swapped": 0,
              "coordinate_mismatch": 0, "allele_mismatch": 0,
              "extreme_z_excluded": 0, "written": 0}
    connection = None
    try:
        connection = sqlite3.connect(database)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=FILE")
        source_counts = [load_table(connection, "first", sources[0]), load_table(connection, "second", sources[1])]
        query = (
            "SELECT first.snp, first.chr, first.bp, first.a1, first.a2, first.z, second.z, "
            "first.p, second.p, first.n, second.n, second.chr, second.bp, second.a1, second.a2 "
            "FROM first JOIN second USING (snp) ORDER BY first.chr, first.bp, first.snp"
        )
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_out:
                    writer = csv.writer(text_out, delimiter="\t", lineterminator="\n")
                    writer.writerow(OUTPUT_FIELDS)
                    for values in connection.execute(query):
                        snp, chromosome, position, a1, a2, z1, z2, p1, p2, n1, n2, chr2, bp2, b1, b2 = values
                        if chromosome != chr2 or position != bp2:
                            counts["coordinate_mismatch"] += 1
                            continue
                        aligned = align(a1, a2, b1, b2, z2)
                        if aligned is None:
                            counts["allele_mismatch"] += 1
                            continue
                        aligned_z2, mode = aligned
                        counts[mode.lower()] += 1
                        if z1 * z1 > z_squared_maximum or aligned_z2 * aligned_z2 > z_squared_maximum:
                            counts["extreme_z_excluded"] += 1
                            continue
                        writer.writerow([
                            snp, chromosome, position, a1, a2,
                            f"{z1:.15g}", f"{aligned_z2:.15g}", f"{p1:.15g}",
                            f"{p2:.15g}", f"{n1:.15g}", f"{n2:.15g}",
                        ])
                        counts["written"] += 1
        connection.close()
        connection = None
        os.replace(temporary, destination)
        provenance = {
            "analysis_id": "atlas-v1.0-pleiotropic-loci",
            "pair_id": args.pair_id,
            "manifest_sha256": sha256(manifest),
            "policy": str(policy_path.relative_to(root)),
            "policy_sha256": sha256(policy_path),
            "placo_z_squared_maximum": z_squared_maximum,
            "source_files": [str(path.relative_to(root)) for path in sources],
            "source_sha256": expected_hashes,
            "source_rows": source_counts,
            "alignment_counts": counts,
            "output": str(destination.relative_to(root)),
            "output_bytes": destination.stat().st_size,
            "output_sha256": sha256(destination),
        }
        temp_provenance = provenance_path.with_suffix(provenance_path.suffix + ".tmp")
        temp_provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temp_provenance.replace(provenance_path)
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)
        database.unlink(missing_ok=True)
    print(f"Published {args.pair_id}: {counts['written']} aligned full-resolution variants")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
