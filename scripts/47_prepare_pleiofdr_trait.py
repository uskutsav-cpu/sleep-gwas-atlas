#!/usr/bin/env python3
"""Align one full GWAS to the pinned pleioFDR template and create a MAT input."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT_SPEC = importlib.util.spec_from_file_location(
    "pleiotropy_preflight", ROOT / "scripts/42_pleiotropy_preflight.py"
)
if PREFLIGHT_SPEC is None or PREFLIGHT_SPEC.loader is None:
    raise RuntimeError("could not load pleiotropy preflight helpers")
preflight = importlib.util.module_from_spec(PREFLIGHT_SPEC)
PREFLIGHT_SPEC.loader.exec_module(preflight)

REQUIRED = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"}
TEMPLATE_FIELDS = ["CHR", "SNP", "GP", "BP", "A1", "A2"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
TRAIT_ID = re.compile(r"[a-z0-9_]+$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows_out(qc: Path) -> int:
    for line in qc.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0] == "rows_out":
            return int(fields[1])
    raise SystemExit(f"ERROR: {qc} lacks rows_out")


def source_rows(path: Path):
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
            except ValueError as exc:
                raise SystemExit(f"ERROR: invalid numeric source row in {path}") from exc
            a1, a2 = row["A1"].upper(), row["A2"].upper()
            if not (
                row["SNP"].lower().startswith("rs")
                and 1 <= chromosome <= 22
                and position > 0
                and len(a1) == len(a2) == 1
                and a1 in "ACGT" and a2 in "ACGT" and a1 != a2
                and math.isfinite(beta) and math.isfinite(se) and se > 0
                and math.isfinite(p_value) and 0 < p_value <= 1
            ):
                raise SystemExit(f"ERROR: invalid source row in {path}: {row.get('SNP', '')}")
            yield row["SNP"].lower(), chromosome, position, a1, a2, beta / se, p_value


def align(ref_a1: str, ref_a2: str, source_a1: str, source_a2: str, z: float) -> float | None:
    if (source_a1, source_a2) == (ref_a1, ref_a2):
        return z
    if (source_a1, source_a2) == (ref_a2, ref_a1):
        return -z
    complement = source_a1.translate(COMPLEMENT), source_a2.translate(COMPLEMENT)
    if complement == (ref_a1, ref_a2):
        return z
    if complement == (ref_a2, ref_a1):
        return -z
    return None


def load_source(connection: sqlite3.Connection, source: Path) -> int:
    connection.execute(
        "CREATE TABLE source (snp TEXT PRIMARY KEY, chr INTEGER NOT NULL, "
        "bp INTEGER NOT NULL, a1 TEXT NOT NULL, a2 TEXT NOT NULL, "
        "z REAL NOT NULL, p REAL NOT NULL) WITHOUT ROWID"
    )
    statement = "INSERT INTO source VALUES (?, ?, ?, ?, ?, ?, ?)"
    count = 0
    batch = []
    try:
        for row in source_rows(source):
            batch.append(row)
            if len(batch) == 50000:
                connection.executemany(statement, batch)
                count += len(batch)
                batch.clear()
        if batch:
            connection.executemany(statement, batch)
            count += len(batch)
        connection.commit()
    except sqlite3.IntegrityError as exc:
        raise SystemExit(f"ERROR: duplicate source SNP in {source}") from exc
    return count


def normalise_mat_header(path: Path) -> None:
    description = b"MATLAB 5.0 MAT-file, Platform: atlas-v1.0, Created deterministically"
    if len(description) > 116:
        raise RuntimeError("deterministic MAT header is too long")
    with path.open("r+b") as handle:
        handle.seek(0)
        handle.write(description.ljust(116, b" "))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trait_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--template", default="ref/pleiofdr/9545380.ref")
    parser.add_argument("--out-dir", default="data/pleiofdr")
    parser.add_argument("--work-dir", default="work/pleiofdr_traits")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    if not TRAIT_ID.fullmatch(args.trait_id):
        raise SystemExit("ERROR: invalid trait_id")
    root = Path(args.root).resolve()
    policy_path = root / "config/pleiotropy_analysis_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    panel = preflight.mixer.read_tsv(root / "config/analysis_panel.tsv")
    matches = [row for row in panel if row["trait_id"] == args.trait_id]
    if len(matches) != 1:
        raise SystemExit("ERROR: trait_id is not in the exact locked panel")
    source, qc, prefilter_strategy = preflight.mixer.choose_harmonized(root, args.trait_id)
    if not prefilter_strategy:
        prefilter_strategy = "not supplied"
    if not source.is_file() or not qc.is_file() or prefilter_strategy != "not supplied":
        raise SystemExit("ERROR: trait lacks full non-HapMap3 post-QC summary statistics")
    template = root / args.template
    if not template.is_file():
        raise SystemExit("ERROR: pinned pleioFDR variant template is absent")
    if (
        template.stat().st_size != policy["pleiofdr_variant_template_bytes"]
        or sha256(template) != policy["pleiofdr_variant_template_sha256"]
    ):
        raise SystemExit("ERROR: pleioFDR variant template differs from its pin")
    expected_source_rows = rows_out(qc)
    estimated_work = 2 * 1024**3 + 5 * source.stat().st_size
    free_bytes = shutil.disk_usage(root).free
    print(f"Estimated trait materialization work space: {estimated_work} bytes")
    print(f"Observed free space: {free_bytes} bytes")
    if not args.materialize:
        print("No MAT materialization requested. Re-run with --materialize on a production host.")
        return 0
    if free_bytes < policy["pleiofdr_trait_materialization_minimum_free_storage_bytes"]:
        raise SystemExit("ERROR: less than the locked 8 GiB free-space floor")
    try:
        import numpy as np
        from scipy.io import savemat
    except ImportError as exc:
        raise SystemExit("ERROR: pinned NumPy/SciPy environment is required") from exc

    work_dir = root / args.work_dir
    work_dir.mkdir(parents=True, exist_ok=True)
    database = work_dir / f"{args.trait_id}.sqlite"
    output = root / args.out_dir / f"{args.trait_id}.mat"
    provenance_path = root / args.out_dir / f"{args.trait_id}.provenance.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    database.unlink(missing_ok=True)
    temporary.unlink(missing_ok=True)
    connection = None
    counts = {
        "template_rows": 0, "matched": 0, "coordinate_mismatch": 0,
        "allele_mismatch": 0, "missing_in_trait": 0,
    }
    try:
        connection = sqlite3.connect(database)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        observed_source_rows = load_source(connection, source)
        if observed_source_rows != expected_source_rows:
            raise SystemExit("ERROR: source row count differs from harmonization QC ledger")
        connection.execute(
            "CREATE TEMP TABLE ref_batch (ord INTEGER PRIMARY KEY, snp TEXT UNIQUE NOT NULL, "
            "chr INTEGER NOT NULL, bp INTEGER NOT NULL, a1 TEXT NOT NULL, a2 TEXT NOT NULL)"
        )
        variant_count = int(policy["pleiofdr_variant_template_variants"])
        logp = np.full(variant_count, np.nan, dtype=np.float64)
        z_score = np.full(variant_count, np.nan, dtype=np.float64)

        def process_batch(batch: list[tuple[int, str, int, int, str, str]]) -> None:
            if not batch:
                return
            connection.execute("DELETE FROM ref_batch")
            try:
                connection.executemany("INSERT INTO ref_batch VALUES (?, ?, ?, ?, ?, ?)", batch)
            except sqlite3.IntegrityError as exc:
                raise SystemExit("ERROR: duplicate SNP within pleioFDR template batch") from exc
            query = (
                "SELECT ref_batch.ord, ref_batch.chr, ref_batch.bp, ref_batch.a1, ref_batch.a2, "
                "source.chr, source.bp, source.a1, source.a2, source.z, source.p "
                "FROM ref_batch LEFT JOIN source USING (snp) ORDER BY ref_batch.ord"
            )
            for values in connection.execute(query):
                order, ref_chr, ref_bp, ref_a1, ref_a2, src_chr, src_bp, src_a1, src_a2, src_z, src_p = values
                if src_chr is None:
                    counts["missing_in_trait"] += 1
                    continue
                if ref_chr != src_chr or ref_bp != src_bp:
                    counts["coordinate_mismatch"] += 1
                    continue
                aligned_z = align(ref_a1, ref_a2, src_a1, src_a2, src_z)
                if aligned_z is None:
                    counts["allele_mismatch"] += 1
                    continue
                z_score[order] = aligned_z
                logp[order] = -math.log10(src_p)
                counts["matched"] += 1

        batch: list[tuple[int, str, int, int, str, str]] = []
        with template.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != TEMPLATE_FIELDS:
                raise SystemExit(f"ERROR: unexpected pleioFDR template header: {reader.fieldnames}")
            for order, row in enumerate(reader):
                try:
                    chromosome, position = int(row["CHR"]), int(row["BP"])
                except ValueError as exc:
                    raise SystemExit("ERROR: invalid pleioFDR template coordinate") from exc
                snp = row["SNP"].lower()
                a1, a2 = row["A1"].upper(), row["A2"].upper()
                if not (
                    snp.startswith("rs") and 1 <= chromosome <= 22 and position > 0
                    and len(a1) == len(a2) == 1 and a1 in "ACGT" and a2 in "ACGT" and a1 != a2
                ):
                    raise SystemExit(f"ERROR: invalid pleioFDR template row: {snp}")
                batch.append((order, snp, chromosome, position, a1, a2))
                if len(batch) == 50000:
                    process_batch(batch)
                    batch.clear()
            process_batch(batch)
            counts["template_rows"] = order + 1 if "order" in locals() else 0
        if counts["template_rows"] != variant_count:
            raise SystemExit("ERROR: pleioFDR template row count differs from its pin")
        if counts["matched"] == 0:
            raise SystemExit("ERROR: no variants aligned to the pleioFDR template")
        savemat(
            temporary,
            {"logpvec": logp.reshape((-1, 1)), "zvec": z_score.reshape((-1, 1))},
            appendmat=False,
            do_compression=True,
            oned_as="column",
        )
        normalise_mat_header(temporary)
        os.replace(temporary, output)
        provenance = {
            "analysis_id": policy["analysis_id"],
            "trait_id": args.trait_id,
            "policy_sha256": sha256(policy_path),
            "source": str(source.relative_to(root)),
            "source_sha256": sha256(source),
            "source_rows": observed_source_rows,
            "template": str(template.relative_to(root)),
            "template_sha256": sha256(template),
            "alignment_counts": counts,
            "output": str(output.relative_to(root)),
            "output_bytes": output.stat().st_size,
            "output_sha256": sha256(output),
            "mat_variables": ["logpvec", "zvec"],
        }
        provenance_temporary = provenance_path.with_suffix(provenance_path.suffix + ".tmp")
        provenance_temporary.write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        provenance_temporary.replace(provenance_path)
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)
        database.unlink(missing_ok=True)
    print(f"pleioFDR trait ready: {args.trait_id}; {counts['matched']} aligned variants")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
