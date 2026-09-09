#!/usr/bin/env python3
"""Checksum-lock the exact PLACO+ pair scans before result access."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_file(row: dict[str, str], path_field: str, checksum_field: str) -> None:
    path = Path(row[path_field])
    if not path.is_file() or len(row[checksum_field]) != 64 or sha256(path) != row[checksum_field]:
        raise SystemExit(f"ERROR: {row['pair_id']} missing or checksum-mismatched {path_field}")


def validate_genomewide_input(path: Path, pair_id: str) -> int:
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    required = {"SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2", "MAF", "INFO1", "INFO2"}
    count = 0
    previous: tuple[int, int, str] | None = None
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"ERROR: {pair_id} merged genome-wide input lacks fields: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            count += 1
            coordinate = (int(row["CHR"]), int(row["BP"]), row["SNP"])
            if not 1 <= coordinate[0] <= 22 or coordinate[1] < 1 or (previous is not None and coordinate <= previous):
                raise SystemExit(f"ERROR: {pair_id} merged genome-wide input is unsorted, duplicated, or non-autosomal")
            previous = coordinate
            if row["A1"] not in {"A", "C", "G", "T"} or row["A2"] not in {"A", "C", "G", "T"} or row["A1"] == row["A2"]:
                raise SystemExit(f"ERROR: {pair_id} merged genome-wide input has invalid alleles")
            values = [float(row[field]) for field in ("Z1", "Z2", "P1", "P2", "MAF", "INFO1", "INFO2")]
            if not all(math.isfinite(value) for value in values):
                raise SystemExit(f"ERROR: {pair_id} merged genome-wide input has non-finite values")
            z1, z2, p1, p2, maf, info1, info2 = values
            if not (0 <= p1 <= 1 and 0 <= p2 <= 1 and 0.01 <= maf <= 0.5 and info1 >= 0.9 and info2 >= 0.9):
                raise SystemExit(f"ERROR: {pair_id} merged genome-wide input violates P/MAF/INFO filters")
            if z1**2 > 80 or z2**2 > 80:
                raise SystemExit(f"ERROR: {pair_id} primary merged input retains Z squared above 80")
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/pleiotropy/pleiotropy_input_queue.tsv"))
    parser.add_argument("--candidate-lock", type=Path, default=Path("discovery_extension/config/pleiotropy_candidate_family.lock.json"))
    parser.add_argument("--readiness", type=Path, default=Path("discovery_extension/results/pleiotropy/pleiotropy_readiness.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/config/pleiotropy_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/pleiotropy_manifest.lock.json"))
    parser.add_argument("--minimum-genomewide-variants", type=int, default=500000)
    args = parser.parse_args()

    fields, rows = read_tsv(args.queue)
    family = json.loads(args.candidate_lock.read_text(encoding="utf-8"))
    readiness_rows = read_tsv(args.readiness)[1]
    if not rows:
        raise SystemExit("ERROR: no replicated high-priority pair is available for PLACO+")
    if [row["pair_id"] for row in rows] != family.get("pair_ids_in_locked_order") or len(rows) != family.get("pair_count"):
        raise SystemExit("ERROR: PLACO+ queue differs from its candidate-family lock")
    if family.get("results_accessed_before_lock") is not False:
        raise SystemExit("ERROR: PLACO+ candidate family was not result-free")
    if len(readiness_rows) != 1 or readiness_rows[0]["code_status"] != "PASS":
        raise SystemExit("ERROR: pinned PLACO+ code is not ready")
    for row in rows:
        if row["selection_status"] != "SELECTED_REPLICATED_PRIORITY_PAIR" or row["replication_class"] != "REPLICATED" or row["priority_tier"] != "B":
            raise SystemExit(f"ERROR: ineligible PLACO+ pair: {row['pair_id']}")
        if row["results_accessed_before_lock"] != "NO":
            raise SystemExit(f"ERROR: PLACO+ results accessed before input lock: {row['pair_id']}")
        if row["input_scope"] != "FULL_GENOME" or row["build"] != "GRCh37" or row["ancestry"] != "EUR":
            raise SystemExit(f"ERROR: PLACO+ full-genome/build/ancestry contract failed: {row['pair_id']}")
        if row["schema_status"] != "VERIFIED" or row["effect_allele_alignment_status"] != "VERIFIED":
            raise SystemExit(f"ERROR: PLACO+ schema/effect alignment failed: {row['pair_id']}")
        if row["maf_info_filter_status"] != "MAF_GE_0.01_INFO_GE_0.9" or row["z_squared_filter_status"] != "PRIMARY_Z2_LE_80_WITH_SEPARATE_LARGE_EFFECT_LEDGER":
            raise SystemExit(f"ERROR: PLACO+ variant filter contract failed: {row['pair_id']}")
        if int(row["eligible_variant_count"]) < args.minimum_genomewide_variants:
            raise SystemExit(f"ERROR: PLACO+ input is not demonstrably genome-wide: {row['pair_id']}")
        if int(row["input_variant_count"]) != int(row["eligible_variant_count"]) + int(row["z2_excluded_count"]) or int(row["z2_excluded_count"]) < 0:
            raise SystemExit(f"ERROR: PLACO+ input/eligible/Z2-exclusion ledger is inconsistent: {row['pair_id']}")
        for field in ("clumping_parameters", "ld_reference_id", "curator", "curation_date"):
            if not row[field].strip() or row[field] == "PENDING":
                raise SystemExit(f"ERROR: PLACO+ pair lacks {field}: {row['pair_id']}")
        for path_field, checksum_field in (
            ("merged_genomewide_path", "merged_genomewide_sha256"),
            ("ld_reference_path", "ld_reference_sha256"),
            ("placo_source_path", "placo_source_sha256"),
        ):
            verify_file(row, path_field, checksum_field)
        if validate_genomewide_input(Path(row["merged_genomewide_path"]), row["pair_id"]) != int(row["eligible_variant_count"]):
            raise SystemExit(f"ERROR: {row['pair_id']} eligible variant count differs from merged input")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    primary_threshold = 5e-8 / len(rows)
    lock = {
        "schema_version": "1.0.0", "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "results_accessed_before_lock": False, "manifest_sha256": sha256(args.out),
        "pair_ids_in_locked_order": [row["pair_id"] for row in rows], "pair_scan_count": len(rows),
        "minimum_genomewide_variants": args.minimum_genomewide_variants,
        "conventional_variant_threshold": 5e-8, "primary_locked_family_threshold": primary_threshold,
        "nuisance_marginal_p_threshold": 1e-4, "primary_z_squared_maximum": 80,
        "candidate_family_lock_sha256": sha256(args.candidate_lock),
        "queue_sha256": sha256(args.queue), "readiness_sha256": sha256(args.readiness),
        "claim_limit": "statistical pleiotropic association; not a shared causal variant or causal mechanism",
    }
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    args.lock.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PLEIOTROPY_MANIFEST_LOCKED pairs={len(rows)} threshold={primary_threshold:.12g} result_free=true")


if __name__ == "__main__":
    main()
