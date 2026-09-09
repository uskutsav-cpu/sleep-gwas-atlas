#!/usr/bin/env python3
"""Fail-closed validation of the frozen discovery-extension panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


REQUIRED_COLUMNS = {
    "extension_trait_id",
    "phenotype_name",
    "phenotype_definition",
    "source",
    "study_accession",
    "PMID",
    "DOI",
    "ancestry",
    "sample_size",
    "cases",
    "controls",
    "build",
    "source_url",
    "checksum",
    "source_tabix_url",
    "source_tabix_checksum",
    "binary_or_continuous",
    "available_beta",
    "available_se",
    "available_effect_allele",
    "available_other_allele",
    "available_frequency",
    "available_info",
    "prior_sleep_screen_coverage",
    "prior_direct_sleep_genetics_evidence",
    "novelty_priority",
    "selection_reason",
}
ALLOWED_NOVELTY = {
    "HEAVILY_STUDIED",
    "PREVIOUSLY_SCREENED",
    "RELATED_EVIDENCE_ONLY",
    "UNDEREXPLORED",
    "HIGH_NOVELTY_PRIORITY",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--lock",
        type=Path,
        default=Path("discovery_extension/config/extension_panel.lock.json"),
    )
    parser.add_argument(
        "--core-verifier",
        type=Path,
        default=Path("discovery_extension/scripts/00_verify_core_checkpoint.py"),
    )
    args = parser.parse_args()

    core = subprocess.run(
        ["python3", str(args.core_verifier)], check=False, text=True, capture_output=True
    )
    if core.returncode:
        raise SystemExit(f"ERROR: immutable core verification failed:\n{core.stdout}{core.stderr}")

    with args.panel.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = set(reader.fieldnames or [])
    lock = json.loads(args.lock.read_text())
    missing_columns = sorted(REQUIRED_COLUMNS - columns)
    if missing_columns:
        raise SystemExit(f"ERROR: panel missing required columns: {missing_columns}")
    if len(rows) != lock["extension_trait_count"] or not 75 <= len(rows) <= 125:
        raise SystemExit("ERROR: panel size does not match the 75-125-trait locked contract")
    ids = [row["extension_trait_id"] for row in rows]
    if len(set(ids)) != len(ids):
        raise SystemExit("ERROR: duplicate extension trait IDs")
    if ids != lock["extension_trait_ids_in_locked_order"]:
        raise SystemExit("ERROR: panel membership/order differs from lock")
    if lock["planned_raw_rg_test_count"] != 12 * len(rows):
        raise SystemExit("ERROR: planned raw rg test count is not 12 x extension traits")
    expected_hash = lock["artifact_hashes"][str(args.panel)]
    if sha256(args.panel) != expected_hash:
        raise SystemExit("ERROR: candidate_traits.tsv byte hash differs from lock")
    for path, expected in lock["artifact_hashes"].items():
        artifact = Path(path)
        if not artifact.is_file() or sha256(artifact) != expected:
            raise SystemExit(f"ERROR: locked input changed or missing: {artifact}")

    for row in rows:
        if float(row["source_h2_z"]) < 4.0:
            raise SystemExit(f"ERROR: source h2 Z below 4 for {row['extension_trait_id']}")
        if float(row["source_LDSC_intercept"]) > 1.2:
            raise SystemExit(f"ERROR: source intercept above 1.2 for {row['extension_trait_id']}")
        if row["novelty_priority"] not in ALLOWED_NOVELTY:
            raise SystemExit(f"ERROR: invalid novelty category for {row['extension_trait_id']}")
        if row["binary_or_continuous"] == "binary":
            if int(row["cases"]) < 1000 or int(row["controls"]) < 1000:
                raise SystemExit(f"ERROR: binary count threshold failed for {row['extension_trait_id']}")
        elif int(row["sample_size"]) < 10000:
            raise SystemExit(f"ERROR: quantitative sample threshold failed for {row['extension_trait_id']}")
        if not row["source_url"].startswith("https://") or not row["checksum"].startswith("md5:"):
            raise SystemExit(f"ERROR: missing persistent URL/checksum for {row['extension_trait_id']}")
        if not row["source_tabix_url"].startswith("https://") or not row[
            "source_tabix_checksum"
        ].startswith("md5:"):
            raise SystemExit(
                f"ERROR: missing persistent tabix URL/checksum for {row['extension_trait_id']}"
            )
        if int(row["source_file_size_bytes"]) <= 0 or int(row["source_tabix_size_bytes"]) <= 0:
            raise SystemExit(f"ERROR: invalid source byte size for {row['extension_trait_id']}")

    print(
        f"EXTENSION_PANEL_OK traits={len(rows)} planned_rg={lock['planned_raw_rg_test_count']} "
        f"panel_sha256={sha256(args.panel)}"
    )


if __name__ == "__main__":
    main()
