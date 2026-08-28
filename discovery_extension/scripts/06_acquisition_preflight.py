#!/usr/bin/env python3
"""Quantify storage required for the locked full-resolution extension sources."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


VARIANT_MANIFEST = {
    "source_filename": "full_variant_qc_metrics.txt.bgz",
    "source_url": (
        "https://pan-ukb-us-east-1.s3.amazonaws.com/"
        "sumstats_release/full_variant_qc_metrics.txt.bgz"
    ),
    "checksum": "etag:e70ebc8289f762dd8d5086f54e766654",
    "source_file_size_bytes": 2701503051,
}
VARIANT_MANIFEST_TABIX = {
    "source_filename": "full_variant_qc_metrics.txt.bgz.tbi",
    "source_url": (
        "https://pan-ukb-us-east-1.s3.amazonaws.com/"
        "sumstats_release/full_variant_qc_metrics.txt.bgz.tbi"
    ),
    "checksum": "etag:555b93fcb299e59c6e75a5a41ec3de40",
    "source_file_size_bytes": 2261519,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--plan-out",
        type=Path,
        default=Path("discovery_extension/results/acquisition_plan.tsv"),
    )
    parser.add_argument(
        "--preflight-out",
        type=Path,
        default=Path("discovery_extension/provenance/acquisition_preflight.json"),
    )
    parser.add_argument("--storage-path", type=Path, default=Path("."))
    parser.add_argument("--safety-factor", type=float, default=1.15)
    args = parser.parse_args()

    with args.panel.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    phenotype_total = sum(int(row["source_file_size_bytes"]) for row in rows)
    tabix_total = sum(int(row["source_tabix_size_bytes"]) for row in rows)
    shared_variant_total = (
        VARIANT_MANIFEST["source_file_size_bytes"]
        + VARIANT_MANIFEST_TABIX["source_file_size_bytes"]
    )
    total = phenotype_total + tabix_total + shared_variant_total
    required = int(total * args.safety_factor)
    usage = shutil.disk_usage(args.storage_path.resolve())
    enough = usage.free >= required
    status = "READY" if enough else "BLOCKED_INSUFFICIENT_STORAGE"

    plan_fields = [
        "extension_trait_id",
        "source_filename",
        "source_url",
        "checksum",
        "source_file_size_bytes",
        "source_tabix_filename",
        "source_tabix_url",
        "source_tabix_checksum",
        "source_tabix_size_bytes",
        "acquisition_status",
        "retention_policy",
    ]
    args.plan_out.parent.mkdir(parents=True, exist_ok=True)
    with args.plan_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=plan_fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "extension_trait_id": row["extension_trait_id"],
                    "source_filename": row["source_filename"],
                    "source_url": row["source_url"],
                    "checksum": row["checksum"],
                    "source_file_size_bytes": row["source_file_size_bytes"],
                    "source_tabix_filename": row["source_tabix_filename"],
                    "source_tabix_url": row["source_tabix_url"],
                    "source_tabix_checksum": row["source_tabix_checksum"],
                    "source_tabix_size_bytes": row["source_tabix_size_bytes"],
                    "acquisition_status": status,
                    "retention_policy": "retain_full_resolution_source_and_checksum",
                }
            )
        for shared in (VARIANT_MANIFEST, VARIANT_MANIFEST_TABIX):
            writer.writerow(
                {
                    "extension_trait_id": "SHARED_VARIANT_REFERENCE",
                    "source_filename": shared["source_filename"],
                    "source_url": shared["source_url"],
                    "checksum": shared["checksum"],
                    "source_file_size_bytes": shared["source_file_size_bytes"],
                    "source_tabix_filename": "NA",
                    "source_tabix_url": "NA",
                    "source_tabix_checksum": "NA",
                    "source_tabix_size_bytes": "0",
                    "acquisition_status": status,
                    "retention_policy": "retain_full_resolution_source_and_checksum",
                }
            )

    report = {
        "schema_version": "1.0.0",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "panel_trait_count": len(rows),
        "storage_path": str(args.storage_path.resolve()),
        "phenotype_sumstats_bytes": phenotype_total,
        "phenotype_tabix_bytes": tabix_total,
        "shared_variant_reference_bytes": shared_variant_total,
        "compressed_source_bytes": total,
        "compressed_source_gib": round(total / 1024**3, 3),
        "safety_factor": args.safety_factor,
        "required_free_bytes": required,
        "required_free_gib": round(required / 1024**3, 3),
        "available_free_bytes": usage.free,
        "available_free_gib": round(usage.free / 1024**3, 3),
        "shortfall_bytes": max(0, required - usage.free),
        "shortfall_gib": round(max(0, required - usage.free) / 1024**3, 3),
        "status": status,
        "download_started": False,
        "reason": (
            "No source download was started because the full-resolution retention contract "
            "cannot be satisfied on the current volume."
            if not enough
            else "The storage preflight passed; acquisition remains a separately authorized step."
        ),
    }
    args.preflight_out.parent.mkdir(parents=True, exist_ok=True)
    args.preflight_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        f"ACQUISITION_PREFLIGHT_{status} source_gib={report['compressed_source_gib']} "
        f"required_gib={report['required_free_gib']} available_gib={report['available_free_gib']} "
        f"shortfall_gib={report['shortfall_gib']}"
    )


if __name__ == "__main__":
    main()
