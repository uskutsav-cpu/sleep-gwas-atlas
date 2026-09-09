#!/usr/bin/env python3
"""Create an executable tiny PLINK-in-ZIP fixture and its exact test policy."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import zipfile
from pathlib import Path


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def add_member(archive: zipfile.ZipFile, name: str, payload: bytes) -> None:
    info = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o100644 << 16
    archive.writestr(info, payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory")
    parser.add_argument(
        "--mode", choices=("valid", "invalid-bim", "missing-member"), default="valid",
    )
    args = parser.parse_args()
    root = Path(args.output_directory).resolve()
    root.mkdir(parents=True, exist_ok=True)
    source = root / "fixture_source"
    source.mkdir()

    bed = b"\x6c\x1b\x01\x00\x00\x00\x00"
    allele = "N" if args.mode == "invalid-bim" else "C"
    bim = (
        f"1\trs1\t0\t10\tA\t{allele}\n"
        "2\trs2\t0\t20\tG\tT\n"
        "22\trs3\t0\t30\tC\tT\n"
        "23\trsX\t0\t40\tA\tG\n"
    ).encode("ascii")
    fam = (
        "F1 I1 0 0 1 -9\n"
        "F2 I2 0 0 2 -9\n"
        "F3 I3 0 0 0 -9\n"
    ).encode("ascii")
    members = {
        "bed": {"name": "g1000_eur.bed", "payload": bed, "records": 4},
        "bim": {"name": "g1000_eur.bim", "payload": bim, "records": 4},
        "fam": {"name": "g1000_eur.fam", "payload": fam, "records": 3},
    }
    archive_path = source / "g1000_eur.zip"
    with zipfile.ZipFile(archive_path, "x") as archive:
        for component in ("bed", "bim", "fam"):
            if args.mode == "missing-member" and component == "fam":
                continue
            member = members[component]
            add_member(archive, member["name"], member["payload"])

    plink_path = source / "plink"
    plink_payload = b"#!/bin/sh\nprintf '%s\\n' 'PLINK v1.9.0-test-fixture'\n"
    plink_path.write_bytes(plink_payload)
    os.chmod(plink_path, 0o755)
    source_members = {
        component: {
            "name": member["name"], "expected_bytes": len(member["payload"]),
            "sha256": digest(member["payload"]), "records": member["records"],
        }
        for component, member in members.items()
    }
    extracted_bytes = sum(member["expected_bytes"] for member in source_members.values())
    archive_payload = archive_path.read_bytes()
    policy = {
        "schema_version": "sleep-atlas-track-b-pleiotropy-ld-fixture.1",
        "analysis_id": "track-b-test-pleiotropy-ld",
        "ld_reference_policy": {
            "purpose": "tiny executable test fixture",
            "ancestry": "EUR",
            "build": "GRCh37/hg19",
            "source_release": "TINY_TEST_FIXTURE",
            "source_archive_path": "fixture_source/g1000_eur.zip",
            "source_archive_expected_bytes": len(archive_payload),
            "source_archive_sha256": digest(archive_payload),
            "source_members": source_members,
            "archive_allowed_members": [
                "g1000_eur.bed", "g1000_eur.bim", "g1000_eur.fam",
            ],
            "materialized_directory": "ref/track_b/pleiotropy/g1000_eur",
            "materialized_prefix": "ref/track_b/pleiotropy/g1000_eur/g1000_eur",
            "materialized_manifest": "ref/track_b/pleiotropy/g1000_eur/g1000_eur.manifest.tsv",
            "materialized_provenance": "ref/track_b/pleiotropy/g1000_eur/g1000_eur.provenance.json",
            "staging_root": "ref/track_b/pleiotropy/.staging",
            "materializer_script": "scripts/127_prepare_track_b_pleiotropy_ld.py",
            "staging_safety_bytes": 1024,
            "minimum_free_storage_bytes": extracted_bytes + 1024,
            "conflicting_activity_markers": ["active_lava_download.tmp"],
            "sample_count": 3,
            "variant_count": 4,
            "reference_chromosomes": [1, 2, 22, 23],
            "analysis_autosomes": [1, 2, 22],
            "nonanalysis_chromosomes_retained_for_source_fidelity": [23],
            "expected_autosomal_variant_count": 3,
            "expected_reference_chromosome_variant_counts": {
                "1": 1, "2": 1, "22": 1, "23": 1,
            },
            "grch37_chromosome_lengths": {
                "1": 1000, "2": 1000, "22": 1000, "23": 1000,
            },
            "full_genome_wide_reference_required": True,
            "hapmap3_subset_forbidden": True,
            "raw_reference_fidelity_rule": "Seal untouched chr1-23 fixture",
            "analysis_chromosome_rule": "Only fixture autosomes are eligible; chr23/X is forbidden",
            "materialization_transform_rule": "Byte-identical extraction only; no LD-block splitting or SNP reduction",
            "materialization_resource_metrics_required": [
                "staging_extraction_validation_wall_seconds", "peak_rss_bytes",
                "streaming_chunk_bytes", "extracted_uncompressed_bytes",
            ],
            "plink_path": "fixture_source/plink",
            "plink_version": "PLINK v1.9.0-test-fixture",
            "plink_sha256": digest(plink_payload),
            "clump_r2_strict_upper_bound": 0.1,
            "clump_window_kb": 1000,
        },
    }
    (root / "fixture_policy.json").write_text(
        json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    print(root / "fixture_policy.json")


if __name__ == "__main__":
    main()
