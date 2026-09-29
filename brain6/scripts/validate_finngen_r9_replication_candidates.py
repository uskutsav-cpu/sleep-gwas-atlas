#!/usr/bin/env python3
"""Validate metadata-only FinnGen R9 candidates without admitting GWAS data."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = Path("brain6/qc/replication_sources/finngen_r9_manifest.tsv")
CANDIDATES = Path("brain6/qc/finngen_r9_replication_candidates.tsv")
PROVENANCE = Path("brain6/qc/finngen_r9_replication_candidates.provenance.json")
EXPECTED = {
    "F5_DEPRESSIO": ("insomnia__mdd", 43280, 329192),
    "F5_BIPO": ("longsleep__bipolar", 7006, 329192),
    "F5_SCHZPHR": ("longsleep__scz", 6515, 364160),
    "G6_PARKINSON": ("longsleep__parkinson", 4235, 373042),
}
FORM_REQUIRED = "OFFICIAL_DOWNLOAD_FORM_REQUIRED; no GWAS bytes acquired"
EXPECTED_MANIFEST_BYTES = 661463
EXPECTED_MANIFEST_SHA256 = "13f3866a224e35471abc0290cba5ed151aad7a4e270ae70b3eacf3a92cf124a2"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate(root: Path = ROOT) -> dict[str, Any]:
    root = root.resolve()
    manifest_path, candidate_path, provenance_path = (
        root / MANIFEST, root / CANDIDATES, root / PROVENANCE
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("status") != "CANDIDATE_SOURCE_METADATA_VERIFIED_ACCESS_FORM_REQUIRED_NO_GWAS_DATA_ACQUIRED":
        raise ValueError("R9 audit is not explicitly metadata-only")
    manifest_meta = provenance.get("sources", {})
    if (manifest_meta.get("manifest_path") != str(MANIFEST) or
            manifest_meta.get("manifest_bytes") != manifest_path.stat().st_size or
            manifest_meta.get("manifest_sha256") != _sha256(manifest_path)):
        raise ValueError("R9 manifest does not match its recorded size/hash")
    if (manifest_path.stat().st_size != EXPECTED_MANIFEST_BYTES or
            _sha256(manifest_path) != EXPECTED_MANIFEST_SHA256 or
            manifest_meta.get("manifest_sha256") != EXPECTED_MANIFEST_SHA256):
        raise ValueError("R9 manifest differs from the frozen, independently verified release manifest")
    table_meta = provenance.get("candidate_table", {})
    if (table_meta.get("path") != str(CANDIDATES) or
            table_meta.get("bytes") != candidate_path.stat().st_size or
            table_meta.get("sha256") != _sha256(candidate_path)):
        raise ValueError("R9 candidate table does not match its recorded size/hash")
    if not provenance.get("sources", {}).get("official_download_instructions"):
        raise ValueError("Official access workflow URL is missing")

    with manifest_path.open(encoding="utf-8", newline="") as stream:
        manifest_rows = {
            row["phenocode"]: row for row in csv.DictReader(stream, delimiter="\t")
        }
    rows = _read_tsv(candidate_path)
    if len(rows) != 4 or table_meta.get("rows") != len(rows):
        raise ValueError("Expected exactly four source-metadata candidates")
    if {row["phenocode"] for row in rows} != set(EXPECTED):
        raise ValueError("Candidate endpoint set differs from the preselected pair mappings")

    head_meta = provenance.get("endpoint_head_metadata", {})
    for row in rows:
        code = row["phenocode"]
        pair, cases, controls = EXPECTED[code]
        manifest_row = manifest_rows.get(code)
        observed_head = head_meta.get(code, {})
        if manifest_row is None:
            raise ValueError(f"Endpoint absent from frozen R9 manifest: {code}")
        if row["pair_id"] != pair or (int(row["cases"]), int(row["controls"])) != (cases, controls):
            raise ValueError(f"Pair mapping or case/control counts changed for {code}")
        if (int(manifest_row["num_cases"]), int(manifest_row["num_controls"])) != (cases, controls):
            raise ValueError(f"Manifest case/control counts disagree for {code}")
        if row["summary_statistics_url"] != manifest_row["path_https"]:
            raise ValueError(f"Candidate URL differs from the official manifest for {code}")
        expected_eff_n = 4 * cases * controls / (cases + controls)
        if abs(float(row["effective_N_4casecontrol"]) - expected_eff_n) > 0.02:
            raise ValueError(f"Effective-N approximation is inconsistent for {code}")
        if row["access_status"] != FORM_REQUIRED or row["analysis_status"] != "CANDIDATE_ONLY_NOT_ADMITTED":
            raise ValueError(f"Unapproved candidate/admission status for {code}")
        if (observed_head.get("downloaded") is not False or
                observed_head.get("url") != manifest_row["path_https"] or
                observed_head.get("http_status") != 200 or
                str(observed_head.get("Content-Length")) != row["compressed_bytes"] or
                str(observed_head.get("x-goog-generation")) != row["gcs_generation"] or
                observed_head.get("ETag", "").strip('"') != row["gcs_etag_md5"] or
                observed_head.get("Last-Modified") != row["last_modified"]):
            raise ValueError(f"Object-header metadata or no-download state changed for {code}")

    return {
        "status": "PASS_R9_METADATA_ONLY_CANDIDATES",
        "candidate_endpoints": len(rows),
        "summary_statistics_acquired": 0,
        "replication_estimates_admitted": 0,
        "official_access_form_required": True,
        "manifest_sha256": manifest_meta["manifest_sha256"],
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
