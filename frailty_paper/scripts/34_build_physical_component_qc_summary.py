#!/usr/bin/env python3
"""Consolidate structural QC and unresolved source gates for five components."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path


EXPECTED = {
    "zenodo_14011550_weight_loss",
    "zenodo_14011550_exhaustion",
    "zenodo_14011550_low_physical_activity",
    "zenodo_14011550_slow_walking_speed",
    "zenodo_14011550_low_grip_strength",
}
REQUIRED_STREAM_FIELDS = ["resource_id", "file", "manifest_bytes", "manifest_sha256", "genome_build", "ancestry", "sample_size", "stream_status", "header_status", "header", "data_rows", "rows_missing_variant", "rows_missing_chromosome", "rows_invalid_position", "rows_missing_alleles", "rows_missing_beta_or_se", "rows_invalid_beta_or_se", "rows_missing_p", "rows_invalid_p_or_log10p", "rows_invalid_eaf"]
REQUIRED_DUP_FIELDS = ["resource_id", "file", "bytes", "sha256", "rows", "coordinate_order", "coordinate_order_descents", "duplicate_key_groups", "duplicate_excess_rows", "global_key_uniqueness", "key_definition"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def unique_index(rows: list[dict[str, str]], key: str, source: str) -> dict[str, dict[str, str]]:
    out = {}
    for row in rows:
        value = row[key].strip()
        if not value or value in out:
            raise ValueError(f"Missing or duplicate {key} in {source}: {value!r}")
        out[value] = row
    return out


def assemble(repo: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    registry_path = repo / "frailty_paper/manifests/all_acquired_resources.tsv"
    stream_path = repo / "frailty_paper/manifests/gwas_source_scan.tsv"
    duplicate_path = repo / "frailty_paper/manifests/gwas_component_duplicate_audit.tsv"
    registry_rows = [r for r in read_rows(registry_path) if r.get("resource_id") in EXPECTED]
    stream_rows = [r for r in read_rows(stream_path) if r.get("resource_id") in EXPECTED]
    duplicate_rows = [r for r in read_rows(duplicate_path) if r.get("resource_id") in EXPECTED]
    registry = unique_index(registry_rows, "resource_id", str(registry_path))
    stream = unique_index(stream_rows, "resource_id", str(stream_path))
    duplicates = unique_index(duplicate_rows, "resource_id", str(duplicate_path))
    for label, data in (("registry", registry), ("stream scan", stream), ("duplicate audit", duplicates)):
        if set(data) != EXPECTED:
            raise ValueError(f"{label} resource IDs mismatch; missing={sorted(EXPECTED-set(data))}; extra={sorted(set(data)-EXPECTED)}")
    result = []
    for rid in sorted(EXPECTED):
        reg, scan, dup = registry[rid], stream[rid], duplicates[rid]
        for field in REQUIRED_STREAM_FIELDS:
            if field not in scan:
                raise ValueError(f"Missing stream field {field} for {rid}")
        for field in REQUIRED_DUP_FIELDS:
            if field not in dup:
                raise ValueError(f"Missing duplicate-audit field {field} for {rid}")
        if scan["manifest_sha256"] != reg["sha256"] or dup["sha256"] != reg["sha256"]:
            raise ValueError(f"Input checksum disagreement for {rid}")
        if scan["manifest_bytes"] != reg["bytes"] or dup["bytes"] != reg["bytes"]:
            raise ValueError(f"Input size disagreement for {rid}")
        if scan["data_rows"] != dup["rows"]:
            raise ValueError(f"Input row-count disagreement for {rid}")
        if scan["stream_status"] != "GZIP_EOF_OK" or scan["header_status"] != "REQUIRED_COLUMNS_PRESENT":
            raise ValueError(f"Structural scan did not pass for {rid}")
        if dup["global_key_uniqueness"] != "UNIQUE" or dup["duplicate_key_groups"] != "0" or dup["duplicate_excess_rows"] != "0":
            raise ValueError(f"Duplicate-key audit did not pass for {rid}")
        trait = reg["trait"]
        result.append({
            "resource_id": rid,
            "trait": trait,
            "file": reg["file"],
            "bytes": reg["bytes"],
            "sha256": reg["sha256"],
            "sample_size": reg["sample_size"],
            "cases": reg["cases"],
            "controls": reg["controls"],
            "genome_build": reg["genome_build"],
            "ancestry": reg["ancestry"],
            "stream_status": scan["stream_status"],
            "header": scan["header"],
            "data_rows": scan["data_rows"],
            "rows_missing_variant": scan["rows_missing_variant"],
            "rows_missing_chromosome": scan["rows_missing_chromosome"],
            "rows_invalid_position": scan["rows_invalid_position"],
            "rows_missing_alleles": scan["rows_missing_alleles"],
            "rows_missing_beta_or_se": scan["rows_missing_beta_or_se"],
            "rows_invalid_beta_or_se": scan["rows_invalid_beta_or_se"],
            "rows_missing_p": scan["rows_missing_p"],
            "rows_invalid_p_or_log10p": scan["rows_invalid_p_or_log10p"],
            "rows_invalid_eaf": scan["rows_invalid_eaf"],
            "duplicate_key_groups": dup["duplicate_key_groups"],
            "duplicate_excess_rows": dup["duplicate_excess_rows"],
            "coordinate_order": dup["coordinate_order"],
            "effect_allele_mapping": "UNVERIFIED (header ALLELE0/ALLELE1 and BETA; coding not specified)",
            "structural_qc": "PASS (file integrity, required fields, value domains, unique coordinate/allele keys)",
            "analysis_eligibility": "BLOCKED",
            "eligibility_reason": "Genome build, ancestry/QC, signed effect model and exact generating-study provenance unresolved; do not harmonize or run LDSC.",
        })
    input_hashes = {p.relative_to(repo).as_posix(): sha256(p) for p in (registry_path, stream_path, duplicate_path)}
    return result, input_hashes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = (args.output or repo / "frailty_paper/manifests/physical_component_qc.tsv").resolve()
    rows, input_hashes = assemble(repo)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "artifact": output.relative_to(repo).as_posix(),
        "row_count": len(rows),
        "resource_ids": [r["resource_id"] for r in rows],
        "operation": "aggregation of existing structural and duplicate-key audits only; raw GWAS records are not changed or reanalyzed",
        "script": Path(__file__).resolve().relative_to(repo).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "command": "python " + Path(__file__).resolve().relative_to(repo).as_posix() + (" " + " ".join(sys.argv[1:]) if sys.argv[1:] else ""),
        "python_version": platform.python_version(),
        "input_sha256": input_hashes,
        "output_sha256": sha256(output),
        "interpretation": "All five files pass recorded structural checks and duplicate-key audits, but none passes the scientific provenance gate for harmonization/LDSC.",
    }
    manifest_path = output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "output": str(output), "manifest": str(manifest_path)}))


if __name__ == "__main__":
    main()
