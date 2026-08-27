#!/usr/bin/env python3
"""Validate and fingerprint the immutable atlas-v1.0 trait panel.

The lock intentionally covers the ordered trait IDs, not every metadata byte.
Source metadata must be curated over time, but changing a phenotype identity or
its analysis order requires an explicit lock-file update and code review.
"""
import argparse
import csv
import hashlib
import json
import os
import re


DEFAULT_MANIFEST = "config/analysis_panel.tsv"
DEFAULT_LOCK = "config/analysis_panel.lock.json"
REQUIRED_COLUMNS = {
    "panel_version", "trait_id", "label", "domain", "type",
    "phenotype_definition", "source_note", "dataset_version", "source_id",
    "raw_file", "ncase", "ncontrol", "n_total", "pop_prev", "build",
    "source_status", "pmid", "doi", "ancestry", "pop_prev_citation",
}
SOURCE_STATES = {"SOURCE_PENDING", "SOURCE_VERIFIED"}
TRAIT_TYPES = {"binary", "continuous"}
TRAIT_ID_RE = re.compile(r"^[a-z][a-z0-9_]*$")


def sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(manifest_path, lock_path):
    with open(lock_path, encoding="utf-8") as handle:
        lock = json.load(handle)
    with open(manifest_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS.difference(columns)
        if missing:
            raise SystemExit(f"ERROR: analysis panel missing columns: {sorted(missing)}")
        rows = list(reader)

    if len(rows) != lock["trait_count"]:
        raise SystemExit(
            f"ERROR: expected {lock['trait_count']} traits, found {len(rows)}"
        )
    ids = [row["trait_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise SystemExit("ERROR: duplicate trait_id in analysis panel")
    invalid_ids = [trait_id for trait_id in ids if not TRAIT_ID_RE.fullmatch(trait_id)]
    if invalid_ids:
        raise SystemExit(f"ERROR: invalid trait_id values: {invalid_ids}")

    ordered_id_payload = "".join(f"{trait_id}\n" for trait_id in ids).encode()
    ids_sha256 = sha256_bytes(ordered_id_payload)
    if ids_sha256 != lock["ordered_trait_ids_sha256"]:
        raise SystemExit(
            "ERROR: ordered 45-trait identity lock changed. Review the scope change "
            "and update config/analysis_panel.lock.json explicitly if intentional."
        )

    versions = {row["panel_version"] for row in rows}
    if versions != {lock["panel_version"]}:
        raise SystemExit(
            f"ERROR: every row must use panel_version={lock['panel_version']}; got {sorted(versions)}"
        )
    counts = {}
    for row in rows:
        counts[row["domain"]] = counts.get(row["domain"], 0) + 1
        if row["type"] not in TRAIT_TYPES:
            raise SystemExit(f"ERROR: {row['trait_id']} has invalid type={row['type']}")
        if row["source_status"] not in SOURCE_STATES:
            raise SystemExit(
                f"ERROR: {row['trait_id']} has invalid source_status={row['source_status']}"
            )
        for column in REQUIRED_COLUMNS:
            if not str(row[column]).strip():
                raise SystemExit(f"ERROR: {row['trait_id']} has a blank {column}")
        if row["source_status"] == "SOURCE_VERIFIED" and row["source_id"] == "UNREGISTERED":
            raise SystemExit(
                f"ERROR: {row['trait_id']} claims SOURCE_VERIFIED without a source_id"
            )
    if counts != lock["domain_counts"]:
        raise SystemExit(
            f"ERROR: domain counts changed: expected {lock['domain_counts']}, got {counts}"
        )
    sleep_count = counts.get("sleep", 0)
    if sleep_count != lock["sleep_trait_count"] or len(rows) - sleep_count != lock["non_sleep_trait_count"]:
        raise SystemExit("ERROR: the panel is no longer exactly 12 sleep + 33 non-sleep traits")

    return {
        "panel_version": lock["panel_version"],
        "manifest": manifest_path,
        "manifest_sha256": file_sha256(manifest_path),
        "ordered_trait_ids_sha256": ids_sha256,
        "trait_count": len(rows),
        "sleep_trait_count": sleep_count,
        "non_sleep_trait_count": len(rows) - sleep_count,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST)
    parser.add_argument("--lock", default=DEFAULT_LOCK)
    parser.add_argument("--write-provenance")
    args = parser.parse_args()

    summary = validate(args.manifest, args.lock)
    if args.write_provenance:
        out_dir = os.path.dirname(args.write_provenance)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.write_provenance, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=summary.keys(), delimiter="\t")
            writer.writeheader()
            writer.writerow(summary)
    print(
        f"PANEL_LOCK_OK {summary['panel_version']} "
        f"{summary['trait_count']} traits "
        f"manifest_sha256={summary['manifest_sha256']}"
    )


if __name__ == "__main__":
    main()
