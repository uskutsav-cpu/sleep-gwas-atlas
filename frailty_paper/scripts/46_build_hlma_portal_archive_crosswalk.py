#!/usr/bin/env python3
"""Build a link-metadata crosswalk from the official HLMA portal API snapshot."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from urllib.parse import urlparse

EXPECTED_SECTION_COUNTS = {
    "snRNA-seq count matrix": 10,
    "snATAC-seq fragments": 46,
    "snATAC-seq peak matrix": 6,
    "scRNA-seq count matrix": 1,
    "scRNA-seq h5ad": 6,
}
ARCHIVE_RE = re.compile(r"/(OMIX\d+)/(OMIX\d+-\d+)(?:\.|$)")
RELEASE_ROWS = {
    "OMIX004305-01": ("ATAC-seq fragments for OM1 and OM2", "5", "9.46 GB"),
    "OMIX004305-02": ("ATAC-seq fragments for OM5 and OM9", "6", "7.91 GB"),
    "OMIX004305-03": ("ATAC-seq fragments for YM2, YM3, YM4", "10", "12.07 GB"),
    "OMIX004305-04": ("ATAC-seq fragments for OM3, OM5, OM6, OM7, OM9", "11", "14.83 GB"),
    "OMIX004305-05": ("snATAC-seq rds object", "7", "2.34 GB"),
    "OMIX006021-01": ("ATAC-seq fragments for P26", "4", "4.39 GB"),
    "OMIX006021-02": ("ATAC-seq fragments for P17", "1", "969.01 MB"),
    "OMIX006021-03": ("ATAC-seq fragments for P29", "1", "1.46 GB"),
    "OMIX006021-04": ("ATAC-seq fragments for P21", "1", "1.37 GB"),
    "OMIX006021-05": ("ATAC-seq fragments for P5", "4", "3.59 GB"),
    "OMIX006021-06": ("ATAC-seq fragments for P23", "3", "2.93 GB"),
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def archive_id(row: dict) -> tuple[str, str]:
    url = row.get("file_path", "")
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.netloc != "download.cncb.ac.cn":
        raise ValueError(f"unexpected official archive URL: {url!r}")
    match = ARCHIVE_RE.search(parsed.path)
    if not match or match.group(1) not in match.group(2):
        raise ValueError(f"cannot derive OMIX file ID from portal path: {url!r}")
    return match.group(2), url


def canonical_note(notes: str, omix_id: str) -> str:
    """Replace prior generated suffixes with one stable, explicit note."""
    prefix = re.sub(r"\s+Official portal API snapshot\b.*$", "", notes)
    return prefix + (
        f" Official portal API snapshot reports archive path {omix_id}; link metadata only, "
        "no archive bytes or member identity verified. Per-file CNP mapping, donor crosswalk, "
        "genome build, and object-specific reuse terms remain unresolved."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo.resolve()
    manifest_dir = root / "frailty_paper/manifests"
    snapshot = manifest_dir / "hlma_portal_download_metadata_2026-09-24.json"
    inventory_path = manifest_dir / "hlma_processed_file_inventory.tsv"
    crosswalk_path = manifest_dir / "hlma_fragment_archive_crosswalk.tsv"

    response = json.loads(snapshot.read_text(encoding="utf-8"))
    sections = response.get("data")
    if not isinstance(sections, list):
        raise ValueError("portal response has no data section list")
    found_counts = {section.get("label"): len(section.get("data", [])) for section in sections}
    if found_counts != EXPECTED_SECTION_COUNTS:
        raise ValueError(f"portal section counts changed: {found_counts!r}")

    all_rows: list[dict] = []
    by_name_size: dict[tuple[str, str], dict] = {}
    section_rows: dict[str, list[dict]] = {}
    for section in sections:
        label = section["label"]
        section_rows[label] = []
        for row in section["data"]:
            if row.get("public_type") != "public":
                raise ValueError(f"non-public portal listing row: {row!r}")
            omix_id, url = archive_id(row)
            item = {
                "section": label,
                "portal_row_id": str(row.get("id", "")),
                "file_name": row.get("file_name", ""),
                "portal_size": row.get("size", ""),
                "omix_archive_id": omix_id,
                "download_url": url,
            }
            if not item["file_name"] or not item["portal_size"]:
                raise ValueError(f"incomplete portal listing row: {row!r}")
            all_rows.append(item)
            section_rows[label].append(item)
            key = (item["file_name"], item["portal_size"])
            if key in by_name_size:
                raise ValueError(f"duplicate portal filename/size key: {key!r}")
            by_name_size[key] = item

    fields, inventory = read_tsv(inventory_path)
    new_inventory_fields = ["portal_archive_id", "portal_archive_url", "portal_link_status"]
    for field in new_inventory_fields:
        if field not in fields:
            fields.append(field)
    required = {"resource_id", "file_name_as_listed", "portal_display_size", "per_file_accession_mapping", "notes"}
    if not required.issubset(fields):
        raise ValueError(f"inventory lacks fields: {sorted(required - set(fields))}")
    targeted = 0
    for item in inventory:
        key = (item["file_name_as_listed"], item["portal_display_size"])
        match = by_name_size.get(key)
        if match is None:
            continue
        targeted += 1
        omix = match["omix_archive_id"]
        item["portal_archive_id"] = omix
        item["portal_archive_url"] = match["download_url"]
        item["portal_link_status"] = "PORTAL_LINK_TARGET_METADATA_ONLY"
        item["per_file_accession_mapping"] = "UNKNOWN; portal archive link does not establish a per-file CNP accession"
        if item["resource_id"] == "HLMA_snRNA_count_matrix_06":
            item["portal_link_status"] = "PORTAL_LINK_TARGET_RELEASE_LABEL_CONFLICT"
        elif item["resource_id"] == "HLMA_snRNA_count_matrix_10":
            item["portal_link_status"] = "PORTAL_LINK_TARGET_RELEASE_TITLE_SIZE_CONFLICT"
        elif item["resource_id"] == "HLMA_scsnRNA_h5ad_03":
            item["portal_link_status"] = "PORTAL_LINK_TARGET_RELEASE_TITLE_SIZE_CONFLICT"
        elif omix == "OMIX004308-01":
            item["portal_link_status"] = "PORTAL_CONTAINER_LINK_GROUPS_MULTIPLE_OBJECT_ROWS"
        item["notes"] = canonical_note(item["notes"], omix)
    if targeted != 23:
        raise ValueError(f"expected 23 targeted objects matched by exact portal filename and size, found {targeted}")

    with inventory_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(inventory)

    fragment_rows = section_rows["snATAC-seq fragments"]
    archive_row_counts: dict[str, int] = {}
    for item in fragment_rows:
        archive_row_counts[item["omix_archive_id"]] = archive_row_counts.get(item["omix_archive_id"], 0) + 1
    fragment_fields = ["portal_row_id", "file_name", "portal_size", "omix_archive_id", "download_url", "official_release_title", "official_release_sample_count", "official_release_size", "portal_rows_sharing_archive_path", "sample_label_comparison", "download_status", "notes"]
    with crosswalk_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fragment_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in fragment_rows:
            omix = item["omix_archive_id"]
            title, sample_count, release_size = RELEASE_ROWS[omix]
            sample = re.match(r"^([A-Z]{1,2}\d+)(?:_|\.)", item["file_name"])
            sample_label = sample.group(1) if sample else "UNKNOWN"
            label_status = "LABEL_MATCHES_RELEASE_TITLE" if sample_label != "UNKNOWN" and sample_label in title else "SAMPLE_LABEL_CONFLICT_WITH_RELEASE_TITLE"
            writer.writerow({
                "portal_row_id": item["portal_row_id"],
                "file_name": item["file_name"],
                "portal_size": item["portal_size"],
                "omix_archive_id": item["omix_archive_id"],
                "download_url": item["download_url"],
                "official_release_title": title,
                "official_release_sample_count": sample_count,
                "official_release_size": release_size,
                "portal_rows_sharing_archive_path": archive_row_counts[omix],
                "sample_label_comparison": label_status + "; portal ID is link metadata, not byte/member validation",
                "download_status": "LISTED_NOT_DOWNLOADED",
                "notes": "Official release title/count/size were read from the linked NGDC OMIX release page. Different counts between portal rows and release sample counts are retained as metadata; they do not establish missing or duplicate files. No archive bytes were downloaded.",
            })

    print(f"HLMA_PORTAL_CROSSWALK_OK targeted_objects={targeted} fragments={len(fragment_rows)} api_rows={len(all_rows)} unique_archive_ids={len({r['omix_archive_id'] for r in all_rows})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
