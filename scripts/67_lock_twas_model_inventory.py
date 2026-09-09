#!/usr/bin/env python3
"""Snapshot and lock the complete result-free Box inventory for phi-enabled GTEx-v8 models."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path


FIELDS = [
    "file_id", "filename", "context", "file_role", "bytes", "download_url",
    "results_accessed_before_lock",
]
FILE_ID = re.compile(r"^f_([0-9]+)$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def table_text(rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


class BoxListingParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.files: list[dict[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "li":
            return
        values = {key: value or "" for key, value in attrs}
        classes = set(values.get("class", "").split())
        if "tbl-list-item" not in classes or values.get("data-item-type") != "file":
            return
        self.files.append({
            "file_id": values.get("data-id", ""), "filename": values.get("data-name", ""),
            "bytes": values.get("data-size", ""), "updated_epoch": values.get("data-lastmodified", ""),
            "parent_folder_id": values.get("data-parent-folder-id", ""),
        })


def fetch(url: str) -> bytes:
    result = subprocess.run(
        ["curl", "-fsSL", "--retry", "3", "--user-agent", "sleep-gwas-atlas/1.0 inventory-lock", url],
        capture_output=True,
    )
    if result.returncode:
        fail(f"Box inventory request failed: {url}: {result.stderr.decode('utf-8', errors='replace').strip()}")
    return result.stdout


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument(
        "--relock-existing-snapshots", action="store_true",
        help="rebuild only the result-free inventory/lock from the existing exact HTML snapshots",
    )
    args = parser.parse_args()
    if not args.execute:
        fail("Box metadata access requires explicit --execute after reviewing the source contract")
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    source = policy["twas"]["phi_model_source"]
    inventory_path = root / source["inventory_path"]
    lock_path = root / source["inventory_lock_path"]
    snapshot_dir = root / source["source_snapshot_dir"]
    existing_lock = None
    if args.relock_existing_snapshots:
        if not inventory_path.is_file() or not lock_path.is_file() or not snapshot_dir.is_dir():
            fail("relocking requires the existing inventory, lock, and exact source snapshots")
        existing_lock = json.loads(lock_path.read_text(encoding="utf-8"))
        if existing_lock.get("results_accessed_before_lock") is not False:
            fail("existing inventory is not a result-free lock")
    elif inventory_path.exists() or lock_path.exists() or snapshot_dir.exists():
        fail("TWAS phi-model inventory outputs already exist; refusing overwrite")
    work_root = root / "work/twas_phi_inventory"
    work_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="box_inventory.", dir=work_root))
    try:
        staged_snapshots = staging / "source_snapshots"; staged_snapshots.mkdir()
        observed: list[dict[str, str]] = []
        page_hashes: dict[str, str] = {}
        page_count = (int(source["expected_file_count"]) + 19) // 20
        for page in range(page_count):
            url = (
                f"{source['box_url']}?sortColumn=name&sortDirection=asc&pageSize=20&pageNumber={page}"
            )
            existing_page = snapshot_dir / f"page_{page + 1:02d}.html"
            payload = existing_page.read_bytes() if args.relock_existing_snapshots else fetch(url)
            page_path = staged_snapshots / f"page_{page + 1:02d}.html"
            page_path.write_bytes(payload)
            listing = BoxListingParser(); listing.feed(payload.decode("utf-8"))
            observed.extend(listing.files)
            page_hashes[page_path.name] = sha256(page_path)
        by_id = {row["file_id"]: row for row in observed}
        if len(observed) != len(by_id) or len(observed) != int(source["expected_file_count"]):
            fail(f"Box inventory is duplicated or incomplete: rows={len(observed)} unique={len(by_id)}")
        prefix = source["expected_db_prefix"]
        rows: list[dict[str, object]] = []
        pairs: dict[str, set[str]] = {}
        total_bytes = 0
        for record in observed:
            match = FILE_ID.fullmatch(record["file_id"])
            try:
                size = int(record["bytes"])
            except ValueError as exc:
                fail(f"Box inventory has a nonnumeric size: {record}")
                raise AssertionError from exc
            if match is None or size <= 0 or record["parent_folder_id"] != source["box_folder_id"]:
                fail(f"Box inventory has an invalid identity: {record}")
            filename = record["filename"]
            if filename.startswith(prefix) and filename.endswith(".db"):
                context, role = filename[len(prefix):-3], "MODEL_DB"
            elif filename.startswith(prefix) and filename.endswith(".txt.gz"):
                context, role = filename[len(prefix):-7], "COVARIANCE"
            else:
                fail(f"unexpected file in phi-model release: {filename}")
            if not context or any(character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-" for character in context):
                fail(f"unsafe or empty phi-model context: {context}")
            pairs.setdefault(context, set()).add(role)
            file_number = match.group(1)
            download_url = (
                "https://uchicago.app.box.com/index.php?rm=box_download_shared_file"
                f"&shared_name={source['box_shared_name']}&file_id=f_{file_number}"
            )
            rows.append({
                "file_id": record["file_id"], "filename": filename, "context": context,
                "file_role": role, "bytes": size, "download_url": download_url,
                "results_accessed_before_lock": "NO",
            })
            total_bytes += size
        if len(pairs) != int(source["expected_context_count"]) or any(roles != {"MODEL_DB", "COVARIANCE"} for roles in pairs.values()):
            fail("Box inventory is not exactly one model database and covariance per context")
        if total_bytes != int(source["expected_total_bytes"]):
            fail(f"Box inventory byte total drifted: expected {source['expected_total_bytes']}, got {total_bytes}")
        rows.sort(key=lambda row: (str(row["context"]), str(row["file_role"])))
        staged_inventory = staging / "inventory.tsv"
        staged_inventory.write_text(table_text(rows), encoding="utf-8")
        lock = {
            "schema_version": "atlas-v1.0-twas-phi-model-inventory.1",
            "locked_utc": (existing_lock or {}).get(
                "locked_utc", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            ),
            "results_accessed_before_lock": False, "source_release": source["release_name"],
            "source_release_updated_utc_date": source["release_updated_utc_date"],
            "box_shared_name": source["box_shared_name"], "box_folder_id": source["box_folder_id"],
            "file_count": len(rows), "context_count": len(pairs), "total_bytes": total_bytes,
            "file_ids_in_locked_order": [row["file_id"] for row in rows],
            "inventory_sha256": sha256(staged_inventory), "policy_sha256": sha256(policy_path),
            "source_snapshot_sha256": page_hashes, "claim_limit": policy["claim_limit"],
        }
        staged_lock = staging / "inventory.lock.json"
        staged_lock.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if args.relock_existing_snapshots and page_hashes != existing_lock.get("source_snapshot_sha256"):
            fail("existing source snapshots differ from their original result-free lock")
        inventory_path.parent.mkdir(parents=True, exist_ok=True); snapshot_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staged_inventory, inventory_path); os.replace(staged_lock, lock_path)
        if args.relock_existing_snapshots:
            shutil.rmtree(staged_snapshots)
        else:
            os.replace(staged_snapshots, snapshot_dir)
        staging.rmdir()
        print(f"TWAS_PHI_MODEL_INVENTORY_OK contexts={len(pairs)} files={len(rows)} bytes={total_bytes} result_free=true")
        return 0
    finally:
        if staging.exists(): shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
