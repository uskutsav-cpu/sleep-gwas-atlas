#!/usr/bin/env python3
"""Pin exact S3 object versions for the locked extension without downloading bodies."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import ssl
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def versioned_url(url: str, version_id: str) -> str:
    parsed = urllib.parse.urlparse(url)
    query = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
    query = [(key, value) for key, value in query if key != "versionId"]
    query.append(("versionId", version_id))
    return urllib.parse.urlunparse(parsed._replace(query=urllib.parse.urlencode(query)))


def expected_etag(checksum: str) -> str | None:
    if checksum.startswith("etag:"):
        return checksum.removeprefix("etag:").lower()
    if checksum.startswith("md5:"):
        return checksum.removeprefix("md5:").lower()
    return None


def snapshot(row: dict[str, str], timeout: int) -> dict[str, str]:
    url = row["url"]
    result = dict(row)
    result.update(
        {
            "http_status": "NA",
            "observed_content_length": "NA",
            "observed_etag": "NA",
            "last_modified": "NA",
            "s3_version_id": "NA",
            "versioned_url": "NA",
            "head_verification_status": "ERROR",
            "head_error": "NA",
        }
    )
    try:
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme != "https" or parsed.hostname != "pan-ukb-us-east-1.s3.amazonaws.com":
            raise ValueError(f"unapproved remote host: {url}")
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "sleep-gwas-atlas-extension/1.0"},
            method="HEAD",
        )
        try:
            import certifi

            context = ssl.create_default_context(cafile=certifi.where())
        except ImportError:
            context = ssl.create_default_context()
        with urllib.request.urlopen(  # nosec B310: host locked above
            request, timeout=timeout, context=context
        ) as response:
            headers = {key.lower(): value for key, value in response.headers.items()}
            status = str(response.status)
        observed_length = headers.get("content-length", "NA")
        observed_etag = headers.get("etag", "NA").strip('"').lower()
        version_id = headers.get("x-amz-version-id", "NA")
        if status != "200":
            raise ValueError(f"HEAD returned HTTP {status}")
        if observed_length == "NA" or int(observed_length) != int(row["expected_size_bytes"]):
            raise ValueError(
                f"content length mismatch: observed {observed_length}, expected {row['expected_size_bytes']}"
            )
        if version_id in {"", "NA", "null"}:
            raise ValueError("S3 version ID is absent")
        checksum_etag = expected_etag(row["expected_checksum"])
        if row["expected_checksum"].startswith("etag:") and checksum_etag != observed_etag:
            raise ValueError(
                f"ETag mismatch: observed {observed_etag}, expected {checksum_etag}"
            )
        if row["expected_checksum"].startswith("md5:") and "-" not in observed_etag and checksum_etag != observed_etag:
            raise ValueError(
                f"single-part ETag differs from manifest MD5: {observed_etag} != {checksum_etag}"
            )
        checksum_status = (
            "HEAD_SIZE_VERSION_AND_ETAG_PASS"
            if row["expected_checksum"].startswith("etag:") or "-" not in observed_etag
            else "HEAD_SIZE_VERSION_PASS_BODY_MD5_PENDING_STREAM"
        )
        result.update(
            {
                "http_status": status,
                "observed_content_length": observed_length,
                "observed_etag": observed_etag,
                "last_modified": headers.get("last-modified", "NA"),
                "s3_version_id": version_id,
                "versioned_url": versioned_url(url, version_id),
                "head_verification_status": checksum_status,
            }
        )
    except Exception as error:  # preserve all source failures as rows
        result["head_error"] = f"{type(error).__name__}: {error}"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plan",
        type=Path,
        default=Path("discovery_extension/results/acquisition_plan.tsv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/remote_object_snapshot.tsv"),
    )
    parser.add_argument(
        "--provenance-out",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/remote_object_snapshot.json"),
    )
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--timeout-seconds", type=int, default=60)
    args = parser.parse_args()

    plan = read_tsv(args.plan)
    requests: list[dict[str, str]] = []
    for row in plan:
        role = (
            "shared_variant_reference"
            if row["extension_trait_id"] == "SHARED_VARIANT_REFERENCE"
            else "phenotype_sumstats"
        )
        requests.append(
            {
                "extension_trait_id": row["extension_trait_id"],
                "object_role": role,
                "filename": row["source_filename"],
                "url": row["source_url"],
                "expected_checksum": row["checksum"],
                "expected_size_bytes": row["source_file_size_bytes"],
            }
        )
        if row["source_tabix_url"] not in {"", "NA"}:
            requests.append(
                {
                    "extension_trait_id": row["extension_trait_id"],
                    "object_role": "phenotype_tabix",
                    "filename": row["source_tabix_filename"],
                    "url": row["source_tabix_url"],
                    "expected_checksum": row["source_tabix_checksum"],
                    "expected_size_bytes": row["source_tabix_size_bytes"],
                }
            )
    identity = [(row["extension_trait_id"], row["object_role"], row["filename"]) for row in requests]
    if len(identity) != len(set(identity)):
        raise SystemExit("ERROR: duplicate remote-object identity in acquisition plan")

    results: list[dict[str, str] | None] = [None] * len(requests)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_index = {
            executor.submit(snapshot, row, args.timeout_seconds): index
            for index, row in enumerate(requests)
        }
        for future in as_completed(future_to_index):
            results[future_to_index[future]] = future.result()
    completed = [row for row in results if row is not None]

    fields = [
        "extension_trait_id", "object_role", "filename", "url", "expected_checksum",
        "expected_size_bytes", "http_status", "observed_content_length", "observed_etag",
        "last_modified", "s3_version_id", "versioned_url", "head_verification_status", "head_error",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(completed)

    failures = [row for row in completed if row["head_verification_status"] == "ERROR"]
    body_md5_pending = sum(
        row["head_verification_status"] == "HEAD_SIZE_VERSION_PASS_BODY_MD5_PENDING_STREAM"
        for row in completed
    )
    provenance = {
        "schema_version": "1.0.0",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "acquisition_plan": str(args.plan),
        "acquisition_plan_sha256": sha256(args.plan),
        "object_count": len(completed),
        "head_pass_count": len(completed) - len(failures),
        "body_md5_pending_stream_count": body_md5_pending,
        "failure_count": len(failures),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "interpretation": (
            "HEAD pins exact S3 versions and validates byte counts. Multipart object MD5 values "
            "remain pending until every compressed byte is consumed by the streaming harmonizer."
        ),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    if failures:
        sample = "; ".join(
            f"{row['extension_trait_id']}:{row['object_role']}={row['head_error']}"
            for row in failures[:5]
        )
        raise SystemExit(f"ERROR: {len(failures)} remote object HEAD checks failed: {sample}")
    print(
        f"REMOTE_SOURCE_SNAPSHOT_OK objects={len(completed)} "
        f"body_md5_pending_stream={body_md5_pending} sha256={provenance['output_sha256']}"
    )


if __name__ == "__main__":
    main()
