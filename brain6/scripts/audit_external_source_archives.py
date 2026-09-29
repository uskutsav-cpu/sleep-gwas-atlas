#!/usr/bin/env python3
"""Audit cached raw GWAS archives against the frozen Brain6 source ledger."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPECTED_TRAITS = {"insomnia", "longsleep", "adhd", "mdd", "scz", "bipolar", "parkinson"}
FIELDS = ("trait", "study", "source_URL", "expected_source_archive_bytes",
          "expected_source_archive_SHA256", "archive_path", "content_transform",
          "container_bytes", "container_SHA256", "observed_source_bytes",
          "observed_source_SHA256", "status")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def audit_archives(manifest: Path, candidates: Path, archive_root: Path) -> list[dict[str, str]]:
    manifest_rows = {row["trait"]: row for row in read_tsv(manifest)
                     if row.get("trait") in EXPECTED_TRAITS}
    candidate_rows = read_tsv(candidates)
    by_trait = {row["trait"]: row for row in candidate_rows}
    if len(by_trait) != len(candidate_rows) or set(by_trait) != EXPECTED_TRAITS:
        raise ValueError("Archive candidate map must contain exactly one row for each of the seven locked traits")
    if set(manifest_rows) != EXPECTED_TRAITS:
        raise ValueError("Frozen GWAS source manifest does not cover the seven selected traits")

    records = []
    for trait in sorted(EXPECTED_TRAITS):
        source = manifest_rows[trait]
        relpath = Path(by_trait[trait]["archive_path"])
        transform = by_trait[trait]["content_transform"]
        if transform not in {"identity", "gzip_decompress"}:
            raise ValueError(f"Unsupported source archive content transform for {trait}: {transform}")
        if relpath.is_absolute() or ".." in relpath.parts:
            raise ValueError(f"Archive candidate must stay under archive root: {relpath}")
        path = (archive_root / relpath).resolve()
        try:
            path.relative_to(archive_root.resolve())
        except ValueError as exc:
            raise ValueError(f"Archive path escapes archive root: {relpath}") from exc
        expected_size = source["expected_source_archive_bytes"]
        expected_sha = source["expected_source_archive_SHA256"]
        if not path.is_file():
            container_size, container_sha = "", ""
            observed_size, observed_sha, status = "", "", "MISSING"
        else:
            container_size = str(path.stat().st_size)
            container_sha = sha256_file(path)
            opener = gzip.open if transform == "gzip_decompress" else open
            digest = hashlib.sha256()
            observed_bytes = 0
            with opener(path, "rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    observed_bytes += len(block)
                    digest.update(block)
            observed_size = str(observed_bytes)
            observed_sha = digest.hexdigest()
            verified_content = observed_size == expected_size and observed_sha == expected_sha
            status = (("VERIFIED_EXACT_ARCHIVE" if transform == "identity" else "VERIFIED_DECOMPRESSED_SOURCE")
                      if verified_content else "MISMATCH")
        records.append({
            "trait": trait,
            "study": source["study"],
            "source_URL": source["source_URL"],
            "expected_source_archive_bytes": expected_size,
            "expected_source_archive_SHA256": expected_sha,
            "archive_path": str(path),
            "content_transform": transform,
            "container_bytes": container_size,
            "container_SHA256": container_sha,
            "observed_source_bytes": observed_size,
            "observed_source_SHA256": observed_sha,
            "status": status,
        })
    return records


def write_tsv(path: Path, records: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "brain6/manifests/gwas_master.tsv")
    parser.add_argument("--candidates", type=Path,
                        default=ROOT / "brain6/manifests/source_archive_candidates.tsv")
    parser.add_argument("--archive-root", type=Path,
                        default=Path("/Volumes/Extreme SSD"))
    parser.add_argument("--output", type=Path,
                        default=ROOT / "brain6/qc/source_archive_audit.tsv")
    args = parser.parse_args()
    records = audit_archives(args.manifest, args.candidates, args.archive_root)
    write_tsv(args.output, records)
    summary = {status: sum(row["status"] == status for row in records)
               for status in ("VERIFIED_EXACT_ARCHIVE", "VERIFIED_DECOMPRESSED_SOURCE", "MISSING", "MISMATCH")}
    print(json.dumps({"audit": "PASS_NO_HASH_MISMATCH" if not summary["MISMATCH"] else "FAIL_HASH_MISMATCH",
                      "counts": summary, "output": str(args.output.resolve())}, indent=2))
    if summary["MISMATCH"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
