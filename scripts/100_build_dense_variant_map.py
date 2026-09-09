#!/usr/bin/env python3
"""Acquire and deterministically build the genome-wide GRCh37 identity map.

The large transfer is opt-in.  The source is a result-free Pan-UK Biobank
variant manifest; only chromosome, position, alleles, and rsID are retained.
Association statistics, allele frequencies, INFO, and Pan-UKBB QC labels are
never used to harmonize another GWAS.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
import urllib.request


SHA256 = re.compile(r"[0-9a-f]{64}$")
MD5 = re.compile(r"[0-9a-f]{32}$")
RSID = re.compile(r"rs[0-9]+$", re.IGNORECASE)
VALID_ALLELES = {"A", "C", "G", "T"}
MAP_HEADER = ["SNP", "CHR", "BP", "A1", "A2"]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def file_hashes(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_policy(root: Path, relative: str) -> tuple[Path, dict[str, object]]:
    path = root / relative
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"dense variant-map policy is unreadable: {exc}")
    required = {
        "schema_version", "analysis_panel", "resource", "source_url", "source_path",
        "source_bytes", "source_md5", "source_s3_version_id",
        "source_header_sha256", "required_columns", "source_build", "map_path",
        "map_provenance_path", "map_schema_version", "map_scope",
        "minimum_source_rows", "minimum_mapped_rows",
        "minimum_free_bytes_before_download", "minimum_free_bytes_before_build",
        "requires_large_download_acknowledgement", "identity_rule", "role_limit",
        "builder_path", "builder_sha256",
    }
    if not required.issubset(policy):
        fail(f"dense variant-map policy lacks: {sorted(required - set(policy))}")
    if (
        policy["schema_version"] != "sleep-atlas-dense-variant-map.1"
        or policy["analysis_panel"] != "atlas-v1.0"
        or policy["source_build"] != "GRCh37"
        or policy["map_schema_version"] != "atlas.grch37-variant-map.v2"
        or policy["map_scope"] != "GENOME_WIDE_IMPUTED_VARIANT_IDENTITY"
        or not str(policy["source_url"]).startswith("https://")
        or not MD5.fullmatch(str(policy["source_md5"]))
        or not SHA256.fullmatch(str(policy["source_header_sha256"]))
        or not SHA256.fullmatch(str(policy["builder_sha256"]))
        or policy["required_columns"] != ["chrom", "pos", "ref", "alt", "rsid"]
        or not policy["requires_large_download_acknowledgement"]
    ):
        fail("dense variant-map policy differs from the reviewed result-free contract")
    for key in ("source_path", "map_path", "map_provenance_path"):
        relative_path = Path(str(policy[key]))
        if relative_path.is_absolute() or ".." in relative_path.parts:
            fail(f"unsafe path in dense variant-map policy: {key}")
    builder_path = Path(str(policy["builder_path"]))
    if (
        builder_path.is_absolute() or ".." in builder_path.parts
        or not (root / builder_path).is_file()
        or sha256(root / builder_path) != policy["builder_sha256"]
    ):
        fail("dense variant-map builder differs from its pre-result policy hash")
    return path, policy


def remote_identity(policy: dict[str, object]) -> dict[str, str | int]:
    request = urllib.request.Request(str(policy["source_url"]), method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        headers = response.headers
        observed = {
            "bytes": int(headers.get("Content-Length", "0")),
            "etag": headers.get("ETag", "").strip('"').lower(),
            "version_id": headers.get("x-amz-version-id", ""),
            "last_modified": headers.get("Last-Modified", ""),
        }
    if (
        observed["bytes"] != int(policy["source_bytes"])
        or observed["etag"] != policy["source_md5"]
        or observed["version_id"] != policy["source_s3_version_id"]
    ):
        fail(f"remote Pan-UKBB variant manifest differs from policy: {observed}")
    return observed


def download_source(root: Path, policy: dict[str, object]) -> None:
    source = root / str(policy["source_path"])
    source.parent.mkdir(parents=True, exist_ok=True)
    if source.is_file():
        verify_source(source, policy)
        print(f"Reusing verified dense-map source: {source}")
        return
    free = shutil.disk_usage(source.parent).free
    if free < int(policy["minimum_free_bytes_before_download"]):
        fail(
            f"dense-map source download requires {policy['minimum_free_bytes_before_download']} "
            f"free bytes; observed {free}"
        )
    remote_identity(policy)
    partial = Path(f"{source}.part")
    start = partial.stat().st_size if partial.is_file() else 0
    expected = int(policy["source_bytes"])
    if start > expected:
        fail(f"partial dense-map source is larger than its exact object: {partial}")
    if start == expected:
        observed_md5, _ = file_hashes(partial)
        if observed_md5 != policy["source_md5"]:
            fail("complete dense-map partial differs from the locked S3 object")
        os.replace(partial, source)
        source.chmod(0o644)
        verify_source(source, policy)
        print(f"Promoted verified complete partial: {source}")
        return
    headers = {}
    mode = "wb"
    if start:
        headers["Range"] = f"bytes={start}-"
        mode = "ab"
    request = urllib.request.Request(str(policy["source_url"]), headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        if start and response.status != 206:
            fail("server did not honor the resume range; partial file was preserved")
        if response.headers.get("x-amz-version-id", "") != policy["source_s3_version_id"]:
            fail("downloaded object version differs from the locked S3 version")
        with partial.open(mode) as handle:
            shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)
    if partial.stat().st_size != expected:
        fail(
            f"partial dense-map source has {partial.stat().st_size} bytes; expected {expected}; "
            "rerun to resume"
        )
    observed_md5, _ = file_hashes(partial)
    if observed_md5 != policy["source_md5"]:
        fail("downloaded dense-map source MD5 differs from the locked S3 ETag")
    os.replace(partial, source)
    source.chmod(0o644)
    verify_source(source, policy)
    print(f"Downloaded and verified {source.stat().st_size:,} bytes: {source}")


def verify_source(source: Path, policy: dict[str, object]) -> tuple[str, str]:
    if not source.is_file() or source.stat().st_size != int(policy["source_bytes"]):
        fail("dense-map source is absent or has the wrong byte count")
    observed_md5, observed_sha256 = file_hashes(source)
    if observed_md5 != policy["source_md5"]:
        fail("dense-map source MD5 differs from the locked S3 object identity")
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        header_line = handle.readline()
    if hashlib.sha256(header_line.encode("utf-8")).hexdigest() != policy["source_header_sha256"]:
        fail("dense-map source header differs from the locked schema")
    header = header_line.rstrip("\r\n").split("\t")
    if any(column not in header for column in policy["required_columns"]):
        fail("dense-map source lacks a required identity column")
    return observed_md5, observed_sha256


def write_map(database: sqlite3.Connection, path: Path) -> int:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    rows_written = 0
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                    writer.writerow(MAP_HEADER)
                    query = """
                        SELECT v.snp, v.chromosome, v.position, v.a1, v.a2
                        FROM variants AS v
                        JOIN (
                            SELECT snp FROM variants GROUP BY snp HAVING COUNT(*) = 1
                        ) AS unique_ids ON unique_ids.snp = v.snp
                        ORDER BY v.chromosome, v.position, v.snp, v.a1, v.a2
                    """
                    for row in database.execute(query):
                        writer.writerow(row)
                        rows_written += 1
        os.replace(temporary, path)
        path.chmod(0o644)
    finally:
        if temporary.exists():
            temporary.unlink()
    return rows_written


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8",
        )
        os.replace(temporary, path)
        path.chmod(0o644)
    finally:
        if temporary.exists():
            temporary.unlink()


def build_map(root: Path, policy_path: Path, policy: dict[str, object]) -> None:
    source = root / str(policy["source_path"])
    output = root / str(policy["map_path"])
    provenance_path = root / str(policy["map_provenance_path"])
    if output.exists() or provenance_path.exists():
        fail("dense variant map or provenance already exists; validate it instead of overwriting")
    output.parent.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(output.parent).free
    if free < int(policy["minimum_free_bytes_before_build"]):
        fail(
            f"dense-map build requires {policy['minimum_free_bytes_before_build']} free bytes; "
            f"observed {free}"
        )
    source_md5, source_sha256 = verify_source(source, policy)
    source_rows = 0
    eligible_rows = 0
    excluded = {
        "wrong_width": 0, "invalid_coordinate": 0, "non_autosomal": 0,
        "non_snp_or_equal_alleles": 0, "missing_or_non_rsid": 0,
        "exact_duplicate": 0,
    }
    with tempfile.TemporaryDirectory(prefix="dense-variant-map-", dir=output.parent) as temporary:
        database_path = Path(temporary) / "identity.sqlite"
        database = sqlite3.connect(database_path)
        try:
            database.execute("PRAGMA journal_mode=OFF")
            database.execute("PRAGMA synchronous=OFF")
            database.execute("PRAGMA temp_store=FILE")
            database.execute("PRAGMA cache_size=-262144")
            database.execute(
                "CREATE TABLE variants ("
                "snp TEXT NOT NULL, chromosome INTEGER NOT NULL, position INTEGER NOT NULL, "
                "a1 TEXT NOT NULL, a2 TEXT NOT NULL, "
                "UNIQUE(snp, chromosome, position, a1, a2))"
            )
            with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
                rows = csv.reader(handle, delimiter="\t")
                header = next(rows, None)
                if header is None:
                    fail("dense-map source is empty")
                indices = {column: header.index(column) for column in policy["required_columns"]}
                width = len(header)
                batch: list[tuple[str, int, int, str, str]] = []
                for values in rows:
                    source_rows += 1
                    if len(values) != width:
                        excluded["wrong_width"] += 1
                        continue
                    try:
                        chromosome = int(values[indices["chrom"]])
                        position = int(values[indices["pos"]])
                    except ValueError:
                        excluded["invalid_coordinate"] += 1
                        continue
                    if chromosome < 1 or chromosome > 22:
                        excluded["non_autosomal"] += 1
                        continue
                    if position <= 0:
                        excluded["invalid_coordinate"] += 1
                        continue
                    a1 = values[indices["ref"]].strip().upper()
                    a2 = values[indices["alt"]].strip().upper()
                    if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES or a1 == a2:
                        excluded["non_snp_or_equal_alleles"] += 1
                        continue
                    snp = values[indices["rsid"]].strip().lower()
                    if not RSID.fullmatch(snp):
                        excluded["missing_or_non_rsid"] += 1
                        continue
                    batch.append((snp, chromosome, position, a1, a2))
                    eligible_rows += 1
                    if len(batch) == 100_000:
                        before = database.total_changes
                        database.executemany("INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?)", batch)
                        excluded["exact_duplicate"] += len(batch) - (database.total_changes - before)
                        database.commit()
                        batch.clear()
                if batch:
                    before = database.total_changes
                    database.executemany("INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?)", batch)
                    excluded["exact_duplicate"] += len(batch) - (database.total_changes - before)
                    database.commit()
            if source_rows < int(policy["minimum_source_rows"]):
                fail(
                    f"dense-map source has {source_rows:,} rows; policy requires at least "
                    f"{policy['minimum_source_rows']:,}"
                )
            database.execute("CREATE INDEX variants_snp ON variants(snp)")
            database.execute(
                "CREATE INDEX variants_coordinate ON variants(chromosome, position, a1, a2)"
            )
            database.commit()
            distinct_rows = database.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
            ambiguous_ids, ambiguous_rows = database.execute(
                "SELECT COUNT(*), COALESCE(SUM(n),0) FROM "
                "(SELECT COUNT(*) AS n FROM variants GROUP BY snp HAVING COUNT(*) > 1)"
            ).fetchone()
            mapped_rows = write_map(database, output)
        finally:
            database.close()
    if mapped_rows < int(policy["minimum_mapped_rows"]):
        output.unlink()
        fail(
            f"dense map has {mapped_rows:,} rows; policy requires at least "
            f"{policy['minimum_mapped_rows']:,}"
        )
    provenance = {
        "schema_version": policy["map_schema_version"],
        "map_scope": policy["map_scope"],
        "resource": policy["resource"],
        "genome_build": "GRCh37/hg19",
        "analysis_panel": policy["analysis_panel"],
        "policy_path": str(policy_path.relative_to(root)),
        "policy_sha256": sha256(policy_path),
        "builder_path": policy["builder_path"],
        "builder_sha256": policy["builder_sha256"],
        "source": {
            "path": str(source.relative_to(root)), "url": policy["source_url"],
            "bytes": source.stat().st_size, "md5": source_md5,
            "sha256": source_sha256, "s3_version_id": policy["source_s3_version_id"],
        },
        "source_rows": source_rows,
        "eligible_source_rows": eligible_rows,
        "distinct_identity_rows": distinct_rows,
        "ambiguous_rsid_count": ambiguous_ids,
        "ambiguous_rsid_rows_excluded": ambiguous_rows,
        "exclusions": excluded,
        "mapped_rows": mapped_rows,
        "identity_rule": policy["identity_rule"],
        "role_limit": policy["role_limit"],
        "map_path": str(output.relative_to(root)),
        "map_bytes": output.stat().st_size,
        "map_sha256": sha256(output),
    }
    atomic_json(provenance_path, provenance)
    print(
        f"DENSE_VARIANT_MAP_BUILT rows={mapped_rows} bytes={output.stat().st_size} "
        f"out={output.relative_to(root)}"
    )


def validate_map(root: Path, policy_path: Path, policy: dict[str, object]) -> None:
    source = root / str(policy["source_path"])
    output = root / str(policy["map_path"])
    provenance_path = root / str(policy["map_provenance_path"])
    source_md5, source_sha256 = verify_source(source, policy)
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"dense-map provenance is unreadable: {exc}")
    if (
        provenance.get("schema_version") != policy["map_schema_version"]
        or provenance.get("map_scope") != policy["map_scope"]
        or provenance.get("genome_build") != "GRCh37/hg19"
        or provenance.get("analysis_panel") != policy["analysis_panel"]
        or provenance.get("policy_sha256") != sha256(policy_path)
        or provenance.get("builder_sha256") != policy["builder_sha256"]
        or provenance.get("source", {}).get("bytes") != source.stat().st_size
        or provenance.get("source", {}).get("md5") != source_md5
        or provenance.get("source", {}).get("sha256") != source_sha256
        or provenance.get("mapped_rows", 0) < int(policy["minimum_mapped_rows"])
        or not output.is_file()
        or provenance.get("map_bytes") != output.stat().st_size
        or provenance.get("map_sha256") != sha256(output)
    ):
        fail("dense variant map differs from its immutable policy/provenance seal")
    print(
        f"DENSE_VARIANT_MAP_VALID rows={provenance['mapped_rows']} "
        f"bytes={provenance['map_bytes']}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/dense_variant_map_policy.json")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--acknowledge-large-download", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, policy = load_policy(root, args.policy)
    if args.report_only and (
        args.download or args.acknowledge_large_download or args.build or args.validate_only
    ):
        fail("--report-only cannot be combined with acquisition, build, or validation modes")
    if args.download and not args.acknowledge_large_download:
        fail("--download requires --acknowledge-large-download for the 2.70 GB source")
    if args.acknowledge_large_download and not args.download:
        fail("--acknowledge-large-download is meaningful only with --download")
    if args.validate_only and (args.download or args.build):
        fail("--validate-only cannot be combined with mutation modes")
    if args.download:
        download_source(root, policy)
    if args.build:
        build_map(root, policy_path, policy)
    if args.validate_only:
        validate_map(root, policy_path, policy)
        return 0
    if args.download or args.build:
        return 0
    source = root / str(policy["source_path"])
    output = root / str(policy["map_path"])
    print(
        f"Dense identity source: {policy['source_bytes']:,} bytes; "
        f"source={'READY' if source.is_file() else 'ABSENT'}; "
        f"map={'READY' if output.is_file() else 'ABSENT'}"
    )
    print(
        "No data changed. Use --download --acknowledge-large-download on the "
        "production host, then --build."
    )
    return 0 if args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
