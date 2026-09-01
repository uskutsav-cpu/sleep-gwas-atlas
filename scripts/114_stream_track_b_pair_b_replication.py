#!/usr/bin/env python3
"""Stream, verify, harmonize, and munge the frozen Pair B FinnGen source."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import shutil
import ssl
import sys
import tempfile
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_AUDIT = ROOT / "results/track_b/replication_source_audit.tsv"
SOURCE_LOCK = ROOT / "results/track_b/replication_source.lock.json"
REFERENCE = ROOT / "discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz"
HM3_ALLELES = ROOT / "ref/w_hm3.snplist"
OUTPUT = ROOT / "data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz"
QC = ROOT / "results/track_b/replication/finngen_r13_F5_ADHD_ingest_qc.tsv"
RECEIPT = ROOT / "results/track_b/replication/finngen_r13_F5_ADHD_ingest_receipt.json"
MIN_FREE_BYTES = 500 * 1024**2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_existing_helpers():
    scripts = ROOT / "discovery_extension/scripts"
    sys.path.insert(0, str(scripts))
    path = scripts / "47_stream_replication_sources.py"
    spec = importlib.util.spec_from_file_location("atlas_replication_stream", path)
    if spec is None or spec.loader is None:
        raise SystemExit("ERROR: cannot load canonical replication streaming helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def selected_source() -> dict[str, str]:
    rows = read_tsv(SOURCE_AUDIT)
    selected = [row for row in rows if row["selection_decision"] == "SELECT_PRIMARY_REPLICATION"]
    if len(selected) != 1 or selected[0]["source_id"] != "finngen_r13_F5_ADHD":
        raise SystemExit("ERROR: frozen Pair B replication source drifted")
    row = selected[0]
    lock = json.loads(SOURCE_LOCK.read_text(encoding="utf-8"))
    if lock["replication_results_accessed_before_selection"] is not False:
        raise SystemExit("ERROR: source selection did not precede results")
    if (
        row["source_generation"] != str(lock["pair_b_source_generation"])
        or row["md5_hex"] != lock["pair_b_source_md5"]
        or int(row["content_length_bytes"]) != lock["pair_b_source_content_length_bytes"]
    ):
        raise SystemExit("ERROR: selected source differs from immutable source lock")
    return row


def head(url: str) -> dict[str, str]:
    request = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": "sleep-gwas-atlas-track-b/1.0"}
    )
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        context = ssl.create_default_context()
    with urllib.request.urlopen(request, context=context, timeout=120) as response:  # nosec B310
        return {key.lower(): value for key, value in response.headers.items()}


def remote_identity(row: dict[str, str]) -> dict[str, str]:
    headers = head(row["source_url"])
    observed_etag = headers.get("etag", "").strip('"')
    if (
        headers.get("content-length") != row["content_length_bytes"]
        or headers.get("x-goog-generation") != row["source_generation"]
        or observed_etag != row["etag"]
    ):
        raise SystemExit("ERROR: frozen FinnGen remote identity changed")
    return {
        "content_length_bytes": headers["content-length"],
        "source_generation": headers["x-goog-generation"],
        "etag": observed_etag,
        "last_modified": headers.get("last-modified", "NA"),
    }


def verified_existing(row: dict[str, str]) -> bool:
    if not RECEIPT.is_file() or not OUTPUT.is_file() or not QC.is_file():
        return False
    try:
        receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        receipt.get("pipeline_status") == "STREAM_HARMONIZE_DIRECT_MUNGE_PASS"
        and receipt.get("source_id") == row["source_id"]
        and receipt.get("source_url") == row["source_url"]
        and receipt.get("source_verification", {}).get("verification_status") == "PASS"
        and receipt.get("munged_output_sha256") == sha256(OUTPUT)
        and receipt.get("qc_sha256") == sha256(QC)
        and receipt.get("source_lock_sha256") == sha256(SOURCE_LOCK)
    )


def write_qc(counts: Counter[str], effective_n: float) -> None:
    QC.parent.mkdir(parents=True, exist_ok=True)
    with QC.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["metric", "value", "notes"], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for metric, value in sorted(counts.items()):
            writer.writerow({"metric": metric, "value": value, "notes": "frozen Track B replication harmonization filter"})
        writer.writerow({"metric": "effective_sample_size", "value": f"{effective_n:.12g}", "notes": "4/(1/cases+1/controls)"})
        writer.writerow({"metric": "source_info_filter", "value": ">0.6", "notes": "FinnGen R13 public-release source-level filter"})
        writer.writerow({"metric": "munged_output_sha256", "value": sha256(OUTPUT), "notes": str(OUTPUT.relative_to(ROOT))})


def stream(row: dict[str, str], helpers) -> None:
    if shutil.disk_usage(ROOT).free < MIN_FREE_BYTES:
        raise SystemExit("ERROR: BLOCKED_BY_COMPUTE less than 500 MiB free for atomic munged output")
    hm3 = helpers.load_hm3(REFERENCE, HM3_ALLELES)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    handle_fd, temporary_name = tempfile.mkstemp(prefix=OUTPUT.name + ".", suffix=".tmp", dir=OUTPUT.parent)
    os.close(handle_fd)
    temporary = Path(temporary_name)
    counts: Counter[str] = Counter()
    seen: set[str] = set()
    cases, controls = int(row["cases"]), int(row["controls"])
    effective_n = 4.0 / (1.0 / cases + 1.0 / controls)
    try:
        with helpers.open_verified_gzip_text(
            url=row["source_url"], expected_md5=f"md5:{row['md5_hex']}",
            expected_size_bytes=int(row["content_length_bytes"]), timeout_seconds=240,
        ) as (source_handle, source_receipt):
            reader = csv.DictReader(source_handle, delimiter="\t")
            required = {"ref", "alt", "rsids", "beta", "sebeta", "af_alt"}
            if not required.issubset(set(reader.fieldnames or [])):
                raise SystemExit(f"ERROR: FinnGen source lacks fields: {sorted(required-set(reader.fieldnames or []))}")
            with gzip.open(temporary, "wt", newline="", encoding="utf-8", compresslevel=1) as destination:
                writer = csv.DictWriter(
                    destination, fieldnames=["SNP", "A1", "A2", "Z", "N"],
                    delimiter="\t", lineterminator="\n",
                )
                writer.writeheader()
                for source_row in reader:
                    counts["input_rows"] += 1
                    rsid = helpers.source_rsid(source_row["rsids"], hm3)
                    if rsid is None:
                        counts["not_nonambiguous_hm3"] += 1
                        continue
                    if rsid in seen:
                        counts["duplicate_hm3_rsid"] += 1
                        continue
                    target_a1, target_a2, chromosome, position = hm3[rsid]
                    if chromosome == 6 and 25_000_000 <= position <= 34_000_000:
                        counts["extended_MHC"] += 1
                        continue
                    frequency = helpers.finite(source_row["af_alt"], "frequency", counts)
                    beta = helpers.finite(source_row["beta"], "beta", counts)
                    se = helpers.finite(source_row["sebeta"], "se", counts)
                    if beta is None or se is None or se <= 0:
                        counts["invalid_beta_or_se"] += 1
                        continue
                    if frequency is None or not 0 <= frequency <= 1 or min(frequency, 1 - frequency) <= 0.01:
                        counts["maf_at_or_below_0_01_or_invalid"] += 1
                        continue
                    sign = helpers.orientation(source_row["alt"], source_row["ref"], target_a1, target_a2)
                    if sign is None:
                        counts["allele_mismatch_or_ambiguous"] += 1
                        continue
                    z = sign * beta / se
                    if not math.isfinite(z):
                        counts["invalid_z"] += 1
                        continue
                    writer.writerow({"SNP": rsid, "A1": target_a1, "A2": target_a2, "Z": f"{z:.12g}", "N": f"{effective_n:.12g}"})
                    seen.add(rsid)
                    counts["output_rows"] += 1
        if counts["output_rows"] < 700_000:
            raise SystemExit(f"ERROR: too few HapMap3 variants survived: {counts['output_rows']}")
        temporary.replace(OUTPUT)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    write_qc(counts, effective_n)
    receipt = {
        "schema_version": 1,
        "sealed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pipeline_status": "STREAM_HARMONIZE_DIRECT_MUNGE_PASS",
        "replication_results_accessed": False,
        "source_id": row["source_id"],
        "source_url": row["source_url"],
        "source_verification": source_receipt,
        "remote_full_resolution_retained_locally": False,
        "munged_output": str(OUTPUT.relative_to(ROOT)),
        "munged_output_sha256": sha256(OUTPUT),
        "qc": str(QC.relative_to(ROOT)),
        "qc_sha256": sha256(QC),
        "source_lock": str(SOURCE_LOCK.relative_to(ROOT)),
        "source_lock_sha256": sha256(SOURCE_LOCK),
        "reference": str(REFERENCE.relative_to(ROOT)),
        "reference_sha256": sha256(REFERENCE),
        "hm3_alleles": str(HM3_ALLELES.relative_to(ROOT)),
        "hm3_alleles_sha256": sha256(HM3_ALLELES),
        "pipeline_code_sha256": {
            "scripts/114_stream_track_b_pair_b_replication.py": sha256(Path(__file__)),
            "discovery_extension/scripts/47_stream_replication_sources.py": sha256(ROOT / "discovery_extension/scripts/47_stream_replication_sources.py"),
            "discovery_extension/scripts/streaming_io.py": sha256(ROOT / "discovery_extension/scripts/streaming_io.py"),
        },
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"TRACK_B_REPLICATION_STREAM_PASS rows={counts['output_rows']} sha256={receipt['munged_output_sha256']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--acknowledge-network-gib", type=float, required=True)
    args = parser.parse_args()
    row = selected_source()
    required_gib = int(row["content_length_bytes"]) / 1024**3
    if args.acknowledge_network_gib + 1e-9 < required_gib:
        raise SystemExit(f"ERROR: source is {required_gib:.6f} GiB compressed")
    if not REFERENCE.is_file() or not HM3_ALLELES.is_file():
        raise SystemExit("ERROR: pinned HapMap3 reference inputs are absent")
    if args.verify:
        if not verified_existing(row):
            raise SystemExit("ERROR: Track B Pair B replication ingest is absent or drifted")
        print(f"verified Track B Pair B replication ingest: {sha256(OUTPUT)}")
        return
    identity = remote_identity(row)
    if args.preflight_only:
        free_bytes = shutil.disk_usage(ROOT).free
        status = "PASS" if free_bytes >= MIN_FREE_BYTES else "BLOCKED_BY_COMPUTE"
        print(
            f"TRACK_B_REPLICATION_PREFLIGHT_{status} "
            f"compressed_gib={required_gib:.6f} free_gib={free_bytes/1024**3:.3f} "
            f"minimum_free_gib={MIN_FREE_BYTES/1024**3:.3f} "
            f"generation={identity['source_generation']} etag={identity['etag']}"
        )
        if status != "PASS":
            raise SystemExit(2)
        return
    if verified_existing(row):
        print(f"TRACK_B_REPLICATION_STREAM_SKIP_VERIFIED sha256={sha256(OUTPUT)}")
        return
    helpers = load_existing_helpers()
    stream(row, helpers)


if __name__ == "__main__":
    main()
