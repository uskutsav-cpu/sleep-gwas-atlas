#!/usr/bin/env python3
"""Stream, body-verify, harmonize, and directly munge frozen replication sources."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import ssl
import tempfile
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from streaming_io import open_verified_gzip_text


ROOT = Path(__file__).resolve().parents[2]
VALID = {"A", "C", "G", "T"}
AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
RSID = re.compile(r"rs[0-9]+")
MIN_FREE_BYTES = 500 * 1024**2


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: str, label: str, counts: Counter[str]) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        counts[f"invalid_{label}"] += 1
        return None
    if not math.isfinite(number):
        counts[f"invalid_{label}"] += 1
        return None
    return number


def load_hm3(reference: Path, allele_list: Path) -> dict[str, tuple[str, str, int, int]]:
    alleles: dict[str, tuple[str, str]] = {}
    with allele_list.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            alleles[row["SNP"]] = (row["A1"].upper(), row["A2"].upper())
    output: dict[str, tuple[str, str, int, int]] = {}
    with gzip.open(reference, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rsid = row["SNP"]
            if rsid not in alleles:
                continue
            a1, a2 = alleles[rsid]
            if a1 not in VALID or a2 not in VALID or frozenset((a1, a2)) in AMBIGUOUS:
                continue
            output[rsid] = (a1, a2, int(row["CHR"]), int(row["BP"]))
    if len(output) < 900_000:
        raise SystemExit(f"ERROR: unexpectedly small nonambiguous HapMap3 reference: {len(output)}")
    return output


def orientation(effect: str, other: str, target_a1: str, target_a2: str) -> int | None:
    effect, other = effect.upper(), other.upper()
    if effect not in VALID or other not in VALID or effect == other:
        return None
    if frozenset((effect, other)) in AMBIGUOUS:
        return None
    if (effect, other) == (target_a1, target_a2):
        return 1
    if (effect, other) == (target_a2, target_a1):
        return -1
    complemented = (effect.translate(COMPLEMENT), other.translate(COMPLEMENT))
    if complemented == (target_a1, target_a2):
        return 1
    if complemented == (target_a2, target_a1):
        return -1
    return None


def source_rsid(value: str, hm3: dict[str, tuple[str, str, int, int]]) -> str | None:
    for rsid in RSID.findall(value):
        if rsid in hm3:
            return rsid
    return None


def head(url: str) -> dict[str, str]:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "sleep-gwas-atlas-extension/1.0"})
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        context = ssl.create_default_context()
    with urllib.request.urlopen(request, context=context, timeout=120) as response:  # nosec B310
        return {key.lower(): value for key, value in response.headers.items()}


def existing_complete(row: dict[str, str], queue_sha: str, reference_sha: str, code_sha: dict[str, str]) -> bool:
    receipt_path = ROOT / row["replication_receipt_path"]
    munged_path = ROOT / row["replication_munged_path"]
    if not receipt_path.is_file() or not munged_path.is_file():
        return False
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return False
    return (
        receipt.get("pipeline_status") == "STREAM_HARMONIZE_DIRECT_MUNGE_PASS"
        and receipt.get("replication_source_id") == row["replication_source_id"]
        and receipt.get("source_url") == row["replication_source_url"]
        and receipt.get("source_verification", {}).get("verification_status") == "PASS"
        and receipt.get("munged_output_sha256") == sha256(munged_path)
        and receipt.get("queue_sha256") == queue_sha
        and receipt.get("reference_sha256") == reference_sha
        and receipt.get("pipeline_code_sha256") == code_sha
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--source-id", action="append")
    selection.add_argument("--all", action="store_true")
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/replication/replication_source_queue.tsv"))
    parser.add_argument("--reference", type=Path, default=Path("discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz"))
    parser.add_argument("--hm3-alleles", type=Path, default=Path("ref/w_hm3.snplist"))
    parser.add_argument("--acknowledge-network-gib", type=float, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    queue = read_tsv(args.queue)
    unique: dict[str, dict[str, str]] = {}
    for row in queue:
        if row["source_curation_status"] == "COMPLETE_BEFORE_RESULTS":
            unique.setdefault(row["replication_source_id"], row)
    if len(unique) != 13:
        raise SystemExit(f"ERROR: expected 13 frozen replication sources, found {len(unique)}")
    selected_ids = list(unique) if args.all else list(args.source_id or [])
    if len(selected_ids) != len(set(selected_ids)) or any(source_id not in unique for source_id in selected_ids):
        raise SystemExit("ERROR: requested source IDs are duplicated or absent from the frozen queue")
    selected = [unique[source_id] for source_id in selected_ids]
    required_bytes = sum(int(row["replication_content_length_bytes"]) for row in selected)
    required_gib = required_bytes / 1024**3
    if args.acknowledge_network_gib + 1e-9 < required_gib:
        raise SystemExit(f"ERROR: selected sources total {required_gib:.6f} GiB; acknowledgement was {args.acknowledge_network_gib:.6f}")
    if not args.reference.is_file() or not args.hm3_alleles.is_file():
        raise SystemExit("ERROR: pinned HapMap3 reference inputs are absent")
    queue_sha, reference_sha = sha256(args.queue), sha256(args.reference)
    code_paths = [
        Path("discovery_extension/scripts/streaming_io.py"),
        Path("discovery_extension/scripts/47_stream_replication_sources.py"),
    ]
    code_sha = {str(path): sha256(path) for path in code_paths}
    print(f"REPLICATION_STREAMING_ACKNOWLEDGED sources={len(selected)} compressed_gib={required_gib:.6f}", flush=True)
    for row in selected:
        headers = head(row["replication_source_url"])
        observed_etag = headers.get("etag", "").strip('"')
        if headers.get("content-length") != row["replication_content_length_bytes"] or observed_etag != row["replication_etag"]:
            raise SystemExit(f"ERROR: remote HEAD identity changed: {row['replication_source_id']}")
        if row["replication_storage_mode"] == "VERSIONED_REMOTE_STREAMING" and headers.get("x-goog-generation") != row["replication_source_generation"]:
            raise SystemExit(f"ERROR: GCS generation changed: {row['replication_source_id']}")
    if args.preflight_only:
        print(f"REPLICATION_STREAMING_PREFLIGHT_PASS sources={len(selected)} free_gib={shutil.disk_usage(ROOT).free / 1024**3:.3f}")
        return

    hm3 = load_hm3(args.reference, args.hm3_alleles)
    print(f"REPLICATION_HM3_REFERENCE_LOADED nonambiguous_snps={len(hm3)}", flush=True)
    for index, row in enumerate(selected, start=1):
        source_id = row["replication_source_id"]
        if existing_complete(row, queue_sha, reference_sha, code_sha):
            print(f"REPLICATION_STREAM_SKIP_VERIFIED index={index}/{len(selected)} source={source_id}", flush=True)
            continue
        if shutil.disk_usage(ROOT).free < MIN_FREE_BYTES:
            raise SystemExit(f"ERROR: less than {MIN_FREE_BYTES / 1024**2:.0f} MiB free before {source_id}")
        output = ROOT / row["replication_munged_path"]
        receipt_path = ROOT / row["replication_receipt_path"]
        qc_path = ROOT / f"discovery_extension/results/replication/qc/{source_id}.tsv"
        output.parent.mkdir(parents=True, exist_ok=True)
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        qc_path.parent.mkdir(parents=True, exist_ok=True)
        handle_fd, temporary_name = tempfile.mkstemp(prefix=output.name + ".", suffix=".tmp", dir=output.parent)
        os.close(handle_fd)
        temporary = Path(temporary_name)
        counts: Counter[str] = Counter()
        seen: set[str] = set()
        cases, controls = int(row["cases"]), int(row["controls"])
        effective_n = 4.0 / (1.0 / cases + 1.0 / controls)
        expected_md5 = row["replication_checksum"]
        print(f"REPLICATION_STREAM_START index={index}/{len(selected)} source={source_id} compressed_gib={int(row['replication_content_length_bytes']) / 1024**3:.6f}", flush=True)
        try:
            with open_verified_gzip_text(
                url=row["replication_source_url"], expected_md5=expected_md5,
                expected_size_bytes=int(row["replication_content_length_bytes"]), timeout_seconds=240,
            ) as (source_handle, source_receipt):
                reader = csv.DictReader(source_handle, delimiter="\t")
                fields = set(reader.fieldnames or [])
                is_finngen = source_id.startswith("finngen_r13_")
                required = (
                    {"ref", "alt", "rsids", "beta", "sebeta", "af_alt"}
                    if is_finngen else
                    {"effect_allele", "other_allele", "odds_ratio", "effect_allele_frequency", "p_value", "rsid", "r2", "ci_upper", "ci_lower"}
                )
                if not required.issubset(fields):
                    raise SystemExit(f"ERROR: {source_id} lacks required fields: {sorted(required-fields)}")
                with gzip.open(temporary, "wt", newline="", encoding="utf-8", compresslevel=1) as destination:
                    writer = csv.DictWriter(destination, delimiter="\t", fieldnames=["SNP", "A1", "A2", "Z", "N"], lineterminator="\n")
                    writer.writeheader()
                    for source_row in reader:
                        counts["input_rows"] += 1
                        rsid = source_rsid(source_row["rsids"] if is_finngen else source_row["rsid"], hm3)
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
                        if is_finngen:
                            effect, other = source_row["alt"], source_row["ref"]
                            frequency = finite(source_row["af_alt"], "frequency", counts)
                            beta = finite(source_row["beta"], "beta", counts)
                            se = finite(source_row["sebeta"], "se", counts)
                            if beta is None or se is None or se <= 0:
                                counts["invalid_beta_or_se"] += 1
                                continue
                            source_z = beta / se
                        else:
                            effect, other = source_row["effect_allele"], source_row["other_allele"]
                            frequency = finite(source_row["effect_allele_frequency"], "frequency", counts)
                            odds_ratio = finite(source_row["odds_ratio"], "odds_ratio", counts)
                            p_value = finite(source_row["p_value"], "p", counts)
                            r2 = finite(source_row["r2"], "r2", counts)
                            ci_upper = finite(source_row["ci_upper"], "ci_upper", counts)
                            ci_lower = finite(source_row["ci_lower"], "ci_lower", counts)
                            if odds_ratio is None or odds_ratio <= 0 or p_value is None or not 0 <= p_value <= 1 or r2 is None or r2 <= 0.9 or ci_upper is None or ci_lower is None or ci_lower <= 0 or ci_upper <= ci_lower:
                                counts["invalid_or_p_ci_or_r2_le_0_9"] += 1
                                continue
                            standard_error = (math.log(ci_upper) - math.log(ci_lower)) / (2 * 1.959963984540054)
                            source_z = math.log(odds_ratio) / standard_error if standard_error > 0 else math.nan
                        if frequency is None or not 0 <= frequency <= 1 or min(frequency, 1-frequency) <= 0.01:
                            counts["maf_at_or_below_0_01_or_invalid"] += 1
                            continue
                        sign = orientation(effect, other, target_a1, target_a2)
                        if sign is None:
                            counts["allele_mismatch_or_ambiguous"] += 1
                            continue
                        if not math.isfinite(source_z):
                            counts["invalid_z"] += 1
                            continue
                        writer.writerow({"SNP": rsid, "A1": target_a1, "A2": target_a2, "Z": f"{sign * source_z:.12g}", "N": f"{effective_n:.12g}"})
                        seen.add(rsid); counts["output_rows"] += 1
            temporary.replace(output)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
        if counts["output_rows"] < 700_000:
            output.unlink(missing_ok=True)
            raise SystemExit(f"ERROR: too few HM3 variants survived for {source_id}: {counts['output_rows']}")
        with qc_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["metric", "value", "notes"], lineterminator="\n")
            writer.writeheader()
            for metric, value in sorted(counts.items()):
                writer.writerow({"metric": metric, "value": value, "notes": "frozen replication harmonization filter"})
            writer.writerow({"metric": "effective_sample_size", "value": f"{effective_n:.12g}", "notes": "4/(1/cases+1/controls)"})
            writer.writerow({"metric": "munged_output_sha256", "value": sha256(output), "notes": str(output.relative_to(ROOT))})
        receipt = {
            "schema_version": "1.0.0", "sealed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "pipeline_status": "STREAM_HARMONIZE_DIRECT_MUNGE_PASS", "replication_results_accessed": False,
            "replication_source_id": source_id, "source_url": row["replication_source_url"],
            "source_verification": source_receipt, "full_resolution_local_retention": False,
            "full_resolution_remote_identity": "generation or study-accession path plus locked MD5, ETag, and byte length",
            "munged_output": row["replication_munged_path"], "munged_output_sha256": sha256(output),
            "harmonization_qc": str(qc_path.relative_to(ROOT)), "harmonization_qc_sha256": sha256(qc_path),
            "queue": str(args.queue), "queue_sha256": queue_sha,
            "reference": str(args.reference), "reference_sha256": reference_sha,
            "hm3_alleles": str(args.hm3_alleles), "hm3_alleles_sha256": sha256(args.hm3_alleles),
            "pipeline_code_sha256": code_sha,
        }
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"REPLICATION_STREAM_PASS index={index}/{len(selected)} source={source_id} rows={counts['output_rows']} sha256={receipt['munged_output_sha256']}", flush=True)
    print(f"REPLICATION_STREAMING_FAMILY_PASS sources={len(selected)} compressed_gib={required_gib:.6f}", flush=True)


if __name__ == "__main__":
    main()
