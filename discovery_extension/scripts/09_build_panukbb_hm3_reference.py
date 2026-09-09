#!/usr/bin/env python3
"""Build the Pan-UKB-to-pinned-HapMap3 identity/INFO reference once."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from streaming_io import open_verified_gzip_text


VALID = {"A", "C", "G", "T"}
RSID = re.compile(r"rs[0-9]+")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def allele_match(ref: str, alt: str, a1: str, a2: str) -> bool:
    source = {ref, alt}
    target = {a1, a2}
    return source == target or {ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)} == target


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--variant-manifest", type=Path)
    source.add_argument("--variant-manifest-url")
    parser.add_argument("--variant-manifest-checksum")
    parser.add_argument("--variant-manifest-size-bytes", type=int)
    parser.add_argument(
        "--hm3-map", type=Path, default=Path("ref/hm3_grch37_variant_map.tsv.gz")
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz"),
    )
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/panukbb/hm3_variant_reference_build.json"),
    )
    args = parser.parse_args()
    if args.variant_manifest is not None and not args.variant_manifest.is_file():
        raise SystemExit(f"ERROR: required input is missing: {args.variant_manifest}")
    if not args.hm3_map.is_file():
        raise SystemExit(f"ERROR: required input is missing: {args.hm3_map}")

    by_coordinate: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
    with gzip.open(args.hm3_map, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or []) != {"SNP", "CHR", "BP", "A1", "A2"}:
            raise SystemExit("ERROR: unexpected pinned HapMap3 map schema")
        for row in reader:
            by_coordinate[(int(row["CHR"]), int(row["BP"]))].append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    seen_keys: set[tuple[int, int, str, str]] = set()
    seen_rsids: set[str] = set()
    try:
        with open_verified_gzip_text(
            local_path=args.variant_manifest,
            url=args.variant_manifest_url,
            expected_md5=args.variant_manifest_checksum,
            expected_size_bytes=args.variant_manifest_size_bytes,
        ) as (source_handle, source_receipt):
            with gzip.open(args.out, "wt", newline="", encoding="utf-8") as destination:
                reader = csv.DictReader(source_handle, delimiter="\t")
                required = {"chrom", "pos", "ref", "alt", "rsid", "info"}
                missing = required - set(reader.fieldnames or [])
                if missing:
                    raise SystemExit(f"ERROR: Pan-UKB variant manifest lacks fields: {sorted(missing)}")
                writer = csv.DictWriter(
                    destination,
                    delimiter="\t",
                    fieldnames=["SNP", "CHR", "BP", "REF", "ALT", "INFO"],
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in reader:
                    counts["input_rows"] += 1
                    try:
                        chromosome = int(row["chrom"])
                        position = int(row["pos"])
                    except ValueError:
                        counts["invalid_coordinate"] += 1
                        continue
                    candidates = by_coordinate.get((chromosome, position))
                    if not candidates:
                        counts["not_hapmap3_coordinate"] += 1
                        continue
                    ref, alt, rsid = row["ref"].upper(), row["alt"].upper(), row["rsid"]
                    if ref not in VALID or alt not in VALID or ref == alt:
                        counts["non_snp_or_invalid_alleles"] += 1
                        continue
                    if not RSID.fullmatch(rsid):
                        counts["missing_or_invalid_rsid"] += 1
                        continue
                    matches = [
                        candidate for candidate in candidates
                        if candidate["SNP"] == rsid
                        and allele_match(ref, alt, candidate["A1"], candidate["A2"])
                    ]
                    if not matches:
                        counts["hm3_identity_or_allele_conflict"] += 1
                        continue
                    if len(matches) != 1:
                        counts["ambiguous_hm3_match"] += 1
                        continue
                    key = (chromosome, position, ref, alt)
                    if key in seen_keys or rsid in seen_rsids:
                        raise SystemExit(f"ERROR: duplicate Pan-UKB/HapMap3 identity: {key} {rsid}")
                    try:
                        info = float(row["info"])
                    except ValueError:
                        counts["missing_or_invalid_info"] += 1
                        continue
                    if not 0 <= info <= 1:
                        counts["missing_or_invalid_info"] += 1
                        continue
                    writer.writerow(
                        {"SNP": rsid, "CHR": chromosome, "BP": position, "REF": ref, "ALT": alt, "INFO": row["info"]}
                    )
                    seen_keys.add(key)
                    seen_rsids.add(rsid)
                    counts["output_rows"] += 1
    except BaseException:
        args.out.unlink(missing_ok=True)
        raise

    provenance = {
        "schema_version": "1.0.0",
        "variant_manifest": str(args.variant_manifest or args.variant_manifest_url),
        "variant_manifest_checksum_expected": args.variant_manifest_checksum or "UNSPECIFIED_LOCAL_TEST",
        "variant_manifest_size_expected": args.variant_manifest_size_bytes,
        "variant_manifest_stream_receipt": source_receipt,
        "variant_manifest_sha256": source_receipt["observed_sha256"],
        "hm3_map": str(args.hm3_map),
        "hm3_map_sha256": sha256(args.hm3_map),
        "matching": "same rsID, GRCh37 coordinate, and allele pair (direct or strand complement); never guess conflicts",
        "info_filtering": "not_applied_here; numeric INFO is retained for per-trait harmonization at >=0.9",
        "counts": dict(sorted(counts.items())),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(
        f"PANUKBB_HM3_REFERENCE_OK rows={counts['output_rows']} "
        f"sha256={provenance['output_sha256']}"
    )


if __name__ == "__main__":
    main()
