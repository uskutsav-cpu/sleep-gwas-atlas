#!/usr/bin/env python3
"""Integrity- and schema-audit one reviewed public substitution candidate.

The scan traverses every decompressed row, which also validates the gzip CRC.
It records exclusions rather than silently treating missing rsIDs, non-SNP
alleles, or invalid statistics as analysis-ready. Malformed row widths fail.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re


SCHEMA_VERSION = "atlas.substitution-source-audit.v1"
FINNGEN_HEADER = [
    "#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval",
    "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls",
]
BURREN_HEADER = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value", "rs_id",
    "info", "ci_lower", "ci_upper", "n",
]
AUTOSOMES = {str(chromosome) for chromosome in range(1, 23)}
DNA = {"A", "C", "G", "T"}


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def digest(path, algorithm):
    value = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def load_candidate(path, candidate_source_id):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["candidate_source_id"] == candidate_source_id
        ]
    if len(rows) != 1:
        fail(
            f"expected exactly one candidate row for {candidate_source_id}; "
            f"found {len(rows)}"
        )
    return rows[0]


def schema_for(candidate):
    if candidate["candidate_source_id"].startswith("finngen_r9_"):
        return "FINNGEN_R9", FINNGEN_HEADER
    if candidate["candidate_source_id"] == "burren_2024_telomere_nfe":
        return "GWAS_SSF_BURREN_NFE", BURREN_HEADER
    fail(f"no reviewed audit schema for {candidate['candidate_source_id']}")


def audit_source(path, candidate):
    schema_name, expected_header = schema_for(candidate)
    archive_bytes = os.path.getsize(path)
    expected_bytes = int(candidate["archive_bytes"])
    if archive_bytes != expected_bytes:
        fail(
            f"archive byte count mismatch: expected {expected_bytes}, "
            f"found {archive_bytes}"
        )
    archive_sha256 = digest(path, "sha256")
    if archive_sha256 != candidate["archive_sha256"]:
        fail("archive SHA-256 does not match the candidate registry")
    archive_md5 = digest(path, "md5")
    expected_md5 = candidate["archive_md5"]
    if re.fullmatch(r"[0-9a-f]{32}", expected_md5) and archive_md5 != expected_md5:
        fail("archive MD5 does not match the upstream candidate registry value")

    counts = {
        "source_rows": 0,
        "malformed_rows": 0,
        "autosomal_rows": 0,
        "explicit_rsid_rows": 0,
        "snp_allele_rows": 0,
        "valid_effect_se_p_rows": 0,
        "valid_frequency_rows": 0,
    }
    first_malformed_lines = []
    n_min = None
    n_max = None
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header != expected_header:
            fail(
                f"literal header mismatch for {schema_name}: expected "
                f"{expected_header!r}, found {header!r}"
            )
        indexes = {name: index for index, name in enumerate(header)}
        if schema_name == "FINNGEN_R9":
            names = {
                "chrom": "#chrom", "a1": "alt", "a2": "ref", "rsid": "rsids",
                "beta": "beta", "se": "sebeta", "p": "pval", "frq": "af_alt",
            }
        else:
            names = {
                "chrom": "chromosome", "a1": "effect_allele",
                "a2": "other_allele", "rsid": "rs_id", "beta": "beta",
                "se": "standard_error", "p": "p_value",
                "frq": "effect_allele_frequency", "n": "n",
            }
        for line_number, line in enumerate(handle, start=2):
            counts["source_rows"] += 1
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                counts["malformed_rows"] += 1
                if len(first_malformed_lines) < 10:
                    first_malformed_lines.append(line_number)
                continue
            chrom = fields[indexes[names["chrom"]]].removeprefix("chr")
            if chrom in AUTOSOMES:
                counts["autosomal_rows"] += 1
            rsid = fields[indexes[names["rsid"]]].strip()
            if re.fullmatch(r"rs[0-9]+", rsid, flags=re.IGNORECASE):
                counts["explicit_rsid_rows"] += 1
            a1 = fields[indexes[names["a1"]]].upper()
            a2 = fields[indexes[names["a2"]]].upper()
            if a1 in DNA and a2 in DNA:
                counts["snp_allele_rows"] += 1
            beta = finite_number(fields[indexes[names["beta"]]])
            se = finite_number(fields[indexes[names["se"]]])
            p_value = finite_number(fields[indexes[names["p"]]])
            if beta is not None and se is not None and se > 0 and p_value is not None and 0 <= p_value <= 1:
                counts["valid_effect_se_p_rows"] += 1
            frequency = finite_number(fields[indexes[names["frq"]]])
            if frequency is not None and 0 <= frequency <= 1:
                counts["valid_frequency_rows"] += 1
            if "n" in names:
                n_value = finite_number(fields[indexes[names["n"]]])
                if n_value is not None and n_value > 0:
                    n_min = n_value if n_min is None else min(n_min, n_value)
                    n_max = n_value if n_max is None else max(n_max, n_value)
            if counts["source_rows"] % 5_000_000 == 0:
                print(
                    f"Scanned {counts['source_rows']:,} rows from "
                    f"{candidate['candidate_source_id']}",
                    flush=True,
                )

    if counts["source_rows"] == 0:
        fail("source contains a header but no data rows")
    if counts["malformed_rows"]:
        fail(
            f"source has {counts['malformed_rows']} malformed rows; first line(s) "
            f"{first_malformed_lines}"
        )
    if counts["valid_effect_se_p_rows"] == 0:
        fail("source has no rows with valid effect, SE, and P")
    return {
        "schema_version": SCHEMA_VERSION,
        "candidate_source_id": candidate["candidate_source_id"],
        "trait_id": candidate["trait_id"],
        "archive_name": candidate["archive_name"],
        "archive_bytes": archive_bytes,
        "archive_md5": archive_md5,
        "archive_sha256": archive_sha256,
        "literal_header": expected_header,
        "schema": schema_name,
        **counts,
        "per_variant_n_min": n_min,
        "per_variant_n_max": n_max,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_source_id")
    parser.add_argument(
        "--candidates", default="config/public_gwas_substitution_candidates.tsv"
    )
    parser.add_argument("--archive-dir", default="data/raw/.archives")
    parser.add_argument("--out-dir", default="data/raw/.audits")
    args = parser.parse_args()
    candidate = load_candidate(args.candidates, args.candidate_source_id)
    path = os.path.join(args.archive_dir, candidate["archive_name"])
    if not os.path.isfile(path):
        fail(f"candidate archive not found: {path}")
    report = audit_source(path, candidate)
    os.makedirs(args.out_dir, exist_ok=True)
    output = os.path.join(args.out_dir, f"{args.candidate_source_id}.json")
    temporary = f"{output}.tmp.{os.getpid()}"
    with open(temporary, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, output)
    print(
        f"Audited {report['source_rows']:,} rows for "
        f"{args.candidate_source_id} -> {output}"
    )


if __name__ == "__main__":
    main()
