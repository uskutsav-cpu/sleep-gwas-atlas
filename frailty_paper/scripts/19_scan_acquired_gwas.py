#!/usr/bin/env python3
"""Stream acquired GWAS files for readability, headers, and basic numeric anomalies."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


ALIASES = {
    "variant": ("variant_id", "rsid", "rs_id", "rsids", "snpid", "snp", "id"),
    "chromosome": ("chromosome", "chr", "chrom", "#chrom"),
    "position": ("base_pair_location", "position", "pos", "bp", "genpos"),
    "effect_allele": ("effect_allele", "ea", "a1", "allele1", "alt"),
    "other_allele": ("other_allele", "nea", "oa", "ra", "a0", "a2", "allele0", "ref"),
    "beta": ("beta", "beta1", "effect", "log_odds", "logor", "beta_sleepduration", "beta_shortsleep", "beta_longsleep"),
    "odds_ratio": ("or", "odds_ratio", "oddsratio"),
    "se": ("standard_error", "se", "stderr", "standard.error", "sebeta", "logor_se", "se_sleepduration", "se_shortsleep", "se_longsleep"),
    "p": ("p_value", "p-value", "pvalue", "p", "pval", "p_bolt_lmm", "p_sleepduration", "p_shortsleep", "p_longsleep"),
    "log10p": ("log10p", "neg_log10_p", "minus_log10_p", "-log10(p-value)"),
    "eaf": ("effect_allele_frequency", "eaf", "freq", "freq1", "a1freq"),
}
OUTPUT_FIELDS = [
    "resource_id", "file", "manifest_bytes", "manifest_sha256", "genome_build",
    "ancestry", "sample_size", "stream_status", "header_status", "header",
    "data_rows", "rows_missing_variant", "rows_missing_chromosome",
    "rows_invalid_position", "rows_missing_alleles", "rows_missing_beta_or_se",
    "rows_invalid_beta_or_se", "rows_missing_p", "rows_invalid_p_or_log10p",
    "rows_invalid_eaf", "duplicates_assessed", "anomaly_examples", "qc_interpretation", "notes",
]


def normalize_header(header: list[str]) -> dict[str, str]:
    lookup = {name.strip().lower(): name for name in header}
    return {
        key: next((lookup[name] for name in names if name in lookup), "")
        for key, names in ALIASES.items()
    }


def number(value: str) -> float | None:
    value = value.strip()
    if not value or value.lower() in {"na", "nan", "null", "."}:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def scan_file(path: Path) -> dict[str, object]:
    counts = {
        "rows_missing_variant": 0,
        "rows_missing_chromosome": 0,
        "rows_invalid_position": 0,
        "rows_missing_alleles": 0,
        "rows_missing_beta_or_se": 0,
        "rows_invalid_beta_or_se": 0,
        "rows_missing_p": 0,
        "rows_invalid_p_or_log10p": 0,
        "rows_invalid_eaf": 0,
    }
    nrows = 0
    examples = []
    try:
        delimiter = "," if path.name.lower().endswith((".csv", ".csv.gz")) else "\t"
        opener = gzip.open if path.suffix.lower() == ".gz" else Path.open
        open_args = (path, "rt") if opener is gzip.open else (path, "r")
        with opener(*open_args, encoding="utf-8-sig", errors="strict", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            header = reader.fieldnames or []
            columns = normalize_header(header)
            required = ("chromosome", "position", "effect_allele", "other_allele", "se")
            missing_columns = [key for key in required if not columns[key]]
            if not columns["beta"] and not columns["odds_ratio"]:
                missing_columns.append("beta_or_odds_ratio")
            if not columns["variant"] and not all(columns[key] for key in ("chromosome", "position", "effect_allele", "other_allele")):
                missing_columns.append("variant_or_coordinate_allele_key")
            if not columns["p"] and not columns["log10p"]:
                missing_columns.append("p_or_log10p")
            for row in reader:
                nrows += 1
                get = lambda key: (row.get(columns[key], "") or "").strip() if columns[key] else ""
                issues = []
                variant = get("variant")
                if not variant and all(columns[key] for key in ("chromosome", "position", "effect_allele", "other_allele")):
                    variant = "|".join(get(key) for key in ("chromosome", "position", "effect_allele", "other_allele"))
                if not variant:
                    counts["rows_missing_variant"] += 1
                if not get("chromosome"):
                    counts["rows_missing_chromosome"] += 1
                pos = number(get("position"))
                if pos is None or pos <= 0 or not pos.is_integer():
                    counts["rows_invalid_position"] += 1
                    issues.append("position")
                if not get("effect_allele") or not get("other_allele"):
                    counts["rows_missing_alleles"] += 1
                beta, se = number(get("beta")), number(get("se"))
                if beta is None and columns["odds_ratio"]:
                    odds_ratio = number(get("odds_ratio"))
                    # Validate the source OR on the log-effect scale used by
                    # harmonization; this is a scan only and does not rewrite
                    # or filter the raw source.
                    beta = math.log(odds_ratio) if odds_ratio is not None and odds_ratio > 0 else None
                if beta is None or se is None:
                    counts["rows_missing_beta_or_se"] += 1
                    issues.append("missing_beta_or_se")
                elif se <= 0:
                    counts["rows_invalid_beta_or_se"] += 1
                    issues.append("nonpositive_se")
                pvalue = number(get("p")) if columns["p"] else None
                log10p = number(get("log10p")) if columns["log10p"] else None
                if pvalue is None and log10p is None:
                    counts["rows_missing_p"] += 1
                    issues.append("missing_p")
                elif (pvalue is not None and not 0 <= pvalue <= 1) or (log10p is not None and log10p < 0):
                    counts["rows_invalid_p_or_log10p"] += 1
                    issues.append("p_or_log10p_domain")
                if columns["eaf"]:
                    eaf = number(get("eaf"))
                    if eaf is None or not 0 <= eaf <= 1:
                        counts["rows_invalid_eaf"] += 1
                        issues.append("eaf_domain")
                if issues and len(examples) < 10:
                    examples.append({
                        "data_row": nrows, "variant": variant,
                        "chromosome": get("chromosome"), "position": get("position"),
                        "beta": get("beta"), "se": get("se"),
                        "p": get("p") or get("log10p"), "issues": issues,
                    })
            status = "GZIP_EOF_OK" if path.suffix.lower() == ".gz" else "TEXT_EOF_OK"
    except (OSError, EOFError, UnicodeError, csv.Error) as exc:
        return {"stream_status": f"ERROR: {type(exc).__name__}: {exc}", "header_status": "NOT_VERIFIED", "header": "", "data_rows": nrows, "anomaly_examples": json.dumps(examples, separators=(",", ":")), **counts}

    header_status = "MISSING_COLUMNS:" + ",".join(missing_columns) if missing_columns else (
        "REQUIRED_COLUMNS_PRESENT;OR_TO_LOG_EFFECT_FOR_SCAN" if not columns["beta"] and columns["odds_ratio"] else
        "REQUIRED_COLUMNS_PRESENT" if columns["variant"] else "REQUIRED_COLUMNS_PRESENT;COMPOSITE_COORDINATE_ALLELE_KEY"
    )
    return {
        "stream_status": status,
        "header_status": header_status,
        "header": "|".join(header),
        "data_rows": nrows,
        "anomaly_examples": json.dumps(examples, separators=(",", ":")),
        **counts,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=Path("frailty_paper/manifests/gwas_source_scan.tsv"))
    parser.add_argument("--resource-id", action="append", help="Scan only this exact resource_id (repeatable); useful for targeted rechecks")
    args = parser.parse_args()
    repo = args.repo.resolve()
    manifest = repo / "frailty_paper/manifests/all_acquired_resources.tsv"
    rows = []
    with manifest.open(newline="", encoding="utf-8") as handle:
        for item in csv.DictReader(handle, delimiter="\t"):
            if item["resource_type"] != "GWAS summary statistics":
                continue
            if args.resource_id and item["resource_id"] not in args.resource_id:
                continue
            path = repo / item["file"]
            result = scan_file(path) if path.is_file() else {
                "stream_status": "MISSING_FILE", "header_status": "NOT_VERIFIED", "header": "", "data_rows": 0,
                **{key: 0 for key in OUTPUT_FIELDS if key.startswith("rows_")},
            }
            numeric_flags = sum(int(result[key]) for key in OUTPUT_FIELDS if key.startswith("rows_") and isinstance(result[key], int))
            if not path.is_file():
                interpretation = "MISSING_FILE"
            elif result["stream_status"] not in {"GZIP_EOF_OK", "TEXT_EOF_OK"} or not str(result["header_status"]).startswith("REQUIRED_COLUMNS_PRESENT"):
                interpretation = "STRUCTURE_OR_STREAM_ERROR"
            elif numeric_flags:
                interpretation = "ANOMALIES_REPORTED_NO_FILTERING_OR_THRESHOLD_APPLIED"
            else:
                interpretation = "STREAM_AND_BASIC_DOMAINS_OK_NOT_FULL_SCIENTIFIC_QC"
            notes = "Duplicate variants, allele compatibility, build, ancestry, and study-specific QC gates are not assessed by this scanner. No records are changed or filtered."
            if "COMPOSITE_COORDINATE_ALLELE_KEY" in str(result["header_status"]):
                notes += " Source has no variant-ID column; basic field checks use a composite chromosome:position:allele key and infer no rsID."
            rows.append({
                **{key: item.get(source, "") for key, source in {
                    "resource_id": "resource_id", "file": "file", "manifest_bytes": "bytes",
                    "manifest_sha256": "sha256", "genome_build": "genome_build", "ancestry": "ancestry", "sample_size": "sample_size",
                }.items()},
                **result,
                "duplicates_assessed": "NO",
                "qc_interpretation": interpretation,
                "notes": notes,
            })
    out = args.output if args.output.is_absolute() else repo / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"GWAS_SOURCE_SCAN_WRITTEN rows={len(rows)} output={out}")
    for row in rows:
        print(f"{row['resource_id']} rows={row['data_rows']} status={row['qc_interpretation']}")


if __name__ == "__main__":
    main()
