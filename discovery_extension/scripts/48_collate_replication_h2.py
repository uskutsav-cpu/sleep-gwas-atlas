#!/usr/bin/env python3
"""Collate the frozen 13-source replication h2 gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


H2 = re.compile(r"Total (Observed|Liability) scale h2:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)")
INTERCEPT = re.compile(r"^Intercept:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)", re.M)
INPUT = re.compile(r"Read summary statistics for ([0-9]+) SNPs\.")
REGRESSION = re.compile(r"After merging with regression SNP LD, ([0-9]+) SNPs remain\.")
MEAN_CHI = re.compile(r"Mean Chi\^2:\s*(-?[\d.eE+-]+)")
LAMBDA = re.compile(r"Lambda GC:\s*(-?[\d.eE+-]+)")
RATIO = re.compile(r"^Ratio:\s*(-?[\d.eE+-]+)", re.M)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric(match: re.Match[str] | None, group: int = 1) -> float:
    return float(match.group(group)) if match else math.nan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/replication_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/replication_manifest.lock.json"))
    parser.add_argument("--logdir", type=Path, default=Path("discovery_extension/logs/replication/h2"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/replication/replication_source_h2.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/replication_h2.json"))
    args = parser.parse_args()
    manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    if sha256(args.manifest) != lock.get("manifest_sha256"):
        raise SystemExit("ERROR: replication manifest differs from lock")
    sources: dict[str, dict[str, str]] = {}
    testable = set(lock["testable_pair_ids_in_locked_order"])
    for row in manifest:
        if row["pair_id"] in testable:
            sources.setdefault(row["replication_source_id"], row)
    if len(sources) != 13:
        raise SystemExit(f"ERROR: expected 13 locked replication sources, found {len(sources)}")
    output: list[dict[str, object]] = []
    for source_id, source in sources.items():
        log = args.logdir / f"h2_{source_id}.log"
        if not log.is_file():
            raise SystemExit(f"ERROR: replication h2 log missing: {log}")
        text = log.read_text(encoding="utf-8")
        h2_match, intercept_match = H2.search(text), INTERCEPT.search(text)
        input_match, regression_match = INPUT.search(text), REGRESSION.search(text)
        if not all((h2_match, intercept_match, input_match, regression_match)):
            raise SystemExit(f"ERROR: required h2 diagnostics missing: {log}")
        assert h2_match and intercept_match and input_match and regression_match
        h2, se = float(h2_match.group(2)), float(h2_match.group(3))
        intercept, intercept_se = float(intercept_match.group(1)), float(intercept_match.group(2))
        z = h2 / se if se > 0 else math.nan
        passes = math.isfinite(z) and z >= 4 and math.isfinite(intercept) and intercept <= 1.2
        reason = "pass" if passes else "h2_z_below_4_or_missing" if not math.isfinite(z) or z < 4 else "intercept_above_1.2_or_missing"
        output.append({
            "replication_source_id": source_id, "study_accession": source["replication_study_accession"],
            "phenotype_definition": source["replication_phenotype_definition"], "ancestry": source["ancestry"],
            "scale": h2_match.group(1).lower(), "h2": h2, "h2_se": se, "h2_z": z,
            "LDSC_intercept": intercept, "LDSC_intercept_se": intercept_se,
            "input_snp_count": int(input_match.group(1)), "ldsc_regression_snp_count": int(regression_match.group(1)),
            "lambda_gc": numeric(LAMBDA.search(text)), "mean_chi2": numeric(MEAN_CHI.search(text)),
            "attenuation_ratio": numeric(RATIO.search(text)) if RATIO.search(text) else "LT_ZERO_OR_NOT_ESTIMATED",
            "primary_status": "PASS" if passes else "FAIL", "qc_reason": reason,
            "analysis_status": "REPLICATION_H2_QC_COMPLETE", "input_log": str(log),
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(output[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0", "manifest_sha256": sha256(args.manifest), "lock_sha256": sha256(args.lock),
        "source_count": len(output), "pass_count": sum(row["primary_status"] == "PASS" for row in output),
        "fail_count": sum(row["primary_status"] == "FAIL" for row in output),
        "thresholds": {"h2_z_min": 4.0, "ldsc_intercept_max": 1.2},
        "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"REPLICATION_H2_COLLATED sources={len(output)} pass={provenance['pass_count']} fail={provenance['fail_count']} sha256={provenance['output_sha256']}")


if __name__ == "__main__":
    main()
