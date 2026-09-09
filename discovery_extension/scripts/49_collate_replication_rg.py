#!/usr/bin/env python3
"""Collate the exact h2-eligible independent-replication rg family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


PAIR_QC = re.compile(
    r"Computing rg for phenotype [0-9]+/[0-9]+\s*\n"
    r"Reading summary statistics from (.+?) \.\.\.\s*\n"
    r"Read summary statistics for ([0-9]+) SNPs\.\s*\n"
    r"After merging with summary statistics, ([0-9]+) SNPs remain\.\s*\n"
    r"([0-9]+) SNPs with valid alleles\.", re.M,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: str | None) -> float:
    try:
        return float(value) if value is not None else math.nan
    except ValueError:
        return math.nan


def parse_log(path: Path) -> list[dict[str, object]]:
    text = path.read_text(encoding="utf-8")
    pair_qc: dict[str, dict[str, int]] = {}
    for match in PAIR_QC.finditer(text):
        source_id = Path(match.group(1)).name.removesuffix(".sumstats.gz")
        pair_qc[source_id] = {
            "replication_input_snp_count": int(match.group(2)),
            "snp_overlap_after_merge": int(match.group(3)),
            "snp_overlap_valid_alleles": int(match.group(4)),
        }
    lines = text.splitlines()
    try:
        start = next(index for index, line in enumerate(lines) if "Summary of Genetic Correlation Results" in line)
    except StopIteration as error:
        raise SystemExit(f"ERROR: rg results table missing from {path}") from error
    header = lines[start + 1].split()
    if not {"p1", "p2", "rg", "se", "p"}.issubset(header):
        raise SystemExit(f"ERROR: unexpected rg header in {path}: {header}")
    rows: list[dict[str, str]] = []
    for line in lines[start + 2:]:
        if not line.strip():
            break
        values = line.split()
        if len(values) != len(header):
            break
        rows.append(dict(zip(header, values)))
    scalar_p = [line.split(":", 1)[1].strip() for line in lines[:start] if line.startswith("P:")]
    if len(scalar_p) == len(rows):
        for row, value in zip(rows, scalar_p):
            row["p"] = value
    output: list[dict[str, object]] = []
    for row in rows:
        sleep = Path(row["p1"]).name.removesuffix(".sumstats.gz")
        source_id = Path(row["p2"]).name.removesuffix(".sumstats.gz")
        if source_id not in pair_qc:
            raise SystemExit(f"ERROR: pair diagnostics missing for {source_id} in {path}")
        z_value = number(row.get("z"))
        p_value = math.erfc(abs(z_value) / math.sqrt(2)) if math.isfinite(z_value) else number(row["p"])
        output.append({
            "sleep_trait": sleep, "replication_source_id": source_id,
            "rg": number(row["rg"]), "se": number(row["se"]), "z": z_value, "p": p_value,
            "replication_h2_observed": number(row.get("h2_obs")),
            "replication_h2_observed_se": number(row.get("h2_obs_se")),
            "replication_h2_intercept": number(row.get("h2_int")),
            "replication_h2_intercept_se": number(row.get("h2_int_se")),
            "cross_trait_LDSC_intercept": number(row.get("gcov_int")),
            "cross_trait_LDSC_intercept_se": number(row.get("gcov_int_se")),
            **pair_qc[source_id], "input_log": str(path),
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/replication_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/replication_manifest.lock.json"))
    parser.add_argument("--h2", type=Path, default=Path("discovery_extension/results/replication/replication_source_h2.tsv"))
    parser.add_argument("--jobs", type=Path, default=Path("discovery_extension/results/replication/replication_rg_jobs.tsv"))
    parser.add_argument("--logdir", type=Path, default=Path("discovery_extension/logs/replication/rg"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/replication/replication_rg.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/replication_rg.json"))
    args = parser.parse_args()
    manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    if sha256(args.manifest) != lock.get("manifest_sha256"):
        raise SystemExit("ERROR: manifest differs from lock")
    h2 = {row["replication_source_id"]: row for row in read_tsv(args.h2)}
    testable = set(lock["testable_pair_ids_in_locked_order"])
    expected: dict[tuple[str, str], dict[str, str]] = {}
    for row in manifest:
        if row["pair_id"] in testable and h2[row["replication_source_id"]]["primary_status"] == "PASS":
            expected[(row["sleep_trait"], row["replication_source_id"])] = row
    observed: list[dict[str, object]] = []
    for job in read_tsv(args.jobs):
        observed.extend(parse_log(args.logdir / f"rg_replication_{job['sleep_trait']}.log"))
    observed_keys = {(str(row["sleep_trait"]), str(row["replication_source_id"])) for row in observed}
    if len(observed_keys) != len(observed) or observed_keys != set(expected):
        raise SystemExit(f"ERROR: replication rg family mismatch expected={len(expected)} observed={len(observed_keys)}")
    output: list[dict[str, object]] = []
    for result in observed:
        key = (str(result["sleep_trait"]), str(result["replication_source_id"]))
        source = expected[key]
        for field in ("rg", "se", "p"):
            if not math.isfinite(float(result[field])):
                raise SystemExit(f"ERROR: nonfinite replication {field}: {source['pair_id']}")
        result.update({
            "pair_id": source["pair_id"], "extension_trait_id": source["extension_trait_id"],
            "ancestry": source["ancestry"], "analysis_status": "REPLICATION_RG_COMPLETE",
        })
        output.append(result)
    order = {pair_id: index for index, pair_id in enumerate(lock["pair_ids_in_locked_order"])}
    output.sort(key=lambda row: order[str(row["pair_id"])])
    fields = [
        "pair_id", "sleep_trait", "extension_trait_id", "replication_source_id", "rg", "se", "z", "p",
        "replication_h2_observed", "replication_h2_observed_se", "replication_h2_intercept",
        "replication_h2_intercept_se", "cross_trait_LDSC_intercept", "cross_trait_LDSC_intercept_se",
        "replication_input_snp_count", "snp_overlap_after_merge", "snp_overlap_valid_alleles",
        "ancestry", "analysis_status", "input_log",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0", "manifest_sha256": sha256(args.manifest), "lock_sha256": sha256(args.lock),
        "h2_sha256": sha256(args.h2), "jobs_sha256": sha256(args.jobs), "pair_count": len(output),
        "multiplicity": "No replication-side FDR; success uses the frozen 0.05/217 Bonferroni threshold.",
        "p_value_precision": "Two-sided normal-tail p recomputed from the LDSC log z statistic because the LDSC scalar display rounds some p values to four decimals.",
        "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"REPLICATION_RG_COLLATED pairs={len(output)} sha256={provenance['output_sha256']}")


if __name__ == "__main__":
    main()
