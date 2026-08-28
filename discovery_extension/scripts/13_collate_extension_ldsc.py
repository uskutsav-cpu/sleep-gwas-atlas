#!/usr/bin/env python3
"""Collate the exact extension LDSC h2 and rg families with isolated BH FDR."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import re
from pathlib import Path


H2_PATTERN = re.compile(r"Total (Observed|Liability) scale h2:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)")
INTERCEPT_PATTERN = re.compile(r"^Intercept:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)", re.M)
MEAN_CHI_PATTERN = re.compile(r"Mean Chi\^2:\s*(-?[\d.eE+-]+)")
LAMBDA_PATTERN = re.compile(r"Lambda GC:\s*(-?[\d.eE+-]+)")
RATIO_PATTERN = re.compile(r"^Ratio:\s*(-?[\d.eE+-]+)", re.M)
RATIO_NEGATIVE_PATTERN = re.compile(r"^Ratio < 0", re.M)
INPUT_SNP_PATTERN = re.compile(r"Read summary statistics for ([0-9]+) SNPs\.")
REGRESSION_SNP_PATTERN = re.compile(r"After merging with regression SNP LD, ([0-9]+) SNPs remain\.")
PAIR_QC_PATTERN = re.compile(
    r"Computing rg for phenotype [0-9]+/[0-9]+\s*\n"
    r"Reading summary statistics from (.+?) \.\.\.\s*\n"
    r"Read summary statistics for ([0-9]+) SNPs\.\s*\n"
    r"After merging with summary statistics, ([0-9]+) SNPs remain\.\s*\n"
    r"([0-9]+) SNPs with valid alleles\.",
    re.M,
)
SYNTHETIC_MARKER = "SYNTHETIC EXTENSION LDSC OUTPUT - NOT REAL RESULTS"


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric(value: str | None) -> float:
    try:
        return float(value) if value is not None else math.nan
    except ValueError:
        return math.nan


def reject_synthetic(logdir: Path) -> None:
    for path_text in glob.glob(str(logdir / "*.log")):
        path = Path(path_text)
        in_synthetic_namespace = any("synthetic" in part.lower() for part in path.parts)
        if SYNTHETIC_MARKER in path.read_text(encoding="utf-8") and not in_synthetic_namespace:
            fail(f"synthetic LDSC log found outside a synthetic namespace: {path}")


def metadata(panel_path: Path, core_path: Path) -> tuple[list[dict[str, str]], list[str]]:
    panel = read_tsv(panel_path)
    if len(panel) != 100 or len({row["extension_trait_id"] for row in panel}) != 100:
        fail("extension panel is not the exact 100-trait lock")
    core = read_tsv(core_path)
    sleeps = [row["trait_id"] for row in core if row["domain"] == "sleep"]
    if len(sleeps) != 12:
        fail("core sleep family is not exactly 12 traits")
    return panel, sleeps


def parse_h2(panel: list[dict[str, str]], logdir: Path) -> list[dict[str, object]]:
    reject_synthetic(logdir)
    rows: list[dict[str, object]] = []
    for trait in panel:
        trait_id = trait["extension_trait_id"]
        path = logdir / f"h2_{trait_id}.log"
        if not path.is_file():
            fail(f"locked extension h2 log is missing: {path}")
        text = path.read_text(encoding="utf-8")
        h2_match = H2_PATTERN.search(text)
        intercept_match = INTERCEPT_PATTERN.search(text)
        if not h2_match or not intercept_match:
            fail(f"required h2/intercept estimate is missing from {path}")
        h2 = numeric(h2_match.group(2))
        se = numeric(h2_match.group(3))
        intercept = numeric(intercept_match.group(1))
        intercept_se = numeric(intercept_match.group(2))
        input_snp_match = INPUT_SNP_PATTERN.search(text)
        regression_snp_match = REGRESSION_SNP_PATTERN.search(text)
        if not input_snp_match or not regression_snp_match:
            fail(f"required h2 SNP-count diagnostics are missing from {path}")
        ratio_match = RATIO_PATTERN.search(text)
        ratio_negative = bool(RATIO_NEGATIVE_PATTERN.search(text))
        if not ratio_match and not ratio_negative:
            fail(f"required attenuation-ratio diagnostic is missing from {path}")
        z = h2 / se if se > 0 else math.nan
        pass_z = math.isfinite(z) and z >= 4.0
        pass_intercept = math.isfinite(intercept) and intercept <= 1.2
        if pass_z and pass_intercept:
            verdict, reason = "PRIMARY_PASS", "pass"
        elif not pass_z:
            verdict, reason = "SENSITIVITY_ONLY", "h2_z_below_4_or_missing"
        else:
            verdict, reason = "SENSITIVITY_ONLY", "intercept_above_1.2_or_missing"
        rows.append({
            "extension_trait_id": trait_id,
            "phenotype_name": trait["phenotype_name"],
            "phenotype_domain": trait["phenotype_domain"],
            "novelty_priority": trait["novelty_priority"],
            "scale": h2_match.group(1).lower(), "h2": h2, "h2_se": se, "h2_z": z,
            "LDSC_intercept": intercept, "LDSC_intercept_se": intercept_se,
            "input_snp_count": int(input_snp_match.group(1)),
            "ldsc_regression_snp_count": int(regression_snp_match.group(1)),
            "lambda_gc": numeric(LAMBDA_PATTERN.search(text).group(1)) if LAMBDA_PATTERN.search(text) else math.nan,
            "mean_chi2": numeric(MEAN_CHI_PATTERN.search(text).group(1)) if MEAN_CHI_PATTERN.search(text) else math.nan,
            "attenuation_ratio": numeric(ratio_match.group(1)) if ratio_match else "LT_ZERO",
            "attenuation_ratio_status": "NUMERIC" if ratio_match else "NEGATIVE_NOT_ESTIMATED",
            "pass_h2_z_ge_4": str(pass_z), "pass_intercept_le_1.2": str(pass_intercept),
            "primary_rg_eligibility": verdict, "analysis_status": "EXTENSION_H2_QC_COMPLETE",
            "qc_reason": reason, "input_log": str(path),
        })
    return rows


def parse_one_rg_log(path: Path) -> list[dict[str, object]]:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    pair_qc: dict[str, dict[str, int]] = {}
    for match in PAIR_QC_PATTERN.finditer(text):
        trait_id = Path(match.group(1)).name.replace(".sumstats.gz", "")
        if trait_id in pair_qc:
            fail(f"duplicate per-pair SNP diagnostics for {trait_id} in {path}")
        pair_qc[trait_id] = {
            "extension_input_snp_count": int(match.group(2)),
            "snp_overlap_after_merge": int(match.group(3)),
            "snp_overlap_valid_alleles": int(match.group(4)),
        }
    try:
        start = next(i for i, line in enumerate(lines) if "Summary of Genetic Correlation Results" in line)
    except StopIteration:
        fail(f"rg results table is missing from {path}")
    header = lines[start + 1].split()
    required = {"p1", "p2", "rg", "se", "p"}
    if not required.issubset(header):
        fail(f"unexpected rg header in {path}: {header}")
    body: list[list[str]] = []
    for line in lines[start + 2:]:
        if not line.strip():
            break
        fields = line.split()
        if len(fields) != len(header):
            break
        body.append(fields)
    scalar_p = [line.split(":", 1)[1].strip() for line in lines[:start] if line.startswith("P:")]
    rows = [dict(zip(header, fields)) for fields in body]
    if len(scalar_p) == len(rows):
        for row, p_value in zip(rows, scalar_p):
            row["p"] = p_value
    output: list[dict[str, object]] = []
    for row in rows:
        sleep = Path(row["p1"]).name.replace(".sumstats.gz", "")
        extension = Path(row["p2"]).name.replace(".sumstats.gz", "")
        if extension not in pair_qc:
            fail(f"per-pair SNP-overlap diagnostics are missing for {extension} in {path}")
        output.append({
            "sleep_trait": sleep, "extension_trait_id": extension,
            "rg": numeric(row["rg"]), "se": numeric(row["se"]),
            "z": numeric(row.get("z")), "p": numeric(row["p"]),
            "extension_h2_observed": numeric(row.get("h2_obs")),
            "extension_h2_observed_se": numeric(row.get("h2_obs_se")),
            "extension_h2_intercept": numeric(row.get("h2_int")),
            "extension_h2_intercept_se": numeric(row.get("h2_int_se")),
            "cross_trait_LDSC_intercept": numeric(row.get("gcov_int")),
            "cross_trait_LDSC_intercept_se": numeric(row.get("gcov_int_se")),
            **pair_qc[extension],
            "input_log": str(path),
        })
    return output


def bh_fdr(rows: list[dict[str, object]]) -> None:
    ordered = sorted(range(len(rows)), key=lambda index: float(rows[index]["p"]))
    m = len(ordered)
    adjusted = [math.nan] * m
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        original_index = ordered[rank_index]
        p_value = float(rows[original_index]["p"])
        running = min(running, p_value * m / (rank_index + 1))
        adjusted[original_index] = min(1.0, running)
    for row, value in zip(rows, adjusted):
        row["extension_fdr"] = value


def parse_rg(
    panel: list[dict[str, str]], sleeps: list[str], h2_rows: list[dict[str, str]],
    logdir: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    reject_synthetic(logdir)
    panel_by_id = {row["extension_trait_id"]: row for row in panel}
    h2_by_id = {row["extension_trait_id"]: row for row in h2_rows}
    if set(h2_by_id) != set(panel_by_id):
        fail("rerun h2 table is not the complete locked extension panel")
    passing = [trait_id for trait_id in panel_by_id if h2_by_id[trait_id]["primary_rg_eligibility"] == "PRIMARY_PASS"]
    expected = {(sleep, trait_id) for sleep in sleeps for trait_id in passing}
    primary: list[dict[str, object]] = []
    for sleep in sleeps:
        path = logdir / f"rg_{sleep}.log"
        if not path.is_file():
            fail(f"required extension rg log is missing: {path}")
        primary.extend(parse_one_rg_log(path))
    observed = {(str(row["sleep_trait"]), str(row["extension_trait_id"])) for row in primary}
    if len(observed) != len(primary):
        fail("duplicate extension rg pair detected")
    if observed != expected:
        fail(
            f"extension rg family mismatch expected={len(expected)} observed={len(observed)} "
            f"missing={sorted(expected-observed)[:5]} extra={sorted(observed-expected)[:5]}"
        )
    if any(not math.isfinite(float(row["p"])) or not 0 <= float(row["p"]) <= 1 for row in primary):
        fail("extension rg result contains invalid p-values")
    bh_fdr(primary)
    for row in primary:
        trait = panel_by_id[str(row["extension_trait_id"])]
        row.update({
            "phenotype_name": trait["phenotype_name"],
            "phenotype_domain": trait["phenotype_domain"],
            "novelty_priority": trait["novelty_priority"],
            "ancestry": trait["ancestry"],
            "analysis_status": "PRIMARY_EXTENSION_RG_COMPLETE",
            "extension_fdr_pass_0.05": str(float(row["extension_fdr"]) < 0.05),
            "abs_rg_ge_0.15": str(abs(float(row["rg"])) >= 0.15),
            "initial_screen_status": (
                "PAIR_NOVELTY_AUDIT_REQUIRED"
                if float(row["extension_fdr"]) < 0.05 and abs(float(row["rg"])) >= 0.15
                else "NOT_TIER_A_SCREEN_THRESHOLD"
            ),
        })
    primary_by_pair = {(str(row["sleep_trait"]), str(row["extension_trait_id"])): row for row in primary}
    universe: list[dict[str, object]] = []
    for sleep in sleeps:
        for trait_id, trait in panel_by_id.items():
            h2 = h2_by_id[trait_id]
            result = primary_by_pair.get((sleep, trait_id), {})
            universe.append({
                "sleep_trait": sleep, "extension_trait_id": trait_id,
                "phenotype_name": trait["phenotype_name"], "phenotype_domain": trait["phenotype_domain"],
                "novelty_priority": trait["novelty_priority"],
                "extension_h2_z": h2["h2_z"], "extension_LDSC_intercept": h2["LDSC_intercept"],
                "extension_h2_verdict": h2["primary_rg_eligibility"],
                "ancestry": trait["ancestry"],
                "pair_status": "PRIMARY_TESTED" if result else "H2_FAILED_PRIMARY_EXCLUSION",
                "rg": result.get("rg", "NA"), "se": result.get("se", "NA"),
                "z": result.get("z", "NA"), "p": result.get("p", "NA"),
                "cross_trait_LDSC_intercept": result.get("cross_trait_LDSC_intercept", "NA"),
                "cross_trait_LDSC_intercept_se": result.get("cross_trait_LDSC_intercept_se", "NA"),
                "snp_overlap_after_merge": result.get("snp_overlap_after_merge", "NA"),
                "snp_overlap_valid_alleles": result.get("snp_overlap_valid_alleles", "NA"),
                "extension_fdr": result.get("extension_fdr", "NA"),
                "initial_screen_status": result.get("initial_screen_status", "SENSITIVITY_ONLY_NOT_RUN"),
            })
    return primary, universe


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("h2", "rg"), required=True)
    parser.add_argument("--logdir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pair-universe-out", type=Path)
    parser.add_argument("--h2", type=Path)
    parser.add_argument("--h2-failed-out", type=Path)
    parser.add_argument(
        "--panel", type=Path, default=Path("discovery_extension/config/candidate_traits.tsv")
    )
    parser.add_argument("--core", type=Path, default=Path("config/analysis_panel.tsv"))
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/ldsc_collation.json"),
    )
    args = parser.parse_args()

    panel, sleeps = metadata(args.panel, args.core)
    if args.mode == "h2":
        rows = parse_h2(panel, args.logdir)
        fields = list(rows[0])
        write_tsv(args.out, fields, rows)
        if args.h2_failed_out is not None:
            failed = [row for row in rows if row["primary_rg_eligibility"] == "SENSITIVITY_ONLY"]
            write_tsv(args.h2_failed_out, fields, failed)
        counts = {
            "rows": len(rows),
            "primary_pass": sum(row["primary_rg_eligibility"] == "PRIMARY_PASS" for row in rows),
            "sensitivity_only": sum(row["primary_rg_eligibility"] == "SENSITIVITY_ONLY" for row in rows),
        }
    else:
        if args.h2 is None or args.pair_universe_out is None:
            fail("rg mode requires --h2 and --pair-universe-out")
        h2_rows = read_tsv(args.h2)
        rows, universe = parse_rg(panel, sleeps, h2_rows, args.logdir)
        write_tsv(args.out, list(rows[0]) if rows else [], rows)
        write_tsv(args.pair_universe_out, list(universe[0]), universe)
        counts = {"primary_pairs": len(rows), "locked_pair_universe": len(universe)}
    provenance = {
        "schema_version": "1.0.0", "mode": args.mode,
        "panel_sha256": sha256(args.panel), "core_panel_sha256": sha256(args.core),
        "logdir": str(args.logdir), "counts": counts,
        "multiple_testing": "Benjamini-Hochberg within primary extension rg family only" if args.mode == "rg" else "not_applicable",
        "output": str(args.out), "output_sha256": sha256(args.out),
    }
    if args.mode == "h2" and args.h2_failed_out is not None:
        provenance["h2_failed_output"] = str(args.h2_failed_out)
        provenance["h2_failed_output_sha256"] = sha256(args.h2_failed_out)
    if args.mode == "rg":
        provenance["pair_universe_output"] = str(args.pair_universe_out)
        provenance["pair_universe_sha256"] = sha256(args.pair_universe_out)
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"EXTENSION_{args.mode.upper()}_COLLATED {counts}")


if __name__ == "__main__":
    main()
