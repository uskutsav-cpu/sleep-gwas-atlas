#!/usr/bin/env python3
"""Validate normalized LAVA/HDL-L rows and emit canonical local evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


RAW_FIELDS = {
    "pair_id", "method_id", "locus_id", "chr", "start", "end",
    "local_h2_trait_a", "local_h2_trait_a_se", "local_h2_trait_a_p",
    "local_h2_trait_b", "local_h2_trait_b_se", "local_h2_trait_b_p",
    "local_rg", "local_rg_se", "local_rg_p", "bivariate_test_status",
    "effect_direction", "convergence_status", "qc_status",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(pairs: list[tuple[int, float]]) -> dict[int, float]:
    ordered = sorted(pairs, key=lambda item: item[1])
    count = len(ordered)
    adjusted: dict[int, float] = {}
    running = 1.0
    for rank, (index, p_value) in reversed(list(enumerate(ordered, start=1))):
        running = min(running, p_value * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def finite_probability(value: str, field: str, identity: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid {field} for {identity}") from exc
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        raise SystemExit(f"ERROR: invalid {field} for {identity}")
    return parsed


def finite_number(value: str, field: str, identity: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid {field} for {identity}") from exc
    if not math.isfinite(parsed):
        raise SystemExit(f"ERROR: invalid {field} for {identity}")
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/local_analysis_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/local_analysis_manifest.lock.json"))
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/local/local_rg_results.tsv"))
    parser.add_argument("--pair-summary-out", type=Path, default=Path("discovery_extension/results/local/local_pair_summary.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/local_architecture_results.json"))
    args = parser.parse_args()

    manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    if sha256(args.manifest) != lock.get("manifest_sha256"):
        raise SystemExit("ERROR: local manifest differs from its pre-result lock")
    if [row["pair_id"] for row in manifest] != lock.get("pair_ids"):
        raise SystemExit("ERROR: local manifest membership/order differs from lock")
    manifest_by_pair = {row["pair_id"]: row for row in manifest}
    raw = read_tsv(args.raw)
    if not raw or not RAW_FIELDS.issubset(raw[0]):
        raise SystemExit(f"ERROR: normalized local table lacks fields: {sorted(RAW_FIELDS - set(raw[0] if raw else []))}")
    if any(row["pair_id"] not in manifest_by_pair for row in raw):
        raise SystemExit("ERROR: normalized local result contains a pair outside the lock")
    keys = [(row["pair_id"], row["method_id"], row["locus_id"]) for row in raw]
    if len(set(keys)) != len(keys):
        raise SystemExit("ERROR: duplicate pair/method/locus local result")

    rows_by_pair_method: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        rows_by_pair_method[(row["pair_id"], row["method_id"])].append(row)
    for pair_id, item in manifest_by_pair.items():
        expected_loci = int(item["locus_count"])
        for path_field, checksum_field in (
            ("trait_a_dense_path", "trait_a_dense_sha256"),
            ("trait_b_dense_path", "trait_b_dense_sha256"),
            ("locus_definition_path", "locus_definition_sha256"),
        ):
            locked_input = Path(item[path_field])
            if not locked_input.is_file() or sha256(locked_input) != item[checksum_field]:
                raise SystemExit(f"ERROR: locked local input drifted for {pair_id}: {path_field}")
        locus_rows = read_tsv(Path(item["locus_definition_path"]))
        required_locus_fields = {"locus_id", "chr", "start", "end"}
        if not locus_rows or not required_locus_fields.issubset(locus_rows[0]):
            raise SystemExit(f"ERROR: locked locus definition lacks required fields for {pair_id}")
        if len(locus_rows) != expected_loci or len({row["locus_id"] for row in locus_rows}) != expected_loci:
            raise SystemExit(f"ERROR: locked locus definition count/identity differs for {pair_id}")
        loci = {row["locus_id"]: row for row in locus_rows}
        for method in item["planned_methods"].split(";"):
            observed = rows_by_pair_method[(pair_id, method)]
            if len(observed) != expected_loci:
                raise SystemExit(f"ERROR: {pair_id}/{method} has {len(observed)} rows; expected {expected_loci}")
            if {row["locus_id"] for row in observed} != set(loci):
                raise SystemExit(f"ERROR: {pair_id}/{method} locus identities differ from the locked definition")
            for row in observed:
                locus = loci[row["locus_id"]]
                if row["chr"] != locus["chr"] or int(row["start"]) != int(locus["start"]) or int(row["end"]) != int(locus["end"]):
                    raise SystemExit(f"ERROR: {pair_id}/{method}/{row['locus_id']} coordinates differ from lock")
        unexpected = {method for pair, method in rows_by_pair_method if pair == pair_id} - set(item["planned_methods"].split(";"))
        if unexpected:
            raise SystemExit(f"ERROR: unplanned method for {pair_id}: {sorted(unexpected)}")

    h2_alpha = float(lock["lava_local_h2_bonferroni_alpha"])
    expected_h2_alpha = 0.05 / (2 * int(lock["lava_pair_locus_family_size"]))
    if abs(h2_alpha - expected_h2_alpha) > 1e-15:
        raise SystemExit("ERROR: locked LAVA local-h2 threshold differs from the exact Bonferroni family")
    eligible: list[tuple[int, float]] = []
    for index, row in enumerate(raw):
        identity = "/".join((row["pair_id"], row["method_id"], row["locus_id"]))
        if row["effect_direction"] not in {"POSITIVE", "NEGATIVE", "ZERO", "NOT_ESTIMATED"}:
            raise SystemExit(f"ERROR: invalid effect_direction for {identity}")
        for trait in ("a", "b"):
            finite_number(row[f"local_h2_trait_{trait}"], f"local_h2_trait_{trait}", identity)
            h2_se = finite_number(row[f"local_h2_trait_{trait}_se"], f"local_h2_trait_{trait}_se", identity)
            finite_probability(row[f"local_h2_trait_{trait}_p"], f"local_h2_trait_{trait}_p", identity)
            if h2_se <= 0:
                raise SystemExit(f"ERROR: non-positive local h2 SE for {identity}")
        both_lava_h2_pass = (
            row["method_id"] == "LAVA_PRIMARY"
            and float(row["local_h2_trait_a_p"]) < h2_alpha
            and float(row["local_h2_trait_b_p"]) < h2_alpha
        )
        row["local_h2_gate_alpha"] = f"{h2_alpha:.12g}" if row["method_id"] == "LAVA_PRIMARY" else "METHOD_SPECIFIC"
        row["local_h2_gate_status"] = (
            "BOTH_PASS" if both_lava_h2_pass else
            ("ONE_OR_BOTH_FAIL" if row["method_id"] == "LAVA_PRIMARY" else "METHOD_SPECIFIC_SENSITIVITY")
        )
        if row["method_id"] == "LAVA_PRIMARY":
            if both_lava_h2_pass and row["bivariate_test_status"] == "NOT_TESTED_LOCAL_H2_GATE":
                raise SystemExit(f"ERROR: LAVA bivariate gate incorrectly excluded {identity}")
            if not both_lava_h2_pass and row["bivariate_test_status"] != "NOT_TESTED_LOCAL_H2_GATE":
                raise SystemExit(f"ERROR: LAVA bivariate test bypassed the locked local-h2 gate: {identity}")
        if row["bivariate_test_status"] == "TESTED":
            finite_probability(row["local_rg_p"], "local_rg_p", identity)
            rg = finite_number(row["local_rg"], "local_rg", identity)
            se = finite_number(row["local_rg_se"], "local_rg_se", identity)
            if se <= 0:
                raise SystemExit(f"ERROR: invalid local rg/SE for {identity}")
            expected_direction = "POSITIVE" if rg > 0 else ("NEGATIVE" if rg < 0 else "ZERO")
            if row["effect_direction"] != expected_direction:
                raise SystemExit(f"ERROR: local effect direction disagrees with estimate for {identity}")
            if row["method_id"] == "LAVA_PRIMARY" and row["convergence_status"] == "PASS" and row["qc_status"] == "PASS":
                eligible.append((index, float(row["local_rg_p"])))
        elif row["bivariate_test_status"] not in {"NOT_TESTED_LOCAL_H2_GATE", "FAILED_ESTIMATION"}:
            raise SystemExit(f"ERROR: invalid bivariate_test_status for {identity}")
    adjusted = bh(eligible)

    significant_by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    for index, row in enumerate(raw):
        row["local_rg_fdr"] = f"{adjusted[index]:.12g}" if index in adjusted else "NA"
        row["primary_fdr_family_member"] = str(index in adjusted)
        row["local_significance_status"] = "FDR_SIGNIFICANT" if index in adjusted and adjusted[index] < 0.05 else (
            "TESTED_NOT_FDR_SIGNIFICANT" if index in adjusted else "NOT_IN_PRIMARY_FDR_FAMILY"
        )
        row["architecture_flags"] = "NONE"
        if row["local_significance_status"] == "FDR_SIGNIFICANT":
            significant_by_pair[row["pair_id"]].append(row)

    pair_summaries: list[dict[str, object]] = []
    for pair_id, item in manifest_by_pair.items():
        significant = significant_by_pair[pair_id]
        lava_sig = [row for row in significant if row["method_id"] == "LAVA_PRIMARY"]
        positive = [row for row in lava_sig if float(row["local_rg"]) > 0]
        negative = [row for row in lava_sig if float(row["local_rg"]) < 0]
        global_null = float(item["global_rg_fdr"]) >= 0.05
        flags: list[str] = []
        if global_null and positive:
            flags.append("GLOBAL_NULL_LOCAL_POSITIVE")
        if global_null and negative:
            flags.append("GLOBAL_NULL_LOCAL_NEGATIVE")
        if positive and negative:
            flags.append("OPPOSING_LOCAL_EFFECTS")
        if len({row["locus_id"] for row in lava_sig}) >= 2:
            flags.append("MULTI_LOCUS_SHARED_ARCHITECTURE")
        hdl_clean = {
            (row["locus_id"], row["effect_direction"]) for row in raw
            if row["pair_id"] == pair_id and row["method_id"] == "HDL_L_ROBUSTNESS"
            and row["bivariate_test_status"] == "TESTED" and row["convergence_status"] == "PASS"
            and row["qc_status"] == "PASS" and float(row["local_rg_p"]) < 0.05
        }
        concordant = sum((row["locus_id"], row["effect_direction"]) in hdl_clean for row in lava_sig)
        if len(lava_sig) >= 2 or (lava_sig and concordant >= 1):
            support = "STRONG_LOCAL_SUPPORT"
        elif lava_sig:
            support = "LOCAL_SUPPORT_LAVA_ONLY"
        else:
            support = "NO_FDR_LOCAL_SUPPORT"
        joined_flags = ";".join(flags) if flags else "NONE"
        for row in raw:
            if row["pair_id"] == pair_id:
                row["architecture_flags"] = joined_flags
        pair_summaries.append({
            "pair_id": pair_id, "sleep_trait": item["sleep_trait"],
            "extension_trait_id": item["extension_trait_id"], "selection_stratum": item["selection_stratum"],
            "global_rg": item["global_rg"], "global_rg_fdr": item["global_rg_fdr"],
            "locked_locus_count": item["locus_count"], "lava_tested_loci": sum(
                row["pair_id"] == pair_id and row["method_id"] == "LAVA_PRIMARY" and row["bivariate_test_status"] == "TESTED"
                for row in raw
            ),
            "lava_fdr_significant_loci": len(lava_sig), "lava_positive_loci": len(positive),
            "lava_negative_loci": len(negative), "hdl_directionally_concordant_loci": concordant,
            "architecture_flags": joined_flags, "local_support_status": support,
            "claim_limit": "LOCAL_GENETIC_SHARING_NOT_SHARED_CAUSAL_VARIANT_OR_CAUSALITY",
        })

    output_fields = list(raw[0])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=output_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(raw)
    summary_fields = list(pair_summaries[0])
    with args.pair_summary_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(pair_summaries)
    provenance = {
        "schema_version": "1.0.0", "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifest_sha256": sha256(args.manifest), "manifest_lock_sha256": sha256(args.lock),
        "normalized_raw_sha256": sha256(args.raw), "primary_method": "LAVA_PRIMARY",
        "primary_fdr_method": "Benjamini-Hochberg", "primary_fdr_threshold": 0.05,
        "lava_local_h2_bonferroni_alpha": h2_alpha,
        "primary_fdr_denominator": len(eligible), "locked_pair_locus_family_size": lock["lava_pair_locus_family_size"],
        "all_locked_tests_preserved": True, "pair_count": len(pair_summaries),
        "local_support_counts": dict(sorted(Counter(str(row["local_support_status"]) for row in pair_summaries).items())),
        "output_sha256": sha256(args.out), "pair_summary_sha256": sha256(args.pair_summary_out),
        "warning": "Local genetic sharing does not establish a shared causal variant or causality.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"LOCAL_ARCHITECTURE_RESULTS_OK pairs={len(pair_summaries)} rows={len(raw)} primary_fdr_n={len(eligible)}")


if __name__ == "__main__":
    main()
