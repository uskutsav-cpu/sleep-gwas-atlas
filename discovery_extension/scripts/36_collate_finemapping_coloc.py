#!/usr/bin/env python3
"""Collate the exact locked fine-mapping family and classify prior-robust coloc evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


VARIANT_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role",
    "dataset_type", "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO",
    "prior_method", "normalized_prior_weight", "PIP", "credible_set_ids", "max_alpha_component",
    "model_converged", "rss_ld_s", "kriging_allele_switch_outlier",
]
CREDIBLE_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role",
    "signal_id", "component_index", "lead_snp", "lead_pip", "credible_set_size",
    "credible_set_snps", "requested_coverage", "achieved_coverage", "min_abs_corr",
    "mean_abs_corr", "median_abs_corr", "cs_log10bf", "model_converged",
]
COLOC_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset1_id", "dataset2_id",
    "molecular_feature_id", "tissue_cell_context", "coloc_method", "p1", "p2", "p12",
    "prior_role", "signal1", "signal2", "hit1", "hit2", "nsnps", "PP_H0", "PP_H1",
    "PP_H2", "PP_H3", "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant",
    "top_shared_variant_PP_H4", "fine_mapping_qc", "analysis_status",
    "single_signal_fallback_justification", "claim_limit",
]
SHARED_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "coloc_method", "p12",
    "signal1", "signal2", "SNP", "SNP_PP_H4",
]
DIAGNOSTIC_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role",
    "variant_count", "model_converged", "niter", "credible_set_count", "max_pip", "rss_ld_s",
    "kriging_allele_switch_outlier_count", "kriging_allele_switch_outliers", "diagnostic_status",
]
CANONICAL_COLOC_FIELDS = COLOC_FIELDS + [
    "source_search_status", "primary_shared_signal_rule_pass", "prior_robust",
    "colocalization_interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def collect(paths: list[Path], expected: list[str], label: str) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for path in paths:
        fields, rows = read_tsv(path)
        if fields != expected:
            raise SystemExit(f"ERROR: {label} fields differ from contract: {path}")
        output.extend(rows)
    return output


def probability(value: str, field: str, identity: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid {field}: {identity}") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise SystemExit(f"ERROR: invalid {field}: {identity}")
    return number


def verify_manifest_inputs(row: dict[str, str]) -> None:
    if row["source_search_status"] == "NO_SUITABLE_DATASET":
        path = Path(row["source_search_log_path"])
        if not path.is_file() or sha256(path) != row["source_search_log_sha256"]:
            raise SystemExit(f"ERROR: unavailable-source evidence drifted: {row['queue_row_id']}")
        return
    for path_field, checksum_field in (
        ("summary1_path", "summary1_sha256"), ("summary2_path", "summary2_sha256"),
        ("ld_path", "ld_sha256"), ("ld_variant_order_path", "ld_variant_order_sha256"),
    ):
        path = Path(row[path_field])
        if not path.is_file() or sha256(path) != row[checksum_field]:
            raise SystemExit(f"ERROR: locked fine-mapping input drifted: {row['queue_row_id']}/{path_field}")
    for prefix in ("dataset1", "dataset2"):
        if row[f"{prefix}_prior_method"] == "POLYFUN_FUNCTIONAL":
            path = Path(row[f"{prefix}_prior_source_path"])
            if not path.is_file() or sha256(path) != row[f"{prefix}_prior_source_sha256"]:
                raise SystemExit(f"ERROR: locked functional prior source drifted: {row['queue_row_id']}/{prefix}")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/fine_mapping_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/fine_mapping_manifest.lock.json"))
    parser.add_argument("--variant-results", type=Path, nargs="+", required=True)
    parser.add_argument("--credible-set-results", type=Path, nargs="+", required=True)
    parser.add_argument("--coloc-results", type=Path, nargs="+", required=True)
    parser.add_argument("--shared-variant-results", type=Path, nargs="+", required=True)
    parser.add_argument("--diagnostic-results", type=Path, nargs="+", required=True)
    parser.add_argument("--variant-out", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_variants.tsv"))
    parser.add_argument("--credible-set-out", type=Path, default=Path("discovery_extension/results/fine_mapping/credible_sets.tsv"))
    parser.add_argument("--shared-variant-out", type=Path, default=Path("discovery_extension/results/fine_mapping/coloc_shared_variant_posteriors.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_colocalization.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/fine_mapping_colocalization_results.json"))
    args = parser.parse_args()

    _, manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    if sha256(args.manifest) != lock.get("manifest_sha256"):
        raise SystemExit("ERROR: fine-mapping manifest differs from its pre-result lock")
    manifest_by_id = {row["queue_row_id"]: row for row in manifest}
    if len(manifest_by_id) != len(manifest):
        raise SystemExit("ERROR: duplicate comparison identifier in fine-mapping manifest")
    for row in manifest:
        verify_manifest_inputs(row)
    analysis_ids = lock["analysis_comparison_ids_in_locked_order"]
    unavailable_ids = lock["unavailable_search_ids_in_locked_order"]
    if [row["queue_row_id"] for row in manifest if row["source_search_status"] == "VERIFIED_ANALYSIS"] != analysis_ids:
        raise SystemExit("ERROR: analyzable comparison family/order differs from lock")
    if [row["queue_row_id"] for row in manifest if row["source_search_status"] == "NO_SUITABLE_DATASET"] != unavailable_ids:
        raise SystemExit("ERROR: unavailable QTL-search family/order differs from lock")

    variants = collect(args.variant_results, VARIANT_FIELDS, "fine-mapping variant")
    credible_sets = collect(args.credible_set_results, CREDIBLE_FIELDS, "credible-set")
    coloc = collect(args.coloc_results, COLOC_FIELDS, "coloc")
    shared = collect(args.shared_variant_results, SHARED_FIELDS, "shared-variant posterior")
    diagnostics = collect(args.diagnostic_results, DIAGNOSTIC_FIELDS, "fine-mapping diagnostic")

    diagnostic_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in diagnostics:
        if row["comparison_id"] not in analysis_ids:
            raise SystemExit(f"ERROR: diagnostic outside locked family: {row['comparison_id']}")
        diagnostic_group[row["comparison_id"]].append(row)
    if set(diagnostic_group) != set(analysis_ids):
        raise SystemExit("ERROR: diagnostics do not cover the exact locked analysis family")
    failed_ids: set[str] = set()
    for comparison_id in analysis_ids:
        rows = diagnostic_group[comparison_id]
        manifest_row = manifest_by_id[comparison_id]
        if len(rows) != 2 or {row["dataset_id"] for row in rows} != {manifest_row["dataset1_id"], manifest_row["dataset2_id"]}:
            raise SystemExit(f"ERROR: {comparison_id} diagnostics do not cover both datasets exactly")
        for row in rows:
            if int(row["variant_count"]) != int(lock["variant_counts"][comparison_id]):
                raise SystemExit(f"ERROR: {comparison_id} diagnostic variant count differs from lock")
            s_value = float(row["rss_ld_s"])
            if not math.isfinite(s_value) or not 0 <= s_value <= 1 or int(row["kriging_allele_switch_outlier_count"]) < 0:
                raise SystemExit(f"ERROR: invalid LD-summary diagnostic: {comparison_id}/{row['dataset_id']}")
            if not row["diagnostic_status"].startswith(("PASS", "LD_SUMSTAT_INCONSISTENCY_FLAGGED")):
                failed_ids.add(comparison_id)

    variant_group: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in variants:
        if row["comparison_id"] not in analysis_ids:
            raise SystemExit(f"ERROR: variant result outside lock: {row['comparison_id']}")
        variant_group[(row["comparison_id"], row["dataset_id"])].append(row)
        probability(row["PIP"], "PIP", f"{row['comparison_id']}/{row['dataset_id']}/{row['SNP']}")
        probability(row["normalized_prior_weight"], "normalized_prior_weight", f"{row['comparison_id']}/{row['dataset_id']}/{row['SNP']}")
    for comparison_id in analysis_ids:
        manifest_row = manifest_by_id[comparison_id]
        for dataset_id in (manifest_row["dataset1_id"], manifest_row["dataset2_id"]):
            rows = variant_group.get((comparison_id, dataset_id), [])
            if comparison_id in failed_ids:
                if rows:
                    raise SystemExit(f"ERROR: failed fine map unexpectedly has variant results: {comparison_id}")
                continue
            expected = int(lock["variant_counts"][comparison_id])
            if len(rows) != expected or len({row["SNP"] for row in rows}) != expected:
                raise SystemExit(f"ERROR: {comparison_id}/{dataset_id} variant PIPs do not cover locked universe")
            if not math.isclose(sum(float(row["normalized_prior_weight"]) for row in rows), 1, rel_tol=0, abs_tol=1e-8):
                raise SystemExit(f"ERROR: {comparison_id}/{dataset_id} normalized prior weights do not sum to one")

    cs_keys: set[tuple[str, str, str]] = set()
    for row in credible_sets:
        comparison_id = row["comparison_id"]
        if comparison_id not in analysis_ids or comparison_id in failed_ids:
            raise SystemExit(f"ERROR: credible set outside successful locked family: {comparison_id}")
        key = (comparison_id, row["dataset_id"], row["signal_id"])
        if key in cs_keys:
            raise SystemExit(f"ERROR: duplicate credible set: {key}")
        cs_keys.add(key)
        if int(row["credible_set_size"]) != len(row["credible_set_snps"].split(";")):
            raise SystemExit(f"ERROR: credible-set size/list mismatch: {key}")
        probability(row["lead_pip"], "lead_pip", "/".join(key))
        if not 0 < float(row["achieved_coverage"]) <= 1 or float(row["min_abs_corr"]) < 0:
            raise SystemExit(f"ERROR: invalid credible-set coverage/purity: {key}")

    p12_grid = [float(value) for value in lock["coloc_p12_sensitivity_grid"]]
    primary_p12 = float(lock["coloc_primary_priors"]["p12"])
    coloc_keys: set[tuple[str, str, str, str, float]] = set()
    coloc_group: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in coloc:
        comparison_id = row["comparison_id"]
        if comparison_id not in analysis_ids or comparison_id in failed_ids:
            raise SystemExit(f"ERROR: coloc row outside successful locked family: {comparison_id}")
        p12 = float(row["p12"])
        if p12 not in p12_grid:
            raise SystemExit(f"ERROR: coloc p12 outside locked grid: {comparison_id}/{p12}")
        key = (comparison_id, row["coloc_method"], row["signal1"], row["signal2"], p12)
        if key in coloc_keys:
            raise SystemExit(f"ERROR: duplicate coloc signal/prior row: {key}")
        coloc_keys.add(key)
        group_key = key[:4]
        coloc_group[group_key].append(row)
        if row["analysis_status"].startswith(("COLOC_SUSIE_COMPLETE", "COLOC_ABF_FALLBACK_COMPLETE")):
            posterior_sum = sum(probability(row[field], field, str(key)) for field in ("PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4"))
            if not math.isclose(posterior_sum, 1, rel_tol=0, abs_tol=1e-5):
                raise SystemExit(f"ERROR: coloc hypotheses do not sum to one: {key}")
    for key, rows in coloc_group.items():
        if sorted(float(row["p12"]) for row in rows) != sorted(p12_grid):
            raise SystemExit(f"ERROR: coloc signal pair lacks the complete prior grid: {key}")
    for comparison_id in set(analysis_ids) - failed_ids:
        if not any(row["comparison_id"] == comparison_id for row in coloc):
            raise SystemExit(f"ERROR: successful fine map lacks coloc disposition: {comparison_id}")

    shared_group: dict[tuple[str, str, str, str, float], list[dict[str, str]]] = defaultdict(list)
    for row in shared:
        key = (row["comparison_id"], row["coloc_method"], row["signal1"], row["signal2"], float(row["p12"]))
        if key not in coloc_keys:
            raise SystemExit(f"ERROR: shared-variant posterior lacks matching coloc row: {key}")
        probability(row["SNP_PP_H4"], "SNP_PP_H4", f"{key}/{row['SNP']}")
        shared_group[key].append(row)
    for key, rows in shared_group.items():
        if len({row["SNP"] for row in rows}) != len(rows) or not math.isclose(sum(float(row["SNP_PP_H4"]) for row in rows), 1, rel_tol=0, abs_tol=1e-5):
            raise SystemExit(f"ERROR: shared-variant posterior is duplicated or does not sum to one: {key}")

    rule_pass: dict[tuple[str, str, str, str, float], bool] = {}
    for row in coloc:
        key = (row["comparison_id"], row["coloc_method"], row["signal1"], row["signal2"], float(row["p12"]))
        complete = row["analysis_status"] in {"COLOC_SUSIE_COMPLETE", "COLOC_ABF_FALLBACK_COMPLETE"}
        pass_qc = row["fine_mapping_qc"] == "PASS"
        if complete:
            pp4, pp3 = float(row["PP_H4"]), float(row["PP_H3"])
            ratio = math.inf if pp3 == 0 and pp4 > 0 else (pp4 / pp3 if pp3 > 0 else math.nan)
            rule_pass[key] = pass_qc and pp4 >= 0.80 and ratio >= 5
        else:
            rule_pass[key] = False
    robust_by_signal = {
        key: all(rule_pass[(key[0], key[1], key[2], key[3], p12)] for p12 in p12_grid)
        for key in coloc_group
    }

    canonical: list[dict[str, object]] = []
    coloc_by_comparison: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in coloc:
        coloc_by_comparison[row["comparison_id"]].append(row)
    for manifest_row in manifest:
        comparison_id = manifest_row["queue_row_id"]
        if manifest_row["source_search_status"] == "NO_SUITABLE_DATASET":
            record = {field: "NA" for field in COLOC_FIELDS}
            record.update({
                "comparison_id": comparison_id, "pair_id": manifest_row["pair_id"],
                "locus_id": manifest_row["locus_id"], "comparison_type": manifest_row["comparison_type"],
                "dataset1_id": manifest_row["dataset1_id"], "dataset2_id": manifest_row["dataset2_id"],
                "molecular_feature_id": manifest_row["molecular_feature_id"],
                "tissue_cell_context": manifest_row["tissue_cell_context"], "coloc_method": "NONE",
                "analysis_status": "NO_SUITABLE_DATASET", "claim_limit": lock["claim_limit"],
                "source_search_status": "NO_SUITABLE_DATASET", "primary_shared_signal_rule_pass": "False",
                "prior_robust": "False", "colocalization_interpretation": "NO_SUITABLE_DATASET",
            })
            canonical.append(record)
            continue
        if comparison_id in failed_ids:
            record = {field: "NA" for field in COLOC_FIELDS}
            record.update({
                "comparison_id": comparison_id, "pair_id": manifest_row["pair_id"],
                "locus_id": manifest_row["locus_id"], "comparison_type": manifest_row["comparison_type"],
                "dataset1_id": manifest_row["dataset1_id"], "dataset2_id": manifest_row["dataset2_id"],
                "molecular_feature_id": manifest_row["molecular_feature_id"],
                "tissue_cell_context": manifest_row["tissue_cell_context"], "coloc_method": "NONE",
                "analysis_status": "FINE_MAPPING_FAILED_OR_NONCONVERGED", "claim_limit": lock["claim_limit"],
                "source_search_status": "VERIFIED_ANALYSIS", "primary_shared_signal_rule_pass": "False",
                "prior_robust": "False", "colocalization_interpretation": "FINE_MAPPING_FAILED",
            })
            canonical.append(record)
            continue
        rows = sorted(coloc_by_comparison[comparison_id], key=lambda row: (float(row["p12"]), row["coloc_method"], row["signal1"], row["signal2"]))
        for row in rows:
            record: dict[str, object] = dict(row)
            key = (comparison_id, row["coloc_method"], row["signal1"], row["signal2"], float(row["p12"]))
            signal_key = key[:4]
            passed = rule_pass[key]
            if row["fine_mapping_qc"] != "PASS":
                interpretation = "QC_FLAGGED_INFERENCE_WITHHELD"
            elif row["analysis_status"] == "NO_CREDIBLE_SET_AND_SINGLE_SIGNAL_FALLBACK_NOT_JUSTIFIED":
                interpretation = "NO_ESTIMABLE_COLOCALIZATION"
            elif row["analysis_status"].startswith(("COLOC_SUSIE_COMPLETE", "COLOC_ABF_FALLBACK_COMPLETE")):
                if passed:
                    interpretation = "SHARED_SIGNAL_MODEL_SUPPORTED" if row["coloc_method"] == "COLOC_SUSIE" else "SINGLE_SIGNAL_FALLBACK_SHARED_MODEL_SUPPORTED"
                elif float(row["PP_H3"]) >= 0.80:
                    interpretation = "DISTINCT_SIGNAL_MODEL_SUPPORTED"
                else:
                    interpretation = "INCONCLUSIVE"
            else:
                interpretation = "NO_ESTIMABLE_COLOCALIZATION"
            record.update({
                "source_search_status": "VERIFIED_ANALYSIS",
                "primary_shared_signal_rule_pass": str(passed),
                "prior_robust": str(robust_by_signal[signal_key]),
                "colocalization_interpretation": interpretation,
            })
            canonical.append(record)

    write_tsv(args.variant_out, VARIANT_FIELDS, variants)
    write_tsv(args.credible_set_out, CREDIBLE_FIELDS, credible_sets)
    write_tsv(args.shared_variant_out, SHARED_FIELDS, shared)
    write_tsv(args.out, CANONICAL_COLOC_FIELDS, canonical)
    counts = Counter(str(row["colocalization_interpretation"]) for row in canonical)
    provenance = {
        "schema_version": "1.0.0",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifest_sha256": sha256(args.manifest),
        "lock_sha256": sha256(args.lock),
        "analysis_comparison_count": len(analysis_ids),
        "unavailable_search_count": len(unavailable_ids),
        "failed_fine_mapping_count": len(failed_ids),
        "variant_pip_row_count": len(variants),
        "credible_set_count": len(credible_sets),
        "coloc_row_count": len(canonical),
        "shared_variant_posterior_row_count": len(shared),
        "interpretation_counts": dict(sorted(counts.items())),
        "primary_p12": primary_p12,
        "p12_sensitivity_grid": p12_grid,
        "susieR_version": "0.14.2",
        "coloc_version": "5.2.3",
        "variant_output_sha256": sha256(args.variant_out),
        "credible_set_output_sha256": sha256(args.credible_set_out),
        "shared_variant_output_sha256": sha256(args.shared_variant_out),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "warning": lock["claim_limit"],
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FINEMAPPING_COLOC_RESULTS_OK analyses={len(analysis_ids)} unavailable={len(unavailable_ids)} rows={len(canonical)} interpretations={dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
