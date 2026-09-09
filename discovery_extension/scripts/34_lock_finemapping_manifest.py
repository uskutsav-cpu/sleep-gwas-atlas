#!/usr/bin/env python3
"""Validate and checksum-lock dense locus, LD, and QTL inputs before result access."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ANALYSIS_TYPES = {
    "TRAIT_TRAIT", "SLEEP_EQTL", "EXTERNAL_EQTL", "SLEEP_SQTL",
    "EXTERNAL_SQTL", "SLEEP_PQTL", "EXTERNAL_PQTL",
}
MODALITIES = {"EQTL", "SQTL", "PQTL"}
SUMMARY_FIELDS = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO", "PRIOR_WEIGHT"}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verified_file(row: dict[str, str], path_field: str, checksum_field: str) -> Path:
    path = Path(row[path_field])
    expected = row[checksum_field]
    if not path.is_file() or len(expected) != 64 or sha256(path) != expected:
        raise SystemExit(f"ERROR: {row['queue_row_id']} missing or checksum-mismatched {path_field}")
    return path


def scalar_metadata(row: dict[str, str], prefix: str) -> None:
    dataset_type = row[f"{prefix}_type"]
    if dataset_type not in {"quant", "cc"}:
        raise SystemExit(f"ERROR: {row['queue_row_id']} invalid {prefix}_type")
    total_n = float(row[f"{prefix}_N"])
    effective_n = float(row[f"{prefix}_effective_N"])
    if not math.isfinite(total_n) or total_n <= 0 or not math.isfinite(effective_n) or effective_n <= 0:
        raise SystemExit(f"ERROR: {row['queue_row_id']} invalid {prefix} sample size")
    if dataset_type == "cc":
        case_fraction = float(row[f"{prefix}_case_fraction"])
        if not 0 < case_fraction < 1 or row[f"{prefix}_sdY"] not in {"NA", "", "."}:
            raise SystemExit(f"ERROR: {row['queue_row_id']} invalid case-control metadata for {prefix}")
    else:
        sd_y = float(row[f"{prefix}_sdY"])
        if not math.isfinite(sd_y) or sd_y <= 0 or row[f"{prefix}_case_fraction"] not in {"NA", "", "."}:
            raise SystemExit(f"ERROR: {row['queue_row_id']} invalid quantitative metadata for {prefix}")
    prior_method = row[f"{prefix}_prior_method"]
    if prior_method not in {"FLAT", "POLYFUN_FUNCTIONAL"}:
        raise SystemExit(f"ERROR: {row['queue_row_id']} invalid {prefix}_prior_method")
    if prior_method == "POLYFUN_FUNCTIONAL":
        verified_file(row, f"{prefix}_prior_source_path", f"{prefix}_prior_source_sha256")


def read_variant_order(path: Path) -> list[str]:
    fields, rows = read_tsv(path)
    if fields != ["SNP"] or not rows:
        raise SystemExit(f"ERROR: variant-order file must contain one nonempty SNP column: {path}")
    snps = [row["SNP"] for row in rows]
    if len(set(snps)) != len(snps) or any(not snp for snp in snps):
        raise SystemExit(f"ERROR: duplicate or empty SNP in variant-order file: {path}")
    return snps


def validate_summary(path: Path, row: dict[str, str], expected_snps: list[str], prefix: str, minimum_variants: int) -> list[tuple[str, int, int, str, str]]:
    fields, rows = read_tsv(path)
    if not SUMMARY_FIELDS.issubset(fields):
        raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} summary lacks fields: {sorted(SUMMARY_FIELDS - set(fields))}")
    if len(rows) != len(expected_snps) or len(rows) < minimum_variants:
        raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} summary has wrong or insufficient variant count")
    identities: list[tuple[str, int, int, str, str]] = []
    prior_weights: list[float] = []
    for expected, variant in zip(expected_snps, rows, strict=True):
        snp = variant["SNP"]
        if snp != expected:
            raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} summary differs from locked SNP order at {snp}")
        chromosome, position = int(variant["CHR"]), int(variant["BP"])
        if chromosome != int(row["CHR"]) or not int(row["start"]) <= position <= int(row["end"]):
            raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} variant outside locked locus: {snp}")
        a1, a2 = variant["A1"], variant["A2"]
        if a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"} or a1 == a2:
            raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} invalid alleles: {snp}")
        beta, se, maf, info = (float(variant[field]) for field in ("BETA", "SE", "MAF", "INFO"))
        if not all(math.isfinite(value) for value in (beta, se, maf, info)) or se <= 0 or not 0.01 <= maf <= 0.5 or info < 0.9:
            raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} violates beta/SE/MAF/INFO contract: {snp}")
        try:
            weight = float(variant["PRIOR_WEIGHT"])
        except ValueError as exc:
            raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} invalid prior weight: {snp}") from exc
        if not math.isfinite(weight) or weight < 0:
            raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} invalid prior weight: {snp}")
        prior_weights.append(weight)
        identities.append((snp, chromosome, position, a1, a2))
    if row[f"{prefix}_prior_method"] == "POLYFUN_FUNCTIONAL" and sum(prior_weights) <= 0:
        raise SystemExit(f"ERROR: {row['queue_row_id']} {prefix} PolyFun weights sum to zero")
    return identities


def validate_ld(path: Path, expected_snps: list[str], row: dict[str, str], contract: dict[str, object]) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if header != ["SNP", *expected_snps]:
            raise SystemExit(f"ERROR: {row['queue_row_id']} LD columns differ from locked SNP order")
        n = len(expected_snps)
        matrix = np.empty((n, n), dtype=float)
        row_count = 0
        for index, values in enumerate(reader):
            if index >= n or len(values) != n + 1 or values[0] != expected_snps[index]:
                raise SystemExit(f"ERROR: {row['queue_row_id']} LD rows differ from locked SNP order")
            try:
                matrix[index, :] = np.asarray(values[1:], dtype=float)
            except ValueError as exc:
                raise SystemExit(f"ERROR: {row['queue_row_id']} LD contains nonnumeric values") from exc
            row_count += 1
    if row_count != len(expected_snps) or not np.isfinite(matrix).all():
        raise SystemExit(f"ERROR: {row['queue_row_id']} LD matrix has wrong dimensions or non-finite entries")
    symmetry = float(np.max(np.abs(matrix - matrix.T)))
    diagonal = float(np.max(np.abs(np.diag(matrix) - 1)))
    if symmetry > float(contract["symmetry_tolerance"]) or diagonal > float(contract["diagonal_tolerance"]):
        raise SystemExit(f"ERROR: {row['queue_row_id']} LD matrix fails symmetry/unit-diagonal contract")
    if float(np.max(np.abs(matrix))) > 1 + float(contract["diagonal_tolerance"]):
        raise SystemExit(f"ERROR: {row['queue_row_id']} LD correlation lies outside [-1,1]")
    minimum_eigenvalue = float(np.linalg.eigvalsh(matrix)[0])
    if minimum_eigenvalue < float(contract["minimum_eigenvalue_tolerance"]):
        raise SystemExit(f"ERROR: {row['queue_row_id']} LD matrix is not positive semidefinite")
    return {"symmetry_max_abs": symmetry, "diagonal_max_abs": diagonal, "minimum_eigenvalue": minimum_eigenvalue}


def unavailable_row(row: dict[str, str]) -> None:
    modality = row["qtl_modality"]
    if modality not in MODALITIES or row["comparison_type"] != f"{modality}_SEARCH":
        raise SystemExit(f"ERROR: {row['queue_row_id']} unavailable row must be an eQTL/sQTL/pQTL search")
    verified_file(row, "source_search_log_path", "source_search_log_sha256")
    for field in ("source_id", "source_url", "source_access_date", "curator", "curation_date", "notes"):
        if row[field] in {"", "PENDING", "NA", "."}:
            raise SystemExit(f"ERROR: {row['queue_row_id']} unavailable search lacks {field}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_input_queue.tsv"))
    parser.add_argument("--candidate-lock", type=Path, default=Path("discovery_extension/config/fine_mapping_candidate_family.lock.json"))
    parser.add_argument("--readiness", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_readiness.tsv"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/fine_mapping_colocalization_contract.json"))
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/fine_mapping_sources.tsv"))
    parser.add_argument("--references", type=Path, default=Path("discovery_extension/config/fine_mapping_method_references.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/config/fine_mapping_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/fine_mapping_manifest.lock.json"))
    args = parser.parse_args()

    fields, rows = read_tsv(args.queue)
    family = json.loads(args.candidate_lock.read_text(encoding="utf-8"))
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    _, readiness = read_tsv(args.readiness)
    if not rows:
        raise SystemExit("ERROR: fine-mapping queue is empty")
    if family.get("results_accessed_before_lock") is not False:
        raise SystemExit("ERROR: fine-mapping candidate locus family was not result-free")
    if len(readiness) != 2 or any(row["code_status"] != "PASS" for row in readiness):
        raise SystemExit("ERROR: pinned SuSiE/coloc code is not ready")
    if len({row["queue_row_id"] for row in rows}) != len(rows):
        raise SystemExit("ERROR: duplicate fine-mapping queue row identifier")
    locus_order = family.get("locus_keys_in_locked_order")
    loci_in_queue: list[str] = []
    for row in rows:
        key = f"{row['pair_id']}__{row['locus_id']}"
        if key not in loci_in_queue:
            loci_in_queue.append(key)
    if loci_in_queue != locus_order:
        raise SystemExit("ERROR: curated fine-mapping queue differs from the locked locus family/order")

    minimum_variants = int(contract["input_contract"]["minimum_locus_variants"])
    diagnostic_by_comparison: dict[str, dict[str, float | int]] = {}
    analysis_ids: list[str] = []
    unavailable_ids: list[str] = []
    variant_counts: dict[str, int] = {}
    for locus_key in locus_order:
        locus_rows = [row for row in rows if f"{row['pair_id']}__{row['locus_id']}" == locus_key]
        trait_trait = [row for row in locus_rows if row["comparison_type"] == "TRAIT_TRAIT"]
        if len(trait_trait) != 1 or trait_trait[0]["source_search_status"] != "VERIFIED_ANALYSIS":
            raise SystemExit(f"ERROR: {locus_key} requires exactly one verified trait-trait comparison")
        for modality in sorted(MODALITIES):
            modality_rows = [row for row in locus_rows if row["qtl_modality"] == modality]
            if not modality_rows or any(row["source_search_status"] not in {"VERIFIED_ANALYSIS", "NO_SUITABLE_DATASET"} for row in modality_rows):
                raise SystemExit(f"ERROR: {locus_key} lacks completed {modality} source search")
            statuses = {row["source_search_status"] for row in modality_rows}
            if statuses == {"VERIFIED_ANALYSIS", "NO_SUITABLE_DATASET"} or (statuses == {"NO_SUITABLE_DATASET"} and len(modality_rows) != 1):
                raise SystemExit(f"ERROR: {locus_key} has contradictory or duplicated {modality} source-search disposition")
        for row in locus_rows:
            identity = row["queue_row_id"]
            if row["results_accessed_before_lock"] != "NO":
                raise SystemExit(f"ERROR: result access preceded fine-mapping lock: {identity}")
            if row["source_search_status"] == "NO_SUITABLE_DATASET":
                unavailable_row(row)
                unavailable_ids.append(identity)
                continue
            if row["source_search_status"] != "VERIFIED_ANALYSIS" or row["comparison_type"] not in ANALYSIS_TYPES:
                raise SystemExit(f"ERROR: unresolved or invalid analysis row: {identity}")
            if row["single_signal_fallback_justification"] in {"", "PENDING", "NA", "."}:
                raise SystemExit(f"ERROR: {identity} lacks a pre-result coloc.abf fallback decision")
            if row["build"] != "GRCh37" or row["ancestry"] != "EUR":
                raise SystemExit(f"ERROR: {identity} build/ancestry differs from primary contract")
            if row["schema_status"] != "VERIFIED" or row["effect_allele_alignment_status"] != "VERIFIED" or row["dense_variant_status"] != "FULL_RESOLUTION":
                raise SystemExit(f"ERROR: {identity} schema/allele/dense-input contract failed")
            if row["phenotype_compatibility"] != "VERIFIED" or row["sample_overlap_status"] not in {"NON_OVERLAPPING", "OVERLAP_DISCLOSED"}:
                raise SystemExit(f"ERROR: {identity} phenotype/sample-overlap status is unresolved")
            if row["ld_source_type"] not in set(contract["ld_contract"]["source_types"]) or row["ld_ancestry"] != "EUR" or row["ld_build"] != "GRCh37" or int(row["ld_sample_size"]) <= 0:
                raise SystemExit(f"ERROR: {identity} LD source metadata failed")
            for field in ("source_id", "source_url", "source_access_date", "source_license", "curator", "curation_date"):
                if row[field] in {"", "PENDING", "NA", "."}:
                    raise SystemExit(f"ERROR: {identity} lacks {field}")
            if row["qtl_modality"] in MODALITIES:
                verified_file(row, "source_search_log_path", "source_search_log_sha256")
                for field in ("molecular_feature_id", "molecular_feature_name", "tissue_cell_context"):
                    if row[field] in {"", "PENDING", "NA", "."}:
                        raise SystemExit(f"ERROR: {identity} molecular-QTL comparison lacks {field}")
            scalar_metadata(row, "dataset1")
            scalar_metadata(row, "dataset2")
            summary1 = verified_file(row, "summary1_path", "summary1_sha256")
            summary2 = verified_file(row, "summary2_path", "summary2_sha256")
            ld_path = verified_file(row, "ld_path", "ld_sha256")
            order_path = verified_file(row, "ld_variant_order_path", "ld_variant_order_sha256")
            snps = read_variant_order(order_path)
            ids1 = validate_summary(summary1, row, snps, "dataset1", minimum_variants)
            ids2 = validate_summary(summary2, row, snps, "dataset2", minimum_variants)
            if ids1 != ids2:
                raise SystemExit(f"ERROR: {identity} summaries are not allele/coordinate identical")
            diagnostics = validate_ld(ld_path, snps, row, contract["ld_contract"])
            diagnostics["variant_count"] = len(snps)
            diagnostic_by_comparison[identity] = diagnostics
            analysis_ids.append(identity)
            variant_counts[identity] = len(snps)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lock = {
        "schema_version": "1.0.0",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "results_accessed_before_lock": False,
        "manifest_sha256": sha256(args.out),
        "locus_keys_in_locked_order": locus_order,
        "analysis_comparison_ids_in_locked_order": analysis_ids,
        "unavailable_search_ids_in_locked_order": unavailable_ids,
        "variant_counts": variant_counts,
        "ld_diagnostics": diagnostic_by_comparison,
        "susie_parameters": {
            "L": contract["fine_mapping"]["maximum_causal_signals"],
            "coverage": contract["fine_mapping"]["credible_set_coverage"],
            "min_abs_corr": contract["fine_mapping"]["minimum_absolute_correlation_for_credible_set"],
            "maxit": contract["fine_mapping"]["maximum_iterations"],
            "estimate_residual_variance": contract["fine_mapping"]["estimate_residual_variance"],
        },
        "coloc_primary_priors": contract["colocalization"]["primary_priors"],
        "coloc_p12_sensitivity_grid": contract["colocalization"]["p12_sensitivity_grid"],
        "primary_shared_signal_rule": contract["colocalization"]["primary_shared_signal_rule"],
        "candidate_family_lock_sha256": sha256(args.candidate_lock),
        "queue_sha256": sha256(args.queue),
        "readiness_sha256": sha256(args.readiness),
        "contract_sha256": sha256(args.contract),
        "sources_sha256": sha256(args.sources),
        "references_sha256": sha256(args.references),
        "claim_limit": contract["claim_limit"],
    }
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    args.lock.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FINEMAPPING_MANIFEST_LOCKED loci={len(locus_order)} analyses={len(analysis_ids)} unavailable={len(unavailable_ids)} result_free=true")


if __name__ == "__main__":
    main()
