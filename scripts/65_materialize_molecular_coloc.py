#!/usr/bin/env python3
"""Materialize one locked trait-molecular comparison on the exact signed-LD universe."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import statistics
import tempfile
from pathlib import Path

import fine_mapping_contract
import molecular_contract


SUMMARY_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO", "PRIOR_WEIGHT"]
TASK_FIELDS = [
    "queue_row_id", "pair_id", "locus_id", "comparison_type", "source_search_status",
    "results_accessed_before_lock", "summary1_path", "summary1_sha256", "summary2_path",
    "summary2_sha256", "ld_path", "ld_sha256", "ld_variant_order_path",
    "ld_variant_order_sha256", "dataset1_id", "dataset1_role", "dataset1_type",
    "dataset1_N", "dataset1_case_fraction", "dataset1_sdY", "dataset1_prior_method",
    "dataset1_prior_source_path", "dataset1_prior_source_sha256", "dataset2_id",
    "dataset2_role", "dataset2_type", "dataset2_N", "dataset2_case_fraction",
    "dataset2_sdY", "dataset2_prior_method", "dataset2_prior_source_path",
    "dataset2_prior_source_sha256", "molecular_feature_id", "tissue_cell_context",
    "single_signal_fallback_justification",
]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def gzip_text(path: Path, text: str) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def aligned_sign(a1: str, a2: str, reference_a1: str, reference_a2: str) -> int | None:
    a1, a2 = a1.upper(), a2.upper()
    reference_a1, reference_a2 = reference_a1.upper(), reference_a2.upper()
    if {a1, a2} in ({"A", "T"}, {"C", "G"}):
        return None
    if (a1, a2) == (reference_a1, reference_a2):
        return 1
    if (a1, a2) == (reference_a2, reference_a1):
        return -1
    complemented = (COMPLEMENT.get(a1, ""), COMPLEMENT.get(a2, ""))
    if complemented == (reference_a1, reference_a2):
        return 1
    if complemented == (reference_a2, reference_a1):
        return -1
    return None


def quant_sd(beta: float, se: float, maf: float, n: float) -> float:
    value = 2 * maf * (1 - maf) * (n * se * se + beta * beta)
    if not math.isfinite(value) or value <= 0:
        fail("cannot estimate a positive molecular phenotype SD")
    return math.sqrt(value)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/molecular_feature_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/molecular_feature_manifest.lock.json")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--task-dir", default="results/molecular/tasks")
    parser.add_argument("--data-dir", default="data/molecular_coloc")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    if not args.materialize:
        fail("materialization requires explicit --materialize after feature-family lock review")
    root = Path(args.root).resolve()
    manifest_path, lock_path, policy_path = root / args.manifest, root / args.manifest_lock, root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    manifest, lock = molecular_contract.validate_feature_family(
        root, manifest_path, lock_path, policy_path,
        root / "results/tables/molecular_preflight.json",
    )
    selected = [row for row in manifest if row["comparison_id"] == args.comparison_id]
    if len(selected) != 1:
        fail("comparison_id must identify exactly one molecular feature comparison")
    row = selected[0]
    for path_field, hash_field in (
        ("normalized_qtl_path", "normalized_qtl_sha256"), ("trait_task_path", "trait_task_sha256"),
        ("trait_task_lock_path", "trait_task_lock_sha256"),
    ):
        path = root / row[path_field]
        if not path.is_file() or sha256(path) != row[hash_field]:
            fail(f"molecular manifest input differs: {path_field}")
    trait_task_path, trait_task_lock_path = root / row["trait_task_path"], root / row["trait_task_lock_path"]
    _, trait_tasks = read_tsv(trait_task_path)
    if len(trait_tasks) != 1:
        fail("fine-mapping trait task is not one row")
    trait_task = trait_tasks[0]
    trait_lock = json.loads(trait_task_lock_path.read_text(encoding="utf-8"))
    if (
        trait_lock.get("schema_version") != "atlas-v1.0-finemapping-task.2"
        or trait_lock.get("task_sha256") != sha256(trait_task_path)
        or trait_lock.get("script_sha256") != fine_mapping_contract.script_hashes(root)
        or trait_lock.get("results_accessed_before_lock") is not False
    ):
        fail("fine-mapping trait task differs from its pre-result lock")
    if row["trait_id"] == trait_task["dataset1_id"]:
        prefix = "dataset1"; trait_summary_path = Path(trait_task["summary1_path"]); trait_summary_hash = trait_task["summary1_sha256"]
    elif row["trait_id"] == trait_task["dataset2_id"]:
        prefix = "dataset2"; trait_summary_path = Path(trait_task["summary2_path"]); trait_summary_hash = trait_task["summary2_sha256"]
    else:
        fail("locked molecular trait is absent from the fine-mapping task")
    if not trait_summary_path.is_absolute():
        trait_summary_path = root / trait_summary_path
    ld_path, order_path = Path(trait_task["ld_path"]), Path(trait_task["ld_variant_order_path"])
    if not ld_path.is_absolute(): ld_path = root / ld_path
    if not order_path.is_absolute(): order_path = root / order_path
    for path, expected in ((trait_summary_path, trait_summary_hash), (ld_path, trait_task["ld_sha256"]), (order_path, trait_task["ld_variant_order_sha256"])):
        if not path.is_file() or sha256(path) != expected:
            fail(f"fine-mapping trait/LD input checksum mismatch: {path}")
    summary_fields, trait_rows = read_tsv(trait_summary_path)
    if not set(SUMMARY_FIELDS).issubset(summary_fields):
        fail("trait summary lacks canonical fine-mapping columns")
    trait_by_snp = {record["SNP"].lower(): record for record in trait_rows}
    if len(trait_by_snp) != len(trait_rows):
        fail("duplicate SNP in trait summary")
    qtl_fields, qtl_rows = read_tsv(root / row["normalized_qtl_path"])
    if qtl_fields != policy["normalized_qtl_schema"]:
        fail("normalized molecular QTL schema drifted")
    qtl_feature = [record for record in qtl_rows if record["feature_id"] == row["feature_id"]]
    qtl_by_snp = {record["rsid"].lower(): record for record in qtl_feature}
    if len(qtl_by_snp) != len(qtl_feature):
        fail("duplicate SNP in normalized molecular feature")
    order_fields, order_rows = read_tsv(order_path)
    if order_fields != ["SNP"] or len({record["SNP"].lower() for record in order_rows}) != len(order_rows):
        fail("invalid fine-mapping LD order")
    order = [record["SNP"].lower() for record in order_rows]
    common = [snp for snp in order if snp in trait_by_snp and snp in qtl_by_snp]
    if not policy["variant_qc"]["minimum_shared_variants"] <= len(common) <= policy["variant_qc"]["maximum_shared_variants"]:
        fail(f"comparison has {len(common)} shared variants outside locked bounds")
    output_trait: list[dict[str, object]] = []
    output_qtl: list[dict[str, object]] = []
    sd_estimates: list[float] = []
    retained: list[str] = []
    for snp in common:
        trait, qtl = trait_by_snp[snp], qtl_by_snp[snp]
        sign = aligned_sign(qtl["effect_allele"], qtl["other_allele"], trait["A1"], trait["A2"])
        if sign is None:
            continue
        beta, se, maf, n = float(qtl["beta"]), float(qtl["se"]), float(qtl["maf"]), float(qtl["n"])
        sd_estimates.append(quant_sd(beta, se, maf, n))
        retained.append(snp)
        output_trait.append({field: trait[field] for field in SUMMARY_FIELDS})
        output_qtl.append({
            "SNP": snp, "CHR": trait["CHR"], "BP": trait["BP"], "A1": trait["A1"], "A2": trait["A2"],
            "BETA": format(sign * beta, ".15g"), "SE": format(se, ".15g"), "MAF": format(maf, ".15g"),
            "INFO": "1", "PRIOR_WEIGHT": "1",
        })
    if len(retained) < policy["variant_qc"]["minimum_shared_variants"]:
        fail(f"only {len(retained)} variants remain after allele alignment")
    molecular_sd = statistics.median(sd_estimates)
    if not math.isfinite(molecular_sd) or molecular_sd <= 0:
        fail("invalid locked molecular phenotype SD")
    ld_fields, ld_rows = read_tsv(ld_path)
    if not ld_fields or ld_fields[0] != "SNP" or len(ld_rows) != len(order):
        fail("signed LD matrix has unexpected shape")
    if ld_fields[1:] != [record["SNP"] for record in order_rows]:
        fail("signed LD columns differ from variant order")
    ld_column = {name.lower(): name for name in ld_fields[1:]}
    ld_by_snp = {record["SNP"].lower(): record for record in ld_rows}
    if set(ld_by_snp) != set(order):
        fail("signed LD rows differ from variant order")
    subset_ld = []
    for snp in retained:
        source = ld_by_snp[snp]
        record: dict[str, object] = {"SNP": snp}
        for target in retained:
            value = float(source[ld_column[target]])
            if not math.isfinite(value) or abs(value) > 1.000001:
                fail("invalid signed LD value")
            record[target] = format(value, ".15g")
        subset_ld.append(record)
    task_path = root / args.task_dir / f"{args.comparison_id}.tsv"
    task_lock_path = root / args.task_dir / f"{args.comparison_id}.lock.json"
    final_data = root / args.data_dir / args.comparison_id
    if task_path.exists() or task_lock_path.exists() or final_data.exists():
        fail("molecular comparison inputs already exist; refusing overwrite")
    staging_root = root / "work/molecular"
    staging_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=args.comparison_id + ".", dir=staging_root))
    try:
        staged_data = staging / "data"; staged_data.mkdir()
        trait_out, molecular_out = staged_data / "trait.tsv.gz", staged_data / "molecular.tsv.gz"
        order_out, ld_out = staged_data / "variants.tsv", staged_data / "ld.tsv.gz"
        gzip_text(trait_out, table_text(SUMMARY_FIELDS, output_trait))
        gzip_text(molecular_out, table_text(SUMMARY_FIELDS, output_qtl))
        order_out.write_text(table_text(["SNP"], [{"SNP": snp} for snp in retained]), encoding="utf-8")
        gzip_text(ld_out, table_text(["SNP", *retained], subset_ld))
        final_trait, final_molecular = final_data / "trait.tsv.gz", final_data / "molecular.tsv.gz"
        final_order, final_ld = final_data / "variants.tsv", final_data / "ld.tsv.gz"
        molecular_dataset_id = f"{row['dataset_id']}__{row['feature_id']}"
        task_record = {
            "queue_row_id": args.comparison_id, "pair_id": row["pair_id"], "locus_id": row["locus_id"],
            "comparison_type": f"TRAIT_{row['modality'].upper()}", "source_search_status": "VERIFIED_ANALYSIS",
            "results_accessed_before_lock": "NO", "summary1_path": str(final_trait.relative_to(root)),
            "summary1_sha256": sha256(trait_out), "summary2_path": str(final_molecular.relative_to(root)),
            "summary2_sha256": sha256(molecular_out), "ld_path": str(final_ld.relative_to(root)), "ld_sha256": sha256(ld_out),
            "ld_variant_order_path": str(final_order.relative_to(root)), "ld_variant_order_sha256": sha256(order_out),
            "dataset1_id": row["trait_id"], "dataset1_role": row["trait_role"],
            "dataset1_type": trait_task[f"{prefix}_type"], "dataset1_N": trait_task[f"{prefix}_N"],
            "dataset1_case_fraction": trait_task[f"{prefix}_case_fraction"], "dataset1_sdY": trait_task[f"{prefix}_sdY"],
            "dataset1_prior_method": "FLAT", "dataset1_prior_source_path": "NA", "dataset1_prior_source_sha256": "NA",
            "dataset2_id": molecular_dataset_id, "dataset2_role": row["modality"], "dataset2_type": "quant",
            "dataset2_N": row["molecular_N"], "dataset2_case_fraction": "NA", "dataset2_sdY": format(molecular_sd, ".15g"),
            "dataset2_prior_method": "FLAT", "dataset2_prior_source_path": "NA", "dataset2_prior_source_sha256": "NA",
            "molecular_feature_id": row["feature_id"], "tissue_cell_context": row["context"],
            "single_signal_fallback_justification": "NOT_JUSTIFIED",
        }
        staged_task, staged_lock = staging / "task.tsv", staging / "task.lock.json"
        staged_task.write_text(table_text(TASK_FIELDS, [task_record]), encoding="utf-8")
        staged_lock.write_text(json.dumps({
            "schema_version": "atlas-v1.0-molecular-coloc-task.2", "comparison_id": args.comparison_id,
            "variant_count": len(retained), "molecular_sdY": molecular_sd,
            "molecular_sdY_method": "MEDIAN_COLOC_SUMMARY_STATISTIC_ESTIMATE",
            "task_sha256": sha256(staged_task), "manifest_sha256": sha256(manifest_path),
            "manifest_lock_sha256": sha256(lock_path), "policy_sha256": sha256(policy_path),
            "source_normalized_qtl_sha256": row["normalized_qtl_sha256"],
            "trait_task_sha256": row["trait_task_sha256"], "trait_task_lock_sha256": row["trait_task_lock_sha256"],
            "script_sha256": molecular_contract.script_hashes(root, "qtl"),
            "outputs": {
                "trait.tsv.gz": sha256(trait_out), "molecular.tsv.gz": sha256(molecular_out),
                "variants.tsv": sha256(order_out), "ld.tsv.gz": sha256(ld_out),
            },
            "trait_molecular_results_accessed_before_lock": False, "claim_limit": policy["claim_limit"],
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        final_data.parent.mkdir(parents=True, exist_ok=True); task_path.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staged_data, final_data); os.replace(staged_task, task_path); os.replace(staged_lock, task_lock_path)
        staging.rmdir()
        print(f"MOLECULAR_COLOC_INPUT_OK comparison={args.comparison_id} variants={len(retained)} feature={row['feature_id']}")
        return 0
    finally:
        if staging.exists(): shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
