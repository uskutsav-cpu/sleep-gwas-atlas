#!/usr/bin/env python3
"""Materialize allele-aligned dense summaries and signed LAVA LD for one locus."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import subprocess
from pathlib import Path

import fine_mapping_contract


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
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str, field: str, identity: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"invalid {field} for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed):
        fail(f"non-finite {field} for {identity}")
    return parsed


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


def reference_rows(path: Path, chromosome: int, start: int, end: int, minimum_maf: float) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"LAVA reference info is absent: {path}")
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        required = {"SNP", "CHR", "POS", "A1", "A2", "FREQ", "NOBS"}
        if not required.issubset(header):
            fail(f"unexpected LAVA reference header in {path}")
        output = []
        seen = set()
        for line_number, line in enumerate(handle, start=2):
            values = line.split()
            if len(values) != len(header):
                fail(f"malformed LAVA reference info line {line_number}")
            row = dict(zip(header, values))
            if int(row["CHR"]) != chromosome:
                fail(f"wrong chromosome in LAVA info line {line_number}")
            position = int(row["POS"])
            if not start <= position <= end:
                continue
            snp = row["SNP"].lower()
            if snp in seen:
                fail(f"duplicate LAVA reference SNP: {snp}")
            seen.add(snp)
            frequency = number(row["FREQ"], "reference frequency", snp)
            maf = min(frequency, 1 - frequency)
            if maf <= minimum_maf:
                continue
            if row["A1"].upper() not in COMPLEMENT or row["A2"].upper() not in COMPLEMENT:
                continue
            row["SNP"] = snp
            row["MAF"] = format(maf, ".12g")
            output.append(row)
    return output


def extract_summary(path: Path, reference: dict[str, dict[str, str]], chromosome: int) -> dict[str, dict[str, str]]:
    if not path.is_file():
        fail(f"full summary statistics are absent: {path}")
    output: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
        if reader.fieldnames != expected:
            fail(f"unexpected canonical summary header: {path}")
        for row in reader:
            snp = row["SNP"].lower()
            ref = reference.get(snp)
            if ref is None:
                continue
            if int(row["CHR"]) != chromosome or int(row["BP"]) != int(ref["POS"]):
                fail(f"build/coordinate conflict for {snp} in {path}")
            sign = aligned_sign(row["A1"], row["A2"], ref["A1"], ref["A2"])
            if sign is None:
                continue
            if snp in output:
                fail(f"duplicate locus SNP {snp} in {path}")
            beta = number(row["BETA"], "BETA", snp) * sign
            se = number(row["SE"], "SE", snp)
            if se <= 0:
                fail(f"non-positive SE for {snp} in {path}")
            output[snp] = {
                "SNP": snp, "CHR": str(chromosome), "BP": ref["POS"],
                "A1": ref["A1"].upper(), "A2": ref["A2"].upper(),
                "BETA": format(beta, ".17g"), "SE": format(se, ".17g"),
                "MAF": ref["MAF"], "INFO": "NA", "PRIOR_WEIGHT": "1",
            }
    return output


def gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    temporary.replace(path)


def plain_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def dataset_metadata(panel: dict[str, str], rows: list[dict[str, str]]) -> tuple[str, str, str, str]:
    if panel["type"] == "binary":
        cases, controls = int(panel["ncase"]), int(panel["ncontrol"])
        total = cases + controls
        return "cc", str(total), format(cases / total, ".17g"), "NA"
    sample_size = int(panel["n_total"])
    x_values = []
    y_values = []
    for row in rows:
        variance = float(row["SE"]) ** 2
        maf = float(row["MAF"])
        x_values.append(1 / variance)
        y_values.append(2 * sample_size * maf * (1 - maf))
    denominator = sum(value * value for value in x_values)
    coefficient = sum(x * y for x, y in zip(x_values, y_values)) / denominator if denominator else -1
    if not math.isfinite(coefficient) or coefficient <= 0:
        fail(f"could not estimate quantitative phenotype SD for {panel['trait_id']}")
    return "quant", str(sample_size), "NA", format(math.sqrt(coefficient), ".17g")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shared_locus_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/fine_mapping_locus_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/fine_mapping_locus_manifest.lock.json")
    parser.add_argument("--policy", default="config/fine_mapping_analysis_policy.json")
    parser.add_argument("--preflight", default="results/tables/fine_mapping_preflight.json")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    if not SAFE_ID.fullmatch(args.shared_locus_id):
        fail("unsafe shared-locus identifier")
    root = Path(args.root).resolve()
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    policy_path, policy = fine_mapping_contract.load_policy(root, args.policy)
    preflight_path = root / args.preflight
    preflight = fine_mapping_contract.validate_preflight(root, policy_path, policy, preflight_path)
    manifest, lock = fine_mapping_contract.validate_manifest(
        root, manifest_path, manifest_lock_path, policy_path, policy, preflight_path,
    )
    selected = [row for row in manifest if row["shared_locus_id"] == args.shared_locus_id]
    if len(selected) != 1:
        fail("shared_locus_id does not identify exactly one locked locus")
    row = selected[0]
    if row["comparison_id"] not in lock["comparison_ids_in_locked_order"]:
        fail("comparison is outside the locked fine-mapping family")
    outputs = [root / row[field] for field in ("summary1_path", "summary2_path", "variant_order_path", "ld_path", "task_path")]
    task_lock = outputs[-1].with_suffix(".lock.json")
    outputs.append(task_lock)
    if not args.materialize:
        print("No materialization requested. Re-run with --materialize after reviewing the locked locus and disk plan.")
        return 0
    existing = [str(path.relative_to(root)) for path in outputs if path.exists()]
    if existing:
        fail("refusing to overwrite existing fine-mapping input/task: " + ", ".join(existing))
    if outputs[0].parent.exists():
        fail(f"refusing non-atomic install into existing locus directory: {outputs[0].parent.relative_to(root)}")
    staging = root / "work/fine_mapping" / f"{args.shared_locus_id}.building"
    if staging.exists():
        fail(f"fine-mapping staging directory already exists: {staging.relative_to(root)}")
    staging_data = staging / "data"
    staging_data.mkdir(parents=True)

    chromosome, start, end = int(row["chromosome"]), int(row["start_bp"]), int(row["end_bp"])
    info_path = root / f"{row['reference_prefix']}_chr{chromosome}.info"
    bcor_path = root / f"{row['reference_prefix']}_chr{chromosome}.bcor"
    reference_provenance_path, reference_provenance = fine_mapping_contract.validate_lava_reference(
        root, policy, rehash_payloads=False,
    )
    reference_records = {
        entry["path"]: entry for entry in reference_provenance["extracted_files"]
    }
    info_relative = str(info_path.relative_to(root))
    bcor_relative = str(bcor_path.relative_to(root))
    info_record = reference_records.get(info_relative)
    bcor_record = reference_records.get(bcor_relative)
    if (
        info_record is None or bcor_record is None or not info_path.is_file() or not bcor_path.is_file()
        or info_path.stat().st_size != int(info_record["bytes"])
        or bcor_path.stat().st_size != int(bcor_record["bytes"])
        or sha256(info_path) != info_record["sha256"]
    ):
        fail("chromosome-specific LAVA reference differs from its sealed provenance")
    reference_list = reference_rows(info_path, chromosome, start, end, float(policy["minimum_maf"]))
    reference = {entry["SNP"]: entry for entry in reference_list}
    sleep = extract_summary(root / row["sleep_full_input"], reference, chromosome)
    non_sleep = extract_summary(root / row["non_sleep_full_input"], reference, chromosome)
    common = [entry["SNP"] for entry in reference_list if entry["SNP"] in sleep and entry["SNP"] in non_sleep]
    if not policy["minimum_locus_variants"] <= len(common) <= policy["maximum_locus_variants"]:
        fail(f"aligned locus has {len(common)} variants outside the locked bounds")
    sleep_rows = [sleep[snp] for snp in common]
    non_sleep_rows = [non_sleep[snp] for snp in common]
    summary1, summary2, order_path, ld_path, task_path = outputs[:5]
    staged_summary1 = staging_data / summary1.name
    staged_summary2 = staging_data / summary2.name
    staged_order = staging_data / order_path.name
    staged_ld = staging_data / ld_path.name
    staged_task = staging / task_path.name
    staged_task_lock = staging / task_lock.name
    gzip_tsv(staged_summary1, SUMMARY_FIELDS, sleep_rows)
    gzip_tsv(staged_summary2, SUMMARY_FIELDS, non_sleep_rows)
    plain_tsv(staged_order, ["SNP"], [{"SNP": snp} for snp in common])
    helper = root / "scripts/58_extract_lava_ld.R"
    result = subprocess.run(
        [str(root / ".r-env/bin/Rscript"), str(helper), str(root / row["reference_prefix"]),
         str(chromosome), str(staged_order), str(staged_ld), str(policy["reference"]["minimum_sample_size"])],
        cwd=root, capture_output=True, text=True, check=False,
    )
    if result.returncode or not staged_ld.is_file():
        fail(f"signed-LD extraction failed: {(result.stdout + result.stderr).strip()}")

    panel = {entry["trait_id"]: entry for entry in read_tsv(root / "config/analysis_panel.tsv")}
    first = dataset_metadata(panel[row["sleep_trait"]], sleep_rows)
    second = dataset_metadata(panel[row["non_sleep_trait"]], non_sleep_rows)
    task = {
        "queue_row_id": row["comparison_id"], "pair_id": row["pair_id"],
        "locus_id": row["shared_locus_id"], "comparison_type": "TRAIT_TRAIT",
        "source_search_status": "VERIFIED_ANALYSIS", "results_accessed_before_lock": "NO",
        "summary1_path": str(summary1), "summary1_sha256": sha256(staged_summary1),
        "summary2_path": str(summary2), "summary2_sha256": sha256(staged_summary2),
        "ld_path": str(ld_path), "ld_sha256": sha256(staged_ld),
        "ld_variant_order_path": str(order_path), "ld_variant_order_sha256": sha256(staged_order),
        "dataset1_id": row["sleep_trait"], "dataset1_role": "SLEEP_TRAIT",
        "dataset1_type": first[0], "dataset1_N": first[1], "dataset1_case_fraction": first[2],
        "dataset1_sdY": first[3], "dataset1_prior_method": "FLAT",
        "dataset1_prior_source_path": "NA", "dataset1_prior_source_sha256": "NA",
        "dataset2_id": row["non_sleep_trait"], "dataset2_role": "NON_SLEEP_TRAIT",
        "dataset2_type": second[0], "dataset2_N": second[1], "dataset2_case_fraction": second[2],
        "dataset2_sdY": second[3], "dataset2_prior_method": "FLAT",
        "dataset2_prior_source_path": "NA", "dataset2_prior_source_sha256": "NA",
        "molecular_feature_id": "NA", "tissue_cell_context": "NA",
        "single_signal_fallback_justification": "NOT_JUSTIFIED",
    }
    plain_tsv(staged_task, TASK_FIELDS, [task])
    staged_hashes = {
        str(summary1.relative_to(root)): sha256(staged_summary1),
        str(summary2.relative_to(root)): sha256(staged_summary2),
        str(order_path.relative_to(root)): sha256(staged_order),
        str(ld_path.relative_to(root)): sha256(staged_ld),
    }
    source_inputs = {}
    for role, field in (("sleep", "sleep_full_input"), ("non_sleep", "non_sleep_full_input")):
        source_path = root / fine_mapping_contract.safe_relative(row[field], field)
        source_inputs[role] = {
            "path": str(source_path.relative_to(root)), "bytes": source_path.stat().st_size,
            "sha256": sha256(source_path),
        }
    staged_task_lock.write_text(json.dumps({
        "schema_version": "atlas-v1.0-finemapping-task.2",
        "comparison_id": row["comparison_id"], "shared_locus_id": row["shared_locus_id"],
        "variant_count": len(common), "task_sha256": sha256(staged_task),
        "manifest_sha256": sha256(manifest_path), "manifest_lock_sha256": sha256(manifest_lock_path),
        "policy_sha256": sha256(policy_path), "preflight_sha256": sha256(preflight_path),
        "reference_provenance_sha256": sha256(reference_provenance_path),
        "reference_content_verified_by_preflight": preflight["lava_reference_content_verified"],
        "reference_info_sha256": info_record["sha256"],
        "reference_bcor_path": str(bcor_path.relative_to(root)),
        "reference_bcor_bytes": bcor_record["bytes"],
        "reference_bcor_sha256": bcor_record["sha256"],
        "source_inputs": source_inputs,
        "outputs": staged_hashes,
        "script_sha256": fine_mapping_contract.script_hashes(root),
        "results_accessed_before_lock": False, "claim_limit": policy["claim_limit"],
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary1.parent.parent.mkdir(parents=True, exist_ok=True)
    task_path.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staging_data, summary1.parent)
    os.replace(staged_task, task_path)
    os.replace(staged_task_lock, task_lock)
    staging.rmdir()
    print(
        f"FINEMAPPING_LOCUS_INPUT_OK locus={row['shared_locus_id']} variants={len(common)} "
        f"sleep={row['sleep_trait']} non_sleep={row['non_sleep_trait']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
