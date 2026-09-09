#!/usr/bin/env python3
"""Run one locked FUMA scRNA MAGMA cell-type task from shared trait gene results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import subprocess
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, str]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def fmt(value: float) -> str:
    return format(value, ".12g")


def parse_gsa(path: Path, expected_variables: list[str]) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            values = line.split()
            if not values:
                continue
            if values[0] == "VARIABLE":
                header = values
                continue
            if header is None or values[0].startswith("#"):
                continue
            if len(values) != len(header):
                fail(f"MAGMA gene-property output contains a malformed row: {path}")
            row = dict(zip(header, values))
            if row.get("TYPE") != "COVAR":
                fail(f"MAGMA gene-property output contains a non-covariate row: {path}")
            rows.append(row)
    if header is None or not {"VARIABLE", "TYPE", "NGENES", "BETA", "SE", "P"}.issubset(header):
        fail(f"MAGMA gene-property output lacks its required result columns: {path}")
    if [row["VARIABLE"] for row in rows] != expected_variables:
        fail(f"MAGMA gene-property variables differ from the frozen matrix header: {path}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--provenance-out", required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["fuma_scrna"]
    task_manifest_path = root / args.manifest
    _, tasks = read_tsv(task_manifest_path)
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    if (
        task["analysis_family"] != "cell_type" or task["method"] != "MAGMA_Celltyping"
        or task["source_id"] != spec["source_id"]
        or task["domain"] not in policy["cell_types"]["required_domains"]
    ):
        fail("task is not supported by the FUMA scRNA adapter")
    source_manifest_path = root / spec["component_manifest"]
    if (
        not source_manifest_path.is_file()
        or sha256(source_manifest_path) != spec["component_manifest_sha256"]
    ):
        fail("FUMA source manifest differs from its policy pin")
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    if task["source_release"] != source_manifest["release"]:
        fail("FUMA task release differs from the source manifest")
    components = {row["component_id"]: row for row in source_manifest["components"]}
    binary = root / spec["magma_binary_path"]
    binary_pin = components["MAGMA_V1_10_ARM64_BINARY"]
    if (
        not binary.is_file() or binary.stat().st_size != binary_pin["bytes"]
        or sha256(binary) != binary_pin["sha256"]
    ):
        fail("MAGMA binary differs from its exact runtime pin")
    matrices_provenance_path = root / spec["matrix_cache_provenance_path"]
    if not matrices_provenance_path.is_file():
        fail("materialized FUMA matrix provenance is absent")
    matrices_provenance = json.loads(matrices_provenance_path.read_text(encoding="utf-8"))
    if (
        matrices_provenance.get("policy_sha256") != sha256(policy_path)
        or matrices_provenance.get("component_manifest_sha256") != sha256(source_manifest_path)
        or matrices_provenance.get("mode") != "MATRICES"
    ):
        fail("materialized FUMA matrices differ from policy-bound provenance")
    gene_prefix = root / spec["gene_results_dir"] / task["trait_id"] / task["trait_id"]
    genes_raw = Path(str(gene_prefix) + ".genes.raw")
    gene_provenance_path = gene_prefix.parent / "provenance.json"
    if not genes_raw.is_file() or not gene_provenance_path.is_file():
        fail("trait MAGMA gene results and provenance are incomplete")
    gene_provenance = json.loads(gene_provenance_path.read_text(encoding="utf-8"))
    if (
        gene_provenance.get("trait_id") != task["trait_id"]
        or gene_provenance.get("policy_sha256") != sha256(policy_path)
        or gene_provenance.get("genes_raw_sha256") != sha256(genes_raw)
    ):
        fail("trait MAGMA gene results differ from policy-bound provenance")
    datasets = [row for row in source_manifest["datasets"] if row["domain"] == task["domain"]]
    if len(datasets) != 2:
        fail("FUMA task domain does not have exactly two prespecified matrices")
    canonical_fields = policy["cell_types"]["canonical_fields"]
    normalized: list[dict[str, str]] = []
    run_records: list[dict[str, object]] = []
    run_dir = root / "work/interpretation_automatic/fuma_scrna/gene_property" / task["trait_id"] / task["domain"]
    run_dir.mkdir(parents=True, exist_ok=True)
    for dataset in datasets:
        dataset_id = dataset["dataset_id"]
        matrix = root / spec["matrix_cache_dir"] / f"{dataset_id}.txt"
        if (
            not matrix.is_file() or matrix.stat().st_size != dataset["text_bytes"]
            or sha256(matrix) != dataset["text_sha256"]
        ):
            fail(f"materialized FUMA matrix differs from its pin: {dataset_id}")
        with matrix.open(encoding="utf-8") as handle:
            matrix_header = handle.readline().split()
        variables = matrix_header[1:-1]
        if len(variables) != dataset["cell_type_count"] or matrix_header[-1] != "Average":
            fail(f"materialized FUMA matrix header differs from its manifest: {dataset_id}")
        prefix = run_dir / dataset_id
        gsa, log = Path(str(prefix) + ".gsa.out"), Path(str(prefix) + ".log")
        if gsa.exists() or log.exists():
            fail(f"FUMA gene-property output already exists; refusing overwrite: {dataset_id}")
        temporary_prefix = prefix.with_name(prefix.name + ".tmp")
        temporary_gsa, temporary_log = (
            Path(str(temporary_prefix) + ".gsa.out"), Path(str(temporary_prefix) + ".log")
        )
        if temporary_gsa.exists() or temporary_log.exists():
            fail(f"stale FUMA task temporary output exists; inspect it before retrying: {dataset_id}")
        command = [
            str(binary), "--gene-results", str(genes_raw), "--gene-covar", str(matrix),
            "--model", "condition-hide=Average", "direction=greater",
            "--settings", "abbreviate=0", "--out", str(temporary_prefix),
        ]
        result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
        if result.returncode or not temporary_gsa.is_file() or not temporary_log.is_file():
            fail(f"MAGMA FUMA task failed for {dataset_id}: {(result.stdout + result.stderr).strip()}")
        result_rows = parse_gsa(temporary_gsa, variables)
        for result_row in result_rows:
            try:
                beta, se, p_value = (
                    float(result_row["BETA"]), float(result_row["SE"]), float(result_row["P"])
                )
                n_genes = int(result_row["NGENES"])
            except ValueError as exc:
                fail(f"MAGMA FUMA result is nonnumeric for {dataset_id}/{result_row['VARIABLE']}")
                raise AssertionError from exc
            if (
                not all(math.isfinite(value) for value in (beta, se, p_value))
                or se < 0 or not 0 <= p_value <= 1
                or n_genes < int(spec["minimum_genes_in_gene_property_model"])
            ):
                fail(f"MAGMA FUMA result is outside policy for {dataset_id}/{result_row['VARIABLE']}")
            variable = result_row["VARIABLE"]
            normalized.append({
                "cell_type_id": f"{dataset_id}::{variable}", "cell_type": variable,
                "tissue": dataset["tissue"], "domain": task["domain"],
                "method": task["method"], "trait_or_locus_id": task["trait_id"],
                "effect": fmt(beta), "se": fmt(se), "p_value": fmt(p_value), "fdr": "NA",
                "evidence_level": "UNADJUSTED_MAGMA_GENE_PROPERTY",
                "dataset": task["source_id"], "version": task["source_release"],
                "provenance_id": "FUMA_SCRNA_TASK_PENDING",
            })
        temporary_gsa.replace(gsa)
        temporary_log.replace(log)
        run_records.append({
            "dataset_id": dataset_id, "domain": dataset["domain"], "species": dataset["species"],
            "tissue": dataset["tissue"], "matrix_sha256": sha256(matrix), "command": command,
            "result_rows": len(result_rows), "gsa_path": str(gsa), "gsa_sha256": sha256(gsa),
            "log_path": str(log), "log_sha256": sha256(log),
        })
    normalized.sort(key=lambda row: (row["cell_type_id"], row["p_value"]))
    if any(list(row) != canonical_fields for row in normalized):
        fail("FUMA normalized row differs from the locked cell-type schema")
    payload = table_text(canonical_fields, normalized)
    out_path = root / args.out
    atomic_text(out_path, payload)
    provenance_path = root / args.provenance_out
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"], "terminal_status": "COMPLETED",
        "terminal_reason": (
            f"MAGMA v1.10 gene-property models completed for both prespecified {task['domain']} "
            f"FUMA scRNA matrices and produced {len(normalized)} unadjusted cell-type rows"
        ),
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(task_manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"], "source_id": task["source_id"],
        "source_release": task["source_release"], "source_manifest_sha256": sha256(source_manifest_path),
        "matrices_provenance_sha256": sha256(matrices_provenance_path),
        "gene_results_provenance_sha256": sha256(gene_provenance_path),
        "genes_raw_sha256": sha256(genes_raw), "runs": run_records,
        "result_rows": len(normalized), "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "multiple_testing_family": spec["multiple_testing_family"],
        "claim_limit": spec["claim_limit"],
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"FUMA_SCRNA_TASK_OK task={task['task_id']} rows={len(normalized)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
