#!/usr/bin/env python3
"""Run one locked CATlas adult scATAC cell-type enrichment task."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

from pathway_sources import hypergeometric_right_tail
from catlas_policy import reference_policy_sha256, trait_policy_sha256


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
    handle = (
        gzip.open(path, "rt", encoding="utf-8", newline="")
        if path.name.endswith(".gz")
        else path.open(encoding="utf-8", newline="")
    )
    with handle:
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
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def fmt(value: float) -> str:
    return format(value, ".12g")


def compute_enrichment_rows(
    policy: dict[str, object], source_manifest: dict[str, object], task: dict[str, str],
    universe_rows: list[dict[str, str]], trait_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    spec = policy["catlas_adult_v4"]
    p_values: dict[str, float] = {}
    for row in trait_rows:
        try:
            p_value = float(row["P"])
        except (KeyError, ValueError):
            fail("CATlas trait cache contains an invalid P value")
        if not math.isfinite(p_value) or not 0 < p_value <= 1 or row["SNP"] in p_values:
            fail("CATlas trait cache contains an invalid or duplicate variant")
        p_values[row["SNP"]] = p_value
    annotations: dict[str, tuple[bool, set[str]]] = {}
    for row in universe_rows:
        snp = row["SNP"]
        if snp in annotations or row["IN_ADULT_CCRE"] not in {"0", "1"}:
            fail("CATlas fixed-universe cache contains an invalid or duplicate variant")
        cells = set() if row["CELL_TYPE_IDS"] == "NA" else set(row["CELL_TYPE_IDS"].split(";"))
        annotations[snp] = (row["IN_ADULT_CCRE"] == "1", cells)
    if not set(p_values).issubset(annotations):
        fail("CATlas trait cache contains a variant outside the fixed universe")
    background = {
        snp for snp in p_values
        if annotations[snp][0]
    }
    if not background:
        fail("trait has no CATlas adult-cCRE background variants")
    signal = {
        snp for snp in background
        if p_values[snp] <= float(spec["instrument_p_threshold"])
    }
    selected_cells = [
        row for row in source_manifest["selected_cells"]
        if row["domain"] == task["domain"]
    ]
    expected_count = int(spec["selected_cells_by_domain"][task["domain"]])
    if len(selected_cells) != expected_count:
        fail("CATlas task domain differs from the frozen cell selection")
    rows: list[dict[str, str]] = []
    for cell in selected_cells:
        cell_id = f"{spec['cell_type_id_prefix']}::{cell['metadata_cell_type']}"
        cell_background = {snp for snp in background if cell_id in annotations[snp][1]}
        overlap = len(signal & cell_background)
        expected = len(signal) * len(cell_background) / len(background)
        effect = overlap / expected if expected > 0 else 0.0
        p_value = hypergeometric_right_tail(
            overlap, len(background), len(cell_background), len(signal),
        )
        rows.append({
            "cell_type_id": cell_id,
            "cell_type": cell["cell_type"],
            "tissue": cell["tissue"],
            "domain": task["domain"],
            "method": task["method"],
            "trait_or_locus_id": task["trait_id"],
            "effect": fmt(effect),
            "se": "NA",
            "p_value": fmt(p_value),
            "fdr": "NA",
            "evidence_level": "UNADJUSTED_CATLAS_HYPERGEOMETRIC",
            "dataset": task["source_id"],
            "version": task["source_release"],
            "provenance_id": "CATLAS_TASK_PENDING",
        })
    return rows, {
        "trait_fixed_universe_variants": len(p_values),
        "adult_ccre_background_variants": len(background),
        "genome_wide_significant_background_variants": len(signal),
        "domain_cell_types": len(selected_cells),
    }


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
    spec = policy["catlas_adult_v4"]
    task_manifest_path = root / args.manifest
    _, tasks = read_tsv(task_manifest_path)
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    if (
        task["analysis_family"] != "cell_type" or task["method"] != "scATAC_enrichment"
        or task["source_id"] != spec["source_id"]
        or task["domain"] not in policy["cell_types"]["required_domains"]
    ):
        fail("task is not supported by the CATlas scATAC adapter")

    source_manifest_path = root / spec["component_manifest"]
    if not source_manifest_path.is_file() or sha256(source_manifest_path) != spec["component_manifest_sha256"]:
        fail("CATlas source manifest differs from its policy pin")
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    if task["source_release"] != source_manifest["source_release"]:
        fail("CATlas task release differs from the source manifest")

    universe_path = root / spec["variant_cache_path"]
    universe_provenance_path = root / spec["variant_cache_provenance_path"]
    if not universe_path.is_file() or not universe_provenance_path.is_file():
        fail("CATlas fixed-universe cache and provenance are incomplete")
    universe_provenance = json.loads(universe_provenance_path.read_text(encoding="utf-8"))
    if (
        universe_provenance.get("reference_policy_sha256") != reference_policy_sha256(policy)
        or universe_provenance.get("component_manifest_sha256") != sha256(source_manifest_path)
        or universe_provenance.get("cache_sha256") != sha256(universe_path)
    ):
        fail("CATlas fixed-universe cache differs from policy-bound provenance")

    trait_path = root / spec["trait_cache_path_template"].format(trait_id=task["trait_id"])
    trait_provenance_path = root / spec["trait_cache_provenance_path_template"].format(trait_id=task["trait_id"])
    if not trait_path.is_file() or not trait_provenance_path.is_file():
        fail("CATlas trait cache and provenance are incomplete")
    trait_provenance = json.loads(trait_provenance_path.read_text(encoding="utf-8"))
    if (
        trait_provenance.get("trait_id") != task["trait_id"]
        or trait_provenance.get("trait_policy_sha256") != trait_policy_sha256(policy)
        or trait_provenance.get("fixed_universe_provenance_sha256") != sha256(universe_provenance_path)
        or trait_provenance.get("cache_sha256") != sha256(trait_path)
    ):
        fail("CATlas trait cache differs from policy-bound provenance")

    universe_fields, universe_rows = read_tsv(universe_path)
    trait_fields, trait_rows = read_tsv(trait_path)
    if universe_fields != spec["variant_cache_fields"] or trait_fields != spec["trait_cache_fields"]:
        fail("CATlas cache schema differs from policy")
    normalized, counts = compute_enrichment_rows(
        policy, source_manifest, task, universe_rows, trait_rows,
    )
    fields = policy["cell_types"]["canonical_fields"]
    if any(list(row) != fields for row in normalized):
        fail("CATlas normalized row differs from the locked cell-type schema")
    payload = table_text(fields, normalized)
    out_path = root / args.out
    atomic_text(out_path, payload)
    provenance_path = root / args.provenance_out
    provenance = {
        "schema_version": policy["schema_version"],
        "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"],
        "terminal_status": "COMPLETED",
        "terminal_reason": (
            f"CATlas adult scATAC hypergeometric enrichment completed for all "
            f"{len(normalized)} prespecified {task['domain']} cell types"
        ),
        "policy_sha256": sha256(policy_path),
        "task_manifest_sha256": sha256(task_manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"],
        "source_id": task["source_id"],
        "source_release": task["source_release"],
        "source_manifest_sha256": sha256(source_manifest_path),
        "fixed_universe_provenance_sha256": sha256(universe_provenance_path),
        "trait_cache_provenance_sha256": sha256(trait_provenance_path),
        "counts": counts,
        "result_rows": len(normalized),
        "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "multiple_testing_family": spec["multiple_testing_family"],
        "claim_limit": spec["claim_limit"],
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(
        f"CATLAS_TASK_OK task={task['task_id']} rows={len(normalized)} "
        f"signals={counts['genome_wide_significant_background_variants']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
