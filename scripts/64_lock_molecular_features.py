#!/usr/bin/env python3
"""Validate complete source-search coverage and lock every analyzable feature before coloc."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

import molecular_contract

COVERAGE_FIELDS = molecular_contract.SEARCH_COVERAGE_FIELDS
MANIFEST_FIELDS = molecular_contract.FEATURE_FIELDS
EXCLUSION_FIELDS = molecular_contract.FEATURE_EXCLUSION_FIELDS


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, compressed: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty input: {path}")
    opener = gzip.open if compressed or path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def write_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def safe_path(root: Path, path: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        fail(f"path escapes repository root: {path}")
    return resolved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--plan", default="results/tables/molecular_search_plan.tsv")
    parser.add_argument("--plan-lock", default="results/tables/molecular_search_plan.lock.json")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--search-dir", default="results/molecular/search")
    parser.add_argument("--run-dir", default="results/molecular/coloc_runs")
    parser.add_argument("--coverage-out", default="results/tables/molecular_search_coverage.tsv")
    parser.add_argument("--manifest-out", default="results/tables/molecular_feature_manifest.tsv")
    parser.add_argument("--exclusions-out", default="results/tables/molecular_feature_exclusions.tsv")
    parser.add_argument("--lock-out", default="results/tables/molecular_feature_manifest.lock.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    plan_path, plan_lock_path, policy_path = root / args.plan, root / args.plan_lock, root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    plan, plan_lock = molecular_contract.validate_search_plan(
        root, plan_path, plan_lock_path, policy_path,
        root / "results/tables/molecular_preflight.json",
    )
    run_dir = root / args.run_dir
    if not args.validate_only and run_dir.exists() and any(run_dir.iterdir()):
        fail("trait-molecular result directory is nonempty before feature-family lock")
    schema = policy["normalized_qtl_schema"]
    terminal = set(policy["search_family"]["terminal_outcomes"])
    minimum = int(policy["variant_qc"]["minimum_shared_variants"])
    maximum = int(policy["variant_qc"]["maximum_shared_variants"])
    coverage: list[dict[str, object]] = []
    manifest: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    input_hashes: dict[str, dict[str, str]] = {}
    comparison_ids: set[str] = set()
    for task in plan:
        task_id = task["search_task_id"]
        directory = safe_path(root, root / args.search_dir / task_id)
        provenance_path = directory / "provenance.json"
        normalized_path = directory / "normalized_qtl.tsv.gz"
        if not provenance_path.is_file() or not normalized_path.is_file():
            fail(f"molecular source search is incomplete: {task_id}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if (
            provenance.get("schema_version") != "atlas-v1.0-molecular-query.2"
            or provenance.get("search_task_id") != task_id
            or provenance.get("search_plan_sha256") != sha256(plan_path)
            or provenance.get("search_plan_lock_sha256") != sha256(plan_lock_path)
            or provenance.get("policy_sha256") != sha256(policy_path)
            or provenance.get("results_accessed_before_lock") is not False
            or provenance.get("script_sha256") != molecular_contract.script_hashes(root, "qtl")
            or "normalized_qtl.tsv.gz" not in provenance.get("outputs", {})
        ):
            fail(f"molecular query provenance differs from locks: {task_id}")
        for name, record in provenance.get("outputs", {}).items():
            path = directory / name
            if not path.is_file() or path.stat().st_size != record.get("bytes") or sha256(path) != record.get("sha256"):
                fail(f"molecular source snapshot checksum mismatch: {task_id}/{name}")
        fields, rows = read_tsv(normalized_path, compressed=True)
        if fields != schema:
            fail(f"normalized molecular schema drifted: {task_id}")
        outcome = provenance.get("search_outcome")
        if outcome != "DATA_READY" and outcome not in terminal:
            fail(f"invalid molecular search outcome: {task_id}/{outcome}")
        if outcome == "DATA_READY" and not rows:
            fail(f"DATA_READY query has no normalized rows: {task_id}")
        if outcome in terminal and rows:
            fail(f"terminal source search contains normalized evidence rows: {task_id}")
        by_feature: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            if row["source_dataset_id"] != task["dataset_id"]:
                fail(f"normalized molecular row has wrong dataset: {task_id}")
            by_feature.setdefault(row["feature_id"], []).append(row)
        analyzable = 0
        for feature_id in sorted(by_feature):
            feature_rows = by_feature[feature_id]
            variant_ids = [row["rsid"] for row in feature_rows]
            if len(variant_ids) != len(set(variant_ids)):
                fail(f"duplicate normalized variant within feature: {task_id}/{feature_id}")
            feature_name = {row["feature_name"] for row in feature_rows}
            gene_id = {row["gene_id"] for row in feature_rows}
            gene_symbol = {row["gene_symbol"] for row in feature_rows}
            context = {row["context"] for row in feature_rows}
            modality = {row["modality"] for row in feature_rows}
            if any(len(values) != 1 for values in (feature_name, gene_id, gene_symbol, context, modality)):
                fail(f"feature annotations are not invariant: {task_id}/{feature_id}")
            count = len(feature_rows)
            if count < minimum:
                exclusions.append({
                    "search_task_id": task_id, "locus_id": task["locus_id"],
                    "source_family_id": task["source_family_id"], "dataset_id": task["dataset_id"],
                    "modality": next(iter(modality)), "feature_id": feature_id,
                    "shared_variant_count": count, "exclusion_reason": f"BELOW_MINIMUM_SHARED_VARIANTS_{minimum}",
                    "query_provenance_path": str(provenance_path.relative_to(root)),
                })
                continue
            if count > maximum:
                fail(f"molecular feature exceeds locked maximum and requires a pre-result boundary rule: {task_id}/{feature_id}")
            n_values = []
            for row in feature_rows:
                try:
                    n = float(row["n"])
                except ValueError as exc:
                    fail(f"nonnumeric molecular N: {task_id}/{feature_id}")
                    raise AssertionError from exc
                if not math.isfinite(n) or n <= 0:
                    fail(f"invalid molecular N: {task_id}/{feature_id}")
                n_values.append(n)
            analyzable += 1
            trait_task_path = root / f"results/fine_mapping/tasks/{task['locus_id']}.tsv"
            trait_task_lock_path = root / f"results/fine_mapping/tasks/{task['locus_id']}.lock.json"
            if not trait_task_path.is_file() or not trait_task_lock_path.is_file():
                fail(f"fine-mapping input task is absent for molecular feature: {task['locus_id']}")
            for trait_id, trait_role in ((task["sleep_trait"], "SLEEP_TRAIT"), (task["non_sleep_trait"], "NON_SLEEP_TRAIT")):
                digest = hashlib.sha256(f"{task_id}\t{trait_id}\t{feature_id}".encode()).hexdigest()[:24]
                comparison_id = f"MQC__{digest}"
                if comparison_id in comparison_ids:
                    fail("molecular comparison hash collision")
                comparison_ids.add(comparison_id)
                manifest.append({
                    "comparison_id": comparison_id, "search_task_id": task_id,
                    "locus_id": task["locus_id"], "pair_id": task["pair_id"],
                    "trait_id": trait_id, "trait_role": trait_role,
                    "source_family_id": task["source_family_id"], "study_id": task["study_id"],
                    "dataset_id": task["dataset_id"], "exact_release": task["exact_release"],
                    "modality": next(iter(modality)), "feature_id": feature_id,
                    "feature_name": next(iter(feature_name)), "gene_id": next(iter(gene_id)),
                    "gene_symbol": next(iter(gene_symbol)), "context": next(iter(context)),
                    "molecular_N": format(statistics.median(n_values), ".15g"),
                    "shared_variant_count": count,
                    "normalized_qtl_path": str(normalized_path.relative_to(root)),
                    "normalized_qtl_sha256": sha256(normalized_path),
                    "trait_task_path": str(trait_task_path.relative_to(root)),
                    "trait_task_sha256": sha256(trait_task_path),
                    "trait_task_lock_path": str(trait_task_lock_path.relative_to(root)),
                    "trait_task_lock_sha256": sha256(trait_task_lock_path),
                    "trait_molecular_results_accessed_before_feature_lock": "NO",
                })
        if outcome in terminal:
            exclusions.append({
                "search_task_id": task_id, "locus_id": task["locus_id"],
                "source_family_id": task["source_family_id"], "dataset_id": task["dataset_id"],
                "modality": task["modality"], "feature_id": "NA", "shared_variant_count": 0,
                "exclusion_reason": outcome, "query_provenance_path": str(provenance_path.relative_to(root)),
            })
        coverage.append({
            "search_task_id": task_id, "locus_id": task["locus_id"],
            "source_family_id": task["source_family_id"], "study_id": task["study_id"],
            "dataset_id": task["dataset_id"], "modality": task["modality"], "context": task["context"],
            "search_outcome": outcome, "normalized_row_count": len(rows), "feature_count": len(by_feature),
            "analyzable_feature_count": analyzable,
            "query_provenance_path": str(provenance_path.relative_to(root)),
            "query_provenance_sha256": sha256(provenance_path),
        })
        input_hashes[task_id] = {"provenance": sha256(provenance_path), "normalized_qtl": sha256(normalized_path)}
    if [row["search_task_id"] for row in coverage] != plan_lock["search_task_ids_in_locked_order"]:
        fail("molecular source coverage differs from the locked query family/order")
    coverage_text = table_text(COVERAGE_FIELDS, coverage)
    manifest_text = table_text(MANIFEST_FIELDS, manifest)
    exclusion_text = table_text(EXCLUSION_FIELDS, exclusions)
    coverage_out, manifest_out = root / args.coverage_out, root / args.manifest_out
    exclusions_out, lock_out = root / args.exclusions_out, root / args.lock_out
    locked_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    if args.validate_only and lock_out.is_file():
        locked_utc = json.loads(lock_out.read_text(encoding="utf-8")).get("locked_utc", locked_utc)
    lock_payload = {
        "schema_version": "atlas-v1.0-molecular-features.2",
        "locked_utc": locked_utc,
        "trait_molecular_results_accessed_before_feature_lock": False,
        "search_task_count": len(coverage), "analyzable_feature_comparison_count": len(manifest),
        "comparison_ids_in_locked_order": [row["comparison_id"] for row in manifest],
        "coverage_sha256": hashlib.sha256(coverage_text.encode()).hexdigest(),
        "manifest_sha256": hashlib.sha256(manifest_text.encode()).hexdigest(),
        "exclusions_sha256": hashlib.sha256(exclusion_text.encode()).hexdigest(),
        "search_plan_sha256": sha256(plan_path),
        "search_plan_lock_sha256": sha256(plan_lock_path), "policy_sha256": sha256(policy_path),
        "search_input_hashes": input_hashes,
        "zero_locus_qtl_not_applicable": len(plan) == 0,
        "script_sha256": molecular_contract.script_hashes(root, "qtl"),
        "claim_limit": policy["claim_limit"],
    }
    lock_text = json.dumps(lock_payload, indent=2, sort_keys=True) + "\n"
    if args.validate_only:
        expected = (
            (coverage_out, coverage_text), (manifest_out, manifest_text),
            (exclusions_out, exclusion_text), (lock_out, lock_text),
        )
        for path, text in expected:
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                fail(f"molecular feature-family artifact drifted: {path}")
    else:
        if any(path.exists() for path in (coverage_out, manifest_out, exclusions_out, lock_out)):
            fail("molecular feature-family outputs already exist; refusing overwrite")
        write_atomic(coverage_out, coverage_text)
        write_atomic(manifest_out, manifest_text)
        write_atomic(exclusions_out, exclusion_text)
        write_atomic(lock_out, lock_text)
    if not args.quiet:
        print(f"MOLECULAR_FEATURE_LOCK_OK searches={len(coverage)} comparisons={len(manifest)} exclusions={len(exclusions)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
