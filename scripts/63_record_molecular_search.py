#!/usr/bin/env python3
"""Record a checksum-bound non-tabix molecular search or terminal coverage outcome."""
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
import tempfile
from datetime import datetime, timezone
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
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def gzip_text(path: Path, text: str) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def table_text(fields: list[str], rows: list[dict[str, str]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("search_task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--plan", default="results/tables/molecular_search_plan.tsv")
    parser.add_argument("--plan-lock", default="results/tables/molecular_search_plan.lock.json")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--variants", default="results/atlas/variants.tsv")
    parser.add_argument("--source-snapshot", type=Path, required=True)
    parser.add_argument("--normalized-qtl", type=Path)
    parser.add_argument("--terminal-outcome", choices=["NO_FEATURES_IN_LOCUS", "NO_ANALYZABLE_FEATURES", "ACCESS_BLOCKED", "NOT_APPLICABLE"])
    parser.add_argument("--note", required=True)
    parser.add_argument("--curator", required=True)
    parser.add_argument("--out-dir", default="results/molecular/search")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        fail("recording source evidence requires explicit --execute")
    if (args.normalized_qtl is None) == (args.terminal_outcome is None):
        fail("supply exactly one of --normalized-qtl or --terminal-outcome")
    root = Path(args.root).resolve()
    plan_path, lock_path = root / args.plan, root / args.plan_lock
    policy_path, variants_path = root / args.policy, root / args.variants
    _, plan = read_tsv(plan_path)
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("plan_sha256") != sha256(plan_path) or lock.get("results_accessed_before_lock") is not False:
        fail("molecular search plan differs from its result-free lock")
    selected = [row for row in plan if row["search_task_id"] == args.search_task_id]
    if len(selected) != 1:
        fail("search_task_id must identify exactly one locked source search")
    task = selected[0]
    if task["query_mode"] == "TABIX_GRCH38_INTERVAL":
        fail("use 63_query_eqtl_catalogue.py for locked tabix tasks")
    if not args.source_snapshot.is_file() or args.source_snapshot.stat().st_size == 0:
        fail("source snapshot is missing or empty")
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    schema = policy["normalized_qtl_schema"]
    rows: list[dict[str, str]] = []
    if args.normalized_qtl is not None:
        if not args.normalized_qtl.is_file() or args.normalized_qtl.stat().st_size == 0:
            fail("normalized QTL snapshot is missing or empty")
        fields, rows = read_tsv(args.normalized_qtl)
        if fields != schema or not rows:
            fail("normalized QTL snapshot has the wrong schema or no rows")
        _, variants = read_tsv(variants_path)
        reference = {(row["variant_id"], row["chromosome"], row["position_bp"]) for row in variants if row["locus_id"] == task["locus_id"]}
        seen: set[tuple[str, str]] = set()
        for row in rows:
            identity = (row["feature_id"], row["rsid"])
            if identity in seen:
                fail(f"duplicate normalized feature/variant row: {identity}")
            seen.add(identity)
            if (row["rsid"], row["chromosome"], row["position"]) not in reference:
                fail(f"normalized QTL variant is outside the locked locus universe: {identity}")
            if row["source_dataset_id"] != task["dataset_id"]:
                fail(f"normalized QTL dataset differs from locked task: {identity}")
            if task["modality"] == "eQTL_or_sQTL":
                if row["modality"] not in {"eQTL", "sQTL"}:
                    fail(f"invalid exact PsychENCODE modality: {identity}")
            elif row["modality"] != task["modality"]:
                fail(f"normalized QTL modality differs from locked task: {identity}")
            try:
                beta, se, maf, n = map(float, (row["beta"], row["se"], row["maf"], row["n"]))
            except ValueError as exc:
                fail(f"nonnumeric normalized QTL statistics: {identity}")
                raise AssertionError from exc
            if not all(math.isfinite(value) for value in (beta, se, maf, n)) or se <= 0 or not policy["variant_qc"]["minimum_maf"] <= maf <= 0.5 or n <= 0:
                fail(f"normalized QTL statistics fail QC: {identity}")
        outcome = "DATA_READY"
    else:
        outcome = args.terminal_outcome
    output_dir = root / args.out_dir / args.search_task_id
    if output_dir.exists():
        fail(f"search output already exists; refusing overwrite: {output_dir}")
    work_root = root / "work/molecular"
    work_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=args.search_task_id + ".", dir=work_root))
    try:
        snapshot_out = staging / "source_snapshot.bin"
        shutil.copyfile(args.source_snapshot, snapshot_out)
        normalized_out = staging / "normalized_qtl.tsv.gz"
        gzip_text(normalized_out, table_text(schema, rows))
        evidence = {
            "source_url": task["source_url"], "source_index_url": task["source_index_url"],
            "exact_release": task["exact_release"], "search_outcome": outcome,
            "note": args.note, "curator": args.curator,
        }
        evidence_out = staging / "search_evidence.json"
        evidence_out.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        provenance = {
            "schema_version": "atlas-v1.0-molecular-query.1",
            "queried_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "search_task_id": args.search_task_id, "search_plan_sha256": sha256(plan_path),
            "search_plan_lock_sha256": sha256(lock_path), "policy_sha256": sha256(policy_path),
            "source_family_id": task["source_family_id"], "study_id": task["study_id"],
            "dataset_id": task["dataset_id"], "exact_release": task["exact_release"],
            "source_url": task["source_url"], "source_index_url": task["source_index_url"],
            "search_outcome": outcome, "raw_source_row_count": "NOT_MACHINE_PARSED",
            "normalized_row_count": len(rows), "normalized_feature_count": len({row["feature_id"] for row in rows}),
            "outputs": {
                "source_snapshot.bin": {"bytes": snapshot_out.stat().st_size, "sha256": sha256(snapshot_out)},
                "normalized_qtl.tsv.gz": {"bytes": normalized_out.stat().st_size, "sha256": sha256(normalized_out)},
                "search_evidence.json": {"bytes": evidence_out.stat().st_size, "sha256": sha256(evidence_out)},
            },
            "results_accessed_before_lock": False, "curator": args.curator,
            "claim_limit": policy["claim_limit"],
        }
        (staging / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staging, output_dir)
        print(f"MOLECULAR_SEARCH_RECORDED task={args.search_task_id} outcome={outcome} rows={len(rows)}")
        return 0
    finally:
        if staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
