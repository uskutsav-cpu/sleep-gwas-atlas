#!/usr/bin/env python3
"""Create one checksum-locked PLACO+ execution task from a materialized pair."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


PAIR_ID = re.compile(r"[a-z0-9_]+__[a-z0-9_]+$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def relative(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except ValueError as exc:
        raise SystemExit(f"ERROR: task path is outside repository root: {path}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/pleiotropy_pair_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/pleiotropy_pair_manifest.lock.json")
    parser.add_argument("--policy", default="config/pleiotropy_analysis_policy.json")
    parser.add_argument("--placo-source", default=".r-env/share/placo/PLACO_v0.2.0.R")
    parser.add_argument("--pair-dir", default="results/pleiotropy/inputs")
    parser.add_argument("--task-dir", default="results/pleiotropy/tasks")
    parser.add_argument("--result-dir", default="results/pleiotropy/placo")
    args = parser.parse_args()

    if not PAIR_ID.fullmatch(args.pair_id):
        raise SystemExit("ERROR: invalid pair_id")
    root = Path(args.root).resolve()
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    policy_path = root / args.policy
    placo_source = root / args.placo_source
    manifest_lock = json.loads(manifest_lock_path.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if (
        manifest_lock.get("schema_version") != "sleep-atlas-pleiotropy-pairs.1"
        or manifest_lock.get("ready_pair_count") != policy["expected_sleep_non_sleep_pairs"]
        or manifest_lock.get("blocked_pair_count") != 0
        or manifest_lock.get("script_sha256")
        != sha256(root / "scripts/43_prepare_pleiotropy_pairs.py")
    ):
        raise SystemExit("ERROR: complete immutable pleiotropy pair family is required")
    if manifest_lock.get("manifest_sha256") != sha256(manifest_path):
        raise SystemExit("ERROR: pleiotropy pair manifest differs from its lock")
    if manifest_lock.get("policy_sha256") != sha256(policy_path):
        raise SystemExit("ERROR: pleiotropy policy differs from the pair-manifest lock")
    if manifest_lock.get("analysis_id") != policy.get("analysis_id"):
        raise SystemExit("ERROR: pleiotropy analysis IDs differ")
    matches = [row for row in read_tsv(manifest_path) if row["pair_id"] == args.pair_id]
    if len(matches) != 1:
        raise SystemExit("ERROR: pair_id must identify exactly one locked pair")
    manifest_row = matches[0]
    if manifest_row["input_status"] != "READY_FULL_SUMSTATS":
        raise SystemExit(f"ERROR: pair is blocked: {manifest_row['blocker']}")

    pair_path = root / args.pair_dir / f"{args.pair_id}.tsv.gz"
    provenance_path = root / args.pair_dir / f"{args.pair_id}.provenance.json"
    if not pair_path.is_file() or not provenance_path.is_file():
        raise SystemExit(
            "ERROR: materialized pair/provenance is absent; run script 44 with --materialize first"
        )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    pair_hash = sha256(pair_path)
    checks = {
        "pair_id": provenance.get("pair_id") == args.pair_id,
        "manifest_sha256": provenance.get("manifest_sha256") == sha256(manifest_path),
        "policy_sha256": provenance.get("policy_sha256") == sha256(policy_path),
        "output_sha256": provenance.get("output_sha256") == pair_hash,
        "output": provenance.get("output") == relative(root, pair_path),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"ERROR: materialized-pair provenance failed: {', '.join(failed)}")
    input_rows = int(provenance.get("alignment_counts", {}).get("written", -1))
    if input_rows <= 0:
        raise SystemExit("ERROR: materialized pair has no locked eligible rows")
    if not placo_source.is_file() or sha256(placo_source) != policy["placo_source_sha256"]:
        raise SystemExit("ERROR: pinned PLACO+ source is absent or checksum-invalid")

    task_path = root / args.task_dir / f"{args.pair_id}.tsv"
    task_lock_path = root / args.task_dir / f"{args.pair_id}.lock.tsv"
    if task_path.exists() or task_lock_path.exists():
        raise SystemExit("ERROR: immutable PLACO+ task or lock already exists")
    result_base = root / args.result_dir / args.pair_id
    task = {
        "analysis_id": policy["analysis_id"],
        "pair_id": args.pair_id,
        "sleep_trait": manifest_row["sleep_trait"],
        "non_sleep_trait": manifest_row["non_sleep_trait"],
        "analysis_tier": manifest_row["analysis_tier"],
        "pair_input": relative(root, pair_path),
        "pair_input_sha256": pair_hash,
        "pair_input_rows": input_rows,
        "pair_provenance_sha256": sha256(provenance_path),
        "manifest_sha256": sha256(manifest_path),
        "policy_sha256": sha256(policy_path),
        "placo_source": relative(root, placo_source),
        "placo_source_sha256": policy["placo_source_sha256"],
        "marginal_p_threshold": policy["placo_marginal_p_threshold"],
        "z_squared_maximum": policy["placo_z_squared_maximum"],
        "abs_tolerance": policy["placo_abs_tolerance"],
        "maximum_failure_fraction": policy["placo_maximum_numerical_failure_fraction"],
        "conventional_threshold": policy["placo_conventional_variant_threshold"],
        "family_threshold": policy["placo_locked_pair_family_threshold"],
        "task_builder_sha256": sha256(Path(__file__)),
        "runner_sha256": sha256(root / "scripts/46_run_placo_pair.R"),
        "variant_hits_out": relative(root, result_base.with_suffix(".hits.tsv")),
        "summary_out": relative(root, result_base.with_suffix(".summary.tsv")),
    }
    atomic_tsv(task_path, [task])
    task_lock = {
        "schema_version": "sleep-atlas-placo-task.1",
        "analysis_id": policy["analysis_id"],
        "pair_id": args.pair_id,
        "task": relative(root, task_path),
        "task_sha256": sha256(task_path),
    }
    atomic_tsv(task_lock_path, [task_lock])
    print(f"PLACO+ task locked: {args.pair_id} ({input_rows} variants)")
    print(f"Run explicitly: .r-env/bin/Rscript scripts/46_run_placo_pair.R {relative(root, task_path)} {relative(root, task_lock_path)} --execute")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
