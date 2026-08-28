#!/usr/bin/env python3
"""Run one checksum-locked conjunction-FDR task with the pinned MATLAB code."""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import pleiotropy_contract


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def one_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(f"ERROR: expected one row in {path}")
    return rows[0]


def atomic_tsv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(path)


def physical_memory() -> int:
    return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task")
    parser.add_argument("lock")
    parser.add_argument("--root", default=".")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        raise SystemExit("ERROR: real conjunction-FDR compute requires the explicit --execute flag")
    root = Path(args.root).resolve()
    task_path, lock_path = root / args.task, root / args.lock
    task, lock = one_row(task_path), one_row(lock_path)
    if lock.get("schema_version") != "sleep-atlas-conjfdr-task.1":
        raise SystemExit("ERROR: conjunction-FDR task lock has the wrong schema")
    if lock.get("task_sha256") != sha256(task_path):
        raise SystemExit("ERROR: conjunction-FDR task differs from its lock")
    if lock.get("analysis_id") != task.get("analysis_id") or lock.get("pair_id") != task.get("pair_id"):
        raise SystemExit("ERROR: conjunction-FDR task identity differs from its lock")
    if (
        task.get("task_builder_sha256") != sha256(root / "scripts/48_prepare_conjfdr_task.py")
        or task.get("runner_sha256") != sha256(Path(__file__))
    ):
        raise SystemExit("ERROR: conjunction-FDR task creation or execution code drifted")
    runtime_provenance = pleiotropy_contract.validate_runtime(root)
    runtime_path = root / task.get("runtime_provenance", "")
    reference_record = runtime_provenance["pleiofdr_reference"]
    if (
        task.get("runtime_provenance") != "ref/pleiofdr/runtime.provenance.json"
        or task.get("runtime_provenance_sha256") != sha256(runtime_path)
        or task.get("reference_sha256") != reference_record["sha256"]
        or int(task.get("reference_bytes", "0")) != reference_record["bytes"]
    ):
        raise SystemExit("ERROR: conjunction-FDR runtime/reference provenance drifted")
    required_files = (
        ("config", "config_sha256"), ("sleep_mat", "sleep_mat_sha256"),
        ("non_sleep_mat", "non_sleep_mat_sha256"), ("template", "template_sha256"),
        ("overlap_patch", "overlap_patch_sha256"),
    )
    for field, hash_field in required_files:
        path = root / task[field]
        if not path.is_file() or sha256(path) != task[hash_field]:
            raise SystemExit(f"ERROR: locked conjunction-FDR input drifted: {field}")
    code = root / task["pleiofdr_code"]
    commit = subprocess.run(
        ["git", "-C", str(code), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(code), "status", "--porcelain"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if commit != task["pleiofdr_commit"] or dirty:
        raise SystemExit("ERROR: pleioFDR checkout drifted")
    if task.get("correct_sample_overlap") != "TRUE":
        raise SystemExit("ERROR: task does not require the locked sample-overlap correction")
    matlab = shutil.which("matlab")
    if not matlab:
        raise SystemExit("ERROR: MATLAB is unavailable; experimental Octave is not accepted")
    if physical_memory() < 17179869184:
        raise SystemExit("ERROR: conjunction-FDR requires at least 16 GiB physical memory")
    if shutil.disk_usage(root).free < 21474836480:
        raise SystemExit("ERROR: conjunction-FDR requires at least 20 GiB free storage")

    outputs = [root / task[field] for field in ("all_results", "locus_results", "result_mat")]
    completion_path = root / task["completion"]
    if completion_path.exists() or any(path.exists() for path in outputs):
        raise SystemExit("ERROR: immutable conjunction-FDR result family already exists")

    runtime = root / "work/pleiotropy/conjfdr_runtime" / task["pair_id"] / "software"
    if runtime.exists():
        raise SystemExit(f"ERROR: clean runtime already exists: {runtime}")
    runtime.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(code, runtime, ignore=shutil.ignore_patterns(".git"))
    patch = root / task["overlap_patch"]
    patched = subprocess.run(
        ["patch", "-p1", "-i", str(patch)], cwd=runtime, capture_output=True, text=True
    )
    if patched.returncode:
        raise SystemExit("ERROR: could not apply locked sample-overlap patch: " + patched.stderr)
    config = root / task["config"]
    config_literal = str(config).replace("'", "''")
    runme_literal = str(runtime / "runme.m").replace("'", "''")
    command = [matlab, "-batch", f"config='{config_literal}'; run('{runme_literal}');"]
    result = subprocess.run(command, cwd=runtime, capture_output=True, text=True)
    result_dir = root / task["result_dir"]
    result_dir.mkdir(parents=True, exist_ok=True)
    log_path = result_dir / "matlab.log"
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise SystemExit(f"ERROR: MATLAB conjunction-FDR failed; see {log_path}")
    missing = [str(path) for path in outputs if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise SystemExit("ERROR: conjunction-FDR omitted expected outputs: " + ", ".join(missing))
    completion = {
        "analysis_id": task["analysis_id"], "pair_id": task["pair_id"],
        "analysis_status": "CONJFDR_COMPLETE", "task_sha256": sha256(task_path),
        "all_results_sha256": sha256(outputs[0]), "locus_results_sha256": sha256(outputs[1]),
        "result_mat_sha256": sha256(outputs[2]), "matlab_log_sha256": sha256(log_path),
        "correct_sample_overlap": "TRUE", "random_prune_iterations": task["random_prune_iterations"],
        "conjfdr_threshold": task["conjfdr_threshold"],
        "runtime_provenance_sha256": sha256(runtime_path),
        "task_builder_sha256": task["task_builder_sha256"],
        "runner_sha256": task["runner_sha256"],
    }
    atomic_tsv(completion_path, completion)
    print(f"CONJFDR_PAIR_OK pair={task['pair_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
