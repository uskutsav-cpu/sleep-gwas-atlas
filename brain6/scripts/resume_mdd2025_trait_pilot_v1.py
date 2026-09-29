#!/usr/bin/env python3
"""Resume the four interrupted, wholly unfinished MDD2025 pilot workers.

The original launch and logs remain immutable. A new attempt gets separate
logs and receipts; the pinned worker configs write only their absent outputs.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from run_mdd2025_trait_pilot_v1 import (
    PILOT, ROOT, RULE, RSCRIPT, RUNNER, save_new, sha, utc, validate_inputs,
)

EMPTY_SHA256 = sha(Path("/dev/null"))


def main() -> None:
    if len(sys.argv) != 2 or not sys.argv[1].isdigit() or int(sys.argv[1]) < 1:
        raise ValueError("Usage: resume_mdd2025_trait_pilot_v1.py ATTEMPT_NUMBER")
    attempt_number = int(sys.argv[1])
    attempt = f"resume_attempt_{attempt_number}"
    config = json.loads(PILOT.read_text())
    output = Path(config["output_root"])
    results = output / "results"
    materialization, workers = validate_inputs(config, output)
    if materialization.get("status") != "PASS_MATERIALIZED_DIAGNOSTIC_INPUTS":
        raise ValueError("Materialization is not valid")
    if config["worker_count"] != 4 or [w["data"]["expected_loci"] for w in workers] != [28, 28, 28, 27]:
        raise ValueError("Pinned four-worker partition changed")
    original_launch = results / "launch.json"
    launch = json.loads(original_launch.read_text())
    expected = {"pilot_config_sha256": sha(PILOT), "rule_sha256": sha(RULE),
                "materialization_receipt_sha256": sha(output / "materialization.receipt.json"),
                "r_runner_sha256": sha(RUNNER), "rscript_sha256": sha(RSCRIPT)}
    for key, value in expected.items():
        if launch.get(key) != value:
            raise ValueError(f"Original launch {key} mismatch")
    if len(launch["workers"]) != 4:
        raise ValueError("Original launch worker count mismatch")
    for item in workers:
        worker = item["data"]["worker_id"]
        original = launch["workers"][worker - 1]
        if original["worker_id"] != worker or original["config_sha256"] != item["sha256"]:
            raise ValueError(f"Original worker {worker} config mismatch")
        original_log = results / f"worker_{worker}.log"
        if not original_log.is_file() or sha(original_log) != EMPTY_SHA256:
            raise ValueError(f"Original worker {worker} log is not empty; inspect it manually")
        for suffix in ("tsv", "summary.json"):
            if (results / f"worker_{worker}.{suffix}").exists():
                raise FileExistsError(f"Worker {worker} {suffix} exists; inspect it manually")
    if (results / "exit.json").exists() or (results / "pilot_decision.json").exists():
        raise FileExistsError("An original terminal receipt or decision exists")
    previous_attempt_sha256 = None
    if attempt_number > 1:
        previous_exit = results / f"resume_attempt_{attempt_number - 1}.exit.json"
        previous = json.loads(previous_exit.read_text())
        if (previous.get("analysis_id") != config["analysis_id"]
                or previous.get("attempt") != f"resume_attempt_{attempt_number - 1}"
                or set(previous.get("exit_codes", {})) != {"1", "2", "3", "4"}
                or not all(code != 0 for code in previous["exit_codes"].values())
                or any(data["tsv_sha256"] is not None or data["summary_sha256"] is not None
                       for data in previous["worker_artifacts"].values())):
            raise ValueError("Previous attempt did not terminate cleanly with zero result files")
        previous_attempt_sha256 = sha(previous_exit)
    runtime_probe = subprocess.run(
        [str(RSCRIPT), "-e", "cat(as.character(packageVersion('LAVA')))"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    if runtime_probe.stdout.strip() != "0.1.5":
        raise ValueError("Pinned LAVA runtime did not load correctly")
    preflight_path = results / f"{attempt}.preflight.json"
    launch_path = results / f"{attempt}.launch.json"
    exit_path = results / f"{attempt}.exit.json"
    planned_logs = {str(i): results / f"worker_{i}.{attempt}.log" for i in range(1, 5)}
    for path in (preflight_path, launch_path, exit_path, *planned_logs.values()):
        if path.exists():
            raise FileExistsError(f"Existing resume attempt artifact: {path}")
    commands = {str(w["data"]["worker_id"]):
                [str(RSCRIPT), str(RUNNER), str(w["path"])] for w in workers}
    save_new(preflight_path, {
        "analysis_id": config["analysis_id"], "attempt": attempt,
        "audited_utc": utc(), "original_launch_sha256": sha(original_launch),
        "previous_attempt_exit_sha256": previous_attempt_sha256,
        "runtime_probe": {"r_version": "4.3.3", "lava_version": runtime_probe.stdout.strip(),
                          "prefix_path": "/private/tmp/brain6ext/r-env-4.3.3-short",
                          "prefix_resolved": str(Path("/private/tmp/brain6ext/r-env-4.3.3-short").resolve())},
        "original_log_sha256": {str(i): EMPTY_SHA256 for i in range(1, 5)},
        "completed_loci_before_resume": {str(i): 0 for i in range(1, 5)},
        "expected_loci": {str(w["data"]["worker_id"]): w["data"]["expected_loci"] for w in workers},
        "pinned_sha256": expected,
        "worker_config_sha256": {str(w["data"]["worker_id"]): w["sha256"] for w in workers},
        "commands": commands,
        "status": "ALL_FOUR_ORIGINAL_WORKERS_UNFINISHED_AND_RESULT_FILES_ABSENT",
        "scientific_scope": "Trait-only diagnostic under a sample-size interpretation hold; no family promotion",
    })
    env = dict(os.environ)
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    handles = []
    children = {}
    try:
        for item in workers:
            worker = str(item["data"]["worker_id"])
            handle = planned_logs[worker].open("xb")
            handles.append(handle)
            child = subprocess.Popen(commands[worker], cwd=ROOT, env=env,
                                     stdout=handle, stderr=subprocess.STDOUT)
            children[worker] = child
        save_new(launch_path, {
            "analysis_id": config["analysis_id"], "attempt": attempt,
            "launched_utc": utc(), "preflight_sha256": sha(preflight_path),
            "workers": {worker: {"pid": child.pid, "command": commands[worker],
                                 "log": str(planned_logs[worker])}
                        for worker, child in children.items()},
        })
        print(json.dumps({"status": "FOUR_UNFINISHED_WORKERS_RESUMED", "attempt": attempt,
                          "pids": {worker: child.pid for worker, child in children.items()},
                          "output": str(results)}, sort_keys=True), flush=True)
        codes = {worker: child.wait() for worker, child in children.items()}
    except BaseException:
        for child in children.values():
            if child.poll() is None:
                child.terminate()
        for child in children.values():
            child.wait()
        raise
    finally:
        for handle in handles:
            handle.close()
    save_new(exit_path, {
        "analysis_id": config["analysis_id"], "attempt": attempt,
        "finished_utc": utc(), "launch_sha256": sha(launch_path),
        "exit_codes": codes,
        "worker_artifacts": {worker: {
            "log_sha256": sha(planned_logs[worker]),
            "tsv_sha256": sha(results / f"worker_{worker}.tsv") if (results / f"worker_{worker}.tsv").is_file() else None,
            "summary_sha256": sha(results / f"worker_{worker}.summary.json") if (results / f"worker_{worker}.summary.json").is_file() else None,
        } for worker in children},
    })
    print(json.dumps({"status": "RESUME_WORKERS_TERMINAL", "exit_codes": codes,
                      "exit_receipt": str(exit_path)}, sort_keys=True), flush=True)
    if any(code != 0 for code in codes.values()):
        raise RuntimeError(f"One or more resumed workers failed: {codes}")


if __name__ == "__main__":
    main()
