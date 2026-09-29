#!/usr/bin/env python3
"""Inspect an interrupted MDD2025 pilot without restarting any worker.

When every pinned shard is terminal and valid, this may write the first
aggregate/decision. It never replaces an existing result or launches work.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from run_mdd2025_trait_pilot_v1 import PILOT, ROOT, RULE, aggregate, sha, table, validate_inputs


SAMPLE_SIZE_AUDIT = (ROOT / "brain6/results/lava_multitrait_feasibility_v1"
                     / "mdd2025_pilot_sample_size_audit_20260927.md")
LOWER_BOUNDS = (ROOT / "brain6/results/lava_multitrait_feasibility_v1"
                / "lower_bounds.json")


def execution_state(config: dict, output: Path, worker_configs: list[dict]) -> dict:
    """Verify the original or latest isolated attempt without changing files."""
    results = output / "results"
    original_launch = results / "launch.json"
    if not original_launch.is_file():
        raise ValueError("Original launch receipt is missing")
    original = json.loads(original_launch.read_text())
    expected = {"pilot_config_sha256": sha(PILOT), "rule_sha256": sha(RULE),
                "materialization_receipt_sha256": sha(output / "materialization.receipt.json"),
                "r_runner_sha256": config["r_runner_sha256"]}
    for key, value in expected.items():
        if original.get(key) != value:
            raise ValueError(f"Original launch {key} mismatch")
    workers = {str(item["data"]["worker_id"]): item for item in worker_configs}
    if len(original.get("workers", [])) != 4:
        raise ValueError("Original launch worker count mismatch")
    for entry in original["workers"]:
        worker = str(entry["worker_id"])
        if worker not in workers or entry["config_sha256"] != workers[worker]["sha256"]:
            raise ValueError(f"Original launch worker {worker} config mismatch")
    original_exit = results / "exit.json"
    if original_exit.is_file():
        receipt = json.loads(original_exit.read_text())
        if receipt["launch_sha256"] != sha(original_launch) or len(receipt["exit_codes"]) != 4:
            raise ValueError("Original exit receipt mismatch")
        for worker, log_hash in receipt["logs"].items():
            if log_hash != sha(results / f"worker_{worker}.log"):
                raise ValueError(f"Original worker {worker} log hash mismatch")
        return {"kind": "original", "attempt": 0,
                "exit_codes": {str(i): code for i, code in enumerate(receipt["exit_codes"], 1)},
                "exit_codes_verified": all(code == 0 for code in receipt["exit_codes"]),
                "logs": {w: results / f"worker_{w}.log" for w in workers}}
    attempts = {}
    for path in results.glob("resume_attempt_*.launch.json"):
        match = re.fullmatch(r"resume_attempt_(\d+)\.launch\.json", path.name)
        if match:
            attempts[int(match.group(1))] = path
    if not attempts:
        return {"kind": "original", "attempt": 0, "exit_codes": None,
                "exit_codes_verified": False,
                "logs": {w: results / f"worker_{w}.log" for w in workers}}
    if set(attempts) != set(range(1, max(attempts) + 1)):
        raise ValueError("Resume attempt sequence has gaps")
    previous_exit_hash = None
    current = None
    for number in sorted(attempts):
        label = f"resume_attempt_{number}"
        preflight_path = results / f"{label}.preflight.json"
        launch_path = attempts[number]
        preflight = json.loads(preflight_path.read_text())
        launch = json.loads(launch_path.read_text())
        if (preflight["original_launch_sha256"] != sha(original_launch)
                or preflight.get("previous_attempt_exit_sha256") != previous_exit_hash
                or preflight["pinned_sha256"]["pilot_config_sha256"] != sha(PILOT)
                or preflight["pinned_sha256"]["rule_sha256"] != sha(RULE)
                or preflight["worker_config_sha256"] != {w: item["sha256"] for w, item in workers.items()}
                or launch["preflight_sha256"] != sha(preflight_path)
                or set(launch["workers"]) != set(workers)):
            raise ValueError(f"Resume attempt {number} preflight/launch mismatch")
        exit_path = results / f"{label}.exit.json"
        logs = {w: results / f"worker_{w}.{label}.log" for w in workers}
        if exit_path.is_file():
            receipt = json.loads(exit_path.read_text())
            if (receipt["launch_sha256"] != sha(launch_path)
                    or set(receipt["exit_codes"]) != set(workers)
                    or set(receipt["worker_artifacts"]) != set(workers)):
                raise ValueError(f"Resume attempt {number} exit receipt mismatch")
            for worker, artifacts in receipt["worker_artifacts"].items():
                if artifacts["log_sha256"] != sha(logs[worker]):
                    raise ValueError(f"Resume attempt {number} worker {worker} log hash mismatch")
                for label_part, filename in (("tsv_sha256", f"worker_{worker}.tsv"),
                                             ("summary_sha256", f"worker_{worker}.summary.json")):
                    recorded = artifacts[label_part]
                    if recorded is not None and recorded != sha(results / filename):
                        raise ValueError(f"Resume attempt {number} worker {worker} output hash mismatch")
            previous_exit_hash = sha(exit_path)
            codes = receipt["exit_codes"]
            verified = all(code == 0 for code in codes.values())
        else:
            if number != max(attempts):
                raise ValueError(f"Resume attempt {number} lacks a terminal receipt")
            codes = None
            verified = False
        current = {"kind": "resume", "attempt": number, "exit_codes": codes,
                   "exit_codes_verified": verified, "logs": logs}
    assert current is not None
    return current


def main() -> None:
    config = json.loads(PILOT.read_text())
    rule = json.loads(RULE.read_text())
    output = Path(config["output_root"])
    if not SAMPLE_SIZE_AUDIT.is_file():
        raise ValueError("Missing prospective sample-size interpretation audit")
    bounds = json.loads(LOWER_BOUNDS.read_text())
    remaining_not_run = (bounds["canonical_not_run"]
                         - bounds["canonical_not_run_by_trait"]["mdd"])
    if (bounds["canonical_family_cells"] != 17465
            or bounds["frozen_maximum_not_run"] != 873):
        raise ValueError("Frozen family size or NOT_RUN ceiling changed")
    if remaining_not_run <= bounds["frozen_maximum_not_run"]:
        raise ValueError("MDD-only feasibility bound no longer supports this hold")
    interpretation = {"status": "SCIENTIFIC_INTERPRETATION_HOLD",
                      "audit_sha256": sha(SAMPLE_SIZE_AUDIT),
                      "full_trait_only_screen_authorized": False,
                      "mdd_only_best_possible_family_not_run": remaining_not_run,
                      "frozen_family_not_run_ceiling": bounds["frozen_maximum_not_run"],
                      "lower_bounds_sha256": sha(LOWER_BOUNDS)}
    if not output.is_dir():
        print(json.dumps({"status": "BLOCKED_EXTERNAL_SSD_UNMOUNTED", "expected_root": str(output),
                          "interpretation": interpretation}, sort_keys=True))
        return
    materialization, worker_configs = validate_inputs(config, output)
    execution = execution_state(config, output, worker_configs)
    worker_state = {}
    complete = execution["exit_codes_verified"]
    for item in worker_configs:
        worker = item["data"]["worker_id"]
        tsv = output / "results" / f"worker_{worker}.tsv"
        summary = output / "results" / f"worker_{worker}.summary.json"
        log = execution["logs"][str(worker)]
        paths = {"tsv": tsv, "summary": summary, "log": log}
        present = {name: path.is_file() for name, path in paths.items()}
        terminal = execution["exit_codes_verified"] and all(present.values())
        rows = []
        if terminal:
            result = json.loads(summary.read_text())
            rows = table(tsv)
            expected_loci = set(item["data"]["locus_ids"])
            terminal = (result.get("rows") == item["data"]["expected_loci"]
                        and result.get("failed") == 0 and len(rows) == len(expected_loci)
                        and {r["locus_id"] for r in rows} == expected_loci
                        and all(r["trait_id"] == "mdd2025_no23andme" and
                                r["status"] in {"TESTED", "NOT_RUN"} for r in rows)
                        and result.get("tested") == sum(r["status"] == "TESTED" for r in rows)
                        and result.get("not_run") == sum(r["status"] == "NOT_RUN" for r in rows)
                        and result.get("strict_gate_pass") == sum(r["strict_gate_pass"].upper() == "TRUE" for r in rows)
                        and "SCREEN_COMPLETE" in log.read_text(errors="replace"))
        if not terminal:
            complete = False
        worker_state[str(worker)] = {"present": present, "terminal": terminal,
                                     "completed_loci": len(rows) if terminal else 0,
                                     "expected_loci": item["data"]["expected_loci"],
                                     "sha256": {name: sha(path) for name, path in paths.items() if path.is_file()}}
    decision_path = output / "results" / "pilot_decision.json"
    if decision_path.exists():
        decision = json.loads(decision_path.read_text())
        if not complete or decision["pilot_config_sha256"] != sha(PILOT) or decision["decision_rule_sha256"] != sha(RULE):
            raise ValueError("Existing pilot decision conflicts with current worker state or frozen rules")
        if decision["aggregate_sha256"] != sha(output / "results" / "pilot_aggregate.tsv"):
            raise ValueError("Existing pilot aggregate checksum mismatch")
        if decision["comparison_sha256"] != sha(output / "results" / "pilot_vs_canonical.tsv"):
            raise ValueError("Existing pilot comparison checksum mismatch")
        for worker, checksums in decision["worker_outputs"].items():
            if (checksums["tsv_sha256"] != worker_state[worker]["sha256"].get("tsv")
                    or checksums["summary_sha256"] != worker_state[worker]["sha256"].get("summary")):
                raise ValueError(f"Existing pilot decision worker {worker} hash mismatch")
        status = "COMPLETE_EXISTING_DECISION_VERIFIED"
    elif complete:
        decision = aggregate(config, rule, output, worker_configs)
        status = "COMPLETE_RECOVERED_FROM_TERMINAL_SHARDS"
    else:
        decision = None
        status = "PARTIAL_SHARDS_NO_RESTART"
    print(json.dumps({"status": status, "workers": worker_state,
                      "decision": decision["status"] if decision else None,
                      "interpretation": interpretation,
                      "materialization_status": materialization["status"],
                      "execution": {k: v for k, v in execution.items() if k != "logs"}}, sort_keys=True))


if __name__ == "__main__":
    main()
