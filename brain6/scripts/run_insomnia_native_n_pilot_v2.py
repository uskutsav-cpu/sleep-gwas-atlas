#!/usr/bin/env python3
"""Run and audit four isolated insomnia native-N LAVA pilot workers."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PILOT = ROOT / "brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_v2.json"
RULE = ROOT / "brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_decision_rule_v2.json"
RUNNER = ROOT / "brain6/scripts/run_insomnia_native_n_pilot_v2.R"
CANONICAL = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_new(path: Path, data: object) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write("\n")


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_inputs(config: dict, output: Path) -> tuple[dict, list[dict]]:
    if sha(Path(__file__)) != config["python_runner_sha256"]:
        raise ValueError("Pinned Python runner changed")
    if sha(RUNNER) != config["r_runner_sha256"]:
        raise ValueError("Pinned R runner changed")
    if len(config["locus_ids"]) != config["expected_loci"] or config["worker_count"] != 4 or config["expected_loci"] != 88:
        raise ValueError("Frozen pilot family or worker count changed")
    if sha(RSCRIPT) != config["rscript_sha256"]:
        raise ValueError("Pinned Rscript changed")
    if sha(CANONICAL) != config["canonical_aggregate_sha256"]:
        raise ValueError("Canonical aggregate changed")
    receipt_path = output / "materialization.receipt.json"
    receipt = json.loads(receipt_path.read_text())
    if receipt["config_sha256"] != sha(PILOT) or set(receipt["records"]) != set(config["locus_ids"]):
        raise ValueError("Materialization receipt does not cover exact frozen pilot")
    for loc, record in receipt["records"].items():
        directory = output / "inputs" / f"locus_{loc}"
        if sha(directory / "insomnia.sumstats.tsv.gz") != record["sumstats_sha256"]:
            raise ValueError(f"Pilot sumstats hash mismatch: locus {loc}")
        if sha(directory / "input_info.tsv") != record["input_info_sha256"]:
            raise ValueError(f"Pilot input info hash mismatch: locus {loc}")
    configs = []
    found = []
    for worker in range(1, 5):
        path = output / f"worker_{worker}.config.json"
        data = json.loads(path.read_text())
        if data["worker_id"] != worker or data["pilot_config_sha256"] != sha(PILOT) or data["materialization_receipt_sha256"] != sha(receipt_path):
            raise ValueError(f"Worker {worker} config identity mismatch")
        if data["expected_loci"] != len(data["locus_ids"]):
            raise ValueError(f"Worker {worker} locus count mismatch")
        found.extend(data["locus_ids"])
        configs.append({"path": path, "sha256": sha(path), "data": data})
    if len(found) != len(set(found)) or set(found) != set(config["locus_ids"]):
        raise ValueError("Worker partitions do not cover pilot exactly once")
    return receipt, configs


def aggregate(config: dict, rule: dict, output: Path, worker_configs: list[dict]) -> dict:
    combined = []
    for item in worker_configs:
        worker = item["data"]["worker_id"]
        tsv = output / "results" / f"worker_{worker}.tsv"
        summary = output / "results" / f"worker_{worker}.summary.json"
        rows = table(tsv)
        stats = json.loads(summary.read_text())
        if len(rows) != len(item["data"]["locus_ids"]) or set(r["locus_id"] for r in rows) != set(item["data"]["locus_ids"]):
            raise ValueError(f"Worker {worker} result partition mismatch")
        if stats["rows"] != len(rows) or stats["failed"] != 0:
            raise ValueError(f"Worker {worker} incomplete or failed")
        combined.extend(rows)
    if len(combined) != config["expected_loci"] or len({r["locus_id"] for r in combined}) != len(combined):
        raise ValueError("Combined pilot results incomplete")
    for row in combined:
        if row["trait_id"] != "insomnia" or row["status"] not in {"TESTED", "NOT_RUN"}:
            raise ValueError("Unexpected pilot trait/status")
        if row["status"] == "TESTED":
            p = float(row["p"])
            if not math.isfinite(p) or not 0 <= p <= 1:
                raise ValueError("Invalid pilot univariate p")
            if (p < config["strict_gate_p"]) != (row["strict_gate_pass"].upper() == "TRUE"):
                raise ValueError("Frozen strict gate mismatch")
    if sha(CANONICAL) != rule["canonical_aggregate_sha256"]:
        raise ValueError("Canonical aggregate hash changed")
    canonical = {r["locus_id"]: r for r in table(CANONICAL) if r["phen"] == "insomnia" and r["locus_id"] in config["locus_ids"]}
    if set(canonical) != set(config["locus_ids"]):
        raise ValueError("Canonical comparison lacks pilot loci")
    pilot = {r["locus_id"]: r for r in combined}
    canon_tested = sum(r["status"] == "TESTED" for r in canonical.values())
    pilot_tested = sum(r["status"] == "TESTED" for r in pilot.values())
    canon_low_h2 = sum(r["reason"] == "LOW_LOCAL_H2_UNDERPOWERED" for r in canonical.values())
    pilot_low_h2 = sum(r["reason"] == "LOW_LOCAL_H2_UNDERPOWERED" for r in pilot.values())
    gain = (pilot_tested - canon_tested) / len(pilot)
    reduction = (canon_low_h2 - pilot_low_h2) / canon_low_h2 if canon_low_h2 else None
    rule_pass = gain >= rule["advance_if_any"]["absolute_tested_fraction_gain_at_least"] or (
        reduction is not None and reduction >= rule["advance_if_any"]["relative_low_local_h2_not_run_reduction_at_least"])
    ordered = [pilot[loc] for loc in config["locus_ids"]]
    result_path = output / "results" / "pilot_aggregate.tsv"
    with result_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ordered[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ordered)
    comparison_path = output / "results" / "pilot_vs_canonical.tsv"
    with comparison_path.open("x", encoding="utf-8", newline="") as stream:
        fields = ("locus_id", "canonical_status", "canonical_reason", "canonical_p", "pilot_status",
                  "pilot_reason", "pilot_p", "pilot_strict_gate_pass")
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for loc in config["locus_ids"]:
            writer.writerow({"locus_id": loc, "canonical_status": canonical[loc]["status"],
                             "canonical_reason": canonical[loc]["reason"], "canonical_p": canonical[loc]["p"],
                             "pilot_status": pilot[loc]["status"], "pilot_reason": pilot[loc]["reason"],
                             "pilot_p": pilot[loc]["p"], "pilot_strict_gate_pass": pilot[loc]["strict_gate_pass"]})
    canonical_strict_pass = sum(r["status"] == "TESTED" and float(r["p"]) < config["strict_gate_p"]
                                for r in canonical.values())
    result = {
        "analysis_id": config["analysis_id"], "status": "PASS_TRAIT_ONLY_DIAGNOSTIC" if rule_pass else "NO_MATERIAL_PILOT_IMPROVEMENT",
        "scope": rule["decision_scope"], "pilot_loci": len(pilot),
        "canonical_tested": canon_tested, "pilot_tested": pilot_tested,
        "canonical_low_local_h2_not_run": canon_low_h2, "pilot_low_local_h2_not_run": pilot_low_h2,
        "absolute_tested_fraction_gain": gain, "relative_low_local_h2_not_run_reduction": reduction,
        "canonical_status_counts": dict(Counter(r["status"] for r in canonical.values())),
        "pilot_status_counts": dict(Counter(r["status"] for r in pilot.values())),
        "canonical_not_run_reasons": dict(Counter(r["reason"] for r in canonical.values() if r["status"] == "NOT_RUN")),
        "pilot_not_run_reasons": dict(Counter(r["reason"] for r in pilot.values() if r["status"] == "NOT_RUN")),
        "canonical_strict_gate_pass": canonical_strict_pass,
        "pilot_strict_gate_pass": sum(r["strict_gate_pass"].upper() == "TRUE" for r in pilot.values()),
        "advance_to_full_trait_only_screen": rule_pass,
        "family_promotion_authorized": False,
        "pilot_config_sha256": sha(PILOT), "decision_rule_sha256": sha(RULE),
        "canonical_aggregate_sha256": sha(CANONICAL),
        "materialization_receipt_sha256": sha(output / "materialization.receipt.json"),
        "aggregate_sha256": sha(result_path),
        "comparison_sha256": sha(comparison_path),
        "worker_outputs": {str(item["data"]["worker_id"]): {
            "config_sha256": item["sha256"],
            "tsv_sha256": sha(output / "results" / f"worker_{item['data']['worker_id']}.tsv"),
            "summary_sha256": sha(output / "results" / f"worker_{item['data']['worker_id']}.summary.json")}
            for item in worker_configs},
        "interpretation_limit": rule["interpretation_limit"],
    }
    save_new(output / "results" / "pilot_decision.json", result)
    return result


def main() -> None:
    config = json.loads(PILOT.read_text())
    rule = json.loads(RULE.read_text())
    if rule["analysis_id"] != config["analysis_id"] or rule["canonical_aggregate_sha256"] != config["canonical_aggregate_sha256"]:
        raise ValueError("Rule/config hash mismatch")
    if sha(RULE) != config["decision_rule_sha256"]:
        raise ValueError("Pinned decision rule changed")
    if rule["analysis_id"] != config["analysis_id"]:
        raise ValueError("Pilot decision rule identity mismatch")
    output = Path(config["output_root"])
    receipt, worker_configs = validate_inputs(config, output)
    if any((output / "results" / f"worker_{w}.tsv").exists() for w in range(1, 5)):
        raise FileExistsError("Some pilot worker results exist; inspect and resume explicitly")
    env = dict(os.environ)
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    handles = []
    children = []
    launch = {"analysis_id": config["analysis_id"], "launched_utc": utc(),
              "pilot_config_sha256": sha(PILOT), "rule_sha256": sha(RULE),
              "materialization_receipt_sha256": sha(output / "materialization.receipt.json"),
              "rscript": str(RSCRIPT), "rscript_sha256": sha(RSCRIPT),
              "r_runner_sha256": sha(RUNNER), "workers": []}
    for item in worker_configs:
        worker = item["data"]["worker_id"]
        log_path = output / "results" / f"worker_{worker}.log"
        log = log_path.open("xb")
        handles.append(log)
        command = [str(RSCRIPT), str(RUNNER), str(item["path"])]
        child = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        children.append(child)
        launch["workers"].append({"worker_id": worker, "pid": child.pid, "command": command,
                                  "config_sha256": item["sha256"], "log": str(log_path)})
    save_new(output / "results" / "launch.json", launch)
    print(json.dumps({"status": "PILOT_RUNNING", "workers": len(children),
                      "pids": [child.pid for child in children], "output": str(output)}, sort_keys=True), flush=True)
    codes = [child.wait() for child in children]
    for handle in handles:
        handle.close()
    exit_receipt = {"analysis_id": config["analysis_id"], "finished_utc": utc(),
                    "launch_sha256": sha(output / "results" / "launch.json"),
                    "exit_codes": codes,
                    "logs": {str(w): sha(output / "results" / f"worker_{w}.log") for w in range(1, 5)}}
    save_new(output / "results" / "exit.json", exit_receipt)
    if any(code != 0 for code in codes):
        raise RuntimeError(f"Pilot worker failure: {codes}")
    result = aggregate(config, rule, output, worker_configs)
    print(json.dumps({"status": result["status"], "canonical_tested": result["canonical_tested"],
                      "pilot_tested": result["pilot_tested"], "advance": result["advance_to_full_trait_only_screen"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
