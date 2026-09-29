#!/usr/bin/env python3
"""Freeze, execute, and adjudicate two predeclared long-sleep LAVA pilot arms."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "brain6/config/lava_confirmatory_pilots_v1/longsleep_linear_88_v1.json"
PROTOCOL = ROOT / "brain6/results/lava_confirmatory_protocol_v1/longsleep_linear_model_pilot_addendum_20260927.md"
SELECTION = ROOT / "brain6/results/lava_confirmatory_pilots_v1/longsleep_linear_88_loci.tsv"
CANONICAL = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
OUTPUT = Path("/Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1/longsleep_linear_88_v1")
MATERIALIZATION = OUTPUT / "materialization.receipt.json"
REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1")
REFERENCE_PROVENANCE = REFERENCE / "reference.provenance.json"
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
RUNNER = ROOT / "brain6/scripts/run_longsleep_linear_pilot_v1.R"
DRIVER = Path(__file__).resolve()
ANALYSIS = "brain6_longsleep_linear_88_pilot_v1"
ARMS = ("linear", "binary")
WORKERS = 4
STRICT_GATE = 2.86286859433152e-06


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def new_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def partitions(selection: list[dict[str, str]]) -> dict[int, list[str]]:
    result = {worker: [] for worker in range(1, WORKERS + 1)}
    for row in selection:
        worker = 1 + (int(row["CHR"]) - 1) % WORKERS
        result[worker].append(row["LOC"])
    if Counter(locus for group in result.values() for locus in group) != Counter(row["LOC"] for row in selection):
        raise ValueError("Worker partition does not cover the exact selected family")
    return result


def materialized() -> tuple[list[dict[str, str]], dict]:
    if sha(SELECTION) != "9182ebca331a7362d96359cc41e130571ce274f7869477864ce7ca0a944a01e7":
        raise ValueError("Selected loci changed")
    selection = read_tsv(SELECTION)
    receipt = read_json(MATERIALIZATION)
    if len(selection) != 88 or receipt["analysis_id"] != ANALYSIS:
        raise ValueError("Materialization identity changed")
    if receipt["loci_file_sha256"] != sha(OUTPUT / "selected.loci"):
        raise ValueError("Selected LAVA locus file changed")
    if receipt["source_sha256"][str(PROTOCOL)] != sha(PROTOCOL):
        raise ValueError("Pre-run protocol changed")
    if receipt["materializer_sha256"] != sha(ROOT / "brain6/scripts/materialize_longsleep_linear_pilot_v1.py"):
        raise ValueError("Materializer changed")
    if len(receipt["selected_loci"]) != 88 or {item["locus_id"] for item in receipt["selected_loci"]} != {row["LOC"] for row in selection}:
        raise ValueError("Materialization locus inventory changed")
    for item in receipt["selected_loci"]:
        locus_id = item["locus_id"]
        for arm in ARMS:
            directory = OUTPUT / "inputs" / arm / f"locus_{locus_id}"
            for name, key in (("longsleep.sumstats.tsv.gz", f"{arm}_shard_sha256"),
                              ("input_info.tsv", f"{arm}_info_sha256")):
                if sha(directory / name) != item[key]:
                    raise ValueError(f"Prepared input changed: {arm}/{locus_id}/{name}")
        if item["linear_shard_sha256"] != item["binary_shard_sha256"]:
            raise ValueError("Source SNP/Z/N differs between arms")
    return selection, receipt


def frozen() -> tuple[dict, list[dict[str, str]], dict]:
    config = read_json(CONFIG)
    selection, receipt = materialized()
    for path, key in ((PROTOCOL, "protocol_sha256"), (SELECTION, "selection_sha256"),
                      (MATERIALIZATION, "materialization_sha256"),
                      (REFERENCE_PROVENANCE, "reference_provenance_sha256"),
                      (RSCRIPT, "rscript_sha256"), (RUNNER, "runner_sha256"),
                      (DRIVER, "driver_sha256"), (CANONICAL, "canonical_aggregate_sha256")):
        if sha(path) != config[key]:
            raise ValueError(f"Frozen dependency changed: {path}")
    if config["analysis_id"] != ANALYSIS or config["worker_count"] != WORKERS or config["arms"] != list(ARMS):
        raise ValueError("Frozen pilot configuration changed")
    if config["locus_ids"] != [row["LOC"] for row in selection]:
        raise ValueError("Frozen pilot locus order changed")
    if config["partitions"] != {str(key): value for key, value in partitions(selection).items()}:
        raise ValueError("Frozen worker partition changed")
    if config["strict_gate_p"] != STRICT_GATE or config["source_N_proxy"] != 339926:
        raise ValueError("Frozen pilot statistics changed")
    return config, selection, receipt


def freeze() -> None:
    if CONFIG.exists():
        raise FileExistsError(CONFIG)
    selection, receipt = materialized()
    config = {
        "analysis_id": ANALYSIS,
        "scope": "PREDECLARED_TRAIT_ONLY_DIAGNOSTIC_NO_PROMOTION",
        "arms": list(ARMS), "worker_count": WORKERS,
        "locus_ids": [row["LOC"] for row in selection],
        "partitions": {str(key): value for key, value in partitions(selection).items()},
        "source_N_proxy": 339926,
        "primary_arm": "linear", "secondary_arm": "binary",
        "strict_gate_p": STRICT_GATE,
        "random_seed": 20260922,
        "execution_policy": {"locus_processing": {"min_K": 2, "prune_thresh": 99,
                           "max_prop_K": 0.75, "drop_failed": True,
                           "max_block_size": 3000, "cap_estimates": True},
                           "univariate": {"cap_estimates": True}},
        "no_promotion": True,
        "protocol_sha256": sha(PROTOCOL),
        "selection_sha256": sha(SELECTION),
        "materialization_sha256": sha(MATERIALIZATION),
        "reference_provenance_sha256": sha(REFERENCE_PROVENANCE),
        "rscript_sha256": sha(RSCRIPT), "runner_sha256": sha(RUNNER),
        "driver_sha256": sha(DRIVER), "canonical_aggregate_sha256": sha(CANONICAL),
        "output_root": str(OUTPUT),
        "retained_rows": {item["locus_id"]: item["retained_rows"] for item in receipt["selected_loci"]},
    }
    new_json(CONFIG, config)
    print(json.dumps({"status": "FROZEN_BEFORE_PILOT", "config_sha256": sha(CONFIG),
                      "loci": len(selection), "workers_per_arm": WORKERS}, sort_keys=True))


def run_arm(config: dict, arm: str) -> None:
    result_dir = OUTPUT / "results" / arm
    exit_path = result_dir / "exit.json"
    if exit_path.exists():
        receipt = read_json(exit_path)
        if receipt.get("exit_codes") != [0] * WORKERS:
            raise ValueError(f"Prior {arm} attempt did not exit cleanly; preserve it for audit")
        if receipt.get("launch_sha256") != sha(result_dir / "launch.json"):
            raise ValueError(f"Prior {arm} launch receipt changed")
        print(json.dumps({"status": "ARM_ALREADY_COMPLETE", "arm": arm}, sort_keys=True), flush=True)
        return
    if result_dir.exists():
        raise ValueError(f"{arm} output directory exists without terminal receipt; inspect live processes first")
    result_dir.mkdir(parents=True, exist_ok=False)
    input_rows = config["retained_rows"]
    env = os.environ.copy()
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        env[key] = "1"
    workers = []
    handles = []
    launch = {"analysis_id": ANALYSIS, "arm_id": arm, "launched_utc": utc(),
              "config_sha256": sha(CONFIG), "runner_sha256": sha(RUNNER),
              "rscript_sha256": sha(RSCRIPT), "workers": []}
    try:
        for worker in range(1, WORKERS + 1):
            ids = config["partitions"][str(worker)]
            worker_cfg = {"analysis_id": ANALYSIS, "arm_id": arm,
                          "worker_id": worker, "trait_id": "longsleep",
                          "scope": config["scope"], "expected_loci": len(ids),
                          "locus_ids": ids, "input_rows": {key: input_rows[key] for key in ids},
                          "input_root": str(OUTPUT / "inputs" / arm),
                          "loci_file": str(OUTPUT / "selected.loci"),
                          "reference_prefix_stem": str(REFERENCE / "lava-ukb-v1.1_chr"),
                          "strict_gate_p": STRICT_GATE,
                          "random_seed": config["random_seed"] + worker,
                          "execution_policy": config["execution_policy"],
                          "output_tsv": str(result_dir / f"worker_{worker}.tsv"),
                          "summary_json": str(result_dir / f"worker_{worker}.summary.json")}
            cfg_path = result_dir / f"worker_{worker}.config.json"
            new_json(cfg_path, worker_cfg)
            log_path = result_dir / f"worker_{worker}.log"
            log = log_path.open("xb")
            handles.append(log)
            command = [str(RSCRIPT), str(RUNNER), str(cfg_path)]
            child = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log,
                                     stderr=subprocess.STDOUT)
            workers.append(child)
            launch["workers"].append({"worker_id": worker, "pid": child.pid,
                                       "config_sha256": sha(cfg_path),
                                       "log_path": str(log_path), "command": command})
        new_json(result_dir / "launch.json", launch)
        print(json.dumps({"status": "ARM_RUNNING", "arm": arm, "workers": WORKERS,
                          "pids": [item.pid for item in workers]}, sort_keys=True), flush=True)
        codes = [child.wait() for child in workers]
    finally:
        for handle in handles:
            handle.close()
    exit_receipt = {"analysis_id": ANALYSIS, "arm_id": arm,
                    "finished_utc": utc(), "launch_sha256": sha(result_dir / "launch.json"),
                    "exit_codes": codes,
                    "logs_sha256": {str(worker): sha(result_dir / f"worker_{worker}.log")
                                    for worker in range(1, WORKERS + 1)}}
    new_json(exit_path, exit_receipt)
    print(json.dumps({"status": "ARM_EXITED", "arm": arm, "exit_codes": codes}, sort_keys=True), flush=True)
    if codes != [0] * WORKERS:
        raise ValueError(f"{arm} workers did not all exit cleanly")


def run() -> None:
    config, _, _ = frozen()
    for arm in ARMS:
        run_arm(config, arm)


def adjudicate() -> None:
    config, selection, materialization = frozen()
    by_old = {row["locus_id"]: row for row in read_tsv(CANONICAL) if row["phen"] == "longsleep"}
    selected_ids = config["locus_ids"]
    arm_rows: dict[str, dict[str, dict[str, str]]] = {}
    output_hashes = {}
    for arm in ARMS:
        directory = OUTPUT / "results" / arm
        terminal = read_json(directory / "exit.json")
        launch = read_json(directory / "launch.json")
        if terminal["exit_codes"] != [0] * WORKERS or terminal["launch_sha256"] != sha(directory / "launch.json"):
            raise ValueError(f"{arm} lacks a valid terminal exit")
        if launch["config_sha256"] != sha(CONFIG):
            raise ValueError(f"{arm} launched under another configuration")
        combined = {}
        worker_hashes = {}
        for worker in range(1, WORKERS + 1):
            cfg_path = directory / f"worker_{worker}.config.json"
            table = directory / f"worker_{worker}.tsv"
            summary = directory / f"worker_{worker}.summary.json"
            log = directory / f"worker_{worker}.log"
            cfg = read_json(cfg_path)
            rows = read_tsv(table)
            counts = read_json(summary)
            if cfg["locus_ids"] != config["partitions"][str(worker)]:
                raise ValueError(f"{arm}/{worker} partition changed")
            if set(row["locus_id"] for row in rows) != set(cfg["locus_ids"]) or len(rows) != len(cfg["locus_ids"]):
                raise ValueError(f"{arm}/{worker} output incomplete")
            if counts["rows"] != len(rows) or counts["failed"] != 0:
                raise ValueError(f"{arm}/{worker} summary failed")
            if terminal["logs_sha256"][str(worker)] != sha(log):
                raise ValueError(f"{arm}/{worker} log changed")
            if launch["workers"][worker - 1]["config_sha256"] != sha(cfg_path):
                raise ValueError(f"{arm}/{worker} worker config changed")
            for row in rows:
                if row["locus_id"] in combined or row["status"] not in {"TESTED", "NOT_RUN"}:
                    raise ValueError(f"{arm} duplicate or invalid locus result")
                combined[row["locus_id"]] = row
            worker_hashes[str(worker)] = {"config_sha256": sha(cfg_path),
                                          "tsv_sha256": sha(table),
                                          "summary_sha256": sha(summary),
                                          "log_sha256": sha(log)}
        if set(combined) != set(selected_ids):
            raise ValueError(f"{arm} does not cover the selected family")
        arm_rows[arm] = combined
        output_hashes[arm] = {"launch_sha256": sha(directory / "launch.json"),
                              "exit_sha256": sha(directory / "exit.json"),
                              "workers": worker_hashes}

    paired = OUTPUT / "results" / "paired_comparison.tsv"
    columns = ("locus_id", "chromosome", "canonical_status", "linear_status", "binary_status",
               "linear_reason", "binary_reason", "linear_n_snps", "binary_n_snps",
               "linear_h2_obs", "binary_h2_obs", "linear_p", "binary_p",
               "linear_strict_gate", "binary_strict_gate")
    with paired.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for locus in selection:
            locus_id = locus["LOC"]
            linear = arm_rows["linear"][locus_id]
            binary = arm_rows["binary"][locus_id]
            writer.writerow({"locus_id": locus_id, "chromosome": locus["CHR"],
                             "canonical_status": by_old[locus_id]["status"],
                             "linear_status": linear["status"], "binary_status": binary["status"],
                             "linear_reason": linear["reason"], "binary_reason": binary["reason"],
                             "linear_n_snps": linear["n_snps"], "binary_n_snps": binary["n_snps"],
                             "linear_h2_obs": linear["h2_obs"], "binary_h2_obs": binary["h2_obs"],
                             "linear_p": linear["p"], "binary_p": binary["p"],
                             "linear_strict_gate": linear["strict_gate_pass"],
                             "binary_strict_gate": binary["strict_gate_pass"]})
    coverage = [item["retained_rows"] / item["original_rows"] if item["original_rows"] else 0.0
                for item in materialization["selected_loci"]]
    coverage_good = sum(value >= 0.8 for value in coverage)
    old_tested = sum(by_old[locus]["status"] == "TESTED" for locus in selected_ids)
    counts = {arm: dict(Counter(row["status"] for row in arm_rows[arm].values())) for arm in ARMS}
    primary_gain = counts["linear"].get("TESTED", 0) - old_tested
    decision = {
        "analysis_id": ANALYSIS, "scope": config["scope"],
        "status": "TRAIT_ONLY_FEASIBILITY_SIGNAL" if primary_gain / 88 >= 0.10 and coverage_good / 88 >= 0.90 else "NO_FULL_SCREEN_JUSTIFICATION",
        "full_family_promotion_authorized": False,
        "full_trait_screen_authorized": False,
        "primary_model": "linear_BOLT_LMM", "secondary_model": "binary_reconstruction",
        "selected_loci": len(selected_ids), "canonical_tested": old_tested,
        "arm_counts": counts,
        "primary_tested_gain_loci": primary_gain,
        "primary_tested_gain_percentage_points": 100 * primary_gain / 88,
        "coverage_loci_retaining_at_least_80_percent": coverage_good,
        "coverage_fraction": coverage_good / 88,
        "linear_binary_status_discordant_loci": [locus for locus in selected_ids
                                                 if arm_rows["linear"][locus]["status"] != arm_rows["binary"][locus]["status"]],
        "linear_binary_strict_gate_discordant_loci": [locus for locus in selected_ids
            if arm_rows["linear"][locus]["strict_gate_pass"] != arm_rows["binary"][locus]["strict_gate_pass"]],
        "config_sha256": sha(CONFIG), "protocol_sha256": sha(PROTOCOL),
        "materialization_sha256": sha(MATERIALIZATION),
        "canonical_aggregate_sha256": sha(CANONICAL),
        "paired_comparison_sha256": sha(paired), "worker_outputs": output_hashes,
        "interpretation": "The prespecified primary linear arm is diagnostic. Neither arm can be selected by testability or association result; source N and mixed-model calibration plus pair overlap remain unresolved for admission.",
    }
    path = OUTPUT / "results" / "pilot_decision.json"
    new_json(path, decision)
    print(json.dumps({"status": decision["status"], "arm_counts": counts,
                      "primary_gain": primary_gain, "coverage_good": coverage_good,
                      "decision_sha256": sha(path)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "run", "adjudicate"))
    action = parser.parse_args().action
    {"freeze": freeze, "run": run, "adjudicate": adjudicate}[action]()
