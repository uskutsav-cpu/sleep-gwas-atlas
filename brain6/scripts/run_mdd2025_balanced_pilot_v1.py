#!/usr/bin/env python3
"""Freeze, run, and receipt-audit a separate balanced-equivalent MDD pilot arm."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OLD_CONFIG = ROOT / "brain6/config/lava_multitrait_feasibility_v1/mdd2025_trait_pilot_v1.json"
CONFIG = ROOT / "brain6/config/lava_multitrait_feasibility_v1/mdd2025_balanced_equivalent_pilot_v1.json"
PROTOCOL = ROOT / "brain6/results/lava_confirmatory_protocol_v1/mdd2025_balanced_equivalent_pilot_addendum_20260927.md"
RUNNER = ROOT / "brain6/scripts/run_mdd2025_balanced_pilot_v1.R"
DRIVER = Path(__file__).resolve()
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
ANALYSIS = "brain6_mdd2025_balanced_equivalent_pilot_v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def new_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def frozen() -> tuple[dict, dict, Path, Path]:
    cfg = json.loads(CONFIG.read_text())
    old = json.loads(OLD_CONFIG.read_text())
    if cfg["analysis_id"] != ANALYSIS or cfg["old_config_sha256"] != sha(OLD_CONFIG):
        raise ValueError("Balanced pilot or original frozen config changed")
    for path, key in ((PROTOCOL, "protocol_sha256"), (RUNNER, "runner_sha256"),
                      (DRIVER, "driver_sha256"),
                      (RSCRIPT, "rscript_sha256")):
        if cfg[key] != sha(path):
            raise ValueError(f"Frozen dependency changed: {path}")
    for key in ("locus_ids", "expected_loci", "worker_count", "seed", "strict_gate_p",
                "locus_processing", "univariate", "loci_file_sha256",
                "reference_provenance_sha256", "source_sha256"):
        if cfg[key] != old[key]:
            raise ValueError(f"Original pilot setting changed: {key}")
    old_root, root = Path(old["output_root"]), Path(cfg["output_root"])
    if not old_root.is_dir() or root == old_root:
        raise ValueError("Original pilot unavailable or output roots overlap")
    receipt = old_root / "materialization.receipt.json"
    if cfg["old_materialization_sha256"] != sha(receipt):
        raise ValueError("Original materialization receipt changed")
    return cfg, old, root, old_root


def freeze() -> None:
    if CONFIG.exists():
        raise FileExistsError(CONFIG)
    old = json.loads(OLD_CONFIG.read_text())
    old_root = Path(old["output_root"])
    config = {key: old[key] for key in ("locus_ids", "expected_loci", "worker_count", "seed",
        "strict_gate_p", "locus_processing", "univariate", "loci_file_sha256",
        "reference_provenance_sha256", "source_sha256")}
    config.update({
        "schema_version": 1, "analysis_id": ANALYSIS,
        "scope": "PREDECLARED_TRAIT_ONLY_SOURCE_MODEL_SENSITIVITY_NO_PROMOTION",
        "sample_size_model": "source_per_variant_NEFF_with_balanced_case_fraction_0.5",
        "encoded_input_info_cases": 1, "encoded_input_info_controls": 1,
        "old_config_sha256": sha(OLD_CONFIG),
        "old_materialization_sha256": sha(old_root / "materialization.receipt.json"),
        "protocol_sha256": sha(PROTOCOL), "runner_sha256": sha(RUNNER),
        "driver_sha256": sha(DRIVER),
        "rscript_sha256": sha(RSCRIPT),
        "output_root": str(old_root.with_name("mdd2025_balanced_equivalent_pilot_v1")),
    })
    new_json(CONFIG, config)
    print(json.dumps({"config_sha256": sha(CONFIG), "loci": len(config["locus_ids"]),
                      "workers": config["worker_count"]}, sort_keys=True))


def verified_inputs(cfg: dict, root: Path, old_root: Path) -> dict:
    receipt = json.loads((old_root / "materialization.receipt.json").read_text())
    if set(receipt["records"]) != set(cfg["locus_ids"]):
        raise ValueError("Original materialization omitted loci")
    for locus in cfg["locus_ids"]:
        old_dir = old_root / "inputs" / f"locus_{locus}"
        record = receipt["records"][locus]
        for name, key in (("mdd2025_no23andme.sumstats.tsv.gz", "sumstats_sha256"),
                          ("input_info.tsv", "input_info_sha256")):
            if sha(old_dir / name) != record[key]:
                raise ValueError(f"Original input hash mismatch: {locus}/{name}")
    if root.exists():
        new_receipt = json.loads((root / "materialization.receipt.json").read_text())
        if new_receipt["config_sha256"] != sha(CONFIG) or set(new_receipt["records"]) != set(cfg["locus_ids"]):
            raise ValueError("Balanced input receipt mismatch")
        for locus, record in new_receipt["records"].items():
            new_dir = root / "inputs" / f"locus_{locus}"
            if sha(new_dir / "mdd2025_no23andme.sumstats.tsv.gz") != record["sumstats_sha256"]:
                raise ValueError(f"Balanced sumstats hash mismatch: {locus}")
            if sha(new_dir / "input_info.tsv") != record["input_info_sha256"]:
                raise ValueError(f"Balanced case-fraction hash mismatch: {locus}")
            expected = "phenotype\tcases\tcontrols\tfilename\n" + (
                "mdd2025_no23andme\t1\t1\tmdd2025_no23andme.sumstats.tsv.gz\n")
            if (new_dir / "input_info.tsv").read_text() != expected:
                raise ValueError(f"Balanced case-fraction content mismatch: {locus}")
            target = new_dir / "mdd2025_no23andme.sumstats.tsv.gz"
            if not target.is_symlink() or target.resolve() != (old_root / "inputs" / f"locus_{locus}" / target.name).resolve():
                raise ValueError(f"Sumstats symlink target mismatch: {locus}")
        return new_receipt
    return receipt


def prepare() -> None:
    cfg, old, root, old_root = frozen()
    old_receipt = verified_inputs(cfg, root, old_root)
    if root.exists():
        raise FileExistsError(f"Balanced pilot already prepared: {root}")
    root.mkdir()
    (root / "inputs").mkdir()
    (root / "results").mkdir()
    records = {}
    for locus in cfg["locus_ids"]:
        old_dir = old_root / "inputs" / f"locus_{locus}"
        new_dir = root / "inputs" / f"locus_{locus}"
        new_dir.mkdir()
        source = old_dir / "mdd2025_no23andme.sumstats.tsv.gz"
        target = new_dir / source.name
        target.symlink_to(source)
        info = new_dir / "input_info.tsv"
        info.write_text("phenotype\tcases\tcontrols\tfilename\n"
            "mdd2025_no23andme\t1\t1\tmdd2025_no23andme.sumstats.tsv.gz\n")
        records[locus] = {
            "source_sumstats_sha256": old_receipt["records"][locus]["sumstats_sha256"],
            "sumstats_sha256": sha(target), "input_info_sha256": sha(info),
        }
    receipt = {
        "analysis_id": ANALYSIS, "status": "BALANCED_EQUIVALENT_INPUTS_PREPARED",
        "config_sha256": sha(CONFIG), "old_materialization_sha256": sha(old_root / "materialization.receipt.json"),
        "records": records,
    }
    new_json(root / "materialization.receipt.json", receipt)
    found = []
    for worker in range(1, cfg["worker_count"] + 1):
        old_worker = json.loads((old_root / f"worker_{worker}.config.json").read_text())
        loci = cfg["locus_ids"][worker - 1::cfg["worker_count"]]
        if old_worker["locus_ids"] != loci:
            raise ValueError(f"Original worker partition mismatch: {worker}")
        found.extend(loci)
        worker_cfg = dict(old_worker)
        worker_cfg.update({
            "analysis_id": ANALYSIS, "scope": cfg["scope"], "input_root": str(root / "inputs"),
            "output_tsv": str(root / "results" / f"worker_{worker}.tsv"),
            "summary_json": str(root / "results" / f"worker_{worker}.summary.json"),
            "pilot_config_sha256": sha(CONFIG),
            "materialization_receipt_sha256": sha(root / "materialization.receipt.json"),
        })
        new_json(root / f"worker_{worker}.config.json", worker_cfg)
    if len(found) != len(set(found)) or set(found) != set(cfg["locus_ids"]):
        raise ValueError("Worker partitions do not cover exact pilot")
    verified_inputs(cfg, root, old_root)
    print(json.dumps({"status": receipt["status"], "loci": len(records), "output": str(root)}, sort_keys=True))


def run() -> None:
    cfg, _, root, old_root = frozen()
    verified_inputs(cfg, root, old_root)
    if (root / "results" / "launch.json").exists():
        raise FileExistsError("Launch exists; audit before any restart")
    workers = []
    for worker in range(1, cfg["worker_count"] + 1):
        path = root / f"worker_{worker}.config.json"
        data = json.loads(path.read_text())
        if data["analysis_id"] != ANALYSIS or data["pilot_config_sha256"] != sha(CONFIG):
            raise ValueError(f"Worker config mismatch: {worker}")
        if (root / "results" / f"worker_{worker}.tsv").exists():
            raise FileExistsError(f"Worker output already exists: {worker}")
        workers.append((worker, path, sha(path)))
    env = dict(os.environ)
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    children = []
    logs = []
    launch = {"analysis_id": ANALYSIS, "launched_utc": utc(), "config_sha256": sha(CONFIG),
              "runner_sha256": sha(RUNNER), "rscript_sha256": sha(RSCRIPT), "workers": []}
    try:
        for worker, path, config_hash in workers:
            log_path = root / "results" / f"worker_{worker}.log"
            log = log_path.open("xb")
            logs.append(log)
            command = [str(RSCRIPT), str(RUNNER), str(path)]
            child = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            children.append(child)
            launch["workers"].append({"worker_id": worker, "pid": child.pid,
                "config_sha256": config_hash, "command": command, "log": str(log_path)})
        new_json(root / "results" / "launch.json", launch)
        print(json.dumps({"status": "RUNNING", "workers": len(children), "pids": [c.pid for c in children]},
                         sort_keys=True), flush=True)
        codes = [child.wait() for child in children]
    finally:
        for log in logs:
            log.close()
    exit_receipt = {"analysis_id": ANALYSIS, "finished_utc": utc(),
        "launch_sha256": sha(root / "results" / "launch.json"),
        "exit_codes": codes, "logs_sha256": {str(w): sha(root / "results" / f"worker_{w}.log")
                                           for w in range(1, cfg["worker_count"] + 1)}}
    new_json(root / "results" / "exit.json", exit_receipt)
    print(json.dumps({"status": "WORKERS_EXITED", "exit_codes": codes}, sort_keys=True))


def adjudicate() -> None:
    cfg, _, root, old_root = frozen()
    verified_inputs(cfg, root, old_root)
    exit_receipt = json.loads((root / "results" / "exit.json").read_text())
    if exit_receipt["launch_sha256"] != sha(root / "results" / "launch.json") or exit_receipt["exit_codes"] != [0] * cfg["worker_count"]:
        raise ValueError("Balanced pilot workers did not all exit successfully")
    prior = {r["locus_id"]: r for r in rows(old_root / "results" / "pilot_aggregate.tsv")}
    if set(prior) != set(cfg["locus_ids"]):
        raise ValueError("Original pilot aggregate locus mismatch")
    combined = {}
    output_hashes = {}
    for worker in range(1, cfg["worker_count"] + 1):
        worker_cfg = json.loads((root / f"worker_{worker}.config.json").read_text())
        table = root / "results" / f"worker_{worker}.tsv"
        summary = root / "results" / f"worker_{worker}.summary.json"
        group = rows(table)
        counts = json.loads(summary.read_text())
        if counts["rows"] != len(group) or set(r["locus_id"] for r in group) != set(worker_cfg["locus_ids"]):
            raise ValueError(f"Worker result mismatch: {worker}")
        for row in group:
            if row["trait_id"] != "mdd2025_no23andme" or row["status"] not in {"TESTED", "NOT_RUN"}:
                raise ValueError(f"Invalid worker result: {worker}/{row['locus_id']}")
            if row["locus_id"] in combined:
                raise ValueError(f"Duplicate locus: {row['locus_id']}")
            combined[row["locus_id"]] = row
        output_hashes[str(worker)] = {"tsv_sha256": sha(table), "summary_sha256": sha(summary),
                                      "config_sha256": sha(root / f"worker_{worker}.config.json")}
    if set(combined) != set(cfg["locus_ids"]):
        raise ValueError("Balanced pilot incomplete")
    comparison = root / "results" / "balanced_vs_original.tsv"
    fields = ("locus_id", "original_status", "balanced_status", "original_reason", "balanced_reason",
              "original_p", "balanced_p", "original_strict_gate_pass", "balanced_strict_gate_pass")
    with comparison.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for locus in cfg["locus_ids"]:
            a, b = prior[locus], combined[locus]
            writer.writerow({"locus_id": locus, "original_status": a["status"],
                "balanced_status": b["status"], "original_reason": a["reason"],
                "balanced_reason": b["reason"], "original_p": a["p"], "balanced_p": b["p"],
                "original_strict_gate_pass": a["strict_gate_pass"],
                "balanced_strict_gate_pass": b["strict_gate_pass"]})
    discordant = [locus for locus in cfg["locus_ids"] if prior[locus]["status"] != combined[locus]["status"]]
    report = {"analysis_id": ANALYSIS, "scope": cfg["scope"], "family_promotion_authorized": False,
        "advance_to_full_family_authorized": False, "n_loci": len(combined),
        "original_status_counts": dict(Counter(r["status"] for r in prior.values())),
        "balanced_status_counts": dict(Counter(r["status"] for r in combined.values())),
        "status_discordant_loci": discordant,
        "strict_gate_discordant_loci": [locus for locus in cfg["locus_ids"]
            if prior[locus]["strict_gate_pass"] != combined[locus]["strict_gate_pass"]],
        "config_sha256": sha(CONFIG), "protocol_sha256": sha(PROTOCOL),
        "old_aggregate_sha256": sha(old_root / "results" / "pilot_aggregate.tsv"),
        "old_decision_sha256": sha(old_root / "results" / "pilot_decision.json"),
        "materialization_receipt_sha256": sha(root / "materialization.receipt.json"),
        "exit_receipt_sha256": sha(root / "results" / "exit.json"),
        "comparison_sha256": sha(comparison), "worker_outputs": output_hashes}
    new_json(root / "results" / "model_sensitivity_decision.json", report)
    print(json.dumps({"status": "TRAIT_ONLY_SENSITIVITY_AUDITED", "original": report["original_status_counts"],
        "balanced": report["balanced_status_counts"], "discordant": len(discordant)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "prepare", "run", "adjudicate"))
    action = parser.parse_args().action
    {"freeze": freeze, "prepare": prepare, "run": run, "adjudicate": adjudicate}[action]()
