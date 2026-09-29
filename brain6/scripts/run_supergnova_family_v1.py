#!/usr/bin/env python3
"""Resume a frozen five-pair SUPERGNOVA family without rewriting receipts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import pathlib
import subprocess
import sys


REPO = pathlib.Path(__file__).resolve().parents[2]
PROTOCOL = REPO / "brain6/results/brain6_alternative_local_validation_v1/protocol_freeze.json"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checked(path: str, expected: str) -> pathlib.Path:
    p = pathlib.Path(path)
    if not p.is_file() or sha256(p) != expected:
        raise RuntimeError(f"Frozen hash mismatch or missing file: {path}")
    return p


def load_protocol() -> dict:
    p = json.loads(PROTOCOL.read_text())
    for item in p["frozen_files"]:
        checked(item["path"], item["sha256"])
    return p


def atomic_json(path: pathlib.Path, object_: dict) -> None:
    if path.exists():
        raise FileExistsError(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(object_, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def prepare_inputs(protocol: dict) -> None:
    root = pathlib.Path(protocol["ssd_root"])
    output_dir = root / "sumstats_v1"
    output_dir.mkdir(parents=True, exist_ok=True)
    ref_bim = pathlib.Path(protocol["adapter_reference_bim"])
    adapter = REPO / "brain6/scripts/prepare_supergnova_sumstats_v1.py"
    for trait, card in protocol["traits"].items():
        output = output_dir / f"{trait}.sumstats.gz"
        receipt = output_dir / f"{trait}.adapter_receipt.json"
        if receipt.exists():
            info = json.loads(receipt.read_text())
            if (not output.is_file() or sha256(output) != info["output_sha256"] or
                    info["input_sha256"] != card["input_sha256"] or
                    info["reference_bim_sha256"] != protocol["reference_bim_sha256"] or
                    info["study_total_n_used_by_supergnova"] != card["n_total"]):
                raise RuntimeError(f"Existing adapter receipt failed: {trait}")
            print(f"SKIP_VERIFIED_ADAPTER {trait}", flush=True)
            continue
        if output.exists():
            raise RuntimeError(f"Unreceipted adapter output exists: {output}")
        command = [sys.executable, str(adapter), "--reference-bim", str(ref_bim),
                   "--source", card["input"], "--output", str(output),
                   "--trait", trait, "--n-total", str(card["n_total"])]
        subprocess.run(command, check=True)


def expected_blocks(protocol: dict) -> set[tuple[int, int, int]]:
    with pathlib.Path(protocol["partition_file"]).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {(int(row["chr"]), int(row["start"]), int(row["end"])) for row in reader}
    if len(rows) != protocol["n_blocks_nonmhc"]:
        raise RuntimeError("Frozen partition cardinality changed")
    return rows


def audit_raw_output(path: pathlib.Path, protocol: dict) -> dict:
    expected = expected_blocks(protocol)
    seen = set()
    counts = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter=" ")
        required = {"chr", "start", "end", "rho", "corr", "h2_1", "h2_2", "var", "p", "m", "analysis_status"}
        if not reader.fieldnames or set(reader.fieldnames) != required:
            raise RuntimeError(f"Unexpected SUPERGNOVA output columns: {reader.fieldnames}")
        for row in reader:
            key = (int(row["chr"]), int(row["start"]), int(row["end"]))
            if key not in expected or key in seen:
                raise RuntimeError(f"Unexpected/duplicate SUPERGNOVA block: {key}")
            seen.add(key)
            status = row["analysis_status"]
            if status not in {"ESTIMATED", "NOT_TESTED_LOW_SNP", "NOT_TESTED_LOW_RANK"}:
                raise RuntimeError(f"Invalid method block status: {status}")
            m = int(row["m"])
            if status == "NOT_TESTED_LOW_SNP" and m >= 120:
                raise RuntimeError(f"Low-SNP status with m={m}: {key}")
            if status == "ESTIMATED":
                p = float(row["p"])
                rho = float(row["rho"])
                variance = float(row["var"])
                if not (math.isfinite(p) and 0 <= p <= 1 and math.isfinite(rho)
                        and math.isfinite(variance) and variance > 0 and m >= 120):
                    raise RuntimeError(f"Nonfinite/invalid estimate: {key}")
            else:
                if row["p"] not in {"NA", "", "nan"}:
                    raise RuntimeError(f"Untested block has p value: {key}")
            counts[status] = counts.get(status, 0) + 1
    if seen != expected:
        raise RuntimeError(f"Missing method blocks: {len(expected - seen)}")
    return counts


def run_pair(protocol: dict, pair: str) -> None:
    if pair not in protocol["pairs"]:
        raise ValueError(pair)
    root = pathlib.Path(protocol["ssd_root"])
    sumstats = root / "sumstats_v1"
    outdir = root / "runs_v1"
    outdir.mkdir(parents=True, exist_ok=True)
    a, b = protocol["pairs"][pair]
    for trait in (a, b):
        adapter = json.loads((sumstats / f"{trait}.adapter_receipt.json").read_text())
        checked(str(sumstats / f"{trait}.sumstats.gz"), adapter["output_sha256"])
    output = outdir / f"{pair}.supergnova.txt"
    log = outdir / f"{pair}.log"
    receipt = outdir / f"{pair}.receipt.json"
    if receipt.exists():
        old = json.loads(receipt.read_text())
        checked(str(output), old["output_sha256"])
        checked(str(log), old["log_sha256"])
        if old["protocol_sha256"] != sha256(PROTOCOL):
            raise RuntimeError(f"Receipt belongs to a different protocol: {pair}")
        audit_raw_output(output, protocol)
        print(f"SKIP_VERIFIED_PAIR {pair}", flush=True)
        return
    if output.exists() or log.exists():
        raise RuntimeError(f"Unreceipted run artifacts exist for {pair}; inspect before retry")
    software = root / "SUPERGNOVA_source/supergnova.py"
    command = [protocol["python_executable"], str(software),
               str(sumstats / f"{a}.sumstats.gz"), str(sumstats / f"{b}.sumstats.gz"),
               "--N1", str(protocol["traits"][a]["n_total"]),
               "--N2", str(protocol["traits"][b]["n_total"]),
               "--bfile", protocol["reference_prefix"],
               "--partition", protocol["partition_file"],
               "--thread", str(protocol["workers"]), "--out", str(output)]
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    with log.open("w") as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
    if completed.returncode:
        raise RuntimeError(f"SUPERGNOVA failed for {pair}: exit={completed.returncode}, log={log}")
    counts = audit_raw_output(output, protocol)
    atomic_json(receipt, {
        "pair_id": pair, "protocol_sha256": sha256(PROTOCOL),
        "command": command, "workers": protocol["workers"],
        "trait_n_total": {a: protocol["traits"][a]["n_total"],
                          b: protocol["traits"][b]["n_total"]},
        "output_sha256": sha256(output), "log_sha256": sha256(log),
        "block_status_counts": counts, "n_blocks": sum(counts.values()),
    })
    print(f"COMPLETE_VERIFIED_PAIR {pair} {counts}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["prepare", "run"], required=True)
    parser.add_argument("--pair", help="Run only this frozen pair; omitted means all pairs")
    args = parser.parse_args()
    protocol = load_protocol()
    if args.stage == "prepare":
        prepare_inputs(protocol)
    else:
        pairs = [args.pair] if args.pair else list(protocol["pairs"])
        for pair in pairs:
            run_pair(protocol, pair)
