#!/usr/bin/env python3
"""Resume an isolated public-input insomnia--ADHD PLACO+ sensitivity run.

This deliberately does not use, amend, or publish to the frozen Brain6 v3
family or the protected legacy Track B result directory.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
from brain6.artifacts import verify_artifact, verify_current_artifact  # noqa: E402
from brain6.executor import run_job  # noqa: E402
from brain6.io import check_hash, file_record, read_json, read_tsv, require, sha256, write_json, write_tsv  # noqa: E402
from brain6.pairs import verify_lock  # noqa: E402
from brain6.pipeline import create_job  # noqa: E402
from brain6.pleiotropy import collate_placo, split_pair  # noqa: E402


OUT = Path("/Volumes/Extreme SSD/brain6-work/placo-track-b-prospective-v1")
PAIR = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1/pair_insomnia__adhd")
FAMILY = ROOT / "brain6/config/placo_family_v3/family_lock.json"
PAIR_LOCK = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
SOURCE = Path("/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R")
RUNTIME = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
CHUNK_ROWS = 20_000
EXPECTED_ROWS = 5_514_402
FIELDS = ["path", "sha256"]


def contract() -> dict:
    family = read_json(FAMILY)
    pair_lock = verify_lock(PAIR_LOCK)
    require(family["protected_legacy_pair"] == "insomnia__adhd", "protected pair changed")
    require(family["n_selected_tracks"] == pair_lock["n_pairs"] == 5, "family changed")
    require(family["pair_selection_lock_sha256"] == sha256(PAIR_LOCK), "pair selection changed")
    require(family["placo_source_sha256"] == sha256(SOURCE), "PLACO source changed")
    require(family["runtime_sha256"] == sha256(RUNTIME), "R runtime changed")
    require(shutil.which("Rscript") == str(RUNTIME), "PATH does not select the frozen R runtime")
    for record in family["extension_code_inventory"]:
        check_hash(ROOT / record["path"], record["sha256"])
    receipt = verify_artifact(PAIR)
    require(receipt["scientific_status"] == "PASS" and not receipt["synthetic"], "pair input failed QC")
    require(read_json(PAIR / "qc.json")["retained_rows"] == EXPECTED_ROWS, "pair row count changed")
    expected_pair_hash = next(row["sha256"] for row in receipt["outputs"] if row["path"] == "pair.tsv.gz")
    require(sha256(PAIR / "pair.tsv.gz") == expected_pair_hash, "pair input hash changed")
    return {
        "schema": "brain6-track-b-prospective-public-placo.1",
        "analysis_id": "insomnia__adhd_public_input_prospective_sensitivity_v1",
        "claim_scope": "SENSITIVITY_ONLY_NOT_FROZEN_V3_OR_PROTECTED_LEGACY_TRACK_B",
        "pair_id": "insomnia__adhd",
        "pair": file_record(PAIR / "pair.tsv.gz"),
        "pair_receipt": file_record(PAIR / "receipt.json"),
        "family_lock": file_record(FAMILY),
        "pair_selection_lock": file_record(PAIR_LOCK),
        "placo_source": file_record(SOURCE),
        "r_runtime": file_record(RUNTIME),
        "runner": file_record(Path(__file__)),
        "method": "PLACO_PLUS",
        "parameters": {
            "extreme_z2": 80,
            "p_threshold": 1e-4,
            "min_variants": 1_000_000,
            "seed": 20260908,
            "absolute_tolerance": 1e-13,
            "max_numerical_failure_rate": 0.001,
            "chunk_rows": CHUNK_ROWS,
            "expected_rows": EXPECTED_ROWS,
            "within_pair_adjustment": "BH_all_eligible_rows",
            "across_pair_adjustment": "Bonferroni_5_tracks",
            "headline_p": 1e-8,
        },
    }


def freeze() -> None:
    proposed = contract()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "contract.json"
    if path.exists():
        require(read_json(path) == proposed, "prospective run contract changed")
    else:
        write_json(path, proposed)


def prepare() -> None:
    freeze()
    chunks = OUT / "chunks"
    if chunks.exists():
        verify_current_artifact(chunks)
    else:
        split_pair(PAIR, OUT, "chunks", chunk_rows=CHUNK_ROWS)
    manifest = read_json(chunks / "chunks.json")
    require(manifest["total_rows"] == EXPECTED_ROWS, "split input is incomplete")
    require(len(manifest["chunks"]) == (EXPECTED_ROWS + CHUNK_ROWS - 1) // CHUNK_ROWS, "chunk count changed")
    print(json.dumps({"stage": "PREPARED", "rows": EXPECTED_ROWS, "chunks": len(manifest["chunks"])}), flush=True)


def settings(**overrides: object) -> dict:
    return {
        "method": "PLACO_PLUS", "placo_source": str(SOURCE),
        "extreme_z2": 80, "p_threshold": 1e-4,
        "absolute_tolerance": 1e-13, "seed": 20260908,
        **overrides,
    }


def job(name: str, cfg: dict, inputs: dict[str, Path], outputs: dict[str, str]) -> tuple[Path, Path]:
    specs = OUT / "specs"
    specs.mkdir(exist_ok=True)
    settings_path = specs / f"{name}.settings.json"
    if settings_path.exists():
        require(read_json(settings_path) == cfg, f"{name} settings drifted")
    else:
        write_json(settings_path, cfg)
    job_path = create_job(
        settings_path, method="placo", job_id=name, inputs=inputs,
        outputs=outputs, root=specs, reviewed=True, synthetic=False,
        minimum_free_bytes=2**30, threads=1,
    )
    return job_path, OUT / "native" / name


def estimate() -> None:
    prepare()
    path, root = job(
        "parameters", settings(mode="estimate", pair_file=str(PAIR / "pair.tsv.gz"), min_variants=1_000_000),
        {"pair": PAIR / "pair.tsv.gz", "pair_lock": PAIR_LOCK, "placo_source": SOURCE},
        {"parameters": "parameters.json", "semantic_status": "status.json"},
    )
    target = root / "parameters"
    if target.exists():
        verify_current_artifact(target)
    else:
        run_job(path, root)
    parameters = read_json(target / "parameters.json")
    require(parameters["method"] == "PLACO_PLUS" and parameters["n_genomewide"] >= 1_000_000, "invalid global fit")
    print(json.dumps({"stage": "PARAMETERS_COMPLETE", "n_genomewide": parameters["n_genomewide"], "receipt_sha256": sha256(target / "receipt.json")}), flush=True)


def chunk_paths() -> list[tuple[str, Path, int]]:
    entries = read_json(OUT / "chunks/chunks.json")["chunks"]
    return [(f"chunk_{index:06d}", OUT / "chunks" / entry["path"], entry["rows"]) for index, entry in enumerate(entries)]


def run_chunks(workers: int) -> None:
    require(1 <= workers <= 4, "workers must be 1-4")
    estimate()
    params = OUT / "native/parameters/parameters/parameters.json"
    queue = []
    for name, chunk, rows in chunk_paths():
        path, root = job(
            name, settings(mode="chunk", chunk_file=str(chunk), parameters=str(params)),
            {"chunk": chunk, "parameters": params, "placo_source": SOURCE},
            {"chunk": "chunk.tsv", "errors": "numerical_errors.json", "semantic_status": "status.json"},
        )
        target = root / name
        if target.exists():
            receipt = verify_current_artifact(target)
            require(receipt["native_status"]["n_rows"] == rows, f"{name} row count changed")
        else:
            queue.append((name, path, root, rows))
    print(json.dumps({"stage": "CHUNKS_RUNNING", "workers": workers, "remaining": len(queue), "total": len(chunk_paths())}), flush=True)

    def execute(item: tuple[str, Path, Path, int]) -> str:
        name, path, root, rows = item
        receipt = run_job(path, root)
        require(receipt["native_status"]["n_rows"] == rows, f"{name} row count changed")
        return name

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(execute, item) for item in queue]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            name = future.result()
            print(json.dumps({"stage": "CHUNK_COMPLETE", "name": name, "new_completed": index, "total": len(chunk_paths())}), flush=True)


def collate() -> None:
    freeze()
    entries = []
    total = 0
    for name, _chunk, rows in chunk_paths():
        target = OUT / "native" / name / name
        receipt = verify_current_artifact(target)
        require(receipt["native_status"]["n_rows"] == rows, f"{name} row count changed")
        output = target / "chunk.tsv"
        entries.append({"path": str(output), "sha256": sha256(output)})
        total += rows
    require(total == EXPECTED_ROWS, "collation denominator changed")
    manifest = OUT / "completed_chunks.tsv"
    if manifest.exists():
        require(list(read_tsv(manifest, FIELDS)) == entries, "completed chunk manifest changed")
    else:
        write_tsv(manifest, FIELDS, entries)
    target = OUT / "collated"
    if target.exists():
        verify_current_artifact(target)
    else:
        collate_placo(manifest, OUT, "collated", expected_rows=EXPECTED_ROWS, n_pairs=5, synthetic=False)
    status = read_json(target / "status.json")
    require(status["status"] == "PASS", "pair QC did not pass")
    print(json.dumps({"stage": "COMPLETE_SENSITIVITY_ONLY", "status": status, "result_sha256": sha256(target / "results.tsv.gz"), "receipt_sha256": sha256(target / "receipt.json")}), flush=True)


def status() -> None:
    root = OUT / "native"
    chunks = chunk_paths() if (OUT / "chunks/chunks.json").exists() else []
    complete = sum((root / name / name / "receipt.json").exists() for name, _, _ in chunks)
    print(json.dumps({"contract": (OUT / "contract.json").exists(), "prepared_chunks": len(chunks), "parameter_complete": (root / "parameters/parameters/receipt.json").exists(), "complete_chunks": complete, "collated": (OUT / "collated/receipt.json").exists()}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "estimate", "run-chunks", "collate", "status"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.action == "prepare": prepare()
    elif args.action == "estimate": estimate()
    elif args.action == "run-chunks": run_chunks(args.workers)
    elif args.action == "collate": collate()
    else: status()


if __name__ == "__main__":
    main()
