#!/usr/bin/env python3
"""Run the canonical 7 x 2,495 LAVA single-trait family in chromosome batches.

One pinned R process loads LAVA and one chromosome reference, then handles a
whole chromosome batch. Results are immutable locus files with checksummed,
per-cell receipts; invalid or absent loci alone are rescheduled on resume.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
from brain6.io import read_json, sha256

TRAITS = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
SCRIPT = ROOT / "brain6/scripts/run_lava_canonical_batch_v3.R"
PATCH_ID = "LAVA015_BLOCK_REDUCED_SYMMETRY_V1"
PROTECTED_RUN_IDS = {
    "13b80e6b64a5179fec70ba510a89c51241c9d58ec05833cc52161b435aee44fa",
    "af43fde06c4c6c0d57a6ab104bea5aa33f24fb2ff2b89b53163529160452d0cf",
}
EXPECTED_THRESHOLD = 0.05 / 17465
CELL_FIELDS = ("phen", "locus_id", "status", "n_snps", "n_components", "h2.obs", "h2.latent", "p", "reason")
VALID_STATUSES = {"TESTED", "NOT_RUN", "FAILED"}


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def isolated_output_root(path: Path) -> bool:
    normalized = path.resolve()
    if any(run_id in str(normalized) for run_id in PROTECTED_RUN_IDS):
        return False
    return normalized.name != "lava-results-v1"


def atomic_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(obj, sort_keys=True, indent=2).encode() + b"\n"
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload); handle.flush(); os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        try: os.unlink(tmp)
        except FileNotFoundError: pass
        raise


def read_loci(path: Path, expected_hash: str, expected_n: int) -> list[dict[str, str]]:
    if sha256(path) != expected_hash:
        raise ValueError("Frozen LAVA locus-definition checksum mismatch")
    with path.open(encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    if not lines or lines[0] != ["LOC", "CHR", "START", "STOP"]:
        raise ValueError("Locus-definition header drifted")
    values = [dict(zip(lines[0], row)) for row in lines[1:]]
    ids = [r["LOC"] for r in values]
    if len(values) != expected_n or len(ids) != len(set(ids)) or any(len(r) != 4 for r in lines[1:]):
        raise ValueError("Frozen LAVA locus inventory is malformed")
    return values


def cell_hash(row: dict[str, str]) -> str:
    normalized = {k: row.get(k, "") for k in CELL_FIELDS}
    return hashlib.sha256(canonical_json(normalized)).hexdigest()


def canonical_manifest_rows(loci: list[dict[str, str]], traits: tuple[str, ...] = TRAITS) -> list[dict[str, Any]]:
    return [{"trait_id": trait, "locus_id": locus["LOC"], "chromosome": locus["CHR"], "cell_index": index}
            for index, (locus, trait) in enumerate(((loc, tr) for loc in loci for tr in traits), 1)]


def validate_worker_count(execution: dict[str, Any], requested: int) -> int:
    frozen = execution.get("worker_count")
    if frozen != 4 or execution.get("execution_policy", {}).get("worker_count") != 4:
        raise ValueError("The frozen canonical execution lock requires four workers")
    if requested != frozen:
        raise ValueError(f"Requested worker count {requested} differs from frozen setting {frozen}")
    return frozen


def chromosome_batches(loci: list[dict[str, str]], batch_size: int) -> list[tuple[int, int, list[dict[str, str]]]]:
    if batch_size < 1:
        raise ValueError("batch_size_loci must be positive")
    grouped: dict[int, list[dict[str, str]]] = {}
    for locus in loci:
        grouped.setdefault(int(locus["CHR"]), []).append(locus)
    return [(chrom, offset // batch_size + 1, chrom_loci[offset:offset + batch_size])
            for chrom, chrom_loci in sorted(grouped.items())
            for offset in range(0, len(chrom_loci), batch_size)]


def pair_gate(status_a: str, p_a: float | None, status_b: str, p_b: float | None,
              threshold: float = EXPECTED_THRESHOLD) -> tuple[bool, str]:
    if status_a != "TESTED" or status_b != "TESTED" or p_a is None or p_b is None:
        return False, "UNIVARIATE_UNDERPOWERED"
    if not (0 <= p_a <= 1 and 0 <= p_b <= 1):
        raise ValueError("Canonical univariate p-values must lie in [0, 1]")
    if not (p_a < threshold and p_b < threshold):
        return False, "UNIVARIATE_UNDERPOWERED"
    return True, "ELIGIBLE"


def untested_fraction(missing_cells: int, status_counts: dict[str, int], intended: int) -> float:
    if intended <= 0 or missing_cells < 0:
        raise ValueError("Family denominator and missing count must be valid")
    return (missing_cells + status_counts.get("NOT_RUN", 0) + status_counts.get("FAILED", 0)) / intended


def load_cells(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def file_identity(path: Path) -> dict[str, Any]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def result_path(run_dir: Path, locus_id: str) -> Path:
    return run_dir / "cells" / f"locus_{locus_id}.tsv"


def receipt_path(run_dir: Path, locus_id: str) -> Path:
    return run_dir / "receipts" / f"locus_{locus_id}.json"


def verify_locus(run_dir: Path, locus: dict[str, str], traits: tuple[str, ...],
                 analysis_id: str, family_sha: str, execution_sha: str,
                 input_root: Path, reference_prov: Path, script_sha: str) -> bool:
    locus_id = locus["LOC"]
    rp, out = receipt_path(run_dir, locus_id), result_path(run_dir, locus_id)
    if not rp.is_file() or not out.is_file():
        return False
    try:
        receipt = read_json(rp)
        if receipt.get("analysis_id") != analysis_id or receipt.get("family_lock_sha256") != family_sha:
            return False
        if receipt.get("execution_lock_sha256") != execution_sha or receipt.get("script_sha256") != script_sha:
            return False
        config_path = Path(receipt.get("worker_config_path", ""))
        if not config_path.is_file() or receipt.get("config_sha256") != sha256(config_path):
            return False
        config = read_json(config_path)
        if config.get("analysis_id") != analysis_id or locus_id not in config.get("locus_ids", []):
            return False
        if receipt.get("output") != file_identity(out):
            return False
        rows = load_cells(out)
        if len(rows) != len(traits) or {r["phen"] for r in rows} != set(traits) or any(r["locus_id"] != locus_id for r in rows):
            return False
        expected = receipt.get("cells", {})
        if set(expected) != set(traits):
            return False
        input_info = input_root / f"locus_{locus_id}/input_info.tsv"
        if receipt.get("input_info") != file_identity(input_info):
            return False
        if set(receipt.get("sumstats", {})) != set(traits):
            return False
        for trait, identity in receipt["sumstats"].items():
            if identity != file_identity(input_root / f"locus_{locus_id}/{trait}.sumstats.tsv.gz"):
                return False
        if receipt.get("reference_provenance") != file_identity(reference_prov):
            return False
        for row in rows:
            if row.get("status") not in VALID_STATUSES or expected.get(row["phen"]) != cell_hash(row):
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def persist_receipt(run_dir: Path, locus: dict[str, str], rows: list[dict[str, str]],
                    traits: tuple[str, ...], analysis_id: str, family_sha: str,
                    execution_sha: str, input_root: Path, reference_prov: Path,
                    script_sha: str, config_path: Path) -> None:
    locus_id = locus["LOC"]
    if len(rows) != len(traits) or {r["phen"] for r in rows} != set(traits):
        raise ValueError(f"Worker output does not cover all seven cells for locus {locus_id}")
    if any(r.get("status") not in VALID_STATUSES for r in rows):
        raise ValueError(f"Worker returned an unknown cell status for locus {locus_id}")
    output = result_path(run_dir, locus_id)
    receipt = {
        "schema_version": 1, "analysis_id": analysis_id,
        "family_lock_sha256": family_sha, "execution_lock_sha256": execution_sha,
        "locus_id": locus_id, "chromosome": int(locus["CHR"]),
        "worker_config_path": str(config_path), "config_sha256": sha256(config_path),
        "script_sha256": script_sha,
        "runtime": {"R": "4.3.3", "LAVA": "0.1.5", "roundoff_patch": PATCH_ID,
                    "blas_threads": 1, "openmp_threads": 1},
        "input_info": file_identity(input_root / f"locus_{locus_id}/input_info.tsv"),
        "sumstats": {trait: file_identity(input_root / f"locus_{locus_id}/{trait}.sumstats.tsv.gz") for trait in traits},
        "reference_provenance": file_identity(reference_prov),
        "output": file_identity(output),
        "cells": {r["phen"]: cell_hash(r) for r in rows}, "state": "COMPLETE",
    }
    atomic_json(receipt_path(run_dir, locus_id), receipt)


def verify_materialized_locus(input_root: Path, records: dict, locus_id: str) -> None:
    for trait in (*TRAITS, "input_info"):
        record = records.get((locus_id, trait))
        if record is None:
            raise ValueError(f"Materialized-input provenance omits {trait}/{locus_id}")
        path = Path(record["path"])
        if path != (input_root / f"locus_{locus_id}/{trait + '.sumstats.tsv.gz' if trait != 'input_info' else 'input_info.tsv'}"):
            raise ValueError(f"Materialized-input path differs from its provenance: {path}")
        if path.stat().st_size != int(record["bytes"]) or sha256(path) != record["sha256"]:
            raise ValueError(f"Materialized input checksum mismatch: {path}")


def run_batch_attempt(chrom: int, batch_index: int, loci: list[dict[str, str]], run_dir: Path,
              input_root: Path, reference_root: Path, rscript: Path,
              cfg_base: dict) -> dict[str, Any]:
    todo = [l for l in loci if not verify_locus(run_dir, l, tuple(TRAITS),
        cfg_base["analysis_id"], cfg_base["family_sha"], cfg_base["execution_sha"],
        input_root, reference_root / "reference.provenance.json", cfg_base["script_sha"])]
    if not todo:
        return {"batch_id": f"chr{chrom}_batch{batch_index:04d}", "chromosome": chrom,
                "loci": 0, "resumed": len(loci), "wall_seconds": 0.0}
    batch_id = f"chr{chrom}_batch{batch_index:04d}"
    output_dir = run_dir / "worker_output" / batch_id
    output_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "analysis_id": cfg_base["analysis_id"], "locus_ids": [l["LOC"] for l in todo],
        "trait_ids": list(TRAITS), "loci_file": str(cfg_base["loci_file"]),
        "input_root": str(input_root), "reference_prefix": str(reference_root / f"lava-ukb-v1.1_chr{chrom}"),
        "output_dir": str(run_dir), "random_seed": 20260922,
        "runtime_validator_passed": True,
        "runtime_validation_output_sha256": cfg_base["runtime_validation_sha"],
        "execution_policy": cfg_base["execution_policy"],
    }
    config_hash = hashlib.sha256(canonical_json(config)).hexdigest()
    config_path = output_dir / f"batch_{config_hash}.json"
    if config_path.exists():
        if read_json(config_path) != config:
            raise ValueError(f"Immutable worker config collision: {config_path}")
    else:
        atomic_json(config_path, config)
    config_sha = sha256(config_path)
    env = os.environ.copy()
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1", "BLIS_NUM_THREADS": "1"})
    command = [str(rscript), str(SCRIPT), str(config_path)]
    started = time.monotonic()
    log_path = output_dir / f"batch_{config_hash}.log"
    proc = subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, bufsize=1)
    assert proc.stdout is not None
    completed: set[str] = set()
    runtime_ready_seconds = None
    reference_load_seconds = None
    input_hash_seconds = 0.0
    input_hash_bytes = 0
    locus_metrics: list[dict[str, Any]] = []
    try:
        with log_path.open("a", encoding="utf-8") as log:
          for line in proc.stdout:
            log.write(line); log.flush()
            sys.stdout.write(f"[lava-v3 {batch_id}] {line}"); sys.stdout.flush()
            if line.startswith("BRAIN6_WORKER_READY "):
                runtime_ready_seconds = time.monotonic() - started
                for field in line.split()[1:]:
                    if field.startswith("reference_load_seconds="):
                        reference_load_seconds = float(field.split("=", 1)[1])
            if line.startswith("BRAIN6_CANONICAL_LOCUS_METRICS "):
                parts = line.split()
                metric = {"locus_id": parts[1]}
                metric.update({k: float(v) for k, v in (part.split("=", 1) for part in parts[2:])})
                locus_metrics.append(metric)
            if line.startswith("BRAIN6_CANONICAL_LOCUS_COMPLETE "):
                locus_id = line.split()[-1]
                locus = next(l for l in todo if l["LOC"] == locus_id)
                hash_started = time.monotonic()
                verify_materialized_locus(input_root, cfg_base["input_records"], locus_id)
                input_hash_seconds += time.monotonic() - hash_started
                input_hash_bytes += sum(int(cfg_base["input_records"][(locus_id, trait)]["bytes"])
                    for trait in (*TRAITS, "input_info"))
                out = result_path(run_dir, locus_id)
                rows = load_cells(out)
                persist_receipt(run_dir, locus, rows, tuple(TRAITS), cfg_base["analysis_id"],
                    cfg_base["family_sha"], cfg_base["execution_sha"], input_root,
                    reference_root / "reference.provenance.json", cfg_base["script_sha"], config_path)
                completed.add(locus_id)
    finally:
        exit_code = proc.wait()
    elapsed = time.monotonic() - started
    if exit_code != 0:
        raise RuntimeError(f"R batch failed {batch_id}; completed {len(completed)}/{len(todo)} loci")
    if completed != {l["LOC"] for l in todo}:
        raise RuntimeError(f"R batch omitted locus completion markers {batch_id}")
    if runtime_ready_seconds is None:
        raise RuntimeError(f"R batch failed to report worker readiness {batch_id}")
    return {"batch_id": batch_id, "chromosome": chrom, "loci": len(todo), "resumed": len(loci)-len(todo),
            "wall_seconds": elapsed, "worker_startup_seconds": max(0.0, runtime_ready_seconds-(reference_load_seconds or 0.0)),
            "reference_load_seconds": reference_load_seconds, "input_checksum_seconds": input_hash_seconds,
            "input_checksum_bytes": input_hash_bytes, "locus_metrics": locus_metrics,
            "worker_config_path": str(config_path), "worker_config_sha256": config_sha,
            "worker_log_path": str(log_path), "worker_log_sha256": sha256(log_path)}


def run_batch(chrom: int, batch_index: int, loci: list[dict[str, str]], run_dir: Path,
              input_root: Path, reference_root: Path, rscript: Path,
              cfg_base: dict) -> dict[str, Any]:
    attempts = int(cfg_base.get("max_attempts_per_locus", 1))
    if attempts < 1: raise ValueError("max_attempts_per_locus must be positive")
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            result = run_batch_attempt(chrom, batch_index, loci, run_dir, input_root,
                reference_root, rscript, cfg_base)
            result["attempts"] = attempt
            return result
        except (RuntimeError, OSError) as exc:
            last_error = exc
    raise RuntimeError(f"Batch chr{chrom}_batch{batch_index:04d} failed after {attempts} attempts: {last_error}") from last_error


def audit_family(run_dir: Path, loci: list[dict[str, str]], traits: tuple[str, ...],
                 analysis_id: str, family_sha: str, execution_sha: str,
                 input_root: Path, reference_prov: Path, script_sha: str) -> dict[str, Any]:
    intended = {(trait, locus["LOC"]) for locus in loci for trait in traits}
    found: set[tuple[str, str]] = set()
    counts = {s: 0 for s in VALID_STATUSES}
    invalid = []
    for locus in loci:
        out = result_path(run_dir, locus["LOC"])
        if not verify_locus(run_dir, locus, traits, analysis_id, family_sha, execution_sha,
                            input_root, reference_prov, script_sha):
            invalid.append(locus["LOC"]); continue
        rows = load_cells(out)
        for row in rows:
            key = (row["phen"], row["locus_id"])
            if key in found: raise ValueError(f"Duplicate canonical cell identity: {key}")
            found.add(key); counts[row["status"]] += 1
    missing = intended - found
    untested = len(missing) + counts["FAILED"] + counts["NOT_RUN"]
    return {"intended_cells": len(intended), "valid_cells": len(found), "missing_cells": len(missing),
        "invalid_loci": invalid, "duplicate_cells": 0, "statuses": counts,
        "untested_fraction": untested / len(intended), "frozen_max_untested_fraction": 0.05,
        "qc_pass": untested / len(intended) <= 0.05}


def write_family_aggregate(run_dir: Path, all_loci: list[dict[str, str]], selected_loci: list[dict[str, str]],
                           traits: tuple[str, ...], analysis_id: str, family_sha: str, execution_sha: str,
                           input_root: Path, reference_prov: Path, script_sha: str) -> dict[str, Any]:
    selected = {l["LOC"]: l for l in selected_loci}
    path = run_dir / "results" / "canonical_family_results.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    fields = ["phen", "locus_id", "chromosome", "status", "n_snps", "n_components",
              "h2.obs", "h2.latent", "p", "reason"]
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for locus in all_loci:
                valid = locus["LOC"] in selected and verify_locus(run_dir, locus, traits, analysis_id,
                    family_sha, execution_sha, input_root, reference_prov, script_sha)
                rows = {r["phen"]: r for r in load_cells(result_path(run_dir, locus["LOC"]))} if valid else {}
                for trait in traits:
                    row = rows.get(trait)
                    writer.writerow({"phen": trait, "locus_id": locus["LOC"], "chromosome": locus["CHR"],
                        "status": row["status"] if row else "MISSING",
                        "n_snps": row.get("n_snps", "") if row else "",
                        "n_components": row.get("n_components", "") if row else "",
                        "h2.obs": row.get("h2.obs", "") if row else "",
                        "h2.latent": row.get("h2.latent", "") if row else "",
                        "p": row.get("p", "") if row else "",
                        "reason": row.get("reason", "") if row else
                            ("MISSING_OR_INVALID_RECEIPT" if locus["LOC"] in selected else "NOT_IN_EXECUTION_SCOPE")})
            handle.flush(); os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        try: os.unlink(tmp)
        except FileNotFoundError: pass
        raise
    return {"path": str(path), "sha256": sha256(path), "rows": len(all_loci) * len(traits)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family-lock", type=Path, required=True)
    parser.add_argument("--execution-lock", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-id", help="Optional assertion against the computed deterministic run hash")
    parser.add_argument("--rscript", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--locus-limit", type=int, help="Diagnostic subset only; never use for family promotion")
    parser.add_argument("--locus-list", type=Path, help="Diagnostic subset ids, one per line; never use for family promotion")
    args = parser.parse_args()
    if args.workers < 1: parser.error("--workers must be positive")
    family, execution = read_json(args.family_lock), read_json(args.execution_lock)
    if not isolated_output_root(args.output_root):
        raise ValueError("Corrected LAVA v3 requires a distinct output root; baseline and roundoff runs are protected")
    family_sha, execution_sha = sha256(args.family_lock), sha256(args.execution_lock)
    if family.get("analysis_id") != "brain6-lava-canonical-v3" or execution.get("analysis_id") != family["analysis_id"]:
        raise ValueError("LAVA v3 family/execution lock identity mismatch")
    if execution.get("family_lock_sha256") != family_sha or execution.get("gate_p_strictly_less_than") != EXPECTED_THRESHOLD:
        raise ValueError("Execution lock does not bind the exact canonical family/gate")
    canonical = family.get("canonical_univariate", {})
    if family.get("trait_ids") != list(TRAITS) or canonical.get("n_tests") != 17465:
        raise ValueError("Canonical family must be the exact ordered 7 × 2,495 family")
    validate_worker_count(execution, args.workers)
    loci_file = ROOT / family["locus_definition"]["path"]
    all_loci = read_loci(loci_file, family["locus_definition"]["sha256"], 2495)
    loci = all_loci
    if args.locus_list:
        requested = [line.strip() for line in args.locus_list.read_text().splitlines() if line.strip()]
        if len(requested) != len(set(requested)):
            raise ValueError("Diagnostic locus list has duplicate ids")
        by_id = {r["LOC"]: r for r in all_loci}
        if any(x not in by_id for x in requested): raise ValueError("Diagnostic locus list contains unknown id")
        loci = [by_id[x] for x in requested]
    if args.locus_limit:
        loci = loci[:args.locus_limit]
    ref_prov = args.reference_root / "reference.provenance.json"
    if not ref_prov.is_file(): raise FileNotFoundError(ref_prov)
    input_prov = read_json(args.input_root / "provenance.json")
    if input_prov.get("n_loci") != 2495: raise ValueError("Materialized input inventory is not the exact 2,495-locus source")
    if input_prov.get("n_traits") != 7 or len(input_prov.get("records", [])) != 2495 * 8:
        raise ValueError("Materialized-input provenance is not the exact 2,495 × (7 traits + input-info) inventory")
    validator = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
    runtime = subprocess.run([str(args.rscript), str(validator)], cwd=ROOT,
        capture_output=True, text=True, check=False)
    if runtime.returncode != 0 or runtime.stdout.count("TRACK_B_LAVA_RUNTIME\tstatus=PASS\t") != 1:
        raise RuntimeError("Pinned LAVA runtime preflight failed: " + (runtime.stdout + runtime.stderr)[-4000:])
    reference_check = subprocess.run([sys.executable, "scripts/lava_contract.py", "--verify-reference"],
        cwd=ROOT, capture_output=True, text=True, check=False)
    if reference_check.returncode != 0 or "LAVA_REFERENCE_VALIDATED" not in reference_check.stdout:
        raise RuntimeError("Official LAVA reference preflight failed: " +
            (reference_check.stdout + reference_check.stderr)[-4000:])
    run_identity = {"analysis_id": family["analysis_id"], "family_lock_sha256": family_sha,
        "execution_lock_sha256": execution_sha, "input_provenance_sha256": sha256(args.input_root / "provenance.json"),
        "reference_provenance_sha256": sha256(ref_prov), "locus_definition_sha256": sha256(loci_file),
        "batch_worker_sha256": sha256(SCRIPT), "coordinator_sha256": sha256(Path(__file__)),
        "runtime_validator_sha256": sha256(validator),
        "runtime_validation_output_sha256": hashlib.sha256(runtime.stdout.encode("utf-8")).hexdigest(),
        "reference_validation_output_sha256": hashlib.sha256(reference_check.stdout.encode("utf-8")).hexdigest(),
        "rscript_path": str(args.rscript.resolve()), "rscript_binary_sha256": sha256(args.rscript.resolve()),
        "R": "4.3.3", "LAVA": "0.1.5", "roundoff_patch": PATCH_ID}
    computed_run_id = hashlib.sha256(canonical_json(run_identity)).hexdigest()
    if args.run_id and args.run_id != computed_run_id:
        raise ValueError("Supplied run id differs from the computed checksum-bound identity")
    args.run_id = computed_run_id
    run_dir = args.output_root / args.run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg_base = {"analysis_id": family["analysis_id"], "family_sha": family_sha,
        "execution_sha": execution_sha, "loci_file": loci_file,
        "script_sha": sha256(SCRIPT),
        "runtime_validation_sha": hashlib.sha256(runtime.stdout.encode("utf-8")).hexdigest(),
        "input_provenance": input_prov,
        "input_records": {(str(r.get("locus_id")), r.get("trait", "input_info")): r
            for r in input_prov.get("records", []) if r.get("kind") in {"sumstats", "input_info"}},
        "execution_policy": execution["execution_policy"],
        "max_attempts_per_locus": execution.get("max_attempts_per_locus", 1)}
    manifest_path = run_dir / "canonical_family_manifest.tsv"
    rows = canonical_manifest_rows(all_loci)
    import io
    manifest_buffer = io.StringIO(newline="")
    manifest_writer = csv.DictWriter(manifest_buffer, fieldnames=["trait_id", "locus_id", "chromosome", "cell_index"],
        delimiter="\t", lineterminator="\n")
    manifest_writer.writeheader(); manifest_writer.writerows(rows)
    manifest_bytes = manifest_buffer.getvalue().encode("utf-8")
    if manifest_path.exists():
        if manifest_path.read_bytes() != manifest_bytes:
            raise ValueError("Existing canonical manifest differs from the deterministic frozen family")
    else:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with manifest_path.open("x", encoding="utf-8", newline="") as handle:
            handle.write(manifest_bytes.decode("utf-8"))
    manifest_sha = sha256(manifest_path)
    family_lock_copy = run_dir / "family_lock.json"
    if family_lock_copy.exists():
        if read_json(family_lock_copy) != family: raise ValueError("Existing family lock copy changed")
    else: atomic_json(family_lock_copy, family)
    execution_copy = run_dir / "execution_lock.json"
    if execution_copy.exists():
        if read_json(execution_copy) != execution: raise ValueError("Existing execution lock copy changed")
    else: atomic_json(execution_copy, execution)
    run_identity_path = run_dir / "run_identity.json"
    identity_doc = {"analysis_id": family["analysis_id"], "run_id": args.run_id,
        "family_lock_sha256": family_sha, "execution_lock_sha256": execution_sha,
        "run_identity_sha256": args.run_id,
        "canonical_manifest_sha256": manifest_sha, "script_sha256": cfg_base["script_sha"],
        "runtime": {"R": "4.3.3", "LAVA": "0.1.5", "BLAS_threads": 1},
        "runtime_validator_sha256": sha256(validator),
        "runtime_validation_output_sha256": cfg_base["runtime_validation_sha"]}
    if run_identity_path.exists():
        if read_json(run_identity_path) != identity_doc: raise ValueError("Existing run identity changed")
    else: atomic_json(run_identity_path, identity_doc)
    batch_size = int(execution.get("batch_size_loci", 10))
    if batch_size < 1:
        raise ValueError("Frozen execution lock batch_size_loci must be positive")
    batches = chromosome_batches(loci, batch_size)
    # Balance heterogeneous chromosome blocks by estimated input bytes; the pool
    # dynamically takes the next largest batch as each worker becomes free.
    records = {(str(r.get("locus_id")), r.get("trait", "input_info")): r
               for r in input_prov.get("records", []) if r.get("kind") in {"sumstats", "input_info"}}
    def estimated_cost(batch: tuple[int, int, list[dict[str, str]]]) -> int:
        return sum(int(records.get((l["LOC"], trait), {}).get("bytes", 0))
                   for l in batch[2] for trait in (*TRAITS, "input_info"))
    batches.sort(key=lambda b: (-estimated_cost(b), b[0], b[1]))
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_batch, chrom, batch_index, group, run_dir, args.input_root,
            args.reference_root, args.rscript, cfg_base)
            for chrom, batch_index, group in batches]
        for future in concurrent.futures.as_completed(futures): results.append(future.result())
    children_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    audit = audit_family(run_dir, loci, TRAITS, family["analysis_id"], family_sha,
        execution_sha, args.input_root, ref_prov, cfg_base["script_sha"])
    aggregate = write_family_aggregate(run_dir, all_loci, loci, TRAITS, family["analysis_id"],
        family_sha, execution_sha, args.input_root, ref_prov, cfg_base["script_sha"])
    if len(loci) != 2495:
        audit["qc_pass"] = False
        audit["promotion_block"] = "diagnostic subset; exact 17,465-cell family not executed"
    report = {"analysis_id": family["analysis_id"], "run_id": args.run_id,
        "run_identity": run_identity,
        "workers": args.workers, "wall_seconds": time.monotonic()-started,
        "child_cpu_seconds": (children_after.ru_utime-children_before.ru_utime)+
            (children_after.ru_stime-children_before.ru_stime),
        "peak_child_rss_bytes": int(children_after.ru_maxrss),
        "chromosome_batches": sorted(results, key=lambda x: (x["chromosome"], x["batch_id"])),
        "manifest_sha256": manifest_sha, "aggregate": aggregate, "audit": audit}
    atomic_json(run_dir / "latest_audit.json", report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if audit["qc_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
