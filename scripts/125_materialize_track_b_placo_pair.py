#!/usr/bin/env python3
"""Materialize, lock, and publish one exact Track B PLACO+ pair family."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import itertools
import json
import math
import os
import re
import resource
import shutil
import sqlite3
import stat
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterator


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_SPEC = importlib.util.spec_from_file_location(
    "track_b_pleiotropy_contract_runtime", ROOT / "scripts/123_track_b_pleiotropy_contract.py"
)
if CONTRACT_SPEC is None or CONTRACT_SPEC.loader is None:
    raise RuntimeError("could not load Track B pleiotropy contract")
contract = importlib.util.module_from_spec(CONTRACT_SPEC)
CONTRACT_SPEC.loader.exec_module(contract)

POLICY = Path("config/track_b_pleiotropy_policy.json")
CONTRACT_LOCK = Path("results/track_b/pleiotropy/contract.lock.json")
INPUT_GATE = Path("results/track_b/pleiotropy/input_gate.tsv")
READINESS_GATE = Path("results/track_b/pleiotropy/readiness_gate.tsv")
INPUT_LOCK = Path("results/track_b/pleiotropy/input_gate.lock.json")
PLACO_WORK_ROOT = Path("work/track_b_pleiotropy/placo")
CANONICAL_ROOT = Path("results/track_b/pleiotropy/results/placo")
RUNNER = Path("scripts/126_run_track_b_placo_pair.R")
MATERIALIZER = Path("scripts/125_materialize_track_b_placo_pair.py")

PAIR_IDS = ("A", "B", "CONTROL")
PAIR_IDENTITIES = {
    "A": ("snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    "B": ("insomnia", "adhd", "PRIMARY_DISCOVERY"),
    "CONTROL": ("insomnia", "frailty", "POSITIVE_CONTROL"),
}
SOURCE_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
ALIGNED_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "FRQ1", "FRQ2", "Z1", "Z2", "P1", "P2"]
LEDGER_FIELDS = [
    "analysis_id", "pair_id", "SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2",
    "T_PLACO_PLUS", "P_PLACO_PLUS", "PLACO_BH_Q", "within_pair_family_n",
    "analysis_status", "numerical_error",
]
SUMMARY_FIELDS = [
    "analysis_id", "pair_id", "run_fingerprint", "aligned_rows", "within_pair_family_n",
    "nuisance_null_rows_variance", "nuisance_null_rows_correlation", "VarZ1", "VarZ2", "CorZ",
    "nuisance_checkpoint_sha256",
    "numerical_failure_count", "numerical_failure_fraction", "minimum_p_placo_plus",
    "primary_headline_count", "pairwise_gws_count", "within_pair_bh_count", "terminal_status",
    "shard_count", "shard_size", "workers", "elapsed_seconds", "ledger", "ledger_bytes",
    "ledger_sha256", "ledger_schema", "task_sha256", "aligned_sha256", "placo_source_sha256",
    "materializer_sha256", "runner_sha256", "r_version", "data_table_version", "runtime_platform",
]
TASK_FIELDS = [
    "analysis_id", "pair_id", "trait1", "trait2", "family_role", "run_fingerprint",
    "aligned_input", "aligned_input_sha256", "aligned_input_rows", "aligned_input_schema",
    "materialization_provenance", "materialization_provenance_sha256", "policy", "policy_sha256",
    "contract_lock", "contract_lock_sha256", "input_gate_lock", "input_gate_lock_sha256",
    "placo_source", "placo_source_sha256", "marginal_p_threshold", "z_squared_maximum",
    "absolute_tolerance", "maximum_numerical_failure_fraction", "primary_headline_threshold",
    "pairwise_gws_threshold", "within_pair_bh_alpha", "minimum_free_bytes", "maximum_workers",
    "pair_concurrency_limit", "shard_size", "checkpoint_dir", "staged_ledger", "run_summary",
    "staged_provenance", "canonical_ledger", "canonical_provenance", "benchmark_raw",
    "benchmark_out", "materializer_sha256", "runner_sha256",
]
BENCHMARK_RAW_FIELDS = [
    "analysis_id", "pair_id", "run_fingerprint", "input_rows", "input_sha256", "task_sha256",
    "policy_sha256", "contract_lock_sha256", "input_gate_lock_sha256", "materializer_sha256",
    "placo_source_sha256", "nuisance_checkpoint_sha256", "nuisance_null_rows_variance",
    "nuisance_null_rows_correlation", "VarZ1", "VarZ2", "CorZ", "global_nuisance_input_rows",
    "global_nuisance_elapsed_seconds", "global_nuisance_checkpoint_reused", "benchmark_variants",
    "workers", "numerical_failures", "single_variant_testing_elapsed_seconds",
    "measured_runner_elapsed_seconds", "variants_per_second", "projected_full_family_seconds",
    "projected_full_family_hours", "scientific_equivalence", "runner_sha256", "r_version",
    "data_table_version", "runtime_platform",
]
REQUIRED_ALIGNMENT_COUNTS = [
    "trait1_rows", "trait2_rows", "exact_rsid_coordinate_matches", "allele_matches", "allele_flips",
    "allele_mismatches_dropped", "duplicate_ids_dropped", "invalid_z_or_p_dropped",
    "z_squared_exclusions", "eligible_written",
]
SHA256_RE = re.compile(r"[0-9a-f]{64}")
RSID_RE = re.compile(r"rs[0-9]+", re.IGNORECASE)
BASES = {"A", "C", "G", "T"}
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
MINIMUM_FREE_BYTES = 8 * 1024**3
MAX_WORKERS = 4
SHARD_SIZE = 20_000


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def safe_relative(root: Path, value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        fail(f"unsafe repository-relative path: {value}")
    candidate = (root / path).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError:
        fail(f"repository-relative path escapes through a symlink: {value}")
    return candidate


def relative(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        fail(f"path is outside the execution root: {path}")


def identity_fields(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def file_identity(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size == 0:
                fail(f"missing non-empty regular artifact: {path}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        fail(f"could not hash artifact {path}: {error}")
    if identity_fields(before) != identity_fields(after) or identity_fields(after) != identity_fields(current):
        fail(f"artifact changed while hashing: {path}")
    return after.st_size, digest.hexdigest()


def sha256(path: Path) -> str:
    return file_identity(path)[1]


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        fail(f"unreadable JSON artifact {path}: {error}")
    if not isinstance(value, dict):
        fail(f"JSON artifact must contain one object: {path}")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty TSV artifact: {path}")
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            if any(None in row or any(value is None for value in row.values()) for row in rows):
                fail(f"ragged TSV artifact: {path}")
            return list(reader.fieldnames or []), rows
    except (OSError, UnicodeError, csv.Error) as error:
        fail(f"unreadable TSV artifact {path}: {error}")


def one_tsv_row(path: Path) -> dict[str, str]:
    _, rows = read_tsv(path)
    if len(rows) != 1:
        fail(f"TSV artifact must contain exactly one row: {path}")
    return rows[0]


def tsv_bytes(fields: list[str], rows: list[dict[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def exclusive_family(sources: list[Path], destinations: list[Path]) -> None:
    if len(sources) != len(destinations):
        fail("internal no-replace publication family length mismatch")
    if any(path.exists() for path in destinations):
        fail("immutable publication destination already exists; overwrite is forbidden")
    published: list[tuple[Path, int, int]] = []
    try:
        for source, destination in zip(sources, destinations, strict=True):
            destination.parent.mkdir(parents=True, exist_ok=True)
            source_stat = source.stat()
            try:
                os.link(source, destination)
            except FileExistsError:
                fail(f"publication destination appeared concurrently: {destination}")
            destination_stat = destination.stat()
            if (destination_stat.st_dev, destination_stat.st_ino) != (source_stat.st_dev, source_stat.st_ino):
                fail(f"no-replace publication identity mismatch: {destination}")
            published.append((destination, source_stat.st_dev, source_stat.st_ino))
    except BaseException:
        for path, device, inode in published:
            try:
                observed = path.stat()
            except FileNotFoundError:
                continue
            if (observed.st_dev, observed.st_ino) == (device, inode):
                path.unlink()
        raise


def physical_memory_bytes() -> int:
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (OSError, TypeError, ValueError) as error:
        fail(f"could not determine physical memory for PLACO+: {error}")


def assert_upstream_ready(row: dict[str, str]) -> None:
    blockers = []
    if row.get("dense_input_gate") != "PASS_ALL_THREE_PAIRS":
        blockers.append("DENSE_INPUT")
    if not row.get("replication_gate", "").startswith("PASS_TERMINAL"):
        blockers.append("PRIMARY_REPLICATION")
    if not row.get("local_analysis_gate", "").startswith("PASS_TERMINAL"):
        blockers.append("LOCAL_ANALYSIS")
    if not row.get("software_gate", "").startswith("READY_PINNED_PLACO_SOURCE"):
        blockers.append("PLACO_SOURCE")
    if blockers:
        fail("BLOCKED_UPSTREAM Track B PLACO+ scan prerequisites are not terminal: " + ",".join(blockers))


def production_context(root: Path, pair_id: str) -> dict[str, Any]:
    if pair_id not in PAIR_IDS:
        fail("pair_id must be exactly A, B, or CONTROL")
    policy = contract.validate_policy(root)
    dense_contract = policy.get("dense_input_contract", {})
    plausibility = dense_contract.get("pair_materialization_plausibility", {})
    if (
        dense_contract.get("required_columns") != SOURCE_FIELDS
        or dense_contract.get("alignment_provenance_counts") != REQUIRED_ALIGNMENT_COUNTS
        or policy.get("required_future_result_schemas", {}).get("placo_full_p_ledger") != LEDGER_FIELDS
        or policy.get("future_result_paths", {}).get("placo_full_ledger_directory") != str(CANONICAL_ROOT)
        or plausibility.get("required_autosomes") != list(range(1, 23))
        or plausibility.get("nonzero_eligible_variants_on_every_required_autosome") is not True
        or plausibility.get("per_autosome_provenance_required") is not True
    ):
        fail("Track B PLACO+ production schema or materialization-floor contract drifted")
    contract_lock = contract.verify_contract(root)
    input_lock_path = root / INPUT_LOCK
    input_gate_path = root / INPUT_GATE
    readiness_path = root / READINESS_GATE
    input_lock = read_json(input_lock_path)
    if (
        input_lock.get("schema_version") != "sleep-atlas-track-b-pleiotropy-input-gate.1"
        or input_lock.get("analysis_id") != policy["analysis_id"]
        or input_lock.get("policy_sha256") != sha256(root / POLICY)
        or input_lock.get("contract_lock_sha256") != sha256(root / CONTRACT_LOCK)
        or input_lock.get("input_gate_sha256") != sha256(input_gate_path)
        or input_lock.get("readiness_gate_sha256") != sha256(readiness_path)
        or input_lock.get("all_pair_dense_input_gates_pass") is not True
        or input_lock.get("scientific_result_count") != 0
    ):
        fail("Track B pleiotropy input-gate lock is incomplete or drifted")
    _, readiness = read_tsv(readiness_path)
    placo_rows = [row for row in readiness if row.get("component_id") == "PLACO_PLUS"]
    if len(placo_rows) != 1:
        fail("Track B readiness gate lacks one PLACO_PLUS row")
    assert_upstream_ready(placo_rows[0])
    _, pair_gate = read_tsv(input_gate_path)
    pair_rows = [row for row in pair_gate if row.get("pair_id") == pair_id]
    if len(pair_rows) != 1 or pair_rows[0].get("input_gate_status") != "READY_DENSE_INPUTS":
        fail(f"Track B pair is not dense-input ready: {pair_id}")
    trait1, trait2, role = PAIR_IDENTITIES[pair_id]
    pair_row = pair_rows[0]
    if (
        (pair_row.get("trait1"), pair_row.get("trait2"), pair_row.get("family_role"))
        != (trait1, trait2, role)
        or pair_row.get("ancestry") != "EUR" or pair_row.get("analysis_build") != "hg19"
        or pair_row.get("full_genome_wide_dense_required") != "TRUE"
        or pair_row.get("hapmap3_only_forbidden") != "TRUE"
        or pair_row.get("required_autosomes") != ",".join(str(value) for value in range(1, 23))
    ):
        fail("Track B PLACO+ pair identity differs from the frozen A/B/CONTROL family")
    scans = input_lock.get("live_dense_inputs", {})
    if not isinstance(scans, dict) or trait1 not in scans or trait2 not in scans:
        fail("Track B input lock lacks the pair's two dense source identities")
    sources = []
    for trait in (trait1, trait2):
        record = scans[trait]
        if (
            not isinstance(record, dict)
            or record.get("schema") != ",".join(SOURCE_FIELDS)
            or record.get("gzip_integrity") != "FULL_DECOMPRESSION_CRC_PASS"
            or not SHA256_RE.fullmatch(str(record.get("sha256", "")))
            or int(record.get("rows", 0)) < int(policy["dense_input_contract"]["minimum_variant_rows_per_trait"])
        ):
            fail(f"frozen dense source identity is incomplete for {trait}")
        path = safe_relative(root, record["path"])
        size, digest = file_identity(path)
        if size != int(record["bytes"]) or digest != record["sha256"]:
            fail(f"live dense source differs from the frozen input gate: {trait}")
        sources.append({**record, "trait_id": trait, "absolute_path": path})
    placo = policy["placo_plus"]
    placo_source = safe_relative(root, placo["source_path"])
    if sha256(placo_source) != placo["source_sha256"]:
        fail("pinned correlated/overlap-aware PLACO+ source differs from policy")
    runner = root / RUNNER
    if not runner.is_file():
        fail("Track B PLACO+ runner is absent")
    fingerprint_payload = {
        "analysis_id": policy["analysis_id"], "pair_id": pair_id, "traits": [trait1, trait2],
        "policy_sha256": sha256(root / POLICY), "contract_lock_sha256": sha256(root / CONTRACT_LOCK),
        "input_gate_lock_sha256": sha256(input_lock_path), "materializer_sha256": sha256(Path(__file__)),
        "runner_sha256": sha256(runner), "placo_source_sha256": placo["source_sha256"],
        "sources": [{key: value for key, value in source.items() if key != "absolute_path"} for source in sources],
        "resource_envelope": {
            "minimum_free_bytes": MINIMUM_FREE_BYTES,
            "maximum_workers": MAX_WORKERS, "pair_concurrency": 1, "shard_size": SHARD_SIZE,
            "ram_rule": "MEASURE_WITH_EXACT_GLOBAL_NUISANCE_BENCHMARK_BEFORE_FREEZING_A_RAM_GATE",
        },
    }
    fingerprint = hashlib.sha256(
        json.dumps(fingerprint_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "root": root, "pair_id": pair_id, "trait1": trait1, "trait2": trait2, "family_role": role,
        "policy": policy, "contract_lock": contract_lock, "input_lock": input_lock,
        "sources": sources, "placo_source": placo_source, "fingerprint": fingerprint,
        "fingerprint_payload": fingerprint_payload,
    }


def parse_source_row(row: dict[str, str]) -> tuple[tuple[object, ...] | None, str | None, int | None]:
    try:
        chromosome = int(row["CHR"])
        position = int(row["BP"])
    except (TypeError, ValueError):
        return None, "identity_or_frequency", None
    if not 1 <= chromosome <= 22 or position <= 0 or not RSID_RE.fullmatch(row.get("SNP", "")):
        return None, "identity_or_frequency", chromosome if 1 <= chromosome <= 22 else None
    a1, a2 = row.get("A1", "").upper(), row.get("A2", "").upper()
    if (
        a1 not in BASES or a2 not in BASES or a1 == a2
        or frozenset((a1, a2)) in PALINDROMIC
    ):
        return None, "identity_or_frequency", chromosome
    try:
        frequency = float(row["FRQ"])
        beta = float(row["BETA"])
        standard_error = float(row["SE"])
        p_value = float(row["P"])
        sample_size = float(row["N"])
    except (TypeError, ValueError):
        return None, "invalid_z_or_p", chromosome
    if not math.isfinite(frequency) or not 0 < frequency < 1 or not math.isfinite(sample_size) or sample_size <= 0:
        return None, "identity_or_frequency", chromosome
    if (
        not math.isfinite(beta) or not math.isfinite(standard_error) or standard_error <= 0
        or not math.isfinite(p_value) or not 0 <= p_value <= 1
    ):
        return None, "invalid_z_or_p", chromosome
    z = beta / standard_error
    if not math.isfinite(z):
        return None, "invalid_z_or_p", chromosome
    return (
        row["SNP"].lower(), chromosome, position, a1, a2, frequency, z, p_value,
    ), None, chromosome


def load_source(
    connection: sqlite3.Connection, table: str, path: Path, expected: dict[str, Any],
) -> dict[str, Any]:
    raw_table = f"{table}_raw"
    connection.execute(
        f"CREATE TABLE {raw_table} (snp TEXT, chr INTEGER, bp INTEGER, a1 TEXT, a2 TEXT, "
        "frq REAL, z REAL, p REAL)"
    )
    statement = f"INSERT INTO {raw_table} VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
    total_rows = 0
    invalid_z_or_p = 0
    invalid_identity_or_frequency = 0
    per_autosome = {chromosome: 0 for chromosome in range(1, 23)}
    batch: list[tuple[object, ...]] = []
    before_identity = file_identity(path)
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or []) != SOURCE_FIELDS:
                fail(f"dense source schema differs from the exact FRQ-bearing contract: {path}")
            for row in reader:
                total_rows += 1
                if None in row or any(row.get(field) is None for field in SOURCE_FIELDS):
                    invalid_identity_or_frequency += 1
                    continue
                parsed, rejection, chromosome = parse_source_row(row)
                if chromosome in per_autosome:
                    per_autosome[chromosome] += 1
                if rejection == "invalid_z_or_p":
                    invalid_z_or_p += 1
                    continue
                if rejection is not None:
                    invalid_identity_or_frequency += 1
                    continue
                batch.append(parsed)
                if len(batch) == 50_000:
                    connection.executemany(statement, batch)
                    batch.clear()
            if batch:
                connection.executemany(statement, batch)
        connection.commit()
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeError, csv.Error) as error:
        fail(f"could not fully scan dense source {path}: {error}")
    after_identity = file_identity(path)
    if before_identity != after_identity:
        fail(f"dense source changed during materialization: {path}")
    if total_rows != int(expected["rows"]) or after_identity != (int(expected["bytes"]), expected["sha256"]):
        fail(f"dense source count/hash/bytes differ from the frozen gate: {path}")
    connection.execute(f"CREATE INDEX {raw_table}_snp ON {raw_table}(snp)")
    duplicate_rows = int(connection.execute(
        f"SELECT COALESCE(SUM(n), 0) FROM (SELECT COUNT(*) AS n FROM {raw_table} GROUP BY snp HAVING n > 1)"
    ).fetchone()[0])
    connection.execute(
        f"CREATE TABLE {table} (snp TEXT PRIMARY KEY, chr INTEGER NOT NULL, bp INTEGER NOT NULL, "
        "a1 TEXT NOT NULL, a2 TEXT NOT NULL, frq REAL NOT NULL, z REAL NOT NULL, p REAL NOT NULL) WITHOUT ROWID"
    )
    connection.execute(
        f"INSERT INTO {table} SELECT snp, chr, bp, a1, a2, frq, z, p FROM {raw_table} "
        "GROUP BY snp HAVING COUNT(*) = 1"
    )
    unique_rows = int(connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
    connection.execute(f"DROP TABLE {raw_table}")
    connection.commit()
    if total_rows != invalid_z_or_p + invalid_identity_or_frequency + duplicate_rows + unique_rows:
        fail(f"dense source rejection accounting does not reconcile: {path}")
    return {
        "trait_rows": total_rows, "invalid_z_or_p_dropped": invalid_z_or_p,
        "invalid_identity_or_frequency_dropped": invalid_identity_or_frequency,
        "duplicate_ids_dropped": duplicate_rows, "valid_unique_rows": unique_rows,
        "per_autosome_rows": per_autosome,
    }


def materialized_paths(context: dict[str, Any]) -> dict[str, Path]:
    root = context["root"]
    base = root / PLACO_WORK_ROOT / context["fingerprint"] / context["pair_id"]
    return {
        "base": base, "aligned": base / "aligned.tsv.gz", "provenance": base / "materialization.provenance.json",
        "task": base / "task.tsv", "database": base / "alignment.sqlite",
        "checkpoint_dir": base / "checkpoints", "staged_ledger": base / "staged.full.tsv.gz",
        "run_summary": base / "run.summary.tsv", "staged_provenance": base / "staged.provenance.json",
        "canonical_ledger": root / CANONICAL_ROOT / f"{context['pair_id']}.full.tsv.gz",
        "canonical_provenance": root / CANONICAL_ROOT / f"{context['pair_id']}.provenance.json",
        "benchmark_raw": base / "benchmark.raw.tsv", "benchmark": base / "benchmark.tsv",
    }


def expected_task(
    context: dict[str, Any], provenance: dict[str, Any], provenance_sha256: str,
) -> dict[str, object]:
    root, policy = context["root"], context["policy"]
    paths = materialized_paths(context)
    aligned = provenance["aligned_output"]
    return {
        "analysis_id": policy["analysis_id"], "pair_id": context["pair_id"],
        "trait1": context["trait1"], "trait2": context["trait2"], "family_role": context["family_role"],
        "run_fingerprint": context["fingerprint"], "aligned_input": relative(root, paths["aligned"]),
        "aligned_input_sha256": aligned["sha256"], "aligned_input_rows": aligned["rows"],
        "aligned_input_schema": ",".join(ALIGNED_FIELDS),
        "materialization_provenance": relative(root, paths["provenance"]),
        "materialization_provenance_sha256": provenance_sha256,
        "policy": str(POLICY), "policy_sha256": sha256(root / POLICY),
        "contract_lock": str(CONTRACT_LOCK), "contract_lock_sha256": sha256(root / CONTRACT_LOCK),
        "input_gate_lock": str(INPUT_LOCK), "input_gate_lock_sha256": sha256(root / INPUT_LOCK),
        "placo_source": relative(root, context["placo_source"]),
        "placo_source_sha256": policy["placo_plus"]["source_sha256"],
        "marginal_p_threshold": policy["placo_plus"]["marginal_p_threshold_for_nuisance_estimation"],
        "z_squared_maximum": policy["placo_plus"]["z_squared_maximum"],
        "absolute_tolerance": policy["placo_plus"]["absolute_tolerance"],
        "maximum_numerical_failure_fraction": policy["placo_plus"]["maximum_numerical_failure_fraction"],
        "primary_headline_threshold": policy["placo_plus"]["primary_pair_family_headline_threshold"],
        "pairwise_gws_threshold": (
            policy["placo_plus"]["control_genome_wide_threshold"] if context["pair_id"] == "CONTROL"
            else policy["placo_plus"]["per_pair_genome_wide_threshold"]
        ),
        "within_pair_bh_alpha": policy["placo_plus"]["within_pair_bh_alpha"],
        "minimum_free_bytes": MINIMUM_FREE_BYTES,
        "maximum_workers": MAX_WORKERS, "pair_concurrency_limit": 1, "shard_size": SHARD_SIZE,
        "checkpoint_dir": relative(root, paths["checkpoint_dir"]),
        "staged_ledger": relative(root, paths["staged_ledger"]),
        "run_summary": relative(root, paths["run_summary"]),
        "staged_provenance": relative(root, paths["staged_provenance"]),
        "canonical_ledger": relative(root, paths["canonical_ledger"]),
        "canonical_provenance": relative(root, paths["canonical_provenance"]),
        "benchmark_raw": relative(root, paths["benchmark_raw"]),
        "benchmark_out": relative(root, paths["benchmark"]),
        "materializer_sha256": sha256(Path(__file__)), "runner_sha256": sha256(root / RUNNER),
    }


def verify_materialized(context: dict[str, Any]) -> dict[str, Any]:
    paths = materialized_paths(context)
    provenance = read_json(paths["provenance"])
    source_records = provenance.get("source_files")
    expected_sources = context["sources"]
    source_identity_ok = (
        isinstance(source_records, list) and len(source_records) == 2
        and all(isinstance(record, dict) for record in source_records)
        and all(
            record.get("trait_id") == source["trait_id"]
            and record.get("path") == relative(context["root"], source["absolute_path"])
            and record.get("bytes") == int(source["bytes"])
            and record.get("rows") == int(source["rows"])
            and record.get("sha256") == source["sha256"]
            and record.get("schema") == ",".join(SOURCE_FIELDS)
            for record, source in zip(source_records, expected_sources, strict=True)
        )
    )
    counts = provenance.get("alignment_counts", {})
    per_autosome = provenance.get("per_autosome")
    floors = context["policy"]["dense_input_contract"]["pair_materialization_plausibility"]
    expected_coverage = (
        int(counts.get("eligible_written", 0))
        / min(int(expected_sources[0]["rows"]), int(expected_sources[1]["rows"]))
    )
    if (
        provenance.get("schema_version") != "sleep-atlas-track-b-placo-materialization.1"
        or provenance.get("analysis_id") != context["policy"]["analysis_id"]
        or provenance.get("pair_id") != context["pair_id"]
        or (provenance.get("trait1"), provenance.get("trait2"), provenance.get("family_role"))
        != (context["trait1"], context["trait2"], context["family_role"])
        or provenance.get("run_fingerprint") != context["fingerprint"]
        or provenance.get("policy_sha256") != sha256(context["root"] / POLICY)
        or provenance.get("contract_lock_sha256") != sha256(context["root"] / CONTRACT_LOCK)
        or provenance.get("input_gate_lock_sha256") != sha256(context["root"] / INPUT_LOCK)
        or provenance.get("materializer_sha256") != sha256(Path(__file__))
        or provenance.get("runner_sha256") != sha256(context["root"] / RUNNER)
        or provenance.get("placo_source_sha256") != context["policy"]["placo_plus"]["source_sha256"]
        or not source_identity_ok
        or set(counts) != set(REQUIRED_ALIGNMENT_COUNTS)
        or provenance.get("plausibility_floors") != floors
        or not isinstance(per_autosome, list) or len(per_autosome) != 22
        or not all(isinstance(row, dict) for row in per_autosome)
        or [row.get("CHR") for row in per_autosome] != list(range(1, 23))
        or any(int(row.get("eligible_written", 0)) <= 0 for row in per_autosome)
        or not math.isclose(
            float(provenance.get("coverage_fraction_of_smaller_frozen_dense_input", -1)),
            expected_coverage, rel_tol=1e-12, abs_tol=1e-15,
        )
        or int(counts.get("eligible_written", 0)) < int(floors["minimum_aligned_eligible_variants"])
        or expected_coverage < float(floors["minimum_fraction_of_smaller_dense_input"])
        or provenance.get("qc_status") != "PASS"
        or provenance.get("aligned_output", {}).get("sha256") != sha256(paths["aligned"])
        or provenance.get("aligned_output", {}).get("bytes") != paths["aligned"].stat().st_size
        or provenance.get("aligned_output", {}).get("rows") != provenance.get("alignment_counts", {}).get("eligible_written")
        or provenance.get("aligned_output", {}).get("schema") != ",".join(ALIGNED_FIELDS)
        or provenance.get("aligned_output", {}).get("path") != relative(context["root"], paths["aligned"])
    ):
        fail("materialized Track B PLACO+ pair differs from its immutable provenance")
    task_fields, task_rows = read_tsv(paths["task"])
    if len(task_rows) != 1 or task_fields != TASK_FIELDS:
        fail("Track B PLACO+ task schema differs from the materialized pair identity")
    task = task_rows[0]
    wanted_task = expected_task(context, provenance, sha256(paths["provenance"]))
    if any(task.get(key) != str(value) for key, value in wanted_task.items()):
        fail("Track B PLACO+ task differs from the materialized pair identity")
    return provenance


def materialize_pair(context: dict[str, Any]) -> dict[str, Any]:
    root, policy = context["root"], context["policy"]
    paths = materialized_paths(context)
    final_family = [paths["aligned"], paths["provenance"], paths["task"]]
    if all(path.exists() for path in final_family):
        return verify_materialized(context)
    if any(path.exists() for path in final_family):
        fail("partial immutable Track B PLACO+ materialization family already exists")
    paths["base"].mkdir(parents=True, exist_ok=True)
    database = paths["database"]
    aligned_temp = paths["base"] / f".aligned.{os.getpid()}.tmp.gz"
    provenance_temp = paths["base"] / f".provenance.{os.getpid()}.tmp.json"
    task_temp = paths["base"] / f".task.{os.getpid()}.tmp.tsv"
    for path in (database, aligned_temp, provenance_temp, task_temp):
        path.unlink(missing_ok=True)
    connection: sqlite3.Connection | None = None
    try:
        connection = sqlite3.connect(database)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=FILE")
        source_stats = [
            load_source(connection, name, source["absolute_path"], source)
            for name, source in zip(("first", "second"), context["sources"], strict=True)
        ]
        rsid_matches = int(connection.execute("SELECT COUNT(*) FROM first JOIN second USING (snp)").fetchone()[0])
        exact_matches = int(connection.execute(
            "SELECT COUNT(*) FROM first JOIN second USING (snp) WHERE first.chr=second.chr AND first.bp=second.bp"
        ).fetchone()[0])
        counts = {
            "trait1_rows": source_stats[0]["trait_rows"], "trait2_rows": source_stats[1]["trait_rows"],
            "exact_rsid_coordinate_matches": exact_matches, "allele_matches": 0, "allele_flips": 0,
            "allele_mismatches_dropped": 0,
            "duplicate_ids_dropped": source_stats[0]["duplicate_ids_dropped"] + source_stats[1]["duplicate_ids_dropped"],
            "invalid_z_or_p_dropped": source_stats[0]["invalid_z_or_p_dropped"] + source_stats[1]["invalid_z_or_p_dropped"],
            "z_squared_exclusions": 0, "eligible_written": 0,
        }
        extras = {
            "exact_rsid_matches": rsid_matches,
            "coordinate_mismatches_dropped": rsid_matches - exact_matches,
            "trait1_valid_unique_unmatched": source_stats[0]["valid_unique_rows"] - rsid_matches,
            "trait2_valid_unique_unmatched": source_stats[1]["valid_unique_rows"] - rsid_matches,
            "invalid_identity_or_frequency_dropped": (
                source_stats[0]["invalid_identity_or_frequency_dropped"]
                + source_stats[1]["invalid_identity_or_frequency_dropped"]
            ),
            "direct_allele_matches": 0,
        }
        per_chromosome = {
            chromosome: {
                "CHR": chromosome,
                "trait1_rows": source_stats[0]["per_autosome_rows"][chromosome],
                "trait2_rows": source_stats[1]["per_autosome_rows"][chromosome],
                "exact_rsid_coordinate_matches": 0, "allele_matches": 0, "eligible_written": 0,
            }
            for chromosome in range(1, 23)
        }
        query = (
            "SELECT first.snp, first.chr, first.bp, first.a1, first.a2, first.frq, first.z, first.p, "
            "second.a1, second.a2, second.frq, second.z, second.p FROM first JOIN second USING (snp) "
            "WHERE first.chr=second.chr AND first.bp=second.bp ORDER BY first.chr, first.bp, first.snp"
        )
        z_max = float(policy["placo_plus"]["z_squared_maximum"])
        with aligned_temp.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                    writer.writerow(ALIGNED_FIELDS)
                    for values in connection.execute(query):
                        snp, chromosome, bp, a1, a2, frq1, z1, p1, b1, b2, frq2, z2, p2 = values
                        per_chromosome[chromosome]["exact_rsid_coordinate_matches"] += 1
                        if (b1, b2) == (a1, a2):
                            aligned_z2, aligned_frq2 = z2, frq2
                            extras["direct_allele_matches"] += 1
                        elif (b1, b2) == (a2, a1):
                            aligned_z2, aligned_frq2 = -z2, 1.0 - frq2
                            counts["allele_flips"] += 1
                        else:
                            counts["allele_mismatches_dropped"] += 1
                            continue
                        counts["allele_matches"] += 1
                        per_chromosome[chromosome]["allele_matches"] += 1
                        if z1 * z1 > z_max or aligned_z2 * aligned_z2 > z_max:
                            counts["z_squared_exclusions"] += 1
                            continue
                        writer.writerow([
                            snp, chromosome, bp, a1, a2, f"{frq1:.15g}", f"{aligned_frq2:.15g}",
                            f"{z1:.15g}", f"{aligned_z2:.15g}", f"{p1:.15g}", f"{p2:.15g}",
                        ])
                        counts["eligible_written"] += 1
                        per_chromosome[chromosome]["eligible_written"] += 1
            raw.flush()
            os.fsync(raw.fileno())
        if (
            exact_matches != counts["allele_matches"] + counts["allele_mismatches_dropped"]
            or counts["allele_matches"] != counts["z_squared_exclusions"] + counts["eligible_written"]
            or extras["direct_allele_matches"] + counts["allele_flips"] != counts["allele_matches"]
            or any(value < 0 for value in (
                extras["coordinate_mismatches_dropped"], extras["trait1_valid_unique_unmatched"],
                extras["trait2_valid_unique_unmatched"],
            ))
        ):
            fail("Track B PLACO+ exact alignment rejection accounting does not reconcile")
        plausibility = policy["dense_input_contract"]["pair_materialization_plausibility"]
        smaller_frozen = min(counts["trait1_rows"], counts["trait2_rows"])
        coverage_fraction = counts["eligible_written"] / smaller_frozen if smaller_frozen else 0.0
        missing_autosomes = [
            chromosome for chromosome in plausibility["required_autosomes"]
            if per_chromosome[int(chromosome)]["eligible_written"] <= 0
        ]
        if (
            counts["eligible_written"] < int(plausibility["minimum_aligned_eligible_variants"])
            or coverage_fraction < float(plausibility["minimum_fraction_of_smaller_dense_input"])
            or missing_autosomes
        ):
            fail(
                "FAILED_QC Track B PLACO+ alignment plausibility floor failed: "
                f"eligible={counts['eligible_written']} coverage={coverage_fraction:.12g} "
                f"missing_autosomes={missing_autosomes}"
            )
        aligned_size, aligned_hash = file_identity(aligned_temp)
        provenance = {
            "schema_version": "sleep-atlas-track-b-placo-materialization.1",
            "analysis_id": policy["analysis_id"], "pair_id": context["pair_id"],
            "trait1": context["trait1"], "trait2": context["trait2"], "family_role": context["family_role"],
            "run_fingerprint": context["fingerprint"], "policy_sha256": sha256(root / POLICY),
            "contract_lock_sha256": sha256(root / CONTRACT_LOCK),
            "input_gate_lock_sha256": sha256(root / INPUT_LOCK),
            "materializer_sha256": sha256(Path(__file__)), "runner_sha256": sha256(root / RUNNER),
            "placo_source_sha256": policy["placo_plus"]["source_sha256"],
            "source_files": [
                {
                    "trait_id": source["trait_id"], "path": relative(root, source["absolute_path"]),
                    "bytes": source["bytes"], "rows": source["rows"], "sha256": source["sha256"],
                    "schema": source["schema"], "scan": source_stats[index],
                }
                for index, source in enumerate(context["sources"])
            ],
            "alignment_rule": policy["dense_input_contract"]["alignment_rule"],
            "alignment_counts": counts, "additional_rejection_counts": extras,
            "allele_match_count_semantics": (
                "allele_matches counts every exact unordered compatible pair; allele_flips is its swapped subset; "
                "direct_allele_matches is recorded in additional_rejection_counts"
            ),
            "coverage_fraction_of_smaller_frozen_dense_input": coverage_fraction,
            "per_autosome": [per_chromosome[chromosome] for chromosome in range(1, 23)],
            "plausibility_floors": plausibility,
            "aligned_output": {
                "path": relative(root, paths["aligned"]), "bytes": aligned_size,
                "rows": counts["eligible_written"], "sha256": aligned_hash,
                "schema": ",".join(ALIGNED_FIELDS),
            },
            "qc_status": "PASS",
        }
        provenance_temp.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        task = expected_task(context, provenance, sha256(provenance_temp))
        if list(task) != TASK_FIELDS:
            fail("internal PLACO+ task schema construction drifted")
        task_temp.write_bytes(tsv_bytes(TASK_FIELDS, [task]))
        exclusive_family(
            [aligned_temp, provenance_temp, task_temp],
            [paths["aligned"], paths["provenance"], paths["task"]],
        )
        return provenance
    finally:
        if connection is not None:
            connection.close()
        for path in (database, aligned_temp, provenance_temp, task_temp):
            path.unlink(missing_ok=True)


def validate_task_against_context(task: dict[str, str], context: dict[str, Any], task_path: Path) -> None:
    paths = materialized_paths(context)
    provenance = verify_materialized(context)
    fields, rows = read_tsv(task_path)
    expected = expected_task(context, provenance, sha256(paths["provenance"]))
    if (
        task_path.resolve() != paths["task"].resolve()
        or fields != TASK_FIELDS or len(rows) != 1 or rows[0] != task
        or any(task.get(key) != str(value) for key, value in expected.items())
    ):
        fail(f"Track B PLACO+ task differs from production context: {task_path}")


def validate_staged_ledger(root: Path, task: dict[str, str], summary: dict[str, str]) -> dict[str, Any]:
    ledger = safe_relative(root, task["staged_ledger"])
    aligned = safe_relative(root, task["aligned_input"])
    expected_rows = int(task["aligned_input_rows"])
    family_n = int(summary["within_pair_family_n"])
    identity = PAIR_IDENTITIES.get(task.get("pair_id", ""))
    if identity != (task.get("trait1"), task.get("trait2"), task.get("family_role")):
        fail("PLACO+ staged family is not exactly frozen pair A, B, or CONTROL")
    if (
        expected_rows <= 0 or family_n != expected_rows or int(summary["aligned_rows"]) != expected_rows
        or task.get("aligned_input_schema") != ",".join(ALIGNED_FIELDS)
        or task.get("aligned_input_sha256") != sha256(aligned)
    ):
        fail("PLACO+ summary family size differs from the complete aligned family")
    if summary.get("ledger_sha256") != sha256(ledger) or int(summary.get("ledger_bytes", -1)) != ledger.stat().st_size:
        fail("staged PLACO+ ledger differs from its run summary")
    expected_summary = {
        "analysis_id": task["analysis_id"], "pair_id": task["pair_id"],
        "run_fingerprint": task["run_fingerprint"], "ledger": task["staged_ledger"],
        "ledger_schema": ",".join(LEDGER_FIELDS), "aligned_sha256": task["aligned_input_sha256"],
        "placo_source_sha256": task["placo_source_sha256"],
        "materializer_sha256": task["materializer_sha256"], "runner_sha256": task["runner_sha256"],
    }
    if any(summary.get(key) != value for key, value in expected_summary.items()):
        fail("PLACO+ run summary identity, schema, or software provenance drifted")
    if any(not summary.get(key, "").strip() for key in ("r_version", "data_table_version", "runtime_platform")):
        fail("PLACO+ run summary lacks R runtime provenance")
    if not SHA256_RE.fullmatch(summary.get("nuisance_checkpoint_sha256", "")):
        fail("PLACO+ run summary lacks its frozen global nuisance-checkpoint hash")
    nuisance_checkpoint = safe_relative(root, task["checkpoint_dir"]) / "nuisance.rds"
    if summary["nuisance_checkpoint_sha256"] != sha256(nuisance_checkpoint):
        fail("PLACO+ frozen global nuisance checkpoint drifted before publication")
    try:
        variance_rows = int(summary["nuisance_null_rows_variance"])
        correlation_rows = int(summary["nuisance_null_rows_correlation"])
        variance1, variance2 = float(summary["VarZ1"]), float(summary["VarZ2"])
        correlation = float(summary["CorZ"])
        shard_count = int(summary["shard_count"])
        shard_size = int(summary["shard_size"])
        workers = int(summary["workers"])
        elapsed = float(summary["elapsed_seconds"])
    except ValueError:
        fail("PLACO+ run summary has invalid nuisance or execution provenance")
    if (
        not 2 <= variance_rows <= family_n or not 2 <= correlation_rows <= family_n
        or not math.isfinite(variance1) or variance1 <= 0
        or not math.isfinite(variance2) or variance2 <= 0
        or not math.isfinite(correlation) or not -1 < correlation < 1
        or shard_size != int(task["shard_size"])
        or shard_count != math.ceil(family_n / shard_size)
        or not 1 <= workers <= int(task["maximum_workers"])
        or not math.isfinite(elapsed) or elapsed < 0
    ):
        fail("PLACO+ run summary has invalid nuisance or execution provenance")
    database = ledger.with_name(f".{ledger.name}.{os.getpid()}.validate.sqlite")
    database.unlink(missing_ok=True)
    connection = sqlite3.connect(database)
    counts = {"rows": 0, "failures": 0, "primary": 0, "pairwise": 0, "bh": 0}
    minimum_p = 1.0
    try:
        connection.execute("CREATE TABLE pvalues (row_index INTEGER PRIMARY KEY, p REAL NOT NULL, q REAL NOT NULL)")
        with (
            gzip.open(aligned, "rt", encoding="utf-8", newline="") as aligned_handle,
            gzip.open(ledger, "rt", encoding="utf-8", newline="") as handle,
        ):
            aligned_reader = csv.DictReader(aligned_handle, delimiter="\t")
            reader = csv.DictReader(handle, delimiter="\t")
            if list(aligned_reader.fieldnames or []) != ALIGNED_FIELDS:
                fail("materialized PLACO+ aligned family schema drifted")
            if list(reader.fieldnames or []) != LEDGER_FIELDS:
                fail("staged PLACO+ ledger schema differs from policy")
            batch = []
            paired_rows = itertools.zip_longest(reader, aligned_reader)
            for index, (row, aligned_row) in enumerate(paired_rows, start=1):
                if row is None or aligned_row is None:
                    fail("staged PLACO+ ledger does not retain every aligned eligible variant exactly once")
                if (
                    None in row or None in aligned_row
                    or any(value is None for value in row.values())
                    or any(value is None for value in aligned_row.values())
                ):
                    fail("staged PLACO+ ledger or aligned family contains a ragged row")
                if row["analysis_id"] != task["analysis_id"] or row["pair_id"] != task["pair_id"]:
                    fail("staged PLACO+ ledger row identity drifted")
                try:
                    p_value = float(row["P_PLACO_PLUS"])
                    q_value = float(row["PLACO_BH_Q"])
                    observed_family = int(row["within_pair_family_n"])
                    integer_identity_ok = all(
                        int(row[field]) == int(aligned_row[field]) for field in ("CHR", "BP")
                    )
                    numeric_identity_ok = all(
                        float(row[field]) == float(aligned_row[field])
                        for field in ("Z1", "Z2", "P1", "P2")
                    )
                except (TypeError, ValueError):
                    fail("staged PLACO+ ledger contains invalid P/q/family provenance")
                if (
                    not math.isfinite(p_value) or not 0 <= p_value <= 1
                    or not math.isfinite(q_value) or not 0 <= q_value <= 1
                    or observed_family != family_n
                    or not integer_identity_ok or not numeric_identity_ok
                    or any(row[field] != aligned_row[field] for field in ("SNP", "A1", "A2"))
                ):
                    fail("staged PLACO+ ledger does not exactly preserve its aligned hypothesis family")
                if row["analysis_status"] == "NUMERICAL_FAILURE_P_SET_TO_ONE":
                    if p_value != 1.0 or row["T_PLACO_PLUS"] not in {"", "NA", "NaN", "nan"} or not row["numerical_error"].strip() or row["numerical_error"] == "NA":
                        fail("numerical failure row does not retain P=1/error semantics")
                    counts["failures"] += 1
                elif row["analysis_status"] == "TESTED":
                    try:
                        statistic = float(row["T_PLACO_PLUS"])
                    except ValueError:
                        fail("tested PLACO+ row lacks a finite statistic")
                    if not math.isfinite(statistic) or row["numerical_error"] not in {"", "NA"}:
                        fail("tested PLACO+ row has invalid statistic/error provenance")
                else:
                    fail("invalid PLACO+ terminal row status")
                minimum_p = min(minimum_p, p_value)
                counts["primary"] += task["pair_id"] != "CONTROL" and p_value <= float(task["primary_headline_threshold"])
                counts["pairwise"] += p_value <= float(task["pairwise_gws_threshold"])
                counts["bh"] += q_value <= float(task["within_pair_bh_alpha"])
                batch.append((index, p_value, q_value))
                if len(batch) == 50_000:
                    connection.executemany("INSERT INTO pvalues VALUES (?, ?, ?)", batch)
                    batch.clear()
            if batch:
                connection.executemany("INSERT INTO pvalues VALUES (?, ?, ?)", batch)
        connection.commit()
        counts["rows"] = int(connection.execute("SELECT COUNT(*) FROM pvalues").fetchone()[0])
        if counts["rows"] != family_n:
            fail("staged PLACO+ ledger is not the complete aligned hypothesis family")
        rank = family_n
        running = 1.0
        cursor = connection.execute("SELECT p, q FROM pvalues ORDER BY p DESC, row_index DESC")
        pending: list[tuple[float, float]] = []
        pending_p: float | None = None

        def check_group(group: list[tuple[float, float]], p_value: float, high_rank: int, running_q: float) -> float:
            wanted = min(running_q, min(1.0, p_value * family_n / high_rank))
            for _, observed_q in group:
                if not math.isclose(observed_q, wanted, rel_tol=1e-12, abs_tol=1e-300):
                    fail("PLACO_BH_Q differs from exact BH over the fixed within-pair family")
            return wanted

        high_rank = rank
        for p_value, q_value in cursor:
            if pending_p is None or p_value == pending_p:
                if pending_p is None:
                    pending_p, high_rank = p_value, rank
                pending.append((p_value, q_value))
            else:
                running = check_group(pending, pending_p, high_rank, running)
                rank -= len(pending)
                pending = [(p_value, q_value)]
                pending_p, high_rank = p_value, rank
        if pending:
            running = check_group(pending, float(pending_p), high_rank, running)
            rank -= len(pending)
        if rank != 0:
            fail("internal BH validation rank mismatch")
    finally:
        connection.close()
        database.unlink(missing_ok=True)
    expected_failure_fraction = counts["failures"] / family_n
    summary_checks = {
        "numerical_failure_count": counts["failures"], "primary_headline_count": counts["primary"],
        "pairwise_gws_count": counts["pairwise"], "within_pair_bh_count": counts["bh"],
    }
    for key, wanted in summary_checks.items():
        if int(summary.get(key, -1)) != wanted:
            fail(f"PLACO+ run summary count drifted: {key}")
    if not math.isclose(float(summary["numerical_failure_fraction"]), expected_failure_fraction, rel_tol=1e-10, abs_tol=1e-12):
        fail("PLACO+ run summary failure fraction drifted")
    if expected_failure_fraction > float(task["maximum_numerical_failure_fraction"]):
        fail("PLACO+ numerical failure fraction exceeds policy")
    if not math.isclose(float(summary["minimum_p_placo_plus"]), minimum_p, rel_tol=1e-12, abs_tol=1e-300):
        fail("PLACO+ run summary minimum P drifted")
    terminal = "COMPLETE_WITH_HITS" if counts["primary"] or counts["pairwise"] or counts["bh"] else "TESTED_NO_HIT"
    if summary.get("terminal_status") != terminal:
        fail("PLACO+ terminal zero-hit/hit status differs from the complete ledger")
    return {**counts, "minimum_p": minimum_p, "terminal_status": terminal}


def publish_run(root: Path, task_path: Path, summary_path: Path) -> dict[str, Any]:
    task_fields, task_rows = read_tsv(task_path)
    summary_fields, summary_rows = read_tsv(summary_path)
    if task_fields != TASK_FIELDS or len(task_rows) != 1:
        fail("PLACO+ publication task schema drifted")
    if summary_fields != SUMMARY_FIELDS or len(summary_rows) != 1:
        fail("PLACO+ publication summary schema drifted")
    task, summary = task_rows[0], summary_rows[0]
    if summary.get("task_sha256") != sha256(task_path):
        fail("PLACO+ run summary differs from its immutable task")
    if summary.get("run_fingerprint") != task.get("run_fingerprint") or summary.get("pair_id") != task.get("pair_id"):
        fail("PLACO+ run summary identity differs from its task")
    locked_artifacts = {
        "policy_sha256": task["policy"], "contract_lock_sha256": task["contract_lock"],
        "input_gate_lock_sha256": task["input_gate_lock"],
        "materialization_provenance_sha256": task["materialization_provenance"],
        "placo_source_sha256": task["placo_source"],
        "materializer_sha256": str(MATERIALIZER), "runner_sha256": str(RUNNER),
    }
    for hash_field, artifact in locked_artifacts.items():
        if task.get(hash_field) != sha256(safe_relative(root, artifact)):
            fail(f"PLACO+ locked publication artifact drifted: {artifact}")
    counts = validate_staged_ledger(root, task, summary)
    ledger = safe_relative(root, task["staged_ledger"])
    staged_provenance = safe_relative(root, task["staged_provenance"])
    canonical_ledger = safe_relative(root, task["canonical_ledger"])
    canonical_provenance = safe_relative(root, task["canonical_provenance"])
    if canonical_ledger.exists() or canonical_provenance.exists() or staged_provenance.exists():
        fail("canonical or staged Track B PLACO+ result already exists; overwrite is forbidden")
    output_size, output_hash = file_identity(ledger)
    provenance = {
        "schema_version": "sleep-atlas-track-b-placo-result.1",
        "analysis_id": task["analysis_id"], "pair_id": task["pair_id"],
        "trait1": task["trait1"], "trait2": task["trait2"], "family_role": task["family_role"],
        "run_fingerprint": task["run_fingerprint"], "policy_sha256": task["policy_sha256"],
        "contract_lock_sha256": task["contract_lock_sha256"],
        "input_gate_lock_sha256": task["input_gate_lock_sha256"],
        "materialization_provenance_sha256": task["materialization_provenance_sha256"],
        "task_sha256": summary["task_sha256"],
        "software_sha256_or_commit": task["placo_source_sha256"],
        "materializer_sha256": task["materializer_sha256"], "runner_sha256": task["runner_sha256"],
        "reference_sha256": "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING",
        "input_sha256": task["aligned_input_sha256"], "input_rows": int(task["aligned_input_rows"]),
        "input_schema": task["aligned_input_schema"], "output": task["canonical_ledger"],
        "output_sha256": output_hash, "output_bytes": output_size, "output_rows": counts["rows"],
        "exact_schema": ",".join(LEDGER_FIELDS), "complete_family_counts": counts,
        "nuisance_estimation": {
            "method": "PINNED_PLACO_PLUS_var.placo_AND_cor.pearson",
            "marginal_p_threshold": float(task["marginal_p_threshold"]),
            "VarZ1": float(summary["VarZ1"]), "VarZ2": float(summary["VarZ2"]),
            "CorZ": float(summary["CorZ"]),
            "variance_null_rows": int(summary["nuisance_null_rows_variance"]),
            "correlation_null_rows": int(summary["nuisance_null_rows_correlation"]),
            "input_rows": counts["rows"],
            "scope": "ALL_VALID_ALIGNED_GENOME_WIDE_VARIANTS_BEFORE_ANY_SINGLE_VARIANT_SHARD",
            "checkpoint_sha256": summary["nuisance_checkpoint_sha256"],
        },
        "numerical_failure_policy": (
            "RETAIN_EVERY_INDIVIDUAL_ERROR_WARNING_NONFINITE_OR_OUT_OF_RANGE_ROW_AND_SET_P_PLACO_PLUS_TO_1"
        ),
        "within_pair_bh_family_n": counts["rows"], "terminal_result_state": counts["terminal_status"],
        "execution_resource_envelope": {
            "pair_concurrency_limit": int(task["pair_concurrency_limit"]),
            "maximum_workers": int(task["maximum_workers"]), "observed_workers": int(summary["workers"]),
            "shard_size": int(task["shard_size"]),
            "ram_gate": "NO_ASSUMED_WHOLE_PIPELINE_GATE;MEASURE_GLOBAL_NUISANCE_WITH_BENCHMARK_HOOK",
        },
        "qc_status": "PASS", "run_summary": {key: summary[key] for key in SUMMARY_FIELDS},
    }
    staged_provenance.parent.mkdir(parents=True, exist_ok=True)
    temporary = staged_provenance.with_name(f".{staged_provenance.name}.{os.getpid()}.tmp")
    temporary.unlink(missing_ok=True)
    try:
        temporary.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if file_identity(ledger) != (output_size, output_hash):
            fail("staged PLACO+ ledger changed after semantic validation")
        exclusive_family(
            [ledger, temporary, temporary],
            [canonical_ledger, canonical_provenance, staged_provenance],
        )
        if file_identity(canonical_ledger) != (output_size, output_hash):
            fail("canonical PLACO+ ledger differs from the validated stage")
        return provenance
    finally:
        temporary.unlink(missing_ok=True)


def resource_preflight(context: dict[str, Any]) -> dict[str, int]:
    memory = physical_memory_bytes()
    source_rows = sum(int(source["rows"]) for source in context["sources"])
    estimated_work = MINIMUM_FREE_BYTES + 600 * source_rows
    free = shutil.disk_usage(context["root"]).free
    if free < estimated_work:
        fail(f"BLOCKED_BY_COMPUTE PLACO+ requires {estimated_work} free bytes for this pair; observed {free}")
    return {"physical_memory_bytes": memory, "free_bytes": free, "estimated_work_bytes": estimated_work}


def process_tree_rss_bytes(root_pid: int) -> int:
    try:
        result = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="], check=False, capture_output=True, text=True,
        )
    except OSError:
        return 0
    if result.returncode != 0:
        return 0
    records: dict[int, tuple[int, int]] = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            pid, parent, rss_kib = map(int, fields)
        except ValueError:
            continue
        records[pid] = (parent, rss_kib)
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, (parent, _) in records.items():
            if parent in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return 1024 * sum(records.get(pid, (0, 0))[1] for pid in descendants)


def benchmark_pair(context: dict[str, Any], variants: int, workers: int) -> dict[str, str]:
    if variants < 1 or workers < 1 or workers > MAX_WORKERS:
        fail(f"benchmark requires variants>0 and workers between 1 and {MAX_WORKERS}")
    paths = materialized_paths(context)
    verify_materialized(context)
    if paths["benchmark_raw"].exists() or paths["benchmark"].exists():
        fail("immutable Track B PLACO+ benchmark already exists")
    command = [
        str(context["root"] / ".r-env/bin/Rscript"), str(context["root"] / RUNNER),
        relative(context["root"], paths["task"]), "--root", str(context["root"]),
        "--execute", "--workers", str(workers), "--benchmark-variants", str(variants),
    ]
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=context["root"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak_rss = 0
    while process.poll() is None:
        peak_rss = max(peak_rss, process_tree_rss_bytes(process.pid))
        time.sleep(0.05)
    stdout, stderr = process.communicate()
    wall_seconds = time.monotonic() - started
    peak_rss = max(peak_rss, process_tree_rss_bytes(process.pid))
    if process.returncode != 0:
        fail("PLACO+ benchmark runner failed: " + (stdout + "\n" + stderr).strip())
    raw_fields, raw_rows = read_tsv(paths["benchmark_raw"])
    if raw_fields != BENCHMARK_RAW_FIELDS or len(raw_rows) != 1:
        fail("PLACO+ raw benchmark schema drifted")
    raw = raw_rows[0]
    task = one_tsv_row(paths["task"])
    nuisance_checkpoint = paths["checkpoint_dir"] / "nuisance.rds"
    expected_raw = {
        "analysis_id": task["analysis_id"], "pair_id": task["pair_id"],
        "run_fingerprint": task["run_fingerprint"], "input_rows": task["aligned_input_rows"],
        "input_sha256": task["aligned_input_sha256"], "task_sha256": sha256(paths["task"]),
        "policy_sha256": task["policy_sha256"], "contract_lock_sha256": task["contract_lock_sha256"],
        "input_gate_lock_sha256": task["input_gate_lock_sha256"],
        "materializer_sha256": task["materializer_sha256"],
        "placo_source_sha256": task["placo_source_sha256"],
        "nuisance_checkpoint_sha256": sha256(nuisance_checkpoint),
        "global_nuisance_input_rows": task["aligned_input_rows"],
        "benchmark_variants": str(min(variants, int(task["aligned_input_rows"]))), "workers": str(workers),
        "runner_sha256": task["runner_sha256"],
        "scientific_equivalence": (
            "GLOBAL_OFFICIAL_NUISANCE_ESTIMATED_ON_ALL_VALID_VARIANTS_"
            "THEN_IDENTICAL_SINGLE_VARIANT_PLACO_PLUS_CALLS_SAMPLED"
        ),
    }
    if any(raw.get(key) != value for key, value in expected_raw.items()):
        fail("PLACO+ raw benchmark identity, hash, or scientific-equivalence provenance drifted")
    try:
        input_rows = int(raw["input_rows"])
        variance_rows = int(raw["nuisance_null_rows_variance"])
        correlation_rows = int(raw["nuisance_null_rows_correlation"])
        failures = int(raw["numerical_failures"])
        benchmark_n = int(raw["benchmark_variants"])
        nuisance_elapsed = float(raw["global_nuisance_elapsed_seconds"])
        test_elapsed = float(raw["single_variant_testing_elapsed_seconds"])
        total_elapsed = float(raw["measured_runner_elapsed_seconds"])
        rate = float(raw["variants_per_second"])
        projected = float(raw["projected_full_family_seconds"])
        variance1, variance2, correlation = float(raw["VarZ1"]), float(raw["VarZ2"]), float(raw["CorZ"])
    except ValueError:
        fail("PLACO+ raw benchmark contains invalid numeric provenance")
    if (
        not 2 <= variance_rows <= input_rows or not 2 <= correlation_rows <= input_rows
        or not 0 <= failures <= benchmark_n
        or not math.isfinite(nuisance_elapsed) or nuisance_elapsed < 0
        or not math.isfinite(test_elapsed) or test_elapsed <= 0
        or not math.isfinite(total_elapsed) or total_elapsed <= 0
        or not math.isfinite(rate) or rate <= 0
        or not math.isclose(projected, input_rows / rate, rel_tol=1e-10, abs_tol=1e-12)
        or not math.isfinite(variance1) or variance1 <= 0
        or not math.isfinite(variance2) or variance2 <= 0
        or not math.isfinite(correlation) or not -1 < correlation < 1
        or raw["global_nuisance_checkpoint_reused"] not in {"TRUE", "FALSE"}
        or any(not raw[key].strip() for key in ("r_version", "data_table_version", "runtime_platform"))
    ):
        fail("PLACO+ raw benchmark contains invalid measured-runtime or nuisance provenance")
    if peak_rss > 0:
        rss_method = "AGGREGATE_PROCESS_TREE_PS_SAMPLED"
    else:
        observed = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        peak_rss = observed if sys.platform == "darwin" else observed * 1024
        rss_method = "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE"
    host_memory = physical_memory_bytes()
    host_free = shutil.disk_usage(context["root"]).free
    fields = list(raw) + [
        "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
        "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
        "peak_rss_fraction_of_physical_memory",
    ]
    row: dict[str, object] = {
        **raw, "peak_process_tree_rss_bytes": peak_rss, "rss_measurement_method": rss_method,
        "wrapper_wall_seconds": f"{wall_seconds:.12g}", "rss_sample_interval_seconds": "0.05",
        "host_physical_memory_bytes": host_memory, "host_free_bytes_at_benchmark": host_free,
        "peak_rss_fraction_of_physical_memory": f"{peak_rss / host_memory:.12g}",
    }
    temporary = paths["benchmark"].with_name(f".{paths['benchmark'].name}.{os.getpid()}.tmp")
    temporary.write_bytes(tsv_bytes(fields, [row]))
    try:
        exclusive_family([temporary], [paths["benchmark"]])
    finally:
        temporary.unlink(missing_ok=True)
    return {key: str(value) for key, value in row.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_id", choices=PAIR_IDS)
    parser.add_argument("--root", default=".")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--materialize", action="store_true")
    action.add_argument("--verify-materialized", action="store_true")
    action.add_argument("--publish-run", action="store_true")
    action.add_argument("--benchmark", action="store_true")
    parser.add_argument("--task")
    parser.add_argument("--summary")
    parser.add_argument("--benchmark-variants", type=int, default=2000)
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    context = production_context(root, args.pair_id)
    if args.preflight:
        resources = resource_preflight(context)
        print(
            f"TRACK_B_PLACO_PREFLIGHT_PASS pair={args.pair_id} fingerprint={context['fingerprint']} "
            f"ram={resources['physical_memory_bytes']} free={resources['free_bytes']}"
        )
    elif args.materialize:
        resource_preflight(context)
        provenance = materialize_pair(context)
        print(
            f"TRACK_B_PLACO_PAIR_MATERIALIZED pair={args.pair_id} "
            f"eligible={provenance['alignment_counts']['eligible_written']} fingerprint={context['fingerprint']}"
        )
        print(f"TASK={relative(root, materialized_paths(context)['task'])}")
    elif args.verify_materialized:
        provenance = verify_materialized(context)
        print(f"TRACK_B_PLACO_PAIR_VERIFIED pair={args.pair_id} eligible={provenance['alignment_counts']['eligible_written']}")
    elif args.benchmark:
        resource_preflight(context)
        row = benchmark_pair(context, args.benchmark_variants, args.workers)
        print(
            f"TRACK_B_PLACO_BENCHMARK_COMPLETE pair={args.pair_id} variants={row['benchmark_variants']} "
            f"workers={row['workers']} peak_rss={row['peak_process_tree_rss_bytes']} "
            f"projected_full_seconds={row['projected_full_family_seconds']}"
        )
    else:
        if not args.task or not args.summary:
            fail("--publish-run requires --task and --summary")
        task_path = safe_relative(root, args.task)
        summary_path = safe_relative(root, args.summary)
        task = one_tsv_row(task_path)
        validate_task_against_context(task, context, task_path)
        provenance = publish_run(root, task_path, summary_path)
        print(
            f"TRACK_B_PLACO_PAIR_PUBLISHED pair={args.pair_id} rows={provenance['output_rows']} "
            f"status={provenance['terminal_result_state']}"
        )


if __name__ == "__main__":
    main()
