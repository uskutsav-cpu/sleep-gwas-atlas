#!/usr/bin/env python3
"""Run Track B LAVA as measured, resumable, whole-locus fresh R processes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from bisect import bisect_right
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "scripts/119_track_b_lava_contract.py"
RUNNER = ROOT / "scripts/120_run_track_b_lava.R"
VALIDATOR = ROOT / "scripts/121_validate_track_b_lava.py"
CHECKPOINT_VALIDATOR = ROOT / "scripts/132_validate_track_b_lava_checkpoint.R"
R_SCRIPT = ROOT / ".r-env/bin/Rscript"
POLICY = ROOT / "config/track_b_local_analysis_policy.json"
INPUT_LOCK = ROOT / "results/track_b/local_analysis_input.lock.json"
REFERENCE_LOCK = ROOT / "ref/lava/ukb_v1.1/reference.provenance.json"
CHROMOSOME_INPUT_LOCK = ROOT / "results/track_b/lava_chromosome_inputs.provenance.json"
LOCUS_FILE = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
REFERENCE_PREFIX = ROOT / "ref/lava/ukb_v1.1/lava-ukb-v1.1"
CHECKPOINT_ROOT = ROOT / "results/track_b/checkpoints/lava"
BENCHMARK = ROOT / "results/track_b/RAM_BENCHMARK.tsv"
BENCHMARK_PROVENANCE = ROOT / "results/track_b/RAM_BENCHMARK.provenance.json"
BENCHMARK_NAMESPACE = ROOT / "results/track_b/RAM_BENCHMARK.namespace.json"
REPORT = ROOT / "results/track_b/RAM_AWARE_EXECUTION_REPORT.md"
RESULT_LOCK = ROOT / "results/track_b/local/lava_results.provenance.json"
EXECUTION_EVIDENCE = ROOT / "results/track_b/local/lava_execution_evidence.provenance.json"

BENCHMARK_FIELDS = [
    "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
    "runtime_sec", "exit_status", "output_hash",
]
SELECTION_FIELDS = [
    "benchmark_phase", "selection_role", "candidate_rank", "locus_index", "locus",
    "chromosome", "start", "stop", "reference_n_snps",
    "reference_chromosome_n_snps", "selection_basis", "execution_fingerprint",
]
CONDITIONAL_CANDIDATE_FIELDS = [
    "locus_index", "locus", "chromosome", "start", "stop", "n_snps",
    "eligible_models", "eligible_pairs",
]
MARKER_PREFIX = "TRACK_B_LAVA_WORKER\t"
SEMANTIC_MARKER_PREFIX = "TRACK_B_LAVA_CHECKPOINT_SEMANTIC\t"
FINGERPRINT_RE = re.compile(r"^[0-9a-f]{64}$")
LOCUS_PHASES = {"discovery", "conditional"}
ALL_PHASES = {"discovery", "aggregate-discovery", "conditional", "finalize"}
ANALYSIS_BY_PHASE = {
    "discovery": "LAVA_DISCOVERY",
    "aggregate-discovery": "LAVA_DISCOVERY_AGGREGATION",
    "conditional": "LAVA_CONDITIONAL",
    "finalize": "LAVA_FINALIZATION",
}
SCIENTIFICALLY_COMPLETE_QC = {
    "discovery": {"PROCESSED"},
    "conditional": {"CONDITIONAL_LOCUS_PROCESSED"},
    "aggregate-discovery": {"FULL_FAMILY_BH_COMPLETE"},
    "finalize": {"STAGED_SEMANTIC_VALIDATION_COMPLETE"},
}


class ExecutionError(RuntimeError):
    pass


def file_identity(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if before.st_size <= 0:
                raise ExecutionError(f"artifact is empty: {path}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        raise ExecutionError(f"could not hash artifact {path}: {error}") from error
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before) != identity(after) or identity(after) != identity(current):
        raise ExecutionError(f"artifact changed while hashing: {path}")
    return after.st_size, digest.hexdigest()


def sha256(path: Path) -> str:
    return file_identity(path)[1]


def stat_fields(value: os.stat_result) -> dict[str, int]:
    return {
        "device": int(value.st_dev), "inode": int(value.st_ino), "bytes": int(value.st_size),
        "mtime_ns": int(value.st_mtime_ns), "ctime_ns": int(value.st_ctime_ns),
    }


def stable_stat(path: Path) -> dict[str, int]:
    try:
        if path.is_symlink():
            raise ExecutionError(f"symbolic-link active inputs are forbidden: {path}")
        before = path.stat()
        if not path.is_file() or before.st_size <= 0:
            raise ExecutionError(f"active input is missing or empty: {path}")
        after = path.stat()
    except OSError as error:
        raise ExecutionError(f"could not stat active input {path}: {error}") from error
    if stat_fields(before) != stat_fields(after):
        raise ExecutionError(f"active input changed while being captured: {path}")
    return stat_fields(after)


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ExecutionError(f"{label} is unreadable: {error}") from error
    if not isinstance(value, dict):
        raise ExecutionError(f"{label} is not a JSON object")
    return value


def physical_memory_bytes() -> int:
    try:
        value = int(os.sysconf("SC_PHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        value = 0
    if value <= 0:
        raise ExecutionError("physical RAM could not be measured")
    return value


def memory_safety_reserve_bytes() -> int:
    try:
        policy = json.loads(POLICY.read_text(encoding="utf-8"))
        reserve = int(policy["ram_aware_execution"]["memory_safety_reserve_bytes"])
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
        raise ExecutionError(f"RAM safety reserve is unreadable: {error}") from error
    if reserve < 1024**3:
        raise ExecutionError("RAM safety reserve must be at least 1 GiB")
    return reserve


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def publish_no_replace(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise ExecutionError(f"refusing to replace existing artifact: {path}") from error
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def tsv_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def load_loci() -> list[dict[str, int]]:
    with LOCUS_FILE.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        if header != ["LOC", "CHR", "START", "STOP"]:
            raise ExecutionError("official LAVA locus header drifted")
        rows = []
        for line in handle:
            values = line.split()
            if len(values) != 4:
                raise ExecutionError("malformed official LAVA locus row")
            row = dict(zip(header, values, strict=True))
            rows.append({key: int(value) for key, value in row.items()})
    if len(rows) != 2495 or len({row["LOC"] for row in rows}) != 2495:
        raise ExecutionError("official LAVA locus family is not the exact 2,495 loci")
    return rows


def reference_snp_burdens(loci: list[dict[str, int]]) -> tuple[list[int], dict[int, int]]:
    """Measure only frozen reference inputs; never inspect local association results."""
    counts = [0] * len(loci)
    chromosome_counts: dict[int, int] = {}
    by_chromosome: dict[int, list[tuple[int, dict[str, int]]]] = {}
    for index, locus in enumerate(loci, start=1):
        by_chromosome.setdefault(locus["CHR"], []).append((index, locus))
    for chromosome in range(1, 23):
        chromosome_loci = by_chromosome.get(chromosome, [])
        starts = [item[1]["START"] for item in chromosome_loci]
        path = Path(f"{REFERENCE_PREFIX}_chr{chromosome}.info")
        if not path.is_file() or path.stat().st_size == 0:
            raise ExecutionError(f"sealed chromosome reference is missing: {path}")
        total = 0
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"CHR", "POS", "SNP"}.issubset(reader.fieldnames or []):
                raise ExecutionError(f"malformed LAVA reference info: {path}")
            for row in reader:
                try:
                    observed_chromosome = int(row["CHR"])
                    position = int(row["POS"])
                except ValueError as error:
                    raise ExecutionError(f"invalid LAVA reference coordinate: {path}") from error
                if observed_chromosome != chromosome:
                    raise ExecutionError(f"wrong chromosome in LAVA reference info: {path}")
                total += 1
                offset = bisect_right(starts, position) - 1
                if offset >= 0:
                    locus_index, locus = chromosome_loci[offset]
                    if position <= locus["STOP"]:
                        counts[locus_index - 1] += 1
        chromosome_counts[chromosome] = total
    return counts, chromosome_counts


def reference_snp_counts(loci: list[dict[str, int]]) -> list[int]:
    return reference_snp_burdens(loci)[0]


def _ranked_candidates(
    loci: list[dict[str, int]], counts: list[int], chromosome_counts: dict[int, int],
    fingerprint: str, *, benchmark_phase: str, eligible_indices: set[int] | None = None,
) -> list[dict[str, object]]:
    eligible = []
    for index, (locus, count) in enumerate(zip(loci, counts, strict=True), start=1):
        if count <= 0 or (eligible_indices is not None and index not in eligible_indices):
            continue
        eligible.append((int(chromosome_counts.get(locus["CHR"], 0)), count, index, locus))
    if not eligible:
        return []
    eligible.sort(key=lambda item: (item[0], item[1], item[2]))
    if len(eligible) >= 3:
        roles_and_targets = [
            ("SMALL", 0), ("MEDIAN", (len(eligible) - 1) // 2), ("LARGE", len(eligible) - 1),
        ]
    elif len(eligible) == 2:
        roles_and_targets = [("SMALL", 0), ("LARGE", 1)]
    else:
        roles_and_targets = [("ONLY", 0)]
    rows: list[dict[str, object]] = []
    for role, target in roles_and_targets:
        ordered_positions = sorted(
            range(len(eligible)),
            key=lambda position: (abs(position - target), position if role != "LARGE" else -position),
        )
        for rank, position in enumerate(ordered_positions, start=1):
            chromosome_count, count, index, locus = eligible[position]
            rows.append({
                "benchmark_phase": benchmark_phase,
                "selection_role": role,
                "candidate_rank": rank,
                "locus_index": index,
                "locus": locus["LOC"],
                "chromosome": locus["CHR"],
                "start": locus["START"],
                "stop": locus["STOP"],
                "reference_n_snps": count,
                "reference_chromosome_n_snps": chromosome_count,
                "selection_basis": (
                    "INPUT_ONLY_COMPLETE_REFERENCE_CHROMOSOME_THEN_LOCUS_SNP_BURDEN_"
                    "MIN_MEDIAN_MAX_WITH_PREFROZEN_FALLBACK_ORDER"
                    if benchmark_phase == "DISCOVERY" else
                    "POST_FULL_FAMILY_BH_CONDITIONAL_ELIGIBILITY_THEN_COMPLETE_REFERENCE_"
                    "CHROMOSOME_AND_LOCUS_SNP_BURDEN_WITH_PREFROZEN_FALLBACK_ORDER"
                ),
                "execution_fingerprint": fingerprint,
            })
    return rows


def benchmark_selection(
    loci: list[dict[str, int]], counts: list[int], fingerprint: str,
    chromosome_counts: dict[int, int] | None = None,
) -> list[dict[str, object]]:
    """Return the three primary input-only representatives (public test helper)."""
    chromosome_counts = chromosome_counts or {chromosome: 0 for chromosome in range(1, 23)}
    plan = _ranked_candidates(
        loci, counts, chromosome_counts, fingerprint, benchmark_phase="DISCOVERY",
    )
    primary = [row for row in plan if int(row["candidate_rank"]) == 1]
    if len(primary) != 3:
        raise ExecutionError("fewer than three loci have reference SNPs")
    return primary


def run_root(fingerprint: str) -> Path:
    if not FINGERPRINT_RE.fullmatch(fingerprint):
        raise ExecutionError("invalid execution fingerprint")
    return CHECKPOINT_ROOT / fingerprint


def input_baseline_path(fingerprint: str) -> Path:
    return run_root(fingerprint) / "active_input_identity_baseline.json"


def _safe_root_path(value: str, label: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ExecutionError(f"unsafe {label} path: {value}")
    path = ROOT / relative
    try:
        path.resolve(strict=False).relative_to(ROOT.resolve())
    except ValueError as error:
        raise ExecutionError(f"{label} path escapes the repository: {value}") from error
    return path


def _expected_baseline_entries() -> list[dict[str, Any]]:
    policy = _load_json(POLICY, "Track B local policy")
    input_lock = _load_json(INPUT_LOCK, "Track B local input lock")
    reference = _load_json(REFERENCE_LOCK, "sealed LAVA reference")
    chromosome = _load_json(CHROMOSOME_INPUT_LOCK, "LAVA chromosome-input lock")
    records: dict[str, dict[str, Any]] = {}

    def add(path_value: str, digest: str, size: int | None, role: str, chromosome_value: int | None = None) -> None:
        path = _safe_root_path(path_value, role)
        if not FINGERPRINT_RE.fullmatch(str(digest)):
            raise ExecutionError(f"invalid expected SHA-256 for active input: {path_value}")
        observed_size = stable_stat(path)["bytes"]
        expected_size = observed_size if size is None else int(size)
        if observed_size != expected_size:
            raise ExecutionError(f"active input size differs from its frozen record: {path_value}")
        record = {
            "path": path_value, "bytes": expected_size, "sha256": str(digest), "role": role,
            "chromosome": chromosome_value,
        }
        prior = records.get(path_value)
        if prior is not None and prior != record:
            raise ExecutionError(f"conflicting active-input provenance: {path_value}")
        records[path_value] = record

    add(str(POLICY.relative_to(ROOT)), sha256(POLICY), POLICY.stat().st_size, "COMMON_INPUT")
    add(str(INPUT_LOCK.relative_to(ROOT)), sha256(INPUT_LOCK), INPUT_LOCK.stat().st_size, "COMMON_INPUT")
    add(str(REFERENCE_LOCK.relative_to(ROOT)), sha256(REFERENCE_LOCK), REFERENCE_LOCK.stat().st_size, "COMMON_INPUT")
    add(str(CHROMOSOME_INPUT_LOCK.relative_to(ROOT)), sha256(CHROMOSOME_INPUT_LOCK), CHROMOSOME_INPUT_LOCK.stat().st_size, "COMMON_INPUT")
    add(str(LOCUS_FILE.relative_to(ROOT)), str(policy["locus_definition"]["sha256"]), LOCUS_FILE.stat().st_size, "COMMON_INPUT")
    locked_artifacts = input_lock.get("input_artifact_sha256")
    if not isinstance(locked_artifacts, dict):
        raise ExecutionError("local input lock lacks its complete artifact hash map")
    for relative, digest in sorted(locked_artifacts.items()):
        path = _safe_root_path(str(relative), "common input")
        add(str(relative), str(digest), path.stat().st_size, "COMMON_INPUT")
    for item in reference.get("extracted_files", []):
        relative = str(item.get("path", ""))
        match = re.search(r"_chr(\d+)\.(?:info|bcor)$", relative)
        if not match:
            raise ExecutionError(f"malformed sealed reference chromosome path: {relative}")
        add(relative, str(item.get("sha256", "")), int(item.get("bytes", -1)), "REFERENCE_CHROMOSOME", int(match.group(1)))
    for item in chromosome.get("shards", []):
        add(str(item.get("path", "")), str(item.get("sha256", "")), int(item.get("bytes", -1)), "CHROMOSOME_SHARD", int(item.get("chromosome", 0)))
    for item in chromosome.get("chromosome_input_info", []):
        add(str(item.get("path", "")), str(item.get("sha256", "")), int(item.get("bytes", -1)), "CHROMOSOME_INPUT_INFO", int(item.get("chromosome", 0)))
    if len([item for item in records.values() if item["role"] == "REFERENCE_CHROMOSOME"]) != 44:
        raise ExecutionError("active-input baseline lacks 44 reference chromosome files")
    if len([item for item in records.values() if item["role"] == "CHROMOSOME_SHARD"]) != 176:
        raise ExecutionError("active-input baseline lacks 176 chromosome shards")
    return [records[key] for key in sorted(records)]


def _baseline_payload(fingerprint: str) -> dict[str, Any]:
    entries = []
    for expected in _expected_baseline_entries():
        path = ROOT / str(expected["path"])
        entries.append(expected | {"verified_stat": stable_stat(path)})
    return {
        "schema_version": 1,
        "analysis_id": "track-b-v1.0-local",
        "execution_fingerprint": fingerprint,
        "content_verification": "FULL_SHA256_PREFLIGHT_BEFORE_BASELINE_PUBLICATION",
        "entries": entries,
    }


def freeze_input_baseline(fingerprint: str) -> dict[str, Any]:
    expected = _baseline_payload(fingerprint)
    path = input_baseline_path(fingerprint)
    encoded = (json.dumps(expected, indent=2, sort_keys=True) + "\n").encode()
    if path.exists():
        if path.read_bytes() != encoded:
            raise ExecutionError("active-input identity baseline differs from fully verified live inputs")
    else:
        publish_no_replace(path, encoded)
    return validate_input_baseline(fingerprint)


def validate_input_baseline(fingerprint: str) -> dict[str, Any]:
    path = input_baseline_path(fingerprint)
    observed = _load_json(path, "active-input identity baseline")
    expected_entries = _expected_baseline_entries()
    observed_entries = observed.get("entries")
    if (
        observed.get("schema_version") != 1
        or observed.get("analysis_id") != "track-b-v1.0-local"
        or observed.get("execution_fingerprint") != fingerprint
        or observed.get("content_verification") != "FULL_SHA256_PREFLIGHT_BEFORE_BASELINE_PUBLICATION"
        or not isinstance(observed_entries, list)
        or len(observed_entries) != len(expected_entries)
    ):
        raise ExecutionError("active-input identity baseline scope drifted")
    for expected, item in zip(expected_entries, observed_entries, strict=True):
        if not isinstance(item, dict) or {key: item.get(key) for key in expected} != expected:
            raise ExecutionError("active-input identity baseline provenance drifted")
        current = stable_stat(ROOT / str(expected["path"]))
        if item.get("verified_stat") != current:
            raise ExecutionError(f"active input changed since full SHA-256 preflight: {expected['path']}")
    return observed


def selection_path(fingerprint: str, phase: str) -> Path:
    label = "discovery" if phase == "discovery" else "conditional"
    return run_root(fingerprint) / f"{label}_benchmark_candidates.tsv"


def freeze_selection(path: Path, rows: list[dict[str, object]]) -> None:
    expected = tsv_text(SELECTION_FIELDS, rows)
    if path.exists():
        if path.read_text(encoding="utf-8") != expected:
            raise ExecutionError(f"frozen RAM benchmark selection drifted: {path}")
        return
    publish_no_replace(path, expected.encode())


def preflight() -> str:
    if not R_SCRIPT.is_file():
        raise ExecutionError(f"pinned R runtime is missing: {R_SCRIPT}")
    result = subprocess.run(
        [sys.executable, str(CONTRACT), "--preflight"], cwd=ROOT,
        text=True, capture_output=True, check=False,
    )
    output = result.stdout + result.stderr
    if result.returncode != 0:
        raise ExecutionError(f"Track B LAVA preflight failed:\n{output}")
    fingerprints = re.findall(r"fingerprint=([0-9a-f]{64})", output)
    if len(fingerprints) != 1:
        raise ExecutionError("Track B LAVA preflight did not return one fingerprint")
    freeze_input_baseline(fingerprints[0])
    return fingerprints[0]


def bundle_name(phase: str, index: int | None) -> str:
    if phase in LOCUS_PHASES:
        if index is None or not 1 <= index <= 2495:
            raise ExecutionError(f"{phase} requires one valid locus index")
        return f"locus_{index:04d}"
    if phase in {"aggregate-discovery", "finalize"} and index is None:
        return "unit_all"
    raise ExecutionError(f"unsupported phase/index: {phase}/{index}")


def bundle_path(phase: str, index: int | None, fingerprint: str) -> Path:
    return run_root(fingerprint) / phase / bundle_name(phase, index)


def logical_result_path(phase: str, index: int | None, fingerprint: str) -> Path:
    return bundle_path(phase, index, fingerprint) / "result.rds"


def _checkpoint_result_spec(
    phase: str, index: int | None, fingerprint: str,
) -> dict[str, Any]:
    receipt = validate_bundle(
        phase, index, fingerprint, _reconcile_active=False, _revalidate_semantic=False,
    )
    if receipt is None:
        raise ExecutionError(f"missing prerequisite checkpoint: {phase}/{index}")
    bundle = bundle_path(phase, index, fingerprint)
    target = _safe_attempt_target(bundle)
    artifact = next(
        (item for item in receipt["artifacts"] if item["path"] == "result.rds"), None,
    )
    if artifact is None:
        raise ExecutionError(f"prerequisite checkpoint has no result artifact: {phase}/{index}")
    path = target / "result.rds"
    return {
        "path": str(path.relative_to(ROOT)), "bytes": int(artifact["bytes"]),
        "sha256": str(artifact["sha256"]), "role": "CHECKPOINT_PREREQUISITE",
        "chromosome": None,
    }


def _active_specs(phase: str, index: int | None, fingerprint: str) -> list[dict[str, Any]]:
    baseline = validate_input_baseline(fingerprint)
    chromosome_value = load_loci()[index - 1]["CHR"] if phase in LOCUS_PHASES and index is not None else None
    specs = [
        {key: item[key] for key in ("path", "bytes", "sha256", "role", "chromosome")}
        for item in baseline["entries"]
        if item["role"] == "COMMON_INPUT"
        or (
            chromosome_value is not None
            and item["role"] in {"REFERENCE_CHROMOSOME", "CHROMOSOME_SHARD", "CHROMOSOME_INPUT_INFO"}
            and int(item["chromosome"]) == chromosome_value
        )
    ]
    if phase == "aggregate-discovery":
        specs.extend(_checkpoint_result_spec("discovery", locus_index, fingerprint) for locus_index in range(1, 2496))
    elif phase == "conditional":
        specs.append(_checkpoint_result_spec("aggregate-discovery", None, fingerprint))
    elif phase == "finalize":
        specs.append(_checkpoint_result_spec("aggregate-discovery", None, fingerprint))
        specs.extend(_checkpoint_result_spec("conditional", locus_index, fingerprint) for locus_index in range(1, 2496))
    specs.sort(key=lambda item: str(item["path"]))
    if len({str(item["path"]) for item in specs}) != len(specs):
        raise ExecutionError(f"duplicate active-input dependency in {phase} worker")
    return specs


def capture_active_inputs(phase: str, index: int | None, fingerprint: str) -> list[dict[str, Any]]:
    captured = []
    for spec in _active_specs(phase, index, fingerprint):
        path = ROOT / str(spec["path"])
        observed = stable_stat(path)
        if observed["bytes"] != int(spec["bytes"]):
            raise ExecutionError(f"active input has the wrong byte count: {spec['path']}")
        if spec["role"] == "CHECKPOINT_PREREQUISITE" and sha256(path) != spec["sha256"]:
            raise ExecutionError(f"checkpoint prerequisite differs from its READY receipt: {spec['path']}")
        captured.append(spec | {"pre_stat": observed})
    return captured


def complete_active_inputs(captured: list[dict[str, Any]]) -> list[dict[str, Any]]:
    completed = []
    for item in captured:
        path = ROOT / str(item["path"])
        observed = stable_stat(path)
        if observed != item["pre_stat"]:
            raise ExecutionError(f"active input changed during worker execution: {item['path']}")
        if item["role"] == "CHECKPOINT_PREREQUISITE" and sha256(path) != item["sha256"]:
            raise ExecutionError(f"checkpoint prerequisite changed during worker execution: {item['path']}")
        completed.append(item | {"post_stat": observed})
    return completed


def parse_marker(log: str, phase: str) -> dict[str, str]:
    lines = [line for line in log.splitlines() if line.startswith(MARKER_PREFIX)]
    if len(lines) != 1:
        raise ExecutionError(f"{phase} worker did not emit exactly one metadata marker")
    values: dict[str, str] = {}
    for field in lines[0].split("\t")[1:]:
        if "=" not in field:
            raise ExecutionError(f"malformed {phase} worker metadata marker")
        key, value = field.split("=", 1)
        values[key] = value
    required = {
        "phase", "index", "locus", "chromosome", "n_snps", "pair", "qc",
        "construction_complete", "output",
    }
    if set(values) != required or values["phase"] != phase:
        raise ExecutionError(f"incomplete or wrong-phase {phase} worker metadata marker")
    return values


def validate_marker(
    marker: dict[str, str], phase: str, index: int | None, expected_output: Path,
) -> None:
    expected_index = 0 if index is None else index
    if marker["index"] != str(expected_index):
        raise ExecutionError(f"{phase} marker has the wrong locus index")
    if phase in LOCUS_PHASES:
        locus = load_loci()[expected_index - 1]
        if marker["locus"] != str(locus["LOC"]) or marker["chromosome"] != str(locus["CHR"]):
            raise ExecutionError(f"{phase} marker has wrong official locus coordinates")
    elif marker["locus"] != "ALL" or marker["chromosome"] != "ALL":
        raise ExecutionError(f"{phase} aggregate marker is malformed")
    try:
        marker_output = (ROOT / marker["output"]).resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExecutionError(f"{phase} marker names an invalid output") from error
    if marker_output != expected_output.resolve(strict=True):
        raise ExecutionError(f"{phase} marker output differs from the supervised attempt")
    if marker["n_snps"] != "NA":
        try:
            if int(float(marker["n_snps"])) < 0:
                raise ValueError
        except ValueError as error:
            raise ExecutionError(f"{phase} marker has invalid SNP count") from error
    if marker["construction_complete"] not in {"TRUE", "FALSE"}:
        raise ExecutionError(f"{phase} marker has invalid construction-complete state")


def quick_execution_fingerprint(expected: str) -> None:
    result = subprocess.run(
        [sys.executable, str(CONTRACT), "--quick-execution-fingerprint"], cwd=ROOT,
        text=True, capture_output=True, check=False,
    )
    hashes = [line for line in (result.stdout + result.stderr).splitlines() if FINGERPRINT_RE.fullmatch(line)]
    if result.returncode != 0 or hashes != [expected]:
        raise ExecutionError("live Track B LAVA inputs/runtime differ from the supervised fingerprint")


def _semantic_validation(
    phase: str, index: int | None, fingerprint: str, result_path: Path,
) -> tuple[dict[str, str], str]:
    command = [
        str(R_SCRIPT), str(CHECKPOINT_VALIDATOR), "--phase", phase,
        "--fingerprint", fingerprint, "--rds", str(result_path),
    ]
    if index is not None:
        command.extend(["--locus-index", str(index)])
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    combined = "# command\n" + " ".join(command) + "\n# stdout\n" + result.stdout + "# stderr\n" + result.stderr
    markers = [line for line in result.stdout.splitlines() if line.startswith(SEMANTIC_MARKER_PREFIX)]
    if result.returncode != 0 or len(markers) != 1:
        raise ExecutionError(f"{phase} checkpoint failed RDS semantic validation:\n{combined}")
    values: dict[str, str] = {}
    for field in markers[0].split("\t")[1:]:
        if "=" not in field:
            raise ExecutionError(f"malformed {phase} checkpoint semantic marker")
        key, value = field.split("=", 1)
        values[key] = value
    required = {
        "phase", "index", "locus", "chromosome", "n_snps", "pair", "qc",
        "construction_complete", "fingerprint", "status",
    }
    if (
        set(values) != required or values["phase"] != phase
        or values["index"] != str(0 if index is None else index)
        or values["fingerprint"] != fingerprint or values["status"] != "PASS"
    ):
        raise ExecutionError(f"{phase} checkpoint semantic marker scope drifted")
    if phase == "finalize":
        staged_command = [
            sys.executable, str(VALIDATOR), "--verify-staged-results", "--staging-dir",
            str(result_path.parent), "--fingerprint", fingerprint,
        ]
        staged = subprocess.run(staged_command, cwd=ROOT, text=True, capture_output=True, check=False)
        combined += "# staged command\n" + " ".join(staged_command) + "\n# staged stdout\n" + staged.stdout + "# staged stderr\n" + staged.stderr
        if staged.returncode != 0 or sum(
            line.startswith("TRACK_B_LAVA_STAGED_RESULTS_READ_ONLY_VERIFIED")
            for line in staged.stdout.splitlines()
        ) != 1:
            raise ExecutionError(f"finalize staged family failed read-only semantic validation:\n{combined}")
    return values, combined


def _require_semantic_marker_agreement(worker: dict[str, str], semantic: dict[str, str]) -> None:
    for key in ("phase", "index", "locus", "chromosome", "n_snps", "pair", "qc", "construction_complete"):
        if worker[key] != semantic[key]:
            raise ExecutionError(f"worker marker differs from deserialized checkpoint semantics: {key}")


def rss_bytes_from_rusage(value: int) -> int:
    return int(value if sys.platform == "darwin" else value * 1024)


def internal_measure(args: argparse.Namespace) -> int:
    command = [
        str(R_SCRIPT), str(RUNNER), "--phase", args.phase,
        "--fingerprint", args.fingerprint, "--worker-output", args.worker_output,
    ]
    if args.locus_index is not None:
        command.extend(["--locus-index", str(args.locus_index)])
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    runtime = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    payload = {
        "command": command, "returncode": result.returncode,
        "runtime_sec": runtime,
        "peak_rss_bytes": rss_bytes_from_rusage(int(usage.ru_maxrss)),
        "stdout": result.stdout, "stderr": result.stderr,
    }
    Path(args.internal_result).write_text(json.dumps(payload), encoding="utf-8")
    return 0


def _safe_attempt_target(bundle: Path) -> Path:
    if not bundle.is_symlink():
        raise ExecutionError(f"published checkpoint is not an atomic ready link: {bundle}")
    phase_dir = bundle.parent.resolve()
    attempts = (bundle.parent / ".attempts").resolve()
    try:
        target = bundle.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExecutionError(f"checkpoint ready link is broken: {bundle}") from error
    if target.parent != attempts or attempts.parent != phase_dir:
        raise ExecutionError(f"checkpoint ready link escapes its fingerprint/phase namespace: {bundle}")
    return target


def validate_bundle(
    phase: str, index: int | None, fingerprint: str, analysis: str | None = None,
    *, _reconcile_active: bool = True, _revalidate_semantic: bool = True,
) -> dict[str, Any] | None:
    bundle = bundle_path(phase, index, fingerprint)
    if not bundle.exists() and not bundle.is_symlink():
        return None
    target = _safe_attempt_target(bundle)
    ready_path, receipt_path = target / "READY", target / "receipt.json"
    log_path, result_path = target / "worker.log", target / "result.rds"
    semantic_path = target / "semantic_validation.txt"
    for path in (ready_path, receipt_path, log_path, result_path, semantic_path):
        if not path.is_file() or path.stat().st_size <= 0:
            raise ExecutionError(f"checkpoint bundle is incomplete: {bundle} ({path.name})")
    try:
        ready = json.loads(ready_path.read_text(encoding="utf-8"))
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ExecutionError(f"checkpoint bundle metadata is unreadable: {bundle}") from error
    if ready != {"schema_version": 1, "receipt_sha256": sha256(receipt_path), "state": "READY"}:
        raise ExecutionError(f"checkpoint READY marker differs from its receipt: {bundle}")
    marker = receipt.get("marker")
    expected = {
        "schema_version": 3, "analysis": analysis or ANALYSIS_BY_PHASE[phase],
        "phase": phase, "locus_index": index, "execution_fingerprint": fingerprint,
        "exit_status": 0, "bundle_path": str(bundle.relative_to(ROOT)),
        "log_bytes": log_path.stat().st_size, "log_sha256": sha256(log_path),
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise ExecutionError(f"checkpoint receipt differs from live bundle: {bundle} ({key})")
    if not isinstance(marker, dict):
        raise ExecutionError(f"checkpoint receipt lacks worker marker: {bundle}")
    validate_marker(marker, phase, index, result_path)
    if receipt.get("semantic_validator_sha256") != sha256(CHECKPOINT_VALIDATOR):
        raise ExecutionError(f"checkpoint receipt is not bound to the live semantic validator: {bundle}")
    baseline = input_baseline_path(fingerprint)
    if receipt.get("input_baseline_sha256") != sha256(baseline):
        raise ExecutionError(f"checkpoint receipt differs from the active-input baseline: {bundle}")
    if not isinstance(receipt.get("peak_rss_bytes"), int) or int(receipt["peak_rss_bytes"]) <= 0:
        raise ExecutionError(f"checkpoint receipt lacks measured peak RSS: {bundle}")
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ExecutionError(f"checkpoint receipt lacks artifact family: {bundle}")
    recorded_files: set[str] = set()
    for item in artifacts:
        if not isinstance(item, dict) or set(item) != {"path", "bytes", "sha256"}:
            raise ExecutionError(f"malformed checkpoint artifact receipt: {bundle}")
        relative = str(item["path"])
        if relative in recorded_files or "/" in relative or relative.startswith("."):
            raise ExecutionError(f"invalid checkpoint artifact path: {bundle}/{relative}")
        size, digest = file_identity(target / relative)
        if (size, digest) != (item["bytes"], item["sha256"]):
            raise ExecutionError(f"checkpoint artifact differs from receipt: {bundle}/{relative}")
        recorded_files.add(relative)
    live_files = {path.name for path in target.iterdir() if path.is_file()}
    if live_files != recorded_files | {"worker.log", "receipt.json", "READY"}:
        raise ExecutionError(f"checkpoint bundle contains unreceipted files: {bundle}")
    if _revalidate_semantic:
        semantic, semantic_text = _semantic_validation(phase, index, fingerprint, result_path)
        if semantic_path.read_text(encoding="utf-8") != semantic_text:
            raise ExecutionError(f"checkpoint semantic attestation differs from read-only revalidation: {bundle}")
        _require_semantic_marker_agreement(marker, semantic)
    active_inputs = receipt.get("active_inputs")
    if not isinstance(active_inputs, list) or not active_inputs:
        raise ExecutionError(f"checkpoint receipt lacks active-input identities: {bundle}")
    expected_keys = {
        "path", "bytes", "sha256", "role", "chromosome", "pre_stat", "post_stat",
    }
    for item in active_inputs:
        if not isinstance(item, dict) or set(item) != expected_keys:
            raise ExecutionError(f"malformed active-input receipt: {bundle}")
        current = stable_stat(ROOT / str(item["path"]))
        if current != item["pre_stat"] or current != item["post_stat"] or current["bytes"] != int(item["bytes"]):
            raise ExecutionError(f"active input differs from the worker receipt: {item['path']}")
        if item["role"] == "CHECKPOINT_PREREQUISITE" and sha256(ROOT / str(item["path"])) != item["sha256"]:
            raise ExecutionError(f"checkpoint dependency differs from the worker receipt: {item['path']}")
    if _reconcile_active:
        expected_specs = _active_specs(phase, index, fingerprint)
        observed_specs = [
            {key: item[key] for key in ("path", "bytes", "sha256", "role", "chromosome")}
            for item in active_inputs
        ]
        if observed_specs != expected_specs:
            raise ExecutionError(f"checkpoint active-input family is incomplete or drifted: {bundle}")
    return receipt


def _benchmark_exit_status(receipt: dict[str, Any]) -> str:
    phase, qc = str(receipt["phase"]), str(receipt["marker"]["qc"])
    return "0" if qc in SCIENTIFICALLY_COMPLETE_QC[phase] else f"0:{qc}"


def benchmark_row(receipt: dict[str, Any]) -> dict[str, object]:
    marker = receipt["marker"]
    return {
        "analysis": receipt["analysis"], "pair": marker["pair"],
        "locus": marker["locus"], "chromosome": marker["chromosome"],
        "n_snps": marker["n_snps"],
        "peak_ram_gb": f"{int(receipt['peak_rss_bytes']) / 1024**3:.6f}",
        "runtime_sec": f"{float(receipt['runtime_sec']):.3f}",
        "exit_status": _benchmark_exit_status(receipt),
        "output_hash": next(
            str(item["sha256"]) for item in receipt["artifacts"] if item["path"] == "result.rds"
        ),
    }


def assert_benchmark_namespace(fingerprint: str) -> None:
    if BENCHMARK_NAMESPACE.is_file():
        payload = _load_json(BENCHMARK_NAMESPACE, "RAM benchmark namespace")
        if payload != {
            "schema_version": 1, "analysis_id": "track-b-v1.0-local",
            "execution_fingerprint": fingerprint,
        }:
            raise ExecutionError("RAM benchmark table belongs to a different execution fingerprint")
        return
    fields, rows = read_tsv(BENCHMARK)
    if fields != BENCHMARK_FIELDS or rows:
        raise ExecutionError("unsealed historical RAM benchmark rows cannot enter this execution")
    publish_no_replace(BENCHMARK_NAMESPACE, (json.dumps({
        "schema_version": 1, "analysis_id": "track-b-v1.0-local",
        "execution_fingerprint": fingerprint,
    }, indent=2, sort_keys=True) + "\n").encode())


def update_benchmark(row: dict[str, object], fingerprint: str) -> None:
    assert_benchmark_namespace(fingerprint)
    rows: list[dict[str, object]] = []
    if BENCHMARK.exists():
        fields, observed = read_tsv(BENCHMARK)
        if fields != BENCHMARK_FIELDS:
            raise ExecutionError("RAM benchmark schema drifted")
        rows.extend(observed)
    key = (
        str(row["analysis"]), str(row["pair"]), str(row["locus"]),
        str(row["chromosome"]), str(row["output_hash"]),
    )
    rows = [
        item for item in rows
        if (
            item["analysis"], item["pair"], item["locus"], item["chromosome"],
            item["output_hash"],
        ) != key
    ]
    rows.append(row)
    rows.sort(key=lambda item: (
        str(item["analysis"]), str(item["pair"]),
        10**9 if str(item["chromosome"]) in {"ALL", "NA"} else int(item["chromosome"]),
        10**12 if str(item["locus"]) in {"ALL", "NA"} else int(item["locus"]),
    ))
    atomic_text(BENCHMARK, tsv_text(BENCHMARK_FIELDS, rows))


def failure_row(
    analysis: str, index: int | None, returncode: int, runtime: float,
    peak: int, failure_log: Path, fingerprint: str,
) -> None:
    if index is None:
        locus, chromosome = "ALL", "ALL"
    else:
        official = load_loci()[index - 1]
        locus, chromosome = str(official["LOC"]), str(official["CHR"])
    update_benchmark({
        "analysis": analysis, "pair": "UNKNOWN", "locus": locus,
        "chromosome": chromosome, "n_snps": "NA",
        "peak_ram_gb": f"{peak / 1024**3:.6f}", "runtime_sec": f"{runtime:.3f}",
        "exit_status": str(returncode), "output_hash": sha256(failure_log),
    }, fingerprint)


def _attempt_directory(phase: str, index: int | None, fingerprint: str) -> Path:
    attempts = run_root(fingerprint) / phase / ".attempts"
    attempts.mkdir(parents=True, exist_ok=True)
    attempt = attempts / f"{bundle_name(phase, index)}.{os.getpid()}.{time.time_ns()}"
    attempt.mkdir()
    fsync_directory(attempts)
    return attempt


def _publish_ready_bundle(attempt: Path, phase: str, index: int | None, fingerprint: str) -> None:
    bundle = bundle_path(phase, index, fingerprint)
    try:
        os.symlink(os.path.relpath(attempt, start=bundle.parent), bundle, target_is_directory=True)
    except FileExistsError as error:
        raise ExecutionError(f"checkpoint ready link appeared during publication: {bundle}") from error
    fsync_directory(bundle.parent)


def validate_phase_prerequisites(phase: str, fingerprint: str) -> None:
    if phase == "aggregate-discovery":
        for locus_index in range(1, 2496):
            if validate_bundle("discovery", locus_index, fingerprint) is None:
                raise ExecutionError(f"discovery checkpoint family is incomplete at locus {locus_index}")
    elif phase == "conditional":
        if validate_bundle("aggregate-discovery", None, fingerprint, _reconcile_active=False) is None:
            raise ExecutionError("conditional worker lacks a verified full-family discovery aggregation")
    elif phase == "finalize":
        validate_phase_prerequisites("aggregate-discovery", fingerprint)
        if validate_bundle("aggregate-discovery", None, fingerprint, _reconcile_active=False) is None:
            raise ExecutionError("finalization lacks a verified discovery aggregation")
        for locus_index in range(1, 2496):
            if validate_bundle("conditional", locus_index, fingerprint) is None:
                raise ExecutionError(f"conditional checkpoint family is incomplete at locus {locus_index}")


def run_measured(
    phase: str, index: int | None, fingerprint: str, analysis: str | None = None,
) -> dict[str, Any]:
    if phase not in ALL_PHASES:
        raise ExecutionError(f"unsupported measured phase: {phase}")
    analysis = analysis or ANALYSIS_BY_PHASE[phase]
    validate_phase_prerequisites(phase, fingerprint)
    existing = validate_bundle(phase, index, fingerprint, analysis)
    if existing is not None:
        update_benchmark(benchmark_row(existing), fingerprint)
        assert_peak_fits_machine(existing)
        return existing
    quick_execution_fingerprint(fingerprint)
    pre_active_inputs = capture_active_inputs(phase, index, fingerprint)
    attempt = _attempt_directory(phase, index, fingerprint)
    output = attempt / "result.rds"
    with tempfile.NamedTemporaryFile(prefix="track-b-lava-measure-", suffix=".json", delete=False) as handle:
        metrics_path = Path(handle.name)
    command = [
        sys.executable, str(Path(__file__).resolve()), "--internal-measure",
        "--phase", phase, "--fingerprint", fingerprint,
        "--worker-output", str(output), "--internal-result", str(metrics_path),
    ]
    if index is not None:
        command.extend(["--locus-index", str(index)])
    try:
        helper = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
        if helper.returncode != 0 or not metrics_path.is_file():
            raise ExecutionError(f"RSS measurement helper failed: {helper.stdout}{helper.stderr}")
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    finally:
        metrics_path.unlink(missing_ok=True)
    log = (
        "# command\n" + " ".join(str(value) for value in metrics["command"]) +
        "\n# stdout\n" + str(metrics["stdout"]) + "\n# stderr\n" + str(metrics["stderr"])
    )
    returncode, runtime = int(metrics["returncode"]), float(metrics["runtime_sec"])
    peak = int(metrics["peak_rss_bytes"])
    quick_execution_fingerprint(fingerprint)
    active_inputs = complete_active_inputs(pre_active_inputs)
    if returncode != 0:
        failure_dir = run_root(fingerprint) / "failed_attempts"
        failure_dir.mkdir(parents=True, exist_ok=True)
        failure_path = failure_dir / f"{phase}_{index if index is not None else 'all'}_{time.time_ns()}.log"
        publish_no_replace(failure_path, log.encode())
        failure_row(analysis, index, returncode, runtime, peak, failure_path, fingerprint)
        build_report(fingerprint)
        raise ExecutionError(f"{phase} worker failed with exit status {returncode}; log: {failure_path}")
    if not output.is_file() or output.stat().st_size <= 0:
        raise ExecutionError(f"successful {phase} worker did not create its attempt result: {output}")
    marker = parse_marker(log, phase)
    validate_marker(marker, phase, index, output)
    try:
        semantic, semantic_text = _semantic_validation(phase, index, fingerprint, output)
        _require_semantic_marker_agreement(marker, semantic)
    except ExecutionError as error:
        failure_dir = run_root(fingerprint) / "failed_attempts"
        failure_dir.mkdir(parents=True, exist_ok=True)
        failure_path = failure_dir / f"{phase}_{index if index is not None else 'all'}_{time.time_ns()}_semantic.log"
        publish_no_replace(failure_path, (log + "\n# semantic failure\n" + str(error) + "\n").encode())
        failure_row(analysis, index, 70, runtime, peak, failure_path, fingerprint)
        build_report(fingerprint)
        raise
    log_path = attempt / "worker.log"
    publish_no_replace(log_path, log.encode())
    publish_no_replace(attempt / "semantic_validation.txt", semantic_text.encode())
    artifact_paths = sorted(
        [path for path in attempt.iterdir() if path.is_file() and path.name not in {"worker.log", "receipt.json", "READY"}],
        key=lambda path: path.name,
    )
    artifacts = []
    for path in artifact_paths:
        size, digest = file_identity(path)
        artifacts.append({"path": path.name, "bytes": size, "sha256": digest})
    receipt: dict[str, Any] = {
        "schema_version": 3, "analysis": analysis, "phase": phase,
        "locus_index": index, "execution_fingerprint": fingerprint,
        "marker": marker, "command": metrics["command"], "runtime_sec": runtime,
        "peak_rss_bytes": peak, "exit_status": 0,
        "bundle_path": str(bundle_path(phase, index, fingerprint).relative_to(ROOT)),
        "log_bytes": log_path.stat().st_size, "log_sha256": sha256(log_path),
        "semantic_validator_sha256": sha256(CHECKPOINT_VALIDATOR),
        "input_baseline_sha256": sha256(input_baseline_path(fingerprint)),
        "active_inputs": active_inputs,
        "artifacts": artifacts,
    }
    receipt_path = attempt / "receipt.json"
    publish_no_replace(receipt_path, (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode())
    ready = {"schema_version": 1, "receipt_sha256": sha256(receipt_path), "state": "READY"}
    publish_no_replace(attempt / "READY", (json.dumps(ready, sort_keys=True) + "\n").encode())
    fsync_directory(attempt)
    _publish_ready_bundle(attempt, phase, index, fingerprint)
    published = validate_bundle(phase, index, fingerprint, analysis)
    if published is None:
        raise ExecutionError("published checkpoint bundle disappeared")
    update_benchmark(benchmark_row(published), fingerprint)
    assert_peak_fits_machine(published)
    return published


def _representative_receipts(fingerprint: str, phase: str) -> dict[str, dict[str, Any]]:
    path = selection_path(fingerprint, phase)
    if not path.is_file():
        return {}
    fields, candidates = read_tsv(path)
    if fields != SELECTION_FIELDS:
        raise ExecutionError(f"RAM benchmark selection schema drifted: {path}")
    selected: dict[str, dict[str, Any]] = {}
    used_loci: set[str] = set()
    for role in dict.fromkeys(row["selection_role"] for row in candidates):
        role_candidates = sorted(
            (row for row in candidates if row["selection_role"] == role),
            key=lambda row: int(row["candidate_rank"]),
        )
        for candidate in role_candidates:
            receipt = validate_bundle(phase, int(candidate["locus_index"]), fingerprint)
            if receipt is None:
                continue
            marker = receipt["marker"]
            if marker["qc"] in SCIENTIFICALLY_COMPLETE_QC[phase] and marker["locus"] not in used_loci:
                selected[role] = receipt
                used_loci.add(marker["locus"])
                break
    return selected


def _expected_roles(fingerprint: str, phase: str) -> list[str]:
    path = selection_path(fingerprint, phase)
    if not path.is_file():
        return []
    fields, rows = read_tsv(path)
    if fields != SELECTION_FIELDS:
        raise ExecutionError("RAM benchmark selection schema drifted")
    return list(dict.fromkeys(row["selection_role"] for row in rows))


def assert_representatives_fit_machine(fingerprint: str, phase: str) -> dict[str, dict[str, Any]]:
    selected, expected_roles = _representative_receipts(fingerprint, phase), _expected_roles(fingerprint, phase)
    if set(selected) != set(expected_roles) or not expected_roles:
        raise ExecutionError(f"representative {phase} RAM benchmark is incomplete")
    measured = measured_admission_receipts(fingerprint)
    _, benchmark_rows = read_tsv(BENCHMARK)
    attempted_peaks = [
        int(round(float(row["peak_ram_gb"]) * 1024**3))
        for row in benchmark_rows if row["analysis"].startswith("LAVA_")
    ]
    maximum = max(
        [int(receipt["peak_rss_bytes"]) for receipt in measured] + attempted_peaks,
    )
    memory, reserve = physical_memory_bytes(), memory_safety_reserve_bytes()
    if maximum + reserve > memory:
        raise ExecutionError(
            "BLOCKED_BY_MEASURED_PER_LOCUS_RAM: an intact representative LAVA "
            f"{phase} worker used {maximum / 1024**3:.3f} GiB and the predeclared "
            f"reserve is {reserve / 1024**3:.3f} GiB versus {memory / 1024**3:.3f} GiB physical RAM"
        )
    return selected


def _ready_phase_receipts(fingerprint: str, phase: str) -> list[dict[str, Any]]:
    directory = run_root(fingerprint) / phase
    if not directory.is_dir():
        return []
    receipts = []
    if phase in LOCUS_PHASES:
        for path in sorted(directory.glob("locus_[0-9][0-9][0-9][0-9]")):
            try:
                index = int(path.name.removeprefix("locus_"))
            except ValueError:
                continue
            receipt = validate_bundle(
                phase, index, fingerprint, _reconcile_active=False, _revalidate_semantic=False,
            )
            if receipt is not None:
                receipts.append(receipt)
    else:
        receipt = validate_bundle(
            phase, None, fingerprint, _reconcile_active=False, _revalidate_semantic=False,
        )
        if receipt is not None:
            receipts.append(receipt)
    return receipts


def measured_admission_receipts(fingerprint: str) -> list[dict[str, Any]]:
    candidate_indices: dict[str, set[int]] = {"discovery": set(), "conditional": set()}
    for phase in LOCUS_PHASES:
        path = selection_path(fingerprint, phase)
        if path.is_file():
            fields, rows = read_tsv(path)
            if fields != SELECTION_FIELDS:
                raise ExecutionError(f"RAM benchmark selection schema drifted: {path}")
            candidate_indices[phase] = {int(row["locus_index"]) for row in rows}
    measured = []
    for phase in ("discovery", "aggregate-discovery", "conditional", "finalize"):
        for receipt in _ready_phase_receipts(fingerprint, phase):
            construction_complete = receipt["marker"]["construction_complete"] == "TRUE"
            candidate_attempt = (
                phase in LOCUS_PHASES and receipt["locus_index"] in candidate_indices[phase]
            )
            if construction_complete or candidate_attempt:
                measured.append(receipt)
    return measured


def assert_peak_fits_machine(receipt: dict[str, Any]) -> None:
    peak = int(receipt["peak_rss_bytes"])
    memory, reserve = physical_memory_bytes(), memory_safety_reserve_bytes()
    if peak + reserve > memory:
        raise ExecutionError(
            "BLOCKED_BY_MEASURED_UNIT_RAM: a scientifically indivisible LAVA unit "
            f"used {peak / 1024**3:.3f} GiB and the predeclared reserve is "
            f"{reserve / 1024**3:.3f} GiB versus {memory / 1024**3:.3f} GiB physical RAM"
        )


def _run_candidate_plan(fingerprint: str, phase: str, rows: list[dict[str, object]]) -> dict[str, dict[str, Any]]:
    roles = list(dict.fromkeys(str(row["selection_role"]) for row in rows))
    selected: dict[str, dict[str, Any]] = {}
    used_loci: set[str] = set()
    for role in roles:
        candidates = sorted(
            (row for row in rows if str(row["selection_role"]) == role),
            key=lambda row: int(row["candidate_rank"]),
        )
        for candidate in candidates:
            index = int(candidate["locus_index"])
            receipt = run_measured(phase, index, fingerprint)
            marker = receipt["marker"]
            if marker["qc"] in SCIENTIFICALLY_COMPLETE_QC[phase] and marker["locus"] not in used_loci:
                selected[role] = receipt
                used_loci.add(marker["locus"])
                print(
                    f"LAVA_{phase.upper()}_BENCHMARK_COMPLETE role={role} index={index} "
                    f"locus={marker['locus']} n_snps={marker['n_snps']} "
                    f"peak_gib={int(receipt['peak_rss_bytes']) / 1024**3:.3f}"
                )
                break
        if role not in selected:
            raise ExecutionError(f"no scientifically completed predeclared fallback remained for {phase}/{role}")
    assert_representatives_fit_machine(fingerprint, phase)
    write_benchmark_provenance(fingerprint)
    return selected


def discovery_candidate_plan(fingerprint: str) -> list[dict[str, object]]:
    loci = load_loci()
    counts, chromosome_counts = reference_snp_burdens(loci)
    rows = _ranked_candidates(loci, counts, chromosome_counts, fingerprint, benchmark_phase="DISCOVERY")
    if len({row["selection_role"] for row in rows}) != 3:
        raise ExecutionError("fewer than three intact LAVA discovery benchmark roles are available")
    freeze_selection(selection_path(fingerprint, "discovery"), rows)
    return rows


def conditional_candidate_plan(fingerprint: str) -> list[dict[str, object]]:
    if validate_bundle("aggregate-discovery", None, fingerprint) is None:
        raise ExecutionError("full discovery aggregation is required before conditional benchmarking")
    candidate_path = bundle_path("aggregate-discovery", None, fingerprint) / "conditional_candidates.tsv"
    fields, candidates = read_tsv(candidate_path)
    if fields != CONDITIONAL_CANDIDATE_FIELDS:
        raise ExecutionError("conditional benchmark candidate schema drifted")
    eligible_indices = {int(row["locus_index"]) for row in candidates}
    loci = load_loci()
    reference_counts, chromosome_counts = reference_snp_burdens(loci)
    rows = _ranked_candidates(
        loci, reference_counts, chromosome_counts, fingerprint,
        benchmark_phase="CONDITIONAL", eligible_indices=eligible_indices,
    )
    freeze_selection(selection_path(fingerprint, "conditional"), rows)
    return rows


def run_discovery_benchmarks(fingerprint: str) -> dict[str, dict[str, Any]]:
    return _run_candidate_plan(fingerprint, "discovery", discovery_candidate_plan(fingerprint))


def _receipt_evidence(receipt: dict[str, Any]) -> dict[str, Any]:
    bundle = ROOT / str(receipt["bundle_path"])
    target = _safe_attempt_target(bundle)
    return {
        "phase": receipt["phase"], "locus_index": receipt["locus_index"],
        "bundle_path": receipt["bundle_path"],
        "receipt_sha256": sha256(target / "receipt.json"),
        "peak_rss_bytes": receipt["peak_rss_bytes"], "runtime_sec": receipt["runtime_sec"],
        "qc": receipt["marker"]["qc"], "locus": receipt["marker"]["locus"],
        "chromosome": receipt["marker"]["chromosome"],
        "construction_complete": receipt["marker"]["construction_complete"],
        "output_sha256": next(
            item["sha256"] for item in receipt["artifacts"] if item["path"] == "result.rds"
        ),
    }


def benchmark_provenance_payload(fingerprint: str) -> dict[str, Any]:
    if not BENCHMARK.is_file():
        raise ExecutionError("RAM benchmark table is missing")
    discovery = _representative_receipts(fingerprint, "discovery")
    conditional = _representative_receipts(fingerprint, "conditional")
    discovery_roles = _expected_roles(fingerprint, "discovery")
    conditional_roles = _expected_roles(fingerprint, "conditional")
    conditional_selection = selection_path(fingerprint, "conditional")
    if set(discovery) == set(discovery_roles) and discovery_roles:
        if conditional_selection.exists() and set(conditional) == set(conditional_roles):
            state = "ALL_REQUIRED_PHASE_REPRESENTATIVES_COMPLETE"
        else:
            state = "DISCOVERY_REPRESENTATIVES_COMPLETE_CONDITIONAL_PENDING"
    else:
        state = "IN_PROGRESS"
    selections = []
    for phase in ("discovery", "conditional"):
        path = selection_path(fingerprint, phase)
        if path.is_file():
            size, digest = file_identity(path)
            selections.append({
                "phase": phase, "path": str(path.relative_to(ROOT)),
                "bytes": size, "sha256": digest,
            })
    measured = measured_admission_receipts(fingerprint)
    fields, benchmark_rows = read_tsv(BENCHMARK)
    if fields != BENCHMARK_FIELDS:
        raise ExecutionError("RAM benchmark schema drifted")
    lava_attempt_rows = [dict(row) for row in benchmark_rows if row["analysis"].startswith("LAVA_")]
    attempted_peaks = [int(round(float(row["peak_ram_gb"]) * 1024**3)) for row in lava_attempt_rows]
    return {
        "schema_version": 2, "analysis_id": "track-b-v1.0-local",
        "execution_fingerprint": fingerprint, "state": state,
        "physical_memory_bytes": physical_memory_bytes(),
        "memory_safety_reserve_bytes": memory_safety_reserve_bytes(),
        "benchmark_path": str(BENCHMARK.relative_to(ROOT)),
        "benchmark_bytes": BENCHMARK.stat().st_size, "benchmark_sha256": sha256(BENCHMARK),
        "selections": selections,
        "selected_representatives": {
            "discovery": {role: _receipt_evidence(receipt) for role, receipt in discovery.items()},
            "conditional": {role: _receipt_evidence(receipt) for role, receipt in conditional.items()},
        },
        "measured_admission_units": [_receipt_evidence(receipt) for receipt in measured],
        "all_lava_attempt_rows": lava_attempt_rows,
        "maximum_measured_admission_peak_rss_bytes": max(
            [int(receipt["peak_rss_bytes"]) for receipt in measured] + attempted_peaks,
            default=None,
        ),
    }


def write_benchmark_provenance(fingerprint: str) -> None:
    payload = benchmark_provenance_payload(fingerprint)
    atomic_text(BENCHMARK_PROVENANCE, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def validate_benchmark_provenance(fingerprint: str) -> dict[str, Any]:
    try:
        observed = json.loads(BENCHMARK_PROVENANCE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ExecutionError(f"RAM benchmark provenance is unreadable: {error}") from error
    expected = benchmark_provenance_payload(fingerprint)
    if observed != expected:
        raise ExecutionError("RAM benchmark provenance differs from live fingerprint-bound evidence")
    return observed


def report_statistics(rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| Analysis | Scientifically completed units | Max peak RAM (GiB) | Median runtime (s) | Total observed runtime (h) |",
        "|---|---:|---:|---:|---:|",
    ]
    for analysis in sorted({row["analysis"] for row in rows}):
        observed = [row for row in rows if row["analysis"] == analysis and row["exit_status"] == "0"]
        peaks = [float(row["peak_ram_gb"]) for row in observed]
        runtimes = [float(row["runtime_sec"]) for row in observed]
        lines.append(
            f"| {analysis} | {len(observed)} | {max(peaks):.3f} | "
            f"{statistics.median(runtimes):.2f} | {sum(runtimes) / 3600:.3f} |"
            if observed else f"| {analysis} | 0 | NA | NA | NA |"
        )
    if len(lines) == 2:
        lines.append("| No measured unit yet | 0 | NA | NA | NA |")
    return lines


def _projected_runtime(rows: list[dict[str, str]]) -> str:
    discovery = [float(row["runtime_sec"]) for row in rows if row["analysis"] == "LAVA_DISCOVERY" and row["exit_status"] == "0"]
    conditional = [float(row["runtime_sec"]) for row in rows if row["analysis"] == "LAVA_CONDITIONAL" and row["exit_status"] == "0"]
    components = []
    if discovery:
        components.append(f"discovery ~{statistics.median(discovery) * 2495 / 3600:.2f} h")
    if conditional:
        components.append(f"conditional upper bound ~{statistics.median(conditional) * 2495 / 3600:.2f} h")
    return "; ".join(components) if components else "unmeasured"


def build_report(fingerprint: str | None = None) -> None:
    rows: list[dict[str, str]] = []
    if BENCHMARK.exists():
        fields, rows = read_tsv(BENCHMARK)
        if fields != BENCHMARK_FIELDS:
            raise ExecutionError("RAM benchmark schema drifted")
    memory, reserve = physical_memory_bytes(), memory_safety_reserve_bytes()
    provenance_state, measured_peaks = "PENDING", []
    if fingerprint and selection_path(fingerprint, "discovery").exists():
        try:
            payload = benchmark_provenance_payload(fingerprint)
            provenance_state = str(payload["state"])
            measured_peaks.extend(
                float(item["peak_rss_bytes"]) / 1024**3
                for item in payload["measured_admission_units"]
            )
            measured_peaks.extend(
                float(item["peak_ram_gb"]) for item in payload["all_lava_attempt_rows"]
            )
        except ExecutionError:
            provenance_state = "INVALID_OR_INCOMPLETE"
    maximum = max(measured_peaks, default=None)
    if provenance_state == "ALL_REQUIRED_PHASE_REPRESENTATIVES_COMPLETE" and maximum is not None:
        feasibility = (
            "PASS — discovery and every available conditional representative fit with the predeclared reserve."
            if maximum * 1024**3 + reserve <= memory else
            "BLOCKED_BY_MEASURED_PER_LOCUS_RAM — an intact representative plus reserve exceeds physical RAM."
        )
    elif provenance_state.startswith("DISCOVERY_REPRESENTATIVES_COMPLETE"):
        feasibility = "PROVISIONAL DISCOVERY PASS — conditional representative measurement remains pending."
    else:
        feasibility = "PENDING — fingerprint-bound representative measurements are incomplete."
    stats = "\n".join(report_statistics(rows))
    maximum_text = f"{maximum:.3f} GiB" if maximum is not None else "NA"
    report = f"""# Track B RAM-aware execution report

This is an operational execution report, not a scientific result. Physical memory visible to the supervisor is {memory} bytes ({memory / 1024**3:.3f} GiB). The predeclared non-worker reserve is {reserve} bytes ({reserve / 1024**3:.3f} GiB). Measurements are OS-reported peak RSS around one fresh analysis process; no 16-GiB monolithic assumption is used.

## Scientifically equivalent decomposition decisions

- **LAVA — safely decomposable by complete predeclared locus.** Each discovery worker loads one exact chromosome summary-statistic shard and the complete reference chromosome, constructs one intact official locus, runs all eight local univariate tests, and runs every locally eligible frozen pair. All 2,495 locus results are collated before the unchanged family-wide BH corrections. The conditional pass starts only after that full discovery family is frozen and runs BMI-only, sleep-apnea-only, and MDD-only models separately. Failed and unstable loci remain explicit; no skipped locus is called null.
- **PLACO+ — the global nuisance/correlation fit is not decomposable.** It must use all valid genome-wide variants. Only single-variant testing after those global parameters are frozen may be chunked, and only under the exact same equations and final global multiple-testing family.
- **Fine-mapping and coloc — safely decomposable by validated whole locus, never within an LD block.** Every unit retains the complete frozen SNP set and full ancestry-matched LD matrix.
- **Indivisible units.** One complete LAVA locus, the PLACO+ genome-wide nuisance fit, and one whole fine-map/coloc LD locus are the smallest permitted memory units.

## Measured resource results

{stats}

Current LAVA assessment: **{feasibility}** Fingerprint-bound benchmark state: `{provenance_state}`. Maximum observed admission-relevant peak (all constructed units plus every attempted frozen benchmark candidate): {maximum_text}. Expected sequential runtime from currently completed units: {_projected_runtime(rows)}.

The benchmark selection is deterministic and frozen before each phase: discovery uses only sealed reference chromosome/locus SNP burden; conditional benchmarking is restricted only by the already-predeclared full-family FDR and conditioner-h2 eligibility gates, then uses the same input-burden ranking. Scientific failures are preserved but cannot serve as a RAM-pass representative. Every worker has an immutable receipt/log/result bundle committed by one atomic ready link.

## Remaining larger-server decisions

- LAVA needs a larger server only if a complete intact locus plus the reserve actually exceeds this machine, or if a non-locus aggregation/finalization unit cannot be made equivalently streaming. It is not blocked by an assumed monolithic requirement.
- PLACO+ remains dependent on the measured global all-variant nuisance/correlation fit; later exact variant shards do not replace that global fit.
- Fine-mapping/coloc requires a larger server only for a complete indivisible locus whose measured peak plus reserve exceeds this machine.

Machine-readable measurements are in `results/track_b/RAM_BENCHMARK.tsv`; their current fingerprint binding is in `results/track_b/RAM_BENCHMARK.provenance.json`.
"""
    atomic_text(REPORT, report)


def ensure_benchmark_file(fingerprint: str | None = None) -> None:
    if not BENCHMARK.exists():
        atomic_text(BENCHMARK, tsv_text(BENCHMARK_FIELDS, []))
    if fingerprint is not None:
        assert_benchmark_namespace(fingerprint)
    build_report(fingerprint)


def _completed_receipts(fingerprint: str) -> list[dict[str, Any]]:
    receipts = []
    for phase in ("discovery", "aggregate-discovery", "conditional", "finalize"):
        indices: list[int | None] = list(range(1, 2496)) if phase in LOCUS_PHASES else [None]
        for index in indices:
            receipt = validate_bundle(phase, index, fingerprint)
            if receipt is not None:
                receipts.append(receipt)
    return receipts


def execution_evidence_payload(fingerprint: str) -> dict[str, Any]:
    validate_benchmark_provenance(fingerprint)
    result_size, result_digest = file_identity(RESULT_LOCK)
    benchmark_size, benchmark_digest = file_identity(BENCHMARK)
    report_size, report_digest = file_identity(REPORT)
    provenance_size, provenance_digest = file_identity(BENCHMARK_PROVENANCE)
    namespace_size, namespace_digest = file_identity(BENCHMARK_NAMESPACE)
    receipts = [_receipt_evidence(receipt) for receipt in _completed_receipts(fingerprint)]
    if len([item for item in receipts if item["phase"] == "discovery"]) != 2495:
        raise ExecutionError("execution evidence lacks the complete 2,495-locus discovery family")
    if len([item for item in receipts if item["phase"] == "conditional"]) != 2495:
        raise ExecutionError("execution evidence lacks the complete 2,495-locus conditional family")
    failures = []
    failure_dir = run_root(fingerprint) / "failed_attempts"
    if failure_dir.is_dir():
        for path in sorted(failure_dir.iterdir()):
            if path.is_file():
                size, digest = file_identity(path)
                failures.append({"path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest})
    selections = []
    for phase in ("discovery", "conditional"):
        path = selection_path(fingerprint, phase)
        if path.is_file():
            size, digest = file_identity(path)
            selections.append({"path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest})
    return {
        "schema_version": 2, "analysis_id": "track-b-v1.0-local",
        "execution_fingerprint": fingerprint,
        "result_lock": {"path": str(RESULT_LOCK.relative_to(ROOT)), "bytes": result_size, "sha256": result_digest},
        "ram_benchmark": {"path": str(BENCHMARK.relative_to(ROOT)), "bytes": benchmark_size, "sha256": benchmark_digest},
        "ram_report": {"path": str(REPORT.relative_to(ROOT)), "bytes": report_size, "sha256": report_digest},
        "ram_benchmark_provenance": {
            "path": str(BENCHMARK_PROVENANCE.relative_to(ROOT)),
            "bytes": provenance_size, "sha256": provenance_digest,
        },
        "ram_benchmark_namespace": {
            "path": str(BENCHMARK_NAMESPACE.relative_to(ROOT)),
            "bytes": namespace_size, "sha256": namespace_digest,
        },
        "selections": selections, "worker_receipts": receipts, "failed_attempts": failures,
    }


def seal_execution_evidence(fingerprint: str) -> None:
    expected = execution_evidence_payload(fingerprint)
    if EXECUTION_EVIDENCE.exists():
        try:
            observed = json.loads(EXECUTION_EVIDENCE.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ExecutionError(f"execution-evidence seal is unreadable: {error}") from error
        if observed != expected:
            raise ExecutionError("execution-evidence seal differs from live completed run")
        return
    publish_no_replace(EXECUTION_EVIDENCE, (json.dumps(expected, indent=2, sort_keys=True) + "\n").encode())


def verify_execution_evidence() -> None:
    try:
        observed = json.loads(EXECUTION_EVIDENCE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ExecutionError(f"execution-evidence seal is unreadable: {error}") from error
    fingerprint = str(observed.get("execution_fingerprint", "")) if isinstance(observed, dict) else ""
    if not FINGERPRINT_RE.fullmatch(fingerprint):
        raise ExecutionError("execution-evidence seal lacks a valid execution fingerprint")
    quick_execution_fingerprint(fingerprint)
    if observed != execution_evidence_payload(fingerprint):
        raise ExecutionError("execution-evidence seal differs from live completed run")


def publish_canonical_results() -> None:
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "--publish-staged-results"],
        cwd=ROOT, text=True, capture_output=True, check=False,
    )
    if result.returncode != 0:
        raise ExecutionError(
            "validated LAVA staging could not be transactionally published:\n"
            + result.stdout + result.stderr
        )
    if not RESULT_LOCK.is_file():
        raise ExecutionError("LAVA publication returned success without a result provenance lock")


def run_all(fingerprint: str) -> None:
    run_discovery_benchmarks(fingerprint)
    for index in range(1, 2496):
        run_measured("discovery", index, fingerprint)
        if index % 25 == 0 or index == 2495:
            print(f"LAVA_DISCOVERY_PROGRESS completed_through={index} total=2495")
    run_measured("aggregate-discovery", None, fingerprint)
    conditional_plan = conditional_candidate_plan(fingerprint)
    if conditional_plan:
        _run_candidate_plan(fingerprint, "conditional", conditional_plan)
    else:
        write_benchmark_provenance(fingerprint)
    for index in range(1, 2496):
        run_measured("conditional", index, fingerprint)
        if index % 100 == 0 or index == 2495:
            print(f"LAVA_CONDITIONAL_PROGRESS completed_through={index} total=2495")
    run_measured("finalize", None, fingerprint)
    publish_canonical_results()
    write_benchmark_provenance(fingerprint)
    build_report(fingerprint)
    seal_execution_evidence(fingerprint)
    print("TRACK_B_LAVA_SEQUENTIAL_COMPLETE loci=2495 phases=discovery,conditional")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--benchmark", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--report-only", action="store_true")
    mode.add_argument("--verify-execution-evidence", action="store_true")
    mode.add_argument("--internal-measure", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--phase", choices=sorted(ALL_PHASES))
    parser.add_argument("--locus-index", type=int)
    parser.add_argument("--fingerprint")
    parser.add_argument("--worker-output", help=argparse.SUPPRESS)
    parser.add_argument("--internal-result", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.internal_measure:
        if (
            not args.phase or not args.fingerprint or not args.internal_result or not args.worker_output
            or not FINGERPRINT_RE.fullmatch(args.fingerprint)
        ):
            parser.error("internal measurement requires phase, fingerprint, worker output, and result path")
        raise SystemExit(internal_measure(args))
    if args.verify_execution_evidence:
        verify_execution_evidence()
        print("TRACK_B_LAVA_EXECUTION_EVIDENCE_VALIDATED")
        return
    if args.report_only:
        ensure_benchmark_file()
        fingerprint = args.fingerprint if args.fingerprint and FINGERPRINT_RE.fullmatch(args.fingerprint) else None
        build_report(fingerprint)
        print(f"wrote {REPORT.relative_to(ROOT)}")
        return
    fingerprint = preflight()
    ensure_benchmark_file(fingerprint)
    if args.benchmark:
        run_discovery_benchmarks(fingerprint)
        write_benchmark_provenance(fingerprint)
        build_report(fingerprint)
        print("TRACK_B_LAVA_DISCOVERY_RAM_BENCHMARK_PASS representatives=3 conditional=PENDING_FULL_DISCOVERY_BH")
    else:
        run_all(fingerprint)


if __name__ == "__main__":
    try:
        main()
    except ExecutionError as error:
        raise SystemExit(f"ERROR: {error}") from error
