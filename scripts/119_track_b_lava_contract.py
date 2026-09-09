#!/usr/bin/env python3
"""Verify and seal the dedicated Track B LAVA execution contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lava_contract


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "config/track_b_local_analysis_policy.json"
INPUT_LOCK = ROOT / "results/track_b/local_analysis_input.lock.json"
PAIR_MANIFEST = ROOT / "results/track_b/pair_manifest.tsv"
INPUT_INFO = ROOT / "results/track_b/lava_input_info.tsv"
INPUT_PROVENANCE = ROOT / "results/track_b/lava_input_provenance.tsv"
LOCUS_FILE = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
REFERENCE = ROOT / "ref/lava/ukb_v1.1/reference.provenance.json"
CHROMOSOME_INPUT_LOCK = ROOT / "results/track_b/lava_chromosome_inputs.provenance.json"
CHROMOSOME_INPUT_BUILDER = ROOT / "scripts/129_prepare_track_b_lava_chromosome_inputs.py"
CHECKPOINT_DIR = ROOT / "results/track_b/checkpoints/lava"
RESULT_DIR = ROOT / "results/track_b/local"
LOW_LEVEL_RESULTS = [
    RESULT_DIR / "lava_locus_status.tsv", RESULT_DIR / "lava_univariate.tsv",
    RESULT_DIR / "lava_bivariate.tsv", RESULT_DIR / "lava_conditional.tsv",
]
SUMMARY_RESULTS = [
    ROOT / "results/track_b/04_lava_local_results.tsv",
    ROOT / "results/track_b/05_local_conditional_results.tsv",
]
RESULTS = LOW_LEVEL_RESULTS + SUMMARY_RESULTS
RESULT_LOCK = RESULT_DIR / "lava_results.provenance.json"
STAGING_ROOT = ROOT / "results/track_b/staging"
R_SCRIPT = ROOT / ".r-env/bin/Rscript"
RUNNER = ROOT / "scripts/120_run_track_b_lava.R"
SUPERVISOR = ROOT / "scripts/131_run_track_b_lava_sequential.py"
VALIDATOR = ROOT / "scripts/121_validate_track_b_lava.py"
CHECKPOINT_VALIDATOR = ROOT / "scripts/132_validate_track_b_lava_checkpoint.R"
RUNTIME_VALIDATOR = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
RUNTIME_PACKAGES = ("LAVA", "matrixsampling", "cpp11", "keep")


def staging_dir(fingerprint: str) -> Path:
    if len(fingerprint) != 64 or any(character not in "0123456789abcdef" for character in fingerprint):
        fail("invalid Track B LAVA staging fingerprint")
    return STAGING_ROOT / f"lava_{fingerprint}"


def staged_results(fingerprint: str) -> list[Path]:
    # The complete staged family is part of the immutable finalize READY bundle.
    # This makes the publication gap crash-resumable: no unreceipted external
    # staging directory can disappear after finalization succeeds.
    directory = CHECKPOINT_DIR / fingerprint / "finalize" / "unit_all"
    return [directory / path.name for path in RESULTS]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def file_identity(path: Path) -> tuple[int, str]:
    try:
        label = path.relative_to(ROOT)
    except ValueError:
        label = path
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size == 0:
                fail(f"missing real non-empty artifact: {label}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        fail(f"could not hash artifact {label}: {error}")

    def identity_fields(value: os.stat_result) -> tuple[int, int, int, int, int]:
        return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns

    if identity_fields(before) != identity_fields(after) or identity_fields(after) != identity_fields(current):
        fail(f"artifact changed while hashing: {label}")
    return after.st_size, digest.hexdigest()


def sha256(path: Path) -> str:
    return file_identity(path)[1]


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def load_policy() -> dict[str, object]:
    try:
        policy = json.loads(POLICY.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B local policy is unreadable: {error}")
    if (
        policy.get("analysis_id") != "track-b-v1.0-local"
        or policy.get("lava_version") != "0.1.5"
        or policy.get("lava_commit") != "e729a245f7b6923967a96804fbf5246eadf2d6c6"
        or policy.get("lava_source_archive_sha256") != "15477c9547c3533d681cbfb6491508a29342fe99b725845054887b9d29dc2bec"
        or policy.get("lava_process_locus_nref_correction") is not False
        or policy.get("expected_loci") != 2495
        or policy.get("expected_analysis_traits") != 8
        or policy.get("expected_discovery_pairs") != 3
        or policy.get("planned_univariate_tests") != 2495 * 8
        or policy.get("planned_bivariate_pair_locus_family_max") != 2495 * 3
        or policy.get("ram_aware_execution", {}).get("execution_unit") != "WHOLE_PREDECLARED_LAVA_LOCUS"
        or policy.get("ram_aware_execution", {}).get("maximum_loci_per_process") != 1
        or policy.get("ram_aware_execution", {}).get("worker_process_rule") != "FRESH_R_PROCESS_PER_LOCUS_WITH_EXPLICIT_OBJECT_REMOVAL_AND_FULL_GARBAGE_COLLECTION"
        or policy.get("locus_definition", {}).get("path") != str(LOCUS_FILE.relative_to(ROOT))
        or policy.get("locus_definition", {}).get("sha256") != "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882"
    ):
        fail("Track B local policy scope drifted")
    return policy


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path.relative_to(ROOT)}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_inputs(
    policy: dict[str, object], *, rehash_sumstats: bool = False,
) -> tuple[dict[str, object], dict[str, object]]:
    try:
        lock = json.loads(INPUT_LOCK.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B local input lock is unreadable: {error}")
    if (
        lock.get("analysis_id") != policy["analysis_id"]
        or lock.get("policy_sha256") != sha256(POLICY)
        or lock.get("local_results_accessed_before_input_freeze") is not False
    ):
        fail("Track B local input lock differs from policy or pre-result state")
    hashes = lock.get("input_artifact_sha256", {})
    expected_artifacts = {
        "results/track_b/lava_input_info.tsv",
        "results/track_b/lava_input_provenance.tsv",
        "results/track_b/lava_pair_manifest.tsv",
        "results/track_b/lava_runtime_policy.tsv",
        "results/track_b/lava_sample_overlap.txt",
        "results/track_b/local_conditional_manifest.tsv",
        "results/track_b/local_method_robustness_plan.tsv",
    }
    if not isinstance(hashes, dict) or set(hashes) != expected_artifacts:
        fail("Track B local input artifact family is incomplete")
    for relative, expected in hashes.items():
        if sha256(ROOT / relative) != expected:
            fail(f"Track B local input differs from lock: {relative}")
    if sha256(PAIR_MANIFEST) != lock.get("pair_manifest_sha256"):
        fail("frozen Track B pair manifest differs from local input lock")
    locus_sha256 = sha256(LOCUS_FILE)
    if (
        locus_sha256 != policy["locus_definition"]["sha256"]
        or lock.get("locus_file_sha256") != policy["locus_definition"]["sha256"]
    ):
        fail("live 2,495-locus definition, policy pin, and local input lock differ")

    input_rows = read_tsv(INPUT_INFO)
    provenance_rows = read_tsv(INPUT_PROVENANCE)
    trait_order = [str(value) for value in policy["trait_order"]]
    if (
        [row.get("phenotype") for row in input_rows] != trait_order
        or [row.get("phenotype") for row in provenance_rows] != trait_order
    ):
        fail("live LAVA sumstats family differs from the frozen trait order")
    info_paths = [row.get("filename", "") for row in input_rows]
    provenance_paths = [row.get("filename", "") for row in provenance_rows]
    if info_paths != provenance_paths or len(set(info_paths)) != len(info_paths):
        fail("LAVA input-info and provenance paths do not identify one exact file per trait")

    live_sumstats: dict[str, str] = {}
    for row in provenance_rows:
        relative = row.get("filename", "")
        path = ROOT / relative
        try:
            expected_bytes = int(row.get("bytes", ""))
        except ValueError:
            fail(f"invalid frozen byte count for LAVA input: {relative}")
        try:
            observed_bytes = path.stat().st_size
        except OSError as error:
            fail(f"live LAVA input is unreadable: {relative}: {error}")
        if not path.is_file() or path.is_symlink() or observed_bytes != expected_bytes:
            fail(f"live LAVA input byte count differs from provenance: {relative}")
        expected_sha256 = str(row.get("sha256", ""))
        if len(expected_sha256) != 64:
            fail(f"invalid frozen LAVA input SHA-256: {relative}")
        if rehash_sumstats and sha256(path) != expected_sha256:
            fail(f"live LAVA input differs from provenance: {relative}")
        live_sumstats[relative] = expected_sha256
    return lock, {"locus_file_sha256": locus_sha256, "sumstats_sha256": live_sumstats}


def validate_reference(*, rehash: bool = False) -> dict[str, object]:
    generic_path, generic_policy = lava_contract.load_policy(ROOT)
    del generic_path
    return lava_contract.validate_reference(ROOT, generic_policy, rehash=rehash)


def validate_chromosome_input_lock(policy: dict[str, object]) -> dict[str, object]:
    try:
        lock = json.loads(CHROMOSOME_INPUT_LOCK.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B LAVA chromosome-input lock is unreadable: {error}")
    bindings = lock.get("bindings", {})
    totals = lock.get("totals", {})
    if (
        lock.get("schema_version") != "track-b-lava-chromosome-inputs.1"
        or lock.get("analysis_id") != policy["analysis_id"]
        or lock.get("chromosomes") != list(range(1, 23))
        or lock.get("trait_order") != policy["trait_order"]
        or lock.get("scientific_filtering") != "NONE_ONLY_EXACT_SEALED_REFERENCE_SNP_INTERSECTION"
        or bindings.get("policy_sha256") != sha256(POLICY)
        or bindings.get("local_input_lock_sha256") != sha256(INPUT_LOCK)
        or bindings.get("reference_provenance_sha256") != sha256(REFERENCE)
        or bindings.get("script_sha256") != sha256(CHROMOSOME_INPUT_BUILDER)
        or totals.get("chromosome_count") != 22
        or totals.get("trait_count") != 8
        or totals.get("shard_count") != 176
        or totals.get("chromosome_input_info_count") != 22
    ):
        fail("Track B LAVA chromosome-input lock scope or bindings drifted")
    shards = lock.get("shards")
    input_info = lock.get("chromosome_input_info")
    if not isinstance(shards, list) or len(shards) != 176 or not isinstance(input_info, list) or len(input_info) != 22:
        fail("Track B LAVA chromosome-input artifact family is incomplete")
    for record in shards + input_info:
        if not isinstance(record, dict):
            fail("malformed Track B LAVA chromosome-input record")
        relative = Path(str(record.get("path", "")))
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            fail("unsafe Track B LAVA chromosome-input path")
        path = ROOT / relative
        try:
            size = int(record.get("bytes", -1))
        except (TypeError, ValueError):
            fail("invalid Track B LAVA chromosome-input byte count")
        if path.is_symlink() or not path.is_file() or path.stat().st_size != size:
            fail(f"Track B LAVA chromosome input is absent or has the wrong size: {relative}")
    return lock


def fully_verify_chromosome_inputs() -> None:
    result = subprocess.run(
        [sys.executable, str(CHROMOSOME_INPUT_BUILDER), "--verify", "--root", str(ROOT)],
        cwd=ROOT, text=True, capture_output=True, check=False,
    )
    if result.returncode != 0:
        fail("Track B LAVA chromosome inputs failed full verification: " + result.stdout + result.stderr)


def load_supervisor():
    path = SUPERVISOR
    specification = importlib.util.spec_from_file_location("track_b_lava_supervisor_contract", path)
    if specification is None or specification.loader is None:
        fail("could not load Track B LAVA sequential supervisor")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def runtime_contract_payload() -> dict[str, object]:
    if not R_SCRIPT.is_file() or not RUNTIME_VALIDATOR.is_file():
        fail("pinned Track B LAVA R runtime contract is missing")
    result = subprocess.run(
        [str(R_SCRIPT), str(RUNTIME_VALIDATOR)], cwd=ROOT,
        text=True, capture_output=True, check=False,
    )
    output = result.stdout + result.stderr
    markers = [line for line in output.splitlines() if line.startswith("TRACK_B_LAVA_RUNTIME\t")]
    if result.returncode != 0 or len(markers) != 1:
        fail("Track B LAVA runtime semantic verification failed: " + output)
    fields: dict[str, str] = {}
    for item in markers[0].split("\t")[1:]:
        if "=" not in item:
            fail("malformed Track B LAVA runtime marker")
        key, value = item.split("=", 1)
        fields[key] = value
    expected = {
        "status": "PASS", "r_version": "4.3.3", "lava_version": "0.1.5",
        "lava_commit": "e729a245f7b6923967a96804fbf5246eadf2d6c6",
        "process_locus_nref_correction": "FALSE",
        "functions": "process.input;read.loci;process.locus;run.univ;run.bivar;run.pcor",
    }
    if any(fields.get(key) != value for key, value in expected.items()) or not fields.get("platform"):
        fail("Track B LAVA runtime marker differs from the pinned implementation")
    package_files: list[dict[str, object]] = []
    library = ROOT / ".r-env/lib/R/library"
    for package in RUNTIME_PACKAGES:
        directory = library / package
        if not directory.is_dir() or directory.is_symlink():
            fail(f"pinned R package directory is missing or unsafe: {package}")
        for path in sorted(item for item in directory.rglob("*") if item.is_file()):
            if path.is_symlink():
                fail(f"symbolic links are forbidden in pinned R package: {path}")
            size, digest = file_identity(path)
            package_files.append({
                "path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest,
            })
    r_size, r_digest = file_identity(R_SCRIPT.resolve(strict=True))
    return {
        "r_marker": fields,
        "rscript": {"path": str(R_SCRIPT.resolve(strict=True)), "bytes": r_size, "sha256": r_digest},
        "package_files": package_files,
    }


def run_fingerprint(*, rehash_sumstats: bool = False) -> str:
    policy = load_policy()
    _, live_inputs = validate_inputs(policy, rehash_sumstats=rehash_sumstats)
    validate_reference(rehash=False)
    validate_chromosome_input_lock(policy)
    payload = {
        "analysis_id": policy["analysis_id"], "policy_sha256": sha256(POLICY),
        "input_lock_sha256": sha256(INPUT_LOCK), "reference_sha256": sha256(REFERENCE),
        "chromosome_input_lock_sha256": sha256(CHROMOSOME_INPUT_LOCK),
        "chromosome_input_builder_sha256": sha256(CHROMOSOME_INPUT_BUILDER),
        "runner_sha256": sha256(RUNNER),
        "sequential_supervisor_sha256": sha256(SUPERVISOR),
        "validator_sha256": sha256(VALIDATOR),
        "checkpoint_validator_sha256": sha256(CHECKPOINT_VALIDATOR),
        "runtime_validator_sha256": sha256(RUNTIME_VALIDATOR),
        "runtime_setup_sha256": sha256(ROOT / "scripts/30_setup_lava.sh"),
        "runtime_contract": runtime_contract_payload(),
        "contract_sha256": sha256(Path(__file__)),
        "reference_contract_sha256": sha256(ROOT / "scripts/lava_contract.py"),
        "live_inputs": live_inputs,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def result_payload() -> dict[str, object]:
    policy = load_policy()
    fingerprint = run_fingerprint()
    supervisor = load_supervisor()
    checkpoints = []
    for phase in ("discovery", "aggregate-discovery", "conditional", "finalize"):
        indices = range(1, int(policy["expected_loci"]) + 1) if phase in {"discovery", "conditional"} else (None,)
        for index in indices:
            receipt = supervisor.validate_bundle(phase, index, fingerprint)
            if receipt is None:
                fail(f"missing ready Track B LAVA checkpoint bundle: {phase}/{index}")
            bundle = supervisor.bundle_path(phase, index, fingerprint)
            target = supervisor._safe_attempt_target(bundle)
            ready_size, ready_digest = file_identity(target / "READY")
            receipt_size, receipt_digest = file_identity(target / "receipt.json")
            checkpoints.append({
                "phase": phase,
                "locus_index": index,
                "bundle_path": str(bundle.relative_to(ROOT)),
                "ready_bytes": ready_size,
                "ready_sha256": ready_digest,
                "receipt_bytes": receipt_size,
                "receipt_sha256": receipt_digest,
                "peak_rss_bytes": receipt["peak_rss_bytes"],
                "runtime_sec": receipt["runtime_sec"],
                "qc": receipt["marker"]["qc"],
            })
    results = []
    for path in RESULTS:
        size, digest = file_identity(path)
        results.append({"path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest})
    return {
        "schema_version": 4, "analysis_id": policy["analysis_id"],
        "lava_version": policy["lava_version"], "run_fingerprint": fingerprint,
        "input_lock_sha256": sha256(INPUT_LOCK), "reference_sha256": sha256(REFERENCE),
        "chromosome_input_lock_sha256": sha256(CHROMOSOME_INPUT_LOCK),
        "checkpoint_locus_count": int(policy["expected_loci"]),
        "checkpoint_artifact_count": len(checkpoints), "checkpoints": checkpoints,
        "results": results,
    }


def exclusive_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            fail("Track B LAVA result provenance already exists; overwrite is forbidden")
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)
        fsync_directory(path.parent)


def seal_results(expected_result_identities: dict[str, tuple[int, str]] | None = None) -> None:
    expected = result_payload()
    if expected_result_identities is not None:
        observed_result_identities = {
            str(item["path"]): (int(item["bytes"]), str(item["sha256"]))
            for item in expected["results"]
        }
        if observed_result_identities != expected_result_identities:
            fail("canonical Track B LAVA results differ from the semantically validated staged identities")
    if RESULT_LOCK.is_file():
        validate_results()
        return
    exclusive_json(RESULT_LOCK, expected)
    print(
        "TRACK_B_LAVA_RESULTS_SEALED "
        f"loci={expected['checkpoint_locus_count']} artifacts={expected['checkpoint_artifact_count']}"
    )


def validate_results() -> None:
    try:
        observed = json.loads(RESULT_LOCK.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"Track B LAVA result provenance is unreadable: {error}")
    if observed != result_payload():
        fail("Track B LAVA results differ from immutable provenance")
    print("TRACK_B_LAVA_RESULTS_VALIDATED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--preflight", action="store_true")
    group.add_argument("--execution-fingerprint", action="store_true")
    group.add_argument("--quick-execution-fingerprint", action="store_true")
    group.add_argument("--run-fingerprint", action="store_true")
    group.add_argument("--verify-results", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        policy = load_policy()
        validate_inputs(policy, rehash_sumstats=True)
        validate_reference(rehash=True)
        fully_verify_chromosome_inputs()
        runtime_contract_payload()
        fingerprint = run_fingerprint()
        print(
            f"TRACK_B_LAVA_PREFLIGHT_PASS fingerprint={fingerprint} "
            "execution=ONE_COMPLETE_LOCUS_PER_FRESH_R_PROCESS ram_gate=MEASURED_PEAK_RSS"
        )
    elif args.execution_fingerprint or args.quick_execution_fingerprint:
        print(run_fingerprint())
    elif args.run_fingerprint:
        print(run_fingerprint())
    else:
        validate_results()


if __name__ == "__main__":
    main()
