#!/usr/bin/env python3
"""Build and verify the live dense-input and method-readiness gates for Track B."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_SPEC = importlib.util.spec_from_file_location(
    "track_b_pleiotropy_contract", ROOT / "scripts/123_track_b_pleiotropy_contract.py"
)
if CONTRACT_SPEC is None or CONTRACT_SPEC.loader is None:
    raise RuntimeError("could not load Track B pleiotropy contract")
contract = importlib.util.module_from_spec(CONTRACT_SPEC)
CONTRACT_SPEC.loader.exec_module(contract)

INPUT_GATE = Path("results/track_b/pleiotropy/input_gate.tsv")
READINESS_GATE = Path("results/track_b/pleiotropy/readiness_gate.tsv")
INPUT_LOCK = Path("results/track_b/pleiotropy/input_gate.lock.json")

INPUT_FIELDS = [
    "pair_id", "family_role", "trait1", "trait2", "ancestry", "analysis_build",
    "trait1_file", "trait1_bytes", "trait1_sha256", "trait1_rows", "trait1_schema",
    "trait1_INFO_status", "trait2_file", "trait2_bytes", "trait2_sha256",
    "trait2_rows", "trait2_schema", "trait2_INFO_status", "minimum_rows_per_trait",
    "full_genome_wide_dense_required", "hapmap3_only_forbidden", "gzip_integrity",
    "minimum_aligned_eligible_variants", "minimum_fraction_of_smaller_dense_input",
    "required_autosomes", "per_autosome_provenance_required",
    "pair_alignment_status", "input_gate_status", "blocker",
]
READINESS_FIELDS = [
    "component_id", "pair_scope", "dense_input_gate", "replication_gate",
    "local_analysis_gate", "implementation_gate", "software_gate", "reference_gate",
    "compute_gate", "scan_gate", "locus_publication_gate", "overall_status",
    "blockers", "blocked_result_semantics", "claim_status",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def tsv_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def full_gzip_scan(path: Path, expected_schema: list[str]) -> dict[str, object]:
    """Checksum the compressed payload and fully decompress it to verify row count/CRC."""
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty dense summary statistics: {path}")
    compressed_sha256 = contract.sha256(path)
    try:
        with gzip.open(path, "rb") as handle:
            header_raw = handle.readline()
            try:
                header = header_raw.rstrip(b"\r\n").decode("utf-8").split("\t")
            except UnicodeDecodeError as error:
                fail(f"dense header is not UTF-8: {path}: {error}")
            if header != expected_schema:
                fail(f"dense schema drifted for {path}: observed={header} expected={expected_schema}")
            rows = 0
            payload_bytes = 0
            last = b""
            while True:
                block = handle.read(8 * 1024 * 1024)
                if not block:
                    break
                payload_bytes += len(block)
                rows += block.count(b"\n")
                last = block[-1:]
            if payload_bytes and last != b"\n":
                rows += 1
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        fail(f"dense gzip integrity scan failed for {path}: {error}")
    if rows <= 0:
        fail(f"dense summary statistics have no data rows: {path}")
    return {
        "path": str(path), "bytes": path.stat().st_size,
        "sha256": compressed_sha256, "rows": rows, "schema": ",".join(header),
        "gzip_integrity": "FULL_DECOMPRESSION_CRC_PASS",
    }


def validate_live_dense_inputs(
    root: Path, policy: dict[str, Any], upstream: dict[str, Any],
) -> dict[str, dict[str, object]]:
    expected_schema = [str(value) for value in policy["dense_input_contract"]["required_columns"]]
    scans: dict[str, dict[str, object]] = {}
    for record in upstream["dense_traits"]:
        trait = str(record["trait_id"])
        relative = Path(str(record["path"]))
        path = contract.root_path(root, relative)
        scan = full_gzip_scan(path, expected_schema)
        if (
            scan["bytes"] != record["bytes"]
            or scan["sha256"] != record["sha256"]
            or scan["rows"] != record["rows"]
            or scan["schema"] != record["schema"]
        ):
            fail(f"live dense identity/count/schema differs from the frozen gate for {trait}")
        scans[trait] = {
            **scan, "path": str(relative),
        }
    if list(scans) != contract.EXPECTED_TRAITS:
        fail("live dense scan family differs from the five frozen traits")
    return scans


def pair_gate_rows(
    root: Path, policy: dict[str, Any], scans: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    _, dense_rows = contract.read_tsv(root / contract.DENSE_QC)
    dense = {row["trait_id"]: row for row in dense_rows}
    rows: list[dict[str, object]] = []
    plausibility = policy["dense_input_contract"]["pair_materialization_plausibility"]
    for pair in policy["pairs"]:
        trait1, trait2 = pair["trait1"], pair["trait2"]
        first, second = scans[trait1], scans[trait2]
        rows.append({
            "pair_id": pair["pair_id"], "family_role": pair["family_role"],
            "trait1": trait1, "trait2": trait2, "ancestry": policy["ancestry"],
            "analysis_build": policy["analysis_build"],
            "trait1_file": first["path"], "trait1_bytes": first["bytes"],
            "trait1_sha256": first["sha256"], "trait1_rows": first["rows"],
            "trait1_schema": first["schema"], "trait1_INFO_status": dense[trait1]["INFO_status"],
            "trait2_file": second["path"], "trait2_bytes": second["bytes"],
            "trait2_sha256": second["sha256"], "trait2_rows": second["rows"],
            "trait2_schema": second["schema"], "trait2_INFO_status": dense[trait2]["INFO_status"],
            "minimum_rows_per_trait": policy["dense_input_contract"]["minimum_variant_rows_per_trait"],
            "full_genome_wide_dense_required": "TRUE", "hapmap3_only_forbidden": "TRUE",
            "gzip_integrity": "FULL_DECOMPRESSION_CRC_PASS_BOTH",
            "minimum_aligned_eligible_variants": plausibility["minimum_aligned_eligible_variants"],
            "minimum_fraction_of_smaller_dense_input": plausibility["minimum_fraction_of_smaller_dense_input"],
            "required_autosomes": ",".join(str(value) for value in plausibility["required_autosomes"]),
            "per_autosome_provenance_required": "TRUE",
            "pair_alignment_status": "NOT_RUN_PRE_RESULT;EXACT_GLOBAL_AND_22_AUTOSOME_COUNTS_REQUIRED_AT_MATERIALIZATION",
            "input_gate_status": "READY_DENSE_INPUTS", "blocker": "NONE",
        })
    return rows


def primary_replication_terminal(root: Path) -> bool:
    path = root / "results/track_b/02_independent_global_replication.tsv"
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        _, rows = contract.read_tsv(path)
    except SystemExit:
        return False
    by_pair = {row.get("pair_id"): row for row in rows}
    return (
        set(by_pair) == {"A", "B"}
        and by_pair["A"].get("replication_status") == "NO_VALID_REPLICATION"
        and by_pair["B"].get("replication_status")
        in {"DIRECTIONAL_REPLICATION", "STRONG_REPLICATION", "FAILED_REPLICATION"}
    )


def local_analysis_terminal(root: Path) -> bool:
    path = root / "results/track_b/local/lava_results.provenance.json"
    return path.is_file() and path.stat().st_size > 0


def exact_file_ready(path: Path, *, size: int | None = None, digest: str | None = None) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if size is not None and path.stat().st_size != size:
        return False
    return digest is None or contract.sha256(path) == digest


def clean_git_commit(path: Path, expected: str) -> bool:
    if not (path / ".git").is_dir():
        return False
    commit = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=False, capture_output=True, text=True,
    )
    status = subprocess.run(
        ["git", "-C", str(path), "status", "--porcelain"],
        check=False, capture_output=True, text=True,
    )
    return (
        commit.returncode == 0 and status.returncode == 0
        and commit.stdout.strip() == expected and not status.stdout.strip()
    )


def physical_memory_bytes() -> int:
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (OSError, ValueError) as error:
        fail(f"could not determine physical memory: {error}")


def materialized_ld_ready(root: Path, policy: dict[str, Any]) -> bool:
    ld = policy["ld_reference_policy"]
    directory = contract.root_path(root, ld["materialized_directory"])
    expected_names = {
        *(ld["source_members"][component]["name"] for component in ("bed", "bim", "fam")),
        Path(ld["materialized_manifest"]).name,
        Path(ld["materialized_provenance"]).name,
    }
    if (
        not directory.is_dir()
        or {path.name for path in directory.iterdir()} != expected_names
    ):
        return False
    specification = importlib.util.spec_from_file_location(
        "track_b_pleiotropy_ld_for_gate", root / contract.LD_MATERIALIZER_SCRIPT,
    )
    if specification is None or specification.loader is None:
        return False
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    try:
        loaded_policy, loaded_ld, policy_path, contract_lock_sha256 = module.load_policy_context(
            root, contract.POLICY, test_fixture=False,
        )
        module.verify_reference(
            root, loaded_policy, loaded_ld, policy_path, contract_lock_sha256,
            test_fixture=False, quiet=True,
        )
    except (OSError, SystemExit, ValueError):
        return False
    return True


def readiness_rows(root: Path, policy: dict[str, Any]) -> list[dict[str, object]]:
    replication_ready = primary_replication_terminal(root)
    local_ready = local_analysis_terminal(root)
    replication_gate = "PASS_TERMINAL_PRIMARY_REPLICATION_FAMILY" if replication_ready else "BLOCKED_PRIMARY_REPLICATION_NOT_TERMINAL"
    local_gate = "PASS_TERMINAL" if local_ready else "BLOCKED_LOCAL_ANALYSIS_NOT_TERMINAL"
    upstream_blockers = []
    if not replication_ready:
        upstream_blockers.append("UPSTREAM_PRIMARY_REPLICATION")
    if not local_ready:
        upstream_blockers.append("UPSTREAM_LOCAL_ANALYSIS")

    placo = policy["placo_plus"]
    placo_source_ready = exact_file_ready(
        contract.root_path(root, placo["source_path"]), digest=placo["source_sha256"],
    )
    ld = policy["ld_reference_policy"]
    archive_ready = exact_file_ready(
        contract.root_path(root, ld["source_archive_path"]),
        size=ld["source_archive_expected_bytes"], digest=ld["source_archive_sha256"],
    )
    plink_path = contract.root_path(root, ld["plink_path"])
    plink_ready = exact_file_ready(plink_path, digest=ld["plink_sha256"])
    if plink_ready:
        version = subprocess.run(
            [str(plink_path), "--version"], check=False, capture_output=True, text=True,
        )
        plink_ready = version.returncode == 0 and (version.stdout + version.stderr).strip().startswith(ld["plink_version"])
    full_ld_ready = materialized_ld_ready(root, policy)
    active_ld_conflicts = [
        marker for marker in ld["conflicting_activity_markers"]
        if contract.root_path(root, marker).exists()
    ]

    conj = policy["conjfdr"]
    code_ready = clean_git_commit(contract.root_path(root, conj["code_path"]), conj["code_commit"])
    template_ready = exact_file_ready(
        contract.root_path(root, conj["variant_template_path"]),
        size=conj["variant_template_expected_bytes"], digest=conj["variant_template_sha256"],
    )
    reference_ready = exact_file_ready(
        contract.root_path(root, conj["reference_path"]), size=conj["reference_expected_bytes"],
    )
    runtime_provenance_ready = exact_file_ready(contract.root_path(root, conj["runtime_provenance_path"]))
    patch_ready = exact_file_ready(
        contract.root_path(root, conj["overlap_patch_path"]), digest=conj["overlap_patch_sha256"],
    )
    matlab_ready = bool(shutil.which("matlab"))
    memory_ready = physical_memory_bytes() >= int(conj["minimum_memory_bytes"])
    storage_ready = shutil.disk_usage(root).free >= int(conj["minimum_free_storage_bytes"])

    blocked_semantics = (
        "NOT_TESTED_METHOD_BLOCKED;NOT_ZERO_OR_NULL;NO_EMPTY_SCIENTIFIC_RESULT;"
        "NO_ONLY_LABEL_IF_OTHER_METHOD_BLOCKED"
    )
    placo_blockers = upstream_blockers + ["TRACK_B_FULL_P_RUNNER_AND_VALIDATOR_NOT_IMPLEMENTED"]
    if not placo_source_ready:
        placo_blockers.append("PLACO_SOURCE")
    if not archive_ready:
        placo_blockers.append("FULL_EUR_LD_SOURCE_ARCHIVE")
    if not plink_ready:
        placo_blockers.append("PLINK")
    if not full_ld_ready:
        placo_blockers.append("FULL_EUR_LD_REFERENCE_MATERIALIZATION")
    if active_ld_conflicts:
        placo_blockers.append("ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION")
    conj_blockers = upstream_blockers + ["TRACK_B_CONJFDR_RUNNER_AND_VALIDATOR_NOT_IMPLEMENTED"]
    if not code_ready or not patch_ready or not matlab_ready:
        conj_blockers.append("CONJFDR_SOFTWARE")
    if not template_ready or not reference_ready or not runtime_provenance_ready:
        conj_blockers.append("CONJFDR_REFERENCE_TEMPLATE_PROVENANCE")
    if not memory_ready or not storage_ready:
        conj_blockers.append("CONJFDR_COMPUTE")
    union_blockers = [
        "UPSTREAM_COMPLETED_METHOD_RESULTS", "TRACK_B_UNION_VALIDATOR_NOT_IMPLEMENTED",
    ]
    if not full_ld_ready:
        union_blockers.append("FULL_EUR_LD_REFERENCE_MATERIALIZATION")
    if active_ld_conflicts:
        union_blockers.append("ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION")

    return [
        {
            "component_id": "PLACO_PLUS", "pair_scope": "A;B;CONTROL",
            "dense_input_gate": "PASS_ALL_THREE_PAIRS", "replication_gate": replication_gate,
            "local_analysis_gate": local_gate,
            "implementation_gate": "BLOCKED_BY_IMPLEMENTATION_TRACK_B_FULL_P_RUNNER_AND_VALIDATOR",
            "software_gate": "READY_PINNED_PLACO_SOURCE" if placo_source_ready else "BLOCKED_BY_SOFTWARE_PLACO_SOURCE",
            "reference_gate": (
                "READY_FULL_EUR_HG19_LD" if full_ld_ready else
                "BLOCKED_BY_DATA_FULL_EUR_HG19_LD_NOT_MATERIALIZED"
            ),
            "compute_gate": "BLOCKED_BY_IMPLEMENTATION_NO_LOCKED_PLACO_RESOURCE_ENVELOPE",
            "scan_gate": "BLOCKED_UPSTREAM_AND_IMPLEMENTATION",
            "locus_publication_gate": "BLOCKED_FULL_LD_AND_VALIDATED_SCAN",
            "overall_status": "BLOCKED_BY_MULTIPLE_GATES", "blockers": ";".join(placo_blockers),
            "blocked_result_semantics": blocked_semantics, "claim_status": "NO_SCIENTIFIC_RESULT",
        },
        {
            "component_id": "CONJFDR", "pair_scope": "A;B;CONTROL",
            "dense_input_gate": "PASS_ALL_THREE_PAIRS", "replication_gate": replication_gate,
            "local_analysis_gate": local_gate,
            "implementation_gate": "BLOCKED_BY_IMPLEMENTATION_TRACK_B_CONJFDR_RUNNER_AND_VALIDATOR",
            "software_gate": (
                "READY_PINNED_CODE_PATCH_AND_MATLAB" if code_ready and patch_ready and matlab_ready
                else "BLOCKED_BY_SOFTWARE_CODE_PATCH_OR_MATLAB"
            ),
            "reference_gate": (
                "READY_SEALED_CONJFDR_REFERENCE_TEMPLATE" if template_ready and reference_ready and runtime_provenance_ready
                else "BLOCKED_BY_DATA_REFERENCE_TEMPLATE_OR_RUNTIME_PROVENANCE"
            ),
            "compute_gate": (
                "PASS_MEMORY_AND_STORAGE" if memory_ready and storage_ready
                else "BLOCKED_BY_COMPUTE_MEMORY_OR_STORAGE"
            ),
            "scan_gate": "BLOCKED_UPSTREAM_IMPLEMENTATION_RUNTIME_AND_COMPUTE",
            "locus_publication_gate": "BLOCKED_NO_VALIDATED_CONJFDR_SCAN",
            "overall_status": "BLOCKED_BY_MULTIPLE_GATES", "blockers": ";".join(conj_blockers),
            "blocked_result_semantics": blocked_semantics, "claim_status": "NO_SCIENTIFIC_RESULT",
        },
        {
            "component_id": "METHOD_UNION_08_09", "pair_scope": "A;B;CONTROL_SEPARATE",
            "dense_input_gate": "PASS_ALL_THREE_PAIRS", "replication_gate": replication_gate,
            "local_analysis_gate": local_gate,
            "implementation_gate": "BLOCKED_BY_IMPLEMENTATION_TRACK_B_UNION_VALIDATOR",
            "software_gate": "NOT_APPLICABLE_UNTIL_METHOD_RESULTS_EXIST",
            "reference_gate": (
                "READY_FULL_EUR_HG19_LD" if full_ld_ready else
                "BLOCKED_BY_DATA_FULL_EUR_HG19_LD_NOT_MATERIALIZED"
            ),
            "compute_gate": "NOT_ASSESSED_UNTIL_METHOD_RESULTS_EXIST",
            "scan_gate": "BLOCKED_UPSTREAM_NO_COMPLETED_METHOD_RESULTS",
            "locus_publication_gate": "BLOCKED_NO_VALIDATED_METHOD_UNION",
            "overall_status": "BLOCKED_BY_MULTIPLE_GATES", "blockers": ";".join(union_blockers),
            "blocked_result_semantics": blocked_semantics, "claim_status": "NO_SCIENTIFIC_RESULT",
        },
    ]


def build(root: Path) -> dict[Path, str]:
    policy = contract.validate_policy(root)
    upstream = contract.validate_upstream(root, policy)
    contract_lock = contract.verify_contract(root)
    scans = validate_live_dense_inputs(root, policy, upstream)
    input_text = tsv_text(INPUT_FIELDS, pair_gate_rows(root, policy, scans))
    readiness_text = tsv_text(READINESS_FIELDS, readiness_rows(root, policy))
    payload = {
        "schema_version": "sleep-atlas-track-b-pleiotropy-input-gate.1",
        "analysis_id": policy["analysis_id"],
        "selection_timing": policy["selection_timing"],
        "pleiotropy_results_accessed_before_input_gate_freeze": False,
        "policy_sha256": contract.sha256(root / contract.POLICY),
        "contract_lock_sha256": contract.sha256(root / contract.CONTRACT_LOCK),
        "contract_script_sha256": contract_lock["script_sha256"][str(contract.CONTRACT_SCRIPT)],
        "input_gate_script_sha256": contract_lock["script_sha256"][str(contract.INPUT_GATE_SCRIPT)],
        "ld_materializer_script_sha256": contract_lock["script_sha256"][str(contract.LD_MATERIALIZER_SCRIPT)],
        "pair_manifest_sha256": upstream["pair_manifest_sha256"],
        "pair_manifest_lock_sha256": upstream["pair_manifest_lock_sha256"],
        "dense_qc_sha256": upstream["dense_qc_sha256"],
        "dense_qc_lock_sha256": upstream["dense_qc_lock_sha256"],
        "input_gate_sha256": hashlib.sha256(input_text.encode("utf-8")).hexdigest(),
        "readiness_gate_sha256": hashlib.sha256(readiness_text.encode("utf-8")).hexdigest(),
        "pair_count": 3,
        "primary_pair_count": 2,
        "control_pair_count": 1,
        "unique_dense_trait_count": 5,
        "all_pair_dense_input_gates_pass": True,
        "full_gzip_integrity_scans_completed": 5,
        "live_dense_inputs": scans,
        "scientific_result_count": 0,
        "scientific_result_substitution_policy": "FORBIDDEN; readiness gates are not pleiotropy results",
    }
    return {
        INPUT_GATE: input_text,
        READINESS_GATE: readiness_text,
        INPUT_LOCK: json.dumps(payload, indent=2, sort_keys=True) + "\n",
    }


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def verify_outputs(root: Path, expected: dict[Path, str]) -> None:
    for relative, payload in expected.items():
        path = root / relative
        if not path.is_file() or path.read_text(encoding="utf-8") != payload:
            fail(f"Track B pleiotropy input/readiness gate missing or drifted: {relative}")
    print("TRACK_B_PLEIOTROPY_INPUT_GATE_VERIFIED pairs=3 traits=5 full_gzip_scans=5")


def write_outputs(root: Path, expected: dict[Path, str], policy: dict[str, Any]) -> None:
    evidence = contract.result_evidence(root, policy)
    if evidence:
        fail(
            "refusing to freeze Track B pleiotropy input gates after result access: "
            + ", ".join(str(path.relative_to(root)) for path in evidence[:3])
        )
    existing = [relative for relative in expected if (root / relative).exists()]
    if existing:
        if len(existing) != len(expected):
            fail("partial immutable Track B pleiotropy input-gate family already exists")
        verify_outputs(root, expected)
        return
    for relative, payload in expected.items():
        atomic_text(root / relative, payload)
    print("TRACK_B_PLEIOTROPY_INPUT_GATE_WRITTEN pairs=3 traits=5 methods_blocked=true science_results=0")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--build", action="store_true")
    action.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    expected = build(root)
    if args.verify:
        verify_outputs(root, expected)
    else:
        write_outputs(root, expected, contract.validate_policy(root))


if __name__ == "__main__":
    main()
