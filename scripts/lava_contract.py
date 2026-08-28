#!/usr/bin/env python3
"""Seal and verify the portable LAVA reference, input, and result contracts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


SOURCE_FIELDS = ["archive_id", "chromosomes", "url", "archive_bytes", "archive_filename"]
DOWNLOAD_FIELDS = SOURCE_FIELDS + ["sha256"]
EXTRACTED_FIELDS = ["chromosome", "file_type", "path", "bytes", "sha256"]
INPUT_PROVENANCE_FIELDS = [
    "phenotype", "filename", "bytes", "sha256", "required_columns", "validation_status",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty checksum-bound artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, expected_fields: list[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty tabular artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = list(reader.fieldnames or []), list(reader)
    if expected_fields is not None and fields != expected_fields:
        fail(f"wrong schema for {path}: {fields}")
    return fields, rows


def safe_relative(value: str, label: str) -> Path:
    path = Path(value)
    if not value or path.is_absolute() or ".." in path.parts:
        fail(f"unsafe {label} path: {value}")
    return path


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def file_record(root: Path, path: Path) -> dict[str, object]:
    digest = sha256(path)
    return {
        "path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": digest,
    }


def load_policy(root: Path) -> tuple[Path, dict[str, object]]:
    path = root / "config/lava_analysis_policy.json"
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"LAVA policy is unreadable: {exc}")
    if (
        policy.get("analysis_id") != "atlas-v1.0-lava"
        or policy.get("lava_version") != "0.1.5"
        or policy.get("expected_traits") != 45
        or policy.get("expected_loci") != 2495
        or policy.get("planned_univariate_tests") != 112275
    ):
        fail("LAVA policy differs from the frozen atlas scope")
    return path, policy


def load_panel(root: Path, policy: dict[str, object]) -> tuple[Path, list[dict[str, str]]]:
    path = root / "config/analysis_panel.tsv"
    fields, rows = read_tsv(path)
    traits = [row.get("trait_id", "") for row in rows]
    if (
        "trait_id" not in fields or len(rows) != int(policy["expected_traits"])
        or len(set(traits)) != len(rows)
        or sum(row.get("domain") == "sleep" for row in rows) != 12
    ):
        fail("LAVA contract requires the exact locked 45-trait panel")
    return path, rows


def load_reference_rows(
    root: Path, policy: dict[str, object], *, require_live: bool,
) -> tuple[Path, Path, Path, list[dict[str, object]], list[dict[str, object]]]:
    source_path = root / str(policy["reference_source_registry"])
    if sha256(source_path) != policy["reference_source_registry_sha256"]:
        fail("LAVA reference source registry differs from policy")
    _, sources = read_tsv(source_path, SOURCE_FIELDS)
    if (
        len(sources) != int(policy["reference_expected_archives"])
        or len({row["archive_id"] for row in sources}) != len(sources)
        or sum(int(row["archive_bytes"]) for row in sources)
        != int(policy["reference_expected_archive_bytes"])
    ):
        fail("LAVA source registry differs from the exact seven-archive family")

    download_path = root / str(policy["reference_download_manifest"])
    extracted_path = root / str(policy["reference_extracted_manifest"])
    _, downloads = read_tsv(download_path, DOWNLOAD_FIELDS)
    _, extracted = read_tsv(extracted_path, EXTRACTED_FIELDS)
    if len(downloads) != len(sources):
        fail("LAVA download manifest does not cover the source registry")
    archive_records: list[dict[str, object]] = []
    payload_root = root / str(policy["reference_payload_root"])
    for source, observed in zip(sources, downloads):
        if any(observed[field] != source[field] for field in SOURCE_FIELDS):
            fail("LAVA download manifest differs from the ordered source registry")
        if not re.fullmatch(r"[0-9a-f]{64}", observed["sha256"]):
            fail(f"invalid LAVA archive SHA-256: {observed['archive_id']}")
        size = int(observed["archive_bytes"])
        path = payload_root / "archives" / safe_relative(observed["archive_filename"], "archive")
        if require_live and (not path.is_file() or path.stat().st_size != size):
            fail(f"LAVA archive is absent or has the wrong size: {path}")
        archive_records.append({
            "path": str(path.relative_to(root)), "bytes": size, "sha256": observed["sha256"],
        })

    expected_extracted = {
        (str(chromosome), suffix): f"lava-ukb-v1.1_chr{chromosome}.{suffix}"
        for chromosome in range(1, 23) for suffix in ("info", "bcor")
    }
    if len(extracted) != int(policy["reference_expected_extracted_files"]):
        fail("LAVA extracted manifest does not contain 44 chromosome payloads")
    extracted_records: list[dict[str, object]] = []
    observed_keys: set[tuple[str, str]] = set()
    for row in extracted:
        key = (row["chromosome"], row["file_type"])
        if key in observed_keys or expected_extracted.get(key) != row["path"]:
            fail("LAVA extracted manifest differs from the exact chromosome family")
        observed_keys.add(key)
        try:
            size = int(row["bytes"])
        except ValueError as exc:
            fail(f"invalid LAVA extracted byte count: {row['path']}")
        if size <= 0 or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]):
            fail(f"invalid LAVA extracted identity: {row['path']}")
        path = payload_root / safe_relative(row["path"], "extracted reference")
        if require_live and (not path.is_file() or path.stat().st_size != size):
            fail(f"LAVA extracted reference is absent or has the wrong size: {path}")
        extracted_records.append({
            "path": str(path.relative_to(root)), "bytes": size, "sha256": row["sha256"],
        })
    if observed_keys != set(expected_extracted):
        fail("LAVA extracted manifest omits a required chromosome payload")
    return source_path, download_path, extracted_path, archive_records, extracted_records


def reference_payload(root: Path, policy: dict[str, object], *, rehash: bool) -> dict[str, object]:
    source, download, extracted, archives, payloads = load_reference_rows(
        root, policy, require_live=True,
    )
    if rehash:
        for record in archives + payloads:
            if sha256(root / str(record["path"])) != record["sha256"]:
                fail(f"LAVA reference SHA-256 mismatch: {record['path']}")
    return {
        "schema_version": "sleep-atlas-lava-reference.1",
        "analysis_id": policy["analysis_id"],
        "reference": policy["reference"],
        "reference_population": policy["reference_population"],
        "source_registry": str(source.relative_to(root)),
        "source_registry_sha256": sha256(source),
        "download_manifest": str(download.relative_to(root)),
        "download_manifest_sha256": sha256(download),
        "extracted_manifest": str(extracted.relative_to(root)),
        "extracted_manifest_sha256": sha256(extracted),
        "archive_count": len(archives),
        "archive_bytes": sum(int(row["bytes"]) for row in archives),
        "extracted_file_count": len(payloads),
        "extracted_bytes": sum(int(row["bytes"]) for row in payloads),
        "archives": archives,
        "extracted_files": payloads,
        "verification": "SHA-256 verified after official HTTPS acquisition",
    }


def seal_reference(root: Path, policy: dict[str, object]) -> Path:
    path = root / str(policy["reference_provenance"])
    if path.exists():
        fail("immutable LAVA reference provenance already exists")
    payload = reference_payload(root, policy, rehash=True)
    payload["contract_script_sha256"] = sha256(Path(__file__))
    atomic_json(path, payload)
    print(
        f"LAVA_REFERENCE_SEALED archives={payload['archive_count']} "
        f"files={payload['extracted_file_count']} out={path.relative_to(root)}"
    )
    return path


def validate_reference(root: Path, policy: dict[str, object], *, rehash: bool = False) -> dict[str, object]:
    path = root / str(policy["reference_provenance"])
    try:
        observed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"LAVA reference provenance is unreadable: {exc}")
    expected = reference_payload(root, policy, rehash=rehash)
    for key, value in expected.items():
        if observed.get(key) != value:
            fail(f"LAVA reference provenance differs from the sealed contract: {key}")
    return observed


def input_lock_payload(
    root: Path, policy_path: Path, policy: dict[str, object], panel_path: Path,
    panel: list[dict[str, str]],
) -> dict[str, object]:
    artifacts = []
    for relative in policy["input_artifacts"]:
        path = root / str(relative)
        artifacts.append(file_record(root, path))
    provenance_path = root / "results/tables/lava_input_provenance.tsv"
    _, provenance = read_tsv(provenance_path, INPUT_PROVENANCE_FIELDS)
    traits = [row["trait_id"] for row in panel]
    if [row["phenotype"] for row in provenance] != traits:
        fail("LAVA input provenance differs from the ordered 45-trait panel")
    sources = []
    for row in provenance:
        relative = safe_relative(row["filename"], "LAVA summary-statistics")
        path = root / relative
        try:
            size = int(row["bytes"])
        except ValueError:
            fail(f"invalid LAVA input byte count: {row['phenotype']}")
        if (
            row["validation_status"] != "VALIDATED" or size <= 0
            or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"])
            or not path.is_file() or path.stat().st_size != size
        ):
            fail(f"LAVA source input differs from provenance: {row['phenotype']}")
        sources.append({
            "phenotype": row["phenotype"], "path": str(relative),
            "bytes": size, "sha256": row["sha256"],
        })
    scripts = [
        "scripts/31_prepare_lava.py", "scripts/33_run_lava.R",
        "scripts/34_validate_lava.py", "scripts/lava_contract.py",
    ]
    return {
        "schema_version": "sleep-atlas-lava-inputs.1",
        "analysis_id": policy["analysis_id"],
        "panel_sha256": sha256(panel_path),
        "policy_sha256": sha256(policy_path),
        "trait_ids_in_locked_order": traits,
        "artifacts": artifacts,
        "source_files": sources,
        "script_sha256": {relative: sha256(root / relative) for relative in scripts},
    }


def write_input_lock(root: Path) -> Path:
    policy_path, policy = load_policy(root)
    panel_path, panel = load_panel(root, policy)
    path = root / str(policy["input_lock"])
    expected = input_lock_payload(root, policy_path, policy, panel_path, panel)
    if path.exists():
        try:
            observed = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            fail(f"LAVA input lock is unreadable: {exc}")
        if observed != expected:
            fail("immutable LAVA input lock differs from current artifacts")
        return path
    atomic_json(path, expected)
    return path


def validate_inputs(root: Path, policy_path: Path, policy: dict[str, object]) -> dict[str, object]:
    panel_path, panel = load_panel(root, policy)
    path = root / str(policy["input_lock"])
    try:
        observed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"LAVA input lock is unreadable: {exc}")
    expected = input_lock_payload(root, policy_path, policy, panel_path, panel)
    if observed != expected:
        fail("LAVA inputs differ from their immutable lock")
    return observed


def run_fingerprint(root: Path) -> str:
    policy_path, policy = load_policy(root)
    validate_reference(root, policy)
    validate_inputs(root, policy_path, policy)
    payload = {
        "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path),
        "input_lock_sha256": sha256(root / str(policy["input_lock"])),
        "reference_provenance_sha256": sha256(root / str(policy["reference_provenance"])),
        "runtime_script_sha256": sha256(root / "scripts/33_run_lava.R"),
        "validator_script_sha256": sha256(root / "scripts/34_validate_lava.py"),
        "contract_script_sha256": sha256(Path(__file__)),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def result_records(root: Path, policy: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    artifacts = []
    for relative in policy["result_artifacts"]:
        path = root / str(relative)
        artifacts.append(file_record(root, path))
    checkpoint_root = root / str(policy["checkpoint_directory"])
    checkpoints = []
    for index in range(1, int(policy["expected_loci"]) + 1):
        path = checkpoint_root / f"locus_{index:04d}.rds"
        checkpoints.append(file_record(root, path))
    return artifacts, checkpoints


def result_payload(root: Path, policy: dict[str, object]) -> dict[str, object]:
    artifacts, checkpoints = result_records(root, policy)
    return {
        "schema_version": "sleep-atlas-lava-results.1",
        "analysis_id": policy["analysis_id"],
        "lava_version": policy["lava_version"],
        "run_fingerprint": run_fingerprint(root),
        "input_lock_sha256": sha256(root / str(policy["input_lock"])),
        "reference_provenance_sha256": sha256(root / str(policy["reference_provenance"])),
        "checkpoint_count": len(checkpoints),
        "checkpoints": checkpoints,
        "result_count": len(artifacts),
        "results": artifacts,
    }


def seal_results(root: Path, policy: dict[str, object]) -> Path:
    path = root / str(policy["result_provenance"])
    expected = result_payload(root, policy)
    if path.exists():
        validate_results(root, policy)
        return path
    atomic_json(path, expected)
    print(f"LAVA_RESULTS_SEALED checkpoints={expected['checkpoint_count']} out={path.relative_to(root)}")
    return path


def validate_results(root: Path, policy: dict[str, object]) -> dict[str, object]:
    path = root / str(policy["result_provenance"])
    try:
        observed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"LAVA result provenance is unreadable: {exc}")
    expected = result_payload(root, policy)
    if observed != expected:
        fail("LAVA results differ from their immutable provenance")
    return observed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--seal-reference", action="store_true")
    group.add_argument("--verify-reference", action="store_true")
    group.add_argument("--run-fingerprint", action="store_true")
    group.add_argument("--seal-results", action="store_true")
    group.add_argument("--verify-results", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    _, policy = load_policy(root)
    if args.seal_reference:
        seal_reference(root, policy)
    elif args.verify_reference:
        validate_reference(root, policy, rehash=True)
        print("LAVA_REFERENCE_VALIDATED")
    elif args.run_fingerprint:
        print(run_fingerprint(root))
    elif args.seal_results:
        seal_results(root, policy)
    else:
        validate_results(root, policy)
        print("LAVA_RESULTS_VALIDATED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
