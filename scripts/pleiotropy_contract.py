#!/usr/bin/env python3
"""Seal and verify the PLACO+/conjunction-FDR production runtime."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path


REFERENCE_FIELDS = [
    "resource_id", "version_or_commit", "source_url", "expected_bytes",
    "verification", "local_path", "status",
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


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_policy(root: Path) -> tuple[Path, dict[str, object]]:
    path = root / "config/pleiotropy_analysis_policy.json"
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"pleiotropy policy is unreadable: {exc}")
    if (
        policy.get("analysis_id") != "atlas-v1.0-pleiotropic-loci"
        or policy.get("expected_traits") != 45
        or policy.get("expected_sleep_non_sleep_pairs") != 396
        or policy.get("placo_version") != "0.2.0"
    ):
        fail("pleiotropy policy differs from the frozen atlas scope")
    return path, policy


def load_registry(root: Path, policy: dict[str, object]) -> tuple[Path, list[dict[str, str]]]:
    path = root / str(policy["reference_source_registry"])
    if sha256(path) != policy["reference_source_registry_sha256"]:
        fail("pleiotropy source registry differs from policy")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = list(reader.fieldnames or []), list(reader)
    if fields != REFERENCE_FIELDS or [row["resource_id"] for row in rows] != [
        "placo_plus", "pleiofdr", "pleiofdr_reference", "pleiofdr_variant_template",
    ]:
        fail("pleiotropy source registry differs from the exact runtime family")
    expected_paths = {
        "placo_plus": policy["placo_source_path"],
        "pleiofdr": policy["pleiofdr_code_path"],
        "pleiofdr_reference": policy["pleiofdr_reference_path"],
        "pleiofdr_variant_template": policy["pleiofdr_variant_template_path"],
    }
    if any(row["local_path"] != expected_paths[row["resource_id"]] for row in rows):
        fail("pleiotropy registry paths differ from policy")
    return path, rows


def git_identity(path: Path) -> str:
    if not (path / ".git").is_dir():
        fail(f"pinned pleioFDR checkout is absent: {path}")
    commit = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    )
    status = subprocess.run(
        ["git", "-C", str(path), "status", "--porcelain"],
        capture_output=True, text=True, check=False,
    )
    if commit.returncode or status.returncode or status.stdout.strip():
        fail("pleioFDR checkout is unreadable or dirty")
    return commit.stdout.strip()


def runtime_payload(
    root: Path, policy_path: Path, policy: dict[str, object], reference_sha256: str,
) -> dict[str, object]:
    if not re.fullmatch(r"[0-9a-f]{64}", reference_sha256):
        fail("pleioFDR reference SHA-256 is invalid")
    registry_path, _ = load_registry(root, policy)
    placo = root / str(policy["placo_source_path"])
    code = root / str(policy["pleiofdr_code_path"])
    reference = root / str(policy["pleiofdr_reference_path"])
    template = root / str(policy["pleiofdr_variant_template_path"])
    patch = root / str(policy["pleiofdr_overlap_patch"])
    if (
        not placo.is_file() or placo.stat().st_size != 7515
        or sha256(placo) != policy["placo_source_sha256"]
    ):
        fail("PLACO+ source differs from its exact pin")
    commit = git_identity(code)
    if commit != policy["pleiofdr_commit"]:
        fail("pleioFDR checkout differs from its exact commit")
    if not reference.is_file() or reference.stat().st_size != policy["pleiofdr_reference_bytes"]:
        fail("pleioFDR LD reference is absent or has the wrong size")
    if (
        not template.is_file()
        or template.stat().st_size != policy["pleiofdr_variant_template_bytes"]
        or sha256(template) != policy["pleiofdr_variant_template_sha256"]
    ):
        fail("pleioFDR variant template differs from its exact pin")
    if sha256(patch) != policy["pleiofdr_overlap_patch_sha256"]:
        fail("pleioFDR overlap-correction patch differs from its exact pin")
    return {
        "schema_version": "sleep-atlas-pleiotropy-runtime.1",
        "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path),
        "source_registry": str(registry_path.relative_to(root)),
        "source_registry_sha256": sha256(registry_path),
        "placo_source": {
            "path": str(placo.relative_to(root)), "bytes": placo.stat().st_size,
            "sha256": sha256(placo),
        },
        "pleiofdr_code": {"path": str(code.relative_to(root)), "commit": commit},
        "pleiofdr_reference": {
            "path": str(reference.relative_to(root)), "bytes": reference.stat().st_size,
            "sha256": reference_sha256, "etag": policy["pleiofdr_reference_etag"],
        },
        "variant_template": {
            "path": str(template.relative_to(root)), "bytes": template.stat().st_size,
            "sha256": sha256(template),
        },
        "overlap_patch": {
            "path": str(patch.relative_to(root)), "bytes": patch.stat().st_size,
            "sha256": sha256(patch),
        },
        "contract_script_sha256": sha256(Path(__file__)),
        "verification": "Official HTTPS/Git acquisition; every file SHA-256 verified at seal time",
    }


def seal_runtime(root: Path) -> Path:
    policy_path, policy = load_policy(root)
    provenance_path = root / str(policy["runtime_provenance"])
    if provenance_path.exists():
        fail("immutable pleiotropy runtime provenance already exists")
    reference = root / str(policy["pleiofdr_reference_path"])
    payload = runtime_payload(root, policy_path, policy, sha256(reference))
    atomic_json(provenance_path, payload)
    print(f"PLEIOTROPY_RUNTIME_SEALED out={provenance_path.relative_to(root)}")
    return provenance_path


def validate_runtime(root: Path, *, rehash_reference: bool = False) -> dict[str, object]:
    policy_path, policy = load_policy(root)
    provenance_path = root / str(policy["runtime_provenance"])
    try:
        observed = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"pleiotropy runtime provenance is unreadable: {exc}")
    reference_record = observed.get("pleiofdr_reference")
    if not isinstance(reference_record, dict):
        fail("pleiotropy runtime provenance lacks a reference identity")
    expected = runtime_payload(
        root, policy_path, policy, str(reference_record.get("sha256", "")),
    )
    for key, value in expected.items():
        if observed.get(key) != value:
            fail(f"pleiotropy runtime provenance differs from policy: {key}")
    if rehash_reference:
        reference = root / str(policy["pleiofdr_reference_path"])
        if sha256(reference) != reference_record["sha256"]:
            fail("pleioFDR reference differs from its sealed SHA-256")
    return observed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--seal-runtime", action="store_true")
    group.add_argument("--verify-runtime", action="store_true")
    group.add_argument("--rehash-reference", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    if args.seal_runtime:
        seal_runtime(root)
    else:
        validate_runtime(root, rehash_reference=args.rehash_reference)
        print("PLEIOTROPY_RUNTIME_VALIDATED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
