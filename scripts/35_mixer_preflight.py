#!/usr/bin/env python3
"""Audit MiXeR inputs, reference files, storage, and supported hardware."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_reference_manifest(
    root: Path, policy: dict[str, object],
) -> tuple[Path, list[dict[str, object]]]:
    manifest_path = root / str(policy["reference_file_manifest"])
    if (
        not manifest_path.is_file()
        or sha256(manifest_path) != policy["reference_file_manifest_sha256"]
    ):
        raise SystemExit("ERROR: MiXeR reference file manifest is absent or differs from policy")
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["path", "bytes", "sha256", "git_blob_sha1"]:
            raise SystemExit("ERROR: MiXeR reference file manifest has the wrong schema")
        source_rows = list(reader)
    rows: list[dict[str, object]] = []
    paths: set[str] = set()
    for source in source_rows:
        relative = Path(source["path"])
        if relative.is_absolute() or ".." in relative.parts or source["path"] in paths:
            raise SystemExit("ERROR: MiXeR reference file manifest has an unsafe or duplicate path")
        try:
            size = int(source["bytes"])
        except ValueError as exc:
            raise SystemExit("ERROR: MiXeR reference file manifest has a nonnumeric byte count") from exc
        if (
            size <= 0 or not re.fullmatch(r"[0-9a-f]{64}", source["sha256"])
            or not re.fullmatch(r"[0-9a-f]{40}", source["git_blob_sha1"])
        ):
            raise SystemExit("ERROR: MiXeR reference file manifest has an invalid object identity")
        paths.add(source["path"])
        rows.append({**source, "bytes": size})
    reference_root = Path(str(policy["reference_root"]))
    payload_root = Path(str(policy["reference_payload_root"]))
    try:
        family_root = reference_root.relative_to(payload_root)
    except ValueError as exc:
        raise SystemExit("ERROR: MiXeR reference roots are not nested as required") from exc
    expected_paths = {
        str(family_root / str(policy["reference_bim_pattern"]).replace("@", str(chromosome)))
        for chromosome in range(1, 23)
    } | {
        str(family_root / str(policy["reference_ld_pattern"]).replace("@", str(chromosome)))
        for chromosome in range(1, 23)
    } | {
        str(family_root / str(policy["reference_extract_pattern"]).replace("@", str(replicate)))
        for replicate in range(1, int(policy["fit_replicates"]) + 1)
    }
    if (
        paths != expected_paths
        or len(rows) != int(policy["reference_expected_files"])
        or sum(int(row["bytes"]) for row in rows) != int(policy["reference_expected_bytes"])
    ):
        raise SystemExit("ERROR: MiXeR reference file manifest differs from the exact runtime family")
    return manifest_path, rows


def seal_reference(
    root: Path, policy: dict[str, object], manifest_path: Path,
    records: list[dict[str, object]],
) -> Path:
    provenance_path = root / str(policy["reference_provenance_path"])
    if provenance_path.exists():
        raise SystemExit("ERROR: MiXeR reference provenance is immutable and already exists")
    payload_root = root / str(policy["reference_payload_root"])
    files = []
    for record in records:
        path = payload_root / str(record["path"])
        if not path.is_file() or path.stat().st_size != record["bytes"]:
            raise SystemExit(f"ERROR: MiXeR reference file is absent or has the wrong size: {path}")
        observed_hash = sha256(path)
        if observed_hash != record["sha256"]:
            raise SystemExit(f"ERROR: MiXeR reference file differs from its exact LFS identity: {path}")
        files.append({
            "path": str(path.relative_to(root)), "bytes": path.stat().st_size,
            "sha256": observed_hash, "git_blob_sha1": record["git_blob_sha1"],
        })
    provenance = {
        "schema_version": "sleep-atlas-mixer-reference.1",
        "analysis_id": policy["analysis_id"],
        "source_repository": policy["mixer_reference_repository"],
        "source_commit": policy["mixer_reference_commit_at_lock"],
        "reference_file_manifest": str(manifest_path.relative_to(root)),
        "reference_file_manifest_sha256": sha256(manifest_path),
        "file_count": len(files),
        "total_bytes": sum(int(row["bytes"]) for row in files),
        "files": files,
        "verification": "Every payload was SHA-256 verified against its Git LFS object identity",
        "script_sha256": sha256(Path(__file__)),
    }
    atomic_json(provenance_path, provenance)
    print(
        f"MIXER_REFERENCE_SEALED files={len(files)} bytes={provenance['total_bytes']} "
        f"out={provenance_path.relative_to(root)}"
    )
    return provenance_path


def reference_check(
    root: Path, policy: dict[str, object], manifest_path: Path,
    records: list[dict[str, object]],
) -> dict[str, object]:
    payload_root = root / str(policy["reference_payload_root"])
    provenance_path = root / str(policy["reference_provenance_path"])
    present = sum(
        (payload_root / str(record["path"])).is_file()
        and (payload_root / str(record["path"])).stat().st_size == record["bytes"]
        for record in records
    )
    result: dict[str, object] = {
        "pass": False,
        "present_files": present,
        "expected_files": len(records),
        "expected_bytes": sum(int(row["bytes"]) for row in records),
        "manifest_sha256": sha256(manifest_path),
        "provenance_path": str(provenance_path.relative_to(root)),
        "verification_mode": "sealed_SHA256_plus_live_size",
    }
    if not provenance_path.is_file():
        result["blocker"] = "exact reference provenance is absent; run --seal-reference after acquisition"
        return result
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        result["blocker"] = f"reference provenance is unreadable: {exc}"
        return result
    expected_files = [
        {
            "path": str(
                Path(str(policy["reference_payload_root"])) / str(record["path"])
            ),
            "bytes": record["bytes"], "sha256": record["sha256"],
            "git_blob_sha1": record["git_blob_sha1"],
        }
        for record in records
    ]
    if (
        provenance.get("schema_version") != "sleep-atlas-mixer-reference.1"
        or provenance.get("analysis_id") != policy["analysis_id"]
        or provenance.get("source_repository") != policy["mixer_reference_repository"]
        or provenance.get("source_commit") != policy["mixer_reference_commit_at_lock"]
        or provenance.get("reference_file_manifest_sha256") != sha256(manifest_path)
        or provenance.get("file_count") != len(records)
        or provenance.get("total_bytes") != sum(int(row["bytes"]) for row in records)
        or provenance.get("files") != expected_files
    ):
        result["blocker"] = "reference provenance differs from the exact 64-file contract"
        return result
    if present != len(records):
        result["blocker"] = "one or more sealed reference files are absent or have the wrong size"
        return result
    result["pass"] = True
    result["blocker"] = ""
    result["provenance_sha256"] = sha256(provenance_path)
    return result


def mixer_input_check(
    root: Path, policy: dict[str, object], panel: list[dict[str, str]],
    panel_path: Path, policy_path: Path,
) -> dict[str, object]:
    manifest_path = root / str(policy["input_manifest_path"])
    lock_path = root / str(policy["input_manifest_lock_path"])
    result: dict[str, object] = {
        "pass": False,
        "ready_traits": 0,
        "expected_traits": int(policy["expected_traits"]),
        "manifest_path": str(manifest_path.relative_to(root)),
        "lock_path": str(lock_path.relative_to(root)),
        "verification_mode": "sealed_SHA256_plus_live_size",
    }
    if not manifest_path.is_file() or not lock_path.is_file():
        result["blocker"] = "checksum-locked converted MiXeR input family is absent"
        return result
    try:
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = list(reader.fieldnames or []), list(reader)
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        result["blocker"] = f"MiXeR input manifest or lock is unreadable: {exc}"
        return result
    trait_ids = [row["trait_id"] for row in panel]
    if (
        fields != policy["input_manifest_fields"]
        or [row.get("trait_id") for row in rows] != trait_ids
        or lock.get("schema_version") != "sleep-atlas-mixer-inputs.1"
        or lock.get("analysis_id") != policy["analysis_id"]
        or lock.get("panel_sha256") != sha256(panel_path)
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("trait_ids_in_locked_order") != trait_ids
        or lock.get("input_manifest") != str(manifest_path.relative_to(root))
        or lock.get("input_manifest_sha256") != sha256(manifest_path)
        or lock.get("input_scope") != "FULL_POST_QC_AUTOSOMAL"
    ):
        result["blocker"] = "MiXeR input manifest or lock differs from the exact panel/policy contract"
        return result
    expected_hashes: dict[str, str] = {}
    ready = 0
    for row in rows:
        trait = row["trait_id"]
        expected_filename = f"data/mixer/{trait}.sumstats.gz"
        try:
            source_bytes = int(row["source_bytes"])
            output_bytes = int(row["output_bytes"])
            row_count = int(row["rows"])
        except ValueError:
            result["blocker"] = f"MiXeR input manifest has nonnumeric counts for {trait}"
            return result
        if (
            row["filename"] != expected_filename
            or source_bytes <= 0 or output_bytes <= 0 or row_count <= 0
            or not re.fullmatch(r"[0-9a-f]{64}", row["source_sha256"])
            or not re.fullmatch(r"[0-9a-f]{64}", row["output_sha256"])
            or row["input_scope"] != "FULL_POST_QC_AUTOSOMAL"
            or row["status"] != "VALIDATED"
        ):
            result["blocker"] = f"MiXeR input manifest row differs from policy for {trait}"
            return result
        path = root / expected_filename
        if not path.is_file() or path.stat().st_size != output_bytes:
            result["blocker"] = f"converted MiXeR input is absent or has the wrong size: {trait}"
            return result
        expected_hashes[trait] = row["output_sha256"]
        ready += 1
    if lock.get("input_file_sha256") != expected_hashes:
        result["blocker"] = "MiXeR input lock differs from the converted file identities"
        return result
    result.update({
        "pass": ready == int(policy["expected_traits"]),
        "ready_traits": ready,
        "blocker": "" if ready == int(policy["expected_traits"]) else "converted input family is incomplete",
        "manifest_sha256": sha256(manifest_path),
        "lock_sha256": sha256(lock_path),
    })
    return result


def qc_value(path: Path, key: str) -> str:
    if not path.is_file():
        return ""
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0] == key:
            return fields[1]
    return ""


def choose_harmonized(root: Path, trait: str) -> tuple[Path, Path, str]:
    candidates = [
        (
            root / "data/harmonized_mixer_full" / f"{trait}.harmonized.tsv.gz",
            root / "data/harmonized_mixer_full" / f"{trait}.qc.txt",
        ),
        (
            root / "data/harmonized" / f"{trait}.harmonized.tsv.gz",
            root / "data/harmonized" / f"{trait}.qc.txt",
        ),
    ]
    for data, qc in candidates:
        if data.is_file() and data.stat().st_size and qc.is_file():
            prefilter = qc_value(qc, "prefilter_strategy")
            variant_map = qc_value(qc, "variant_map_strategy")
            if variant_map not in {"not supplied"}:
                prefilter = f"HAPMAP3_VARIANT_MAP_{variant_map or 'UNKNOWN'}"
            return data, qc, prefilter
    return candidates[0][0], candidates[0][1], ""


def physical_core_count() -> tuple[int, str]:
    """Return physical cores when the host exposes them, else a labelled fallback."""
    if platform.system() == "Darwin" and shutil.which("sysctl"):
        result = subprocess.run(
            ["sysctl", "-n", "hw.physicalcpu"], capture_output=True, text=True
        )
        if result.returncode == 0 and result.stdout.strip().isdigit():
            return int(result.stdout.strip()), "sysctl hw.physicalcpu"
    cpuinfo = Path("/proc/cpuinfo")
    if platform.system() == "Linux" and cpuinfo.is_file():
        pairs = set()
        physical_id = core_id = None
        for line in cpuinfo.read_text(encoding="utf-8", errors="ignore").splitlines() + [""]:
            if not line:
                if physical_id is not None and core_id is not None:
                    pairs.add((physical_id, core_id))
                physical_id = core_id = None
            elif line.startswith("physical id"):
                physical_id = line.split(":", 1)[1].strip()
            elif line.startswith("core id"):
                core_id = line.split(":", 1)[1].strip()
        if pairs:
            return len(pairs), "/proc/cpuinfo physical/core IDs"
    return os.cpu_count() or 0, "logical-core fallback"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--seal-reference", action="store_true")
    parser.add_argument("--json-out", default="results/tables/mixer_preflight.json")
    parser.add_argument("--trait-out", default="results/tables/mixer_input_readiness.tsv")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / "config/mixer_analysis_policy.json"
    panel_path = root / "config/analysis_panel.tsv"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    reference_manifest_path, reference_records = load_reference_manifest(root, policy)
    if args.seal_reference:
        seal_reference(root, policy, reference_manifest_path, reference_records)
    panel = read_tsv(panel_path)
    if len(panel) != 45 or len({row["trait_id"] for row in panel}) != 45:
        raise SystemExit("ERROR: expected the exact locked 45-trait panel")

    planned_prefilters = {
        row["trait_id"] for row in read_tsv(root / "config/hm3_prefilter_plans.tsv")
    }
    trait_rows = []
    for row in panel:
        trait = row["trait_id"]
        data, qc, prefilter = choose_harmonized(root, trait)
        if not prefilter and trait not in planned_prefilters:
            prefilter = "not supplied"
        raw = root / "data/raw" / row["raw_file"]
        full = data.is_file() and qc.is_file() and prefilter == "not supplied"
        if full:
            status = "READY_FULL_SUMSTATS"
            blocker = ""
        elif data.is_file() and prefilter.startswith("HAPMAP3_"):
            status = "BLOCKED_HAPMAP3_CONSTRAINED"
            blocker = "regenerate full post-QC summary statistics without a HapMap3 allowlist or HapMap3-only variant map"
        else:
            status = "BLOCKED_FULL_HARMONIZED_MISSING"
            blocker = "materialize full post-QC harmonized summary statistics"
        trait_rows.append({
            "trait_id": trait,
            "harmonized_file": str(data.relative_to(root)),
            "harmonized_bytes": data.stat().st_size if data.is_file() else 0,
            "qc_file": str(qc.relative_to(root)),
            "prefilter_strategy": prefilter or "UNKNOWN",
            "raw_file": str(raw.relative_to(root)),
            "raw_available": str(raw.is_file()).upper(),
            "input_status": status,
            "blocker": blocker,
        })

    machine = platform.machine().lower()
    memory_bytes = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    logical_cores = os.cpu_count() or 0
    physical_cores, core_source = physical_core_count()
    free_bytes = shutil.disk_usage(root).free
    engine = "docker" if shutil.which("docker") else ""
    checks = {
        "architecture": {
            "pass": machine in policy["supported_architectures"],
            "observed": machine,
            "required": policy["supported_architectures"],
        },
        "memory": {
            "pass": memory_bytes >= policy["minimum_memory_bytes"],
            "observed_bytes": memory_bytes,
            "minimum_bytes": policy["minimum_memory_bytes"],
        },
        "cores": {
            "pass": physical_cores >= policy["recommended_physical_cores"],
            "observed_physical_cores": physical_cores,
            "observed_logical_cores": logical_cores,
            "measurement": core_source,
            "recommended_physical_cores": policy["recommended_physical_cores"],
        },
        "storage": {
            "pass": free_bytes >= policy["minimum_free_storage_bytes"],
            "observed_free_bytes": free_bytes,
            "minimum_free_bytes": policy["minimum_free_storage_bytes"],
        },
        "container_engine": {
            "pass": bool(engine),
            "observed": engine or "NONE",
        },
        "reference": reference_check(
            root, policy, reference_manifest_path, reference_records,
        ),
        "full_trait_inputs": {
            "pass": all(row["input_status"] == "READY_FULL_SUMSTATS" for row in trait_rows),
            "ready_traits": sum(row["input_status"] == "READY_FULL_SUMSTATS" for row in trait_rows),
            "expected_traits": 45,
            "blocked_traits": [row["trait_id"] for row in trait_rows if row["input_status"] != "READY_FULL_SUMSTATS"],
        },
        "converted_mixer_inputs": mixer_input_check(
            root, policy, panel, panel_path, policy_path,
        ),
    }
    overall = all(check["pass"] for check in checks.values())
    report = {
        "analysis_id": policy["analysis_id"],
        "mixer_release": policy["mixer_release"],
        "container_amd64_manifest_digest": policy["container_amd64_manifest_digest"],
        "ready": overall,
        "checks": checks,
    }

    trait_out = root / args.trait_out
    trait_out.parent.mkdir(parents=True, exist_ok=True)
    with trait_out.open("w", encoding="utf-8", newline="") as handle:
        fields = list(trait_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(trait_rows)
    json_out = root / args.json_out
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    for name, check in checks.items():
        print(f"{'PASS' if check['pass'] else 'BLOCKED':7} {name}: {check}")
    print(f"MiXeR preflight: {'READY' if overall else 'BLOCKED'}")
    return 0 if overall or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
