#!/usr/bin/env python3
"""Fail-closed readiness audit for the pinned LAVA and HDL-L layer."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_METHOD_FIELDS = {
    "method_id", "analysis_role", "software", "package_version", "source_repository",
    "source_commit", "source_commit_date", "entrypoint", "paper_title", "paper_DOI",
    "paper_PMID", "reference_id", "reference_ancestry", "reference_description", "use_policy",
}
REQUIRED_REFERENCE_FIELDS = {
    "reference_id", "component_id", "chromosomes", "official_source_url",
    "expected_archive_name", "published_unzipped_size_gib", "published_checksum",
    "local_reference_root", "lock_file", "acquisition_status", "verification_requirement",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = required - fields
        if missing:
            raise SystemExit(f"ERROR: {path} lacks fields: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"ERROR: {path} is empty")
    return rows


def inspect_r_runtime(rscript: Path) -> tuple[dict[str, dict[str, str]], str]:
    if not rscript.is_file():
        return {}, f"Rscript missing: {rscript}"
    expression = (
        "for (x in list(c('LAVA','run.bivar'),c('HDL','HDL.L'))) {"
        "p<-x[[1]]; f<-x[[2]]; ok<-requireNamespace(p,quietly=TRUE);"
        "v<-if(ok) as.character(packageVersion(p)) else 'MISSING';"
        "e<-if(ok) exists(f,where=asNamespace(p),inherits=FALSE) else FALSE;"
        "cat(p,v,f,e,sep='\\t');cat('\\n')}"
    )
    result = subprocess.run(
        [str(rscript), "-e", expression], check=False, text=True, capture_output=True,
    )
    if result.returncode != 0:
        return {}, (result.stdout + result.stderr).strip()
    observed: dict[str, dict[str, str]] = {}
    for line in result.stdout.splitlines():
        values = line.split("\t")
        if len(values) == 4:
            package, version, entrypoint, available = values
            observed[package] = {
                "version": version, "entrypoint": entrypoint, "entrypoint_available": available,
            }
    return observed, result.stderr.strip()


def verify_reference_lock(reference_id: str, lock_path: Path) -> tuple[str, str, int]:
    if not lock_path.is_file():
        return "MISSING_UNVERIFIED", f"reference lock absent: {lock_path}", 0
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return "INVALID_LOCK", f"reference lock cannot be read: {exc}", 0
    if lock.get("reference_id") != reference_id or lock.get("status") != "VERIFIED":
        return "INVALID_LOCK", "reference_id/status is not the required VERIFIED identity", 0
    files = lock.get("files")
    if not isinstance(files, list) or not files:
        return "INVALID_LOCK", "verified reference lock has no files", 0
    total = 0
    for item in files:
        path = Path(str(item.get("path", "")))
        expected = str(item.get("sha256", ""))
        if not path.is_file() or len(expected) != 64 or sha256(path) != expected:
            return "CHECKSUM_FAILURE", f"missing or checksum-mismatched reference file: {path}", total
        total += path.stat().st_size
    return "VERIFIED", f"{len(files)} locked files passed SHA256", total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--methods", type=Path, default=Path("discovery_extension/config/local_architecture_sources.tsv"))
    parser.add_argument("--references", type=Path, default=Path("discovery_extension/config/local_reference_sources.tsv"))
    parser.add_argument("--rscript", type=Path, default=Path(".r-env/bin/Rscript"))
    parser.add_argument("--discovery-rg", type=Path, default=Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"))
    parser.add_argument("--dense-manifest", type=Path, default=Path("discovery_extension/provenance/dense_local_inputs.lock.json"))
    parser.add_argument("--storage-path", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/local/local_architecture_readiness.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/local_architecture_dependencies.json"))
    args = parser.parse_args()

    methods = read_tsv(args.methods, REQUIRED_METHOD_FIELDS)
    references = read_tsv(args.references, REQUIRED_REFERENCE_FIELDS)
    if len({row["method_id"] for row in methods}) != len(methods):
        raise SystemExit("ERROR: duplicate method_id")
    if {row["method_id"] for row in methods} != {"LAVA_PRIMARY", "HDL_L_ROBUSTNESS"}:
        raise SystemExit("ERROR: local method contract must contain exact LAVA primary and HDL-L robustness rows")
    ref_ids = {row["reference_id"] for row in references}
    if any(row["reference_id"] not in ref_ids for row in methods):
        raise SystemExit("ERROR: method references an unregistered LD source")
    if len([row for row in references if row["reference_id"] == "LAVA_UKB_V1.1"]) != 7:
        raise SystemExit("ERROR: LAVA UKB v1.1 source manifest must contain seven official download components")

    observed, runtime_stderr = inspect_r_runtime(args.rscript)
    usage = shutil.disk_usage(args.storage_path.resolve())
    available_gib = usage.free / 1024**3
    lava_published_unzipped_gib = 15.0
    discovery_ready = args.discovery_rg.is_file()
    dense_ready = False
    dense_reason = f"dense input lock absent: {args.dense_manifest}"
    if args.dense_manifest.is_file():
        try:
            dense_lock = json.loads(args.dense_manifest.read_text(encoding="utf-8"))
            dense_ready = dense_lock.get("status") == "VERIFIED" and dense_lock.get("full_resolution") is True
            dense_reason = "full-resolution dense input lock verified" if dense_ready else "dense input lock is not VERIFIED/full_resolution"
        except (json.JSONDecodeError, OSError) as exc:
            dense_reason = f"dense input lock invalid: {exc}"

    lock_paths: dict[str, Path] = {}
    for reference_id in ref_ids:
        paths = {Path(row["lock_file"]) for row in references if row["reference_id"] == reference_id}
        if len(paths) != 1:
            raise SystemExit(f"ERROR: inconsistent lock_file for {reference_id}")
        lock_paths[reference_id] = next(iter(paths))

    rows: list[dict[str, object]] = []
    all_code_ready = True
    all_refs_ready = True
    for method in methods:
        package = observed.get(method["software"], {})
        version_ok = package.get("version") == method["package_version"]
        function_ok = package.get("entrypoint") == method["entrypoint"] and package.get("entrypoint_available") == "TRUE"
        code_ready = version_ok and function_ok
        all_code_ready &= code_ready
        reference_status, reference_detail, reference_bytes = verify_reference_lock(
            method["reference_id"], lock_paths[method["reference_id"]],
        )
        all_refs_ready &= reference_status == "VERIFIED"
        if method["reference_id"] == "LAVA_UKB_V1.1" and reference_status != "VERIFIED" and available_gib < lava_published_unzipped_gib:
            reference_status = "MISSING_AND_INSUFFICIENT_SPACE_FOR_PUBLISHED_UNZIPPED_SIZE"
            reference_detail += f"; published unzipped size=15 GiB, available={available_gib:.3f} GiB"
        if not code_ready:
            readiness = "BLOCKED_CODE"
        elif reference_status != "VERIFIED":
            readiness = "CODE_READY_REFERENCE_BLOCKED"
        elif not discovery_ready or not dense_ready:
            readiness = "CODE_AND_REFERENCE_READY_DENSE_SUMSTATS_BLOCKED"
        else:
            readiness = "READY"
        rows.append({
            "method_id": method["method_id"], "analysis_role": method["analysis_role"],
            "software": method["software"], "expected_version": method["package_version"],
            "observed_version": package.get("version", "MISSING"),
            "entrypoint": method["entrypoint"], "entrypoint_available": package.get("entrypoint_available", "FALSE"),
            "source_commit": method["source_commit"], "code_status": "PASS" if code_ready else "FAIL",
            "reference_id": method["reference_id"], "reference_status": reference_status,
            "reference_detail": reference_detail, "verified_reference_bytes": reference_bytes,
            "discovery_rg_status": "PRESENT" if discovery_ready else "BLOCKED_UPSTREAM_DISCOVERY_RG_ABSENT",
            "dense_summary_statistics_status": "VERIFIED" if dense_ready else "BLOCKED_UPSTREAM_FULL_RESOLUTION_ACQUISITION",
            "dense_summary_statistics_detail": dense_reason,
            "available_free_gib": f"{available_gib:.3f}", "readiness": readiness,
        })

    fields = list(rows[0])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if not all_code_ready:
        overall = "BLOCKED_CODE"
    elif not all_refs_ready or not discovery_ready or not dense_ready:
        overall = "CODE_READY_REFERENCE_AND_DENSE_SUMSTATS_BLOCKED"
    else:
        overall = "READY"
    provenance = {
        "schema_version": "1.0.0",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall_status": overall,
        "download_started": False,
        "methods_sha256": sha256(args.methods),
        "references_sha256": sha256(args.references),
        "rscript": str(args.rscript),
        "runtime_observations": observed,
        "runtime_stderr": runtime_stderr,
        "storage_path": str(args.storage_path.resolve()),
        "available_free_bytes": usage.free,
        "available_free_gib": round(available_gib, 3),
        "lava_reference_published_unzipped_gib": lava_published_unzipped_gib,
        "discovery_rg_present": discovery_ready,
        "dense_input_lock_present": args.dense_manifest.is_file(),
        "dense_input_verified": dense_ready,
        "readiness_output": str(args.out),
        "readiness_output_sha256": sha256(args.out),
        "warning": "Code readiness is not analysis completion; no local genetic result exists without verified references and dense inputs.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"LOCAL_ARCHITECTURE_PREFLIGHT_{overall} methods={len(rows)} "
        f"code_ready={str(all_code_ready).lower()} refs_ready={str(all_refs_ready).lower()} "
        f"dense_ready={str(dense_ready).lower()} available_gib={available_gib:.3f}"
    )


if __name__ == "__main__":
    main()
