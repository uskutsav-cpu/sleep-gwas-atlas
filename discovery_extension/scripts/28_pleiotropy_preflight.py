#!/usr/bin/env python3
"""Audit pinned PLACO+ code and fail closed on missing real inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/pleiotropy_sources.tsv"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/pleiotropy_contract.json"))
    parser.add_argument("--rscript", type=Path, default=Path(".r-env/bin/Rscript"))
    parser.add_argument("--priority", type=Path, default=Path("discovery_extension/results/prioritization/novel_hit_priority.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("discovery_extension/results/replication/replication_results.tsv"))
    parser.add_argument("--dense-lock", type=Path, default=Path("discovery_extension/provenance/dense_local_inputs.lock.json"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/pleiotropy/pleiotropy_readiness.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/pleiotropy_dependencies.json"))
    args = parser.parse_args()

    rows = read_tsv(args.sources)
    if len(rows) != 1 or rows[0].get("method_id") != "PLACO_PLUS_PRIMARY":
        raise SystemExit("ERROR: pleiotropy source contract must contain the exact PLACO+ primary row")
    source = rows[0]
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("primary_method") != source["method_id"] or contract.get("software_version") != source["software_version"]:
        raise SystemExit("ERROR: PLACO+ method/version differs between source and analysis contracts")
    installed = Path(source["installed_source_path"])
    source_hash_ok = installed.is_file() and sha256(installed) == source["source_sha256"]
    entrypoints_ok = False
    smoke_detail = "installed source missing or SHA256 mismatch"
    if source_hash_ok and args.rscript.is_file():
        expression = (
            f"source('{installed.as_posix()}');"
            "cat(all(vapply(c('var.placo','cor.pearson','placo.plus'),exists,logical(1),inherits=TRUE)))"
        )
        smoke = subprocess.run([str(args.rscript), "-e", expression], check=False, text=True, capture_output=True)
        entrypoints_ok = smoke.returncode == 0 and smoke.stdout.strip().endswith("TRUE")
        smoke_detail = (smoke.stdout + smoke.stderr).strip()
    code_ready = source_hash_ok and entrypoints_ok
    dense_ready = False
    dense_detail = f"dense input lock absent: {args.dense_lock}"
    if args.dense_lock.is_file():
        try:
            dense = json.loads(args.dense_lock.read_text(encoding="utf-8"))
            dense_ready = dense.get("status") == "VERIFIED" and dense.get("full_resolution") is True
            dense_detail = "verified full-resolution dense input lock" if dense_ready else "dense lock is not VERIFIED/full_resolution"
        except (json.JSONDecodeError, OSError) as exc:
            dense_detail = f"invalid dense input lock: {exc}"
    priority_ready, replication_ready = args.priority.is_file(), args.replication.is_file()
    if not code_ready:
        readiness = "BLOCKED_CODE"
    elif not priority_ready or not replication_ready or not dense_ready:
        readiness = "CODE_READY_UPSTREAM_DISCOVERY_REPLICATION_AND_DENSE_INPUTS_BLOCKED"
    else:
        readiness = "READY_FOR_PAIR_MANIFEST_CURATION"
    output = [{
        "method_id": source["method_id"], "software_version": source["software_version"],
        "source_commit": source["source_commit"], "installed_source_path": str(installed),
        "source_sha256": source["source_sha256"], "source_hash_status": "PASS" if source_hash_ok else "FAIL",
        "entrypoints_status": "PASS" if entrypoints_ok else "FAIL", "code_status": "PASS" if code_ready else "FAIL",
        "priority_status": "PRESENT" if priority_ready else "BLOCKED_UPSTREAM",
        "replication_status": "PRESENT" if replication_ready else "BLOCKED_UPSTREAM",
        "dense_input_status": "VERIFIED" if dense_ready else "BLOCKED_UPSTREAM",
        "dense_input_detail": dense_detail, "readiness": readiness,
    }]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(output[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0", "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall_status": readiness, "download_started": False, "analysis_started": False,
        "sources_sha256": sha256(args.sources), "contract_sha256": sha256(args.contract),
        "installed_source_sha256": sha256(installed) if installed.is_file() else None,
        "runtime_smoke_detail": smoke_detail, "priority_present": priority_ready,
        "replication_present": replication_ready, "dense_input_verified": dense_ready,
        "output": str(args.out), "output_sha256": sha256(args.out),
        "warning": "PLACO+ code readiness is not a pleiotropy result; genome-wide inputs and the replicated pair family are absent.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PLEIOTROPY_PREFLIGHT_{readiness} code_ready={str(code_ready).lower()} dense_ready={str(dense_ready).lower()}")


if __name__ == "__main__":
    main()
