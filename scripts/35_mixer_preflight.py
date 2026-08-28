#!/usr/bin/env python3
"""Audit MiXeR inputs, reference files, storage, and supported hardware."""
from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
            return data, qc, qc_value(qc, "prefilter_strategy")
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
    parser.add_argument("--json-out", default="results/tables/mixer_preflight.json")
    parser.add_argument("--trait-out", default="results/tables/mixer_input_readiness.tsv")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads((root / "config/mixer_analysis_policy.json").read_text(encoding="utf-8"))
    panel = read_tsv(root / "config/analysis_panel.tsv")
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
        elif data.is_file() and prefilter == "HAPMAP3_RSID_ALLOWLIST":
            status = "BLOCKED_HAPMAP3_PREFILTERED"
            blocker = "regenerate full post-QC summary statistics without a HapMap3 allowlist"
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

    reference_root = root / policy["reference_root"]
    reference_files = []
    for chromosome in range(1, 23):
        reference_files.extend([
            reference_root / policy["reference_bim_pattern"].replace("@", str(chromosome)),
            reference_root / policy["reference_ld_pattern"].replace("@", str(chromosome)),
        ])
    for replicate in range(1, policy["fit_replicates"] + 1):
        reference_files.append(reference_root / policy["reference_extract_pattern"].replace("@", str(replicate)))

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
        "reference": {
            "pass": all(path.is_file() and path.stat().st_size for path in reference_files),
            "present_files": sum(path.is_file() and path.stat().st_size > 0 for path in reference_files),
            "expected_files": len(reference_files),
            "first_missing": next((str(path.relative_to(root)) for path in reference_files if not path.is_file() or not path.stat().st_size), ""),
        },
        "full_trait_inputs": {
            "pass": all(row["input_status"] == "READY_FULL_SUMSTATS" for row in trait_rows),
            "ready_traits": sum(row["input_status"] == "READY_FULL_SUMSTATS" for row in trait_rows),
            "expected_traits": 45,
            "blocked_traits": [row["trait_id"] for row in trait_rows if row["input_status"] != "READY_FULL_SUMSTATS"],
        },
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
