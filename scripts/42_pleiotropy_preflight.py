#!/usr/bin/env python3
"""Audit PLACO+/conjFDR inputs, software, reference, memory, and storage."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MIXER_SPEC = importlib.util.spec_from_file_location(
    "mixer_preflight", ROOT / "scripts/35_mixer_preflight.py"
)
if MIXER_SPEC is None or MIXER_SPEC.loader is None:
    raise RuntimeError("could not load full-summary-statistics helpers")
mixer = importlib.util.module_from_spec(MIXER_SPEC)
MIXER_SPEC.loader.exec_module(mixer)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--json-out", default="results/tables/pleiotropy_preflight.json")
    parser.add_argument("--trait-out", default="results/tables/pleiotropy_input_readiness.tsv")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads(
        (root / "config/pleiotropy_analysis_policy.json").read_text(encoding="utf-8")
    )
    panel = mixer.read_tsv(root / "config/analysis_panel.tsv")
    if len(panel) != 45 or sum(row["domain"] == "sleep" for row in panel) != 12:
        raise SystemExit("ERROR: pleiotropy workflow requires the exact locked panel")
    planned_prefilters = {
        row["trait_id"] for row in mixer.read_tsv(root / "config/hm3_prefilter_plans.tsv")
    }
    trait_rows = []
    for row in panel:
        trait = row["trait_id"]
        data, qc, prefilter = mixer.choose_harmonized(root, trait)
        if not prefilter and trait not in planned_prefilters:
            prefilter = "not supplied"
        full = data.is_file() and qc.is_file() and prefilter == "not supplied"
        trait_rows.append({
            "trait_id": trait,
            "domain": row["domain"],
            "harmonized_file": str(data.relative_to(root)),
            "harmonized_bytes": data.stat().st_size if data.is_file() else 0,
            "prefilter_strategy": prefilter or "UNKNOWN",
            "input_status": "READY_FULL_SUMSTATS" if full else "BLOCKED_FULL_SUMSTATS",
            "blocker": "" if full else "regenerate full post-QC harmonized summary statistics without HapMap3 prefiltering",
        })

    placo = root / ".r-env/share/placo/PLACO_v0.2.0.R"
    pleiofdr = root / "work/pleiofdr/pleiotropy_analysis.m"
    pleiofdr_git = root / "work/pleiofdr/.git"
    pleiofdr_commit = ""
    if pleiofdr_git.is_dir():
        result = subprocess.run(
            ["git", "-C", str(pleiofdr.parent), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            pleiofdr_commit = result.stdout.strip()
    reference = root / "ref/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat"
    template = root / "ref/pleiofdr/9545380.ref"
    memory_bytes = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    free_bytes = shutil.disk_usage(root).free
    matlab = shutil.which("matlab") or ""
    checks = {
        "full_trait_inputs": {
            "pass": all(row["input_status"] == "READY_FULL_SUMSTATS" for row in trait_rows),
            "ready_traits": sum(row["input_status"] == "READY_FULL_SUMSTATS" for row in trait_rows),
            "expected_traits": 45,
            "blocked_traits": [row["trait_id"] for row in trait_rows if row["input_status"] != "READY_FULL_SUMSTATS"],
        },
        "placo_source": {
            "pass": placo.is_file() and sha256(placo) == policy["placo_source_sha256"],
            "path": str(placo.relative_to(root)),
        },
        "pleiofdr_code": {
            "pass": pleiofdr.is_file() and pleiofdr_commit == policy["pleiofdr_commit"],
            "path": str(pleiofdr.relative_to(root)),
            "observed_commit": pleiofdr_commit or "ABSENT",
            "required_commit": policy["pleiofdr_commit"],
        },
        "matlab": {
            "pass": bool(matlab),
            "observed": matlab or "NONE",
            "note": "Octave is experimental upstream and is not accepted for the primary atlas run",
        },
        "pleiofdr_reference": {
            "pass": reference.is_file() and reference.stat().st_size == policy["pleiofdr_reference_bytes"],
            "path": str(reference.relative_to(root)),
            "observed_bytes": reference.stat().st_size if reference.is_file() else 0,
            "expected_bytes": policy["pleiofdr_reference_bytes"],
        },
        "pleiofdr_variant_template": {
            "pass": template.is_file()
            and template.stat().st_size == policy["pleiofdr_variant_template_bytes"]
            and sha256(template) == policy["pleiofdr_variant_template_sha256"],
            "path": str(template.relative_to(root)),
            "observed_bytes": template.stat().st_size if template.is_file() else 0,
            "expected_bytes": policy["pleiofdr_variant_template_bytes"],
        },
        "scipy_mat_writer": {
            "pass": importlib.util.find_spec("scipy") is not None,
            "note": "SciPy is required to materialize deterministic MATLAB trait vectors",
        },
        "memory": {
            "pass": memory_bytes >= policy["pleiofdr_minimum_memory_bytes"],
            "observed_bytes": memory_bytes,
            "minimum_bytes": policy["pleiofdr_minimum_memory_bytes"],
        },
        "storage": {
            "pass": free_bytes >= policy["minimum_free_storage_bytes"],
            "observed_free_bytes": free_bytes,
            "minimum_free_bytes": policy["minimum_free_storage_bytes"],
        },
    }
    overall = all(check["pass"] for check in checks.values())
    report = {
        "analysis_id": policy["analysis_id"],
        "ready": overall,
        "planned_pair_scans": policy["expected_sleep_non_sleep_pairs"],
        "ready_pair_scans": 12 * sum(
            row["domain"] != "sleep" and row["input_status"] == "READY_FULL_SUMSTATS"
            for row in trait_rows
        ),
        "checks": checks,
    }
    write_tsv(root / args.trait_out, trait_rows)
    output = root / args.json_out
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    for name, check in checks.items():
        print(f"{'PASS' if check['pass'] else 'BLOCKED':7} {name}: {check}")
    print(f"Pleiotropy preflight: {'READY' if overall else 'BLOCKED'}")
    return 0 if overall or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
