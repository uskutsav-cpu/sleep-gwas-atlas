#!/usr/bin/env python3
"""Audit exact fine-mapping code, full inputs, and LAVA signed-LD readiness."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_METHODS = {"SUSIE_RSS_PRIMARY", "COLOC_SUSIE_PRIMARY"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def qc_prefiltered(path: Path) -> bool:
    return any(line.startswith("prefilter_strategy\t") for line in path.read_text(encoding="utf-8").splitlines())


def choose_full_input(root: Path, trait: str) -> tuple[Path, Path, bool]:
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
    for harmonized, qc in candidates:
        if harmonized.is_file() and qc.is_file() and not qc_prefiltered(qc):
            return harmonized, qc, True
    harmonized, qc = candidates[-1]
    return harmonized, qc, False


def memory_bytes() -> int:
    try:
        result = subprocess.run(
            ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=False,
        )
        if result.returncode == 0:
            return int(result.stdout.strip())
    except (OSError, ValueError):
        pass
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        size = os.sysconf("SC_PAGE_SIZE")
        return int(pages * size)
    except (AttributeError, OSError, ValueError):
        return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/fine_mapping_analysis_policy.json")
    parser.add_argument("--sources", default="config/fine_mapping_sources.tsv")
    parser.add_argument("--report", default="results/tables/fine_mapping_preflight.json")
    parser.add_argument("--traits-out", default="results/tables/fine_mapping_input_readiness.tsv")
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    source_path = root / args.sources
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    downstream = json.loads((root / "config/downstream_analysis_policy.json").read_text(encoding="utf-8"))
    if (
        policy["analysis_panel"] != "atlas-v1.0"
        or policy["fine_mapping"]["version"] != downstream["fine_mapping"]["susieR_version"]
        or policy["colocalization"]["version"] != downstream["colocalization"]["coloc_version"]
        or policy["minimum_locus_variants"] != downstream["fine_mapping"]["minimum_locus_variants"]
    ):
        raise SystemExit("ERROR: fine-mapping policy differs from the locked downstream contract")
    sources = read_tsv(source_path)
    if len(sources) != 2 or {row["method_id"] for row in sources} != EXPECTED_METHODS:
        raise SystemExit("ERROR: fine-mapping source registry is not the exact two-method family")

    code_checks = []
    for source in sources:
        archive = root / source["installed_source_archive_path"]
        description = root / source["installed_description_path"]
        archive_ok = archive.is_file() and sha256(archive) == source["source_archive_sha256"]
        description_ok = description.is_file() and f"Version: {source['software_version']}" in description.read_text(
            encoding="utf-8", errors="ignore"
        )
        code_checks.append({
            "method_id": source["method_id"], "archive_present_and_pinned": archive_ok,
            "installed_description_version_pinned": description_ok,
        })
    engine = root / policy["shared_engine"]["path"]
    engine_ok = engine.is_file() and sha256(engine) == policy["shared_engine"]["sha256"]
    rscript = root / ".r-env/bin/Rscript"
    runtime_ok = False
    runtime_detail = "Rscript absent"
    if rscript.is_file():
        expression = (
            "suppressPackageStartupMessages(library(susieR));"
            "suppressPackageStartupMessages(library(coloc));"
            "ok<-as.character(packageVersion('susieR'))=='0.14.2' && "
            "as.character(packageVersion('coloc'))=='5.2.3' && "
            "all(vapply(c('susie_rss','susie_get_cs','estimate_s_rss','kriging_rss'),exists,logical(1),where=asNamespace('susieR'),inherits=FALSE)) && "
            "all(vapply(c('coloc.susie','coloc.abf','runsusie','sensitivity'),exists,logical(1),where=asNamespace('coloc'),inherits=FALSE));cat(ok)"
        )
        smoke = subprocess.run([str(rscript), "-e", expression], capture_output=True, text=True, check=False)
        runtime_ok = smoke.returncode == 0 and smoke.stdout.strip().endswith("TRUE")
        runtime_detail = (smoke.stdout + smoke.stderr).strip()

    panel = read_tsv(root / "config/analysis_panel.tsv")
    trait_rows = []
    for row in panel:
        trait = row["trait_id"]
        harmonized, qc, full = choose_full_input(root, trait)
        trait_rows.append({
            "trait_id": trait,
            "harmonized_path": str(harmonized.relative_to(root)),
            "harmonized_present": str(harmonized.is_file()).upper(),
            "qc_present": str(qc.is_file()).upper(),
            "input_scope": "FULL_AUTOSOMAL_POST_QC" if full else "HAPMAP3_PREFILTERED_BLOCKED",
            "ready": str(full).upper(),
        })
    full_count = sum(row["ready"] == "TRUE" for row in trait_rows)
    if len(trait_rows) != 45:
        raise SystemExit("ERROR: fine-mapping preflight does not cover the exact panel")

    prefix = root / policy["reference"]["prefix"]
    reference_paths = [
        Path(f"{prefix}_chr{chromosome}.{suffix}")
        for chromosome in range(1, 23) for suffix in ("info", "bcor")
    ]
    extracted_manifest = prefix.parent / "extracted_manifest.tsv"
    reference_count = 0
    if extracted_manifest.is_file():
        extracted_rows = read_tsv(extracted_manifest)
        expected_keys = {(str(chromosome), suffix) for chromosome in range(1, 23) for suffix in ("info", "bcor")}
        observed_keys = {(row.get("chromosome", ""), row.get("file_type", "")) for row in extracted_rows}
        if len(extracted_rows) == 44 and observed_keys == expected_keys:
            for row in extracted_rows:
                path = root / row["path"]
                if (
                    path.is_file() and row["bytes"].isdigit()
                    and path.stat().st_size == int(row["bytes"])
                    and sha256(path) == row["sha256"]
                ):
                    reference_count += 1
    shared = root / policy["entry_source"]
    disk = os.statvfs(root)
    free_bytes = disk.f_bavail * disk.f_frsize
    observed_memory = memory_bytes()
    code_ready = all(all(check.values()) for check in code_checks) and engine_ok and runtime_ok
    ready = (
        code_ready and full_count == 45 and reference_count == 44 and shared.is_file()
        and shared.stat().st_size > 0 and free_bytes >= policy["minimum_free_storage_bytes"]
        and observed_memory >= policy["minimum_memory_bytes"]
    )
    blockers = []
    if not code_ready:
        blockers.append("pinned SuSiE/coloc code or runtime absent")
    if full_count != 45:
        blockers.append(f"full-resolution inputs {full_count}/45")
    if reference_count != 44:
        blockers.append(f"LAVA signed-LD files {reference_count}/44")
    if not shared.is_file() or shared.stat().st_size == 0:
        blockers.append("cross-method shared_loci.tsv absent")
    if free_bytes < policy["minimum_free_storage_bytes"]:
        blockers.append(f"free storage {free_bytes} < {policy['minimum_free_storage_bytes']}")
    if observed_memory < policy["minimum_memory_bytes"]:
        blockers.append(f"memory {observed_memory} < {policy['minimum_memory_bytes']}")

    traits_out = root / args.traits_out
    traits_out.parent.mkdir(parents=True, exist_ok=True)
    with traits_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(trait_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(trait_rows)
    report = {
        "analysis_id": policy["analysis_id"],
        "checked_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "ready": ready,
        "blockers": blockers,
        "machine": {"architecture": platform.machine(), "memory_bytes": observed_memory, "free_storage_bytes": free_bytes},
        "code_checks": code_checks,
        "shared_engine": {"path": policy["shared_engine"]["path"], "checksum_pass": engine_ok},
        "runtime": {"pass": runtime_ok, "detail": runtime_detail},
        "full_input_traits": full_count,
        "lava_reference_files": reference_count,
        "lava_extracted_manifest": str(extracted_manifest.relative_to(root)),
        "lava_extracted_manifest_sha256": sha256(extracted_manifest) if extracted_manifest.is_file() else "ABSENT",
        "shared_loci_present": shared.is_file() and shared.stat().st_size > 0 if shared.exists() else False,
        "policy_sha256": sha256(policy_path),
        "sources_sha256": sha256(source_path),
        "trait_readiness_sha256": sha256(traits_out),
        "analysis_started": False,
    }
    report_path = root / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"FINEMAPPING_PREFLIGHT_{'PASS' if ready else 'BLOCKED'} "
        f"code={str(code_ready).lower()} full_inputs={full_count}/45 ld={reference_count}/44"
    )
    if blockers:
        for blocker in blockers:
            print(f"BLOCKER\t{blocker}")
    return 0 if ready or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
