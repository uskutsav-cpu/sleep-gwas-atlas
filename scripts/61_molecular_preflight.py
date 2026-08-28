#!/usr/bin/env python3
"""Fail-closed readiness audit for molecular-QTL colocalization and TWAS."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


SOURCE_FIELDS = [
    "source_family_id", "modality", "study_or_model", "release_or_record",
    "accession_or_dataset", "build", "ancestry_or_population", "context_scope",
    "access_mode", "official_url", "asset_url", "asset_bytes", "asset_checksum",
    "access_terms", "primary_citation", "role", "availability_policy",
]
EXPECTED_SOURCES = {
    "EQTL_CATALOGUE_R7_GE", "EQTL_CATALOGUE_R7_LEAFCUTTER",
    "EQTL_CATALOGUE_R7_SUN2018", "EQTL_CATALOGUE_GTEX_V8_IMPORTED",
    "FENLAND_CIS_PQTL_V1", "PSYCHENCODE_R3_BRAIN",
    "PREDICTDB_GTEX_V8_ELASTIC_NET_PHI_EQTL",
    "METAXCAN_V081",
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def memory_gib() -> float:
    try:
        output = subprocess.check_output(
            ["sysctl", "-n", "hw.memsize"], text=True, stderr=subprocess.DEVNULL,
        ).strip()
        return int(output) / 1024**3
    except (OSError, ValueError, subprocess.CalledProcessError):
        try:
            pages = os.sysconf("SC_PHYS_PAGES")
            page_size = os.sysconf("SC_PAGE_SIZE")
            return pages * page_size / 1024**3
        except (AttributeError, OSError, ValueError):
            return 0.0


def r_package_version(root: Path, package: str) -> str:
    rscript = root / ".r-env/bin/Rscript"
    if not rscript.is_file():
        return "MISSING"
    result = subprocess.run(
        [str(rscript), "-e", f"cat(as.character(packageVersion('{package}')))"],
        capture_output=True, text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "MISSING"


def runtime_versions(python: Path, packages: dict[str, str]) -> tuple[bool, str]:
    code = """
import importlib.metadata
import json
import platform
import sys
names = json.loads(sys.argv[1])
observed = {"python": platform.python_version()}
for name in names:
    if name == "python":
        continue
    try:
        observed[name] = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        observed[name] = "MISSING"
print(json.dumps(observed, sort_keys=True))
"""
    result = subprocess.run(
        [str(python), "-c", code, json.dumps(packages, sort_keys=True)],
        capture_output=True, text=True,
    )
    if result.returncode:
        return False, (result.stdout + result.stderr).strip()[-500:]
    try:
        observed = json.loads(result.stdout)
    except json.JSONDecodeError:
        return False, result.stdout.strip()[-500:]
    matches = all(
        observed.get(name) == expected or observed.get(name, "").startswith(expected + ".")
        for name, expected in packages.items()
    )
    return matches, json.dumps(observed, sort_keys=True, separators=(",", ":"))


def qc_prefiltered(path: Path) -> bool:
    return any(line.startswith("prefilter_strategy\t") for line in path.read_text(encoding="utf-8").splitlines())


def choose_full_input(root: Path, trait: str) -> Path | None:
    candidates = [
        (
            root / f"data/harmonized_mixer_full/{trait}.harmonized.tsv.gz",
            root / f"data/harmonized_mixer_full/{trait}.qc.txt",
        ),
        (
            root / f"data/harmonized/{trait}.harmonized.tsv.gz",
            root / f"data/harmonized/{trait}.qc.txt",
        ),
    ]
    for harmonized, qc in candidates:
        if harmonized.is_file() and qc.is_file() and not qc_prefiltered(qc):
            return harmonized
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--sources", default="config/molecular_source_registry.tsv")
    parser.add_argument("--readiness-out", default="results/tables/molecular_input_readiness.tsv")
    parser.add_argument("--out", default="results/tables/molecular_preflight.json")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, sources_path = root / args.policy, root / args.sources
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    source_fields, sources = read_tsv(sources_path)
    if source_fields != SOURCE_FIELDS:
        fail("molecular source registry header differs from the locked schema")
    source_ids = [row["source_family_id"] for row in sources]
    if set(source_ids) != EXPECTED_SOURCES or len(source_ids) != len(set(source_ids)):
        fail("molecular source registry is incomplete or duplicated")
    if policy["source_release_policy"]["eqtl_catalogue_release"] != "r7_June_2024":
        fail("stable eQTL Catalogue r7 is no longer selected")
    if "r8_pre-release_January_2026" not in policy["source_release_policy"]["reject_prerelease"]:
        fail("r8 pre-release exclusion is absent")
    for row in sources:
        if not row["official_url"].startswith("https://") or not row["asset_url"].startswith("https://"):
            fail(f"non-HTTPS molecular source route: {row['source_family_id']}")
        if not row["primary_citation"].startswith("10.") or not row["availability_policy"].strip():
            fail(f"molecular source lacks citation or availability policy: {row['source_family_id']}")

    rows: list[dict[str, str]] = []
    metadata_ready = True
    for asset in policy["metadata_assets"]:
        path = root / asset["path"]
        ok = path.is_file() and path.stat().st_size == asset["bytes"] and sha256(path) == asset["sha256"]
        metadata_ready &= ok
        rows.append({
            "component": asset["id"], "required": "YES", "status": "READY" if ok else "MISSING_OR_HASH_MISMATCH",
            "observed": str(path.relative_to(root)) if path.is_file() else "ABSENT",
            "requirement": f"bytes={asset['bytes']};sha256={asset['sha256']}",
        })

    chain = policy["reference_build"]
    chain_path = root / chain["chain_path"]
    chain_ok = (
        chain_path.is_file() and chain_path.stat().st_size == chain["chain_bytes"]
        and sha256(chain_path) == chain["chain_sha256"]
    )
    rows.append({
        "component": "GRCH37_TO_GRCH38_CHAIN", "required": "YES",
        "status": "READY" if chain_ok else "MISSING_OR_HASH_MISMATCH",
        "observed": str(chain_path.relative_to(root)) if chain_path.is_file() else "ABSENT",
        "requirement": f"bytes={chain['chain_bytes']};sha256={chain['chain_sha256']}",
    })
    reverse_chain = policy["source_to_analysis_chain"]
    reverse_chain_path = root / reverse_chain["path"]
    reverse_chain_ok = (
        reverse_chain_path.is_file() and reverse_chain_path.stat().st_size == reverse_chain["bytes"]
        and sha256(reverse_chain_path) == reverse_chain["sha256"]
    )
    rows.append({
        "component": "GRCH38_TO_GRCH37_CHAIN", "required": "YES",
        "status": "READY" if reverse_chain_ok else "MISSING_OR_HASH_MISMATCH",
        "observed": str(reverse_chain_path.relative_to(root)) if reverse_chain_path.is_file() else "ABSENT",
        "requirement": f"bytes={reverse_chain['bytes']};sha256={reverse_chain['sha256']}",
    })
    engine = policy["generic_engine"]
    engine_path = root / engine["path"]
    engine_ok = engine_path.is_file() and sha256(engine_path) == engine["sha256"]
    rows.append({
        "component": "SUSIE_COLOC_ENGINE", "required": "YES",
        "status": "READY" if engine_ok else "MISSING_OR_HASH_MISMATCH",
        "observed": str(engine_path.relative_to(root)) if engine_path.is_file() else "ABSENT",
        "requirement": engine["sha256"],
    })
    susie, coloc = r_package_version(root, "susieR"), r_package_version(root, "coloc")
    r_ok = susie == policy["fine_mapping"]["susieR_version"] and coloc == policy["colocalization"]["coloc_version"]
    rows.append({
        "component": "R_PACKAGES", "required": "YES", "status": "READY" if r_ok else "MISSING_OR_VERSION_MISMATCH",
        "observed": f"susieR={susie};coloc={coloc}",
        "requirement": f"susieR={policy['fine_mapping']['susieR_version']};coloc={policy['colocalization']['coloc_version']}",
    })
    tabix = shutil.which("tabix")
    rows.append({
        "component": "TABIX", "required": "YES", "status": "READY" if tabix else "MISSING",
        "observed": tabix or "ABSENT", "requirement": "tabix executable",
    })
    twas = policy["twas"]
    metaxcan = root / twas["entrypoint_path"]
    metaxcan_archive = root / twas["source_archive_path"]
    metaxcan_ok = (
        metaxcan.is_file() and sha256(metaxcan) == twas["entrypoint_sha256"]
        and metaxcan_archive.is_file() and metaxcan_archive.stat().st_size == twas["source_archive_bytes"]
        and sha256(metaxcan_archive) == twas["source_archive_sha256"]
    )
    rows.append({
        "component": "METAXCAN_V081", "required": "YES", "status": "READY" if metaxcan_ok else "MISSING",
        "observed": str(metaxcan.relative_to(root)) if metaxcan_ok else "ABSENT",
        "requirement": policy["twas"]["metaxcan_commit"],
    })
    metaxcan_python = root / twas["runtime_environment"] / "bin/python"
    runtime_ok = False
    runtime_detail = "ABSENT"
    if metaxcan_python.is_file() and metaxcan_ok:
        versions_ok, version_detail = runtime_versions(metaxcan_python, twas["runtime_versions"])
        smoke = subprocess.run(
            [str(metaxcan_python), str(metaxcan), "--help"], capture_output=True, text=True,
        )
        runtime_ok = versions_ok and smoke.returncode == 0
        runtime_detail = f"versions={version_detail};help_returncode={smoke.returncode}"
    rows.append({
        "component": "METAXCAN_RUNTIME", "required": "YES", "status": "READY" if runtime_ok else "MISSING_OR_VERSION_MISMATCH",
        "observed": str(metaxcan_python.relative_to(root)) if metaxcan_python.is_file() else "ABSENT",
        "requirement": json.dumps(twas["runtime_versions"], sort_keys=True, separators=(",", ":")),
    })
    model_source = twas["phi_model_source"]
    inventory_path = root / model_source["inventory_path"]
    inventory_lock_path = root / model_source["inventory_lock_path"]
    inventory_ok = False
    inventory_rows: list[dict[str, str]] = []
    if inventory_path.is_file() and inventory_lock_path.is_file():
        try:
            _, inventory_rows = read_tsv(inventory_path)
            inventory_lock = json.loads(inventory_lock_path.read_text(encoding="utf-8"))
            inventory_ok = (
                inventory_lock.get("inventory_sha256") == sha256(inventory_path)
                and inventory_lock.get("policy_sha256") == sha256(policy_path)
                and inventory_lock.get("file_ids_in_locked_order") == [row["file_id"] for row in inventory_rows]
                and len(inventory_rows) == model_source["expected_file_count"]
                and sum(int(row["bytes"]) for row in inventory_rows) == model_source["expected_total_bytes"]
            )
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            inventory_ok = False
    rows.append({
        "component": "PREDICTDB_PHI_MODEL_INVENTORY", "required": "YES",
        "status": "READY" if inventory_ok else "MISSING_OR_HASH_MISMATCH",
        "observed": f"{len(inventory_rows)}/{model_source['expected_file_count']}",
        "requirement": f"result-free Box inventory files={model_source['expected_file_count']};bytes={model_source['expected_total_bytes']}",
    })
    download_lock_path = root / model_source["download_lock_path"]
    model_ready = False
    ready_model_files = 0
    if inventory_ok and download_lock_path.is_file():
        try:
            download_lock = json.loads(download_lock_path.read_text(encoding="utf-8"))
            model_ready = (
                download_lock.get("inventory_sha256") == sha256(inventory_path)
                and download_lock.get("inventory_lock_sha256") == sha256(inventory_lock_path)
                and download_lock.get("policy_sha256") == sha256(policy_path)
            )
            for record in inventory_rows:
                path = root / model_source["install_dir"] / record["filename"]
                evidence = download_lock.get("files", {}).get(record["filename"], {})
                ok = (
                    path.is_file() and path.stat().st_size == int(record["bytes"])
                    and evidence.get("file_id") == record["file_id"]
                    and evidence.get("bytes") == int(record["bytes"])
                    and evidence.get("sha256") == sha256(path)
                )
                ready_model_files += int(ok)
                model_ready &= ok
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            model_ready = False
    rows.append({
        "component": "PREDICTDB_PHI_MODELS", "required": "YES",
        "status": "READY" if model_ready else "MISSING_OR_HASH_MISMATCH",
        "observed": f"{ready_model_files}/{model_source['expected_file_count']}",
        "requirement": "all phi-enabled model databases and covariances byte- and SHA-256-locked",
    })
    full = []
    for trait in read_tsv(root / "config/analysis_panel.tsv")[1]:
        found = choose_full_input(root, trait["trait_id"])
        full.append((trait["trait_id"], found))
    full_count = sum(path is not None for _, path in full)
    rows.append({
        "component": "FULL_DENSE_GWAS", "required": "YES", "status": "READY" if full_count == 45 else "INCOMPLETE",
        "observed": f"{full_count}/45", "requirement": "45/45 full dense canonical inputs",
    })
    upstream = [
        root / "results/atlas/loci.tsv", root / "results/atlas/variants.tsv",
        root / "results/tables/trait_trait_colocalization.tsv", root / "results/atlas/fine_mapping.provenance.json",
    ]
    upstream_ok = all(path.is_file() and path.stat().st_size > 0 for path in upstream)
    rows.append({
        "component": "UPSTREAM_FINE_MAPPING", "required": "YES", "status": "READY" if upstream_ok else "BLOCKED_UPSTREAM",
        "observed": ";".join(str(path.relative_to(root)) for path in upstream if path.is_file()) or "ABSENT",
        "requirement": "locked loci variants trait-trait colocalization and provenance",
    })
    ram, free = memory_gib(), shutil.disk_usage(root).free / 1024**3
    resources_ok = ram >= policy["runtime"]["minimum_memory_gib"] and free >= policy["runtime"]["minimum_free_disk_gib_for_production"]
    rows.append({
        "component": "LOCAL_RESOURCES", "required": "YES", "status": "READY" if resources_ok else "INSUFFICIENT",
        "observed": f"memory_gib={ram:.2f};free_disk_gib={free:.2f}",
        "requirement": f"memory_gib>={policy['runtime']['minimum_memory_gib']};free_disk_gib>={policy['runtime']['minimum_free_disk_gib_for_production']}",
    })

    required_ready = all(row["status"] == "READY" for row in rows if row["required"] == "YES")
    code_ready = metadata_ready and chain_ok and reverse_chain_ok and engine_ok and r_ok
    status = "READY_FOR_PRODUCTION" if required_ready else "PROTOCOL_READY_INPUTS_OR_RUNTIME_BLOCKED" if code_ready else "PROTOCOL_DEPENDENCIES_BLOCKED"
    readiness_out = root / args.readiness_out
    readiness_out.parent.mkdir(parents=True, exist_ok=True)
    with readiness_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["component", "required", "status", "observed", "requirement"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    provenance = {
        "analysis_id": policy["analysis_id"], "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall_status": status, "code_ready": code_ready, "production_ready": required_ready,
        "policy_sha256": sha256(policy_path), "sources_sha256": sha256(sources_path),
        "readiness_sha256": sha256(readiness_out), "full_dense_gwas_count": full_count,
        "metadata_ready": metadata_ready, "chain_ready": chain_ok and reverse_chain_ok, "engine_ready": engine_ok,
        "tabix_ready": bool(tabix), "metaxcan_source_ready": metaxcan_ok,
        "metaxcan_runtime_ready": runtime_ok, "metaxcan_runtime_detail": runtime_detail,
        "predictdb_models_ready": model_ready,
        "upstream_fine_mapping_ready": upstream_ok, "memory_gib": round(ram, 3), "free_disk_gib": round(free, 3),
        "results_accessed": False, "claim_limit": policy["claim_limit"],
    }
    out = root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not args.quiet:
        print(f"MOLECULAR_PREFLIGHT_{status} full_dense={full_count}/45 free_gib={free:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
