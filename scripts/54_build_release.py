#!/usr/bin/env python3
"""Build atlas-v1.0 atomically after every non-release acceptance gate passes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import downstream_contract

RELEASE_NAME = "atlas-v1.0"
RELEASE_SCHEMA = "atlas-v1.0-release.1"
SOURCE_FIELDS = [
    "trait_id", "source_id", "dataset_version", "source_page_url", "download_url",
    "access", "archive_name", "archive_bytes", "archive_sha256", "archive_member",
    "raw_file", "pmid", "doi", "ancestry", "build", "phenotype_definition",
    "source_status", "acquisition_status",
]
RESULT_SUFFIXES = {".tsv", ".json", ".txt", ".md", ".html", ".gz", ".rds", ".png", ".pdf"}
FORBIDDEN_RESULT_MARKERS = ("SYNTHETIC SMOKE-TEST OUTPUT", "PLACEHOLDER", "FAKE_RESULT")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def git(root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", *arguments], cwd=root, capture_output=True, text=True, check=False,
    )
    if result.returncode:
        fail(f"git {' '.join(arguments)} failed: {(result.stdout + result.stderr).strip()}")
    return result.stdout.strip()


def relative_files(root: Path, directory: str) -> list[Path]:
    base = root / directory
    if not base.is_dir():
        return []
    return [
        path.relative_to(root) for path in sorted(base.rglob("*"))
        if path.is_file() and path.suffix.lower() in RESULT_SUFFIXES and path.name != ".gitkeep"
    ]


def source_provenance(root: Path) -> list[dict[str, str]]:
    panel = read_tsv(root / "config/analysis_panel.tsv")
    sources = read_tsv(root / "config/public_gwas_sources.tsv")
    by_trait: dict[str, dict[str, str]] = {}
    for source in sources:
        for trait in source["trait_ids"].split(","):
            if trait in by_trait:
                fail(f"multiple public source rows map to {trait}")
            by_trait[trait] = source
    output = []
    for trait in panel:
        identity = trait["trait_id"]
        source = by_trait.get(identity)
        if source is None or source["source_id"] != trait["source_id"]:
            fail(f"source registry does not map exactly to {identity}")
        output.append({
            "trait_id": identity,
            "source_id": trait["source_id"],
            "dataset_version": trait["dataset_version"],
            "source_page_url": source["source_page_url"],
            "download_url": source["download_url"],
            "access": source["access"],
            "archive_name": source["archive_name"],
            "archive_bytes": source["archive_bytes"],
            "archive_sha256": source["archive_sha256"],
            "archive_member": source["archive_member"],
            "raw_file": trait["raw_file"],
            "pmid": trait["pmid"],
            "doi": trait["doi"],
            "ancestry": trait["ancestry"],
            "build": trait["build"],
            "phenotype_definition": trait["phenotype_definition"],
            "source_status": trait["source_status"],
            "acquisition_status": source["acquisition_status"],
        })
    if len(output) != 45:
        fail("source provenance is not the exact 45-trait family")
    return output


def require_non_release_gates(root: Path) -> list[dict[str, str]]:
    with tempfile.TemporaryDirectory() as directory:
        audit_path = Path(directory) / "acceptance.json"
        result = subprocess.run(
            [sys.executable, str(root / "scripts/99_atlas_acceptance.py"),
             "--root", str(root), "--json-out", str(audit_path), "--report-only"],
            cwd=root, capture_output=True, text=True, check=False,
        )
        if result.returncode or not audit_path.is_file():
            fail(f"acceptance audit could not run: {(result.stdout + result.stderr).strip()}")
        gates = json.loads(audit_path.read_text(encoding="utf-8"))
    blockers = [gate for gate in gates if gate["gate"] != "atlas_v1_release" and gate["status"] != "PASS"]
    if blockers:
        fail("non-release gates remain blocked: " + ", ".join(gate["gate"] for gate in blockers))
    if len(gates) != 23 or len(gates) - 1 != 22:
        fail("acceptance gate family has drifted")
    return gates


def copy_payload(root: Path, stage: Path, relative: Path) -> None:
    source = root / relative
    if not source.is_file() or source.stat().st_size == 0:
        fail(f"release payload is absent or empty: {relative}")
    target = stage / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--execute", action="store_true", help="create the immutable release directory")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    target = root / "releases" / RELEASE_NAME
    stage = root / "releases" / f".{RELEASE_NAME}.building"
    if target.exists() or stage.exists():
        fail("release target/staging path already exists; releases are never overwritten")
    if git(root, "status", "--porcelain", "--untracked-files=no"):
        fail("tracked worktree must be clean before freezing a release")
    code_commit = git(root, "rev-parse", "HEAD")
    if len(code_commit) != 40:
        fail("could not resolve an exact producing commit")
    gates = require_non_release_gates(root)
    if not args.execute:
        print(f"RELEASE_READY name={RELEASE_NAME} code_commit={code_commit}; rerun with --execute")
        return 0

    stage.mkdir(parents=True)
    try:
        tracked = git(root, "ls-files").splitlines()
        code_prefixes = ("config/", "environment/", "patches/", "scripts/", "docs/", "tests/", ".github/")
        payload = {
            Path(path) for path in tracked
            if path in {
                "README.md", "Snakefile", ".gitignore",
                "discovery_extension/scripts/35_run_susie_coloc.R",
            } or path.startswith(code_prefixes)
        }
        for directory in (
            "results/atlas", "results/tables", "results/logs", "results/figures",
            "results/sources", "results/fine_mapping", "results/molecular",
            "results/twas", "results/pleiotropy", "results/interpretation",
            "results/robustness",
        ):
            payload.update(relative_files(root, directory))
        for relative in sorted(payload):
            copy_payload(root, stage, relative)

        shutil.copy2(root / "config/analysis_panel.tsv", stage / "analysis_panel.tsv")
        shutil.copy2(root / "environment/tool_versions.tsv", stage / "tool_versions.tsv")
        (stage / "code_commit.txt").write_text(code_commit + "\n", encoding="utf-8")
        write_tsv(stage / "source_provenance.tsv", SOURCE_FIELDS, source_provenance(root))
        (stage / "acceptance_gates.json").write_text(
            json.dumps(gates, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        manifest_files = []
        for path in sorted(stage.rglob("*")):
            if not path.is_file() or path.name in {"checksums.sha256", "release_manifest.json"}:
                continue
            relative = path.relative_to(stage).as_posix()
            if relative.startswith("results/") and path.suffix.lower() in {".tsv", ".txt", ".md", ".json"}:
                content = path.read_text(encoding="utf-8", errors="ignore").upper()
                if any(marker in content for marker in FORBIDDEN_RESULT_MARKERS):
                    fail(f"forbidden non-real marker in release result: {relative}")
            role = "code_or_configuration"
            if relative.startswith("results/"):
                role = "scientific_result_or_qc"
            elif relative == "source_provenance.tsv":
                role = "source_provenance"
            elif relative == "tool_versions.tsv":
                role = "tool_versions"
            elif relative in {"analysis_panel.tsv", "code_commit.txt", "acceptance_gates.json"}:
                role = "release_metadata"
            manifest_files.append({
                "path": relative, "bytes": path.stat().st_size, "sha256": sha256(path), "role": role,
            })
        checksum_text = "".join(f"{item['sha256']}  {item['path']}\n" for item in manifest_files)
        (stage / "checksums.sha256").write_text(checksum_text, encoding="utf-8")
        manifest = {
            "schema_version": RELEASE_SCHEMA,
            "release_name": RELEASE_NAME,
            "analysis_panel": "atlas-v1.0",
            "code_commit": code_commit,
            "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "panel_sha256": sha256(stage / "analysis_panel.tsv"),
            "downstream_policy_sha256": sha256(stage / "config/downstream_analysis_policy.json"),
            "atlas_schema_sha256": sha256(stage / "config/atlas_table_schema.json"),
            "checksums_sha256": sha256(stage / "checksums.sha256"),
            "non_release_acceptance_gates": 22,
            "release_script_sha256": downstream_contract.script_hashes(root, "release"),
            "file_count": len(manifest_files),
            "files": manifest_files,
            "immutability_policy": "The release directory is created atomically and is never overwritten in place.",
        }
        (stage / "release_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        validation = subprocess.run(
            [sys.executable, str(root / "scripts/55_validate_release.py"), "--root", str(root),
             "--release", str(stage), "--quiet"],
            cwd=root, capture_output=True, text=True, check=False,
        )
        if validation.returncode:
            fail(f"staged release validation failed: {(validation.stdout + validation.stderr).strip()}")
        os.replace(stage, target)
    except BaseException:
        # Keep a failed stage for forensic inspection; never silently delete evidence.
        raise
    print(f"RELEASE_FROZEN path={target.relative_to(root)} code_commit={code_commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
