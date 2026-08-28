#!/usr/bin/env python3
"""Validate atlas-v1.0 payload identity, checksums, provenance, and schemas."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import downstream_contract

RELEASE_SCHEMA = "atlas-v1.0-release.1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_TOP_LEVEL = {
    "analysis_panel.tsv", "checksums.sha256", "source_provenance.tsv",
    "tool_versions.tsv", "code_commit.txt", "release_manifest.json",
}
REQUIRED_TOOLS = {
    "analysis_panel", "python_workflow", "snakemake", "python_ldsc", "ldsc_CBIIT",
    "plink", "R", "GenomicSEM", "LAVA", "MiXeR", "PLACO+", "pleioFDR",
    "susieR", "coloc", "MetaXcan", "MAGMA", "TwoSampleMR", "MR-PRESSO",
    "CAUSE_or_LHC_MR",
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
        rows = list(reader)
        return reader.fieldnames or [], rows


def validate_command(command: list[str], cwd: Path, label: str) -> None:
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    if result.returncode:
        fail(f"{label} failed: {(result.stdout + result.stderr).strip()}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository containing validator code")
    parser.add_argument("--release", default="releases/atlas-v1.0")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    release_arg = Path(args.release)
    release = release_arg.resolve() if release_arg.is_absolute() else (root / release_arg).resolve()
    if not release.is_dir():
        fail(f"release directory is absent: {release}")
    top_level = {path.name for path in release.iterdir() if path.is_file()}
    missing = sorted(REQUIRED_TOP_LEVEL - top_level)
    if missing:
        fail("release metadata is incomplete: " + ", ".join(missing))
    manifest_path = release / "release_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"invalid release manifest: {exc}")
    if manifest.get("schema_version") != RELEASE_SCHEMA or manifest.get("release_name") != "atlas-v1.0":
        fail("release name/schema differs from atlas-v1.0")
    code_commit = (release / "code_commit.txt").read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"[0-9a-f]{40}", code_commit) or manifest.get("code_commit") != code_commit:
        fail("producing commit identity is invalid or inconsistent")
    git_check = subprocess.run(
        ["git", "cat-file", "-e", f"{code_commit}^{{commit}}"], cwd=root,
        capture_output=True, text=True, check=False,
    )
    if git_check.returncode:
        fail("producing commit is not present in this repository")
    for field, relative in (
        ("panel_sha256", "analysis_panel.tsv"),
        ("downstream_policy_sha256", "config/downstream_analysis_policy.json"),
        ("atlas_schema_sha256", "config/atlas_table_schema.json"),
        ("checksums_sha256", "checksums.sha256"),
    ):
        expected = manifest.get(field, "")
        if not SHA256.fullmatch(expected) or sha256(release / relative) != expected:
            fail(f"manifest {field} does not match {relative}")
    if (release / "analysis_panel.tsv").read_bytes() != (release / "config/analysis_panel.tsv").read_bytes():
        fail("top-level and configured analysis panels differ")

    entries = manifest.get("files")
    if not isinstance(entries, list) or manifest.get("file_count") != len(entries) or not entries:
        fail("release file manifest is empty or count-inconsistent")
    by_path: dict[str, dict[str, object]] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            fail("release file entry is not an object")
        relative = entry.get("path")
        if not isinstance(relative, str) or relative.startswith("/") or ".." in Path(relative).parts:
            fail("unsafe release path in manifest")
        if relative in by_path:
            fail(f"duplicate release manifest path: {relative}")
        path = release / relative
        if not path.is_file() or path.stat().st_size != entry.get("bytes"):
            fail(f"release file missing or byte count differs: {relative}")
        expected = entry.get("sha256", "")
        if not isinstance(expected, str) or not SHA256.fullmatch(expected) or sha256(path) != expected:
            fail(f"release file checksum differs: {relative}")
        by_path[relative] = entry
    actual = {
        path.relative_to(release).as_posix() for path in release.rglob("*")
        if path.is_file() and path.name not in {"checksums.sha256", "release_manifest.json"}
    }
    if actual != set(by_path):
        fail("release contains unmanifested files or omits payload files")
    checksum_rows = []
    for line in (release / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        if "  " not in line:
            fail("malformed checksums.sha256 line")
        checksum, relative = line.split("  ", 1)
        checksum_rows.append((relative, checksum))
    if checksum_rows != [(relative, by_path[relative]["sha256"]) for relative in sorted(by_path)]:
        fail("checksums.sha256 differs from the ordered release manifest")

    panel_fields, panel = read_tsv(release / "analysis_panel.tsv")
    provenance_fields, provenance = read_tsv(release / "source_provenance.tsv")
    del panel_fields, provenance_fields
    if len(panel) != 45 or len(provenance) != 45:
        fail("release panel/source provenance is not the exact 45-trait family")
    if [row["trait_id"] for row in panel] != [row["trait_id"] for row in provenance]:
        fail("source provenance order differs from the release panel")
    for panel_row, source_row in zip(panel, provenance):
        if (
            panel_row["source_id"] != source_row["source_id"]
            or panel_row["dataset_version"] != source_row["dataset_version"]
            or panel_row["source_status"] != "SOURCE_VERIFIED"
            or source_row["source_status"] != "SOURCE_VERIFIED"
        ):
            fail(f"source provenance drift for {panel_row['trait_id']}")
        if (
            not source_row["source_page_url"].startswith("https://")
            or not source_row["download_url"].startswith("https://")
            or not source_row["archive_bytes"].isdigit()
            or int(source_row["archive_bytes"]) <= 0
            or not SHA256.fullmatch(source_row["archive_sha256"])
        ):
            fail(f"incomplete source URL/size/checksum provenance for {panel_row['trait_id']}")
    tool_fields, tools = read_tsv(release / "tool_versions.tsv")
    if tool_fields != ["component", "version_or_commit", "status", "purpose"] or not tools:
        fail("tool version registry has an unexpected schema")
    if any(row["status"] not in {"PINNED", "ACTIVE"} or not row["version_or_commit"] for row in tools):
        fail("tool version registry contains an unpinned component")
    tool_names = {row["component"] for row in tools}
    if not REQUIRED_TOOLS.issubset(tool_names):
        fail("tool version registry omits required atlas components")

    acceptance = json.loads((release / "acceptance_gates.json").read_text(encoding="utf-8"))
    non_release = [gate for gate in acceptance if gate.get("gate") != "atlas_v1_release"]
    if len(acceptance) != 23 or len(non_release) != 22 or any(gate.get("status") != "PASS" for gate in non_release):
        fail("frozen pre-release acceptance evidence is incomplete")
    if manifest.get("non_release_acceptance_gates") != 22:
        fail("release manifest has the wrong non-release gate count")
    if manifest.get("release_script_sha256") != downstream_contract.script_hashes(release, "release"):
        fail("release validator/builder script family differs from the manifest")

    validate_command(
        [sys.executable, str(release / "scripts/52_validate_integrated_atlas.py"),
         "--root", str(release), "--quiet"], root, "integrated atlas validation",
    )
    validate_command(
        [sys.executable, str(release / "scripts/53_validate_robustness.py"),
         "--root", str(release), "--quiet"], root, "robustness validation",
    )
    if not args.quiet:
        print(
            f"RELEASE_OK name=atlas-v1.0 files={len(entries)} "
            f"code_commit={code_commit} manifest_sha256={sha256(manifest_path)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
