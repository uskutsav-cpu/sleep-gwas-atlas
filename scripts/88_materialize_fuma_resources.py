#!/usr/bin/env python3
"""Materialize checksum-pinned FUMA matrices or the disk-guarded MAGMA LD reference."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import tarfile
import zipfile
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_bundle(root: Path, policy: dict[str, object]) -> tuple[dict[str, object], Path, dict[str, Path]]:
    spec = policy["fuma_scrna"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("FUMA scRNA component manifest is absent or differs from the policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components: dict[str, Path] = {}
    for component in manifest.get("components", []):
        relative = Path(component["path"])
        if relative.is_absolute() or ".." in relative.parts:
            fail(f"unsafe FUMA component path: {component['component_id']}")
        path = root / relative
        if (
            not path.is_file() or path.stat().st_size != component["bytes"]
            or sha256(path) != component["sha256"]
        ):
            fail(f"FUMA component differs from its exact pin: {component['component_id']}")
        components[component["component_id"]] = path
    return manifest, manifest_path, components


def materialize_matrices(
    archive_path: Path, prefix: str, datasets: list[dict[str, object]], output_dir: Path,
) -> list[dict[str, object]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    with tarfile.open(archive_path, "r:gz") as archive:
        for dataset in datasets:
            identity = str(dataset["dataset_id"])
            member_name = prefix + str(dataset["member"])
            try:
                member = archive.getmember(member_name)
            except KeyError:
                fail(f"selected FUMA matrix is absent from commit archive: {identity}")
            handle = archive.extractfile(member)
            if handle is None:
                fail(f"selected FUMA matrix cannot be read: {identity}")
            raw = handle.read()
            if (
                len(raw) != dataset["compressed_bytes"]
                or hashlib.sha256(raw).hexdigest() != dataset["compressed_sha256"]
            ):
                fail(f"selected FUMA matrix differs from its compressed pin: {identity}")
            try:
                payload = gzip.decompress(raw)
            except gzip.BadGzipFile as exc:
                fail(f"selected FUMA matrix is not valid gzip: {identity}")
                raise AssertionError from exc
            if (
                len(payload) != dataset["text_bytes"]
                or hashlib.sha256(payload).hexdigest() != dataset["text_sha256"]
            ):
                fail(f"selected FUMA matrix differs from its text pin: {identity}")
            output = output_dir / f"{identity}.txt"
            if output.exists():
                if output.stat().st_size != len(payload) or sha256(output) != dataset["text_sha256"]:
                    fail(f"existing materialized FUMA matrix differs from its pin: {identity}")
            else:
                temporary = output.with_suffix(output.suffix + ".tmp")
                temporary.write_bytes(payload)
                if temporary.stat().st_size != len(payload) or sha256(temporary) != dataset["text_sha256"]:
                    fail(f"materialized FUMA matrix failed checksum validation: {identity}")
                temporary.replace(output)
            records.append({
                "dataset_id": identity, "domain": dataset["domain"], "species": dataset["species"],
                "tissue": dataset["tissue"], "path": str(output), "bytes": output.stat().st_size,
                "sha256": sha256(output), "gene_rows": dataset["gene_rows"],
                "cell_type_count": dataset["cell_type_count"],
            })
    return records


def materialize_reference(
    archive_path: Path, members: list[dict[str, object]], output_dir: Path,
    minimum_free_bytes: int, acknowledged: bool,
) -> list[dict[str, object]]:
    if not acknowledged:
        fail(
            "MAGMA EUR reference extraction expands to 3,600,697,827 bytes; warn the user first "
            "and rerun with --acknowledge-large-extract after approval"
        )
    expected = {str(row["name"]): row for row in members}
    present_bytes = 0
    for name, row in expected.items():
        path = output_dir / name
        if path.is_file() and path.stat().st_size == row["bytes"] and sha256(path) == row["sha256"]:
            present_bytes += int(row["bytes"])
        elif path.exists():
            fail(f"existing MAGMA reference member differs from its pin: {name}")
    required_missing = sum(int(row["bytes"]) for row in expected.values()) - present_bytes
    safety_margin = minimum_free_bytes - sum(int(row["bytes"]) for row in expected.values())
    free_bytes = shutil.disk_usage(output_dir.parent).free
    if free_bytes < required_missing + safety_margin:
        fail(
            f"insufficient disk for guarded MAGMA reference extraction: free={free_bytes} "
            f"required={required_missing + safety_margin}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    try:
        with zipfile.ZipFile(archive_path) as archive:
            infos = {info.filename: info for info in archive.infolist() if not info.is_dir()}
            if set(infos) != set(expected):
                fail("MAGMA EUR reference archive member family differs from its pin")
            for name, row in expected.items():
                output = output_dir / name
                if not output.exists():
                    temporary = output.with_suffix(output.suffix + ".tmp")
                    digest = hashlib.sha256()
                    with archive.open(infos[name]) as source, temporary.open("wb") as target:
                        for block in iter(lambda: source.read(1024 * 1024), b""):
                            target.write(block)
                            digest.update(block)
                    if temporary.stat().st_size != row["bytes"] or digest.hexdigest() != row["sha256"]:
                        fail(f"extracted MAGMA reference member differs from its pin: {name}")
                    temporary.replace(output)
                records.append({"name": name, "path": str(output), "bytes": output.stat().st_size, "sha256": sha256(output)})
    except (OSError, zipfile.BadZipFile, RuntimeError) as exc:
        fail(f"MAGMA EUR reference extraction failed: {exc}")
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--matrices", action="store_true")
    mode.add_argument("--reference", action="store_true")
    parser.add_argument("--acknowledge-large-extract", action="store_true")
    parser.add_argument("--provenance-out")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["fuma_scrna"]
    manifest, manifest_path, components = load_bundle(root, policy)
    if args.matrices:
        records = materialize_matrices(
            components["FUMA_SCRNA_COMMIT_ARCHIVE"], spec["archive_member_prefix"],
            manifest["datasets"], root / spec["matrix_cache_dir"],
        )
        default_provenance = spec["matrix_cache_provenance_path"]
        mode_name = "MATRICES"
    else:
        records = materialize_reference(
            components["MAGMA_1000G_PHASE3_EUR_ARCHIVE"], manifest["reference_members"],
            root / spec["reference_dir"], int(spec["minimum_free_bytes_before_reference_extract"]),
            args.acknowledge_large_extract,
        )
        default_provenance = spec["reference_extract_provenance_path"]
        mode_name = "REFERENCE"
    provenance_path = root / (args.provenance_out or default_provenance)
    atomic_json(provenance_path, {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "mode": mode_name, "policy_sha256": sha256(policy_path),
        "component_manifest_sha256": sha256(manifest_path), "source_release": manifest["release"],
        "source_component_sha256": {
            identity: sha256(path) for identity, path in sorted(components.items())
        },
        "records": records,
    })
    print(f"FUMA_RESOURCE_{mode_name}_OK records={len(records)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
