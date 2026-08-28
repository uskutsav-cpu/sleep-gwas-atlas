#!/usr/bin/env python3
"""Plan, materialize, and verify the frozen LDSC-SEG reference family."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import quote


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha1(path: Path) -> str:
    digest = hashlib.sha1()
    digest.update(f"blob {path.stat().st_size}\0".encode("ascii"))
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def fetch_json(url: str) -> object:
    result = subprocess.run(
        ["curl", "--fail", "--location", "--silent", "--show-error", "--max-time", "60", url],
        capture_output=True, text=True, check=False,
    )
    if result.returncode:
        fail(f"reference metadata request failed: {url}: {(result.stdout + result.stderr).strip()}")
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        fail(f"reference metadata response is not JSON: {url}: {exc}")
    raise AssertionError


def family_records(config: dict[str, object], family: dict[str, object]) -> list[dict[str, object]]:
    url = str(config["mirror_api_template"]).format(
        commit=config["mirror_commit"], path=quote(str(family["remote_path"]), safe="/"),
    )
    payload = fetch_json(url)
    if not isinstance(payload, list):
        fail(f"reference metadata family is not a list: {family['family_id']}")
    pattern = re.compile(str(family["selection_regex"]))
    rows: list[dict[str, object]] = []
    for item in payload:
        if not isinstance(item, dict) or item.get("type") != "file":
            continue
        path = str(item.get("path", ""))
        if not pattern.fullmatch(path):
            continue
        lfs = item.get("lfs")
        if isinstance(lfs, dict):
            content_oid, oid_kind = str(lfs.get("oid", "")), "sha256"
        else:
            content_oid, oid_kind = str(item.get("oid", "")), "git_blob_sha1"
        row = {
            "content_oid": content_oid,
            "oid_kind": oid_kind,
            "path": path,
            "size": int(item.get("size", -1)),
        }
        if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", content_oid) or row["size"] <= 0:
            fail(f"reference metadata has an invalid object identity: {path}")
        rows.append(row)
    rows.sort(key=lambda row: str(row["path"]))
    if (
        len(rows) != int(family["expected_files"])
        or sum(int(row["size"]) for row in rows) != int(family["expected_bytes"])
        or canonical_sha256(rows) != family["metadata_sha256"]
    ):
        fail(f"remote reference family differs from its pre-result pin: {family['family_id']}")
    return rows


def remote_plan(config: dict[str, object]) -> dict[str, list[dict[str, object]]]:
    families = {
        str(family["family_id"]): family_records(config, family)
        for family in config["reference_families"]
    }
    observed = sum(int(row["size"]) for rows in families.values() for row in rows)
    observed += int(config["hm3_print_snps"]["expected_bytes"])
    if observed != int(config["expected_download_bytes"]):
        fail("reference family byte total differs from the frozen download plan")
    return families


def download_url(config: dict[str, object], repository_path: str) -> str:
    return str(config["mirror_download_template"]).format(
        commit=config["mirror_commit"], path=quote(repository_path, safe="/"),
    )


def validate_content(path: Path, size: int, oid: str, oid_kind: str) -> None:
    if not path.is_file() or path.stat().st_size != size:
        fail(f"reference file is absent or has the wrong byte count: {path}")
    observed = sha256(path) if oid_kind == "sha256" else git_blob_sha1(path)
    if observed != oid:
        fail(f"reference file differs from its immutable content identity: {path}")


def download_record(
    config: dict[str, object], root: Path, destination_root: Path, record: dict[str, object],
) -> Path:
    repository_path = str(record["path"])
    path = destination_root / repository_path
    if path.exists():
        validate_content(
            path, int(record["size"]), str(record["content_oid"]), str(record["oid_kind"]),
        )
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    if temporary.exists():
        fail(f"stale reference download exists; inspect it before retrying: {temporary}")
    result = subprocess.run(
        [
            "curl", "--fail", "--location", "--silent", "--show-error",
            "--retry", "3", "--output", str(temporary), download_url(config, repository_path),
        ],
        cwd=root, capture_output=True, text=True, check=False,
    )
    if result.returncode:
        fail(f"reference download failed: {repository_path}: {(result.stdout + result.stderr).strip()}")
    validate_content(
        temporary, int(record["size"]), str(record["content_oid"]), str(record["oid_kind"]),
    )
    temporary.replace(path)
    return path


def hm3_record(config: dict[str, object]) -> dict[str, object]:
    spec = config["hm3_print_snps"]
    return {
        "path": spec["remote_path"], "size": int(spec["expected_bytes"]),
        "content_oid": spec["content_sha256"], "oid_kind": "sha256",
    }


def runtime_versions(root: Path, config: dict[str, object]) -> dict[str, str]:
    ldsc = root / config["ldsc_script"]
    python = root / config["ldsc_python"]
    if not ldsc.is_file() or sha256(ldsc) != config["ldsc_script_sha256"]:
        fail("LDSC script differs from its exact pin")
    commit = subprocess.run(
        ["git", "-C", str(root / "ldsc"), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    ).stdout.strip()
    if commit != config["ldsc_commit"]:
        fail("LDSC checkout differs from its exact commit pin")
    code = (
        "import numpy,scipy,pandas,sys;"
        "print('.'.join(map(str,sys.version_info[:3])),numpy.__version__,scipy.__version__,pandas.__version__)"
    )
    result = subprocess.run([str(python), "-c", code], capture_output=True, text=True, check=False)
    if result.returncode:
        fail(f"LDSC runtime is unavailable: {(result.stdout + result.stderr).strip()}")
    values = result.stdout.split()
    expected = [
        config["python_version"], config["numpy_version"],
        config["scipy_version"], config["pandas_version"],
    ]
    if values != expected:
        fail(f"LDSC runtime versions differ from their pins: observed={values} expected={expected}")
    return dict(zip(("python", "numpy", "scipy", "pandas"), values))


def static_reference(
    root: Path, config_path: Path, config: dict[str, object],
    families: dict[str, list[dict[str, object]]],
) -> Path:
    static_root = root / config["static_reference_dir"]
    provenance_path = root / config["static_reference_provenance_path"]
    if provenance_path.exists():
        fail("static LDSC-SEG reference provenance is immutable and already exists")
    records = [
        *families["BASELINE_LD_V2_2"], *families["WEIGHTS_HM3_NO_MHC"], hm3_record(config),
    ]
    paths = [download_record(config, root, static_root, record) for record in records]
    provenance = {
        "schema_version": "sleep-atlas-ldsc-seg-static-reference.1",
        "source_release": config["source_release"],
        "mirror_commit": config["mirror_commit"],
        "reference_config_sha256": sha256(config_path),
        "files": [
            {
                "path": str(path.relative_to(root)), "bytes": path.stat().st_size,
                "sha256": sha256(path), "repository_path": record["path"],
                "content_oid": record["content_oid"], "oid_kind": record["oid_kind"],
            }
            for path, record in zip(paths, records)
        ],
        "file_count": len(paths),
        "total_bytes": sum(path.stat().st_size for path in paths),
    }
    atomic_json(provenance_path, provenance)
    print(
        f"LDSC_SEG_STATIC_REFERENCE_OK files={len(paths)} bytes={provenance['total_bytes']} "
        f"out={provenance_path.relative_to(root)}"
    )
    return provenance_path


def validate_static(
    root: Path, config_path: Path, config: dict[str, object],
    families: dict[str, list[dict[str, object]]],
) -> dict[str, object]:
    path = root / config["static_reference_provenance_path"]
    if not path.is_file():
        fail("static LDSC-SEG reference is absent; run --download-static after reviewing the plan")
    provenance = json.loads(path.read_text(encoding="utf-8"))
    expected_records = [
        *families["BASELINE_LD_V2_2"], *families["WEIGHTS_HM3_NO_MHC"], hm3_record(config),
    ]
    expected_by_repository_path = {
        str(record["path"]): record for record in expected_records
    }
    observed_records = provenance.get("files", [])
    if not isinstance(observed_records, list):
        fail("static LDSC-SEG reference provenance has an invalid file family")
    observed_by_repository_path = {
        str(record.get("repository_path", "")): record for record in observed_records
        if isinstance(record, dict)
    }
    if (
        provenance.get("reference_config_sha256") != sha256(config_path)
        or provenance.get("mirror_commit") != config["mirror_commit"]
        or provenance.get("file_count") != len(expected_records)
        or len(observed_by_repository_path) != len(observed_records)
        or set(observed_by_repository_path) != set(expected_by_repository_path)
    ):
        fail("static LDSC-SEG reference provenance differs from the frozen plan")
    total = 0
    static_root = root / config["static_reference_dir"]
    for repository_path, expected in expected_by_repository_path.items():
        record = observed_by_repository_path[repository_path]
        relative = Path(str(record.get("path", "")))
        expected_path = static_root / repository_path
        if (
            relative.is_absolute() or ".." in relative.parts
            or root / relative != expected_path
            or record.get("bytes") != expected["size"]
            or record.get("content_oid") != expected["content_oid"]
            or record.get("oid_kind") != expected["oid_kind"]
        ):
            fail("static LDSC-SEG reference provenance contains an unsafe path")
        file_path = root / relative
        if (
            not file_path.is_file() or file_path.stat().st_size != record.get("bytes")
            or sha256(file_path) != record.get("sha256")
        ):
            fail(f"static LDSC-SEG reference file differs from provenance: {relative}")
        total += file_path.stat().st_size
    if total != provenance.get("total_bytes"):
        fail("static LDSC-SEG reference byte total differs from provenance")
    return provenance


def read_source_manifest(
    root: Path, policy: dict[str, object], selection_path: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    spec = policy["ldsc_seg_gtex"]
    manifest_path = root / spec["source_manifest"]
    if (
        not manifest_path.is_file() or manifest_path.stat().st_size != spec["source_manifest_bytes"]
        or sha256(manifest_path) != spec["source_manifest_sha256"]
    ):
        fail("LDSC-SEG annotation source manifest differs from its policy pin")
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("config_sha256") != sha256(selection_path):
        fail("LDSC-SEG annotation manifest differs from its selection config")
    return selection, manifest


def chromosome_genotype_records(
    records: list[dict[str, object]], chromosome: int,
) -> list[dict[str, object]]:
    suffixes = {f".{chromosome}.bed", f".{chromosome}.bim", f".{chromosome}.fam"}
    selected = [row for row in records if any(str(row["path"]).endswith(suffix) for suffix in suffixes)]
    if len(selected) != 3:
        fail(f"genotype metadata is incomplete for chromosome {chromosome}")
    return selected


def derive_reference(
    root: Path, policy_path: Path, policy: dict[str, object], selection_path: Path,
    reference_path: Path, config: dict[str, object], families: dict[str, list[dict[str, object]]],
    *, acknowledge_large: bool, acknowledge_deletion: bool,
) -> None:
    if not acknowledge_large or not acknowledge_deletion:
        fail(
            "the full reference plan transfers 1.88 GB, including 1.59 GB of streamed "
            "genotypes, and deletes only the 66 checksum-verified temporary chromosome PLINK "
            "files; rerun with both acknowledgement flags after approval"
        )
    free = shutil.disk_usage(root).free
    if free < int(config["minimum_free_bytes_before_streaming"]):
        fail(
            f"streamed derivation requires at least {config['minimum_free_bytes_before_streaming']} "
            f"free bytes; observed {free}"
        )
    static = validate_static(root, reference_path, config, families)
    selection, source_manifest = read_source_manifest(root, policy, selection_path)
    runtime = runtime_versions(root, config)
    output_root = root / config["ldscore_cache_dir"]
    provenance_path = root / config["reference_provenance_path"]
    ldcts_path = root / config["selected_ldcts_path"]
    if provenance_path.exists() or ldcts_path.exists():
        fail("derived LDSC-SEG reference is immutable and already exists")
    existing = [path for path in output_root.glob("GTEx.*") if path.is_file()]
    if existing:
        fail(f"derived LDSC-SEG output exists without immutable provenance: {existing[0]}")
    output_root.mkdir(parents=True, exist_ok=True)
    temporary_root = root / config["temporary_genotype_dir"]
    if temporary_root.exists() and any(temporary_root.rglob("*")):
        fail("LDSC-SEG temporary genotype directory is not empty; inspect it before retrying")
    source_root = root / selection["local_root"]
    print_snps = root / config["static_reference_dir"] / config["hm3_print_snps"]["remote_path"]
    output_records: list[dict[str, object]] = []
    deletion_records: list[dict[str, object]] = []
    selected_indices = [int(row["source_index"]) for row in selection["selected_tissues"]]
    for chromosome in range(1, 23):
        build_dir = output_root / f".build_chr{chromosome}"
        genotype_dir = temporary_root / f"chr{chromosome}"
        if build_dir.exists() or genotype_dir.exists():
            fail(f"stale LDSC-SEG chromosome workspace exists: chromosome {chromosome}")
        build_dir.mkdir(parents=True)
        genotype_records = chromosome_genotype_records(families["GENOTYPE"], chromosome)
        genotype_paths = [
            download_record(config, root, genotype_dir, record) for record in genotype_records
        ]
        genotype_prefix = genotype_dir / "1000G_EUR_Phase3_plink" / f"1000G.EUR.QC.{chromosome}"
        annotations = [("control", "All_Genes"), *[(str(index), "ANNOT") for index in selected_indices]]
        chromosome_outputs: list[tuple[Path, Path]] = []
        for label, expected_header in annotations:
            annotation = source_root / selection["source_prefix"] / f"GTEx.{label}.{chromosome}.annot.gz"
            with gzip.open(annotation, "rt", encoding="ascii") as handle:
                if handle.readline().rstrip("\n") != expected_header:
                    fail(f"LDSC-SEG source annotation header differs on chromosome {chromosome}: {label}")
            annotation_prefix = str(annotation)[:-len(".annot.gz")]
            temporary_prefix = build_dir / f"GTEx.{label}.{chromosome}"
            command = [
                str(root / config["ldsc_python"]), str(root / config["ldsc_script"]),
                "--l2", "--bfile", str(genotype_prefix),
                "--ld-wind-cm", str(config["ld_window_cm"]),
                "--annot", annotation_prefix, "--thin-annot",
                "--print-snps", str(print_snps), "--out", str(temporary_prefix),
            ]
            result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
            if result.returncode:
                fail(
                    f"LDSC-SEG LD-score derivation failed for chromosome {chromosome}/{label}: "
                    f"{(result.stdout + result.stderr).strip()}"
                )
            for suffix in (".l2.ldscore.gz", ".l2.M", ".l2.M_5_50", ".log"):
                temporary_output = Path(str(temporary_prefix) + suffix)
                final_output = output_root / (temporary_prefix.name + suffix)
                if not temporary_output.is_file() or temporary_output.stat().st_size == 0:
                    fail(f"LDSC-SEG derivation omitted output: {temporary_output}")
                chromosome_outputs.append((temporary_output, final_output))
        for temporary_output, final_output in chromosome_outputs:
            temporary_output.replace(final_output)
            output_records.append({
                "path": str(final_output.relative_to(root)),
                "bytes": final_output.stat().st_size,
                "sha256": sha256(final_output),
            })
        build_dir.rmdir()
        for path, record in zip(genotype_paths, genotype_records):
            deletion_records.append({
                "path": str(path.relative_to(root)), "bytes": path.stat().st_size,
                "content_oid": record["content_oid"], "oid_kind": record["oid_kind"],
            })
            path.unlink()
        for directory in sorted(
            [path for path in genotype_dir.rglob("*") if path.is_dir()], reverse=True,
        ):
            directory.rmdir()
        genotype_dir.rmdir()
        print(f"LDSC_SEG_CHROMOSOME_OK chromosome={chromosome} temporary_files_deleted=3")
    if temporary_root.exists():
        temporary_root.rmdir()
    ldcts_lines = []
    for tissue in selection["selected_tissues"]:
        target = f"{config['ldscore_cache_dir']}/GTEx.{tissue['source_index']}."
        control = f"{config['ldscore_cache_dir']}/GTEx.control."
        ldcts_lines.append(f"{tissue['source_label']}\t{target},{control}")
    ldcts_payload = "\n".join(ldcts_lines) + "\n"
    ldcts_path.parent.mkdir(parents=True, exist_ok=True)
    ldcts_temporary = ldcts_path.with_name(ldcts_path.name + ".tmp")
    ldcts_temporary.write_text(ldcts_payload, encoding="utf-8")
    ldcts_temporary.replace(ldcts_path)
    provenance = {
        "schema_version": "sleep-atlas-ldsc-seg-derived-reference.1",
        "source_release": config["source_release"],
        "mirror_commit": config["mirror_commit"],
        "reference_config_sha256": sha256(reference_path),
        "selection_config_sha256": sha256(selection_path),
        "annotation_source_manifest_sha256": policy["ldsc_seg_gtex"]["source_manifest_sha256"],
        "static_reference_provenance_sha256": sha256(root / config["static_reference_provenance_path"]),
        "runtime": runtime,
        "ldsc_script_sha256": sha256(root / config["ldsc_script"]),
        "ld_window_cm": config["ld_window_cm"],
        "thin_annotation": config["thin_annotation"],
        "print_snps_sha256": sha256(print_snps),
        "selected_tissues": selection["selected_tissues"],
        "selected_ldcts_path": str(ldcts_path.relative_to(root)),
        "selected_ldcts_sha256": sha256(ldcts_path),
        "derived_file_count": len(output_records),
        "derived_files": sorted(output_records, key=lambda row: str(row["path"])),
        "temporary_deletions": deletion_records,
        "temporary_deleted_bytes": sum(int(row["bytes"]) for row in deletion_records),
        "static_reference_bytes": static["total_bytes"],
        "analysis_rule": config["analysis_rule"],
        "claim_limit": config["claim_limit"],
        "script_sha256": sha256(Path(__file__)),
    }
    genotype_family = next(
        row for row in config["reference_families"] if row["family_id"] == "GENOTYPE"
    )
    if (
        len(output_records) != 22 * 17 * 4
        or len(deletion_records) != int(genotype_family["expected_files"])
        or provenance["temporary_deleted_bytes"] != int(genotype_family["expected_bytes"])
    ):
        fail("derived LDSC-SEG output or temporary deletion family is incomplete")
    atomic_json(provenance_path, provenance)
    print(
        f"LDSC_SEG_REFERENCE_OK files={len(output_records)} "
        f"temporary_deleted_bytes={provenance['temporary_deleted_bytes']} "
        f"out={provenance_path.relative_to(root)}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--selection", default="config/interpretation_ldsc_seg_gtex.json")
    parser.add_argument("--reference", default="config/interpretation_ldsc_seg_reference.json")
    parser.add_argument("--download-static", action="store_true")
    parser.add_argument("--derive", action="store_true")
    parser.add_argument("--acknowledge-large-download", action="store_true")
    parser.add_argument("--acknowledge-temporary-deletion", action="store_true")
    args = parser.parse_args()
    if args.download_static and args.derive:
        fail("choose one materialization action at a time")
    root = Path(args.root).resolve()
    policy_path, selection_path, reference_path = (
        root / args.policy, root / args.selection, root / args.reference,
    )
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    config = json.loads(reference_path.read_text(encoding="utf-8"))
    spec = policy["ldsc_seg_gtex"]
    if (
        spec["reference_config"] != str(reference_path.relative_to(root))
        or spec["reference_config_sha256"] != sha256(reference_path)
        or spec["selection_config"] != str(selection_path.relative_to(root))
        or spec["selection_config_sha256"] != sha256(selection_path)
        or config["mirror_commit"] != json.loads(selection_path.read_text(encoding="utf-8"))["mirror_commit"]
        or config["ldscore_cache_dir"] != spec["ldscore_cache_dir"]
        or config["reference_provenance_path"] != spec["ldscore_cache_provenance_path"]
    ):
        fail("LDSC-SEG source, reference, and interpretation policies are not aligned")
    families = remote_plan(config)
    free = shutil.disk_usage(root).free
    genotype = next(row for row in config["reference_families"] if row["family_id"] == "GENOTYPE")
    persistent = int(config["expected_download_bytes"]) - int(genotype["expected_bytes"])
    print(
        f"LDSC_SEG_REFERENCE_PLAN transfer_bytes={config['expected_download_bytes']} "
        f"persistent_source_bytes={persistent} temporary_genotype_bytes={genotype['expected_bytes']} "
        f"free_bytes={free} delete_exactly=66_temporary_plink_files_after_each_chromosome"
    )
    if args.download_static:
        static_reference(root, reference_path, config, families)
    elif args.derive:
        derive_reference(
            root, policy_path, policy, selection_path, reference_path, config, families,
            acknowledge_large=args.acknowledge_large_download,
            acknowledge_deletion=args.acknowledge_temporary_deletion,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
