#!/usr/bin/env python3
"""Create the checksum-bound Ensembl-v92 ±1 kb MAGMA SNP-to-gene annotation."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--out-prefix")
    parser.add_argument("--provenance-out")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["fuma_scrna"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("FUMA component manifest differs from its policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    component = {row["component_id"]: row for row in manifest["components"]}
    binary = root / spec["magma_binary_path"]
    binary_pin = component["MAGMA_V1_10_ARM64_BINARY"]
    if (
        not binary.is_file() or binary.stat().st_size != binary_pin["bytes"]
        or sha256(binary) != binary_pin["sha256"]
    ):
        fail("MAGMA v1.10 binary differs from its exact ARM64 source-build pin")
    version = subprocess.run([binary, "--version"], capture_output=True, text=True, check=False)
    if version.returncode or (version.stdout + version.stderr).strip() != "MAGMA version: v1.10 (custom)":
        fail("MAGMA binary does not report v1.10")
    gene_location = root / spec["gene_location_path"]
    gene_pin = component["FUMA_ENSEMBL_V92_GENE_BOUNDARIES"]
    if (
        not gene_location.is_file() or gene_location.stat().st_size != gene_pin["bytes"]
        or sha256(gene_location) != gene_pin["sha256"]
    ):
        fail("Ensembl-v92 gene boundaries differ from their exact pin")
    reference_dir = root / spec["reference_dir"]
    reference_pins = {row["name"]: row for row in manifest["reference_members"]}
    reference_paths: dict[str, Path] = {}
    for name, pin in reference_pins.items():
        path = reference_dir / name
        if not path.is_file() or path.stat().st_size != pin["bytes"] or sha256(path) != pin["sha256"]:
            fail(
                f"extracted MAGMA reference differs from its pin: {name}; run the disk-guarded "
                "FUMA reference materializer after confirming sufficient free space"
            )
        reference_paths[name] = path
    prefix = root / (args.out_prefix or spec["gene_annotation_prefix"])
    output, log = Path(str(prefix) + ".genes.annot"), Path(str(prefix) + ".log")
    provenance_path = root / (args.provenance_out or spec["gene_annotation_provenance_path"])
    if output.exists() or log.exists() or provenance_path.exists():
        fail("MAGMA annotation output already exists; refusing overwrite")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    temporary_prefix = prefix.with_name(prefix.name + ".tmp")
    temporary_output = Path(str(temporary_prefix) + ".genes.annot")
    temporary_log = Path(str(temporary_prefix) + ".log")
    if temporary_output.exists() or temporary_log.exists():
        fail("stale MAGMA annotation temporary output exists; inspect it before retrying")
    upstream, downstream = map(int, spec["gene_window_kb"])
    command = [
        str(binary), "--annotate", f"window={upstream},{downstream}",
        "--snp-loc", str(reference_paths["g1000_eur.bim"]),
        "--gene-loc", str(gene_location), "--out", str(temporary_prefix),
    ]
    result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
    if result.returncode or not temporary_output.is_file() or not temporary_log.is_file():
        fail(f"MAGMA annotation failed: {(result.stdout + result.stderr).strip()}")
    annotation_rows = sum(bool(line.strip()) for line in temporary_output.open("rb"))
    if not int(spec["minimum_genes_in_gene_property_model"]) <= annotation_rows <= manifest["gene_count"]:
        fail(f"MAGMA annotation produced an implausible gene family: {annotation_rows}")
    temporary_output.replace(output)
    temporary_log.replace(log)
    atomic_json(provenance_path, {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "component_manifest_sha256": sha256(manifest_path),
        "source_release": manifest["release"], "magma_version": spec["magma_version"],
        "magma_binary_sha256": sha256(binary), "gene_location_sha256": sha256(gene_location),
        "reference_member_sha256": {
            name: sha256(path) for name, path in sorted(reference_paths.items())
        },
        "gene_window_kb": spec["gene_window_kb"], "command": command,
        "annotation_rows": annotation_rows, "annotation_path": str(output),
        "annotation_sha256": sha256(output), "log_path": str(log), "log_sha256": sha256(log),
    })
    print(f"MAGMA_ANNOTATION_OK genes={annotation_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
