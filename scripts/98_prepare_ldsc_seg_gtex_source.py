#!/usr/bin/env python3
"""Download and verify the frozen GTEx LDSC-SEG annotation subset."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.parse import quote


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(root: Path, value: str, label: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        fail(f"unsafe {label} path: {value}")
    return root / relative


def source_url(config: dict[str, object], repository_path: str) -> str:
    return (
        f"{config['mirror_base_url']}/{config['mirror_commit']}/"
        + quote(repository_path, safe="/")
        + "?download=true"
    )


def download(url: str, path: Path) -> None:
    if path.is_file() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    result = subprocess.run(
        [
            "curl", "--fail", "--location", "--silent", "--show-error",
            "--max-time", "60", "--output", str(temporary), url,
        ],
        capture_output=True, text=True, check=False,
    )
    if result.returncode:
        fail(f"source download failed: {url}: {(result.stdout + result.stderr).strip()}")
    if not temporary.is_file() or temporary.stat().st_size == 0:
        fail(f"source download is empty: {url}")
    temporary.replace(path)


def annotation_rows(path: Path, expected_header: str) -> int:
    rows = 0
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        if handle.readline().rstrip("\n") != expected_header:
            fail(f"LDSC-SEG annotation header differs from {expected_header}: {path}")
        for line in handle:
            if line.rstrip("\n") not in {"0", "1"}:
                fail(f"LDSC-SEG annotation contains a nonbinary value: {path}")
            rows += 1
    if not rows:
        fail(f"LDSC-SEG annotation is empty: {path}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--config", default="config/interpretation_ldsc_seg_gtex.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    config_path = root / args.config
    config = json.loads(config_path.read_text(encoding="utf-8"))
    local_root = safe_relative(root, config["local_root"], "LDSC-SEG local root")
    manifest_path = safe_relative(root, config["source_manifest_path"], "LDSC-SEG source manifest")
    if manifest_path.exists():
        fail("LDSC-SEG source manifest is immutable and already exists")
    selected = config["selected_tissues"]
    domains = {domain: 0 for domain in ("brain", "immune", "metabolic", "vascular")}
    labels: set[str] = set()
    indices: set[int] = set()
    for tissue in selected:
        domain, label, index = tissue["domain"], tissue["source_label"], int(tissue["source_index"])
        if domain not in domains or label in labels or index in indices or not 1 <= index <= 53:
            fail("LDSC-SEG tissue selection contains an invalid domain, label, or index")
        domains[domain] += 1
        labels.add(label)
        indices.add(index)
    if domains != {"brain": 4, "immune": 3, "metabolic": 4, "vascular": 5}:
        fail("LDSC-SEG tissue selection differs from the frozen four-domain family")

    repository_paths = [config["source_ldcts"]]
    for chromosome in range(1, 23):
        repository_paths.append(f"{config['source_prefix']}/GTEx.control.{chromosome}.annot.gz")
        for index in sorted(indices):
            repository_paths.append(f"{config['source_prefix']}/GTEx.{index}.{chromosome}.annot.gz")
    if len(repository_paths) != int(config["expected_source_files"]):
        fail("LDSC-SEG expected source-file count differs from the selected family")

    records: list[dict[str, object]] = []
    for repository_path in repository_paths:
        local_path = local_root / repository_path
        url = source_url(config, repository_path)
        download(url, local_path)
        records.append({
            "repository_path": repository_path,
            "local_path": str(local_path.relative_to(root)),
            "bytes": local_path.stat().st_size,
            "sha256": sha256(local_path),
            "url": url,
        })

    ldcts_path = local_root / config["source_ldcts"]
    with ldcts_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        ldcts = {row[0]: row[1] for row in reader if len(row) == 2}
    if len(ldcts) != 205:
        fail("Multi-tissue expression ldcts does not contain the exact 205-entry release family")
    for tissue in selected:
        expected = (
            f"Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.{tissue['source_index']}.,"
            "Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control."
        )
        if ldcts.get(tissue["source_label"]) != expected:
            fail(f"selected GTEx label/prefix differs from the source ldcts: {tissue['source_label']}")

    row_counts: dict[str, int] = {}
    for chromosome in range(1, 23):
        chromosome_paths = [
            local_root / f"{config['source_prefix']}/GTEx.control.{chromosome}.annot.gz",
            *[
                local_root / f"{config['source_prefix']}/GTEx.{index}.{chromosome}.annot.gz"
                for index in sorted(indices)
            ],
        ]
        counts = [
            annotation_rows(path, "All_Genes" if index == 0 else "ANNOT")
            for index, path in enumerate(chromosome_paths)
        ]
        if len(set(counts)) != 1:
            fail(f"GTEx target/control annotation row counts differ on chromosome {chromosome}")
        row_counts[str(chromosome)] = counts[0]
    manifest = {
        "schema_version": "sleep-atlas-ldsc-seg-source.1",
        "source_release": config["source_release"],
        "mirror_repository": config["mirror_repository"],
        "mirror_commit": config["mirror_commit"],
        "config_sha256": sha256(config_path),
        "selected_tissues": selected,
        "selected_tissues_by_domain": domains,
        "source_ldcts_entries": len(ldcts),
        "variant_rows_by_chromosome": row_counts,
        "source_file_count": len(records),
        "source_bytes": sum(int(row["bytes"]) for row in records),
        "files": records,
        "selection_rule": config["selection_rule"],
        "claim_limit": config["claim_limit"],
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = manifest_path.with_name(manifest_path.name + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(manifest_path)
    print(
        f"LDSC_SEG_GTEX_SOURCE_OK files={len(records)} bytes={manifest['source_bytes']} "
        f"tissues={len(selected)} out={manifest_path.relative_to(root)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
