#!/usr/bin/env python3
"""Stream full harmonized GWAS into the official MiXeR input schema."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT_SPEC = importlib.util.spec_from_file_location(
    "mixer_preflight", ROOT / "scripts/35_mixer_preflight.py"
)
if PREFLIGHT_SPEC is None or PREFLIGHT_SPEC.loader is None:
    raise RuntimeError("could not load MiXeR preflight helpers")
preflight_module = importlib.util.module_from_spec(PREFLIGHT_SPEC)
PREFLIGHT_SPEC.loader.exec_module(preflight_module)
OUTPUT_COLUMNS = ["SNP", "CHR", "BP", "A1", "A2", "N", "Z"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def convert(source: Path, destination: Path) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    rows = 0
    try:
        with gzip.open(source, "rt", encoding="utf-8", newline="") as incoming, temporary.open("wb") as raw_out:
            reader = csv.DictReader(incoming, delimiter="\t")
            required = {"SNP", "CHR", "BP", "A1", "A2", "N", "BETA", "SE"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise SystemExit(f"ERROR: {source} lacks columns: {sorted(missing)}")
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_out, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as outgoing:
                    writer = csv.DictWriter(outgoing, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    for row in reader:
                        try:
                            beta = float(row["BETA"])
                            se = float(row["SE"])
                            sample_size = float(row["N"])
                        except ValueError as exc:
                            raise SystemExit(f"ERROR: non-numeric BETA/SE/N in {source}") from exc
                        if not (
                            math.isfinite(beta) and math.isfinite(se) and se > 0
                            and math.isfinite(sample_size) and sample_size > 0
                        ):
                            raise SystemExit(f"ERROR: invalid BETA/SE/N in {source}")
                        z = beta / se
                        writer.writerow({
                            "SNP": row["SNP"], "CHR": row["CHR"], "BP": row["BP"],
                            "A1": row["A1"], "A2": row["A2"], "N": row["N"], "Z": f"{z:.15g}",
                        })
                        rows += 1
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return {
        "source_bytes": source.stat().st_size,
        "source_sha256": sha256(source),
        "output_bytes": destination.stat().st_size,
        "output_sha256": sha256(destination),
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--out-dir", default="data/mixer")
    parser.add_argument("--manifest", default="results/tables/mixer_input_manifest.tsv")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    panel = preflight_module.read_tsv(root / "config/analysis_panel.tsv")
    planned_prefilters = {
        row["trait_id"]
        for row in preflight_module.read_tsv(root / "config/hm3_prefilter_plans.tsv")
    }
    sources = []
    blockers = []
    for row in panel:
        data, qc, prefilter = preflight_module.choose_harmonized(root, row["trait_id"])
        if not prefilter and row["trait_id"] not in planned_prefilters:
            prefilter = "not supplied"
        if not data.is_file() or not qc.is_file() or prefilter != "not supplied":
            blockers.append(row["trait_id"])
        sources.append((row["trait_id"], data))
    if blockers:
        raise SystemExit(
            "ERROR: full non-HapMap3 harmonized inputs are required for all traits; blocked: "
            + ", ".join(blockers)
        )

    source_bytes = sum(path.stat().st_size for _, path in sources)
    required_free = math.ceil(source_bytes * 1.5) + 2 * 1024**3
    free = shutil.disk_usage(root).free
    print(f"Full harmonized source bytes: {source_bytes}")
    print(f"Required free-space preflight: {required_free}; observed: {free}")
    if not args.materialize:
        print("No conversion requested. Re-run with --materialize on a preflight-ready x86 analysis host.")
        return 0
    if free < required_free:
        raise SystemExit("ERROR: insufficient free space to materialize MiXeR inputs safely")

    provenance = []
    for index, (trait, source) in enumerate(sources, start=1):
        destination = root / args.out_dir / f"{trait}.sumstats.gz"
        metrics = convert(source, destination)
        provenance.append({
            "trait_id": trait,
            "filename": str(destination.relative_to(root)),
            **metrics,
            "input_scope": "FULL_POST_QC_AUTOSOMAL",
            "status": "VALIDATED",
        })
        print(f"[{index}/45] {trait}: {metrics['rows']} rows")
    manifest = root / args.manifest
    manifest.parent.mkdir(parents=True, exist_ok=True)
    temporary = manifest.with_suffix(manifest.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        fields = list(provenance[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(provenance)
    temporary.replace(manifest)
    print("Published 45 checksum-bound full-summary-statistics MiXeR inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
