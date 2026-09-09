#!/usr/bin/env python3
"""Materialize one trait's P values on the fixed CATlas variant universe."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

from catlas_policy import reference_policy_sha256, trait_policy_sha256


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


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def atomic_gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    temporary.replace(path)


def load_universe(path: Path, expected_fields: list[str]) -> tuple[list[str], dict[str, tuple[int, int]]]:
    order: list[str] = []
    coordinates: dict[str, tuple[int, int]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected_fields:
            fail("CATlas fixed-universe cache header differs from policy")
        for row in reader:
            snp = row["SNP"]
            try:
                coordinate = (int(row["CHR37"]), int(row["BP37"]))
            except ValueError:
                fail(f"CATlas fixed-universe cache contains an invalid coordinate: {snp}")
            if snp in coordinates:
                fail(f"CATlas fixed-universe cache contains a duplicate variant: {snp}")
            coordinates[snp] = coordinate
            order.append(snp)
    if not order:
        fail("CATlas fixed-universe cache is empty")
    return order, coordinates


def materialize_trait(
    gwas_path: Path, universe_order: list[str], universe_coordinates: dict[str, tuple[int, int]],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    selected = set(universe_order)
    values: dict[str, str] = {}
    rows_read = 0
    coordinate_mismatches = 0
    with gzip.open(gwas_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "CHR", "BP", "P"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            fail("harmonized GWAS lacks SNP, CHR, BP, or P")
        for row in reader:
            rows_read += 1
            snp = row["SNP"]
            if snp not in selected:
                continue
            if snp in values:
                fail(f"harmonized GWAS contains a duplicate fixed-universe variant: {snp}")
            try:
                chromosome, position = int(row["CHR"]), int(row["BP"])
                p_value = float(row["P"])
            except ValueError:
                fail(f"harmonized GWAS contains invalid CATlas fields for {snp}")
            if (chromosome, position) != universe_coordinates[snp]:
                coordinate_mismatches += 1
                continue
            if not math.isfinite(p_value) or not 0 < p_value <= 1:
                fail(f"harmonized GWAS contains an invalid P value for {snp}")
            values[snp] = format(p_value, ".12g")
    rows = [{"SNP": snp, "P": values[snp]} for snp in universe_order if snp in values]
    return rows, {
        "harmonized_rows": rows_read,
        "fixed_universe_variants": len(universe_order),
        "matched_variants": len(rows),
        "missing_variants": len(universe_order) - len(rows),
        "coordinate_mismatch_variants": coordinate_mismatches,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trait_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["catlas_adult_v4"]
    with (root / "config/analysis_panel.tsv").open(encoding="utf-8", newline="") as handle:
        traits = {row["trait_id"] for row in csv.DictReader(handle, delimiter="\t")}
    if args.trait_id not in traits:
        fail(f"trait is outside the locked panel: {args.trait_id}")

    universe_path = safe_relative(root, spec["variant_cache_path"], "CATlas fixed-universe cache")
    universe_provenance_path = safe_relative(
        root, spec["variant_cache_provenance_path"], "CATlas fixed-universe provenance",
    )
    if not universe_path.is_file() or not universe_provenance_path.is_file():
        fail("CATlas fixed-universe cache and provenance are incomplete")
    universe_provenance = json.loads(universe_provenance_path.read_text(encoding="utf-8"))
    if (
        universe_provenance.get("reference_policy_sha256") != reference_policy_sha256(policy)
        or universe_provenance.get("component_manifest_sha256") != spec["component_manifest_sha256"]
        or universe_provenance.get("cache_sha256") != sha256(universe_path)
    ):
        fail("CATlas fixed-universe cache differs from current policy-bound provenance")

    gwas_path = safe_relative(
        root, spec["gwas_path_template"].format(trait_id=args.trait_id), "harmonized GWAS",
    )
    if not gwas_path.is_file() or gwas_path.stat().st_size == 0:
        fail(f"harmonized GWAS is absent: {args.trait_id}")
    out_path = safe_relative(
        root, spec["trait_cache_path_template"].format(trait_id=args.trait_id), "CATlas trait cache",
    )
    provenance_path = safe_relative(
        root, spec["trait_cache_provenance_path_template"].format(trait_id=args.trait_id),
        "CATlas trait provenance",
    )
    if out_path.exists() or provenance_path.exists():
        fail("CATlas trait cache is immutable and already exists")

    universe_order, universe_coordinates = load_universe(universe_path, spec["variant_cache_fields"])
    rows, counts = materialize_trait(gwas_path, universe_order, universe_coordinates)
    coverage = counts["matched_variants"] / counts["fixed_universe_variants"]
    if coverage < float(spec["minimum_trait_variant_coverage"]):
        fail(
            f"trait covers only {coverage:.3%} of the fixed CATlas universe; "
            f"minimum is {float(spec['minimum_trait_variant_coverage']):.3%}"
        )
    atomic_gzip_tsv(out_path, spec["trait_cache_fields"], rows)
    provenance = {
        "schema_version": "sleep-atlas-catlas-trait.2",
        "trait_id": args.trait_id,
        "trait_policy_sha256": trait_policy_sha256(policy),
        "fixed_universe_provenance_sha256": sha256(universe_provenance_path),
        "fixed_universe_sha256": sha256(universe_path),
        "harmonized_gwas_path": str(gwas_path.relative_to(root)),
        "harmonized_gwas_bytes": gwas_path.stat().st_size,
        "harmonized_gwas_sha256": sha256(gwas_path),
        "counts": counts,
        "coverage_fraction": coverage,
        "cache_path": str(out_path.relative_to(root)),
        "cache_bytes": out_path.stat().st_size,
        "cache_sha256": sha256(out_path),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    atomic_json(provenance_path, provenance)
    print(
        f"CATLAS_TRAIT_OK trait={args.trait_id} matched={counts['matched_variants']} "
        f"coverage={coverage:.6f} out={out_path.relative_to(root)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
