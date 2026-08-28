#!/usr/bin/env python3
"""Fail closed unless the full locked-panel GenomicSEM covariance export is valid."""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def expected_elements(traits: list[str]) -> list[tuple[str, str, str]]:
    return [
        (f"{traits[i]}__{traits[j]}", traits[i], traits[j])
        for i in range(len(traits))
        for j in range(i, len(traits))
    ]


def validate_square_matrix(
    path: Path, traits: list[str], first_column: str, errors: list[str]
) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        errors.append(f"missing or empty matrix: {path}")
        return
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        require(header == [first_column, *traits], f"wrong header/order: {path}", errors)
        rows = list(reader)
    require(len(rows) == len(traits), f"wrong row count: {path}", errors)
    if len(rows) != len(traits):
        return
    require([row[0] for row in rows] == traits, f"wrong row trait order: {path}", errors)
    values: list[list[float]] = []
    for index, row in enumerate(rows):
        require(len(row) == len(traits) + 1, f"wrong width at row {index + 2}: {path}", errors)
        if len(row) != len(traits) + 1:
            return
        require(all(finite(value) for value in row[1:]), f"non-finite value at row {index + 2}: {path}", errors)
        values.append([float(value) for value in row[1:]])
    if all(len(row) == len(traits) for row in values):
        max_error = max(
            abs(values[i][j] - values[j][i])
            for i in range(len(traits))
            for j in range(len(traits))
        )
        require(max_error <= 1e-10, f"matrix is not symmetric: {path} (max error {max_error})", errors)


def validate_sampling_matrix(
    path: Path, elements: list[tuple[str, str, str]], errors: list[str]
) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        errors.append(f"missing or empty sampling matrix: {path}")
        return
    ids = [row[0] for row in elements]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        require(header == ["element_id", *ids], f"wrong sampling-matrix header/order: {path}", errors)
        row_count = 0
        for row_count, row in enumerate(reader, start=1):
            if row_count <= len(ids):
                require(row and row[0] == ids[row_count - 1], f"wrong sampling row ID at row {row_count + 1}", errors)
            require(len(row) == len(ids) + 1, f"wrong sampling width at row {row_count + 1}", errors)
            if len(row) == len(ids) + 1:
                require(all(finite(value) for value in row[1:]), f"non-finite sampling value at row {row_count + 1}", errors)
        require(row_count == len(ids), f"wrong sampling-matrix row count: {row_count}", errors)


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    panel_path = root / "config/analysis_panel.tsv"
    panel = read_tsv(panel_path) if panel_path.is_file() else []
    traits = [row.get("trait_id", "") for row in panel]
    require(len(traits) == 45 and len(set(traits)) == 45, "panel is not 45 unique traits", errors)
    require(sum(row.get("domain") == "sleep" for row in panel) == 12, "panel is not 12 sleep traits", errors)
    if errors:
        return errors
    table_dir = root / "results/tables"
    for name in (
        "ldsc_covariance_45x45.tsv",
        "ldsc_genetic_correlation_45x45.tsv",
        "ldsc_intercept_45x45.tsv",
    ):
        validate_square_matrix(table_dir / name, traits, "trait_id", errors)
    elements = expected_elements(traits)
    validate_sampling_matrix(
        table_dir / "ldsc_sampling_covariance_1035x1035.tsv.gz", elements, errors
    )

    pair_path = table_dir / "ldsc_covariance_pairs.tsv"
    pairs = read_tsv(pair_path) if pair_path.is_file() else []
    required_pair_columns = {
        "element_id", "trait_1", "trait_2", "genetic_covariance",
        "genetic_covariance_se", "genetic_covariance_z", "genetic_covariance_p",
        "genetic_correlation", "genetic_correlation_se", "genetic_correlation_z",
        "genetic_correlation_p", "cross_trait_intercept", "effective_n",
        "scale_1", "scale_2", "genetic_covariance_fdr_off_diagonal",
        "genetic_correlation_fdr_off_diagonal",
    }
    observed_columns = set(pairs[0]) if pairs else set()
    require(required_pair_columns <= observed_columns, f"pair table missing columns: {sorted(required_pair_columns - observed_columns)}", errors)
    require(len(pairs) == len(elements), f"pair table has {len(pairs)} rows, expected {len(elements)}", errors)
    if len(pairs) == len(elements):
        observed = [(row.get("element_id"), row.get("trait_1"), row.get("trait_2")) for row in pairs]
        require(observed == elements, "pair table order or identities differ from locked panel", errors)
        numeric = required_pair_columns - {
            "element_id", "trait_1", "trait_2", "scale_1", "scale_2",
            "genetic_covariance_fdr_off_diagonal", "genetic_correlation_fdr_off_diagonal",
        }
        for row_number, row in enumerate(pairs, start=2):
            require(all(finite(row.get(column, "")) for column in numeric), f"non-finite pair estimate at row {row_number}", errors)
            if row.get("trait_1") == row.get("trait_2"):
                require(abs(float(row.get("genetic_correlation", "nan")) - 1.0) <= 1e-10, f"diagonal rg differs from 1 at row {row_number}", errors)
            else:
                require(finite(row.get("genetic_covariance_fdr_off_diagonal", "")), f"missing covariance FDR at row {row_number}", errors)
                require(finite(row.get("genetic_correlation_fdr_off_diagonal", "")), f"missing correlation FDR at row {row_number}", errors)

    scales_path = table_dir / "ldsc_covariance_trait_scales.tsv"
    scales = read_tsv(scales_path) if scales_path.is_file() else []
    require([row.get("trait_id") for row in scales] == traits, "trait-scale table differs from locked order", errors)
    metadata_path = table_dir / "ldsc_covariance_metadata.tsv"
    metadata_rows = read_tsv(metadata_path) if metadata_path.is_file() else []
    metadata = {row.get("field", ""): row.get("value", "") for row in metadata_rows}
    for field, value in {
        "panel_version": "atlas-v1.0", "trait_count": "45",
        "covariance_elements": "1035", "jackknife_blocks": "1082",
        "r_version": "4.3.3", "genomicsem_version": "0.0.5",
        "genomicsem_commit": "6b65ca5db39fdade08b0d811477be1cdd57b5039",
    }.items():
        require(metadata.get(field) == value, f"metadata {field!r} is {metadata.get(field)!r}, expected {value!r}", errors)
    diagnostics_path = table_dir / "ldsc_covariance_diagnostics.tsv"
    diagnostics = read_tsv(diagnostics_path) if diagnostics_path.is_file() else []
    require(len(diagnostics) == 4, "covariance diagnostics must contain four matrices", errors)
    require(all(row.get("all_finite") == "TRUE" for row in diagnostics), "diagnostics report non-finite matrices", errors)
    rds_path = table_dir / "ldsc_covariance_structure.rds"
    require(rds_path.is_file() and rds_path.stat().st_size > 0, "missing covariance RDS", errors)
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    errors = validate(Path(args.root).resolve())
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    if not args.quiet:
        print("PASS: complete locked 45-trait covariance export is structurally valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
