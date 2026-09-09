#!/usr/bin/env python3
"""Prepare and validate the locked 45-trait LAVA input contracts."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import lava_contract


REQUIRED_SUMSTATS_COLUMNS = {"SNP", "A1", "A2", "Z", "N"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def load_panel(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    traits = [row.get("trait_id", "") for row in rows]
    domains = [row.get("domain", "") for row in rows]
    if len(rows) != 45 or len(set(traits)) != 45:
        raise SystemExit("ERROR: LAVA requires the exact locked 45-trait panel")
    if domains.count("sleep") != 12 or len(domains) - domains.count("sleep") != 33:
        raise SystemExit("ERROR: expected 12 sleep and 33 non-sleep traits")
    return rows


def validate_sumstats(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing LAVA summary statistics: {path}")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    missing = REQUIRED_SUMSTATS_COLUMNS.difference(header)
    if missing:
        raise SystemExit(f"ERROR: {path} lacks LAVA columns: {sorted(missing)}")


def prepare_info(panel: list[dict[str, str]], munged_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    info: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    for row in panel:
        trait = row["trait_id"]
        path = munged_dir / f"{trait}.sumstats.gz"
        validate_sumstats(path)
        if row["type"] == "binary":
            try:
                cases = int(row["ncase"])
                controls = int(row["ncontrol"])
                prevalence = float(row["pop_prev"])
            except ValueError as exc:
                raise SystemExit(f"ERROR: unresolved binary LAVA metadata for {trait}") from exc
            if cases <= 0 or controls <= 0 or not 0 < prevalence < 1:
                raise SystemExit(f"ERROR: invalid binary LAVA metadata for {trait}")
        elif row["type"] == "continuous":
            cases = controls = prevalence = "NA"
        else:
            raise SystemExit(f"ERROR: unsupported trait type for {trait}: {row['type']}")
        info.append({
            "phenotype": trait,
            "cases": cases,
            "controls": controls,
            "prevalence": prevalence,
            "filename": path.as_posix(),
        })
        provenance.append({
            "phenotype": trait,
            "filename": path.as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "required_columns": ",".join(sorted(REQUIRED_SUMSTATS_COLUMNS)),
            "validation_status": "VALIDATED",
        })
    return info, provenance


def prepare_overlap(path: Path, traits: list[str]) -> tuple[np.ndarray, dict[str, object]]:
    rows = read_tsv(path)
    if len(rows) != 45 or [row.get("trait_id") for row in rows] != traits:
        raise SystemExit("ERROR: LDSC intercept matrix rows do not match locked panel order")
    if not rows or list(rows[0].keys())[1:] != traits:
        raise SystemExit("ERROR: LDSC intercept matrix columns do not match locked panel order")
    matrix = np.array([[float(row[trait]) for trait in traits] for row in rows], dtype=float)
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, atol=1e-10):
        raise SystemExit("ERROR: LDSC intercept matrix must be finite and symmetric")
    diagonal = np.diag(matrix)
    if np.any(diagonal <= 0):
        raise SystemExit("ERROR: LDSC intercept diagonal must be positive")
    correlation = matrix / np.sqrt(np.outer(diagonal, diagonal))
    correlation = np.round(correlation, 5)
    correlation = (correlation + correlation.T) / 2
    np.fill_diagonal(correlation, 1.0)
    eigenvalues = np.linalg.eigvalsh(correlation)
    if eigenvalues[0] <= 0:
        raise SystemExit(f"ERROR: sample-overlap correlation is not positive definite: {eigenvalues[0]}")
    diagnostics = {
        "traits": 45,
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "condition_number": float(eigenvalues[-1] / eigenvalues[0]),
        "minimum_off_diagonal": float(correlation[np.triu_indices(45, 1)].min()),
        "maximum_off_diagonal": float(correlation[np.triu_indices(45, 1)].max()),
    }
    return correlation, diagnostics


def write_overlap(path: Path, traits: list[str], matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        handle.write(" " + " ".join(traits) + "\n")
        for trait, values in zip(traits, matrix):
            handle.write(trait + " " + " ".join(f"{value:.5f}" for value in values) + "\n")
    temporary.replace(path)


def prepare_pairs(path: Path, panel: list[dict[str, str]]) -> list[dict[str, object]]:
    rg = read_tsv(path)
    sleep = {row["trait_id"] for row in panel if row["domain"] == "sleep"}
    non_sleep = {row["trait_id"] for row in panel if row["domain"] != "sleep"}
    expected = {(first, second) for first in sleep for second in non_sleep}
    observed = {(row.get("sleep_trait", ""), row.get("disease_trait", "")) for row in rg}
    if len(rg) != 396 or observed != expected:
        raise SystemExit("ERROR: rg table is not the exact locked 396-pair family")
    rows = []
    for row in rg:
        rows.append({
            "sleep_trait": row["sleep_trait"],
            "non_sleep_trait": row["disease_trait"],
            "analysis_tier": row["analysis_tier"],
            "interpretation_status": row["interpretation_status"],
            "global_rg": row["rg"],
            "global_rg_se": row["se"],
            "global_rg_p": row["p"],
            "global_rg_fdr": row["fdr"],
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", default="config/analysis_panel.tsv")
    parser.add_argument("--munged-dir", default="data/munged")
    parser.add_argument("--intercepts", default="results/tables/ldsc_intercept_45x45.tsv")
    parser.add_argument("--rg", default="results/tables/rg_matrix.tsv")
    parser.add_argument("--policy", default="config/lava_analysis_policy.json")
    parser.add_argument("--input-info", default="results/tables/lava_input_info.tsv")
    parser.add_argument("--sample-overlap", default="results/tables/lava_sample_overlap.txt")
    parser.add_argument("--pairs", default="results/tables/lava_pair_manifest.tsv")
    parser.add_argument("--provenance", default="results/tables/lava_input_provenance.tsv")
    parser.add_argument("--runtime-policy", default="results/tables/lava_runtime_policy.tsv")
    parser.add_argument("--diagnostics", default="results/tables/lava_input_diagnostics.json")
    parser.add_argument("--lock", default="results/tables/lava_input.lock.json")
    args = parser.parse_args()

    root = Path(".").resolve()
    policy_path = Path(args.policy).resolve()
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if policy.get("expected_traits") != 45 or policy.get("expected_sleep_non_sleep_pairs") != 396:
        raise SystemExit("ERROR: LAVA policy differs from the locked atlas scope")
    requested_artifacts = [
        args.input_info, args.sample_overlap, args.pairs, args.provenance,
        args.runtime_policy, args.diagnostics,
    ]
    if (
        [str(Path(path)) for path in requested_artifacts] != policy["input_artifacts"]
        or str(Path(args.lock)) != policy["input_lock"]
        or policy_path != root / "config/lava_analysis_policy.json"
    ):
        raise SystemExit("ERROR: LAVA preparation paths differ from the immutable policy")
    lock_path = root / str(policy["input_lock"])
    if lock_path.exists():
        lava_contract.validate_inputs(root, policy_path, policy)
        print("Validated existing immutable LAVA input family; no artifacts were replaced")
        return 0
    expected_threshold = policy["univariate_alpha"] / policy["planned_univariate_tests"]
    if not math.isclose(expected_threshold, policy["univariate_p_threshold"], rel_tol=1e-12):
        raise SystemExit("ERROR: inconsistent LAVA univariate threshold")

    panel = load_panel(Path(args.panel))
    traits = [row["trait_id"] for row in panel]
    info, provenance = prepare_info(panel, Path(args.munged_dir))
    overlap, overlap_diagnostics = prepare_overlap(Path(args.intercepts), traits)
    pairs = prepare_pairs(Path(args.rg), panel)

    write_tsv(Path(args.input_info), ["phenotype", "cases", "controls", "prevalence", "filename"], info)
    write_overlap(Path(args.sample_overlap), traits, overlap)
    write_tsv(
        Path(args.pairs),
        ["sleep_trait", "non_sleep_trait", "analysis_tier", "interpretation_status", "global_rg", "global_rg_se", "global_rg_p", "global_rg_fdr"],
        pairs,
    )
    write_tsv(Path(args.provenance), ["phenotype", "filename", "bytes", "sha256", "required_columns", "validation_status"], provenance)
    runtime_keys = [
        "analysis_id", "lava_version", "reference_prefix", "locus_definition",
        "expected_loci", "expected_traits", "expected_sleep_non_sleep_pairs",
        "planned_univariate_tests", "univariate_p_threshold", "bivariate_fdr_alpha",
        "max_locus_failure_fraction", "max_univariate_untested_fraction", "random_seed",
        "max_bivariate_failure_fraction",
        "min_K", "prune_threshold", "max_proportion_K", "max_block_size", "cap_estimates",
    ]
    write_tsv(
        Path(args.runtime_policy),
        ["key", "value"],
        [{"key": key, "value": str(policy[key]).lower() if isinstance(policy[key], bool) else policy[key]} for key in runtime_keys],
    )
    diagnostics = {
        "analysis_id": policy["analysis_id"],
        "input_traits": len(info),
        "input_pairs": len(pairs),
        "planned_loci": policy["expected_loci"],
        "planned_univariate_tests": policy["planned_univariate_tests"],
        "sample_overlap": overlap_diagnostics,
        "input_info_sha256": sha256(Path(args.input_info)),
        "sample_overlap_sha256": sha256(Path(args.sample_overlap)),
        "pair_manifest_sha256": sha256(Path(args.pairs)),
        "input_provenance_sha256": sha256(Path(args.provenance)),
        "runtime_policy_sha256": sha256(Path(args.runtime_policy)),
    }
    diagnostics_path = Path(args.diagnostics)
    diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = diagnostics_path.with_suffix(diagnostics_path.suffix + ".tmp")
    temporary.write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(diagnostics_path)
    lava_contract.write_input_lock(root)
    print("Prepared validated LAVA inputs for 45 traits, 2,495 loci, and 396 locked pairs")
    print(f"Sample-overlap minimum eigenvalue: {overlap_diagnostics['minimum_eigenvalue']:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
