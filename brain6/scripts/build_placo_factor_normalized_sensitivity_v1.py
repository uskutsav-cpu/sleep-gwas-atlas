"""Diagnostic PLACO clumping sensitivity from normalized LAVA factor Gram matrices.

This creates an isolated, non-promoted sensitivity. Within each official UKB
LAVA format-14 primary block, it reconstructs the stored low-rank Gram matrix
and converts that approximation to a correlation matrix by normalizing its
factor-derived diagonal. Cross-block edges remain the pinned LAVA reader
output. Only floating-point excursions up to 1e-12 beyond the mathematical
correlation boundary are normalized to exactly +/-1; material excursions fail.
No existing PLACO/LAVA result is replaced.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import clump_placo_candidates as clump
from audit_placo_ld_factor_normalization_v1 import load_primary_blocks
from rebuild_partial_placo_clumps_valid_ld import extract_windowed_ld, read_tsv

RULE_PATH = ROOT / "brain6/config/shared_locus_rule_v1.json"
FAMILY_PATH = ROOT / "brain6/config/placo_family_v3/family_lock.json"
RUN_PATH = ROOT / "brain6/manifests/placo_family_v3_run.json"
MASTER_PATH = ROOT / "brain6/results/placo/placo_master.tsv"
REFERENCE_PREFIX = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
VALID_LD_PROVENANCE = ROOT / "brain6/results/loci/placo_candidate_valid_ld_partial_v2.provenance.json"
OUT_DEFAULT = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clip_roundoff_correlation_boundary(correlation: np.ndarray,
                                      tolerance: float = 1e-12) -> tuple[np.ndarray, int]:
    """Clamp only finite correlation values within tolerance of +/-1."""
    values = np.asarray(correlation, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("normalized correlation matrix is nonfinite")
    excess = np.abs(values) - 1.0
    if np.max(excess) > tolerance:
        raise ValueError("normalized factor correlation materially exceeds [-1, 1]")
    roundoff = np.abs(values) > 1.0
    normalized = values.copy()
    normalized[roundoff] = np.sign(normalized[roundoff])
    return normalized, int(np.count_nonzero(roundoff))


def normalized_correlation(factors: np.ndarray, eigenvalues: np.ndarray,
                           *, report_roundoff: bool = False) -> np.ndarray | tuple[np.ndarray, int]:
    """Convert a low-rank factor Gram matrix to a unit-diagonal correlation matrix.

    Correlations outside [-1, 1] by at most 1e-12 are mapped to the exact
    mathematical boundary. This only resolves floating-point representation
    error in an otherwise valid Gram matrix; material violations fail closed.
    """
    factors64 = np.asarray(factors, dtype=np.float64)
    eig64 = np.asarray(eigenvalues, dtype=np.float64)
    if factors64.ndim != 2 or eig64.ndim != 1 or factors64.shape[0] != eig64.size:
        raise ValueError("factor/eigenvalue dimensions disagree")
    if not np.all(np.isfinite(factors64)) or not np.all(np.isfinite(eig64)) or np.any(eig64 < 0):
        raise ValueError("factor block contains nonfinite or negative eigenvalues")
    gram = factors64.T @ (eig64[:, None] * factors64)
    diagonal = np.diag(gram)
    if np.any(~np.isfinite(diagonal)) or np.any(diagonal <= 0):
        raise ValueError("factor-derived diagonal must be finite and positive")
    scale = np.sqrt(diagonal)
    correlation = gram / np.outer(scale, scale)
    correlation = (correlation + correlation.T) * 0.5
    if not np.all(np.isfinite(correlation)):
        raise ValueError("normalized correlation matrix is nonfinite")
    if np.max(np.abs(np.diag(correlation) - 1.0)) > 1e-12:
        raise ValueError("normalized correlation matrix does not have unit diagonal")
    np.fill_diagonal(correlation, 1.0)
    # Unit-diagonal factor Gram matrices are PSD and have correlations bounded
    # by one. Values a few ULPs beyond the boundary are arithmetic roundoff,
    # not observed LD to clamp. Preserve all in-range values exactly.
    correlation, roundoff_count = clip_roundoff_correlation_boundary(correlation)
    smallest = float(np.linalg.eigvalsh(correlation)[0])
    if smallest < -1e-9 * max(1.0, float(np.linalg.norm(correlation, ord=2))):
        raise ValueError(f"normalized factor matrix is not PSD: min eigenvalue={smallest}")
    if report_roundoff:
        return correlation, roundoff_count
    return correlation


def edge_chromosome(snp1: str, snp2: str,
                    chromosomes_by_snp: dict[str, set[int]]) -> int | None:
    """Resolve an edge's chromosome in constant time from a prebuilt index."""
    common = chromosomes_by_snp.get(snp1, set()) & chromosomes_by_snp.get(snp2, set())
    if len(common) > 1:
        raise ValueError(f"candidate edge has ambiguous chromosome by SNP ID: {snp1}/{snp2}")
    return next(iter(common)) if common else None


def normalized_candidate_edges(candidates: list[dict[str, Any]], reference: Path,
                               ref_info: dict, raw_ld: dict[tuple[str, str], float],
                               window_bp: int) -> tuple[dict[tuple[str, str], float], list[dict[str, Any]], dict[str, Any]]:
    """Replace same-primary-block edges with normalized factor correlations."""
    exact: dict[int, set[str]] = {}
    candidate_bp: dict[tuple[int, str], int] = {}
    for row in candidates:
        chr_id, snp, bp = int(row["CHR"]), str(row["SNP"]).lower(), int(row["BP"])
        candidate_bp[(chr_id, snp)] = bp
        match = ref_info.get((chr_id, snp))
        if match is not None and int(match["BP"]) == bp:
            exact.setdefault(chr_id, set()).add(snp)

    block_for: dict[tuple[int, str], tuple[int, int]] = {}
    matrices: dict[tuple[int, int], np.ndarray] = {}
    roundoff_boundary_values_normalized = 0
    for chr_id, wanted in sorted(exact.items()):
        snp_to_block, factors_by_block = load_primary_blocks(reference, chr_id, wanted)
        for snp, loc in snp_to_block.items():
            block_for[(chr_id, snp)] = loc
        snps_by_block: dict[int, list[str]] = {}
        for snp, (block_id, _) in snp_to_block.items():
            snps_by_block.setdefault(block_id, []).append(snp)
        for block_id in snps_by_block:
            factors, eigenvalues, _start = factors_by_block[block_id]
            corr, n_roundoff = normalized_correlation(factors, eigenvalues, report_roundoff=True)
            roundoff_boundary_values_normalized += n_roundoff
            matrices[(chr_id, block_id)] = corr
    chromosomes_by_snp: dict[str, set[int]] = {}
    for chr_id, snp in block_for:
        chromosomes_by_snp.setdefault(snp, set()).add(chr_id)

    out = dict(raw_ld)
    replacements: list[dict[str, Any]] = []
    cross_block_invalid = []
    same_block_edges_replaced = 0
    for (snp1, snp2), raw_r in raw_ld.items():
        chr_id = edge_chromosome(snp1, snp2, chromosomes_by_snp)
        if chr_id is None:
            continue
        left = block_for.get((chr_id, snp1))
        right = block_for.get((chr_id, snp2))
        if left is None or right is None or left[0] != right[0]:
            if abs(raw_r) > 1.0:
                cross_block_invalid.append((chr_id, snp1, snp2, raw_r))
            continue
        if abs(candidate_bp[(chr_id, snp1)] - candidate_bp[(chr_id, snp2)]) > window_bp:
            continue
        corr = matrices[(chr_id, left[0])]
        # The diagonal of a correlation matrix is defined to be exactly one;
        # evaluate this identity directly instead of carrying roundoff from
        # the same factor norm through a division.
        factor_r = 1.0 if snp1 == snp2 else float(corr[left[1], right[1]])
        if abs(factor_r) > 1.0 + 1e-12:
            raise ValueError(f"normalized correlation outside range for {snp1}/{snp2}")
        if abs(factor_r) > 1.0:
            factor_r = float(np.sign(factor_r))
        key = tuple(sorted((snp1, snp2)))
        out[key] = factor_r
        same_block_edges_replaced += 1
        if abs(raw_r) > 1.0:
            replacements.append({"CHR": chr_id, "SNP1": snp1, "SNP2": snp2,
                                 "BP1": candidate_bp[(chr_id, snp1)], "BP2": candidate_bp[(chr_id, snp2)],
                                 "lava_r_unscaled": raw_r, "factor_r_normalized": factor_r,
                                 "factor_r2_normalized": factor_r * factor_r,
                                 "diagnostic_status": "NORMALIZED_FACTOR_SENSITIVITY_ONLY_NOT_GENOTYPE_VALIDATED"})
    if cross_block_invalid:
        raise ValueError(f"found {len(cross_block_invalid)} invalid cross-block LD edges; cannot normalize by primary block")
    summary = {"exact_candidate_snps": sum(map(len, exact.values())),
               "primary_blocks_touched": len(matrices),
               "same_block_edges_replaced": same_block_edges_replaced,
               "roundoff_boundary_values_normalized": roundoff_boundary_values_normalized,
               "invalid_cross_block_edges": 0}
    return out, replacements, summary


def atomic_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(reference_prefix: Path = REFERENCE_PREFIX, rscript: Path = RSCRIPT,
          output_dir: Path = OUT_DEFAULT,
          analysis_id: str = "brain6_placo_factor_normalized_sensitivity_v1") -> dict[str, Any]:
    rule = json.loads(RULE_PATH.read_text(encoding="utf-8"))
    family = json.loads(FAMILY_PATH.read_text(encoding="utf-8"))
    candidates, input_hashes = clump.read_candidates(MASTER_PATH, float(family["headline_threshold"]))
    all_rows = [row for pair in sorted(candidates) for row in candidates[pair]]
    window_bp = int(rule["clumping"]["window_kb"]) * 1000
    ref_info, raw_ld, ref_n = extract_windowed_ld(all_rows, reference_prefix, rscript, window_bp)
    if ref_n != int(rule["clumping"]["reference_sample_size"]):
        raise ValueError(f"reference sample-size mismatch: {ref_n}")
    normalized_ld, replacements, matrix_summary = normalized_candidate_edges(
        all_rows, reference_prefix.parent, ref_info, raw_ld, window_bp)
    variants, loci = clump.build_outputs(candidates, ref_info, normalized_ld, rule)
    out = output_dir
    paths = {"variants": out / "placo_candidate_variants.tsv",
             "loci": out / "placo_candidate_loci.tsv",
             "replaced_edges": out / "normalized_range_exceptions.tsv",
             "summary": out / "pair_summary.tsv",
             "provenance": out / "provenance.json"}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError(f"refusing to overwrite sensitivity artifacts in {out}")
    original_variants = read_tsv(ROOT / "brain6/results/loci/placo_candidate_variants_partial.tsv")
    original_loci = read_tsv(ROOT / "brain6/results/loci/placo_candidate_loci_partial.tsv")
    valid_variants = read_tsv(ROOT / "brain6/results/loci/placo_candidate_variants_valid_ld_partial_v2.tsv")
    valid_loci = read_tsv(ROOT / "brain6/results/loci/placo_candidate_loci_valid_ld_partial_v2.tsv")
    summary = []
    for pair in sorted(candidates):
        pair_variants = [r for r in variants if r["pair_id"] == pair]
        pair_original_variants = [r for r in original_variants if r["pair_id"] == pair]
        pair_valid_variants = [r for r in valid_variants if r["pair_id"] == pair]
        pair_original_loci = [r for r in original_loci if r["pair_id"] == pair]
        pair_valid_loci = [r for r in valid_loci if r["pair_id"] == pair]
        pair_loci = [r for r in loci if r["pair_id"] == pair]
        pair_replacements = [r for r in replacements if any(
            c["SNP"].lower() == r["SNP1"] for c in candidates[pair]) and any(
            c["SNP"].lower() == r["SNP2"] for c in candidates[pair])]
        summary.append({"pair_id": pair,
                        "candidate_rows": len(pair_variants),
                        "original_leads": sum(r["candidate_status"] == "LEAD" for r in pair_original_variants),
                        "exclude_invalid_ld_leads": sum(r["candidate_status"] == "LEAD" for r in pair_valid_variants),
                        "factor_normalized_sensitivity_leads": sum(r["candidate_status"] == "LEAD" for r in pair_variants),
                        "original_intervals": len(pair_original_loci),
                        "exclude_invalid_ld_intervals": len(pair_valid_loci),
                        "factor_normalized_sensitivity_intervals": len(pair_loci),
                        "normalized_exception_edges": len(pair_replacements),
                        "interpretation": "Diagnostic clumping sensitivity only; normalized low-rank factors are not genotype-level validated LD."})
    atomic_tsv(paths["variants"], clump.VARIANT_COLUMNS, variants)
    atomic_tsv(paths["loci"], clump.LOCUS_COLUMNS, loci)
    edge_cols = ["CHR", "SNP1", "BP1", "SNP2", "BP2", "lava_r_unscaled", "factor_r_normalized", "factor_r2_normalized", "diagnostic_status"]
    atomic_tsv(paths["replaced_edges"], edge_cols, replacements)
    atomic_tsv(paths["summary"], list(summary[0]), summary)
    ref_provenance = reference_prefix.parent / "reference.provenance.json"
    ref_files = {}
    chromosomes = sorted({int(row["CHR"]) for row in all_rows})
    for chr_id in chromosomes:
        for suffix in ("info", "bcor"):
            path = reference_prefix.parent / f"{reference_prefix.name}_chr{chr_id}.{suffix}"
            ref_files[str(path)] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    result = {"analysis_id": analysis_id,
              "status": "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED",
              "interpretation": "Within primary format-14 blocks, the LAVA reader's unit-diagonal/off-diagonal factor-product matrix is replaced by the correlation normalization of the same low-rank factor Gram matrix. Floating-point excursions <=1e-12 outside [-1,1] are normalized to the exact mathematical boundary; larger excursions fail closed. This is an internal matrix sensitivity, not exact genotype-level LD validation. No prior or canonical files are modified.",
              "scope": {"pairs": sorted(candidates), "candidate_rows": len(all_rows), "candidate_variants": len(variants),
                        "candidate_intervals": len(loci), "reference_N": ref_n, **matrix_summary,
                        "range_exception_edges_normalized": len(replacements),
                        "full_five_pair_family": False},
              "method": {"source": "Pinned LAVA 0.1.5 UKB v1.1 format-14 eigen factors", "window_bp": window_bp,
                         "normalization": "G = V diag(lambda) V^T; R_ij = G_ij/sqrt(G_ii*G_jj); unit diagonal; only <=1e-12 floating-point excursions are normalized to +/-1, larger violations fail",
                         "within_block": "factor-normalized correlation", "cross_block": "pinned LAVA read.ld output unchanged",
                         "external_genotype_validation": "NOT_AVAILABLE", "promotion": "NONE"},
              "sources": {"master": {"path": str(MASTER_PATH.relative_to(ROOT)), "sha256": sha256(MASTER_PATH)},
                          "builder": {"path": "brain6/scripts/build_placo_factor_normalized_sensitivity_v1.py", "sha256": sha256(Path(__file__).resolve())},
                          "family_lock": {"path": str(FAMILY_PATH.relative_to(ROOT)), "sha256": sha256(FAMILY_PATH)},
                          "run_manifest": {"path": str(RUN_PATH.relative_to(ROOT)), "sha256": sha256(RUN_PATH)},
                          "locus_rule": {"path": str(RULE_PATH.relative_to(ROOT)), "sha256": sha256(RULE_PATH)},
                          "candidate_input_hashes": input_hashes,
                          "prior_valid_ld_provenance": {"path": str(VALID_LD_PROVENANCE.relative_to(ROOT)), "sha256": sha256(VALID_LD_PROVENANCE)},
                          "reference_provenance": {"path": str(ref_provenance), "sha256": sha256(ref_provenance)},
                          "reference_files": ref_files,
                          "factor_reader_audit_script": {"path": "brain6/scripts/audit_placo_ld_factor_normalization_v1.py", "sha256": sha256(ROOT / "brain6/scripts/audit_placo_ld_factor_normalization_v1.py")},
                          "factor_psd_audit": {"path": "brain6/qc/placo_ld_primary_block_psd_v2.provenance.json", "sha256": sha256(ROOT / "brain6/qc/placo_ld_primary_block_psd_v2.provenance.json")},
                          "lava_reference_docs": "https://github.com/josefin-werme/LAVA/blob/main/REFERENCE.md",
                          "lava_reader_source": "https://github.com/josefin-werme/LAVA/blob/main/src/load_ld.cpp"},
              "outputs": []}
    for name, path in paths.items():
        if name != "provenance":
            result["outputs"].append({"name": name, "path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size})
    paths["provenance"].parent.mkdir(parents=True, exist_ok=True)
    with paths["provenance"].open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True); stream.write("\n")
    result["provenance_sha256"] = sha256(paths["provenance"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-prefix", type=Path, default=REFERENCE_PREFIX)
    parser.add_argument("--rscript", type=Path, default=RSCRIPT)
    parser.add_argument("--output-dir", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--analysis-id", default="brain6_placo_factor_normalized_sensitivity_v1")
    args = parser.parse_args()
    print(json.dumps(build(args.reference_prefix.resolve(), args.rscript.resolve(),
                           args.output_dir.resolve(), args.analysis_id), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
