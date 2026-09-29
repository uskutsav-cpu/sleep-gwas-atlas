#!/usr/bin/env python3
"""Compare stored factor Gram blocks with LAVA's unit-diagonal block output.

Diagnostic only. It does not replace LD values, rebuild candidates, or promote
any locus. The normalized Gram matrix is a representation sensitivity only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1")
EDGE_PROV = ROOT / "brain6/qc/placo_ld_factor_normalization_v3.provenance.json"
EXCEPTIONS = ROOT / "brain6/results/loci/placo_candidate_ld_range_exceptions_partial_v2.tsv"
SCRIPT = ROOT / "brain6/scripts/audit_placo_ld_factor_normalization_v1.py"

spec = importlib.util.spec_from_file_location("lava_ld_factor_audit", SCRIPT)
factor_audit = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(factor_audit)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def block_matrices(factors: np.ndarray, eigenvalues: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return factor Gram, LAVA-style unit-diagonal, and normalized Gram matrices."""
    v = factors.astype(np.float64)
    eigen = eigenvalues.astype(np.float64)
    gram = v.T @ (eigen[:, None] * v)
    diagonal = np.diag(gram).copy()
    lava_style = gram.copy()
    np.fill_diagonal(lava_style, 1.0)
    normalized = gram / np.sqrt(np.outer(diagonal, diagonal))
    np.fill_diagonal(normalized, 1.0)
    return gram, lava_style, normalized


def read_edges(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(reference: Path, output: Path, provenance_output: Path) -> dict[str, object]:
    edge_prov = json.loads(EDGE_PROV.read_text(encoding="utf-8"))
    edge_table = ROOT / edge_prov["output"]["path"]
    if sha256(edge_table) != edge_prov["output"]["sha256"]:
        raise ValueError("factor-edge diagnostic changed since its checksum-bound audit")
    for path_text, identity in edge_prov["reference_files"].items():
        path = Path(path_text)
        if not path.is_file() or path.stat().st_size != identity["bytes"] or sha256(path) != identity["sha256"]:
            raise ValueError(f"reference checksum mismatch: {path}")
    for key in ("exception_table", "audit_script"):
        record = edge_prov["sources"][key]
        path = ROOT / record["path"]
        if sha256(path) != record["sha256"]:
            raise ValueError(f"factor-edge audit source changed: {path}")

    edges = read_edges(EXCEPTIONS)
    wanted_by_chr: dict[int, set[str]] = {}
    for row in edges:
        chrn = int(row["CHR"])
        wanted_by_chr.setdefault(chrn, set()).update((row["SNP1"].lower(), row["SNP2"].lower()))

    records: list[dict[str, object]] = []
    for chromosome, snps in sorted(wanted_by_chr.items()):
        snp_map, factor_blocks = factor_audit.load_primary_blocks(reference, chromosome, snps)
        block_edge_counts: dict[int, int] = {}
        for row in edges:
            if int(row["CHR"]) != chromosome:
                continue
            bid1, _ = snp_map[row["SNP1"].lower()]
            bid2, _ = snp_map[row["SNP2"].lower()]
            if bid1 != bid2:
                raise ValueError(f"exception edge crosses primary blocks: {row}")
            block_edge_counts[bid1] = block_edge_counts.get(bid1, 0) + 1
        for bid, (factors, eigenvalues, _) in sorted(factor_blocks.items()):
            gram, lava_style, normalized = block_matrices(factors, eigenvalues)
            eval_gram = np.linalg.eigvalsh(gram)
            eval_lava = np.linalg.eigvalsh(lava_style)
            eval_norm = np.linalg.eigvalsh(normalized)
            diag = np.diag(gram)
            mask = ~np.eye(gram.shape[0], dtype=bool)
            tol = 1e-9 * max(1.0, float(np.max(np.abs(eval_norm))))
            records.append({
                "CHR": chromosome, "primary_block": bid,
                "n_snps": factors.shape[1], "n_pcs": factors.shape[0],
                "exception_edges_in_block": block_edge_counts.get(bid, 0),
                "stored_eigenvalue_min": float(np.min(eigenvalues)),
                "stored_eigenvalue_max": float(np.max(eigenvalues)),
                "factor_diagonal_min": float(np.min(diag)),
                "factor_diagonal_median": float(np.median(diag)),
                "factor_diagonal_max": float(np.max(diag)),
                "factor_raw_offdiag_abs_r_gt_1": int(np.sum(np.abs(gram[mask]) > 1.0)),
                "unit_scaled_offdiag_abs_r_gt_1": int(np.sum(np.abs(normalized[mask]) > 1.0)),
                "factor_gram_min_eigenvalue": float(eval_gram[0]),
                "lava_unit_diagonal_min_eigenvalue": float(eval_lava[0]),
                "unit_scaled_gram_min_eigenvalue": float(eval_norm[0]),
                "unit_scaled_psd_with_numeric_tolerance": bool(eval_norm[0] >= -tol),
                "promoted_for_analysis": "NO",
            })
    if not records:
        raise ValueError("no primary blocks were selected")
    output, provenance_output = output.resolve(), provenance_output.resolve()
    if output.exists() or provenance_output.exists():
        raise FileExistsError("refusing to overwrite primary-block diagnostic outputs")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    psd_bad = sum(not row["unit_scaled_psd_with_numeric_tolerance"] for row in records)
    prov = {
        "analysis_id": "brain6-placo-ld-primary-block-psd-v1",
        "interpretation": "Diagnostic comparison of stored factor Gram blocks, LAVA-style off-diagonals with diagonal fixed to one, and factor-diagonal-normalized Gram blocks. Normalized blocks are a sensitivity representation, not independently validated genotype LD and not eligible for candidate processing.",
        "method": {"format_id": 14, "matrix_form": "V.T @ diag(stored_eigenvalues) @ V using decoded float32 factors promoted to float64", "lava_style": "same Gram off-diagonals with diagonal replaced by one", "normalized_sensitivity": "divide each Gram entry by sqrt(d_i*d_j), then set diagonal to one", "PSD_numeric_tolerance": "1e-9 * max(1, max(abs(eigenvalues(normalized_block))))", "no_clamping_or_promotion": True},
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "sources": {"factor_edge_output": {"path": str(edge_table.relative_to(ROOT)), "sha256": sha256(edge_table)},
                    "factor_edge_provenance": {"path": str(EDGE_PROV.relative_to(ROOT)), "sha256": sha256(EDGE_PROV)},
                    "exception_table": {"path": str(EXCEPTIONS.relative_to(ROOT)), "sha256": sha256(EXCEPTIONS)},
                    "factor_parser": {"path": str(SCRIPT.relative_to(ROOT)), "sha256": sha256(SCRIPT)},
                    "audit_script": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256(Path(__file__).resolve())}},
        "reference_files": edge_prov["reference_files"],
        "summary": {"blocks": len(records), "chromosomes": sorted(wanted_by_chr),
                    "n_snp_values": sorted({int(r["n_snps"]) for r in records}),
                    "rank_min": min(int(r["n_pcs"]) for r in records),
                    "rank_max": max(int(r["n_pcs"]) for r in records),
                    "blocks_with_positive_stored_eigenvalues": sum(float(r["stored_eigenvalue_min"]) > 0 for r in records),
                    "blocks_with_lava_style_negative_min_eigenvalue": sum(float(r["lava_unit_diagonal_min_eigenvalue"]) < 0 for r in records),
                    "blocks_normalized_psd_with_numeric_tolerance": len(records) - psd_bad,
                    "promoted_blocks": 0},
        "output": {"path": str(output.relative_to(ROOT)), "sha256": sha256(output), "rows": len(records)},
    }
    with provenance_output.open("x", encoding="utf-8") as stream:
        json.dump(prov, stream, indent=2, sort_keys=True); stream.write("\n")
    return prov


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.reference.resolve(), args.output, args.provenance)["summary"]
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
