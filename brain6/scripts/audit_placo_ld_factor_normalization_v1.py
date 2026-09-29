#!/usr/bin/env python3
"""Read-only diagnostic of LAVA eigen-decomposed LD range exceptions.

This audit compares native LAVA output with the stored float factor products
and reports a unit-diagonal rescaling as a diagnostic only. It never edits LD,
candidate tables, or any production result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1")
EXCEPTIONS = ROOT / "brain6/results/loci/placo_candidate_ld_range_exceptions_partial_v2.tsv"
REPLAY = ROOT / "brain6/qc/placo_ld_range_exception_replay_v1.tsv"
REPLAY_PROVENANCE = ROOT / "brain6/qc/placo_ld_range_exception_replay_v1.provenance.json"
OUTPUT = ROOT / "brain6/qc/placo_ld_factor_normalization_v1.tsv"
PROVENANCE = ROOT / "brain6/qc/placo_ld_factor_normalization_v1.provenance.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bcor_layout(path: Path) -> dict[str, int]:
    """Decode the public LAVA v1.1 format-14 header and fail closed."""
    size = path.stat().st_size
    with path.open("rb") as stream:
        header = stream.read(4 * 19)
    if len(header) != 76:
        raise ValueError(f"short BCOR header: {path}")
    fields = struct.unpack("<19I", header)
    check, header_words, n_snps, value_type, value_size, fmt = fields[:6]
    p = fields[6:]
    index0 = p[0] | (p[1] << 32)
    index1 = p[2] | (p[3] << 32)
    n_blocks, n_entries = p[4], p[5]
    if (check, fmt, value_type, value_size) != (21775, 14, 3, 4):
        raise ValueError(f"unexpected LAVA BCOR format/header in {path}")
    if header_words * 4 != 76 or index1 - index0 != 8 * n_blocks:
        raise ValueError(f"inconsistent LAVA BCOR indexes in {path}")
    if size != index1 + 20 * n_entries or index0 < 76:
        raise ValueError(f"inconsistent LAVA BCOR file size in {path}")
    return {"n_snps": n_snps, "data_offset": 76, "block_index_offset": index0,
            "entry_index_offset": index1, "n_blocks": n_blocks,
            "n_entries": n_entries, "file_size": size}


def factor_correlation(factors: np.ndarray, eigenvalues: np.ndarray, i: int, j: int) -> tuple[float, float, float, float]:
    """Return raw product, its two factor diagonals, and diagnostic unit scaling."""
    vi = factors[:, i].astype(np.float64)
    vj = factors[:, j].astype(np.float64)
    eig = eigenvalues.astype(np.float64)
    raw = float(np.sum((vi * eig) * vj, dtype=np.float64))
    di = float(np.sum((vi * eig) * vi, dtype=np.float64))
    dj = float(np.sum((vj * eig) * vj, dtype=np.float64))
    scaled = raw / math.sqrt(di * dj)
    return raw, di, dj, scaled


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_primary_blocks(reference: Path, chromosome: int, wanted: set[str]) -> tuple[dict[str, tuple[int, int]], dict[int, tuple[np.ndarray, np.ndarray, int]]]:
    prefix = reference / f"lava-ukb-v1.1_chr{chromosome}"
    info_path, bcor_path = Path(str(prefix) + ".info"), Path(str(prefix) + ".bcor")
    with info_path.open(newline="", encoding="utf-8") as stream:
        snps = [row["SNP"].lower() for row in csv.DictReader(stream, delimiter="\t")]
    if len(snps) != len(set(snps)):
        raise ValueError(f"duplicate SNP identifiers in {info_path}")
    positions = {snp: idx for idx, snp in enumerate(snps)}
    missing = wanted - positions.keys()
    if missing:
        raise ValueError(f"{len(missing)} exception SNPs absent on chr{chromosome}")
    layout = bcor_layout(bcor_path)
    if len(snps) != layout["n_snps"]:
        raise ValueError(f".info/.bcor SNP count mismatch on chr{chromosome}")

    with bcor_path.open("rb") as stream:
        stream.seek(layout["block_index_offset"])
        block_raw = np.fromfile(stream, dtype="<u4", count=2 * layout["n_blocks"])
        stream.seek(layout["entry_index_offset"])
        entries_raw = np.fromfile(stream, dtype="<u4", count=5 * layout["n_entries"])
    if block_raw.size != 2 * layout["n_blocks"] or entries_raw.size != 5 * layout["n_entries"]:
        raise ValueError(f"truncated BCOR index on chr{chromosome}")
    blocks = block_raw.reshape((2, layout["n_blocks"]), order="F")
    entries = entries_raw.reshape((5, layout["n_entries"]), order="F")
    wanted_indices = {positions[snp] for snp in wanted}
    block_of = np.empty(layout["n_snps"], dtype=np.int64)
    block_sizes: dict[int, tuple[int, int, int]] = {}
    cursor = 0
    for bid, (n_snp, n_pc) in enumerate(blocks.T):
        block_of[cursor:cursor + int(n_snp)] = bid
        block_sizes[bid] = (int(n_snp), int(n_pc), cursor)
        cursor += int(n_snp)
    if cursor != layout["n_snps"]:
        raise ValueError("block SNP counts do not sum to chromosome SNP count")
    requested_blocks = {int(block_of[idx]) for idx in wanted_indices}
    payload_offsets = np.cumsum(np.r_[0, entries[3].astype(np.uint64)])
    if int(payload_offsets[-1]) != layout["block_index_offset"] - layout["data_offset"]:
        raise ValueError("primary/product payload size disagrees with index")
    result: dict[int, tuple[np.ndarray, np.ndarray, int]] = {}
    for entry_id in range(layout["n_entries"]):
        bid1, bid2 = int(entries[0, entry_id]), int(entries[1, entry_id])
        if bid1 != bid2 or bid1 not in requested_blocks:
            continue
        nsnp, npc, start = block_sizes[bid1]
        if int(entries[2, entry_id]) != npc * (nsnp + 1) or int(entries[3, entry_id]) != 4 * npc * (nsnp + 1):
            raise ValueError(f"invalid primary block dimensions for chr{chromosome} block {bid1}")
        with bcor_path.open("rb") as stream:
            stream.seek(layout["data_offset"] + int(payload_offsets[entry_id]))
            packed = np.fromfile(stream, dtype="<f4", count=npc * (nsnp + 1))
        data = packed.reshape((npc, nsnp + 1), order="F")
        result[bid1] = (data[:, :nsnp].copy(), data[:, nsnp].copy(), start)
    if set(result) != requested_blocks:
        raise ValueError(f"missing primary factor blocks on chr{chromosome}")
    snp_block_local = {snps[idx]: (int(block_of[idx]), idx - block_sizes[int(block_of[idx])][2]) for idx in wanted_indices}
    return snp_block_local, result


def run(reference: Path, output: Path, provenance_out: Path) -> dict[str, object]:
    replay_prov = json.loads(REPLAY_PROVENANCE.read_text())
    expected_exceptions = replay_prov["inputs"]["exception_table"]["sha256"]
    if sha256(EXCEPTIONS) != expected_exceptions or sha256(REPLAY) != replay_prov["output"]["sha256"]:
        raise ValueError("archived exception or LAVA replay checksum changed")
    exceptions = read_tsv(EXCEPTIONS)
    replay_rows = read_tsv(REPLAY)
    replay_map: dict[tuple[int, str, str], float] = {}
    for row in exceptions:
        key = (int(row["CHR"]), row["SNP1"].lower(), row["SNP2"].lower())
        if key in replay_map and not math.isclose(replay_map[key], float(row["lava_r_raw"]), abs_tol=0.0, rel_tol=0.0):
            raise ValueError(f"duplicate edge has inconsistent raw values: {key}")
        replay_map[key] = float(row["lava_r_raw"])
    if len(replay_map) != int(replay_prov["scope"]["unique_edges"]):
        raise ValueError("exception edge count differs from pinned replay provenance")
    if sum(int(row["exception_edges"]) for row in replay_rows) != len(replay_map):
        raise ValueError("replay table does not account for every unique exception edge")

    outputs: list[dict[str, object]] = []
    ref_hashes: dict[str, dict[str, object]] = {}
    for chromosome in sorted({k[0] for k in replay_map}):
        keys = [k for k in replay_map if k[0] == chromosome]
        wanted = {s for _, a, b in keys for s in (a, b)}
        snp_to_block, factors_by_block = load_primary_blocks(reference, chromosome, wanted)
        for key in keys:
            _, snp1, snp2 = key
            b1, i = snp_to_block[snp1]
            b2, j = snp_to_block[snp2]
            if b1 != b2:
                raise ValueError(f"exception edge spans factor blocks: {key}")
            factors, eigenvalues, _ = factors_by_block[b1]
            raw, d1, d2, scaled = factor_correlation(factors, eigenvalues, i, j)
            if abs(raw - replay_map[key]) > 1e-12:
                raise ValueError(f"factor product disagrees with LAVA replay for {key}: {raw} vs {replay_map[key]}")
            outputs.append({"CHR": chromosome, "SNP1": snp1, "SNP2": snp2,
                            "factor_block": b1, "factor_rank": factors.shape[0],
                            "lava_r_raw": replay_map[key], "factor_product_r": raw,
                            "factor_diagonal_snp1": d1, "factor_diagonal_snp2": d2,
                            "unit_diagonal_scaled_r_diagnostic_only": scaled,
                            "factor_raw_delta": abs(raw - replay_map[key]),
                            "promoted_for_analysis": "NO"})
        for suffix in ("info", "bcor"):
            ref_path = Path(f"{reference}/lava-ukb-v1.1_chr{chromosome}.{suffix}")
            prior = replay_prov["inputs"]["reference_files"].get(str(ref_path))
            actual_hash = sha256(ref_path)
            if prior and actual_hash != prior["sha256"]:
                raise ValueError(f"reference file checksum differs from pinned replay: {ref_path}")
            ref_hashes[str(ref_path)] = {"bytes": ref_path.stat().st_size, "sha256": actual_hash}

    outputs.sort(key=lambda r: (int(r["CHR"]), str(r["SNP1"]), str(r["SNP2"])))
    if output.exists() or provenance_out.exists():
        raise FileExistsError("refusing to overwrite diagnostic outputs")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(outputs[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(outputs)
    raw_values = np.array([float(row["lava_r_raw"]) for row in outputs])
    scaled_values = np.array([float(row["unit_diagonal_scaled_r_diagnostic_only"]) for row in outputs])
    diag_values = np.array([[float(row["factor_diagonal_snp1"]), float(row["factor_diagonal_snp2"])] for row in outputs])
    prov: dict[str, object] = {
        "analysis_id": "brain6-placo-ld-factor-normalization-v1",
        "interpretation": "Diagnostic-only reconstruction of LAVA format-14 factor products. Unit-diagonal scaling is a sensitivity diagnostic, not corrected LD; all original raw edges remain excluded and no candidate/locus result is promoted.",
        "sources": {"exception_table": {"path": str(EXCEPTIONS.relative_to(ROOT)), "sha256": sha256(EXCEPTIONS)},
                    "native_replay": {"path": str(REPLAY.relative_to(ROOT)), "sha256": sha256(REPLAY)},
                    "replay_provenance": {"path": str(REPLAY_PROVENANCE.relative_to(ROOT)), "sha256": sha256(REPLAY_PROVENANCE)},
                    "audit_script": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256(Path(__file__).resolve())},
                    "lava_load_ld_source": "https://github.com/josefin-werme/LAVA/blob/main/src/load_ld.cpp",
                    "lava_load_ld_source_scope": "Current upstream format-level implementation; pinned runtime replay remains authoritative for the historical run."},
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "method": {"format_id": 14, "stored_values": "float32 eigenvectors and eigenvalues", "reconstruction": "sum((v_i * eigenvalue) * v_j) in float64; compare with native LAVA read.ld output", "diagnostic_scaling": "raw / sqrt(factor_diagonal_i * factor_diagonal_j); not used downstream", "overwrite_prior_outputs": False},
        "reference_files": ref_hashes,
        "summary": {"unique_edges": len(outputs), "chromosomes": sorted({int(r["CHR"]) for r in outputs}),
                    "all_within_primary_blocks": True, "max_factor_vs_lava_abs_delta": max(float(r["factor_raw_delta"]) for r in outputs),
                    "raw_abs_r_min": float(np.min(np.abs(raw_values))), "raw_abs_r_max": float(np.max(np.abs(raw_values))),
                    "raw_abs_r_median": float(np.median(np.abs(raw_values))),
                    "factor_diagonal_min": float(np.min(diag_values)), "factor_diagonal_median": float(np.median(diag_values)),
                    "factor_diagonal_max": float(np.max(diag_values)),
                    "scaled_abs_r_min": float(np.min(np.abs(scaled_values))), "scaled_abs_r_max": float(np.max(np.abs(scaled_values))),
                    "scaled_abs_r_gt_1": int(np.sum(np.abs(scaled_values) > 1.0)),
                    "diagnostic_near_perfect_ld_threshold_abs_r_ge_0999_count": int(np.sum(np.abs(raw_values) >= 0.999)),
                    "promoted_edges": 0},
        "output": {"path": str(output.relative_to(ROOT)), "sha256": sha256(output), "rows": len(outputs)},
    }
    with provenance_out.open("x", encoding="utf-8") as stream:
        json.dump(prov, stream, indent=2, sort_keys=True); stream.write("\n")
    return prov


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE)
    args = parser.parse_args()
    args.reference = args.reference.resolve()
    args.output = args.output.resolve()
    args.provenance = args.provenance.resolve()
    print(json.dumps(run(args.reference, args.output, args.provenance)["summary"], indent=2))


if __name__ == "__main__":
    main()
