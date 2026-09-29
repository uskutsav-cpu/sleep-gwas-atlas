#!/usr/bin/env python3
"""Audit whether source-paper Q_SNP regions would alter latent-factor inputs.

The Foote et al. workflow uses noQ_postGWAS files for LDSC after removing
Bonferroni-significant Q_SNP regions (P <= 5e-8/7) plus a 1 Mb flank. This audit
checks the Q_metric_pvalue field in the seven GWAS Catalog exports without
rewriting or filtering any analysis input.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import TextIO


REPO = Path(__file__).resolve().parents[2]
FACTORS = [
    ("GCST90624046", "frailty_general"),
    ("GCST90624047", "frailty_factor_1"),
    ("GCST90624048", "frailty_factor_2"),
    ("GCST90624049", "frailty_factor_3"),
    ("GCST90624050", "frailty_factor_4"),
    ("GCST90624051", "frailty_factor_5"),
    ("GCST90624052", "frailty_factor_6"),
]
Q_THRESHOLD = 5e-8 / 7
FLANK_BP = 1_000_000
OUTPUT = REPO / "frailty_paper/analysis/latent_factor_qsnp_pruning_audit.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")


def is_q_significant(pvalue: float) -> bool:
    """Match the author pipeline's inclusive Q_SNP Bonferroni comparison."""
    return pvalue <= Q_THRESHOLD


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge_windows(points: list[tuple[str, int]]) -> dict[str, list[tuple[int, int]]]:
    windows: dict[str, list[tuple[int, int]]] = {}
    for chrom, position in points:
        chrom = chrom.removeprefix("chr")
        windows.setdefault(chrom, []).append(
            (max(1, position - FLANK_BP), position + FLANK_BP)
        )
    for chrom, values in windows.items():
        values.sort()
        merged: list[tuple[int, int]] = []
        for start, end in values:
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        windows[chrom] = merged
    return windows


def in_windows(chrom: str, position: int,
               windows: dict[str, list[tuple[int, int]]]) -> bool:
    # Small number of merged Q_SNP regions makes a direct scan clear and safe.
    chrom = chrom.removeprefix("chr")
    return any(start <= position <= end for start, end in windows.get(chrom, ()))


def audit_file(path: Path, accession: str, trait: str) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    required = {"chromosome", "base_pair_location", "Q_metric_pvalue"}
    n_rows = n_qp_valid = n_qp_invalid = 0
    points: list[tuple[str, int]] = []
    with open_text(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing columns {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            n_rows += 1
            try:
                pvalue = float(row["Q_metric_pvalue"])
                if not math.isfinite(pvalue) or not 0 <= pvalue <= 1:
                    raise ValueError
                n_qp_valid += 1
                if is_q_significant(pvalue):
                    chrom = row["chromosome"].removeprefix("chr")
                    position = int(row["base_pair_location"])
                    if not chrom or position < 1:
                        raise ValueError
                    points.append((chrom, position))
            except (TypeError, ValueError, KeyError):
                n_qp_invalid += 1
    windows = merge_windows(points)
    n_excluded = 0
    if windows:
        with open_text(path) as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            for row in reader:
                try:
                    chrom = row["chromosome"].removeprefix("chr")
                    position = int(row["base_pair_location"])
                    n_excluded += int(in_windows(chrom, position, windows))
                except (TypeError, ValueError, KeyError):
                    continue
    return {
        "accession": accession,
        "trait": trait,
        "source_file": str(path),
        "source_bytes": str(path.stat().st_size),
        "source_sha256": sha256(path),
        "rows_total": str(n_rows),
        "Q_metric_pvalue_valid": str(n_qp_valid),
        "Q_metric_pvalue_invalid_or_missing": str(n_qp_invalid),
        "threshold": f"P <= {Q_THRESHOLD:.12g} (5e-8/7)",
        "significant_Q_SNP_variants": str(len(points)),
        "merged_1Mb_windows": str(sum(map(len, windows.values()))),
        "rows_in_excluded_windows": str(n_excluded),
        "rows_retained_after_exclusion": str(n_rows - n_excluded),
        "result": "NO_SOURCE_ROWS_EXCLUDED" if n_excluded == 0 else "Q_SNP_REGIONS_PRESENT",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = args.output if args.output.is_absolute() else repo / args.output
    provenance = args.provenance if args.provenance.is_absolute() else repo / args.provenance
    output = output.resolve()
    provenance = provenance.resolve()
    rows = []
    for accession, trait in FACTORS:
        source = repo / "frailty_paper/data/gwas/latent_frailty_catalog" / accession / f"{accession}.tsv"
        rows.append(audit_file(source, accession, trait))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "schema_version": "latent_factor_qsnp_pruning_audit.v1",
        "status": "PASS" if all(r["result"] == "NO_SOURCE_ROWS_EXCLUDED" for r in rows) else "REGIONS_IDENTIFIED",
        "operation": "audit only; no source or analysis inputs modified",
        "command": "python frailty_paper/scripts/43_audit_latent_qsnp_pruning.py",
        "source_paper": "Foote et al. 2025, https://www.nature.com/articles/s41588-025-02269-0",
        "source_workflow": "https://github.com/IsyFoote/Frailty-Multivariate-GWAS/blob/main/12_LDSC_Frailty_latent_factors.R",
        "q_pruning_workflow": "https://github.com/IsyFoote/Frailty-Multivariate-GWAS/blob/main/6_Calculate_Q_and_Nhat.R",
        "archived_code_release": "https://doi.org/10.5281/zenodo.15654249",
        "q_threshold": Q_THRESHOLD,
        "q_threshold_comparison": "<=",
        "flank_bp_each_side": FLANK_BP,
        "rows": len(rows),
        "total_source_rows": sum(int(r["rows_total"]) for r in rows),
        "significant_q_variants": sum(int(r["significant_Q_SNP_variants"]) for r in rows),
        "excluded_rows": sum(int(r["rows_in_excluded_windows"]) for r in rows),
        "output": str(output),
        "output_sha256": sha256(output),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "python_version": __import__("platform").python_version(),
        "sources": [{"accession": r["accession"], "path": r["source_file"],
                     "sha256": r["source_sha256"], "rows": int(r["rows_total"])} for r in rows],
    }
    provenance.parent.mkdir(parents=True, exist_ok=True)
    provenance.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"LATENT_QSNP_AUDIT_OK factors={len(rows)} threshold={Q_THRESHOLD:.12g} output={output}")
    print(f"source_rows={sum(int(r['rows_total']) for r in rows)} "
          f"significant_q_variants={sum(int(r['significant_Q_SNP_variants']) for r in rows)} "
          f"excluded_rows={sum(int(r['rows_in_excluded_windows']) for r in rows)}")


if __name__ == "__main__":
    main()
