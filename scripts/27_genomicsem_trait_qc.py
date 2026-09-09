#!/usr/bin/env python3
"""Create the locked SEM-input QC ledger from the multivariable LDSC diagonal."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


Z_MIN = 4.0
INTERCEPT_MAX = 1.20


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", default="config/analysis_panel.tsv")
    parser.add_argument("--phase1-h2", default="results/tables/h2_summary.tsv")
    parser.add_argument("--covariance-pairs", default="results/tables/ldsc_covariance_pairs.tsv")
    parser.add_argument("--out", default="results/tables/genomicsem_trait_inclusion.tsv")
    args = parser.parse_args()

    panel = read_tsv(Path(args.panel))
    if len(panel) != 45 or len({row["trait_id"] for row in panel}) != 45:
        raise SystemExit("ERROR: expected the locked 45-trait panel")
    phase1 = {row["trait"]: row for row in read_tsv(Path(args.phase1_h2))}
    pairs = read_tsv(Path(args.covariance_pairs))
    diagonals = {
        row["trait_1"]: row
        for row in pairs
        if row["trait_1"] == row["trait_2"]
    }
    expected = {row["trait_id"] for row in panel}
    if set(phase1) != expected or set(diagonals) != expected:
        raise SystemExit("ERROR: Phase-1 h2 or covariance diagonal differs from the locked panel")

    rows = []
    for order, trait in enumerate(panel, start=1):
        trait_id = trait["trait_id"]
        old = phase1[trait_id]
        new = diagonals[trait_id]
        h2 = as_float(new["genetic_covariance"])
        se = as_float(new["genetic_covariance_se"])
        z = h2 / se if se > 0 else math.nan
        intercept = as_float(new["cross_trait_intercept"])
        pass_z = math.isfinite(z) and z >= Z_MIN
        pass_intercept = math.isfinite(intercept) and intercept <= INTERCEPT_MAX
        include = pass_z and pass_intercept
        if not math.isfinite(z):
            reason = "h2_or_se_missing"
        elif not pass_z:
            reason = f"h2_z_below_{Z_MIN:g}"
        elif not math.isfinite(intercept):
            reason = "intercept_missing"
        elif not pass_intercept:
            reason = f"intercept_above_{INTERCEPT_MAX:g}"
        else:
            reason = "pass"
        warnings = []
        if trait["type"] == "binary" and h2 > 1:
            warnings.append("liability_h2_above_1")
        phase1_pass = old["verdict"] == "PASS"
        if phase1_pass != include:
            warnings.append("phase1_vs_genomicsem_implementation_sensitivity")
        rows.append({
            "order": order,
            "trait_id": trait_id,
            "domain": trait["domain"],
            "type": trait["type"],
            "scale": new["scale_1"],
            "phase1_h2": old["h2"],
            "phase1_h2_se": old["se"],
            "phase1_intercept": old["intercept"],
            "phase1_verdict": old["verdict"],
            "genomicsem_h2": new["genetic_covariance"],
            "genomicsem_h2_se": new["genetic_covariance_se"],
            "genomicsem_h2_z": f"{z:.12g}",
            "genomicsem_intercept": new["cross_trait_intercept"],
            "pass_h2_z_ge_4": str(pass_z).upper(),
            "pass_intercept_le_1_2": str(pass_intercept).upper(),
            "include_genomic_sem": str(include).upper(),
            "qc_reason": reason,
            "warning": ";".join(warnings) if warnings else "none",
        })

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    included = sum(row["include_genomic_sem"] == "TRUE" for row in rows)
    excluded = [row["trait_id"] for row in rows if row["include_genomic_sem"] != "TRUE"]
    print(f"Published GenomicSEM input QC: {included}/45 included; excluded: {', '.join(excluded)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
