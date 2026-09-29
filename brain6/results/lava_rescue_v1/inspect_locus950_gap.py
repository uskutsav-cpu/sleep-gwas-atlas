#!/usr/bin/env python3
"""Read-only source-versus-reference check for the chr6 locus-950 MHC gap."""
from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "pilot_locus_950_reference_source_gap.json"
BASE = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1")
TRAITS = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
INTERVALS = {"left_949": (24950380, 25684629),
             "locus_950": (25684630, 26396200),
             "MHC_950_969": (25684630, 33864262),
             "right_970": (33864263, 34979270)}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ref = ROOT / "ref/lava/ukb_v1.1/lava-ukb-v1.1_chr6.info"
    loc = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
    ref_counts = {k: 0 for k in INTERVALS}
    with ref.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pos = int(row["POS"])
            for k, (start, stop) in INTERVALS.items():
                if start <= pos <= stop:
                    ref_counts[k] += 1
    traits = {}
    for trait in TRAITS:
        normalized = BASE / f"normalize_{trait}"
        db_path = normalized / "variants.sqlite"
        receipt_path = normalized / "receipt.json"
        receipt = json.loads(receipt_path.read_text())
        db = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        try:
            counts = {k: db.execute("SELECT COUNT(*) FROM variants WHERE chr=6 AND bp BETWEEN ? AND ?", (a, b)).fetchone()[0]
                      for k, (a, b) in INTERVALS.items()}
        finally:
            db.close()
        traits[trait] = {"normalized_variant_counts": counts,
                         "normalization_receipt_sha256": sha(receipt_path),
                         "normalization_exclude_regions": receipt["parameters"]["exclude_regions"],
                         "input_rows": receipt["qc"]["input_rows"],
                         "retained_rows": receipt["qc"]["retained_rows"]}
    result = {"analysis_id": "brain6-lava-rescue-v1-locus950-gap-check", "schema_version": 1,
              "reference_chr6_info_path": str(ref.relative_to(ROOT)), "reference_chr6_info_sha256": sha(ref),
              "locus_definition_path": str(loc.relative_to(ROOT)), "locus_definition_sha256": sha(loc),
              "intervals_grch37_inclusive": INTERVALS, "reference_variant_counts": ref_counts,
              "traits": traits,
              "interpretation": "The sealed reference contains 1907 variants in locus 950 while every pre-materialization normalized source has zero there and zero across loci 950-969. Each normalization receipt records no excluded regions and all input rows retained. This localizes the gap upstream of LAVA materialization, within the shared harmonized source inputs. It does not prove whether the upstream MHC removal was scientifically intentional; no source or threshold is changed.",
              "decision_effect": "At most the bounded min-K category is implicated; resolving all 156 such cells would still leave 3564 low-local-h2 NOT_RUN and fail the 873-cell cap."}
    payload = (json.dumps(result, sort_keys=True, indent=2) + "\n").encode()
    if OUT.exists():
        if OUT.read_bytes() != payload:
            raise FileExistsError("Refusing to rewrite gap receipt")
    else:
        OUT.write_bytes(payload)
    print(json.dumps({"reference": ref_counts, "trait_locus950": {t: v["normalized_variant_counts"]["locus_950"] for t, v in traits.items()}}, indent=2))


if __name__ == "__main__":
    main()
