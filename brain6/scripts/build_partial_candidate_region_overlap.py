#!/usr/bin/env python3
"""Describe cross-pair overlap among gated partial PLACO candidate intervals.

Interval overlap and repeated lead IDs are diagnostics only. They do not
establish a shared causal signal or an independent locus.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "brain6/results/annotation/partial_candidate_region_gate_v1.tsv"
GATE_PROVENANCE = ROOT / "brain6/results/annotation/partial_candidate_region_gate_v1.provenance.json"
SENSITIVITY = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/provenance.json"
OUT = ROOT / "brain6/results/annotation/partial_candidate_region_overlap_v1.tsv"
PROVENANCE = ROOT / "brain6/results/annotation/partial_candidate_region_overlap_v1.provenance.json"
FIELDS = [
    "pair_a", "candidate_region_a", "pair_b", "candidate_region_b", "CHR",
    "start_a", "stop_a", "start_b", "stop_b", "overlap_start", "overlap_stop",
    "overlap_bp", "lead_variants_a", "lead_variants_b", "exact_shared_lead_variants",
    "relationship_class", "evidence_tier", "interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build(output: Path = OUT, provenance_path: Path = PROVENANCE) -> dict[str, Any]:
    gate_prov = json.loads(GATE_PROVENANCE.read_text(encoding="utf-8"))
    sensitivity = json.loads(SENSITIVITY.read_text(encoding="utf-8"))
    if (gate_prov.get("status") != "PASS_NOT_TIERED_PREREQUISITES_BLOCKED" or
            gate_prov.get("output", {}).get("sha256") != sha256(GATE) or
            sensitivity.get("analysis_id") != "brain6_placo_factor_normalized_sensitivity_v2" or
            sensitivity.get("scope", {}).get("full_five_pair_family") is not False):
        raise ValueError("Expected hash-valid NOT_TIERED gate on the partial sensitivity regions")
    regions = read_tsv(GATE)
    if len(regions) != 19 or any(row["evidence_tier"] != "NOT_TIERED" for row in regions):
        raise ValueError("Expected all 19 partial candidate regions to remain NOT_TIERED")
    rows: list[dict[str, str]] = []
    cross_pair_comparisons = 0
    for a, b in itertools.combinations(regions, 2):
        if a["pair_id"] == b["pair_id"] or a["CHR"] != b["CHR"]:
            continue
        cross_pair_comparisons += 1
        overlap_start = max(int(a["START"]), int(b["START"]))
        overlap_stop = min(int(a["STOP"]), int(b["STOP"]))
        if overlap_start > overlap_stop:
            continue
        leads_a, leads_b = set(a["lead_variants"].split(";")), set(b["lead_variants"].split(";"))
        shared = sorted(leads_a & leads_b)
        rows.append({
            "pair_a": a["pair_id"], "candidate_region_a": a["candidate_region_id"],
            "pair_b": b["pair_id"], "candidate_region_b": b["candidate_region_id"],
            "CHR": a["CHR"], "start_a": a["START"], "stop_a": a["STOP"],
            "start_b": b["START"], "stop_b": b["STOP"],
            "overlap_start": str(overlap_start), "overlap_stop": str(overlap_stop),
            "overlap_bp": str(overlap_stop - overlap_start + 1),
            "lead_variants_a": a["lead_variants"], "lead_variants_b": b["lead_variants"],
            "exact_shared_lead_variants": ";".join(shared) if shared else "NONE",
            "relationship_class": "EXACT_LEAD_REUSE_AND_INTERVAL_OVERLAP" if shared else "INTERVAL_OVERLAP_ONLY",
            "evidence_tier": "NOT_TIERED",
            "interpretation": "Coordinate diagnostic only; no shared causal signal or independent-locus claim",
        })
    rows.sort(key=lambda r: (int(r["CHR"]), int(r["overlap_start"]), r["pair_a"], r["pair_b"]))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    repeated_leads = sum(row["exact_shared_lead_variants"] != "NONE" for row in rows)
    provenance: dict[str, Any] = {
        "schema_version": 1,
        "analysis_id": "brain6_partial_placo_cross_pair_region_overlap_v1",
        "status": "PASS_DESCRIPTIVE_NOT_SHARED_LOCUS_EVIDENCE",
        "scope": {
            "candidate_regions": len(regions), "pairs": sorted({r["pair_id"] for r in regions}),
            "cross_pair_region_comparisons_same_chromosome": cross_pair_comparisons,
            "overlapping_region_pairs": len(rows), "pairs_with_exact_shared_lead_variant": repeated_leads,
            "exact_shared_lead_variant_ids": sorted({
                lead for row in rows if row["exact_shared_lead_variants"] != "NONE"
                for lead in row["exact_shared_lead_variants"].split(";")
            }),
            "all_regions_remain_not_tiered": True,
            "shared_locus_claims": False,
        },
        "method": "Inclusive genomic interval overlap across different PLACO pair outputs; exact overlap of lead-variant IDs tracked separately.",
        "limitations": [
            "Source is a four-pair diagnostic PLACO family; protected Track B is absent.",
            "Intervals are 500-kb lead-flank regions after factor-normalized reference-LD grouping, not genotype-validated independent loci.",
            "Shared lead IDs or broad coordinate overlap do not demonstrate shared causal variants, colocalization, or independent replication.",
        ],
        "sources": {
            "candidate_region_gate": {"path": str(GATE.relative_to(ROOT)), "sha256": sha256(GATE)},
            "candidate_region_gate_provenance": {
                "path": str(GATE_PROVENANCE.relative_to(ROOT)), "sha256": sha256(GATE_PROVENANCE),
            },
            "factor_normalized_sensitivity_provenance": {
                "path": str(SENSITIVITY.relative_to(ROOT)), "sha256": sha256(SENSITIVITY),
            },
        },
        "builder": {"path": str(Path(__file__).resolve().relative_to(ROOT)),
                    "sha256": sha256(Path(__file__).resolve())},
        "output": {"path": str(output.relative_to(ROOT)), "rows": len(rows), "sha256": sha256(output)},
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.provenance), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
