#!/usr/bin/env python3
"""Reconcile five geographic overlap edges with existing 1000G genotype LD."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCI = ROOT / "brain6/results/loci"
OUT = ROOT / "brain6/results/exploratory_five_track_v1/independent_cross_pair_ld"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned LD output exists; use a new analysis version")
    reconciliation = LOCI / "five_track_cross_pair_region_reconciliation_v1"
    b = LOCI / "admitted_B_1000g_ld_v1"
    other = LOCI / "independent_1000g_phase3_eur_ld_v1"
    edges = read(reconciliation / "cross_pair_overlap_edges.tsv")
    b_assignments = {x["SNP"]: x for x in read(b / "variant_assignments.tsv")}
    b_matches = {x["SNP"]: x for x in read(b / "variant_matches.tsv")}
    external = read(other / "cross_region_lead_ld.tsv")
    if len(edges) != 5:
        raise ValueError("Expected five frozen geographic overlap edges")
    result = []
    for edge in edges:
        lead_a, lead_b, ukb_r2 = edge["observed_frozen_ukb_lead_r2"].split("|")
        if edge["pair_a"] == "insomnia__adhd":
            candidate = b_assignments[lead_b]
            if (candidate["lead_SNP"] != lead_a or
                    b_matches[lead_a]["status"] != "EXACT_ALLELE_GENOTYPE" or
                    b_matches[lead_b]["status"] != "EXACT_ALLELE_GENOTYPE"):
                raise ValueError("B cross-pair lead does not match exact genotype assignment")
            independent_r2 = float(candidate["r2_to_lead"])
            if b_matches[lead_a]["reference"] != b_matches[lead_b]["reference"]:
                raise ValueError("B cross-pair leads came from different genotype panels")
            reference = b_matches[lead_a]["reference"] + "_EUR_503"
        else:
            hits = [x for x in external if
                    {x["locus_a"], x["locus_b"]} == {edge["locus_a"], edge["locus_b"]} and
                    {x["lead_a"], x["lead_b"]} == {lead_a, lead_b}]
            if len(hits) != 1 or hits[0]["both_vcf_match"] != "1" or hits[0]["eur_n"] != "503":
                raise ValueError("Non-B cross-pair lead lacks a unique exact-genotype LD record")
            independent_r2 = float(hits[0]["eur_r2"])
            reference = "1000G_PHASE3_EUR_503"
        result.append({
            "locus_a": edge["locus_a"], "locus_b": edge["locus_b"],
            "pair_a": edge["pair_a"], "pair_b": edge["pair_b"],
            "lead_a": lead_a, "lead_b": lead_b,
            "frozen_ukb_r2": ukb_r2,
            "independent_1000g_r2": format(independent_r2, ".12g"),
            "absolute_r2_difference": format(abs(independent_r2 - float(ukb_r2)), ".12g"),
            "independent_genotype_reference": reference,
            "analysis_label": "EXPLORATORY",
            "interpretation": "GEOGRAPHIC_OVERLAP_WITH_DIAGNOSTIC_LD_NOT_INDEPENDENT_SIGNAL_COUNT",
            "lava_local_rg": "BLOCKED_LAVA",
        })
    OUT.mkdir(parents=True, exist_ok=False)
    target = OUT / "overlap_edge_ld.tsv"
    with target.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(result[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    receipt = {
        "analysis_id": "brain6_exploratory_cross_pair_ld_v1",
        "analysis_label": "EXPLORATORY", "overlap_edges": len(result),
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in (
            reconciliation / "cross_pair_overlap_edges.tsv", reconciliation / "provenance.json",
            b / "variant_assignments.tsv", b / "variant_matches.tsv", b / "provenance.json",
            other / "cross_region_lead_ld.tsv", other / "provenance.json")},
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {target.name: sha(target)},
        "limitations": [
            "The high-coverage and Phase 3 calls reuse the same 503 EUR individuals, so they are not independent replication cohorts.",
            "R-squared is an orientation-free LD diagnostic; signed direction is a separate analysis.",
            "No result changes frozen UKB clumping or establishes independent causal signals or LAVA local rg.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"overlap_edges": len(result),
                      "min_independent_r2": min(float(x["independent_1000g_r2"]) for x in result)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
