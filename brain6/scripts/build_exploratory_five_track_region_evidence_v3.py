#!/usr/bin/env python3
"""Complete the descriptive five-track 20-region ledger with current signed LD."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / "brain6/results/exploratory_five_track_v1"
V2 = ROOT / "brain6/results/exploratory_five_track_v2"
OUT = ROOT / "brain6/results/exploratory_five_track_v3"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check_receipt(folder: Path, receipt_name: str = "provenance.json") -> None:
    d = json.loads((folder / receipt_name).read_text())
    for name, expected in d["output_sha256"].items():
        if sha(folder / name) != expected:
            raise ValueError(f"Hash mismatch: {folder / name}")


def write(path: Path, data: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned exploratory v3 output exists")
    check_receipt(V1, "geographic_region_evidence.provenance.json")
    check_receipt(V1 / "b_signed_ld")
    check_receipt(V1 / "fourpair_signed_ld")
    check_receipt(V1 / "independent_cross_pair_ld")
    check_receipt(V2)
    previous = read(V2 / "geographic_region_evidence.tsv")
    b = read(V1 / "b_signed_ld/candidate_signed_ld.tsv")
    four = read(V1 / "fourpair_signed_ld/candidate_signed_ld.tsv")
    signed = b + four
    overlap = read(V1 / "independent_cross_pair_ld/overlap_edge_ld.tsv")
    unassigned = read(V1 / "unassigned_reference_unmatched_candidates.tsv")
    if len(previous) != 20 or len(signed) != 2686 or len(b) != 389 or len(four) != 2297 or len(unassigned) != 24 or len(overlap) != 5:
        raise ValueError("20/25/2,686/five-edge evidence coverage failed")
    if len({(r["pair_id"], r["SNP"]) for r in signed}) != 2686:
        raise ValueError("Signed-LD table has duplicate or missing candidate identities")
    counts = Counter(r["ld_qc_status"] for r in signed)
    if counts != {"PASS_SIGNED_LD": 2656, "OUT_OF_RANGE_RAW_R": 6, "NO_EXACT_REFERENCE_MATCH": 24}:
        raise ValueError(f"Unexpected signed-LD status accounting: {counts}")
    if any(r["ld_qc_status"] != "PASS_SIGNED_LD" for r in signed if r["candidate_status"] == "LEAD"):
        raise ValueError("One or more five-track leads lacks valid signed LD")
    if sum(r["candidate_status"] == "LEAD" for r in signed) != 27:
        raise ValueError("Expected 27 lead rows across 25 pair-specific loci")
    if {(r["pair_id"], r["SNP"]) for r in signed if r["ld_qc_status"] == "NO_EXACT_REFERENCE_MATCH"} != {
            (r["pair_id"], r["SNP"]) for r in unassigned}:
        raise ValueError("Reference-unmatched candidate identities changed")

    by_locus = defaultdict(list)
    for record in signed:
        if record["locus_id"] != "NA":
            by_locus[record["locus_id"]].append(record)
    updated = []
    used_edges = set()
    for item in previous:
        entry = dict(item)
        ids = set(item["candidate_locus_ids"].split(";"))
        members = [r for locus_id in ids for r in by_locus[locus_id]]
        if len(members) != int(item["n_candidate_variant_rows"]):
            raise ValueError(f"Region signed-LD candidate count mismatch: {item['geographic_region']}")
        cross = [r for r in overlap if r["locus_a"] in ids and r["locus_b"] in ids]
        used_edges.update((r["locus_a"], r["locus_b"]) for r in cross)
        entry["n_valid_signed_ld_rows"] = str(sum(r["ld_qc_status"] == "PASS_SIGNED_LD" for r in members))
        entry["n_out_of_range_raw_signed_r"] = str(sum(r["ld_qc_status"] == "OUT_OF_RANGE_RAW_R" for r in members))
        entry["n_effect_direction_opposite_to_signed_ld"] = str(sum(
            r["ld_qc_status"] == "PASS_SIGNED_LD" and
            (r["sleep_effect_vs_lead_ld"] == "OPPOSITE_TO_LD_SIGN" or
             (r.get("adhd_effect_vs_lead_ld") or r.get("disorder_effect_vs_lead_ld")) == "OPPOSITE_TO_LD_SIGN")
            for r in members))
        entry["n_independent_cross_pair_1000g_ld_edges"] = str(len(cross))
        entry["independent_cross_pair_1000g_r2"] = ";".join(
            r["lead_a"] + "|" + r["lead_b"] + "|" + r["independent_1000g_r2"] for r in cross
        ) or "NOT_APPLICABLE"
        entry["signed_ld_scope_current_candidates"] = "COMPLETE_EXPLORATORY_WITH_EXPLICIT_QC_EXCEPTIONS"
        updated.append(entry)
    if len(used_edges) != 5 or sum(int(r["n_valid_signed_ld_rows"]) for r in updated) != 2656 or sum(
            int(r["n_out_of_range_raw_signed_r"]) for r in updated) != 6:
        raise ValueError("Current signed-LD evidence did not reconcile to all regions")

    OUT.mkdir(parents=True, exist_ok=False)
    table = OUT / "geographic_region_evidence.tsv"
    write(table, updated)
    report = OUT / "region_evidence_report.md"
    lines = [
        "# Five-track PLACO exploratory 20-region evidence — current diagnostic ledger",
        "",
        "**EXPLORATORY. Zero LAVA-confirmed regions and zero final shared-locus tiers.**",
        "",
        "All 25 pair-specific candidates remain in the 20 coordinate-union",
        "groups below. Across 2,686 candidate rows, 2,656 have valid signed",
        "UKB-reference LD/effect orientation, six retain raw |r| > 1 as QC",
        "exceptions, and 24 lack an exact frozen-reference match. All 27 lead",
        "rows pass signed LD. Direction agreement is a post-selection",
        "descriptive check and cannot establish a shared causal effect.",
        "",
        "| Geographic region | Pairs | Leads | Candidate rows | Valid signed LD | Invalid raw r | Pair-level global replication | Positional genes (four-pair / B) |",
        "|---|---|---|---:|---:|---:|---|---|",
    ]
    for r in updated:
        lines.append(
            f"| {r['geographic_region']} | {r['pairs']} | {r['lead_variants']} | "
            f"{r['n_candidate_variant_rows']} | {r['n_valid_signed_ld_rows']} | "
            f"{r['n_out_of_range_raw_signed_r']} | {r['pair_level_global_replication']} | "
            f"{r['nearest_gene_context']} / {r['b_nearest_gene_positional_context']} |"
        )
    lines.extend([
        "",
        "The detailed table adds exact 1000G genotype coverage, five measured",
        "cross-pair overlap edges, source-specific GWAS Catalog context, and",
        "explicit empty evidence states. Those five independent-genotype LD",
        "comparisons have r² from 0.755 to 1.000; this corroborates LD among",
        "overlapping leads without establishing independent signals.",
        "The 24 reference-unmatched rows are preserved in the v1 unassigned",
        "table, outside the 20 coordinate groups.",
        "",
        "Only insomnia-ADHD has archived **global pair-level** directional",
        "replication in FinnGen; no region has locus-specific replication.",
        "Neither a GWAS Catalog exact-rsID record nor a nearby gene is",
        "independent biological support. No region has validated SuSiE,",
        "trait-trait coloc, eQTL/sQTL coloc, regulatory, tissue/cell-type, or",
        "pathway support. Thus no biologically supported region can be ranked",
        "as a validated shared locus. The multi-pair groups on chromosomes",
        "5, 7, 15, and 17 remain descriptive hypotheses, without a formal",
        "order or evidence-tier promotion.",
        "",
        "Canonical LAVA v3 remains `FAILED_QC_NOT_PROMOTED`: 3,720/17,465",
        "NOT_RUN versus a frozen maximum of 873. Local-rg claims, final",
        "region tiers, and confirmatory cross-layer claims are `BLOCKED_LAVA`.",
        "Any further molecular result in this branch stays exploratory and",
        "cannot enter the frozen tier decision.",
        "",
    ])
    report.write_text("\n".join(lines))
    sources = [
        V2 / "geographic_region_evidence.tsv", V2 / "provenance.json",
        V1 / "b_signed_ld/candidate_signed_ld.tsv", V1 / "b_signed_ld/provenance.json",
        V1 / "fourpair_signed_ld/candidate_signed_ld.tsv", V1 / "fourpair_signed_ld/provenance.json",
        V1 / "independent_cross_pair_ld/overlap_edge_ld.tsv",
        V1 / "independent_cross_pair_ld/provenance.json",
        V1 / "unassigned_reference_unmatched_candidates.tsv",
    ]
    receipt = {
        "analysis_id": "brain6_exploratory_five_track_region_evidence_v3",
        "analysis_label": "EXPLORATORY", "geographic_regions": 20,
        "pair_specific_candidates": 25, "candidate_variant_rows": 2686,
        "signed_ld_status_counts": dict(counts),
        "cross_pair_genotype_ld_edges": 5,
        "locus_specific_replication": 0,
        "fine_mapping_coloc_qtl_or_regulatory_support_established": 0,
        "lava_confirmed_regions": 0,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {table.name: sha(table), report.name: sha(report)},
        "interpretation": "Descriptive five-track PLACO/LD/annotation context only; all confirmatory tiers blocked by immutable LAVA QC failure.",
    }
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"regions": len(updated), "signed_ld": dict(counts),
                      "cross_pair_edges": len(used_edges)}, sort_keys=True))


if __name__ == "__main__":
    main()
