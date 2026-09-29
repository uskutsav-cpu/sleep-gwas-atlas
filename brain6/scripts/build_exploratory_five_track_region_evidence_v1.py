#!/usr/bin/env python3
"""Join existing five-track candidate evidence without changing any promotion gate."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/exploratory_five_track_v1"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, items: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(items)


def verified_outputs(folder: Path) -> None:
    receipt = json.loads((folder / "provenance.json").read_text())
    expected = receipt.get("output_sha256", receipt.get("outputs_sha256", {}))
    for name, digest in expected.items():
        if sha(folder / name) != digest:
            raise ValueError(f"Source receipt hash mismatch: {folder / name}")


def join_unique(values: list[str]) -> str:
    return ";".join(sorted(set(v for v in values if v and v != "NA")))


def main() -> None:
    base = ROOT / "brain6/results"
    locus_root = base / "loci"
    reconciliation = locus_root / "five_track_cross_pair_region_reconciliation_v1"
    candidates = locus_root / "placo_five_track_candidate_loci_v1"
    phase3 = locus_root / "independent_1000g_phase3_eur_ld_v1"
    admitted_b = locus_root / "admitted_B_1000g_ld_v1"
    catalog = locus_root / "gwas_catalog_candidate_leads_v1"
    genes = base / "annotation/partial_candidate_lead_nearest_gene_grch37_v1.tsv"
    signed = locus_root / "placo_candidate_signed_effect_ld_partial_v1.tsv"
    replication = base / "replication/replication_master.tsv"
    lava_decision = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json"
    if sha(lava_decision) != "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22":
        raise ValueError("Immutable canonical LAVA decision changed")

    for folder in (reconciliation, candidates, phase3, admitted_b, catalog):
        verified_outputs(folder)
    gene_receipt = json.loads((genes.with_suffix(".provenance.json")).read_text())
    if sha(genes) != gene_receipt["output"]["sha256"]:
        raise ValueError("Nearest-gene source hash mismatch")
    signed_receipt = json.loads((signed.with_suffix(".provenance.json")).read_text())
    if sha(signed) != signed_receipt["outputs"][0]["sha256"]:
        raise ValueError("Signed-LD source hash mismatch")
    candidate_receipt = json.loads((candidates / "provenance.json").read_text())
    for key, name in (("candidate_loci_sha256", "candidate_loci.tsv"),
                      ("candidate_variants_sha256", "candidate_variants.tsv")):
        if sha(candidates / name) != candidate_receipt[key]:
            raise ValueError(f"Five-track source hash mismatch: {name}")

    groups = rows(reconciliation / "geographic_region_groups.tsv")
    members = rows(reconciliation / "pair_candidate_members.tsv")
    loci = rows(candidates / "candidate_loci.tsv")
    variants = rows(candidates / "candidate_variants.tsv")
    phase3_rows = {r["locus_id"]: r for r in rows(phase3 / "region_summary.tsv")}
    overlap_edges = rows(reconciliation / "cross_pair_overlap_edges.tsv")
    b_matches = {r["SNP"]: r for r in rows(admitted_b / "variant_matches.tsv")}
    catalog_rows = defaultdict(list)
    for record in rows(catalog / "trait_relevant_context.tsv"):
        catalog_rows[record["locus_id"]].append(record)
    gene_rows = defaultdict(list)
    for record in rows(genes):
        gene_rows[record["candidate_region_id"]].append(record)
    signed_rows = defaultdict(list)
    for record in rows(signed):
        signed_rows[record["locus_id"]].append(record)
    replication_rows = {
        r["sleep_trait"] + "__" + r["brain_disorder"]: r
        for r in rows(replication) if r["sleep_trait"] != "NA"
    }

    if len(groups) != 20 or len(members) != 25 or len(loci) != 25 or len(variants) != 2686:
        raise ValueError("Five-track 20/25/2,686 coverage invariant failed")
    by_locus = {r["locus_id"]: r for r in loci}
    if len(by_locus) != 25 or {r["candidate_locus_id"] for r in members} != set(by_locus):
        raise ValueError("Pair-specific candidate membership mismatch")
    by_group = defaultdict(list)
    for record in members:
        by_group[record["region_group"]].append(record)
    if set(by_group) != {r["region_group"] for r in groups}:
        raise ValueError("Geographic group membership mismatch")
    by_variants = defaultdict(list)
    unassigned: list[dict[str, str]] = []
    for record in variants:
        if record["locus_id"] not in by_locus:
            if (record["locus_id"] == "NA" and
                    record["candidate_status"] == "NOT_CLUMPED_NO_EXACT_REFERENCE_MATCH"):
                unassigned.append(record)
                continue
            raise ValueError("Candidate variant has unknown locus")
        by_variants[record["locus_id"]].append(record)
    if len(unassigned) != 24 or sum(len(v) for v in by_variants.values()) != 2662:
        raise ValueError("Candidate variant accounting mismatch")

    per_locus: dict[str, dict[str, object]] = {}
    for locus_id, locus in by_locus.items():
        pair = locus["pair_id"]
        variant_rows = by_variants[locus_id]
        if len(variant_rows) != int(locus["n_candidate_variants"]):
            raise ValueError(f"Candidate count mismatch: {locus_id}")
        if pair == "insomnia__adhd":
            exact = sum(b_matches[v["SNP"]]["status"] == "EXACT_ALLELE_GENOTYPE"
                        for v in variant_rows)
            ld_measured = "NOT_RECORDED_FOR_B_ASSIGNMENTS"
            ld_discordant = "NOT_RECORDED_FOR_B_ASSIGNMENTS"
            ld_panel = "1000G EUR Phase3/30x diagnostic; six lead identities reproduce"
            signed_status = "NOT_RUN_FOR_B"
            signed_valid = 0
            signed_inconsistent = 0
        else:
            ld = phase3_rows[locus_id]
            exact = int(ld["n_exact_vcf"])
            ld_measured = ld["n_edges_measured"]
            ld_discordant = ld["n_r2_threshold_discordant"]
            ld_panel = "1000G Phase3 EUR diagnostic; see targeted 30x sensitivity"
            signal_rows = signed_rows[locus_id]
            if ({x["SNP"] for x in signal_rows} == {x["SNP"] for x in variant_rows} and
                    all(x["lead_SNP"] == y["lead_SNP"] for x, y in zip(
                        sorted(signal_rows, key=lambda z: z["SNP"]),
                        sorted(variant_rows, key=lambda z: z["SNP"])))):
                signed_valid = sum(x["ld_qc_status"] == "PASS_VALID_SIGNED_R" for x in signal_rows)
                signed_inconsistent = sum(
                    x["ld_qc_status"] == "PASS_VALID_SIGNED_R" and
                    (x["sleep_effect_vs_lead_ld"] != "CONSISTENT_WITH_LD_SIGN" or
                     x["disorder_effect_vs_lead_ld"] != "CONSISTENT_WITH_LD_SIGN")
                    for x in signal_rows
                )
                signed_status = "DESCRIPTIVE_UKB_SIGNED_LD_SAME_CANDIDATE_ASSIGNMENTS"
            else:
                signed_valid = 0
                signed_inconsistent = 0
                signed_status = "NOT_REUSABLE_DIFFERENT_CLUMP_VERSION"
        per_locus[locus_id] = {
            "pair_id": pair, "candidate_variants": len(variant_rows),
            "exact_1000g_allele_genotypes": exact,
            "independent_ld_measured_edges": ld_measured,
            "independent_ld_r2_threshold_discordances": ld_discordant,
            "independent_ld_panel": ld_panel,
            "signed_ld_status": signed_status,
            "signed_ld_valid_rows": signed_valid,
            "signed_ld_inconsistent_rows": signed_inconsistent,
        }

    output: list[dict[str, object]] = []
    for group in groups:
        region = group["region_group"]
        group_members = by_group[region]
        ids = [m["candidate_locus_id"] for m in group_members]
        pairs = [m["pair_id"] for m in group_members]
        if len(ids) != int(group["n_pair_specific_candidates"]):
            raise ValueError(f"Geographic group count mismatch: {region}")
        catalog_here = [r for locus_id in ids for r in catalog_rows[locus_id]]
        genes_here = [r for locus_id in ids for r in gene_rows[locus_id]]
        global_rep = [p for p in pairs if replication_rows[p]["replication_class"] ==
                      "DIRECTIONAL_REPLICATION"]
        statuses = [str(per_locus[x]["signed_ld_status"]) for x in ids]
        output.append({
            "geographic_region": region,
            "chromosome": group["chromosome"], "start": group["start"], "stop": group["stop"],
            "n_pair_specific_candidates": len(ids), "pairs": join_unique(pairs),
            "candidate_locus_ids": ";".join(ids),
            "lead_variants": group["lead_variants"],
            "best_lead_p_placo": format(min(float(m["lead_p_placo"]) for m in group_members), ".12g"),
            "n_candidate_variant_rows": sum(int(per_locus[x]["candidate_variants"]) for x in ids),
            "n_exact_1000g_allele_genotypes": sum(int(per_locus[x]["exact_1000g_allele_genotypes"]) for x in ids),
            "n_independent_ld_measured_edges_four_pairs": sum(
                int(per_locus[x]["independent_ld_measured_edges"]) for x in ids
                if str(per_locus[x]["independent_ld_measured_edges"]).isdigit()),
            "n_independent_ld_r2_threshold_discordances_four_pairs": sum(
                int(per_locus[x]["independent_ld_r2_threshold_discordances"]) for x in ids
                if str(per_locus[x]["independent_ld_r2_threshold_discordances"]).isdigit()),
            "overlap_ukb_lead_r2": ";".join(
                e["observed_frozen_ukb_lead_r2"] for e in overlap_edges
                if e["locus_a"] in ids and e["locus_b"] in ids) or "NOT_APPLICABLE",
            "independent_ld_scope": join_unique([str(per_locus[x]["independent_ld_panel"]) for x in ids]),
            "signed_ld_scope": join_unique(statuses),
            "signed_ld_valid_rows_four_pairs": sum(int(per_locus[x]["signed_ld_valid_rows"]) for x in ids),
            "signed_ld_inconsistent_rows_four_pairs": sum(int(per_locus[x]["signed_ld_inconsistent_rows"]) for x in ids),
            "pair_level_global_replication": join_unique(global_rep) or "NONE",
            "locus_specific_replication": "NOT_ESTABLISHED",
            "prior_catalog_trait_relevant_records": len(catalog_here),
            "prior_catalog_studies": join_unique([r["study_accession"] for r in catalog_here]) or "NONE_IN_EXISTING_SNAPSHOT",
            "nearest_gene_context": join_unique([r["gene_symbol"] for r in genes_here]) or "NOT_ANNOTATED",
            "fine_mapping": "NOT_RUN_METHOD_PREREQUISITES_UNREVIEWED",
            "trait_trait_coloc": "NOT_RUN_METHOD_PREREQUISITES_UNREVIEWED",
            "eqtl_sqtl_coloc": "NOT_RUN_NO_VERIFIED_QTL_INPUT",
            "regulatory_support": "NOT_ESTABLISHED",
            "tissue_cell_type_support": "NOT_ESTABLISHED",
            "pathway_support": "NOT_ESTABLISHED",
            "lava_local_rg": "BLOCKED_LAVA",
            "final_region_tier": "BLOCKED_LAVA",
            "analysis_label": "EXPLORATORY",
        })
    output.sort(key=lambda r: (int(r["chromosome"]), int(r["start"])))
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "geographic_region_evidence.tsv"
    if target.exists():
        raise FileExistsError("Exploratory output already exists; create a new version")
    write_rows(target, output, list(output[0]))
    unassigned_target = OUT / "unassigned_reference_unmatched_candidates.tsv"
    write_rows(unassigned_target, unassigned, list(unassigned[0]))
    report_target = OUT / "region_evidence_report.md"
    report_lines = [
        "# Five-track PLACO exploratory region evidence",
        "",
        "**EXPLORATORY; no LAVA-confirmed region or final shared-locus tier.**",
        "",
        "The five-track PLACO family has 25 pair-specific candidate intervals. Their",
        "coordinate union comprises 20 geographic groups; neither count is an",
        "independent causal-signal count. Of 2,686 candidate variant rows, 2,662",
        "are assigned to these intervals and 24 lack an exact frozen-reference",
        "match and remain in `unassigned_reference_unmatched_candidates.tsv`.",
        "",
        "| Geographic region | Pairs | Leads | Candidate rows | Exact 1000G genotypes | Prior catalog records | Positional gene context |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for item in output:
        report_lines.append(
            f"| {item['geographic_region']} | {item['pairs']} | {item['lead_variants']} | "
            f"{item['n_candidate_variant_rows']} | {item['n_exact_1000g_allele_genotypes']} | "
            f"{item['prior_catalog_trait_relevant_records']} | {item['nearest_gene_context']} |"
        )
    report_lines.extend([
        "",
        "`geographic_region_evidence.tsv` records the full pair membership,",
        "candidate count, diagnostic genotype/LD coverage, signed-LD reuse status,",
        "pair-level global replication, prior GWAS Catalog study IDs, and the",
        "method-specific empty evidence states for every region. The six",
        "insomnia-ADHD candidate regions inherit a **pair-level** FinnGen global",
        "directional replication label; no region has locus-specific replication.",
        "The existing signed-LD diagnostic applies only where its exact",
        "four-pair candidate-to-lead assignments match this five-track set.",
        "It is not reusable for three loci with changed clumps or any of the",
        "six admitted B loci. Independent 1000G LD is diagnostic and does not",
        "replace the frozen UKB locus definition.",
        "",
        "No region has validated SuSiE, trait-trait coloc, eQTL/sQTL coloc,",
        "regulatory, tissue/cell-type, or pathway support in this ledger.",
        "The nearest-gene names are descriptive positions. Curated catalog",
        "matches may reuse discovery cohorts and do not establish independent",
        "replication. Accordingly, there is **no biologically supported region",
        "that can be ranked as a validated shared locus**. The four multi-pair",
        "groups on chromosomes 5, 7, 15, and 17 are descriptive overlap",
        "hypotheses, retained with all of their pair-specific members.",
        "",
        "The canonical LAVA v3 family remains `FAILED_QC_NOT_PROMOTED`",
        "(3,720/17,465 NOT_RUN, frozen maximum 873). All local-rg claims,",
        "confirmatory region tiers, and dependent cross-layer categories are",
        "`BLOCKED_LAVA`. Biological annotations in this exploratory branch",
        "cannot enter the frozen tier rule.",
        "",
    ])
    report_target.write_text("\n".join(report_lines))
    sources = [
        reconciliation / "geographic_region_groups.tsv",
        reconciliation / "pair_candidate_members.tsv",
        reconciliation / "cross_pair_overlap_edges.tsv",
        reconciliation / "provenance.json",
        candidates / "candidate_loci.tsv", candidates / "candidate_variants.tsv",
        candidates / "provenance.json", phase3 / "region_summary.tsv", phase3 / "provenance.json",
        admitted_b / "variant_matches.tsv", admitted_b / "provenance.json",
        catalog / "trait_relevant_context.tsv", catalog / "provenance.json",
        genes, signed, replication, lava_decision,
    ]
    receipt = {
        "analysis_id": "brain6_exploratory_five_track_region_evidence_v1",
        "analysis_label": "EXPLORATORY",
        "geographic_regions": len(output), "pair_specific_candidates": len(members),
        "candidate_variant_rows": len(variants),
        "assigned_candidate_variant_rows": 2662,
        "unassigned_reference_unmatched_candidate_rows": len(unassigned),
        "locus_specific_replication_established": 0,
        "fine_mapping_coloc_qtl_or_regulatory_support_established": 0,
        "lava_confirmed_regions": 0,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {target.name: sha(target), unassigned_target.name: sha(unassigned_target),
                          report_target.name: sha(report_target)},
        "limitations": [
            "Twenty geographic groups are coordinate unions, not independent genetic signals.",
            "FinnGen ADHD evidence is pair-level global directional replication, not locus replication.",
            "GWAS Catalog is a dated curated association snapshot, not independent replication.",
            "Nearest gene is positional context, not gene prioritization.",
            "Partial signed-LD analysis covers the four original pairs; B is pending signed orientation.",
            "All final tiers remain blocked by immutable canonical LAVA family QC failure.",
        ],
    }
    (OUT / "geographic_region_evidence.provenance.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    counts = Counter(int(r["n_pair_specific_candidates"]) for r in output)
    print(json.dumps({"regions": len(output), "candidates": len(members),
                      "candidate_variants": len(variants), "unassigned": len(unassigned),
                      "group_sizes": dict(counts)}, sort_keys=True))


if __name__ == "__main__":
    main()
