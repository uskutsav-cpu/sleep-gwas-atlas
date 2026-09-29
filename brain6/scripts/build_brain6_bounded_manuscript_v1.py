#!/usr/bin/env python3
"""Reconcile frozen Brain6 candidates into bounded, source-bound manuscript tables.

This is a deterministic descriptive join. It never changes a protected candidate,
canonical LAVA receipt, or biological tier, and never estimates missing methods.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results"
OUT = BASE / "brain6_bounded_manuscript_v1"
INPUTS = {
    "regions": BASE / "exploratory_five_track_v6/geographic_region_evidence.tsv",
    "members": BASE / "loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv",
    "canonical": BASE / "lava_longsleep_source_rescue_v1/BRAIN6_FINAL_REGION_EVIDENCE.tsv",
    "yale": BASE / "lava_longsleep_source_rescue_v1/exploratory_long_sleep_replication_v1/pair_candidate_replication.tsv",
    "qtl": BASE / "exploratory_five_track_v5/gtex_brain_qtl_credible_sets/candidate_molecular_credible_set_memberships.tsv",
    "nominal_eqtl": BASE / "exploratory_five_track_v4/lead_eqtl_context/lead_eqtl_associations.tsv",
    "regulatory": BASE / "exploratory_five_track_v6/lead_regulatory_context/lead_regulatory_features.tsv",
    "nearest": BASE / "exploratory_five_track_v1/b_annotation/nearest_gene.tsv",
    "catalog": BASE / "loci/gwas_catalog_candidate_leads_v1/trait_relevant_context.tsv",
    "global_replication": BASE / "replication/replication_master.tsv",
    "method_readiness": BASE / "exploratory_five_track_v3/method_readiness_20260927.json",
    "canonical_decision": ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(key: str) -> list[dict[str, str]]:
    with INPUTS[key].open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def stable_write(name: str, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Empty table: {name}")
    columns = list(rows[0])
    from io import StringIO

    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path = OUT / name
    content = buffer.getvalue()
    if path.exists():
        if path.read_text() != content:
            raise FileExistsError(f"Refusing to overwrite changed output: {path}")
    else:
        path.write_text(content)


def separated(items: set[str] | list[str]) -> str:
    return ";".join(sorted(x for x in items if x)) or "NONE"


def main() -> None:
    for key, path in INPUTS.items():
        if not path.is_file():
            raise FileNotFoundError(f"{key}: {path}")
    regions = read("regions")
    members = read("members")
    canonical = read("canonical")
    yale = read("yale")
    qtl = read("qtl")
    nominal = read("nominal_eqtl")
    regulatory = read("regulatory")
    nearest = read("nearest")
    catalog = read("catalog")
    global_replication = read("global_replication")
    if len(regions) != 20 or len(members) != 25 or len(canonical) != 25 or len(yale) != 25:
        raise ValueError("Frozen 20-region/25-pair cardinality changed")
    region_by_id = {r["geographic_region"]: r for r in regions}
    if len(region_by_id) != 20:
        raise ValueError("Duplicate geographic region")
    member_by_id = {r["candidate_locus_id"]: r for r in members}
    if len(member_by_id) != 25:
        raise ValueError("Duplicate pair-specific candidate")
    canonical_by_id = {r["candidate_locus_id"]: r for r in canonical}
    yale_by_id = {r["candidate_locus_id"]: r for r in yale}
    if set(member_by_id) != set(canonical_by_id) or set(member_by_id) != set(yale_by_id):
        raise ValueError("Candidate IDs differ across frozen/current tables")
    qtl_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    nominal_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    reg_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    near_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    catalog_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for source, target, key in (
        (qtl, qtl_by_id, "locus_id"), (nominal, nominal_by_id, "locus_id"),
        (regulatory, reg_by_id, "locus_id"), (nearest, near_by_id, "locus_id"),
        (catalog, catalog_by_id, "locus_id"),
    ):
        for row in source:
            if row[key] not in member_by_id:
                raise ValueError(f"Unknown candidate in {key}: {row[key]}")
            target[row[key]].append(row)
    global_by_pair = {r["sleep_trait"] + "__" + r["brain_disorder"]: r for r in global_replication}
    if len(global_by_pair) != len(global_replication):
        raise ValueError("Duplicate global replication pair")
    decision = json.loads(INPUTS["canonical_decision"].read_text())
    if decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED":
        raise ValueError("Canonical family status changed")

    region_table = []
    for r in regions:
        region_table.append({
            "geographic_region_grch37": r["geographic_region"],
            "pair_specific_candidates": r["n_pair_specific_candidates"],
            "pairs": r["pairs"], "candidate_locus_ids": r["candidate_locus_ids"],
            "lead_variants": r["lead_variants"],
            "best_placo_p": r["best_lead_p_placo"],
            "candidate_variant_rows": r["n_candidate_variant_rows"],
            "valid_signed_ld_rows": r["n_valid_signed_ld_rows"],
            "out_of_range_raw_signed_r": r["n_out_of_range_raw_signed_r"],
            "independent_cross_pair_1000g_edges": r["n_independent_cross_pair_1000g_ld_edges"],
            "pair_level_global_replication": r["pair_level_global_replication"],
            "independent_locus_replication": "NOT_ESTABLISHED",
            "lead_nominal_eqtl_rows": r["gtex_brain_lead_eqtl_distinct_nominal_rows"],
            "candidate_molecular_cs_variants": r["gtex_brain_qtl_cs_candidate_variants"],
            "max_molecular_eqtl_pip": r["gtex_brain_eqtl_cs_max_molecular_pip"],
            "max_molecular_sqtl_pip": r["gtex_brain_sqtl_cs_max_molecular_pip"],
            "exact_lead_regulatory_features": r["ensembl_exact_lead_regulatory_feature_count"],
            "prior_trait_relevant_catalog_rows": r["prior_catalog_trait_relevant_records"],
            "gwas_fine_mapping": "NOT_ESTIMABLE_METHOD_PREREQUISITES",
            "trait_trait_coloc": "NOT_ESTIMABLE_METHOD_PREREQUISITES",
            "gwas_qtl_coloc": "NOT_ESTIMABLE_FULL_CIS_QTL_MISSING",
            "tissue_cell_type_enrichment": "NOT_ESTIMABLE_REVIEWED_BACKGROUND_MISSING",
            "pathway_enrichment": "NOT_ESTIMABLE_REVIEWED_GENE_MAPPING_MISSING",
            "lava_local_rg": "BLOCKED_FROZEN_FAMILY_QC",
            "final_shared_locus_tier": "NONE_NOT_PROMOTED",
            "analysis_label": "EXPLORATORY_DESCRIPTIVE",
        })
    stable_write("region_evidence_20.tsv", region_table)

    candidate_table = []
    replication_table = []
    qtl_table = []
    regulatory_table = []
    tissue_table = []
    method_table = []
    for locus_id in sorted(member_by_id):
        m, c, y = member_by_id[locus_id], canonical_by_id[locus_id], yale_by_id[locus_id]
        r = region_by_id[m["region_group"]]
        if m["region_group"] != c["geographic_region"] or m["region_group"] != y["region_group"]:
            raise ValueError(f"Region mismatch: {locus_id}")
        if m["pair_id"] != c["pair_id"] or m["pair_id"] != y["pair_id"]:
            raise ValueError(f"Pair mismatch: {locus_id}")
        if c["confirmatory_promotion"] != "NO" or c["confirmatory_family_status"] != "FAILED_QC_NOT_PROMOTED":
            raise ValueError(f"Unexpected promotion/family state: {locus_id}")
        qrows, nrows, frows = qtl_by_id[locus_id], nominal_by_id[locus_id], reg_by_id[locus_id]
        genes = near_by_id[locus_id]
        publications = catalog_by_id[locus_id]
        expression = [x for x in qrows if x["qtl_method"] == "gene_expression"]
        splicing = [x for x in qrows if x["qtl_method"] == "leafcutter_splicing"]
        global_row = global_by_pair.get(m["pair_id"])
        common = {"candidate_locus_id": locus_id, "geographic_region_grch37": m["region_group"],
                  "pair_id": m["pair_id"], "lead_variants": m["lead_variants"]}
        candidate_table.append({**common, "placo_lead_p": m["lead_p_placo"],
                                "canonical_lava_locus_id": c["canonical_lava_locus_id"],
                                "canonical_trait1_status": c["canonical_trait1_status"],
                                "canonical_trait2_status": c["canonical_trait2_status"],
                                "both_strict_univariate_gates_pass": c["canonical_both_strict_univariate_gates_pass"],
                                "independent_locus_replication": "NOT_ESTABLISHED",
                                "exploratory_adjacent_long_sleep_lookup": y["classification"] if m["pair_id"].startswith("longsleep__") else "NOT_APPLICABLE",
                                "molecular_qtl_candidate_variants": len({x["snp"] for x in qrows}),
                                "exact_lead_regulatory_features": len(frows),
                                "prior_catalog_rows": len(publications),
                                "lava_family_status": c["confirmatory_family_status"],
                                "final_tier": "NONE_NOT_PROMOTED", "analysis_label": "EXPLORATORY"})
        replication_table.append({**common,
                                  "global_pair_replication_source": global_row["replication_source"] if global_row else "NONE",
                                  "global_pair_replication_class": global_row["replication_class"] if global_row else "NONE",
                                  "global_result_is_locus_specific": "NO",
                                  "adjacent_long_sleep_source": c["exploratory_source"] if m["pair_id"].startswith("longsleep__") else "NOT_APPLICABLE",
                                  "adjacent_phenotype_relation": y["phenotype_relation"],
                                  "adjacent_exact_lead": y["best_exact_lead"],
                                  "adjacent_exact_p": y["best_exact_p"],
                                  "adjacent_exact_direction": y["best_exact_direction"],
                                  "adjacent_exact_classification": y["classification"],
                                  "independent_two_trait_locus_replication": "NOT_ESTABLISHED"})
        qtl_table.append({**common, "candidate_variants_in_published_molecular_sets": len({x["snp"] for x in qrows}),
                          "published_molecular_credible_sets": len({(x["dataset_id"], x["cs_id"]) for x in qrows}),
                          "expression_membership_rows": len(expression),
                          "splicing_membership_rows": len(splicing),
                          "max_published_expression_pip": max((float(x["molecular_pip"]) for x in expression), default="NA"),
                          "max_published_splicing_pip": max((float(x["molecular_pip"]) for x in splicing), default="NA"),
                          "source_sample_groups": separated({x["sample_group"] for x in qrows}),
                          "nominal_exact_lead_eqtl_rows": len(nrows),
                          "minimum_nominal_eqtl_p": min((float(x["nominal_pvalue"]) for x in nrows), default="NA"),
                          "gwas_pip": "NOT_ESTIMATED", "gwas_qtl_coloc": "NOT_ESTIMATED",
                          "interpretation": "EXPLORATORY_MOLECULAR_CONTEXT_ONLY"})
        regulatory_table.append({**common,
                                 "nearest_positional_gene_symbols": separated({x["gene_symbol"] for x in genes}),
                                 "nearest_positional_gene_ids": separated({x["gene_id"] for x in genes}),
                                 "exact_lead_regulatory_feature_ids": separated({x["feature_id"] for x in frows}),
                                 "exact_lead_regulatory_feature_types": separated({x["description"] for x in frows}),
                                 "prior_trait_relevant_catalog_studies": separated({x["study_accession"] for x in publications}),
                                 "prior_trait_relevant_catalog_rows": len(publications),
                                 "causal_gene": "NOT_ESTABLISHED", "analysis_label": "POSITIONAL_CONTEXT_ONLY"})
        tissue_table.append({**common,
                             "nominal_eqtl_source_tissues": separated({x["sample_group"] for x in nrows}),
                             "published_molecular_set_source_tissues": separated({x["sample_group"] for x in qrows}),
                             "published_molecular_set_gene_ids": separated({x["gene_id"] for x in qrows}),
                             "cell_type_source": "NONE_VERIFIED_FOR_THIS_CANDIDATE",
                             "tissue_enrichment": "NOT_ESTIMATED",
                             "cell_type_enrichment": "NOT_ESTIMATED",
                             "interpretation": "SOURCE_TISSUE_DESCRIPTIVE_ONLY"})
        method_table.append({**common,
                             "dense_gwas_available": "YES_SOURCE_HASHED_SEMANTICS_UNREVIEWED",
                             "full_regional_allele_aligned_ld": "NOT_VERIFIED",
                             "trait_specific_n_case_fraction_contract": "NOT_VERIFIED",
                             "pinned_susieR": "UNAVAILABLE",
                             "gwas_fine_mapping": "NOT_ESTIMABLE",
                             "pinned_coloc": "UNAVAILABLE",
                             "trait_trait_coloc": "NOT_ESTIMABLE",
                             "full_cis_qtl": "NOT_BOUND",
                             "eqtl_sqtl_coloc": "NOT_ESTIMABLE",
                             "reviewed_gene_universe_weights": "NOT_DEFINED",
                             "pathway_enrichment": "NOT_ESTIMABLE",
                             "lava_local_rg": "BLOCKED_FROZEN_FAMILY_QC"})
    for name, table in (
        ("candidate_evidence_25.tsv", candidate_table), ("replication_25.tsv", replication_table),
        ("qtl_context_25.tsv", qtl_table), ("regulatory_gene_25.tsv", regulatory_table),
        ("tissue_context_25.tsv", tissue_table), ("method_feasibility_25.tsv", method_table),
    ):
        stable_write(name, table)
    stable_write("prior_literature_records.tsv", [{**row, "interpretation": "PRIOR_TRAIT_CONTEXT_NOT_LOCUS_REPLICATION"}
                                                   for row in catalog])
    provenance = {
        "schema_version": 1, "analysis_label": "BOUNDED_EXPLORATORY_MANUSCRIPT_PACKAGE",
        "canonical_family_status": "FAILED_QC_NOT_PROMOTED", "candidate_count": 25, "geographic_region_count": 20,
        "source_sha256": {str(path.relative_to(ROOT)): digest(path) for path in INPUTS.values()},
        "builder_sha256": digest(Path(__file__)),
        "output_sha256": {p.name: digest(p) for p in sorted(OUT.glob("*.tsv"))},
        "method_note": "Deterministic descriptive joins only; missing method states are not negative biological results.",
    }
    provenance_path = OUT / "table_provenance.json"
    content = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    if provenance_path.exists() and provenance_path.read_text() != content:
        raise FileExistsError(f"Refusing to overwrite changed provenance: {provenance_path}")
    provenance_path.write_text(content)


if __name__ == "__main__":
    main()
