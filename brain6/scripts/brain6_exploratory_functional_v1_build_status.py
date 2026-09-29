#!/usr/bin/env python3
"""Build bounded Phase 11–14 status from frozen candidate and exploratory receipts."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/brain6_exploratory_functional_v1"
SSD = Path("/Volumes/Extreme SSD/brain6-work/brain6-exploratory-functional-v1")
BASE = ROOT / "brain6/results/brain6_bounded_manuscript_v1"
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
REGULATORY = ROOT / "brain6/results/exploratory_five_track_v6/lead_regulatory_context/lead_regulatory_summary.tsv"
QTL_CONTEXT = ROOT / "brain6/results/exploratory_five_track_v5/gtex_brain_qtl_credible_sets/candidate_molecular_credible_set_memberships.tsv"
GWAS_MANIFEST = ROOT / "brain6/manifests/gwas_external_availability.tsv"
DATASETS = {"QTD000171": ("brain_cortex", "expression", "all"),
            "QTD000175": ("brain_cortex", "splicing", "cc"),
            "QTD000176": ("brain_frontal_cortex", "expression", "all"),
            "QTD000180": ("brain_frontal_cortex", "splicing", "cc")}
GATE = OUT / "coloc_gate.tsv"
COLOC = OUT / "coloc_results.tsv"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]) -> None:
    with path.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    candidates = read(MEMBERS)
    eqtl = {r["candidate_locus_id"]: r for r in read(BASE / "eqtl_sqtl_evidence_25.tsv")}
    repl = {r["candidate_locus_id"]: r for r in read(BASE / "replication_25.tsv")}
    reg = {r["candidate_locus_id"]: r for r in read(BASE / "regulatory_gene_25.tsv")}
    lead_reg = defaultdict(list)
    for r in read(REGULATORY):
        lead_reg[r["locus_id"]].append(r)
    memberships = read(QTL_CONTEXT)
    gwas_manifest = {r["trait"]: r for r in read(GWAS_MANIFEST)}
    gate = read(GATE)
    coloc = read(COLOC)
    if (len(candidates), len(eqtl), len(repl), len(reg), len({r["region_group"] for r in candidates})) != (25, 25, 25, 25, 20):
        raise ValueError("Frozen 25-candidate/20-region identity changed")
    if len(coloc) != 4 or len({r["candidate_locus_id"] for r in coloc}) != 1:
        raise ValueError("Unexpected numerical coloc result family")
    by_candidate_coloc = defaultdict(list)
    for r in coloc:
        by_candidate_coloc[r["candidate_locus_id"]].append(r)
    candidate_rows = []
    for c in candidates:
        cid = c["candidate_locus_id"]
        pair = c["pair_id"]
        cs = by_candidate_coloc[cid]
        gstatus = ("ADHD_COMPONENT_EXPLORATORY_ABF_ONLY" if cs else
                   "NOT_ESTIMATED_NO_ELIGIBLE_MOLECULAR_CONTEXT" if pair == "insomnia__adhd" else
                   "NOT_ESTIMATED_SOURCE_N_MODEL_HOLD")
        regcount = sum(int(x["n_regulatory_features_overlapping_exact_lead"]) for x in lead_reg[cid])
        candidate_rows.append({
            "candidate_locus_id": cid, "pair_id": pair, "region_grch37": c["region_group"],
            "lead_variants": c["lead_variants"], "canonical_lava": "FAILED_QC_NOT_PROMOTED",
            "published_eqtl_cs_memberships": eqtl[cid]["published_eqtl_membership_rows"],
            "published_sqtl_cs_memberships": eqtl[cid]["published_sqtl_membership_rows"],
            "nominal_exact_lead_eqtl_rows": eqtl[cid]["nominal_exact_lead_eqtl_rows"],
            "exact_lead_regulatory_overlaps": regcount,
            "gwas_qtl_coloc_status": gstatus, "gwas_qtl_coloc_component": "adhd" if cs else "NA",
            "n_eligible_coloc_tests": len(cs), "max_exploratory_pp_h4": max((float(x["pp_h4"]) for x in cs), default="NA"),
            "sleep_component_coloc": "NOT_ESTIMATED", "independent_two_trait_locus_replication": "NOT_ESTABLISHED",
            "global_pair_replication_class": repl[cid]["global_pair_replication_class"],
            "regulatory_inference": "DESCRIPTIVE_COORDINATE_OVERLAP_ONLY",
            "tissue_cell_type_enrichment": "NOT_ESTIMATED_MISSING_WEIGHTED_SET_AND_BACKGROUND",
            "pathway_enrichment": "NOT_ESTIMATED_MISSING_CAUSAL_GENE_WEIGHTS_UNIVERSE",
            "final_shared_locus_tier": "BLOCKED_LAVA"})
    candidate_rows.sort(key=lambda x: (x["pair_id"], x["region_grch37"], x["candidate_locus_id"]))
    candidate_columns = tuple(candidate_rows[0])
    write(OUT / "candidate_functional_status_25.tsv", candidate_rows, candidate_columns)
    region_map = defaultdict(list)
    for c in candidate_rows:
        region_map[c["region_grch37"]].append(c)
    region_rows = []
    for region, rows in sorted(region_map.items()):
        region_rows.append({"region_grch37": region, "n_pair_candidates": len(rows),
                            "pairs": ";".join(sorted({r["pair_id"] for r in rows})),
                            "candidate_locus_ids": ";".join(sorted(r["candidate_locus_id"] for r in rows)),
                            "published_eqtl_cs_membership_rows_pair_counted": sum(int(r["published_eqtl_cs_memberships"]) for r in rows),
                            "published_sqtl_cs_membership_rows_pair_counted": sum(int(r["published_sqtl_cs_memberships"]) for r in rows),
                            "exact_lead_regulatory_feature_rows_pair_counted": sum(int(r["exact_lead_regulatory_overlaps"]) for r in rows),
                            "n_component_coloc_tests": sum(int(r["n_eligible_coloc_tests"]) for r in rows),
                            "independent_two_trait_locus_replication": "NOT_ESTABLISHED",
                            "tissue_cell_type_pathway_inference": "NOT_ESTIMATED",
                            "final_shared_locus_tier": "BLOCKED_LAVA"})
    write(OUT / "geographic_region_status_20.tsv", region_rows, tuple(region_rows[0]))
    source_rows = []
    for dataset, (sample, method, suffix) in DATASETS.items():
        receipts = sorted((SSD / "qtl_regional").glob(dataset + "_*.receipt.json"))
        if len(receipts) != 6:
            raise ValueError(f"Incomplete QTL regional source: {dataset}")
        rec = [json.loads(p.read_text()) for p in receipts]
        source_rows.append({"source_id": dataset, "source_kind": "GTEX_QTL_REGIONAL",
                            "phenotype_or_tissue": sample + "_" + method,
                            "url": f"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{dataset}/{dataset}.{suffix}.tsv.gz",
                            "original_file_hash": "NOT_DOWNLOADED_REMOTE_TABIX",
                            "local_receipt_count": len(rec), "local_rows": sum(x["n_rows"] for x in rec),
                            "n_semantics": "AN_DIVIDED_BY_2_AUTOSOMAL",
                            "model": "QTL_CATALOGUE_LINEAR_REGRESSION", "decision": "EXPLORATORY_REGIONAL_INPUT"})
    for trait, nmodel in (("insomnia", "LITERAL_PER_SNP_N"),
                          ("adhd", "NCA_PLUS_NCO_PER_SNP")):
        receipts = sorted((SSD / "gwas_regional").glob(trait + "_*.receipt.json"))
        if len(receipts) != 6:
            raise ValueError(f"Incomplete GWAS regional source: {trait}")
        rec = [json.loads(p.read_text()) for p in receipts]
        source_rows.append({"source_id": trait, "source_kind": "RAW_GWAS_REGIONAL",
                            "phenotype_or_tissue": trait, "url": gwas_manifest[trait]["source_URL"],
                            "original_file_hash": rec[0]["source_sha256"],
                            "local_receipt_count": len(rec), "local_rows": sum(x["n_rows"] for x in rec),
                            "n_semantics": nmodel, "model": "SOURCE_LOG_OR_AND_SE",
                            "decision": "EXPLORATORY_COLOC_GATE_ONLY"})
    for trait in ("longsleep", "mdd", "scz", "bipolar", "parkinson"):
        source_rows.append({"source_id": trait, "source_kind": "BRAIN6_GWAS", "phenotype_or_tissue": trait,
                            "url": "SEE_BRAIN6_SOURCE_MANIFEST", "original_file_hash": "SEE_BRAIN6_SOURCE_MANIFEST",
                            "local_receipt_count": 0, "local_rows": 0, "n_semantics": "NOT_REVIEWED_FOR_COLOC",
                            "model": "SOURCE_SPECIFIC", "decision": "NO_NEW_COLOC_N_MODEL_HOLD"})
    source_rows.extend([
        {"source_id": "FinnGen_R13_F5_ADHD", "source_kind": "REPLICATION", "phenotype_or_tissue": "adhd",
         "url": "https://www.finngen.fi/en/access_results", "original_file_hash": "PAIR_LEVEL_RECEIPT_ONLY",
         "local_receipt_count": 0, "local_rows": 0, "n_semantics": "NOT_APPLICABLE",
         "model": "GLOBAL_GENETIC_CORRELATION", "decision": "NOT_LOCUS_SPECIFIC_TWO_TRAIT_REPLICATION"},
        {"source_id": "Ensembl_GRCh37_exact_lead_regulatory", "source_kind": "REGULATORY",
         "phenotype_or_tissue": "unspecified_regulatory_features", "url": "https://grch37.rest.ensembl.org/overlap/region/human/",
         "original_file_hash": sha(REGULATORY), "local_receipt_count": 27, "local_rows": sum(int(x["exact_lead_regulatory_overlaps"]) for x in candidate_rows),
         "n_semantics": "NOT_APPLICABLE", "model": "COORDINATE_OVERLAP",
         "decision": "DESCRIPTIVE_ONLY"},
        {"source_id": "tissue_cell_type_pathway", "source_kind": "ENRICHMENT",
         "phenotype_or_tissue": "not_defined", "url": "NA", "original_file_hash": "NA",
         "local_receipt_count": 0, "local_rows": 0, "n_semantics": "NOT_APPLICABLE",
         "model": "NO_WEIGHTED_GENE_SET_OR_BACKGROUND", "decision": "NOT_ESTIMATED"},
    ])
    write(OUT / "source_gate.tsv", source_rows, tuple(source_rows[0]))
    tracked = [MEMBERS, REGULATORY, QTL_CONTEXT, GWAS_MANIFEST, GATE, COLOC,
               BASE / "eqtl_sqtl_evidence_25.tsv", BASE / "replication_25.tsv", BASE / "regulatory_gene_25.tsv"]
    outputs = [OUT / x for x in ("candidate_functional_status_25.tsv", "geographic_region_status_20.tsv", "source_gate.tsv")]
    receipt = {"analysis_label": "EXPLORATORY_NOT_PROMOTED", "candidate_count": 25, "geographic_region_count": 20,
               "qtl_regional_slice_count": 24, "gwas_regional_slice_count": 12,
               "coloc_gate_rows": len(gate), "coloc_results": len(coloc),
               "candidate_coloc_status_counts": dict(Counter(x["gwas_qtl_coloc_status"] for x in candidate_rows)),
               "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in tracked},
               "output_sha256": {p.name: sha(p) for p in outputs},
               "ssd_qtl_receipt_sha256": {p.name: sha(p) for p in sorted((SSD / "qtl_regional").glob("*.receipt.json"))},
               "ssd_gwas_receipt_sha256": {p.name: sha(p) for p in sorted((SSD / "gwas_regional").glob("*.receipt.json"))},
               "coloc_package": "coloc_5.2.3", "coloc_priors": {"p1": 1e-4, "p2": 1e-4, "p12": 1e-5}}
    with (OUT / "provenance.json").open("x") as f:
        json.dump(receipt, f, indent=2, sort_keys=True)
        f.write("\n")
    print(json.dumps({"candidates": 25, "regions": 20, "coloc_tests": len(coloc),
                      "status": receipt["candidate_coloc_status_counts"]}, indent=2))


if __name__ == "__main__":
    main()
