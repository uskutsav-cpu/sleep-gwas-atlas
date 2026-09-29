#!/usr/bin/env python3
"""Materialize explicit not-estimable method tables from source-bound candidate rows."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/brain6_bounded_manuscript_v1"
SOURCES = {
    "candidate": OUT / "candidate_evidence_25.tsv",
    "region": OUT / "region_evidence_20.tsv",
    "qtl": OUT / "qtl_context_25.tsv",
    "method": OUT / "method_feasibility_25.tsv",
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_new(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Empty table: {path}")
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    candidate = read(SOURCES["candidate"])
    regions = read(SOURCES["region"])
    qtl = read(SOURCES["qtl"])
    method = read(SOURCES["method"])
    if not (len(candidate) == len(qtl) == len(method) == 25 and len(regions) == 20):
        raise ValueError("Frozen candidate/region counts changed")
    by_method = {r["candidate_locus_id"]: r for r in method}
    by_qtl = {r["candidate_locus_id"]: r for r in qtl}
    if set(by_method) != set(by_qtl) or set(by_method) != {r["candidate_locus_id"] for r in candidate}:
        raise ValueError("Candidate IDs differ between method inputs")
    fine, coloc, molecular = [], [], []
    for c in candidate:
        locus = c["candidate_locus_id"]
        m, q = by_method[locus], by_qtl[locus]
        key = {k: c[k] for k in ("candidate_locus_id", "geographic_region_grch37", "pair_id", "lead_variants")}
        fine.append({**key, "dense_gwas_available": m["dense_gwas_available"],
                     "trait_specific_n_case_fraction_contract": m["trait_specific_n_case_fraction_contract"],
                     "full_regional_signed_allele_aligned_ld": m["full_regional_allele_aligned_ld"],
                     "pinned_susieR": m["pinned_susieR"], "gwas_pip": "NOT_ESTIMATED",
                     "status": "NOT_ESTIMABLE_METHOD_PREREQUISITES",
                     "interpretation": "NO_FINE_MAP_INFERENCE"})
        coloc.append({**key, "paired_dense_gwas": "FILES_PRESENT_BUT_N_AND_ALLELES_UNREVIEWED",
                      "regional_ld_and_multisignal_inputs": "NOT_VERIFIED",
                      "sample_overlap_and_prior_contract": "NOT_REVIEWED_FOR_THIS_METHOD",
                      "pinned_coloc": m["pinned_coloc"], "pp_h4": "NOT_ESTIMATED",
                      "status": "NOT_ESTIMABLE_METHOD_PREREQUISITES",
                      "interpretation": "NO_COLOCALIZATION_INFERENCE"})
        molecular.append({**key,
                          "published_eqtl_membership_rows": q["expression_membership_rows"],
                          "published_sqtl_membership_rows": q["splicing_membership_rows"],
                          "max_molecular_eqtl_pip": q["max_published_expression_pip"],
                          "max_molecular_sqtl_pip": q["max_published_splicing_pip"],
                          "nominal_exact_lead_eqtl_rows": q["nominal_exact_lead_eqtl_rows"],
                          "full_cis_qtl": m["full_cis_qtl"],
                          "gwas_eqtl_coloc_pp_h4": "NOT_ESTIMATED",
                          "gwas_sqtl_coloc_pp_h4": "NOT_ESTIMATED",
                          "interpretation": "EXPLORATORY_MOLECULAR_OVERLAP_ONLY"})
    pathway = [{"geographic_region_grch37": r["geographic_region_grch37"],
                "pairs": r["pairs"], "reviewed_causal_gene_weights": "NOT_AVAILABLE",
                "reviewed_gene_universe": "NOT_DEFINED", "pathway_collection": "NOT_FROZEN",
                "multiplicity_family": "NOT_FROZEN", "pathway_p": "NOT_ESTIMATED",
                "status": "NOT_ESTIMABLE_METHOD_PREREQUISITES",
                "interpretation": "POSITIONAL_GENES_NOT_USED_AS_ENRICHMENT_INPUT"} for r in regions]
    outputs = {
        "fine_mapping_25.tsv": fine,
        "trait_trait_coloc_25.tsv": coloc,
        "eqtl_sqtl_evidence_25.tsv": molecular,
        "pathway_20.tsv": pathway,
    }
    for name, rows in outputs.items():
        write_new(OUT / name, rows)
    provenance = {
        "schema_version": 1, "analysis_label": "BOUNDED_METHOD_STATUS_ADDENDUM",
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in SOURCES.values()},
        "builder_sha256": sha(Path(__file__)),
        "output_sha256": {name: sha(OUT / name) for name in outputs},
        "status": "NO_NEW_FINE_MAPPING_COLOCALIZATION_OR_PATHWAY_ESTIMATES",
    }
    (OUT / "method_tables_addendum.provenance.json").open("x").write(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
