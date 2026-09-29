#!/usr/bin/env python3
"""Record all predeclared v1 SuSiE units as explicit input/LD-QC NOT_RUN."""

import csv
import json
from pathlib import Path

from audit_brain6_exploratory_finemap_coloc_v1 import sha256, write_tsv


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/brain6_exploratory_finemap_coloc_v1"


def read(name):
    return list(csv.DictReader((OUT / name).open(), delimiter="\t"))


def main():
    trait = read("fine_mapping_input_manifest.tsv")
    pair = read("trait_coloc_input_manifest.tsv")
    ld = {r["analysis_id"]: r for r in read("ld_numeric_qc.tsv")}
    if len(trait) != 50 or len(pair) != 25 or len(ld) != 32:
        raise RuntimeError("Incomplete v1 analysis family")
    if any(r["valid_for_v1_susie"] != "NO" for r in ld.values()):
        raise RuntimeError("This summary applies only to the observed 0/32 LD-QC failure")
    fine_rows = []
    for row in trait:
        status, reason = row["status"], row["reason"]
        if status == "READY":
            if row["analysis_id"] not in ld:
                raise RuntimeError("Missing LD audit for admitted trait")
            status, reason = "NOT_RUN_LD_NUMERIC_QC", "0.005_REGULARIZED_LD_FAILED_CORRELATION_OR_CHOLESKY"
        fine_rows.append({
            "candidate_locus_id": row["candidate_locus_id"], "pair_id": row["pair_id"],
            "trait_id": row["trait_id"], "status": status, "reason": reason,
            "n_input_snps": row["n_after_n_gate"], "credible_set_count": "",
            "maximum_pip": "", "analysis_label": "EXPLORATORY_SUSIE_NOT_ESTIMATED",
        })
    pair_rows = []
    for row in pair:
        status, reason = row["status"], row["reason"]
        if status == "READY":
            if row["candidate_locus_id"] + "__shared" not in ld:
                raise RuntimeError("Missing LD audit for admitted pair")
            status, reason = "NOT_RUN_LD_NUMERIC_QC", "0.005_REGULARIZED_LD_FAILED_CORRELATION_OR_CHOLESKY"
        pair_rows.append({
            "candidate_locus_id": row["candidate_locus_id"], "pair_id": row["pair_id"],
            "status": status, "reason": reason, "n_shared": row["n_shared"],
            "PP.H0": "", "PP.H1": "", "PP.H2": "", "PP.H3": "", "PP.H4": "",
            "analysis_label": "EXPLORATORY_SUSIE_COLOC_NOT_ESTIMATED",
        })
    write_tsv(OUT / "fine_mapping_50.tsv", fine_rows, list(fine_rows[0]))
    write_tsv(OUT / "trait_coloc_25.tsv", pair_rows, list(pair_rows[0]))
    provenance = {
        "status": "NOT_RUN_NUMERIC_LD_QC",
        "planned_trait_units": 50,
        "planned_pair_units": 25,
        "admitted_ld_matrices": 32,
        "ld_matrices_passing_frozen_regularization": 0,
        "method_lock_sha256": sha256(OUT / "PRE_OUTCOME_METHOD_LOCK.md"),
        "ld_numeric_qc_sha256": sha256(OUT / "ld_numeric_qc.tsv"),
        "script_sha256": sha256(Path(__file__)),
        "fine_mapping_sha256": sha256(OUT / "fine_mapping_50.tsv"),
        "trait_coloc_sha256": sha256(OUT / "trait_coloc_25.tsv"),
        "interpretation": "No SuSiE PIP/credible set or SuSiE-coloc H0-H4 generated.",
    }
    (OUT / "v1_failure_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print("v1: 50 fine-map and 25 pair rows; zero valid SuSiE outcomes")


if __name__ == "__main__":
    main()
