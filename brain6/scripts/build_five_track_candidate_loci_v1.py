#!/usr/bin/env python3
"""Build five-track PLACO candidate loci under the frozen UKB/LAVA LD rule."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

from clump_placo_candidates import (
    LOCUS_COLUMNS, VARIANT_COLUMNS, build_outputs, extract_reference_ld, read_candidates,
)

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "brain6/results/track_b_admission_v1/five_track_family_audit.json"
MASTER = ROOT / "brain6/results/track_b_admission_v1/placo_master_five_track_v1.tsv"
FAMILY = ROOT / "brain6/config/placo_family_v3/family_lock.json"
RULE = ROOT / "brain6/config/shared_locus_rule_v1.json"
RECEIPT = ROOT / "brain6/results/placo/insomnia__adhd/admission_receipt.json"
B_FULL = ROOT / "brain6/results/placo/insomnia__adhd/B.full.tsv.gz"
REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
OUT = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def compatible_alleles(a1: str, a2: str, b1: str, b2: str) -> bool:
    bases = {"A", "C", "G", "T"}
    if any(x not in bases for x in (a1,a2,b1,b2)) or a1 == a2 or b1 == b2:
        return False
    if {a1,a2} in ({"A","T"},{"C","G"}):
        return False
    complement = str.maketrans("ACGT","TGCA")
    return (a1,a2) in ((b1,b2),(b2,b1),(b1.translate(complement),b2.translate(complement)),
                       (b2.translate(complement),b1.translate(complement)))


def main() -> None:
    audit, family, rule, receipt = (json.loads(p.read_text()) for p in (AUDIT,FAMILY,RULE,RECEIPT))
    if audit["status"] != "PASS_FIVE_PAIR_QC_LOCUS_DEDUP_PENDING" or receipt["status"] != "PROMOTED_AFTER_ALL_CHECKS_PASS":
        raise ValueError("full five-track PLACO pair family is not admitted")
    if audit["five_pair_master_sha256"] != sha(MASTER) or audit["family_lock_sha256"] != sha(FAMILY):
        raise ValueError("admitted family identity drifted")
    if sha(B_FULL) != receipt["files"]["B.full.tsv.gz"]["sha256"]:
        raise ValueError("archived B ledger in protected slot changed")
    if family["headline_threshold"] != 1e-8 or family["n_selected_tracks"] != 5:
        raise ValueError("frozen five-track headline rule changed")
    if rule["clumping"]["r2_threshold"] != 0.1 or rule["clumping"]["window_kb"] != 500 or rule["region_definition"]["lead_flank_kb"] != 500:
        raise ValueError("frozen locus rule changed")
    if OUT.exists():
        raise FileExistsError("five-track candidate locus directory already exists")
    candidates, output_hashes = read_candidates(MASTER, float(family["headline_threshold"]))
    if set(candidates) != set(audit["pairs"]) or output_hashes != audit["output_sha256"]:
        raise ValueError("five-track candidate source set differs from family audit")
    if {pair:len(rows) for pair,rows in candidates.items()} != audit["headline_variants_by_pair"]:
        raise ValueError("candidate count differs from frozen headline audit")
    b_hits = {row["SNP"]:row for row in candidates["insomnia__adhd"]}
    b_alleles = {}
    with gzip.open(B_FULL, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            if row["SNP"] in b_hits:
                if row["SNP"] in b_alleles or int(row["CHR"]) != b_hits[row["SNP"]]["CHR"] or int(row["BP"]) != b_hits[row["SNP"]]["BP"]:
                    raise ValueError(f"B hit identity differs: {row['SNP']}")
                b_alleles[row["SNP"]] = (row["A1"],row["A2"])
    if set(b_alleles) != set(b_hits):
        raise ValueError("B source allele coverage is incomplete")
    all_candidates = [r for pair in sorted(candidates) for r in candidates[pair]]
    ref_info, ld, ref_n = extract_reference_ld(all_candidates, REFERENCE, RSCRIPT)
    if ref_n != rule["clumping"]["reference_sample_size"]:
        raise ValueError("pinned LAVA reference sample size changed")
    b_allele_match = b_allele_mismatch = b_missing = 0
    for row in candidates["insomnia__adhd"]:
        key = (int(row["CHR"]),row["SNP"])
        info = ref_info.get(key)
        if info is None or int(info["BP"]) != int(row["BP"]):
            b_missing += 1
        elif compatible_alleles(*b_alleles[row["SNP"]], info["A1"], info["A2"]):
            b_allele_match += 1
        else:
            b_allele_mismatch += 1
            # The frozen rule forbids non-exact reference matching.
            ref_info.pop(key)
    variants, loci = build_outputs(candidates, ref_info, ld, rule)
    if len(variants) != len(all_candidates):
        raise ValueError("candidate variant accounting incomplete")
    OUT.mkdir(parents=True,exist_ok=False)
    variant_path, locus_path = OUT / "candidate_variants.tsv", OUT / "candidate_loci.tsv"
    write_tsv(variant_path,VARIANT_COLUMNS,variants)
    write_tsv(locus_path,LOCUS_COLUMNS,loci)
    report = {
        "analysis_id":"brain6-five-track-placo-candidate-loci-v1",
        "status":"PASS_COMPLETE_PAIR_FAMILY_CANDIDATE_LOCI_ONLY",
        "five_track_audit_sha256":sha(AUDIT),"five_track_master_sha256":sha(MASTER),
        "family_lock_sha256":sha(FAMILY),"frozen_locus_rule_sha256":sha(RULE),
        "B_admission_receipt_sha256":sha(RECEIPT),"B_full_ledger_sha256":sha(B_FULL),
        "reference_prefix":str(REFERENCE),"reference_n":ref_n,
        "rscript_sha256":sha(RSCRIPT),
        "reference_provenance_sha256":sha(REFERENCE.parent / "reference.provenance.json"),
        "extractor_sha256":sha(ROOT / "brain6/scripts/extract_placo_ld.R"),
        "clumper_library_sha256":sha(ROOT / "brain6/scripts/clump_placo_candidates.py"),
        "builder_sha256":sha(Path(__file__)),
        "pair_output_sha256":output_hashes,
        "n_candidate_variants":len(variants),"n_candidate_loci":len(loci),
        "n_exact_reference_matched_variants":sum(r["reference_status"] == "EXACT_MATCH" for r in variants),
        "B_exact_reference_allele_matches":b_allele_match,
        "B_reference_allele_mismatches_excluded":b_allele_mismatch,
        "B_missing_or_position_mismatch":b_missing,
        "candidate_variants_sha256":sha(variant_path),"candidate_loci_sha256":sha(locus_path),
        "interpretation":"Complete five-track PLACO candidate loci under the frozen UKB LD rule. LAVA local-rg QC is a separate unmet prerequisite for final evidence tiers or downstream causal claims."
    }
    (OUT / "provenance.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:report[k] for k in ("status","n_candidate_variants","n_candidate_loci","B_exact_reference_allele_matches","B_reference_allele_mismatches_excluded")}))


if __name__ == "__main__":
    main()
