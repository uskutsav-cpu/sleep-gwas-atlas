#!/usr/bin/env python3
"""Descriptive lead-Z versus local-covariance sign audit; no status changes."""

import csv
import gzip
import hashlib
import json
import pathlib


REPO = pathlib.Path(__file__).resolve().parents[2]
RESULTS = REPO / "brain6/results/brain6_alternative_local_validation_v1"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sign(value):
    return "POSITIVE" if value > 0 else "NEGATIVE" if value < 0 else "ZERO"


def main():
    protocol = json.loads((RESULTS / "protocol_freeze.json").read_text())
    lock = json.loads((RESULTS / "source_eligibility_lock_v1.json").read_text())
    if sha256(RESULTS / "protocol_freeze.json") != lock["parent_protocol_sha256"]:
        raise RuntimeError("Eligibility lock/protocol mismatch")
    candidate_path = RESULTS / "candidate_crosswalk.tsv"
    with candidate_path.open() as stream:
        candidates = list(csv.DictReader(stream, delimiter="\t"))
    if len(candidates) != 25:
        raise RuntimeError("Candidate cardinality mismatch")
    leads = {rs for row in candidates for rs in row["lead_variants"].split(";")}
    reference = {}
    with pathlib.Path(protocol["adapter_reference_bim"]).open() as stream:
        for line in stream:
            chrom, rs, cm, bp, a1, a2 = line.split()
            if rs in leads:
                reference[rs] = (a1, a2)
    traits = {trait for pair, names in protocol["pairs"].items()
              if lock["pair_status"][pair] == "ELIGIBLE_SECONDARY" for trait in names}
    z_by_trait = {}
    adapter_receipts = {}
    root = pathlib.Path(protocol["ssd_root"]) / "sumstats_v1"
    for trait in traits:
        source = root / f"{trait}.sumstats.gz"
        receipt = root / f"{trait}.adapter_receipt.json"
        card = json.loads(receipt.read_text())
        if sha256(source) != card["output_sha256"]:
            raise RuntimeError(f"Adapter receipt mismatch: {trait}")
        adapter_receipts[trait] = sha256(receipt)
        z_by_trait[trait] = {}
        with gzip.open(source, "rt") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                rs = row["SNP"]
                if rs not in leads or rs not in reference:
                    continue
                a1, a2 = reference[rs]
                if (row["A1"], row["A2"]) == (a1, a2):
                    factor = 1
                elif (row["A1"], row["A2"]) == (a2, a1):
                    factor = -1
                else:
                    continue
                z_by_trait[trait][rs] = factor * float(row["Z"])
    rows = []
    for candidate in candidates:
        pair = candidate["pair_id"]
        eligibility = lock["pair_status"][pair]
        a, b = protocol["pairs"][pair]
        lead_signs = []
        if eligibility == "ELIGIBLE_SECONDARY":
            for rs in candidate["lead_variants"].split(";"):
                if rs in z_by_trait[a] and rs in z_by_trait[b]:
                    lead_signs.append(f"{rs}:{sign(z_by_trait[a][rs] * z_by_trait[b][rs])}")
                else:
                    lead_signs.append(f"{rs}:UNAVAILABLE")
        direction_set = {item.split(":", 1)[1] for item in lead_signs}
        lead_direction = (next(iter(direction_set)) if len(direction_set) == 1 else
                          "MIXED" if direction_set else "UNAVAILABLE")
        covariance = candidate["covariance_sign"] or "UNAVAILABLE"
        comparison = ("UNAVAILABLE" if eligibility == "METHOD_INAPPLICABLE" or
                      lead_direction == "UNAVAILABLE" or covariance == "UNAVAILABLE" else
                      "MIXED_LEADS" if lead_direction == "MIXED" else
                      "CONCORDANT" if lead_direction == covariance else "DISCORDANT")
        rows.append({"candidate_locus_id": candidate["candidate_locus_id"],
                     "pair_id": pair, "source_eligibility": eligibility,
                     "lead_variant_reference_aligned_z_product_signs": ";".join(lead_signs),
                     "lead_cross_trait_sign": lead_direction,
                     "best_local_covariance_sign": covariance,
                     "descriptive_sign_comparison": comparison,
                     "candidate_status": candidate["status"]})
    output = RESULTS / "candidate_direction_audit.tsv"
    receipt = RESULTS / "candidate_direction_audit_receipt.json"
    if output.exists() or receipt.exists():
        raise FileExistsError("Direction audit already exists")
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt.write_text(json.dumps({
        "protocol_sha256": sha256(RESULTS / "protocol_freeze.json"),
        "source_eligibility_lock_sha256": sha256(RESULTS / "source_eligibility_lock_v1.json"),
        "candidate_crosswalk_sha256": sha256(candidate_path),
        "reference_bim_sha256": sha256(pathlib.Path(protocol["adapter_reference_bim"])),
        "adapter_receipt_sha256": adapter_receipts,
        "output_sha256": sha256(output),
        "n_candidates": len(rows),
        "purpose": "Descriptive sign concordance only; no candidate status changes"
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
