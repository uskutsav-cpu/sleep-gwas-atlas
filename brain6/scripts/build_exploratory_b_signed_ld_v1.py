#!/usr/bin/env python3
"""Exploratory signed UKB LD and effect directions for admitted Track B only."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from build_partial_placo_signed_effect_ld import (  # noqa: E402
    allele_sign, ld_direction_consistency, pair_effects, read_tsv,
)

PAIR = "insomnia__adhd"
SOURCE = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_variants.tsv"
B_LEDGER = ROOT / "brain6/results/placo/insomnia__adhd/B.full.tsv.gz"
B_RECEIPT = ROOT / "brain6/results/placo/insomnia__adhd/admission_receipt.json"
REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
REFERENCE_RECEIPT = REFERENCE.parent / "reference.provenance.json"
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
EXTRACTOR = ROOT / "brain6/scripts/extract_placo_candidate_lead_ld.R"
OUT = ROOT / "brain6/results/exploratory_five_track_v1/b_signed_ld"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(path: Path, data: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def conservative_sign(ref_a1: str, ref_a2: str, a1: str, a2: str) -> int | None:
    if {ref_a1.upper(), ref_a2.upper()} in ({"A", "T"}, {"C", "G"}):
        return None
    return allele_sign(ref_a1, ref_a2, a1, a2)


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned B signed-LD output already exists")
    admission = json.loads(B_RECEIPT.read_text())
    if admission["status"] != "PROMOTED_AFTER_ALL_CHECKS_PASS" or sha(B_LEDGER) != admission["files"]["B.full.tsv.gz"]["sha256"]:
        raise ValueError("Protected B ledger/receipt mismatch")
    candidate_receipt = json.loads((SOURCE.parent / "provenance.json").read_text())
    if sha(SOURCE) != candidate_receipt["candidate_variants_sha256"]:
        raise ValueError("Five-track candidate table hash changed")
    candidates = [r for r in read_tsv(SOURCE) if r["pair_id"] == PAIR]
    if (len(candidates) != 389 or len({r["SNP"] for r in candidates}) != 389 or
            sum(r["candidate_status"] == "LEAD" for r in candidates) != 6 or
            any(r["reference_status"] != "EXACT_MATCH" for r in candidates)):
        raise ValueError("Admitted B candidate invariant failed")
    wanted = {r["SNP"].lower() for r in candidates}
    wanted.update(r["lead_SNP"].lower() for r in candidates)
    effects, pair_source = pair_effects(PAIR, wanted)
    if set(effects) != wanted:
        raise ValueError("B candidate missing from receipt-verified joined pair")
    archived = {}
    with gzip.open(B_LEDGER, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row["SNP"].lower()
            if key in wanted:
                if key in archived:
                    raise ValueError(f"Duplicate protected B SNP: {key}")
                archived[key] = row
    if set(archived) != wanted:
        raise ValueError("Protected B ledger does not contain every candidate and lead")
    for snp in wanted:
        a, b = effects[snp], archived[snp]
        if ((a["CHR"], a["BP"], a["A1"].upper(), a["A2"].upper()) !=
                (b["CHR"], b["BP"], b["A1"].upper(), b["A2"].upper())):
            raise ValueError(f"Protected B allele/coordinate disagreement: {snp}")
        for key in ("Z1", "Z2"):
            if not math.isclose(float(a[key]), float(b[key]), rel_tol=1e-11, abs_tol=1e-11):
                raise ValueError(f"Protected B Z disagreement: {snp}/{key}")

    with tempfile.TemporaryDirectory(prefix="brain6-b-signed-ld-") as temporary:
        temp = Path(temporary)
        request = temp / "candidates.tsv"
        write(request, [{k: row[k] for k in ("pair_id", "SNP", "CHR", "BP", "lead_SNP")}
                        for row in candidates], ["pair_id", "SNP", "CHR", "BP", "lead_SNP"])
        result = subprocess.run([str(RSCRIPT), str(EXTRACTOR), str(REFERENCE), str(request),
                                 str(temp / "ld")], capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError("Pinned LAVA signed-LD extraction failed: " +
                               (result.stdout + result.stderr)[-3000:])
        info = {(int(r["CHR"]), r["SNP"].lower()): r for r in
                read_tsv(temp / "ld/reference_alleles.tsv")}
        edges = {(r["pair_id"], r["SNP"].lower(), r["lead_SNP"].lower()): float(r["R"])
                 for r in read_tsv(temp / "ld/candidate_lead_ld.tsv")}
        sample_size = (temp / "ld/reference_sample_size.txt").read_text().strip()
    if sample_size != "100000":
        raise ValueError("Unexpected official UKB LAVA reference N")

    output = []
    for row in candidates:
        snp, lead, chrom, bp = row["SNP"].lower(), row["lead_SNP"].lower(), int(row["CHR"]), int(row["BP"])
        candidate_ref, lead_ref = info.get((chrom, snp)), info.get((chrom, lead))
        candidate_effect, lead_effect = effects[snp], effects[lead]
        status = "NO_EXACT_REFERENCE_OR_LEAD"
        r = edges.get((PAIR, snp, lead))
        beta1 = beta2 = beta1_lead = beta2_lead = None
        if candidate_ref and lead_ref:
            if int(candidate_ref["POS"]) != bp or int(candidate_effect["BP"]) != bp:
                raise ValueError(f"Reference/effect coordinate mismatch: {snp}")
            sign = conservative_sign(candidate_ref["A1"], candidate_ref["A2"],
                                     candidate_effect["A1"], candidate_effect["A2"])
            lead_sign = conservative_sign(lead_ref["A1"], lead_ref["A2"],
                                          lead_effect["A1"], lead_effect["A2"])
            if sign is None or lead_sign is None:
                status = "AMBIGUOUS_OR_UNALIGNABLE_ALLELES"
            elif r is None:
                status = "MISSING_SIGNED_LD_EDGE"
            else:
                beta1, beta2 = float(candidate_effect["BETA1"]) * sign, float(candidate_effect["BETA2"]) * sign
                beta1_lead, beta2_lead = float(lead_effect["BETA1"]) * lead_sign, float(lead_effect["BETA2"]) * lead_sign
                status = "OUT_OF_RANGE_RAW_R" if not math.isfinite(r) or abs(r) > 1 else "PASS_SIGNED_LD"
        valid = status == "PASS_SIGNED_LD"
        output.append({
            "pair_id": PAIR, "locus_id": row["locus_id"], "SNP": row["SNP"], "lead_SNP": row["lead_SNP"],
            "CHR": chrom, "BP": bp, "candidate_status": row["candidate_status"],
            "reference_A1": candidate_ref["A1"] if candidate_ref else "NA",
            "reference_A2": candidate_ref["A2"] if candidate_ref else "NA",
            "beta_sleep_reference_A1": format(beta1, ".12g") if beta1 is not None else "NA",
            "beta_adhd_reference_A1": format(beta2, ".12g") if beta2 is not None else "NA",
            "raw_r_reference_A1": format(r, ".12g") if r is not None else "NA",
            "signed_r_reference_A1": format(r, ".12g") if valid else "NA",
            "r2_reference": format(r * r, ".12g") if valid else "NA",
            "sleep_effect_vs_lead_ld": ld_direction_consistency(beta1, beta1_lead, r) if valid else "NOT_ASSESSED",
            "adhd_effect_vs_lead_ld": ld_direction_consistency(beta2, beta2_lead, r) if valid else "NOT_ASSESSED",
            "ld_qc_status": status,
            "analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA",
        })
    OUT.mkdir(parents=True, exist_ok=False)
    target = OUT / "candidate_signed_ld.tsv"
    write(target, output, list(output[0]))
    statuses = Counter(r["ld_qc_status"] for r in output)
    provenance = {
        "analysis_id": "brain6_admitted_b_exploratory_signed_ld_v1",
        "analysis_label": "EXPLORATORY", "candidate_rows": len(output),
        "lead_rows": sum(r["candidate_status"] == "LEAD" for r in output),
        "status_counts": dict(statuses),
        "source_sha256": {
            str(SOURCE.relative_to(ROOT)): sha(SOURCE),
            str(B_LEDGER.relative_to(ROOT)): sha(B_LEDGER),
            str(B_RECEIPT.relative_to(ROOT)): sha(B_RECEIPT),
            str(EXTRACTOR.relative_to(ROOT)): sha(EXTRACTOR),
            "joined_pair_receipt": pair_source["pair_receipt_sha256"],
            "joined_pair_tsv": pair_source["pair_tsv_sha256"],
            "reference_provenance": sha(REFERENCE_RECEIPT),
            "pinned_rscript": sha(RSCRIPT),
        },
        "reference_sample_size": 100000,
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {target.name: sha(target)},
        "limitations": [
            "Post-selection descriptive signed-LD/effect-direction comparison, not independent replication or fine-mapping.",
            "Palindromic SNPs are omitted from orientation-based comparison.",
            "Raw signed R outside [-1,1] is retained only as a QC failure and never clipped.",
            "No LAVA local-rg or final region tier is promoted.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"rows": len(output), "statuses": dict(statuses)}, sort_keys=True))


if __name__ == "__main__":
    main()
