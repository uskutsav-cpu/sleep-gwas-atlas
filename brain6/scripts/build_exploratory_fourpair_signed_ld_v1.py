#!/usr/bin/env python3
"""Signed effect/LD diagnostic for current four-pair five-track candidates."""

from __future__ import annotations

import csv
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
    ld_direction_consistency, pair_effects, read_tsv,
)
from build_exploratory_b_signed_ld_v1 import conservative_sign, sha, write  # noqa: E402

PAIRS = ("insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz")
SOURCE = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_variants.tsv"
REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
REFERENCE_RECEIPT = REFERENCE.parent / "reference.provenance.json"
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
EXTRACTOR = ROOT / "brain6/scripts/extract_placo_candidate_lead_ld.R"
OUT = ROOT / "brain6/results/exploratory_five_track_v1/fourpair_signed_ld"


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned four-pair signed-LD output already exists")
    source_receipt = json.loads((SOURCE.parent / "provenance.json").read_text())
    if sha(SOURCE) != source_receipt["candidate_variants_sha256"]:
        raise ValueError("Five-track candidate source changed")
    candidates = [r for r in read_tsv(SOURCE) if r["pair_id"] in PAIRS]
    if len(candidates) != 2297 or sum(r["reference_status"] != "EXACT_MATCH" for r in candidates) != 24:
        raise ValueError("Four-pair 2,297/24 candidate invariant failed")
    exact = [r for r in candidates if r["reference_status"] == "EXACT_MATCH"]
    wanted_by_pair = {pair: set() for pair in PAIRS}
    for row in exact:
        wanted_by_pair[row["pair_id"]].update((row["SNP"].lower(), row["lead_SNP"].lower()))
    effects = {}
    pair_sources = {}
    for pair in PAIRS:
        effects[pair], pair_sources[pair] = pair_effects(pair, wanted_by_pair[pair])
        if set(effects[pair]) != wanted_by_pair[pair]:
            raise ValueError(f"Missing receipt-verified pair effects for {pair}")

    with tempfile.TemporaryDirectory(prefix="brain6-fourpair-signed-ld-") as temporary:
        temp = Path(temporary)
        request = temp / "candidates.tsv"
        write(request, [{k: row[k] for k in ("pair_id", "SNP", "CHR", "BP", "lead_SNP")}
                        for row in exact], ["pair_id", "SNP", "CHR", "BP", "lead_SNP"])
        proc = subprocess.run([str(RSCRIPT), str(EXTRACTOR), str(REFERENCE), str(request),
                               str(temp / "ld")], capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError("Pinned LAVA signed-LD extraction failed: " +
                               (proc.stdout + proc.stderr)[-3000:])
        info = {(int(r["CHR"]), r["SNP"].lower()): r for r in
                read_tsv(temp / "ld/reference_alleles.tsv")}
        edges = {(r["pair_id"], r["SNP"].lower(), r["lead_SNP"].lower()): float(r["R"])
                 for r in read_tsv(temp / "ld/candidate_lead_ld.tsv")}
        sample_size = (temp / "ld/reference_sample_size.txt").read_text().strip()
    if sample_size != "100000":
        raise ValueError("Unexpected UKB LAVA reference sample size")

    output = []
    for row in candidates:
        pair = row["pair_id"]
        snp, lead, chrom, bp = row["SNP"].lower(), row["lead_SNP"].lower(), int(row["CHR"]), int(row["BP"])
        ref, lead_ref = info.get((chrom, snp)), info.get((chrom, lead))
        effect, lead_effect = effects[pair].get(snp), effects[pair].get(lead)
        r = edges.get((pair, snp, lead))
        beta1 = beta2 = beta1_lead = beta2_lead = None
        status = "NO_EXACT_REFERENCE_MATCH" if row["reference_status"] != "EXACT_MATCH" else "NO_EXACT_REFERENCE_OR_LEAD"
        if row["reference_status"] == "EXACT_MATCH" and ref and lead_ref and effect and lead_effect:
            if int(ref["POS"]) != bp or int(effect["BP"]) != bp:
                raise ValueError(f"Reference/effect coordinate mismatch: {pair}/{snp}")
            sign = conservative_sign(ref["A1"], ref["A2"], effect["A1"], effect["A2"])
            lead_sign = conservative_sign(lead_ref["A1"], lead_ref["A2"],
                                          lead_effect["A1"], lead_effect["A2"])
            if sign is None or lead_sign is None:
                status = "AMBIGUOUS_OR_UNALIGNABLE_ALLELES"
            elif r is None:
                status = "MISSING_SIGNED_LD_EDGE"
            else:
                beta1, beta2 = float(effect["BETA1"]) * sign, float(effect["BETA2"]) * sign
                beta1_lead, beta2_lead = float(lead_effect["BETA1"]) * lead_sign, float(lead_effect["BETA2"]) * lead_sign
                status = "OUT_OF_RANGE_RAW_R" if not math.isfinite(r) or abs(r) > 1 else "PASS_SIGNED_LD"
        valid = status == "PASS_SIGNED_LD"
        output.append({
            "pair_id": pair, "locus_id": row["locus_id"], "SNP": row["SNP"], "lead_SNP": row["lead_SNP"],
            "CHR": chrom, "BP": bp, "candidate_status": row["candidate_status"],
            "reference_A1": ref["A1"] if ref else "NA", "reference_A2": ref["A2"] if ref else "NA",
            "beta_sleep_reference_A1": format(beta1, ".12g") if beta1 is not None else "NA",
            "beta_disorder_reference_A1": format(beta2, ".12g") if beta2 is not None else "NA",
            "raw_r_reference_A1": format(r, ".12g") if r is not None else "NA",
            "signed_r_reference_A1": format(r, ".12g") if valid else "NA",
            "r2_reference": format(r * r, ".12g") if valid else "NA",
            "sleep_effect_vs_lead_ld": ld_direction_consistency(beta1, beta1_lead, r) if valid else "NOT_ASSESSED",
            "disorder_effect_vs_lead_ld": ld_direction_consistency(beta2, beta2_lead, r) if valid else "NOT_ASSESSED",
            "ld_qc_status": status,
            "analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA",
        })
    OUT.mkdir(parents=True, exist_ok=False)
    target = OUT / "candidate_signed_ld.tsv"
    write(target, output, list(output[0]))
    counts = Counter(r["ld_qc_status"] for r in output)
    receipt = {
        "analysis_id": "brain6_fourpair_current_candidates_exploratory_signed_ld_v1",
        "analysis_label": "EXPLORATORY", "candidate_rows": len(output),
        "exact_reference_rows": len(exact), "status_counts": dict(counts),
        "source_sha256": {
            str(SOURCE.relative_to(ROOT)): sha(SOURCE),
            str(EXTRACTOR.relative_to(ROOT)): sha(EXTRACTOR),
            "reference_provenance": sha(REFERENCE_RECEIPT),
            "pinned_rscript": sha(RSCRIPT),
            **{f"{pair}_joined_pair_receipt": pair_sources[pair]["pair_receipt_sha256"] for pair in PAIRS},
            **{f"{pair}_joined_pair_tsv": pair_sources[pair]["pair_tsv_sha256"] for pair in PAIRS},
        },
        "reference_sample_size": 100000,
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {target.name: sha(target)},
        "limitations": [
            "Exploratory post-selection effect/LD sign comparison, not independent replication or fine-mapping.",
            "Palindromic or unalignable alleles are not assigned direction.",
            "Raw R outside [-1,1] is retained as QC failure and never clipped.",
            "No LAVA local-rg or final region tier is promoted.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"rows": len(output), "statuses": dict(counts)}, sort_keys=True))


if __name__ == "__main__":
    main()
