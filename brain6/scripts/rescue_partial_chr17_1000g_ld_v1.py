#!/usr/bin/env python3
"""Test Phase 3 rescue of chr17 candidates absent from 1000G 30x genotypes.

This is a partial-family LD diagnostic. It neither changes the v1 clumping
artifact nor assigns final Brain6 loci or tiers.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pysam

from audit_1000g_phase3_eur_ld_v1 import BASE_URL, VCF_NAME, dosage, pearson_ld, rows, sha256, write_rows


ROOT = Path(__file__).resolve().parents[2]
ASSIGNMENTS = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v1/variant_assignments.tsv"
LEADS = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v1/diagnostic_leads.tsv"
SIGNED = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
LIFTOVER = ROOT / "brain6/results/loci/independent_1000g_highcov_targeted_ld_v1/target_variant_liftover.tsv"
PANEL = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/integrated_call_samples_v3.20130502.ALL.panel")
OUT = ROOT / "brain6/results/loci/partial_chr17_phase3_rescue_v1"
PHASE3_URL = f"{BASE_URL}/{VCF_NAME.format(chrom='17')}"
HIGHCOV_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def capture(vcf: pysam.VariantFile, chrom: str, targets: dict[int, list[dict]],
            signed: dict[tuple[str, str], dict], samples: list[str]) -> tuple[dict[str, np.ndarray], int]:
    found: dict[str, np.ndarray] = {}
    scanned = 0
    for bp in sorted(targets):
        for record in vcf.fetch(chrom, bp - 1, bp):
            scanned += 1
            if record.pos != bp or len(record.alts or ()) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                continue
            for target in targets[bp]:
                key = (target["pair_id"], target["SNP"])
                expected = {signed[key]["reference_A1"], signed[key]["reference_A2"]}
                if {record.ref, record.alts[0]} != expected:
                    continue
                if target["SNP"] in found:
                    raise ValueError(f"duplicate exact allele record for {target['SNP']}")
                gt = dosage(record, samples)
                if signed[key]["reference_A1"] == record.ref:
                    gt = 2 - gt
                found[target["SNP"]] = gt
    return found, scanned


def main() -> None:
    assignments = rows(ASSIGNMENTS)
    missing = [x for x in assignments if x["CHR"] == "17" and x["status"] == "NOT_CLUMPED_NO_EXACT_GENOTYPE"]
    leads = [x for x in rows(LEADS) if x["pair_id"] == "longsleep__scz" and x["CHR"] == "17"]
    if len(missing) != 31 or len(leads) != 3:
        raise ValueError(f"unexpected fixed partial-family chr17 scope: {len(missing)} missing, {len(leads)} leads")
    if any(x["pair_id"] != "longsleep__scz" for x in missing):
        raise ValueError("chr17 missing candidate pair changed")
    signed = {(x["pair_id"], x["SNP"]): x for x in rows(SIGNED)}
    samples = [x["sample"] for x in rows(PANEL) if x["super_pop"] == "EUR"]
    if len(samples) != 503 or len(set(samples)) != 503:
        raise ValueError("EUR panel changed")
    targets: dict[int, list[dict]] = {}
    for x in missing:
        targets.setdefault(int(x["BP"]), []).append(x)
    for x in leads:
        targets.setdefault(int(x["lead_BP"]), []).append({"pair_id": x["pair_id"], "SNP": x["lead_SNP"]})
    with pysam.VariantFile(PHASE3_URL) as vcf:
        if set(samples) - set(vcf.header.samples):
            raise ValueError("Phase 3 VCF lacks EUR samples")
        vcf.subset_samples(samples)
        phase3, phase3_scanned = capture(vcf, "17", targets, signed, samples)
    if not all(lead["lead_SNP"] in phase3 for lead in leads):
        raise ValueError("one or more Phase 3 lead genotypes are absent")

    lifted = {(x["pair_id"], x["SNP"]): x for x in rows(LIFTOVER)}
    high_targets: dict[int, list[dict]] = {}
    for lead in leads:
        x = lifted[(lead["pair_id"], lead["lead_SNP"])]
        if x["liftover_status"] != "UNIQUE_PLUS" or x["CHR38"] != "chr17":
            raise ValueError("high coverage lead liftOver changed")
        high_targets.setdefault(int(x["BP38"]), []).append({"pair_id": lead["pair_id"], "SNP": lead["lead_SNP"]})
    with pysam.VariantFile(HIGHCOV_URL) as vcf:
        if set(samples) - set(vcf.header.samples):
            raise ValueError("30x VCF lacks EUR samples")
        vcf.subset_samples(samples)
        highcov, highcov_scanned = capture(vcf, "chr17", high_targets, signed, samples)
    if len(highcov) != len(leads):
        raise ValueError("30x lead genotype absent")

    concordance = []
    for lead in leads:
        snp = lead["lead_SNP"]
        a, b = phase3[snp], highcov[snp]
        valid = np.isfinite(a) & np.isfinite(b)
        r, r2, n, _, _ = pearson_ld(a, b)
        concordance.append({"lead_SNP": snp, "n_same_sample_nonmissing": n,
                            "discordant_genotype_calls": int(np.sum(a[valid] != b[valid])),
                            "cross_panel_r2": r2})
    ordered_leads = sorted(leads, key=lambda x: float(x["lead_P_PLACO"]))
    rescue = []
    for x in sorted(missing, key=lambda y: float(y["P_PLACO"])):
        if x["SNP"] not in phase3:
            rescue.append({"pair_id": x["pair_id"], "CHR": x["CHR"], "BP": x["BP"],
                           "SNP": x["SNP"], "P_PLACO": x["P_PLACO"],
                           "phase3_exact_genotype": "NO", "earlier_leads_within_500kb": "",
                           "max_r2_to_earlier_30x_lead": "", "best_earlier_lead": "",
                           "n_samples": "", "status": "UNRESOLVED_NO_EXACT_GENOTYPE"})
            continue
        eligible = [lead for lead in ordered_leads if float(lead["lead_P_PLACO"]) < float(x["P_PLACO"])
                    and abs(int(lead["lead_BP"]) - int(x["BP"])) <= 500_000]
        edges = []
        for lead in eligible:
            r, r2, n, _, _ = pearson_ld(phase3[x["SNP"]], highcov[lead["lead_SNP"]])
            edges.append((r2, lead["lead_SNP"], n))
        best = max((e for e in edges if math.isfinite(e[0])), default=None)
        rescue.append({"pair_id": x["pair_id"], "CHR": x["CHR"], "BP": x["BP"],
                       "SNP": x["SNP"], "P_PLACO": x["P_PLACO"],
                       "phase3_exact_genotype": "YES", "earlier_leads_within_500kb": len(eligible),
                       "max_r2_to_earlier_30x_lead": best[0] if best else "",
                       "best_earlier_lead": best[1] if best else "",
                       "n_samples": best[2] if best else "",
                       "status": "CLUMPED_TO_EXISTING_LEAD" if best and best[0] >= 0.1 else "NEW_LEAD_POTENTIAL_REQUIRES_FULL_GREEDY_REBUILD"})
    OUT.mkdir(parents=True, exist_ok=False)
    write_rows(OUT / "rescue_assignments.tsv", rescue, list(rescue[0]))
    write_rows(OUT / "lead_cross_panel_concordance.tsv", concordance, list(concordance[0]))
    provenance = {"analysis_id": "brain6_partial_chr17_phase3_rescue_v1",
                  "created_utc": datetime.now(timezone.utc).isoformat(),
                  "scope": "Partial-family diagnostic only; no frozen PLACO, LAVA, locus, or tier change",
                  "phase3_url": PHASE3_URL, "highcov_url": HIGHCOV_URL,
                  "reference_records_scanned": {"phase3": phase3_scanned, "highcov": highcov_scanned},
                  "counts": {"missing_30x_candidates": len(rescue),
                             "phase3_exact_genotypes": sum(x["phase3_exact_genotype"] == "YES" for x in rescue),
                             "unresolved_no_exact_genotype": sum(x["status"] == "UNRESOLVED_NO_EXACT_GENOTYPE" for x in rescue),
                             "clumped_to_existing_lead": sum(x["status"] == "CLUMPED_TO_EXISTING_LEAD" for x in rescue),
                             "new_lead_potential": sum(x["status"] == "NEW_LEAD_POTENTIAL_REQUIRES_FULL_GREEDY_REBUILD" for x in rescue)},
                  "inputs_sha256": {str(p): sha256(p) for p in (ASSIGNMENTS, LEADS, SIGNED, LIFTOVER, PANEL, Path(__file__))},
                  "outputs_sha256": {p.name: sha256(p) for p in (OUT / "rescue_assignments.tsv", OUT / "lead_cross_panel_concordance.tsv")}}
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance["counts"]), flush=True)


if __name__ == "__main__":
    main()
