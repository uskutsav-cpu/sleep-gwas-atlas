#!/usr/bin/env python3
"""Diagnostic 1000G EUR clumping of four-pair PLACO candidates.

This intentionally cannot assign a final locus tier while protected Track B is
absent. It applies the frozen 500 kb/r² 0.1 rule to exact-allele public genotypes.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pysam
from pyliftover import LiftOver

from audit_1000g_phase3_eur_ld_v1 import BASE_URL, VCF_NAME, ROOT, dosage, pearson_ld, rows, sha256, write_rows


SOURCE = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv"
SIGNED = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
FROZEN = ROOT / "brain6/config/shared_locus_rule_v1.json"
PANEL = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/integrated_call_samples_v3.20130502.ALL.panel")
CHAIN = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz")
OUT = ROOT / "brain6/results/loci/partial_candidates_1000g_clumping_v1"
HIGHCOV_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def main() -> None:
    frozen = json.loads(FROZEN.read_text())
    assert frozen["clumping"]["r2_threshold"] == 0.1
    assert frozen["clumping"]["window_kb"] == 500
    assert frozen["region_definition"]["lead_flank_kb"] == 500
    candidates = rows(SOURCE)
    signed = {(x["pair_id"], x["SNP"]): x for x in rows(SIGNED)}
    samples = [x["sample"] for x in rows(PANEL) if x["super_pop"] == "EUR"]
    assert len(samples) == 503
    if len(candidates) != 2297:
        raise ValueError(f"Unexpected partial candidate count: {len(candidates)}")
    lifter = LiftOver(str(CHAIN))
    target = defaultdict(lambda: defaultdict(list))
    mapping = {}
    for row in candidates:
        key = row["pair_id"], row["SNP"]
        if row["reference_status"] != "EXACT_MATCH":
            mapping[key] = ("UNMATCHED_LOCKED_REFERENCE", "", "")
            continue
        chrom, bp = row["CHR"], int(row["BP"])
        if chrom in ("5", "17"):
            lifted = lifter.convert_coordinate(f"chr{chrom}", bp - 1)
            lifted = [(c, p + 1) for c, p, strand, _ in lifted if c == f"chr{chrom}" and strand == "+"]
            if len(lifted) != 1:
                mapping[key] = ("LIFTOVER_UNRESOLVED", "", "")
                continue
            target[chrom][lifted[0][1]].append(row)
            mapping[key] = ("HIGHCOV_GRCH38", lifted[0][1], "")
        else:
            target[chrom][bp].append(row)
            mapping[key] = ("PHASE3_GRCH37", bp, "")

    genotypes, scans = {}, []
    for chrom in sorted(target, key=int):
        is_highcov = chrom in ("5", "17")
        url = HIGHCOV_URL.format(chrom=chrom) if is_highcov else f"{BASE_URL}/{VCF_NAME.format(chrom=chrom)}"
        vcf = pysam.VariantFile(url)
        if set(samples) - set(vcf.header.samples):
            raise ValueError(f"EUR samples absent from chr{chrom} reference")
        vcf.subset_samples(samples)
        positions = target[chrom]
        intervals = []
        for position in sorted(positions):
            if intervals and position - intervals[-1][1] <= 1_000_000:
                intervals[-1][1] = position
            else:
                intervals.append([position, position])
        for start, stop in intervals:
            n_records = 0
            for record in vcf.fetch(f"chr{chrom}" if is_highcov else chrom, start - 1, stop):
                n_records += 1
                if record.pos not in positions or len(record.alts or ()) != 1:
                    continue
                if len(record.ref) != 1 or len(record.alts[0]) != 1:
                    continue
                for row in positions[record.pos]:
                    key = row["pair_id"], row["SNP"]
                    ref = signed.get(key)
                    if ref is None or {record.ref, record.alts[0]} != {ref["reference_A1"], ref["reference_A2"]}:
                        continue
                    if key in genotypes:
                        raise ValueError(f"Duplicate exact allele reference record for {key}")
                    gt = dosage(record, samples)
                    if ref["reference_A1"] == record.ref:
                        gt = 2 - gt
                    genotypes[key] = gt
            scans.append({"chrom": chrom, "start": start, "stop": stop,
                          "reference": "30x GRCh38" if is_highcov else "Phase3 GRCh37",
                          "url": url, "vcf_records_scanned": n_records})
            print(f"chr{chrom}:{start}-{stop} {scans[-1]['reference']}: {n_records} VCF records", flush=True)
        vcf.close()

    by_pair_chr = defaultdict(list)
    for row in candidates:
        by_pair_chr[(row["pair_id"], row["CHR"])].append(row)
    assignments, leads, regions = [], [], []
    for (pair, chrom), group in sorted(by_pair_chr.items(), key=lambda x: (x[0][0], int(x[0][1]))):
        ordered = sorted(group, key=lambda x: (float(x["P_PLACO"]), int(x["BP"]), x["SNP"]))
        unassigned = {(x["pair_id"], x["SNP"]) for x in ordered if (x["pair_id"], x["SNP"]) in genotypes}
        for row in ordered:
            key = row["pair_id"], row["SNP"]
            if key not in unassigned:
                continue
            unassigned.remove(key)
            leads.append({"pair_id": pair, "CHR": chrom, "lead_SNP": row["SNP"],
                          "lead_BP": row["BP"], "lead_P_PLACO": row["P_PLACO"],
                          "reference": mapping[key][0]})
            assignments.append({"pair_id": pair, "CHR": chrom, "SNP": row["SNP"],
                                "BP": row["BP"], "P_PLACO": row["P_PLACO"],
                                "status": "LEAD", "lead_SNP": row["SNP"], "r2_to_lead": 1.0,
                                "reference": mapping[key][0]})
            for other in ordered:
                other_key = other["pair_id"], other["SNP"]
                if other_key not in unassigned or abs(int(other["BP"]) - int(row["BP"])) > 500_000:
                    continue
                r, r2, n, _, _ = pearson_ld(genotypes[key], genotypes[other_key])
                if math.isfinite(r2) and r2 >= 0.1:
                    unassigned.remove(other_key)
                    assignments.append({"pair_id": pair, "CHR": chrom, "SNP": other["SNP"],
                                        "BP": other["BP"], "P_PLACO": other["P_PLACO"],
                                        "status": "LD_CLUMPED", "lead_SNP": row["SNP"],
                                        "r2_to_lead": r2, "reference": mapping[other_key][0]})
        for row in ordered:
            key = row["pair_id"], row["SNP"]
            if key not in genotypes:
                assignments.append({"pair_id": pair, "CHR": chrom, "SNP": row["SNP"],
                                    "BP": row["BP"], "P_PLACO": row["P_PLACO"],
                                    "status": "NOT_CLUMPED_NO_EXACT_GENOTYPE", "lead_SNP": "",
                                    "r2_to_lead": "", "reference": mapping[key][0]})
        lead_group = sorted([x for x in leads if x["pair_id"] == pair and x["CHR"] == chrom], key=lambda x: int(x["lead_BP"]))
        for lead in lead_group:
            start, stop = int(lead["lead_BP"]) - 500_000, int(lead["lead_BP"]) + 500_000
            if regions and regions[-1]["pair_id"] == pair and regions[-1]["CHR"] == chrom and start <= int(regions[-1]["STOP"]):
                regions[-1]["STOP"] = max(int(regions[-1]["STOP"]), stop)
                regions[-1]["lead_variants"] += ";" + lead["lead_SNP"]
                regions[-1]["n_leads"] += 1
            else:
                regions.append({"pair_id": pair, "CHR": chrom, "START": start, "STOP": stop,
                                "lead_variants": lead["lead_SNP"], "n_leads": 1})

    assert len(assignments) == len(candidates)
    OUT.mkdir(parents=True, exist_ok=False)
    for name, values in (("variant_assignments.tsv", assignments), ("diagnostic_leads.tsv", leads),
                         ("diagnostic_regions.tsv", regions)):
        write_rows(OUT / name, values, list(values[0]))
    receipt = {"analysis_id": "brain6_partial_candidates_1000g_clumping_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "scope": "Diagnostic four-pair candidates only; protected Track B missing; no final loci or tiers",
               "frozen_rule": "Greedy P order, 500 kb, r2>=0.1 clumped, 500 kb lead flanks, overlapping intervals merged within pair",
               "reference_selection": "1000 Genomes 30x GRCh38 for chr5/chr17 with unique plus-strand liftOver; Phase3 GRCh37 for other chromosomes; same 503 EUR IDs",
               "scans": scans,
               "counts": {"input_variants": len(candidates), "exact_genotypes": len(genotypes),
                          "unmatched": len(candidates) - len(genotypes), "diagnostic_leads": len(leads),
                          "diagnostic_regions": len(regions)},
               "inputs_sha256": {str(x): sha256(x) for x in (SOURCE, SIGNED, FROZEN, PANEL, CHAIN, Path(__file__))},
               "outputs_sha256": {name: sha256(OUT / name) for name in
                                  ("variant_assignments.tsv", "diagnostic_leads.tsv", "diagnostic_regions.tsv")}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt["counts"], indent=2))


if __name__ == "__main__":
    main()
