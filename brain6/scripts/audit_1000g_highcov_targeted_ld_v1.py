#!/usr/bin/env python3
"""Targeted 30x GRCh38 genotype follow-up for Phase 3 gaps and chr17 LD conflict."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pysam
from pyliftover import LiftOver

from audit_1000g_phase3_eur_ld_v1 import ROOT, dosage, pearson_ld, rows, sha256, write_rows


REFDIR = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1")
CHAIN = REFDIR / "hg19ToHg38.over.chain.gz"
PANEL = REFDIR / "integrated_call_samples_v3.20130502.ALL.panel"
CANDIDATES = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv"
SIGNED = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
PHASE3 = ROOT / "brain6/results/loci/independent_1000g_phase3_eur_ld_v1/candidate_lead_ld.tsv"
OUT = ROOT / "brain6/results/loci/independent_1000g_highcov_targeted_ld_v1"
URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def main() -> None:
    selected = [x for x in rows(CANDIDATES) if x["locus_id"] in {
        "placo_insomnia__mdd_chr5_103481726_104481726",
        "placo_longsleep__scz_chr17_43325576_45343136"}]
    signed = {(x["pair_id"], x["SNP"]): x for x in rows(SIGNED)}
    previous = {(x["pair_id"], x["SNP"]): x for x in rows(PHASE3)}
    samples = [x["sample"] for x in rows(PANEL) if x["super_pop"] == "EUR"]
    lifter = LiftOver(str(CHAIN))
    targets: dict[str, dict[int, list[dict]]] = defaultdict(lambda: defaultdict(list))
    liftover_rows = []
    for candidate in selected:
        mapped = lifter.convert_coordinate(f"chr{candidate['CHR']}", int(candidate["BP"]) - 1)
        valid = [(chrom, pos + 1, strand) for chrom, pos, strand, _ in mapped
                 if chrom == f"chr{candidate['CHR']}" and strand == "+"]
        if len(valid) == 1:
            targets[candidate["CHR"]][valid[0][1]].append(candidate)
        liftover_rows.append({"pair_id": candidate["pair_id"], "SNP": candidate["SNP"],
                              "CHR37": candidate["CHR"], "BP37": candidate["BP"],
                              "CHR38": valid[0][0] if len(valid) == 1 else "",
                              "BP38": valid[0][1] if len(valid) == 1 else "",
                              "liftover_status": "UNIQUE_PLUS" if len(valid) == 1 else "UNRESOLVED"})
    matches, genotypes, scans = {}, {}, []
    for chrom in sorted(targets, key=int):
        url = URL.format(chrom=chrom)
        vcf = pysam.VariantFile(url)
        missing_samples = set(samples) - set(vcf.header.samples)
        if missing_samples:
            raise ValueError(f"Missing {len(missing_samples)} Phase3 EUR samples in high coverage VCF")
        vcf.subset_samples(samples)
        positions = targets[chrom]
        start, stop = min(positions), max(positions)
        n_records = 0
        for record in vcf.fetch(f"chr{chrom}", start - 1, stop):
            n_records += 1
            if record.pos not in positions or len(record.alts or ()) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                continue
            for candidate in positions[record.pos]:
                key = candidate["pair_id"], candidate["SNP"]
                ref = signed.get(key)
                if ref is None or {record.ref, record.alts[0]} != {ref["reference_A1"], ref["reference_A2"]}:
                    continue
                if key in matches:
                    raise ValueError(f"Duplicate exact allele high coverage record: {key}")
                gt = dosage(record, samples)
                if ref["reference_A1"] == record.ref:
                    gt = 2 - gt
                genotypes[key] = gt
                matches[key] = {"vcf_id": record.id or "", "vcf_ref": record.ref,
                                "vcf_alt": record.alts[0], "n_nonmissing": int(np.isfinite(gt).sum())}
        scans.append({"chrom": chrom, "start": start, "stop": stop,
                      "records_scanned": n_records, "vcf_url": url})
        vcf.close()
        print(f"chr{chrom}:{start}-{stop} scanned {n_records}, matched {sum(k[0].startswith('insomnia') if chrom == '5' else k[0].startswith('longsleep') for k in matches)}", flush=True)

    output = []
    for candidate in selected:
        key, lead = (candidate["pair_id"], candidate["SNP"]), (candidate["pair_id"], candidate["lead_SNP"])
        lift = next(x for x in liftover_rows if x["pair_id"] == key[0] and x["SNP"] == key[1])
        ld = (math.nan,) * 5
        if key in genotypes and lead in genotypes:
            ld = pearson_ld(genotypes[key], genotypes[lead])
        prior = previous.get(key)
        output.append({"locus_id": candidate["locus_id"], "pair_id": key[0], "SNP": key[1],
                       "candidate_status": candidate["candidate_status"], "lead_SNP": lead[1],
                       "BP37": candidate["BP"], "BP38": lift["BP38"],
                       "liftover_status": lift["liftover_status"],
                       "highcov_exact_alleles": int(key in genotypes),
                       "highcov_lead_exact_alleles": int(lead in genotypes),
                       "highcov_vcf_id": matches.get(key, {}).get("vcf_id", ""),
                       "highcov_n": ld[2], "highcov_r": ld[0], "highcov_r2": ld[1],
                       "phase3_r2": prior["eur_r2"] if prior else "",
                       "lava_r2": candidate["r2_to_lead"]})
    OUT.mkdir(parents=True, exist_ok=False)
    write_rows(OUT / "target_variant_liftover.tsv", liftover_rows, list(liftover_rows[0]))
    write_rows(OUT / "targeted_highcov_ld.tsv", output, list(output[0]))
    receipt = {"analysis_id": "brain6_1000g_highcov_targeted_ld_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "sample_selection": "503 Phase 3 EUR participants, matched by sample ID within the 3202 high coverage call set",
               "method": "Unique plus-strand UCSC hg19ToHg38 chain mapping; exact unordered SNP alleles; unphased dosage Pearson r",
               "limitations": ["The high coverage and Phase 3 panels reuse the same 503 individuals; they are not independent replication cohorts",
                               "Targeted to chr5 missing lead and chr17 LD discrepancy, not all 19 regions"],
               "scans": scans,
               "inputs_sha256": {str(x): sha256(x) for x in (CHAIN, PANEL, CANDIDATES, SIGNED, PHASE3, Path(__file__))},
               "outputs_sha256": {name: sha256(OUT / name) for name in ("target_variant_liftover.tsv", "targeted_highcov_ld.tsv")}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for lid in sorted({x["locus_id"] for x in output}):
        region = [x for x in output if x["locus_id"] == lid]
        measured = [x for x in region if math.isfinite(x["highcov_r2"])]
        print(f"{lid}: {sum(x['highcov_exact_alleles'] for x in region)}/{len(region)} variants, "
              f"{sum(x['candidate_status']=='LEAD' and x['highcov_exact_alleles'] for x in region)} lead(s), "
              f"{len(measured)} measured LD rows", flush=True)


if __name__ == "__main__":
    main()
