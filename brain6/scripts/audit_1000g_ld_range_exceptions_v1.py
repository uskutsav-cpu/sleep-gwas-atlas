#!/usr/bin/env python3
"""Test LAVA's out-of-range candidate LD edges against independent EUR genotypes."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pysam

from audit_1000g_phase3_eur_ld_v1 import BASE_URL, VCF_NAME, ROOT, dosage, merged_intervals, pearson_ld, rows, sha256, write_rows


SOURCE = ROOT / "brain6/results/loci/placo_candidate_ld_range_exceptions_partial_v2.tsv"
SIGNED = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
LOCI = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_loci.tsv"
PANEL = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/integrated_call_samples_v3.20130502.ALL.panel")
OUT = ROOT / "brain6/results/loci/independent_1000g_ld_range_exceptions_v1"


def main() -> None:
    exceptions = rows(SOURCE)
    signed = {(x["pair_id"], x["SNP"]): x for x in rows(SIGNED)}
    panel = rows(PANEL)
    samples = [x["sample"] for x in panel if x["super_pop"] == "EUR"]
    regions = rows(LOCI)
    by_chr = defaultdict(list)
    for region in regions:
        by_chr[region["CHR"]].append(region)
    targets = defaultdict(list)
    for edge in exceptions:
        for suffix in ("1", "2"):
            key = edge["pair_id"], edge[f"SNP{suffix}"]
            targets[(edge["CHR"], int(edge[f"BP{suffix}"]))].append(key)
    genotypes = {}
    scans = []
    for chrom in sorted(by_chr, key=int):
        if not any(key[0] == chrom for key in targets):
            continue
        url = f"{BASE_URL}/{VCF_NAME.format(chrom=chrom)}"
        vcf = pysam.VariantFile(url)
        vcf.subset_samples(samples)
        for start, stop in merged_intervals(by_chr[chrom]):
            count = 0
            for record in vcf.fetch(chrom, start - 1, stop):
                count += 1
                keys = targets.get((chrom, record.pos), [])
                if not keys or len(record.alts or ()) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                    continue
                for key in keys:
                    if key in genotypes:
                        continue
                    ref = signed[key]
                    if {record.ref, record.alts[0]} != {ref["reference_A1"], ref["reference_A2"]}:
                        continue
                    value = dosage(record, samples)
                    if ref["reference_A1"] == record.ref:
                        value = 2 - value
                    genotypes[key] = value
            scans.append({"chrom": chrom, "start": start, "stop": stop, "records_scanned": count, "vcf_url": url})
            print(f"chr{chrom}:{start}-{stop}: {count} VCF records", flush=True)
        vcf.close()
    output = []
    for edge in exceptions:
        a, b = (edge["pair_id"], edge["SNP1"]), (edge["pair_id"], edge["SNP2"])
        ld = (math.nan,) * 5
        if a in genotypes and b in genotypes:
            ld = pearson_ld(genotypes[a], genotypes[b])
        output.append({**edge, "snp1_exact_alleles": int(a in genotypes),
                       "snp2_exact_alleles": int(b in genotypes),
                       "eur_r": ld[0], "eur_r2": ld[1], "eur_n": ld[2]})
    measured = [x for x in output if math.isfinite(x["eur_r2"])]
    OUT.mkdir(parents=True, exist_ok=False)
    write_rows(OUT / "invalid_raw_edges_1000g_ld.tsv", output, list(output[0]))
    stats = {"n_raw_invalid": len(output), "n_1000g_measured": len(measured),
             "n_r2_below_0_1": sum(x["eur_r2"] < 0.1 for x in measured),
             "n_r2_below_0_9": sum(x["eur_r2"] < 0.9 for x in measured),
             "median_1000g_r2": statistics.median(x["eur_r2"] for x in measured),
             "min_1000g_r2": min(x["eur_r2"] for x in measured),
             "max_1000g_r2": max(x["eur_r2"] for x in measured)}
    receipt = {"analysis_id": "brain6_1000g_ld_range_exceptions_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(), "stats": stats,
               "method": "Unphased allele-dosage Pearson r in 503 1000 Genomes Phase 3 EUR participants; exact GRCh37 position and unordered alleles",
               "limitations": ["This independent panel does not repair the frozen LAVA UKB factor matrices",
                               "Out-of-range LAVA factor products are invalid as literal correlations"],
               "vcf_scans": scans,
               "inputs_sha256": {str(x): sha256(x) for x in (SOURCE, SIGNED, LOCI, PANEL, Path(__file__))},
               "output_sha256": sha256(OUT / "invalid_raw_edges_1000g_ld.tsv")}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
