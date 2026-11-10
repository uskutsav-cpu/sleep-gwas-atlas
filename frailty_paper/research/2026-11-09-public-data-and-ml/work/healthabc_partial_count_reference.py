#!/usr/bin/env python3
"""Compare typed-array subset vs complete count in public reference genotypes.

The public reference panel is not Health ABC and contains no project outcomes.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pysam

root = Path(__file__).resolve().parent.parent
out = root / "outputs"
score_path = out / "unvalidated_22_shared_direction_score.csv"
coverage_path = out / "healthabc_geo_platform_coordinate_screen.json"
panel_path = root / "work/1kg_phase3_panel.tsv"
score = list(csv.DictReader(score_path.open()))
screen = json.loads(coverage_path.read_text())
assert screen["complete_platform_table"] and len(score) == 22
subset = {x["SNP_ID"] for x in screen["matches"]} & {x["variant"] for x in score}
assert len(subset) == 8
panel = list(csv.DictReader(panel_path.open(), delimiter="\t"))
groups = {group: sorted(r["sample"] for r in panel if r["super_pop"] == group) for group in ("EUR", "AFR")}
assert {k: len(v) for k, v in groups.items()} == {"EUR": 503, "AFR": 661}
handles, columns = {}, {group: [] for group in groups}
for item in score:
    chrom, pos = item["chr_GRCh37_assumed"], int(item["position_GRCh37_assumed"])
    if chrom not in handles:
        handles[chrom] = pysam.VariantFile(f"https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
    f = handles[chrom]
    hits = [v for v in f.fetch(chrom, pos-1, pos) if v.pos == pos and v.id == item["variant"] and len(v.alleles) == 2 and set(v.alleles) == {item["counted_allele"], item["other_allele"]}]
    if len(hits) != 1:
        raise ValueError("Reference record not uniquely biallelic: " + item["variant"])
    v = hits[0]; index = v.alleles.index(item["counted_allele"])
    for group, samples in groups.items():
        genotypes = []
        for sample in samples:
            gt = v.samples[sample].get("GT")
            if gt is None or len(gt) != 2 or any(a is None for a in gt):
                raise ValueError("Missing or non-diploid reference genotype")
            genotypes.append(sum(a == index for a in gt))
        columns[group].append(genotypes)
mask = np.array([r["variant"] in subset for r in score])
results = {}
for group, vectors in columns.items():
    x = np.asarray(vectors, dtype=float).T
    complete, partial, omitted = x.sum(axis=1), x[:,mask].sum(axis=1), x[:,~mask].sum(axis=1)
    results[group] = {"n":len(complete), "correlation":float(np.corrcoef(complete,partial)[0,1]),
                      "full_sd":float(complete.std(ddof=1)), "subset_sd":float(partial.std(ddof=1)),
                      "omitted_sd":float(omitted.std(ddof=1))}
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {"groups":results, "subset_rsids":sorted(subset),
           "input_sha256":{"score":digest(score_path),"coverage":digest(coverage_path),"panel":digest(panel_path)},
           "limitations":"Genotype-only reference correlations; no Health ABC data, SNP QC, imputation, outcome, independent score validation or clinical prediction. The remote VCF full files were not checksum verified."}
(out / "healthabc_partial_count_reference.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
