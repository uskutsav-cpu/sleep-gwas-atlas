#!/usr/bin/env python3
"""Describe postselected BBJ count subsets in public reference samples only."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pysam

root = Path(__file__).resolve().parent.parent
out = root / 'outputs'
score_path = out / 'unvalidated_22_shared_direction_score.csv'
subset_path = out / 'stage6_bbj_insomnia_marker_coverage_sensitivity.json'
panel_path = root / 'work/1kg_phase3_panel.tsv'
with score_path.open(newline='') as f:
    score = list(csv.DictReader(f))
with panel_path.open(newline='') as f:
    panel = list(csv.DictReader(f, delimiter='\t'))
groups = {
    'JPT': sorted(r['sample'] for r in panel if r['pop'] == 'JPT'),
    'EAS': sorted(r['sample'] for r in panel if r['super_pop'] == 'EAS'),
    'EUR': sorted(r['sample'] for r in panel if r['super_pop'] == 'EUR'),
}
assert {k: len(v) for k, v in groups.items()} == {'JPT': 104, 'EAS': 504, 'EUR': 503}
prior = json.loads(subset_path.read_text())['result']
subsets = {name: set(prior['JPT_n104_' + key]['markers']) for name, key in [
    ('matched14', 'matched14'), ('selected15', 'expected_10_or_more_cases'),
    ('selected17', 'all_17_after_low_total_count_exclusion')]}
assert len(score) == 22 and [len(subsets[k]) for k in subsets] == [14, 15, 17]
handles, arrays = {}, {g: [] for g in groups}
for item in score:
    chrom, pos = item['chr_GRCh37_assumed'], int(item['position_GRCh37_assumed'])
    if chrom not in handles:
        url = f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        handles[chrom] = pysam.VariantFile(url)
    matches = [r for r in handles[chrom].fetch(chrom, pos - 1, pos)
               if r.pos == pos and r.id == item['variant']]
    if len(matches) != 1 or len(matches[0].alleles) != 2 or set(matches[0].alleles) != {item['counted_allele'], item['other_allele']}:
        raise ValueError('Ambiguous reference record: ' + item['variant'])
    rec = matches[0]
    allele_index = rec.alleles.index(item['counted_allele'])
    for group, names in groups.items():
        values = []
        for name in names:
            gt = rec.samples[name].get('GT')
            if gt is None or len(gt) != 2 or any(a is None for a in gt):
                raise ValueError('Missing or nondiploid genotype: ' + item['variant'])
            values.append(sum(a == allele_index for a in gt))
        arrays[group].append(values)

result = {}
for group, columns in arrays.items():
    a = np.asarray(columns, dtype=float).T
    full = a.sum(axis=1)
    result[group] = {'n': len(full), 'full_mean': float(full.mean()), 'full_sd': float(full.std(ddof=1)), 'subsets': {}}
    for name, ids in subsets.items():
        mask = np.array([item['variant'] in ids for item in score])
        part, omitted = a[:, mask].sum(axis=1), a[:, ~mask].sum(axis=1)
        result[group]['subsets'][name] = {
            'n_markers': int(mask.sum()), 'correlation_with_full': float(np.corrcoef(full, part)[0, 1]),
            'omitted_count_sd': float(omitted.std(ddof=1)), 'omitted_markers': [item['variant'] for item, keep in zip(score, mask) if not keep]}
receipt = {'inputs_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (score_path, subset_path, panel_path)},
           'source': '1000 Genomes phase 3 via UCSC hg19 indexed remote VCF',
           'limitations': 'Reference sample genotype correlations only, not target LD or outcome validation; JPT is part of EAS; VCF full-file checksum not verified.', 'groups': result}
(out / 'stage6_bbj_subset_reference_correlation.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(result, indent=2))
