#!/usr/bin/env python3
"""Check stored 22-marker FFS weights against cached source rows, not source semantics."""
import csv
from pathlib import Path

out = Path(__file__).resolve().parents[1] / 'outputs'
with (out / 'unvalidated_22_variant_score_prototype.csv').open(newline='') as f:
    scores = list(csv.DictReader(f))
with (out / 'published_candidate_marginal_lookup.csv').open(newline='') as f:
    source_rows = list(csv.DictReader(f))
source = {x['candidate']: x for x in source_rows}
assert len(source) == len(source_rows)
rows = []
for s in scores:
    x = source[s['variant']]
    assert x['FFS_rows'] == '1'
    assert s['chr'] == x['FFS_CHROM'] and s['position_GRCh37_assumed'] == x['FFS_POS']
    pair = (s['counted_allele'], s['other_allele'])
    reported = (x['FFS_A1'], x['FFS_A2'])
    if pair == reported:
        relation, expected = 'same', float(x['FFS_BETA'])
    elif pair == reported[::-1]:
        relation, expected = 'swapped', -float(x['FFS_BETA'])
    else:
        raise ValueError(f"Unmatched alleles: {s['variant']}")
    assert abs(float(s['FFS_beta_weight_inferred']) - expected) < 1e-10
    multi = len(set(s['reference_alleles'].split('/'))) > 2
    rows.append([s['variant'], *pair, *reported, x['FFS_BETA'],
                 s['FFS_beta_weight_inferred'], relation, s['reference_alleles'], multi])
assert len(rows) == 22 and len({r[0] for r in rows}) == 22
with (out / 'score22_ffs_alignment_audit.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['variant', 'score_counted', 'score_other', 'source_A1', 'source_A2',
                'source_beta', 'score_beta', 'alignment', 'reference_alleles', 'reference_multiallelic'])
    w.writerows(rows)
print(f"{len(rows)} markers; {sum(r[7] == 'swapped' for r in rows)} swapped; "
      f"{sum(r[7] == 'swapped' and r[9] for r in rows)} swapped with multiallelic reference annotations")
