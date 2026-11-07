#!/usr/bin/env python3
"""Join frozen, previously extracted candidate tables; no source inference."""
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parent
ukb = Path('/Volumes/Extreme SSD/brain6-hfrs-candidate-results-20261028.tsv')
fin = Path('/Volumes/Extreme SSD/brain6-mvp-finngen-hfrs-results-20261028.tsv')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

u = [r for r in rows(ukb) if r['score'] == 'original']
f = rows(fin)
assert len(u) == len(f) == 425
assert len({r['rsid'] for r in u}) == len({r['rsid'] for r in f}) == 425
f_by_id = {r['rsid']: r for r in f}
out = []
for a in u:
    b = f_by_id[a['rsid']]
    assert a['grch38'] == b['grch38'] and a['region'] == b['region']
    if a['status'] != 'matched' or b['hfrs_status'] != 'matched':
        continue
    p1, p2 = float(a['p']), float(b['hfrs_p'])
    beta1, beta2 = float(a['beta_candidate_alt']), float(b['hfrs_beta_alt'])
    assert all(map(math.isfinite, (p1, p2, beta1, beta2))) and 0 <= p1 <= 1 and 0 <= p2 <= 1
    out.append({'rsid': a['rsid'], 'region': a['region'], 'grch38': a['grch38'],
                'ukb_p': p1, 'finnish_p': p2, 'ukb_beta_alt': beta1, 'finnish_beta_alt': beta2,
                'same_direction': beta1 * beta2 > 0, 'max_p': max(p1, p2),
                'bonferroni_425': min(1, 425 * max(p1, p2))})
out.sort(key=lambda r: (r['max_p'], r['rsid']))
with (root / 'results.tsv').open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter='\t')
    w.writeheader(); w.writerows(out)
hits = [r for r in out if r['bonferroni_425'] < .05]
summary = {'input_sha256': {str(ukb): digest(ukb), str(fin): digest(fin)},
           'results_sha256': digest(root / 'results.tsv'), 'candidate_slots': 425,
           'matched': len(out), 'same_direction': sum(r['same_direction'] for r in out),
           'hits': len(hits), 'hits_same_direction': sum(r['same_direction'] for r in hits),
           'hit_counts_by_region': dict(Counter(r['region'] for r in hits)),
           'top': out[:5], 'selected': [r for r in out if r['rsid'] in
               {'rs62264764', 'rs62264767', 'rs2431108', 'rs77960', 'rs4916723', 'rs1476535'}]}
(root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
