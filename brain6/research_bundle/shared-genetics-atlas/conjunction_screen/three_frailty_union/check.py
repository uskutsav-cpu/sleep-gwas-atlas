#!/usr/bin/env python3
"""Derive a post-hoc union-family no-rejection result from saved screen receipts."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
out = Path(__file__).resolve().parent
fi_path = root / 'conjunction_screen/fdr_summary.json'
hfrs_path = root / 'hfrs_conjunction/summary.json'
fi = json.loads(fi_path.read_text())
hfrs = json.loads(hfrs_path.read_text())
assert fi['level'] == 0.05
assert fi['BH_rejections'] == 0 and hfrs['BH_rejections_padded_family'] == 0
m1, m2 = fi['number_hypotheses'], hfrs['family_size_insomnia_rows']
assert m1 == 2276524 and m2 == 1158908
minimum = min(fi['minimum_p'], hfrs['min_conjunction_p'])
assert minimum > 0.05 / (m1 + m2)
result = {
    'family_sizes': {'FI_and_FFS': m1, 'HFRS': m2, 'union': m1 + m2},
    'level': 0.05, 'bonferroni_cutoff': 0.05 / (m1 + m2),
    'minimum_reported_max_p': minimum,
    'bonferroni_rejections': 0, 'BH_rejections': 0, 'BY_rejections': 0,
    'derivation': 'BH union no-rejection theorem; BY rejection implies BH rejection',
    'input_receipt_sha256': {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (fi_path, hfrs_path)},
}
(out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
