#!/usr/bin/env python3
"""Check a completed stream and summarize LDSC diagnostics, without altering output."""
import json, math, re, sys
from pathlib import Path
if len(sys.argv) != 3:
    raise SystemExit('usage: python audit_result.py STREAM_RECEIPT LDSC_LOG')
receipt = json.loads(Path(sys.argv[1]).read_text())
log = Path(sys.argv[2]).read_text()
assert receipt['source_md5'] == 'c10d4dca6efaa4dba2e39ae29add61ae'
assert receipt['counts']['retained'] > 200000
assert 'Analysis finished at ' in log
blocks = re.split(r'Computing rg for phenotype \d+/5', log)[1:]
assert len(blocks) == 4, f'expected 4 comparisons, got {len(blocks)}'
def find(pattern, text):
    match = re.search(pattern, text, re.MULTILINE)
    if not match: raise ValueError(f'missing diagnostic: {pattern}')
    result = [float(x) for x in match.groups()]
    assert all(math.isfinite(x) for x in result)
    return result
num = r'([-+0-9.eE]+)'
pair = num + r' \(' + num + r'\)'
out=[]
for b, target in zip(blocks, ['MVP', 'FI', 'FFS', 'HFRS']):
    assert f'/{target}.sumstats.gz' in b, f'wrong target order: {target}'
    counts = re.findall(r'(\d+) SNPs with valid alleles', b)
    assert len(counts) == 1
    rg, se = find(r'Genetic Correlation: '+pair, b)
    p, = find(r'^P: '+num, b)
    h = re.findall(r'Total Observed scale h2: '+pair, b)
    assert h and all(float(v[0]) > 0 for v in h)
    intercepts = re.findall(r'^Intercept: '+pair, b, re.MULTILINE)
    assert len(intercepts) >= 2
    out.append({'target': target, 'valid_snps': int(counts[0]), 'rg': rg, 'se': se, 'p': p,
                'h2_reports': h, 'intercepts_in_log_order': intercepts})
assert len({x['valid_snps'] for x in out}) == 1, 'marker counts differ across pairs; inspect sets and logs'
print(json.dumps({'stream_receipt': sys.argv[1], 'ldsc_log': sys.argv[2], 'comparisons': out,
                  'note': 'Marker-count equality does not itself prove marker-set equality; review log warnings and source receipt.'}, indent=2))
