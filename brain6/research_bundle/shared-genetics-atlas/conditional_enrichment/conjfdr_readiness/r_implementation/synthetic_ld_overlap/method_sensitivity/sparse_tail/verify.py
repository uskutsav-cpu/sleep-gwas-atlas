#!/usr/bin/env python3
import hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
d=json.loads((p/'results.json').read_text());t=d['exact_tail_counts']
assert d['snp']=='s021341' and len(t)==5 and all(x['retained']==4400 for x in t)
assert all(x['both']<=min(x['trait1_tail'],x['trait2_tail']) for x in t)
assert math.isclose(d['outputs']['TRUE']['conj'],max(d['outputs']['TRUE']['cfdr12'],d['outputs']['TRUE']['cfdr21']),rel_tol=1e-14)
assert math.isclose(d['outputs']['FALSE']['conj'],max(d['outputs']['FALSE']['cfdr12'],d['outputs']['FALSE']['cfdr21']),rel_tol=1e-14)
assert 'SPARSE_TAIL_COMPLETED' in (p/'run.log').read_text()
files=[p/x for x in ('PROTOCOL.md','run.R','run.log','results.json')]
receipt={'sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
 'trait2_tail_range':[min(x['trait2_tail'] for x in t),max(x['trait2_tail'] for x in t)],
 'joint_tail_counts':[x['both'] for x in t], 'target_retained':[x['target_retained'] for x in t],
 'note':'Exact corrected p thresholds, not package grid bin counts; selected case, not calibrated inference.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
