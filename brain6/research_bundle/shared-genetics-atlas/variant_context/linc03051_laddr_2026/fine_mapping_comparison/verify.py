#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;receipt=json.loads((p/'receipt.json').read_text())
with (p/'comparison.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==12 and len({x['snp'] for x in r})==6
assert {(x['snp'],x['trait']) for x in r}=={(s,t) for s in {x['snp'] for x in r} for t in ('insomnia','adhd')}
for x in r:
 assert math.isclose(float(x['wald_p_approx']),math.erfc(abs(float(x['beta_a1'])/float(x['se']))/math.sqrt(2)),rel_tol=1e-12)
 assert 0<=float(x['conditional_pip'])<=1
assert hashlib.sha256((p/'comparison.tsv').read_bytes()).hexdigest()==receipt['comparison_sha256']
out={'rows':len(r),'unique_markers':6,'wald_arithmetic_checked':True,
 'receipt_sha256':hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest(),
 'note':'Arithmetic check only; not independent fine-mapping validation.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
