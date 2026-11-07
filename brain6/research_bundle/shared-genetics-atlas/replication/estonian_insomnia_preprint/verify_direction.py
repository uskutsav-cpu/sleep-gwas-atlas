#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'direction_comparison.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
r=json.loads((p/'direction_receipt.json').read_text());assert len(rows)==2
for x in rows:
 for group in ('discovery','estonian'):
  b=float(x[group+'_beta_alt']);s=float(x[group+'_se'])
  assert math.isclose(float(x[group+'_ci_low']),b-1.96*s,abs_tol=1e-12)
  assert math.isclose(float(x[group+'_ci_high']),b+1.96*s,abs_tol=1e-12)
 assert float(x['discovery_beta_alt'])*float(x['estonian_beta_alt'])>0
 assert float(x['estonian_ci_low'])<0<float(x['estonian_ci_high'])
assert hashlib.sha256((p/'direction_comparison.tsv').read_bytes()).hexdigest()==r['output_sha256']
out={'two_intervals_recalculated':True,'both_estonian_intervals_cross_zero':True,'output_sha256':r['output_sha256'],'scope':'Recalculates interval arithmetic only; does not establish independence, phenotype equivalence, replication or effective sample size.'}
(p/'direction_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
