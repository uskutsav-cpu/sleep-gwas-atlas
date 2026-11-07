#!/usr/bin/env python3
import csv, hashlib, json, math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'repair_ldz_comparison.tsv').open() as f: sums=list(csv.DictReader(f,delimiter='\t'))
with (p/'repair_ldz_details.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert len(sums)==12 and len(rows)==12*1547
assert {(x['reference'],x['trait']) for x in sums}=={(r,t) for r in ('eigenclip','ridge_0p001','1000g_eur') for t in ('oli','exc','insomnia','adhd')}
maxerr=0
for s in sums:
 assert not s['error']
 d=[x for x in rows if x['reference']==s['reference'] and x['trait']==s['trait']]
 assert len(d)==1547 and len({x['variant'] for x in d})==1547
 for x in d:
  err=abs((float(x['z'])-float(x['condmean']))/math.sqrt(float(x['condvar']))-float(x['z_std_diff']))
  maxerr=max(maxerr,err)
 assert math.isclose(max(abs(float(x['z_std_diff'])) for x in d),float(s['max_abs_z_std_diff']),abs_tol=1e-10)
 assert sum(abs(float(x['z_std_diff']))>3 for x in d)==int(s['count_gt3'])
 assert sum(abs(float(x['z_std_diff']))>5 for x in d)==int(s['count_gt5'])
assert maxerr<1e-9
r={'combinations_verified':len(sums),'rows_verified':len(rows),'maximum_residual_arithmetic_error':maxerr,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('repair_ldz_comparison.tsv','repair_ldz_details.tsv','repair_sensitivity_receipt.json','diagnose_repaired.R')},'limitation':'Checks summary and residual arithmetic, not SuSiE calibration, assumptions, or study LD compatibility.'}
(p/'repair_verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
