#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
sums=read(p/'adhd_min_n_summary.tsv');rows=read(p/'adhd_min_n_details.tsv')
assert len(sums)==3 and len(rows)==3*1547
maxerr=0
for s in sums:
 assert not s['error'] and int(s['n'])==64699
 d=[x for x in rows if x['reference']==s['reference']];assert len(d)==1547 and len({x['variant'] for x in d})==1547
 for x in d:maxerr=max(maxerr,abs((float(x['z'])-float(x['condmean']))/math.sqrt(float(x['condvar']))-float(x['z_std_diff'])))
 assert math.isclose(max(abs(float(x['z_std_diff'])) for x in d),float(s['max_abs_z_std_diff']),abs_tol=1e-10)
 assert sum(abs(float(x['z_std_diff']))>3 for x in d)==int(s['count_gt3'])
 assert sum(abs(float(x['z_std_diff']))>5 for x in d)==int(s['count_gt5'])
assert maxerr<1e-9
receipt={'references_verified':3,'rows_verified':len(rows),'max_residual_arithmetic_error':maxerr,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('adhd_min_n_summary.tsv','adhd_min_n_details.tsv','adhd_source_audit.json','adhd_coverage_strata.json')},'limitations':'Verifies arithmetic and file integrity, not validity of sample sizes, genotype correlations, statistical calibration or causal hypotheses.'}
(p/'adhd_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
