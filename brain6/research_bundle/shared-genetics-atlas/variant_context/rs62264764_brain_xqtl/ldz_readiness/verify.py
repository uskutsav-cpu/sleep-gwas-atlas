import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'ldz_details.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with (p/'ldz_summary.tsv').open() as f:sums=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==3096 and len(sums)==2
maxerr=0
for s in sums:
 d=[x for x in rows if x['context']==s['context']];assert len(d)==1548 and len({x['variant'] for x in d})==1548
 for r in d:
  expected=(float(r['z'])-float(r['condmean']))/math.sqrt(float(r['condvar']));maxerr=max(maxerr,abs(expected-float(r['z_std_diff'])))
 assert math.isclose(max(abs(float(x['z_std_diff'])) for x in d),float(s['max_abs_z_std_diff']),abs_tol=1e-12)
 assert sum(abs(float(x['z_std_diff']))>3 for x in d)==int(s['count_abs_z_std_diff_gt3'])
 assert sum(abs(float(x['z_std_diff']))>5 for x in d)==int(s['count_abs_z_std_diff_gt5'])
assert maxerr<1e-9
r={'counts_reproduced':True,'n_rows':len(rows),'max_standardized_residual_error':maxerr,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('ldz_details.tsv','ldz_summary.tsv','diagnose.R','matching_receipt.json')},'caveat':'Verifies summary arithmetic, not calibration or study-matched LD.'}
(p/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
