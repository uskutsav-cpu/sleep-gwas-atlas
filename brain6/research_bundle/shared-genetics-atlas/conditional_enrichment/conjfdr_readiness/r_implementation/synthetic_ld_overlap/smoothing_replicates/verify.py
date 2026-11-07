#!/usr/bin/env python3
import csv,hashlib,json,statistics
from pathlib import Path
p=Path(__file__).resolve().parent
def read(f):
 with f.open() as h:return list(csv.DictReader(h,delimiter='\t'))
rows=read(p/'runs.tsv');orig=read(p.parent/'runs.tsv')
assert len(rows)==20 and {(int(r['seed']),r['smooth']) for r in rows}=={(i,s) for i in range(69001,69011) for s in ('TRUE','FALSE')}
assert all(r['status']=='OK' and int(r['finite'])==44000 and int(r['nonfinite'])==0 for r in rows)
for r in rows:
 if r['smooth']=='TRUE':
  q=[x for x in orig if x['seed']==r['seed'] and float(x['rho'])==0];assert len(q)==1
  assert int(r['rejections'])==int(q[0]['rejections']) and abs(float(r['min_conj'])-float(q[0]['min_conj']))<1e-12
log=(p/'retry.log').read_text();assert 'ALL_COMPLETED' in log
res={}
for mode in ('TRUE','FALSE'):
 group=[r for r in rows if r['smooth']==mode]
 res[mode]={'runs':len(group),'runs_with_crossing':sum(int(r['rejections'])>0 for r in group),
  'total_crossings':sum(int(r['rejections']) for r in group),'min_run_minimum':min(float(r['min_conj']) for r in group),
  'median_run_minimum':statistics.median(float(r['min_conj']) for r in group)}
files=[p/x for x in ('PROTOCOL.md','run.R','run.log','retry.log','runs.tsv')]+[p.parent/'runs.tsv']
hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
receipt={'summary':res,'sha256':hashes,'initial_attempt':'Stopped on erroneous data.table row selection; no setting results from that attempt','note':'Paired post-hoc reanalysis; not independent calibration.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(res,indent=2))
