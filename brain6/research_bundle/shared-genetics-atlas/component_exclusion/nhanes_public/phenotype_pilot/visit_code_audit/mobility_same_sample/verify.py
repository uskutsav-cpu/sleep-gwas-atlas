#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
rows=read(p/'correlations.tsv');pairs=read(p/'paired_mobility.tsv');counts=read(p/'counts.tsv');prior=read(p.parent/'corrected/correlations.tsv')
assert len(rows)==8 and len(pairs)==4 and all(int(x['eligible_n'])==2986 and int(x['df'])==32 for x in counts)
for x in pairs:
 a=next(r for r in rows if r['mode']=='mobility_zero' and r['x']==x['x'] and r['y']==x['y']);b=next(r for r in rows if r['mode']=='mobility_one' and r['x']==x['x'] and r['y']==x['y'])
 prev=next(r for r in prior if r['mode']=='mobility_proxy' and r['x']==x['x'] and r['y']==x['y'])
 assert math.isclose(float(b['estimate'])-float(a['estimate']),float(x['difference_one_minus_zero']),abs_tol=1e-12)
 for key in ('estimate','se','low','high'):assert math.isclose(float(b[key]),float(prev[key]),abs_tol=1e-12)
 assert float(x['low'])<float(x['difference_one_minus_zero'])<float(x['high'])
sleep=read(p/'paired_contrasts.tsv');change=read(p/'paired_sleep_contrast_change.tsv');assert len(sleep)==2 and len(change)==1
z=next(x for x in sleep if x['mode']=='mobility_zero');o=next(x for x in sleep if x['mode']=='mobility_one')
assert math.isclose(float(o['estimate'])-float(z['estimate']),float(change[0]['difference_one_minus_zero']),abs_tol=1e-12)
assert int(change[0]['eligible_n'])==2986 and float(change[0]['low'])<float(change[0]['difference_one_minus_zero'])<float(change[0]['high'])
r={'identical_eligible_counts':2986,'comparisons':len(pairs),'paired_differences_reproduced':True,'sleep_contrast_change_reproduced':True,'mobility_one_matches_corrected_prior':True,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('correlations.tsv','paired_mobility.tsv','paired_sleep_contrast_change.tsv','counts.tsv','run_generated.R')},'limitations':'Arithmetic checks and prior-output reproduction; do not validate missing-data assumptions, original instruments or causal interpretations.'}
(p/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
