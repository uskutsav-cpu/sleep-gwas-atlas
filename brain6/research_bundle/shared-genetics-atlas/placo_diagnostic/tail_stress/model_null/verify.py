#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'model_probabilities.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with (p.parent/'tail_probabilities.tsv').open() as f:old=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==12
maximum_error=0
for pair in ('FI','FFS'):
 for level in ('nominal','screen'):
  q={r['subnull']:r for r in rows if r['pair']==pair and r['level']==level};assert set(q)=={'H00','H01','H02'}
  a=float(q['H00']['alpha']);reconstructed=float(q['H01']['uncapped_probability'])+float(q['H02']['uncapped_probability'])-float(q['H00']['uncapped_probability'])
  assert math.isclose(reconstructed,a,rel_tol=1e-6,abs_tol=1e-15)
  original=[r for r in old if r['pair']==pair and r['level']==level and float(r['mean_nonnull'])==0];assert len(original)==1
  assert math.isclose(float(q['H00']['capped_probability']),float(original[0]['probability']),rel_tol=1e-6,abs_tol=1e-15)
  for r in q.values():
   u=float(r['uncapped_probability']);v=float(r['capped_probability']);assert 0<=v<=u+1e-15 and u<=a+1e-15
   maximum_error=max(maximum_error,abs(u-float(r['alternate_uncapped'])),abs(v-float(r['alternate_capped'])))
receipt={'groups':4,'subnull_rows':len(rows),'reconstructed_official_tail':True,'h00_matches_prior_zero_mean':True,'all_stipulated_subnull_rates_le_alpha':True,'maximum_conditioning_difference':maximum_error,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('calculate.R','model_probabilities.tsv','runtime.txt')},'limits':'Arithmetic and model-conditional integrals only; no proof that the real GWAS obey the assumed null model or that genome-wide candidates are valid.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
