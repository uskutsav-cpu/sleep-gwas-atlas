import csv,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
r=read(p/'results.tsv'); old=read(p.parent/'finngen_adhd_raw_lookup/indexed_results.tsv'); by={x['rsid']:x for x in r};assert len(r)==len(by)==425
matches=[x for x in r if x['max_p']!='']; assert len(matches)==419
for x in matches:
 m=max(float(x['mvp_p']),float(x['hfrs_p']),float(x['adhd_p']));assert math.isclose(m,float(x['max_p']),abs_tol=1e-15) and math.isclose(min(1,425*m),float(x['max_p_times_425']),abs_tol=1e-15)
 assert set(x['grch38'].split(':')[-1].split('>')) not in ({'A','T'},{'C','G'})
checked=0
for x in old:
 if x['variant'] not in by or x['status']!='MATCHED': continue
 y=by[x['variant']]; assert y['adhd_status']=='matched' and math.isclose(float(y['adhd_p']),float(x['source_p']),abs_tol=1e-15)
 alt=y['grch38'].split('>')[-1]; expected=float(x['beta_aligned'])*(1 if x['effect_allele']==alt else -1);assert math.isclose(expected,float(y['adhd_beta_alt']),abs_tol=1e-12);checked+=1
assert checked==5
(p/'verification.json').write_text(json.dumps({'candidate_rows':len(r),'complete':len(matches),'old_lead_comparisons':checked,'crossings':sum(float(x['max_p_times_425'])<.05 for x in matches),'scope':'arithmetic and earlier selected lead comparison, not independent verification of full remote queries or cohort overlap'},indent=2)+'\n')
print('Verified 425 slots, 419 complete rows, five previous lead records, and all max-p values')
