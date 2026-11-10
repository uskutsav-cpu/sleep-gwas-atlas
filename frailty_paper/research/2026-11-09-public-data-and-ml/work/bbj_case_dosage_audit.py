import csv,json,hashlib
from pathlib import Path
p=Path('outputs/stage6_bbj_insomnia_oriented_candidates.csv');q=Path('outputs/stage6_bbj_insomnia_matched14_z.json');prior=json.loads(q.read_text());ids=set(prior['markers']);ncase=389;nctrl=177864
with p.open(newline='') as f: rows=[r for r in csv.DictReader(f) if r['score']=='shared22' and r['status']=='exact position and pair']
assert len(rows)==18 and len(ids)==14 and ids.issubset({r['variant'] for r in rows})
out=[]
for r in rows:
 af=float(r['allele2_dosage_count'])/(2*(ncase+nctrl));p_case=float(r['allele2_frequency_in_cases']);p_ctrl=float(r['allele2_frequency_in_controls']); assert all(0<=x<=1 for x in (af,p_case,p_ctrl))
 if af>0.5:p_case,p_ctrl=1-p_case,1-p_ctrl
 out.append({'variant':r['variant'],'in_matched14':r['variant'] in ids,'estimated_case_minor_dosage_sum':2*ncase*p_case,'estimated_control_minor_dosage_sum':2*nctrl*p_ctrl})
def summary(items):
 v=[x['estimated_case_minor_dosage_sum'] for x in items];return {'n_variants':len(v),'minimum_case_minor_dosage_sum':min(v),'variants_below_5':sum(x<5 for x in v),'variants_below_20':sum(x<20 for x in v),'variants_below_50':sum(x<50 for x in v)}
res={'ncase_assumed':ncase,'ncontrols_assumed':nctrl,'matched18':summary(out),'matched14':summary([x for x in out if x['in_matched14']]),'input_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [p,q]},'limits':'Approximations from source case/control frequencies and advertised total N, assuming all case and control samples contribute at every variant. Minor allele defined using overall allele2 dosage frequency. Imputed dosage sums are not counts of people or hard genotype counts; sample missingness and uncertainty unknown. Postselected and not score validation.'}
with Path('outputs/stage6_bbj_case_dosage_components.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
Path('outputs/stage6_bbj_case_dosage_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
