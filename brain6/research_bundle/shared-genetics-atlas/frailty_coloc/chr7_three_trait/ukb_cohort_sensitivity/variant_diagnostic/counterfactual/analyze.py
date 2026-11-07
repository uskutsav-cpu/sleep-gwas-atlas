import csv, math, json, hashlib
from pathlib import Path
p=Path(__file__).resolve().parent;f=Path('/Volumes/Brain6Work/chr7-three-trait-ukb-sensitivity.tsv');x=list(csv.DictReader(f.open(),delimiter='\t'));assert len(x)==1055 and len({r['rsid'] for r in x})==len(x)
def ls(v):
 m=max(v);return m+math.log(math.fsum(math.exp(a-m) for a in v))
def abf(beta,se,sd):
 b=float(beta);s=float(se);r=sd**2/(sd**2+s*s);return .5*(math.log1p(-r)+r*(b/s)**2)
result=[]
for trait in ['adhd','insomnia']:
 q=[abf(a[trait+'_beta'],a[trait+'_se'],.2) for a in x];sq=ls(q)
 for bsrc,ssrc in [('fi','fi'),('ukb_fi','fi'),('fi','ukb_fi'),('ukb_fi','ukb_fi')]:
  z=[abf(a[bsrc+'_beta'],a[ssrc+'_se'],.15) for a in x];sz=ls(z);s=ls([i+j for i,j in zip(z,q)]);d=sz+sq+math.log1p(-math.exp(s-sz-sq))
  for prior in [1e-5,1e-6]:
   weights=[0,math.log(1e-4)+sz,math.log(1e-4)+sq,math.log(1e-8)+d,math.log(prior)+s];total=ls(weights);h=[math.exp(v-total) for v in weights]
   result.append(dict(trait=trait,beta_source=bsrc,se_source=ssrc,p12=prior,H3=h[3],H4=h[4]))
old=list(csv.DictReader((p.parents[1]/'pairwise.tsv').open(),delimiter='\t'))
for v in result:
 if v['beta_source']!=v['se_source']:continue
 source='meta' if v['beta_source']=='fi' else 'UKB'
 row=next(z for z in old if z['source']==source and z['pair'].lower()==v['trait'] and float(z['p12'])==v['p12'])
 assert abs(float(row['H4'])-v['H4'])<1e-10 and abs(float(row['H3'])-v['H3'])<1e-10
with (p/'counterfactual.tsv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=result[0],delimiter='\t');w.writeheader();w.writerows(result)
variants=[]
for rid in ['rs9969232','rs936145','rs2396766']:
 a=next(z for z in x if z['rsid']==rid);variants.append({'rsid':rid,'effect_allele':a['effect_allele'],'meta_fi_z':float(a['fi_beta'])/float(a['fi_se']),'ukb_fi_z':float(a['ukb_fi_beta'])/float(a['ukb_fi_se'])})
out={'input_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'n':len(x),'natural_model_h3_h4_reproduced':True,'selected_fi_z':variants,'warning':'Hybrids are mathematical perturbations; beta and SE from different fits have no sampling interpretation.'};(p/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('\n'.join(str(v) for v in result if v['trait']=='adhd' and v['p12']==1e-5))
