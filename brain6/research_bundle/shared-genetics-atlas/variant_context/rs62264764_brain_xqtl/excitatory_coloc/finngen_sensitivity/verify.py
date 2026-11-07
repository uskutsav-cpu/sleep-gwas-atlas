#!/usr/bin/env python3
import csv,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def read(x):
 with x.open() as f:return list(csv.DictReader(f,delimiter='\t'))
d=read(p/'matched.tsv');outputs=read(p/'posteriors.tsv')
def lse(v):
 m=max(v);return m+math.log(sum(math.exp(x-m) for x in v))
def logbf(beta,se,sd):
 v=se*se;r=sd*sd/(sd*sd+v);return .5*(math.log1p(-r)+r*(beta/se)**2)
x=[logbf(float(t['qtl_beta']),float(t['qtl_se']),.15) for t in d];a=lse(x);maxdiff=0
for r in outputs:
 bkey,skey=('gwas_beta_alt','gwas_se') if r['trait']=='brain6' else ('fg_beta_alt','fg_se')
 y=[logbf(float(t[bkey]),float(t[skey]),.2) for t in d];b=lse(y);c=lse([u+v for u,v in zip(x,y)]);ab=a+b
 logs=[0,math.log(1e-4)+a,math.log(1e-4)+b,2*math.log(1e-4)+ab+math.log1p(-math.exp(c-ab)),math.log(float(r['p12']))+c]
 den=lse(logs);post=[math.exp(z-den) for z in logs]
 maxdiff=max(maxdiff,max(abs(post[k]-float(r['H'+str(k)])) for k in range(5)))
assert len(outputs)==4 and maxdiff<1e-10
v=[x for x in d if x['variant']=='chr3:117920728:G:A']
(p/'verification.json').write_text(json.dumps({'fits':len(outputs),'n':len(d),'maximum_absolute_posterior_difference':maxdiff,'target_rows':len(v),'target':v[0] if v else None,'limitations':'Arithmetic only; no source, LD, model or causal calibration.'},indent=2)+'\n')
print('fits',len(outputs),'n',len(d),'max diff',maxdiff,'target',len(v))
