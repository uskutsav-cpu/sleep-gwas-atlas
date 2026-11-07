#!/usr/bin/env python3
import csv,hashlib,json,math,statistics
from pathlib import Path
from statistics import NormalDist
p=Path(__file__).resolve().parent;source=p.parent/'matched.tsv';base=p.parent/'posteriors.tsv'
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
d=read(source);baseline=read(base);assert len(d)==2441
normal=NormalDist()
def vals(row,trait):
 b,se,pv=map(float,(row[trait+'_beta_alt'] if trait=='fg' else row['qtl_beta'], row[trait+'_se'],row[trait+'_p']))
 assert math.isfinite(b) and se>0 and 0<pv<=1
 z=b/se;zm=0 if pv==1 else -normal.inv_cdf(pv/2)
 return b,se,pv,z,math.copysign(zm,b)
def logbf(z,se,sd):
 v=se*se;r=sd*sd/(sd*sd+v);return .5*(math.log1p(-r)+r*z*z)
def lse(v):
 m=max(v);return m+math.log(sum(math.exp(x-m) for x in v))
def posterior(x,y,p12):
 a=lse(x);b=lse(y);c=lse([u+v for u,v in zip(x,y)]);ab=a+b
 h3=ab+math.log1p(-math.exp(c-ab))+2*math.log(1e-4)
 logs=[0,math.log(1e-4)+a,math.log(1e-4)+b,h3,math.log(p12)+c];den=lse(logs)
 return [math.exp(z-den) for z in logs]
q=[vals(t,'qtl') for t in d];g=[vals(t,'fg') for t in d]
out=[]
def compute(label,qz,gz,sd):
 x=[logbf(v[qz],v[1],.15) for v in q];y=[logbf(v[gz],v[1],sd) for v in g]
 for p12 in (1e-5,1e-6):
  h=posterior(x,y,p12);out.append({'analysis':label,'prior_sd_disease':sd,'p12':p12,**{f'H{i}':v for i,v in enumerate(h)}})
for sd in (.1,.2,.3):compute('beta_over_se',3,3,sd)
for qz,gz,name in [(4,3,'qtl_reported_p_hybrid'),(3,4,'finngen_reported_p_hybrid'),(4,4,'both_reported_p_hybrid')]:compute(name,qz,gz,.2)
with (p/'posteriors.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
metrics={}
for name,vs in [('qtl',q),('finngen',g)]:
 dif=sorted(abs(v[3]-v[4]) for v in vs);p_dif=sorted(abs(math.erfc(abs(v[3])/math.sqrt(2))-v[2]) for v in vs)
 metrics[name]={'median_abs_z_difference':statistics.median(dif),'p95_abs_z_difference':dif[int(.95*(len(dif)-1))],'max_abs_z_difference':max(dif),'max_abs_p_difference':max(p_dif)}
check=next(r for r in baseline if r['trait']=='finngen' and float(r['p12'])==1e-5);new=next(r for r in out if r['analysis']=='beta_over_se' and r['prior_sd_disease']==.2 and r['p12']==1e-5)
assert max(abs(float(check[f'H{i}'])-new[f'H{i}']) for i in range(5))<1e-10
(p/'receipt.json').write_text(json.dumps({'n':len(d),'statistic_discrepancy':metrics,'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline_verified':True,'limitation':'p-derived Z with SE variance is a numerical hybrid, not a validated likelihood'},indent=2)+'\n')
print(json.dumps(metrics,indent=2));print((p/'posteriors.tsv').read_text())
