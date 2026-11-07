import csv, hashlib, itertools, json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; source=p.parents[3]/'source_matrix.tsv'; prior=p.parent/'predictions.tsv'
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
r=read(source); old=read(prior); names=list(dict.fromkeys(x['sleep_trait'] for x in r)); feats=['adhd','mdd','scz','bipolar','alz','parkinson']; targets=['FI','FFS','HFRS']; lookup={(x['sleep_trait'],x['outcome']):x for x in r};assert len(lookup)==108 and len(names)==12
X=np.array([[float(lookup[n,k]['rg']) for k in feats] for n in names]);S=np.array([[float(lookup[n,k]['se']) for k in feats] for n in names]);Y=np.array([[float(lookup[n,k]['rg']) for k in targets] for n in names]); assert np.isfinite(S).all() and np.all(S>=0) and np.all(np.abs(Y)<1)
groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names}; oldby={(x['sleep_trait'],x['target']):x for x in old};assert len(oldby)==36
rows=[];maxerr=0
for i,n in enumerate(names):
 tr=np.array([j for j,m in enumerate(names) if groups[m]!=groups[n]]); mu=X[tr].mean(0); sd=X[tr].std(0);sd=np.where(sd>1e-12,sd,1); a=(X[tr]-mu)/sd
 for j,t in enumerate(targets):
  o=oldby[n,t];lam=float(o['lambda']); z=np.arctanh(Y[tr,j]);b=np.linalg.solve(a.T@a+lam*np.eye(6),a.T@(z-z.mean())); center=z.mean()+((X[i]-mu)/sd)@b;radius=np.sum(np.abs(b*S[i]/sd)); lo,hi=np.tanh([center-radius,center+radius]);original=float(np.tanh(center));assert abs(original-float(o['bounded_ridge']))<1e-10
  corners=[float(np.tanh(z.mean()+((X[i]+np.array(v)*S[i]-mu)/sd)@b)) for v in itertools.product([-1,1],repeat=6)]; maxerr=max(maxerr,abs(min(corners)-lo),abs(max(corners)-hi));
  rows.append({'sleep_trait':n,'target':t,'held_group':groups[n],'observed':Y[i,j],'prediction':original,'min_1se_box':lo,'max_1se_box':hi,'box_width':hi-lo,'crosses_zero':bool(lo<=0<=hi),'input_outside_correlation_range':bool(np.any(np.abs(X[i])+S[i]>1)),'lambda':lam})
assert len(rows)==36 and maxerr<1e-12
with (p/'ranges.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(z):hashlib.sha256(z.read_bytes()).hexdigest() for z in (source,prior)},'rows':len(rows),'corners_per_row':64,'max_analytic_corner_difference':maxerr},indent=2)+'\n')
print(json.dumps({'widest':sorted(rows,key=lambda x:-x['box_width'])[:4],'cross_zero':sum(x['crosses_zero'] for x in rows),'invalid_input_bounds':sum(x['input_outside_correlation_range'] for x in rows),'corner_max_error':maxerr},indent=2))
