import itertools,json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;y=np.loadtxt(p/'joint_fit/moments.txt');R=np.eye(3)
for (i,j),v in zip(((0,1),(0,2),(1,2)),y[:3]):R[i,j]=R[j,i]=v
assert np.linalg.eigvalsh(R).min()>0
def value(w,r):return float(w@r/np.sqrt(w@R@w))
out={}
for source,r in [('registry',y[3:6]),('questionnaire',y[6:9])]:
 candidates=[]
 for n in (1,2,3):
  for support in itertools.combinations(range(3),n):
   indices=list(support);direction=np.linalg.solve(R[np.ix_(indices,indices)],r[indices])
   for sign in (-1,1):
    a=sign*direction
    if (a>=0).all() and a.sum()>0:
     w=np.zeros(3);w[indices]=a/a.sum();candidates.append({'weights_walking_grip_activity':w.tolist(),'correlation':value(w,r)})
 candidates.sort(key=lambda x:x['correlation']);assert candidates
 lo,hi=candidates[0],candidates[-1]
 grid=[]
 for i in range(101):
  for j in range(101-i):
   w=np.array([i,j,100-i-j])/100;grid.append(value(w,r))
 assert min(grid)>=lo['correlation']-1e-10 and max(grid)<=hi['correlation']+1e-10
 assert min(grid)-lo['correlation']<.001 and hi['correlation']-max(grid)<.001
 path=[]
 for t in (0,.25,1/3,.5,.75,1):
  w=np.array([(1-t)/2,(1-t)/2,t]);path.append({'activity_weight':t,'correlation':value(w,r)})
 out[source]={'minimum':lo,'maximum':hi,'fixed_activity_weight_path':path,'grid_minimum':min(grid),'grid_maximum':max(grid),'qualification':'Conditional point-estimate range across outcome-optimized nonnegative weights; not a confidence interval or validated weight choice.'}
(p/'joint_fit/weight_sensitivity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
