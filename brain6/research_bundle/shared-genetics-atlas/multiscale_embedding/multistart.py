#!/usr/bin/env python3
import csv,collections,json,math,pathlib,numpy as np
p=pathlib.Path('outputs/shared-genetics-atlas/multiscale_embedding')
def read(n):
 with (p/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
s=read('source_matrix.tsv');a=read('masked_predictions.tsv');sleeps=list(dict.fromkeys(x['sleep_trait'] for x in s));outs=list(dict.fromkeys(x['outcome'] for x in s));assert len(sleeps)==12 and len(outs)==8 and len(s)==96 and len(a)==960
M=np.empty((12,8));seen=set()
for x in s:
 i=sleeps.index(x['sleep_trait']);j=outs.index(x['outcome']);assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(x['rg'])
groups=collections.defaultdict(list)
for x in a:groups[(int(x['repeat']),int(x['fold']))].append(x)
assert len(groups)==50
def fit(mask,seed):
 rng=np.random.RandomState(seed);U=rng.normal(scale=.15,size=(12,2));V=rng.normal(scale=.15,size=(8,2));I=np.eye(2)*.02
 for _ in range(100):
  for i in range(12):
   z=np.flatnonzero(mask[i]);A=V[z];U[i]=np.linalg.solve(A.T@A+I,A.T@M[i,z])
  for j in range(8):
   z=np.flatnonzero(mask[:,j]);A=U[z];V[j]=np.linalg.solve(A.T@A+I,A.T@M[z,j])
 pred=U@V.T;objective=float(np.sum((M[mask]-pred[mask])**2)+.02*(np.sum(U*U)+np.sum(V*V)));return objective,pred
results=[];trains=[]
for (rep,fold),rows in sorted(groups.items()):
 mask=np.ones((12,8),dtype=bool)
 for x in rows:mask[sleeps.index(x['sleep_trait']),outs.index(x['outcome'])]=False
 fits=[fit(mask,7000+rep*50+fold*10+i) for i in range(5)];best=min(range(5),key=lambda i:fits[i][0]);trains.extend({'repeat':rep,'fold':fold,'initialization':i,'seed':7000+rep*50+fold*10+i,'training_objective':fits[i][0],'selected':i==best} for i in range(5));mat=fits[best][1]
 for x in rows:
  i=sleeps.index(x['sleep_trait']);j=outs.index(x['outcome']);assert float(x['observed'])==M[i,j];results.append({'repeat':rep,'fold':fold,'sleep_trait':x['sleep_trait'],'outcome':x['outcome'],'observed':M[i,j],'column_mean':x['column_mean'],'original_rank2':x['rank2'],'additive':x['additive'],'multistart_rank2':mat[i,j]})
for name,data in [('multistart_predictions.tsv',results),('multistart_training_objectives.tsv',trains)]:
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
summary={'n_predictions':len(results),'n_fits':len(trains),'selection':'minimum observed-training SSE plus 0.02*(sum(U^2)+sum(V^2))','rmse':{},'out_of_range':{}}
for k in ('column_mean','original_rank2','additive','multistart_rank2'):
 summary['rmse'][k]=math.sqrt(sum((float(r[k])-r['observed'])**2 for r in results)/len(results));summary['out_of_range'][k]=sum(abs(float(r[k]))>1 for r in results)
(p/'multistart_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
