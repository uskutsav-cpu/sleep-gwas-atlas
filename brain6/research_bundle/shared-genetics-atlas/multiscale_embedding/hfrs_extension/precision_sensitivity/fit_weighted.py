#!/usr/bin/env python3
import csv,collections,hashlib,json,math,pathlib,numpy as np
p=pathlib.Path('outputs/shared-genetics-atlas/multiscale_embedding/hfrs_extension');out=p/'precision_sensitivity'
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
s=read(p/'source_matrix.tsv');prior=read(p/'masked_predictions.tsv');sleeps=list(dict.fromkeys(x['sleep_trait'] for x in s));traits=list(dict.fromkeys(x['outcome'] for x in s));assert len(s)==108 and len(sleeps)==12 and len(traits)==9 and len(prior)==1080
M=np.empty((12,9));SE=np.empty_like(M);seen=set()
for x in s:
 i=sleeps.index(x['sleep_trait']);j=traits.index(x['outcome']);assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(x['rg']);SE[i,j]=float(x['se'])
assert len(seen)==108 and np.all(np.isfinite(SE)) and np.all(SE>0)
groups=collections.defaultdict(list)
for r in prior:groups[(int(r['repeat']),int(r['fold']))].append(r)
assert len(groups)==50
design=np.zeros((108,20));design[:,0]=1
for i in range(12):
 for j in range(9):
  if i:design[9*i+j,i]=1
  if j:design[9*i+j,11+j]=1
def fit(mask,w,seed):
 rng=np.random.RandomState(seed);U=rng.normal(scale=.15,size=(12,2));V=rng.normal(scale=.15,size=(9,2));I=np.eye(2)*.02
 for _ in range(100):
  for i in range(12):
   z=np.flatnonzero(mask[i]);A=V[z];wz=w[i,z];U[i]=np.linalg.solve(A.T@(A*wz[:,None])+I,A.T@(M[i,z]*wz))
  for j in range(9):
   z=np.flatnonzero(mask[:,j]);A=U[z];wz=w[z,j];V[j]=np.linalg.solve(A.T@(A*wz[:,None])+I,A.T@(M[z,j]*wz))
 pred=U@V.T;obj=float(np.sum(w[mask]*(M[mask]-pred[mask])**2)+.02*(np.sum(U*U)+np.sum(V*V)));return obj,pred
rows=[];trials=[];foldinfo=[]
for (rep,fold),grp in sorted(groups.items()):
 mask=np.ones((12,9),dtype=bool)
 for x in grp:
  i=sleeps.index(x['sleep_trait']);j=traits.index(x['outcome']);assert mask[i,j] and float(x['observed'])==M[i,j];mask[i,j]=False
 med=float(np.median(SE[mask]));w=np.clip((med/SE)**2,.25,4);foldinfo.append({'repeat':rep,'fold':fold,'training_median_se':med,'n_train':int(mask.sum()),'n_clipped_low':int(np.sum(w[mask]==.25)),'n_clipped_high':int(np.sum(w[mask]==4))})
 column=np.array([np.average(M[mask[:,j],j],weights=w[mask[:,j],j]) for j in range(9)])
 wm=w.ravel()[mask.ravel()];X=design[mask.ravel()];coef=np.linalg.lstsq(X*np.sqrt(wm)[:,None],M.ravel()[mask.ravel()]*np.sqrt(wm),rcond=None)[0];add=(design@coef).reshape(12,9)
 test=[fit(mask,w,9000+rep*50+fold*10+i) for i in range(5)];best=min(range(5),key=lambda i:test[i][0]);pred=test[best][1]
 for j,(obj,_) in enumerate(test):trials.append({'repeat':rep,'fold':fold,'initialization':j,'seed':9000+rep*50+fold*10+j,'objective':obj,'selected':j==best})
 for x in grp:
  i=sleeps.index(x['sleep_trait']);j=traits.index(x['outcome']);rows.append({'repeat':rep,'fold':fold,'sleep_trait':x['sleep_trait'],'outcome':x['outcome'],'observed':M[i,j],'SE':SE[i,j],'unweighted_mean':float(x['column_mean']),'unweighted_additive':float(x['additive']),'unweighted_rank2':float(x['rank2_multistart']),'weighted_mean':column[j],'weighted_additive':add[i,j],'weighted_rank2':pred[i,j]})
def write(name,r):
 with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,list(r[0]),delimiter='\t');w.writeheader();w.writerows(r)
write('predictions.tsv',rows);write('training_objectives.tsv',trials);write('fold_weights.tsv',foldinfo)
summary=[]
for group in ['ALL']+traits:
 vals=[x for x in rows if group=='ALL' or x['outcome']==group]
 for model in ['unweighted_mean','weighted_mean','unweighted_additive','weighted_additive','unweighted_rank2','weighted_rank2']:
  err=np.array([x[model]-x['observed'] for x in vals]);se=np.array([x['SE'] for x in vals]);summary.append({'group':group,'model':model,'n':len(vals),'rmse':float(np.sqrt(np.mean(err**2))),'se_scaled_rmse':float(np.sqrt(np.mean((err/se)**2))),'out_of_range':sum(abs(x[model])>1 for x in vals)})
write('summary.tsv',summary)
meta={'source_matrix_sha256':hashlib.sha256((p/'source_matrix.tsv').read_bytes()).hexdigest(),'source_masks_sha256':hashlib.sha256((p/'masked_predictions.tsv').read_bytes()).hexdigest(),'n_predictions':len(rows),'folds':len(groups),'fits':len(trials),'weights':'clip((median_training_SE/SE)^2,0.25,4)','selection':'minimum weighted training SSE + 0.02 ridge','caveat':'SE-scaled errors are descriptive; covariance and external validation missing'}
(out/'receipt.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2));print([x for x in summary if x['group'] in ('ALL','FI','FFS','HFRS')])
