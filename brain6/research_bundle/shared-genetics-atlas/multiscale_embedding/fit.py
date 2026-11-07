#!/usr/bin/env python3
import csv,hashlib,json,math,pathlib,platform,numpy as np
root=pathlib.Path('outputs/shared-genetics-atlas');out=root/'multiscale_embedding';b=root/'trait_embedding/source_brain6_72.tsv';f=root/'evidence/ALL12_FI_FFS_CONTRAST_RECOMPUTED.tsv'
with b.open() as x:brain=list(csv.DictReader(x,delimiter='\t'))
with f.open() as x:frail=list(csv.DictReader(x,delimiter='\t'))
sleeps=list(dict.fromkeys(r['sleep_trait'] for r in brain));disorders=list(dict.fromkeys(r['brain_disorder'] for r in brain));assert len(sleeps)==12 and len(disorders)==6 and len(brain)==72 and len(frail)==12;assert len({r['sleep_trait'] for r in frail})==12 and set(sleeps)=={r['sleep_trait'] for r in frail};traits=disorders+['FI','FFS'];M=np.empty((12,8));S=np.empty_like(M);seen=set();sources=[]
for r in brain:
 i=sleeps.index(r['sleep_trait']);j=traits.index(r['brain_disorder']);assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(r['rg']);S[i,j]=float(r['se']);sources.append((sleeps[i],traits[j],M[i,j],S[i,j],'Brain6 frozen 72; '+r['source_GWAS_ID_sleep']))
for r in frail:
 assert r['status']=='ESTIMATED_SENSITIVITY' and r['h2_gate_pass']=='True'
 i=sleeps.index(r['sleep_trait'])
 for name in ('FI','FFS'):
  j=traits.index(name);assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(r[name+'_rg']);S[i,j]=float(r[name+'_SE']);sources.append((sleeps[i],name,M[i,j],S[i,j],'paired FI-FFS archive; source by sleep label, version identity not independently established'))
assert len(seen)==96 and np.all(np.isfinite(M)) and np.all(np.isfinite(S)) and np.all(S>0)
u,d,v=np.linalg.svd(M,full_matrices=False)
for k in range(len(d)):
 sign=1 if v[k,np.argmax(abs(v[k]))]>=0 else -1;u[:,k]*=sign;v[k]*=sign
sleep=u[:,:2]*np.sqrt(d[:2]);trait=v[:2,:].T*np.sqrt(d[:2]);rec=sleep@trait.T
with (out/'source_matrix.tsv').open('w') as x:
 w=csv.writer(x,delimiter='\t');w.writerow(['sleep_trait','outcome','rg','se','source']);w.writerows(sources)
with (out/'coordinates.tsv').open('w') as x:
 w=csv.writer(x,delimiter='\t');w.writerow(['type','trait','dimension_1','dimension_2']);
 for label,names,coords in [('sleep',sleeps,sleep),('outcome',traits,trait)]:
  for name,c in zip(names,coords):w.writerow([label,name,*c])
with (out/'reconstruction.tsv').open('w') as x:
 w=csv.writer(x,delimiter='\t');w.writerow(['sleep_trait','outcome','rg','se','rank2_reconstruction','residual']);
 for i,s in enumerate(sleeps):
  for j,t in enumerate(traits):w.writerow([s,t,M[i,j],S[i,j],rec[i,j],M[i,j]-rec[i,j]])
def als(obs,mask,seed,rank=2,lam=.02,steps=100):
 rng=np.random.RandomState(seed);U=rng.normal(scale=.15,size=(12,rank));V=rng.normal(scale=.15,size=(8,rank));I=np.eye(rank)*lam
 for _ in range(steps):
  for i in range(12):
   z=np.flatnonzero(mask[i]);A=V[z];U[i]=np.linalg.solve(A.T@A+I,A.T@obs[i,z])
  for j in range(8):
   z=np.flatnonzero(mask[:,j]);A=U[z];V[j]=np.linalg.solve(A.T@A+I,A.T@obs[z,j])
 return U@V.T
rng=np.random.RandomState(20261027);cv=[];design=np.zeros((96,19));design[:,0]=1
for i in range(12):
 for j in range(8):
  if i:design[8*i+j,i]=1
  if j:design[8*i+j,11+j]=1
for rep in range(10):
 while True:
  perm=rng.permutation(96);fold=np.empty(96,dtype=int);fold[perm]=np.arange(96)%5;fold=fold.reshape(12,8)
  if all(np.all((fold!=k).sum(axis=0)>=3) and np.all((fold!=k).sum(axis=1)>=3) for k in range(5)):break
 for k in range(5):
  train=fold!=k;mean=np.array([np.mean(M[train[:,j],j]) for j in range(8)]);pred=als(M,train,3000+rep*5+k);coeff=np.linalg.lstsq(design[train.ravel()],M.ravel()[train.ravel()],rcond=None)[0];add=(design@coeff).reshape(12,8)
  for i,j in zip(*np.where(~train)):cv.append([rep,k,sleeps[i],traits[j],M[i,j],mean[j],pred[i,j],add[i,j]])
with (out/'masked_predictions.tsv').open('w') as x:
 w=csv.writer(x,delimiter='\t');w.writerow(['repeat','fold','sleep_trait','outcome','observed','column_mean','rank2','additive']);w.writerows(cv)
cv=np.array([[float(z) for z in row[4:]] for row in cv]);rmse=lambda a,b:float(np.sqrt(np.mean((a-b)**2)));meta={'sources_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (b,f)},'n':96,'rank2_in_sample_squared_magnitude_fraction':float(sum(d[:2]**2)/sum(d**2)),'in_sample_rmse':rmse(M,rec),'masked_cell_rmse':{name:rmse(cv[:,0],cv[:,j]) for j,name in [(1,'column_mean'),(2,'rank2_ALS'),(3,'additive')]},'n_predictions':len(cv),'folds':5,'repeats':10,'seed':20261027,'rank2_lambda':.02,'ALS_steps':100,'python':platform.python_version(),'numpy':np.__version__}
(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
