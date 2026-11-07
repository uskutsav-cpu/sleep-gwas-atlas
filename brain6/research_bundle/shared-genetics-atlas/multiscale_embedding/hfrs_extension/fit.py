#!/usr/bin/env python3
import csv,hashlib,json,math,pathlib,platform,numpy as np
root=pathlib.Path('outputs/shared-genetics-atlas');out=root/'multiscale_embedding/hfrs_extension'
inputs={'matrix':root/'multiscale_embedding/source_matrix.tsv','hfrs':root/'evidence/NEW_VALIDATION_FAMILY_RESULTS.tsv','crosswalk':root/'multiscale_embedding/sleep_source_crosswalk.tsv'}
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
a,b,c=(read(inputs[n]) for n in ('matrix','hfrs','crosswalk'));h=[x for x in b if x['trait2']=='HFRS_GCST90707291'];assert len(a)==96 and len(h)==12 and len(c)==12
sleeps=list(dict.fromkeys(x['sleep_trait'] for x in a));traits=list(dict.fromkeys(x['outcome'] for x in a))+['HFRS'];assert len(sleeps)==12 and len(traits)==9 and len(set(sleeps))==12
cross={x['sleep_trait']:x for x in c};hm={x['trait1']:x for x in h};assert set(cross)==set(hm)==set(sleeps)
M=np.empty((12,9));records=[];seen=set()
for x in a:
 i=sleeps.index(x['sleep_trait']);j=traits.index(x['outcome']);assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(x['rg']);records.append(dict(x,metadata_match='existing_12x8'))
for name in sleeps:
 x=hm[name];z=cross[name];assert x['status']=='ADMITTED_EXPLORATORY_GLOBAL_RESULT' and x['source1']==z['core_dataset_version'] and int(float(x['n1']))==int(float(z['core_sleep_n'])) and float(x['SE'])>0
 i=sleeps.index(name);j=8;assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(x['rg_or_contrast']);records.append({'sleep_trait':name,'outcome':'HFRS','rg':M[i,j],'se':float(x['SE']),'source':x['source1']+'; '+x['source2'],'metadata_match':'sleep_dataset_and_N_match;paired_bytes_not_established'})
assert len(seen)==108 and np.all(np.isfinite(M));u,d,v=np.linalg.svd(M,full_matrices=False);svd=(u[:,:2]*d[:2])@v[:2,:]
with (out/'source_matrix.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,['sleep_trait','outcome','rg','se','source','metadata_match'],delimiter='\t');w.writeheader();w.writerows(records)
for k in range(2):
 if v[k,np.argmax(abs(v[k]))]<0:u[:,k]*=-1;v[k]*=-1
with (out/'coordinates.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['type','trait','dimension_1','dimension_2']);
 for kind,names,coords in [('sleep',sleeps,u[:,:2]*np.sqrt(d[:2])),('outcome',traits,v[:2,:].T*np.sqrt(d[:2]))]:
  for name,c in zip(names,coords):w.writerow([kind,name,*c])
with (out/'reconstruction.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sleep_trait','outcome','observed','rank2_reconstruction','residual']);
 for i,name in enumerate(sleeps):
  for j,trait in enumerate(traits):w.writerow([name,trait,M[i,j],svd[i,j],M[i,j]-svd[i,j]])
def fit(mask,seed):
 rng=np.random.RandomState(seed);U=rng.normal(scale=.15,size=(12,2));V=rng.normal(scale=.15,size=(9,2));I=.02*np.eye(2)
 for _ in range(100):
  for i in range(12):
   z=np.flatnonzero(mask[i]);A=V[z];U[i]=np.linalg.solve(A.T@A+I,A.T@M[i,z])
  for j in range(9):
   z=np.flatnonzero(mask[:,j]);A=U[z];V[j]=np.linalg.solve(A.T@A+I,A.T@M[z,j])
 pred=U@V.T;objective=float(np.sum((M[mask]-pred[mask])**2)+.02*(np.sum(U*U)+np.sum(V*V)));return objective,pred
rng=np.random.RandomState(20261029);records=[];fits=[];design=np.zeros((108,20));design[:,0]=1
for i in range(12):
 for j in range(9):
  if i:design[9*i+j,i]=1
  if j:design[9*i+j,11+j]=1
for rep in range(10):
 while True:
  perm=rng.permutation(108);fold=np.empty(108,dtype=int);fold[perm]=np.arange(108)%5;fold=fold.reshape(12,9)
  if all(np.all((fold!=k).sum(axis=0)>=3) and np.all((fold!=k).sum(axis=1)>=3) for k in range(5)):break
 for k in range(5):
  mask=fold!=k;mean=np.array([np.mean(M[mask[:,j],j]) for j in range(9)]);coef=np.linalg.lstsq(design[mask.ravel()],M.ravel()[mask.ravel()],rcond=None)[0];add=(design@coef).reshape(12,9)
  trials=[fit(mask,9000+rep*50+k*10+i) for i in range(5)];best=min(range(5),key=lambda i:trials[i][0]);pred=trials[best][1]
  for index,(obj,_) in enumerate(trials):fits.append({'repeat':rep,'fold':k,'initialization':index,'seed':9000+rep*50+k*10+index,'objective':obj,'selected':index==best})
  for i,j in zip(*np.where(~mask)):records.append({'repeat':rep,'fold':k,'sleep_trait':sleeps[i],'outcome':traits[j],'observed':M[i,j],'column_mean':mean[j],'additive':add[i,j],'rank2_multistart':pred[i,j]})
def write(name,rows):
 with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('masked_predictions.tsv',records);write('training_objectives.tsv',fits)
summary=[]
for group in ['ALL']+traits:
 rows=[x for x in records if group=='ALL' or x['outcome']==group]
 for model in ['column_mean','additive','rank2_multistart']:
  summary.append({'group':group,'model':model,'n_predictions':len(rows),'rmse':math.sqrt(sum((x[model]-x['observed'])**2 for x in rows)/len(rows)),'out_of_range':sum(abs(x[model])>1 for x in rows)})
write('model_summary.tsv',summary)
meta={'n_cells':108,'n_masks':len(records),'n_fits':len(fits),'seed':20261029,'rank2_in_sample_squared_magnitude_fraction':float(sum(d[:2]**2)/sum(d**2)),'rank2_in_sample_rmse':float(np.sqrt(np.mean((M-svd)**2))),'input_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs.values()},'source_version_and_sample_matches':12,'paired_input_bytes':'not established','python':platform.python_version(),'numpy':np.__version__,'caution':'post-hoc cell reconstruction only; not independent cohort or new-phenotype validation'}
(out/'receipt.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2));print('summary', [x for x in summary if x['group'] in ('ALL','FI','FFS','HFRS')])
