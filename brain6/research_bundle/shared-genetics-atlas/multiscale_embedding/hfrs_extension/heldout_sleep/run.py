import csv,hashlib,json,math,pathlib,numpy as np
p=pathlib.Path(__file__).resolve().parent;src=p.parent/'source_matrix.tsv'
with src.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
sleeps=list(dict.fromkeys(r['sleep_trait'] for r in rows));traits=list(dict.fromkeys(r['outcome'] for r in rows));features=[t for t in traits if t not in ('FI','FFS','HFRS')];targets=['FI','FFS','HFRS'];assert len(sleeps)==12 and len(features)==6 and len(rows)==108
x=np.empty((12,6));y=np.empty((12,3));seen=set()
for r in rows:
 i=sleeps.index(r['sleep_trait']);o=r['outcome'];assert (i,o) not in seen;seen.add((i,o));v=float(r['rg']);
 if o in features:x[i,features.index(o)]=v
 else:y[i,targets.index(o)]=v
assert len(seen)==108 and np.all(np.isfinite(x)) and np.all(np.isfinite(y))
lambdas=[.1,1.,10.,100.]
def pred(train_x,train_y,test_x,lam):
 mu=train_x.mean(axis=0);scale=train_x.std(axis=0);scale=np.where(scale>1e-12,scale,1.);a=(train_x-mu)/scale;b=(test_x-mu)/scale;ym=train_y.mean();coef=np.linalg.solve(a.T@a+lam*np.eye(a.shape[1]),a.T@(train_y-ym));return ym+b@coef
out=[];inner=[]
for held in range(12):
 train=np.array([i for i in range(12) if i!=held]);baseline=y[train].mean(axis=0)
 mu=x[train].mean(axis=0);scale=x[train].std(axis=0);scale=np.where(scale>1e-12,scale,1.);distance=np.sum(((x[train]-x[held])/scale)**2,axis=1);neighbors=train[np.argsort(distance,kind='stable')[:3]];knn=y[neighbors].mean(axis=0)
 for target in range(3):
  scores=[]
  for lam in lambdas:
   errors=[]
   for validation in train:
    fitrows=train[train!=validation];estimate=float(pred(x[fitrows],y[fitrows,target],x[validation],lam));errors.append((estimate-y[validation,target])**2)
   score=float(np.mean(errors));scores.append(score);inner.append({'held_sleep':sleeps[held],'target':targets[target],'lambda':lam,'inner_mse':score})
  best=min(range(len(lambdas)),key=lambda j:(scores[j],-lambdas[j]));prediction=float(pred(x[train],y[train,target],x[held],lambdas[best]));out.append({'held_sleep':sleeps[held],'target':targets[target],'observed':y[held,target],'baseline':baseline[target],'knn3':knn[target],'knn3_neighbors':','.join(sleeps[i] for i in neighbors),'ridge_prediction':prediction,'selected_lambda':lambdas[best],'inner_mse':scores[best]})
def write(name,z):
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,list(z[0]),delimiter='\t');w.writeheader();w.writerows(z)
write('predictions.tsv',out);write('nested_scores.tsv',inner)
summary=[]
for target in ['ALL']+targets:
 rr=[z for z in out if target=='ALL' or z['target']==target]
 for name in ['baseline','knn3','ridge_prediction']:
  err=np.array([r[name]-r['observed'] for r in rr]);summary.append({'target':target,'model':name,'n':len(rr),'rmse':float(np.sqrt(np.mean(err**2))),'mae':float(np.mean(abs(err))),'out_of_range':int(sum(abs(r[name])>1 for r in rr))})
write('summary.tsv',summary)
receipt={'matrix_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'n_sleep_traits':len(sleeps),'features':features,'targets':targets,'heldout_predictions':len(out),'additional_posthoc_baseline':'3 nearest sleep traits by standardized brain-correlation features, equal weights','lambda_grid':lambdas,'normalization':'training-fold mean and population SD; zero SD -> 1','selection':'inner leave-one-training-row-out MSE, ties select larger lambda','caveat':'six observed brain correlations required for held-out sleep; source cohorts reused; not external-cohort validation'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summary,indent=2))
