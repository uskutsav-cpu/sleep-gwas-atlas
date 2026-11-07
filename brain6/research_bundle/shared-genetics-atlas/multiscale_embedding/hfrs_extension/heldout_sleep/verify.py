import csv,pathlib,collections,math,hashlib,json,numpy as np
p=pathlib.Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
s=read(p.parent/'source_matrix.tsv');pred=read(p/'predictions.tsv');scores=read(p/'nested_scores.tsv');summ=read(p/'summary.tsv')
sleeps=list(dict.fromkeys(x['sleep_trait'] for x in s));outcomes=list(dict.fromkeys(x['outcome'] for x in s));features=[x for x in outcomes if x not in ('FI','FFS','HFRS')];targets=['FI','FFS','HFRS'];m={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in s};X=np.array([[m[(z,t)] for t in features] for z in sleeps]);Y=np.array([[m[(z,t)] for t in targets] for z in sleeps]);assert len(pred)==36 and len(scores)==144
# Independent augmented least-squares formulation of ridge, intercept unpenalized.
def estimate(Xt,Yt,x,lam):
 mu=Xt.mean(axis=0);sd=Xt.std(axis=0);sd=np.where(sd>1e-12,sd,1.);xx=(Xt-mu)/sd;A=np.column_stack([np.ones(len(xx)),xx]);pen=np.column_stack([np.zeros(6),math.sqrt(lam)*np.eye(6)]);b=np.r_[Yt,np.zeros(6)];coef=np.linalg.lstsq(np.vstack([A,pen]),b,rcond=None)[0];return float(np.r_[1,(x-mu)/sd]@coef)
for row in pred:
 held=sleeps.index(row['held_sleep']);j=targets.index(row['target']);train=np.array([i for i in range(12) if i!=held]);assert math.isclose(float(row['observed']),Y[held,j],rel_tol=1e-12);assert math.isclose(float(row['baseline']),Y[train,j].mean(),rel_tol=1e-12)
 val=[z for z in scores if z['held_sleep']==row['held_sleep'] and z['target']==row['target']];assert len(val)==4
 for z in val:
  lam=float(z['lambda']);errs=[]
  for q in train:
   fit=train[train!=q];errs.append((estimate(X[fit],Y[fit,j],X[q],lam)-Y[q,j])**2)
  assert math.isclose(sum(errs)/len(errs),float(z['inner_mse']),rel_tol=1e-10,abs_tol=1e-12)
 best=min(val,key=lambda z:(float(z['inner_mse']),-float(z['lambda'])));assert float(row['selected_lambda'])==float(best['lambda']);assert math.isclose(estimate(X[train],Y[train,j],X[held],float(best['lambda'])),float(row['ridge_prediction']),rel_tol=1e-10,abs_tol=1e-12)
 mu=X[train].mean(axis=0);sd=X[train].std(axis=0);sd=np.where(sd>1e-12,sd,1.);idx=train[np.argsort(np.sum(((X[train]-X[held])/sd)**2,axis=1),kind='stable')[:3]];assert row['knn3_neighbors']==','.join(sleeps[i] for i in idx);assert math.isclose(float(row['knn3']),Y[idx,j].mean(),rel_tol=1e-12)
for x in summ:
 z=[r for r in pred if x['target']=='ALL' or r['target']==x['target']];rmse=math.sqrt(sum((float(r[x['model']])-float(r['observed']))**2 for r in z)/len(z));assert math.isclose(rmse,float(x['rmse']),rel_tol=1e-12)
report={'status':'PASS','outer_predictions':len(pred),'nested_scores_checked':len(scores),'independent_solver':'augmented least squares','checks':'nested folds, fitted predictions, neighbor selection and summary RMSE','input_files_sha256':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ['run.py','predictions.tsv','nested_scores.tsv','summary.tsv']},'caveat':'No validation of independence, source compatibility, or causal relevance'};(p/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS: all 36 held-out predictions and 144 nested scores independently recomputed')
