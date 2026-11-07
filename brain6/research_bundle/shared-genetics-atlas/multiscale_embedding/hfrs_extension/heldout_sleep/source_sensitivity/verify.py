import csv,pathlib,json,hashlib,math,numpy as np
p=pathlib.Path(__file__).resolve().parent
def read(f):
 with f.open() as z:return list(csv.DictReader(z,delimiter='\t'))
s=read(p.parent.parent/'source_matrix.tsv');pred=read(p/'predictions.tsv');scores=read(p/'nested_scores.tsv');summary=read(p/'summary.tsv');prev={(x['held_sleep'],x['target']):x for x in read(p.parent/'predictions.tsv')};names=list(dict.fromkeys(x['sleep_trait'] for x in s));outcomes=list(dict.fromkeys(x['outcome'] for x in s));features=[x for x in outcomes if x not in ('FI','FFS','HFRS')];targets=['FI','FFS','HFRS'];M={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in s};X=np.array([[M[(i,j)] for j in features] for i in names]);Y=np.array([[M[(i,j)] for j in targets] for i in names]);assert len(pred)==66 and len(scores)==264 and len(summary)==24
# Use independently expressed augmented least-squares solver.
def est(x,y,v,l):
 mu=x.mean(0);sd=x.std(0);sd=np.where(sd>1e-12,sd,1.);z=(x-mu)/sd;A=np.column_stack([np.ones(len(z)),z]);pen=np.column_stack([np.zeros(6),math.sqrt(l)*np.eye(6)]);coef=np.linalg.lstsq(np.vstack([A,pen]),np.r_[y,np.zeros(6)],rcond=None)[0];return float(np.r_[1,(v-mu)/sd]@coef)
for r in pred:
 ex=names.index(r['excluded']);held=names.index(r['held_sleep']);j=targets.index(r['target']);assert ex!=held;train=np.array([i for i in range(12) if i not in (ex,held)]);assert len(train)==10
 assert math.isclose(float(r['observed']),Y[held,j],rel_tol=1e-12);assert math.isclose(float(r['old_full_training_prediction']),float(prev[(r['held_sleep'],r['target'])]['ridge_prediction']),rel_tol=1e-12);assert math.isclose(float(r['excluded_training_mean']),Y[train,j].mean(),rel_tol=1e-12)
 vals=[x for x in scores if (x['excluded'],x['held_sleep'],x['target'])==(r['excluded'],r['held_sleep'],r['target'])];assert len(vals)==4
 for x in vals:
  lam=float(x['lambda']);errs=[]
  for q in train:
   fit=train[train!=q];errs.append((est(X[fit],Y[fit,j],X[q],lam)-Y[q,j])**2)
  assert math.isclose(sum(errs)/len(errs),float(x['inner_mse']),rel_tol=1e-10,abs_tol=1e-12)
 best=min(vals,key=lambda z:(float(z['inner_mse']),-float(z['lambda'])));assert float(best['lambda'])==float(r['selected_lambda']);assert math.isclose(est(X[train],Y[train,j],X[held],float(best['lambda'])),float(r['excluded_training_prediction']),rel_tol=1e-10,abs_tol=1e-12)
for r in summary:
 z=[x for x in pred if x['excluded']==r['excluded'] and (r['target']=='ALL' or x['target']==r['target'])];v=math.sqrt(sum((float(x[r['model']])-float(x['observed']))**2 for x in z)/len(z));assert len(z)==int(r['n']) and math.isclose(v,float(r['rmse']),rel_tol=1e-12)
meta={'status':'PASS','n_outer':len(pred),'n_inner_scores':len(scores),'method':'independent augmented least squares; nested folds, original predictions and RMSE checked','hashes':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['run.py','predictions.tsv','nested_scores.tsv','summary.tsv']},'caveat':'Arithmetic validation does not establish external generalization or remove source biases.'};(p/'verification.json').write_text(json.dumps(meta,indent=2)+'\n');print('PASS 66 predictions and 264 nested scores')
