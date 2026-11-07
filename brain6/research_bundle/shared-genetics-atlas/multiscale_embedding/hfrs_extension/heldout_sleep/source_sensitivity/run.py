import csv,json,hashlib,math,pathlib,numpy as np
p=pathlib.Path(__file__).resolve().parent;source=p.parent.parent/'source_matrix.tsv';prev=p.parent/'predictions.tsv'
def read(f):
 with f.open() as z:return list(csv.DictReader(z,delimiter='\t'))
rows=read(source);original=read(prev);sleeps=list(dict.fromkeys(x['sleep_trait'] for x in rows));outs=list(dict.fromkeys(x['outcome'] for x in rows));features=[x for x in outs if x not in ('FI','FFS','HFRS')];targets=['FI','FFS','HFRS'];m={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in rows};X=np.array([[m[(s,t)] for t in features] for s in sleeps]);Y=np.array([[m[(s,t)] for t in targets] for s in sleeps]);original={(r['held_sleep'],r['target']):r for r in original};assert len(sleeps)==12 and len(features)==6 and len(original)==36
lams=[.1,1.,10.,100.]
def predict(x,y,test,lam):
 mu=x.mean(0);scale=x.std(0);scale=np.where(scale>1e-12,scale,1.);a=(x-mu)/scale;return float(y.mean()+((test-mu)/scale)@np.linalg.solve(a.T@a+lam*np.eye(6),a.T@(y-y.mean())))
results=[];scores=[]
for excluded in ('longsleep','sleep_apnea'):
 available=np.array([i for i,s in enumerate(sleeps) if s!=excluded]);assert len(available)==11
 for held in available:
  train=available[available!=held];assert len(train)==10
  for j,trait in enumerate(targets):
   mse=[]
   for lam in lams:
    errors=[]
    for val in train:
     fit=train[train!=val];errors.append((predict(X[fit],Y[fit,j],X[val],lam)-Y[val,j])**2)
    score=float(np.mean(errors));mse.append(score);scores.append({'excluded':excluded,'held_sleep':sleeps[held],'target':trait,'lambda':lam,'inner_mse':score})
   best=min(range(4),key=lambda k:(mse[k],-lams[k]));old=original[(sleeps[held],trait)];assert float(old['observed'])==Y[held,j]
   results.append({'excluded':excluded,'held_sleep':sleeps[held],'target':trait,'observed':Y[held,j],'old_full_training_prediction':float(old['ridge_prediction']),'excluded_training_prediction':predict(X[train],Y[train,j],X[held],lams[best]),'excluded_training_mean':float(Y[train,j].mean()),'selected_lambda':lams[best]})
def write(name,rows):
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('predictions.tsv',results);write('nested_scores.tsv',scores);summary=[]
for ex in ('longsleep','sleep_apnea'):
 for trait in ['ALL']+targets:
  z=[r for r in results if r['excluded']==ex and (trait=='ALL' or r['target']==trait)]
  for model in ('old_full_training_prediction','excluded_training_prediction','excluded_training_mean'):
   e=[r[model]-r['observed'] for r in z];summary.append({'excluded':ex,'target':trait,'model':model,'n':len(z),'rmse':math.sqrt(sum(v*v for v in e)/len(e)),'mae':sum(abs(v) for v in e)/len(e)})
write('summary.tsv',summary);receipt={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'previous_predictions_sha256':hashlib.sha256(prev.read_bytes()).hexdigest(),'heldout_predictions':len(results),'nested_scores':len(scores),'features':features,'excluded_traits':['longsleep','sleep_apnea'],'caveat':'post-hoc row sensitivity, not independent validation or resolution of source QC'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('outputs',len(results),len(scores))
