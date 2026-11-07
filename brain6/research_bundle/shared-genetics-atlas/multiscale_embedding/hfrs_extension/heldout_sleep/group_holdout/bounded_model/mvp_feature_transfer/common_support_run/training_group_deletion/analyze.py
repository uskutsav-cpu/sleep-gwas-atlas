import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas');repo=atlas.parent.parent/'work/sleep-gwas-atlas/brain6/research_prototypes/frailty_prediction';sys.path.insert(0,str(repo));import predict
matrix=repo/'source_matrix.tsv';features=repo/'mvp_common_support_features.json';prior=p.parent/'predictions.tsv';names,x,y=predict.load_matrix(matrix);v=json.loads(features.read_text());z=np.array([v[k] for k in predict.FEATURES]);groups=list(dict.fromkeys(predict.group(n) for n in names if predict.group(n)!='insomnia'));assert len(groups)==7
with prior.open(newline='') as f:original={r['target']:r for r in csv.DictReader(f,delimiter='\t')}
base=predict.predict(names,x,y,z,'insomnia');assert set(original)==set(predict.TARGETS)
for t in predict.TARGETS:assert abs(base[t]['bounded_ridge']-float(original[t]['mvp_feature_bounded_prediction']))<1e-10
out=[]
for deleted in groups:
 idx=np.array([i for i,n in enumerate(names) if predict.group(n)!=deleted]);subset=[names[i] for i in idx];fit=predict.predict(subset,x[idx],y[idx],z,'insomnia')
 for target in predict.TARGETS:
  prior_prediction=base[target]['bounded_ridge'];pred=fit[target]['bounded_ridge'];observed=float(original[target]['observed_mvp_rg']);out.append({'deleted_training_group':deleted,'target':target,'prediction':pred,'change_from_full_training':pred-prior_prediction,'observed_rg':observed,'absolute_difference_from_observed':abs(pred-observed),'lambda':fit[target]['lambda'],'n_training_rows':len([n for n in subset if predict.group(n)!='insomnia'])})
assert len(out)==21
with (p/'predictions.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
summary={}
for t in predict.TARGETS:
 v=[r['prediction'] for r in out if r['target']==t];summary[t]={'full_training':base[t]['bounded_ridge'],'minimum':min(v),'maximum':max(v),'maximum_absolute_shift':max(abs(n-base[t]['bounded_ridge']) for n in v)}
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(f):sha(f) for f in (matrix,features,prior,repo/'predict.py')},'removed_groups':groups,'rows':len(out),'summary':summary,'scope':'post-hoc training sensitivity, not an independent evaluation or prediction interval'},indent=2)+'\n');print(json.dumps(summary,indent=2))
