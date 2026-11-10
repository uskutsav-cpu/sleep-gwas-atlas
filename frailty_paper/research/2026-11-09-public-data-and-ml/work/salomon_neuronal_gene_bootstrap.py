import runpy, numpy as np, json, hashlib
from pathlib import Path
m=runpy.run_path('work/salomon_neuronal_invariant_transfer.py')
g=np.array(m['group_labels']); scores=m['scores']; y=m['labels']; assert len(g)==len(scores)==len(y) and all(g)
names=np.unique(g); groups=[np.flatnonzero(g==name) for name in names];rng=np.random.default_rng(20316)
def auc(z,v):
 order=np.argsort(-z,kind='stable');s=z[order];t=v[order];n=int(t.sum());neg=len(t)-n
 if not n or not neg:return np.nan
 ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1];tp=np.cumsum(t)[ends];fp=ends+1-tp
 return float(np.trapezoid(np.r_[0,tp/n],np.r_[0,fp/neg]))
point=auc(scores,y); assert abs(point-m['receipt']['auroc'])<1e-12
vals=[]
for i in range(2000):
 draw=rng.integers(0,len(groups),len(groups));ix=np.concatenate([groups[k] for k in draw]);vals.append(auc(scores[ix],y[ix]))
valid=np.array(vals);valid=valid[np.isfinite(valid)]
result={'groups':len(names),'draws':2000,'valid_draws':len(valid),'seed':20316,'original_n':len(y),'point_auroc':point,'percentile_95':np.quantile(valid,[.025,.975]).tolist(),'group_field':'associated_gene_name','limitations':'Assumes reported gene groups are exchangeable independent clusters; does not capture source selection, training uncertainty, label error, shared donor or inter-gene dependence. A percentile range is not proof of no effect.'}
Path('outputs/stage6_salomon_gene_cluster_bootstrap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
