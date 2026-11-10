import runpy,json
import numpy as np
v=runpy.run_path('work/matoba_invariant_baseline.py'); t=runpy.run_path('work/matoba_locus_bootstrap.py')
pairs=v['pairs'];comp=v['comp'];b=v['b'];y=t['y'];auc=t['auc'];clusters=t['clusters'];keys=list(clusters)
if len(y)!=len(pairs) or not np.array_equal(y,np.array([p[2] for p in pairs])):raise ValueError('test set order differs')
new=v['predictions'](v['invariant'](v['raw']))
raw_plus=np.stack([np.abs(b['kcounts'](a.translate(comp)[::-1])-b['kcounts'](c.translate(comp)[::-1])) for a,c,_ in pairs])
x=np.clip((raw_plus-b['mean'])/b['sd'],-20,20); old=1/(1+np.exp(-np.clip(x@b['wgt']+b['intercept'],-30,30)))
assert abs(auc(old,y)-t['point'])<1e-12
rng=np.random.default_rng(20261028);draw=[]
for _ in range(5000):
    idx=np.array([i for k in rng.choice(keys,len(keys),replace=True) for i in clusters[k]],int)
    if len(np.unique(y[idx]))==2:draw.append(auc(new[idx],y[idx])-auc(old[idx],y[idx]))
r={'n':len(y),'reported_loci':len(keys),'fixed_baseline_reference_plus_auroc':auc(old,y),'new_invariant_auroc':auc(new,y),'paired_difference':auc(new,y)-auc(old,y),'paired_locus_bootstrap_percentile_95':list(map(float,np.quantile(draw,[.025,.975]))),'valid_draws':len(draw),'seed':20261028,'warning':'Post hoc redesign after test inspection; fixed-model conditional bootstrap assumes exchangeable independent locus labels and does not address training, selection, label or donor uncertainty. Not independent validation or patient score evaluation.'}
open('outputs/stage6_matoba_invariant_paired_bootstrap.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
