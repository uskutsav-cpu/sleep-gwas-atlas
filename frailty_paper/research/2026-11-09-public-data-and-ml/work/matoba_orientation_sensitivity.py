import runpy, json
import numpy as np
m=runpy.run_path('work/matoba_locus_bootstrap.py')
b=m['b']; pairs=m['pairs']; y=m['y']; comp=m['comp']; auc=m['auc']; clusters=m['clusters']
def score(reverse):
    features=[]
    for a,c,_ in pairs:
        if reverse:
            a,c=a.translate(comp)[::-1],c.translate(comp)[::-1]
        features.append(np.abs(b['kcounts'](a)-b['kcounts'](c)))
    x=np.clip((np.stack(features)-b['mean'])/b['sd'],-20,20)
    return 1/(1+np.exp(-np.clip(x@b['wgt']+b['intercept'],-30,30)))
a=score(False); z=score(True)
def average_ranks(x):
    order=np.argsort(x,kind='stable'); ranks=np.empty(len(x)); i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]: j+=1
        ranks[order[i:j]]=(i+j-1)/2; i=j
    return ranks
ra,rz=average_ranks(a),average_ranks(z)
rng=np.random.default_rng(20261028);keys=list(clusters);d=[]
for _ in range(5000):
    idx=np.array([i for k in rng.choice(keys,len(keys),replace=True) for i in clusters[k]],int)
    if len(np.unique(y[idx]))==2:d.append(auc(z[idx],y[idx])-auc(a[idx],y[idx]))
r={'n':len(y),'auroc_as_released':auc(a,y),'auroc_reference_plus':auc(z,y),'paired_auroc_difference_reference_minus_released':auc(z,y)-auc(a,y),'paired_locus_bootstrap_95_percentile':list(map(float,np.quantile(d,[.025,.975]))),'valid_draws':len(d),'seed':20261028,'median_absolute_score_difference':float(np.median(abs(z-a))),'maximum_absolute_score_difference':float(np.max(abs(z-a))),'mean_absolute_rank_shift':float(np.mean(abs(rz-ra))),'maximum_absolute_rank_shift':float(np.max(abs(rz-ra))),'spearman_rank_correlation':float(np.corrcoef(ra,rz)[0,1]),'caveats':'Post hoc orientation sensitivity on same selected tests. Locus cluster bootstrap conditions on fixed trained model and assumes exchangeable reported loci; labels are detection calls, not ground truth. Scores are model outputs, not calibrated probabilities or clinical risk.'}
open('outputs/stage6_matoba_orientation_sensitivity.json','w').write(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
