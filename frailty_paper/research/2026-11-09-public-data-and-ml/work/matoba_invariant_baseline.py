import runpy,json
import numpy as np
m=runpy.run_path('work/matoba_fixed_baseline_transfer.py'); b=m['b']; pairs=m['pairs']; comp=m['comp']
kmers=b['kmers']; idx={k:i for i,k in enumerate(kmers)}
permutation=np.array([idx[k.translate(comp)[::-1]] for k in kmers])
assert np.array_equal(permutation[permutation],np.arange(len(kmers)))
def invariant(x):return (x+x[...,permutation])/2
X=invariant(np.stack(b['xs'])); y=np.array(b['ys']); train=~np.isin(np.array(b['chrs']),[1,2]); test=~train
mean=X[train].mean(0);sd=X[train].std(0);sd[sd==0]=1
X=np.clip((X-mean)/sd,-20,20);weights=np.zeros(X.shape[1]);intercept=0.
for _ in range(4000):
    err=1/(1+np.exp(-np.clip(X[train]@weights+intercept,-30,30)))-y[train]
    weights-=.05*(X[train].T@err/train.sum()+.1*weights);intercept-=.05*err.mean()
def predictions(x):return 1/(1+np.exp(-np.clip(np.clip((x-mean)/sd,-20,20)@weights+intercept,-30,30)))
def metrics(scores,labels):
    pos=scores[labels==1];neg=scores[labels==0]
    auc=float(((pos[:,None]>neg).sum()+.5*(pos[:,None]==neg).sum())/(len(pos)*len(neg)))
    order=np.argsort(-scores,kind='stable'); lab=labels[order]; s=scores[order]; ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1]; tp=np.cumsum(lab)[ends]; prec=tp/(ends+1);rec=tp/len(pos)
    return {'n':len(labels),'detected':int(labels.sum()),'auroc':auc,'average_precision':float(np.sum(np.diff(np.r_[0,rec])*prec))}
raw=np.stack([np.abs(b['kcounts'](a)-b['kcounts'](c)) for a,c,t in pairs]); rev=np.stack([np.abs(b['kcounts'](a.translate(comp)[::-1])-b['kcounts'](c.translate(comp)[::-1])) for a,c,t in pairs]);
if not np.array_equal(invariant(raw),invariant(rev)):raise AssertionError('strand invariance failed')
mat=metrics(predictions(invariant(raw)),np.array([t for a,c,t in pairs])); ret=metrics(1/(1+np.exp(-np.clip(X[test]@weights+intercept,-30,30))),y[test])
r={'approach':'absolute allele k-mer count differences averaged with reverse-complement permutation before standardization','train_n':int(train.sum()),'retallick_chromosome_1_2_test':ret,'matoba_chromosome_1_2_test':mat,'matoba_exact_invariant_feature_check':True,'retallick_original_test_auroc':b['res']['test_auroc'],'matoba_original_released_auroc':m['r']['as_released_orientation']['auroc'],'matoba_original_plus_auroc':m['r']['reference_plus_orientation']['auroc'],'caveat':'Post hoc model redesign after seeing both test sets and orientation sensitivity; no untouched independent test. Detection labels are not functional truth; no clinical genetic score validation.'}
open('outputs/stage6_matoba_invariant_baseline.json','w').write(json.dumps(r,indent=2)+'\n'); print(json.dumps(r,indent=2))
