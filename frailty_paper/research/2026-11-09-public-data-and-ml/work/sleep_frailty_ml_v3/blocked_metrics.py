"""Fast, tied-score-correct weighted ranking and paired genomic bootstrap."""
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
class RankedMetrics:
    def __init__(self,y,p):
        y=np.asarray(y,dtype='float64');p=np.asarray(p,dtype='float64');assert y.shape==p.shape and np.isfinite(p).all() and set(np.unique(y))<={0.,1.};self.order=np.argsort(p,kind='stable');self.y=y[self.order];scores=p[self.order];self.starts=np.r_[0,np.flatnonzero(np.diff(scores)!=0)+1]
    def evaluate(self,weights):
        w=np.asarray(weights,dtype='float64')[self.order];tp=np.add.reduceat(w*self.y,self.starts);fp=np.add.reduceat(w*(1-self.y),self.starts);P=tp.sum();N=fp.sum()
        if P==0 or N==0:return np.array([np.nan,np.nan])
        auc=np.sum(tp*(np.cumsum(fp)-.5*fp))/(P*N);pos=np.cumsum(tp[::-1]);total=np.cumsum(tp[::-1]+fp[::-1]);precision=np.divide(pos,total,out=np.zeros_like(pos),where=total>0);ap=np.sum(tp[::-1]/P*precision);return np.array([auc,ap])
def paired_block_bootstrap(y,predictions,blocks,repetitions=1000,seed=2311):
    y=np.asarray(y,dtype='int32');blocks=np.asarray(blocks);_,group=np.unique(blocks,return_inverse=True);nblocks=int(group.max())+1;prepared={k:RankedMetrics(y,p) for k,p in predictions.items()};draws={k:np.empty((repetitions,2)) for k in predictions};rng=np.random.default_rng(seed)
    for b in range(repetitions):
        counts=np.bincount(rng.integers(0,nblocks,size=nblocks),minlength=nblocks);w=counts[group]
        for name,ranking in prepared.items():
            value=ranking.evaluate(w);draws[name][b]=value
            if b<3 and np.isfinite(value).all():
                expected=[roc_auc_score(y,predictions[name],sample_weight=w),average_precision_score(y,predictions[name],sample_weight=w)];assert np.max(np.abs(value-expected))<1e-12,(name,value,expected)
    out={}
    for name,values in draws.items():
        finite=np.isfinite(values).all(1);assert finite.sum()>=.95*repetitions
        out[name]={'valid_repetitions':int(finite.sum()),'auroc_CI95':np.quantile(values[finite,0],[.025,.975]).tolist(),'AP_CI95':np.quantile(values[finite,1],[.025,.975]).tolist()}
    diffs={}
    for name in predictions:
        diffs[name]={}
        for comparator in predictions:
            if name==comparator:continue
            d=draws[name]-draws[comparator];ok=np.isfinite(d).all(1);diffs[name][comparator]={'paired_auroc_difference_CI95':np.quantile(d[ok,0],[.025,.975]).tolist(),'paired_AP_difference_CI95':np.quantile(d[ok,1],[.025,.975]).tolist()}
    return {'block_count':nblocks,'repetitions':repetitions,'seed':seed,'models':out,'paired_differences':diffs,'qualification':'Within-assay genomic-block sampling uncertainty. Does not establish independent donors, disease outcomes or causal effects.'}
def check():
    y=np.array([0,1,1,0,1,0,0,1]);p=np.array([.2,.2,.4,.4,.6,.6,.8,.8]);w=np.array([0,1,2,1,3,0,4,1]);r=RankedMetrics(y,p);actual=r.evaluate(w);expected=np.array([roc_auc_score(y,p,sample_weight=w),average_precision_score(y,p,sample_weight=w)]);assert np.max(np.abs(actual-expected))<1e-12
    result=paired_block_bootstrap(y,{'ties':p,'reverse':1-p},np.array([0,0,1,1,2,2,3,3]),repetitions=200,seed=27);assert result['models']['ties']['valid_repetitions']>=190;return {'weighted_fixture_error':float(np.max(np.abs(actual-expected))),'paired_block_fixture_repetitions':200,'reserved_test_evaluated':False}
if __name__=='__main__':
    import json
    from common import OUT
    r=check();(OUT/'blocked_metrics_implementation_check.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
