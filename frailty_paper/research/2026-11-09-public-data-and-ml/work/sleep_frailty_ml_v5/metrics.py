"""Balanced development measures; all final-study outcomes remain excluded."""
from common import *
from sklearn.metrics import roc_auc_score,average_precision_score
def point(y,p):
    if len(y)==0 or len(np.unique(y))<2:return None
    ap=float(average_precision_score(y,p));rate=float(np.mean(y))
    return {'n':len(y),'calls':int(np.sum(y)),'AUROC':float(roc_auc_score(y,p)),'AP':ap,'prevalence':rate,'APskill':(ap-rate)/(1-rate)}
def evaluate(d,ids,p,pairs,neural_only=False):
    assert len(ids)==len(p) and np.isfinite(p).all();out={};neur=[];other=[];conditional=[]
    for key in np.unique(d['context'][ids]):
        use=d['context'][ids]==key;y=d['y'][ids[use]]
        if len(y)<30 or y.sum()<20 or len(y)-y.sum()<20:continue
        m=point(y,p[use]);out[str(key)]=m
        if 'neuronal_' in str(key):
            neur.append(m['APskill']);a=use&(d['active'][ids]==1);q=point(d['y'][ids[a]],p[a])
            if q is not None and q['calls']>=20:conditional.append(q['APskill']);out[str(key)+'|active']=q
        else:other.append(m['APskill'])
    global_p=np.full(d['n'],np.nan);global_p[ids]=p;valid=np.isfinite(global_p[pairs[:,0]])&np.isfinite(global_p[pairs[:,1]]);hp=pairs[valid];hard=[]
    if len(hp):
        groups=d['pair_group'][hp[:,0]]
        for group in np.unique(groups):
            use=groups==group;a=global_p[hp[use,0]];b=global_p[hp[use,1]];hard.append(float(np.mean((a>b)+.5*(a==b))))
    assert neur and conditional
    if neural_only:selection=.65*np.mean(neur)+.35*np.mean(conditional)
    else:
        assert other and hard
        selection=.4*np.mean(neur)+.25*np.mean(conditional)+.2*np.mean(other)+.15*(2*np.mean(hard)-1)
    return {'selection_score':float(selection),'contexts':out,'mean_neural_APskill':float(np.mean(neur)),'mean_active_neural_APskill':float(np.mean(conditional)),'mean_other_cell_APskill':float(np.mean(other)) if other else None,'hard_pairs':len(hp),'hard_experiment_groups':len(hard),'mean_hard_group_accuracy':float(np.mean(hard)) if hard else None,'scope':'Previously observed and new development sources, locus-held-out; not final-study validation.'}
