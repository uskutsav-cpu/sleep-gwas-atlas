"""Main task: raw magnitude error, local ordering and between-cell differences."""
from common import *
from scipy.stats import spearmanr
def rho(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>2 and np.ptp(a)>0 and np.ptp(b)>0 else 0.
def baseline(d,fit,ids):
    main=fit[d['unit'][fit]==0];global_median=float(np.median(d['effect'][main]));means={int(c):float(np.median(d['effect'][main[d['cell'][main]==c]])) for c in np.unique(d['cell'][main])};return np.array([means.get(int(c),global_median) for c in d['cell'][ids]])
def evaluate(d,fit,ids,prediction):
    assert len(ids)==len(prediction) and np.isfinite(prediction).all() and np.all(prediction>=0);main=d['unit'][ids]==0;truth=d['effect'][ids];base=baseline(d,fit,ids);contexts={};skills=[];rhos=[];local=[]
    for key in np.unique(d['context'][ids[main]]):
        use=main&(d['context'][ids]==key);ix=ids[use]
        if len(ix)<30:continue
        y=truth[use];p=prediction[use];b=base[use];mae=float(np.mean(abs(p-y)));bm=float(np.mean(abs(b-y)));r=rho(y,p);skill=1-mae/max(bm,1e-6);skills.append(skill);rhos.append(r);contexts[d['context_names'][int(key)]]={'n':len(ix),'MAE':mae,'RMSE':float(np.sqrt(np.mean((p-y)**2))),'baseline_MAE':bm,'MAE_skill':skill,'effect_Spearman':r}
        g=d['group'][ix];order=np.argsort(g,kind='stable');edges=np.r_[0,np.flatnonzero(np.diff(g[order]))+1,len(g)]
        for lo,hi in zip(edges[:-1],edges[1:]):
            if hi-lo<5:continue
            use=order[lo:hi]
            if np.ptp(y[use])>0:local.append(rho(y[use],p[use]))
    # Packed pairs preserve the same locus and role. Full-length vectors
    # let the paired evaluator check actual differences, never source calls.
    full=np.full(d['n'],np.nan);full[ids]=prediction;bfull=np.full(d['n'],np.nan);bfull[ids]=base;pairs=d['cross_cell_pairs'];valid=np.isfinite(full[pairs[:,0]])&np.isfinite(full[pairs[:,1]]);pairs=pairs[valid];pair_skill=[];pair_results={}
    for study in np.unique(d['study'][pairs[:,0]]) if len(pairs) else []:
        pp=pairs[d['study'][pairs[:,0]]==study];actual=d['effect'][pp[:,0]]-d['effect'][pp[:,1]];estimated=full[pp[:,0]]-full[pp[:,1]];reference=bfull[pp[:,0]]-bfull[pp[:,1]];mae=float(np.mean(abs(actual-estimated)));bm=float(np.mean(abs(actual-reference)));skill=1-mae/max(bm,1e-6);pair_skill.append(skill);pair_results[d['study_names'][int(study)]]={'pairs':len(pp),'magnitude_difference_MAE':mae,'cell_median_difference_MAE':bm,'MAE_skill':skill,'ordering_accuracy':float(np.mean((estimated>0)+.5*(estimated==0)))}
    assert skills;mae_skill=float(np.mean(skills));global_rho=float(np.mean(rhos));local_rho=float(np.mean(local)) if local else 0.;cell_skill=float(np.mean(pair_skill)) if pair_skill else 0.;score=.5*mae_skill+.25*global_rho+.15*local_rho+.1*cell_skill
    return {'selection_score':score,'macro_MAE_skill':mae_skill,'macro_effect_Spearman':global_rho,'mean_within_locus_Spearman':local_rho,'within_locus_groups':len(local),'cross_cell_MAE_skill':cell_skill,'contexts':contexts,'cross_cell_comparisons':pair_results,'primary_observations':int(main.sum()),'scope':'Previously observed/new development experiments with held-out loci. Primary measured log2effect magnitude; auxiliary saturation coefficient units excluded from primary error. Not a new-study or clinical validation.'}
