"""Paired genomic-block uncertainty for effect ranks, local ranks and abstention."""
import subprocess
from common import *
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
class WeightedRanks:
    def __init__(self,x):
        self.order=np.argsort(x,kind='stable');s=np.asarray(x)[self.order];self.starts=np.r_[0,np.flatnonzero(np.diff(s)!=0)+1];self.sizes=np.diff(np.r_[self.starts,len(x)])
    def ranks(self,w):
        weight=np.add.reduceat(w[self.order],self.starts);middle=np.cumsum(weight)-.5*weight;out=np.empty(len(w));out[self.order]=np.repeat(middle,self.sizes);return out
def correlation(a,b,w):
    total=w.sum();aa=a-np.sum(w*a)/total;bb=b-np.sum(w*b)/total;den=np.sqrt(np.sum(w*aa*aa)*np.sum(w*bb*bb));return float(np.sum(w*aa*bb)/den) if den>0 else np.nan
def summarize(values):
    a=np.asarray(values);ok=np.isfinite(a);return {'valid_repetitions':int(ok.sum()),'repetitions':len(a),'CI95':np.quantile(a[ok],[.025,.975]).tolist() if ok.sum()>=.95*len(a) else None,'qualification':None if ok.sum()>=.95*len(a) else 'Too few valid block resamples; interval not reported.'}
def main():
    report=json.loads((OUT/'reserved_study_results.json').read_text());assert sha(OUT/'reserved_study_predictions.csv')==report['predictions_sha256']
    with np.load(OUT/'reserved_analysis_arrays.npz') as z:d={k:z[k] for k in z.files}
    y=d['label'];effect=d['effect'];models={k.removeprefix('model::'):v for k,v in d.items() if k.startswith('model::')};unique,group=np.unique(d['blocks'],return_inverse=True);nblocks=len(unique);er=WeightedRanks(effect);pr={n:WeightedRanks(p) for n,p in models.items()};cor_draw={n:[] for n in models};local_draw={n:{'Spearman':[],'AUROC':[]} for n in models};tiles={}
    for tile in np.unique(d['tiles']):
        ix=np.flatnonzero(d['tiles']==tile);tileblock=int(np.bincount(group[ix],minlength=nblocks).argmax());tiles[tile]=(ix,tileblock)
    local={}
    for name,p in models.items():
        rhos=[];aucs=[];rg=[];ag=[]
        for tile,(ix,block) in tiles.items():
            if len(ix)>=10 and len(np.unique(effect[ix]))>1 and len(np.unique(p[ix]))>1:rhos.append(float(spearmanr(effect[ix],p[ix]).statistic));rg.append(block)
            if len(ix)>=10 and y[ix].sum()>=2 and (len(ix)-y[ix].sum())>=2:aucs.append(float(roc_auc_score(y[ix],p[ix])));ag.append(block)
        local[name]={'Spearman':(np.array(rhos),np.array(rg,dtype='int32')),'AUROC':(np.array(aucs),np.array(ag,dtype='int32'))}
    ranked=d['status']=='ranked';positive=ranked&d['set_includes1'];negative=ranked&d['set_includes0'];covered=np.where(y==1,d['set_includes1'],d['set_includes0']);correct=(positive&(y==1))|(negative&(y==0));prob=d['calibrated_call_probability']
    def mean(a,w,use=None):
        if use is not None:a=a[use];w=w[use]
        return float(np.sum(w*a)/w.sum()) if w.sum()>0 else np.nan
    def uncertainty(w):return {'coverage':mean(covered,w),'detected_coverage':mean(covered,w,y==1),'non_detected_coverage':mean(covered,w,y==0),'ranked_fraction':mean(ranked,w),'singleton_accuracy':mean(correct,w,ranked),'positive_precision':mean(y,w,positive),'positive_recall':mean(positive,w,y==1),'Brier':mean((prob-y)**2,w)}
    u_draw={n:[] for n in uncertainty(np.ones(len(y)))};rng=np.random.default_rng(20261107)
    for b in range(1000):
        counts=np.bincount(rng.integers(0,nblocks,size=nblocks),minlength=nblocks);w=counts[group].astype('float64');rank_effect=er.ranks(w)
        for name in models:
            value=correlation(rank_effect,pr[name].ranks(w),w);cor_draw[name].append(value)
            if b<3 and name=='primary':
                replicated=np.repeat(np.arange(len(y)),w.astype('int32'));expected=float(spearmanr(effect[replicated],models[name][replicated]).statistic);assert abs(value-expected)<1e-12
            for metric,(values,blocks) in local[name].items():weights=counts[blocks];local_draw[name][metric].append(float(np.sum(weights*values)/weights.sum()) if weights.sum()>0 else np.nan)
        for name,value in uncertainty(w).items():u_draw[name].append(value)
    paired={}
    for name in models:
        if name=='primary':continue
        paired[name]={'effect_Spearman_primary_minus_comparator':summarize(np.array(cor_draw['primary'])-np.array(cor_draw[name])),'within_tile_Spearman_primary_minus_comparator':summarize(np.array(local_draw['primary']['Spearman'])-np.array(local_draw[name]['Spearman'])),'within_tile_AUROC_primary_minus_comparator':summarize(np.array(local_draw['primary']['AUROC'])-np.array(local_draw[name]['AUROC']))}
    src=OUT/'reserved_R_effect_input.csv';dst=OUT/'reserved_R_effect_metrics.csv';usrc=OUT/'reserved_R_uncertainty_input.csv';udst=OUT/'reserved_R_uncertainty_metrics.csv'
    with src.open('w',newline='') as f:w=csv.writer(f);w.writerow(['row_id','effect',*models]);w.writerows(zip(d['design_index'],effect,*models.values()))
    with usrc.open('w',newline='') as f:w=csv.writer(f);w.writerow(['row_id','label','probability','set0','set1','ranked']);w.writerows(zip(d['design_index'],y,prob,d['set_includes0'].astype('int32'),d['set_includes1'].astype('int32'),ranked.astype('int32')))
    subprocess.run(['Rscript',str(WORK/'independent_reserved_additional.R'),str(src),str(dst),str(usrc),str(udst)],check=True,capture_output=True,text=True,env=r_environment());rchecks=[]
    for r in csv.DictReader(dst.open()):
        value=float(spearmanr(effect,models[r['model']]).statistic);error=abs(float(r['Spearman'])-value);assert error<1e-12;rchecks.append({'model':r['model'],'metric':'Spearman','error':error})
    expected=uncertainty(np.ones(len(y)))
    for r in csv.DictReader(udst.open()):
        value=expected[r['metric']]
        if not np.isfinite(value):assert r['value'] in ['NA','NaN',''];continue
        error=abs(float(r['value'])-value);assert error<1e-12;rchecks.append({'metric':r['metric'],'error':error})
    result={'genomic_blocks':nblocks,'repetitions':1000,'seed':20261107,'effect_Spearman_intervals':{n:summarize(v) for n,v in cor_draw.items()},'within_tile_intervals':{n:{m:summarize(v) for m,v in values.items()} for n,values in local_draw.items()},'paired_differences':paired,'uncertainty_intervals':{n:summarize(v) for n,v in u_draw.items()},'independent_R_checks':rchecks,'weighted_rank_bootstrap_materialized_replication_checks':3,'scope':'Within-study genomic-block uncertainty only. No independent donors or disease outcomes, and no model re-fitting. Tile statistics use equal tile weight and a modal genomic block per tile.'}
    write_json(OUT/'reserved_additional_checks.json',result);print('FINAL EFFECT/LOCAL-RANK/ABSTENTION BLOCK CHECKS AND INDEPENDENT R CHECKS PASSED',len(rchecks),flush=True)
if __name__=='__main__':main()
