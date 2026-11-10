"""Finite, continuous-outcome experiment. Final-study outcomes are never loaded."""
import argparse,gc
from data import *
from network import initialize,forward,loss,mx
from metrics import evaluate,baseline

def freeze_plan():
    path=OUT/'model_search_plan.json'
    if path.exists():return json.loads(path.read_text())
    assert not (OUT/'model_search_started.json').exists()
    assert json.loads((OUT/'feature_and_gradient_audit.json').read_text())['passed']
    plan={'created_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
      'seeds':[42,73,101],'recipes':{'cell_effect':{'pair_weight':0.,'no_cell':False},'cell_effect_pairs':{'pair_weight':.15,'no_cell':False},'no_cell_pairs':{'pair_weight':.15,'no_cell':True}},
      'max_epochs':24,'patience':5,'batch_size':1024,'pair_batch':128,'learning_rate':.0003,'weight_decay':.01,'gradient_clip':5.,'dropout':.1,
      'ridge_penalties':[100.,1000.,10000.],'tree_depths':[3,5],'tree_iterations':160,'tree_min_leaf':50,'tree_l2':10.,
      'compact_columns':'native context first32; pattern51 all168; variant motif155; local motif mean155; all cell flags',
      'selection':'.5 macro raw-MAE skill + .25 macro magnitude Spearman + .15 local Spearman + .1 paired-cell raw-MAE skill; no_cell is ablation only; compare all other families and top2 mixtures .5/.75',
      'AVI':'Matched added-input version of selected family/families with identical seeds, rows and selected epoch counts; no added-input-specific search. Permissive AVI composite only.',
      'AVI_retention':{'minimum_relative_macro_MAE_improvement':.01,'paired_locus_bootstrap_MAE_improvement_95CI_lower_above':0.,'maximum_per_cell_relative_MAE_harm':.05,'maximum_local_Spearman_drop':.01,'scope':'Conditional recommendation after a prespecified whole-study comparison; no test-dependent refitting. AVI includes other predictors and cannot isolate AlphaGenome contribution.'},
      'uncertainty':{'method':'Leave-one-primary-study-out residuals on calibration loci. Hash-partition calibration loci into scale-fitting A and quantile B; scale is group median absolute error on A, fallback global median. B90% absolute standardized-error quantile separately by omitted study; use worst study quantile. No formal exchangeability guarantee across studies.','minimum_cell_group_training_studies':2,'domain_max_abs_standardized_feature':12.,'maximum_interval_width':'Twice training90%quantile of primary effect magnitude; zero-truncated intervals.','unknown_cell_group':'Unranked','one_study_cell_group':'Unranked','risk_reporting':'Outside-study coverage, width, accepted fraction and accepted error; never call all-abstention a solved system.'},
      'study_folds':'Omit each of five primary sources; retain auxiliary coefficient head and omit no held-out source rows; fixed final recipe and epochs. Calibration groups absent from all fits.',
      'outside_test':'McAfee2023 both cells, all9602source-verified effect measurements;4820design pairs. Absolute original fitted log2effect; raw counts not independently reconciled.',
      'hashes':{name:sha(OUT/name) for name in ['feature_receipt.json','development_data_audit.json','compact_data_and_pair_receipt.json','reserved_design_and_original_source_audit.json','feature_and_gradient_audit.json']},
      'previous_V6_lock':assert_previous_intact(),'final_test_performance_seen':False}
    write_json(path,plan);return plan

def transform_save(path,t):
    arrays={k:v for k,v in t.items() if isinstance(v,np.ndarray)};np.savez_compressed(path,**arrays)
    write_json(path.with_suffix('.json'),{k:v for k,v in t.items() if not isinstance(v,np.ndarray)})
def transform_load(path):
    with np.load(path) as z:t={k:z[k] for k in z.files}
    t.update(json.loads(path.with_suffix('.json').read_text()));return t
def predict(d,ids,t,p,no_cell=False,avi=None):
    parts=[];pm={k:mx.array(v) for k,v in p.items()}
    for start in range(0,len(ids),2048):
        ix=ids[start:start+2048];c,v=features(d,ix,t,no_cell,avi);z=np.array(forward(pm,mx.array(c),mx.array(v)));a=z[np.arange(len(ix)),d['unit'][ix]]*t['effect_scales'][d['unit'][ix]];parts.append(np.expm1(np.clip(a,0,12)))
    return np.concatenate(parts)

def fit_nn(d,fit,va,recipe,seed,stage='development',epochs=None,avi=None,omit=None):
    plan=freeze_plan();name=f'{stage}_{recipe}_{seed}'+('_avi' if avi is not None else '')+(f'_omit{omit}' if omit is not None else '');record=OUT/(name+'_record.json');path=MODELS/(name+'.npz')
    if record.exists():
        r=json.loads(record.read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256'];return r
    config=plan['recipes'][recipe];t=fit_transform(d,fit)
    if avi is not None:t=add_avi_transform(d,fit,t,avi)
    tx=MODELS/(name+'_transform.npz');transform_save(tx,t)
    c,v=features(d,fit[:1],t,config['no_cell'],avi);p={k:mx.array(a) for k,a in initialize(c.shape[1],v.shape[1],seed,avi is not None).items()};moment={k:mx.zeros_like(a) for k,a in p.items()};variance={k:mx.zeros_like(a) for k,a in p.items()}
    selected=np.zeros(d['n'],dtype=bool);selected[fit]=True
    lp=d['local_pairs'];lp=lp[selected[lp[:,0]]&selected[lp[:,1]]];cp=d['cross_cell_pairs'];cp=cp[selected[cp[:,0]]&selected[cp[:,1]]]
    rng=np.random.default_rng(seed);prng=np.random.default_rng(seed+700);step=0;curve=[];best=-np.inf;best_epoch=0;last=time.time();vag=mx.value_and_grad(loss)
    limit=epochs if epochs is not None else plan['max_epochs'];start_time=time.time()
    for epoch in range(1,limit+1):
        order=rng.permutation(fit);losses=[]
        for start in range(0,len(order),plan['batch_size']):
            main=order[start:start+plan['batch_size']];nl=nc=0;parts=[main]
            if config['pair_weight']:
                if len(lp):a=lp[prng.integers(len(lp),size=plan['pair_batch'])];parts.append(a.ravel());nl=len(a)
                if len(cp):a=cp[prng.integers(len(cp),size=plan['pair_batch'])];parts.append(a.ravel());nc=len(a)
            ix=np.concatenate(parts);c,v=features(d,ix,t,config['no_cell'],avi);u=d['unit'][ix].astype('int32');target=d['log_effect'][ix]/t['effect_scales'][u];weights=t['weights'][ix];masks=[]
            # Identical main-row dropout streams in matched experiments.
            for width in [128,64,128,64]:
                mask=np.ones((len(ix),width),dtype='float32');mask[:len(main)]=(rng.random((len(main),width))>=plan['dropout'])/(1-plan['dropout']);masks.append(mx.array(mask))
            val,grads=vag(p,mx.array(c),mx.array(v),mx.array(u),mx.array(target),mx.array(weights),masks,len(main),nl,nc,config['pair_weight'],config.get('raw_mae',False),mx.array(t['effect_scales']));norm=mx.sqrt(sum(mx.sum(g*g) for g in grads.values()));factor=mx.minimum(1.,plan['gradient_clip']/(norm+1e-8));step+=1
            for k in p:
                g=grads[k]*factor;moment[k]=.9*moment[k]+.1*g;variance[k]=.999*variance[k]+.001*g*g;p[k]=p[k]*(1-plan['learning_rate']*plan['weight_decay'])-plan['learning_rate']*(moment[k]/(1-.9**step))/(mx.sqrt(variance[k]/(1-.999**step))+1e-8)
            mx.eval(p,moment,variance,val);losses.append(float(val))
        cpu={k:np.array(a) for k,a in p.items()}
        row={'epoch':epoch,'training_objective':float(np.mean(losses)),'seconds':time.time()-last};last=time.time()
        if epochs is None:
            prediction=predict(d,va,t,cpu,config['no_cell'],avi);m=evaluate(d,fit,va,prediction);row.update({'selection_score':m['selection_score'],'macro_MAE_skill':m['macro_MAE_skill'],'macro_effect_Spearman':m['macro_effect_Spearman'],'mean_within_locus_Spearman':m['mean_within_locus_Spearman'],'cross_cell_MAE_skill':m['cross_cell_MAE_skill']})
            if m['selection_score']>best+1e-5:
                best=m['selection_score'];best_epoch=epoch;np.savez_compressed(path,**cpu);np.savez_compressed(MODELS/(name+'_validation.npz'),row_index=va,prediction=prediction);write_json(OUT/(name+'_metrics.json'),m)
        else:best_epoch=epoch;np.savez_compressed(path,**cpu)
        curve.append(row);write_json(OUT/(name+'_curve.json'),curve);print(name,json.dumps(row),flush=True)
        if epochs is None and epoch-best_epoch>=plan['patience']:break
    r={'name':name,'recipe':recipe,'seed':seed,'stage':stage,'selected_epoch':best_epoch,'total_epochs':epoch,'seconds':time.time()-start_time,'fit_rows':len(fit),'fit_loci':len(np.unique(d['group'][fit])),'training_loss':'raw_magnitude_MAE' if config.get('raw_mae') else 'scaled_log1p_Huber','plan_sha256':sha(OUT/'model_search_plan.json'),'AVI_added_input':avi is not None,'omitted_study':omit,'checkpoint':str(path),'checkpoint_sha256':sha(path),'transform':str(tx),'transform_sha256':sha(tx),'curve_sha256':sha(OUT/(name+'_curve.json')),'actual_backprop_updates':step,'final_study_loaded':False};write_json(record,r);mx.metal.clear_cache();return r

def main():
    plan=freeze_plan();write_json(OUT/'model_search_started.json',{'plan_sha256':sha(OUT/'model_search_plan.json'),'time':time.time(),'final_study_outcomes_loaded':False})
    d=load();fit=np.flatnonzero(d['split']==0);va=np.flatnonzero(d['split']==1);mx.metal.set_memory_limit(1024*1024**2);mx.metal.set_cache_limit(64*1024**2)
    for recipe in plan['recipes']:
        for seed in plan['seeds']:fit_nn(d,fit,va,recipe,seed)
    print('ALL NINE DEVELOPMENT NEURAL FITS COMPLETE',flush=True)
if __name__=='__main__':main()
