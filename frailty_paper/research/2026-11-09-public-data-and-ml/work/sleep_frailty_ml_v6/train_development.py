"""Finite preregistered measured-data search, with real neural backpropagation."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import sys,argparse,subprocess,pickle
from data import *
from metrics import evaluate
from network import initialize,forward,branches,mx
import mlx.optimizers as optim
from sklearn.ensemble import HistGradientBoostingClassifier

RECIPES=[
 {'name':'local_mlp','architecture':'mlp','view':'local','pair_loss':0.,'neural_only':False},
 {'name':'biological_mlp','architecture':'mlp','view':'biological','pair_loss':0.,'neural_only':False},
 {'name':'contrast_standard','architecture':'contrast','view':'biological','pair_loss':0.,'neural_only':False},
 {'name':'contrast_hard','architecture':'contrast','view':'biological','pair_loss':.2,'neural_only':False},
 {'name':'contrast_neural_only','architecture':'contrast','view':'biological','pair_loss':0.,'neural_only':True},
]
def freeze_plan():
    assert not (OUT/'final_test_evaluation_started.json').exists()
    plan={'recipes':RECIPES,'seeds':[42,73,101],'architecture':'MLP128/32 GELU or context128/32 and bias-free contrast128/32, cell-conditioned context modulation; scalar detected-call head and source-specific effect/activity auxiliary heads. Original DNA encoder fixed in this first finite comparison.','optimizer':{'name':'AdamW','learning_rate':.0003,'weight_decay':.01,'gradient_norm_cap':5.,'batch_size':512,'dropout':.1,'max_epochs':40,'patience':6,'minimum_selection_improvement':1e-5},'hard_loss':'0.2 mean softplus(-(positive_logit-negative_logit)),128 training pairs per minibatch with a separate RNG; actual measured strata declared in partition_audit.json.','auxiliary_loss':'0.1 train-context-standardized log1p(abs measured effect) Huber plus0.1 train-context-standardized activity Huber where observed. No outcome is a prediction input.','weights':'Sqrt inverse source/cell-context count, clipped0.2–5 and mean1; sqrt negative/positive class weight per source/cell fitted on training only.','selection':'0.40 mean neural APskill +0.25 mean within-active neural APskill +0.20 mean other-cell APskill +0.15(2*mean hard-experiment ranking accuracy-1). Context metrics require>=20calls and>=20other observations. Neural-only ablation stops on0.65neuralAPskill+0.35within-activeAPskill; excluded from primary-family selection.','trees':[{'depth':3,'iterations':250},{'depth':5,'iterations':250}],'tree_features':'Local sequence1–3mer and declared motif features, plus cell-group indicators. HistGradientBoosting learning_rate0.05,l2=1,min_leaf40,max_bins64,early_stoppingFalse,seed42; same train-derived weights/normalization.','ensemble_rule':'Arithmetic three-seed means for each neural recipe; compare fixed50/50 and75/25 mixtures of the top two primary eligible family ensembles (including tree families) on selection only. Retain all candidates; no final-study-based choices.','final_fit':'Refit selected family/families on train+validation with each member fixed at its own development-selected epoch; trees keep fixed tree count. Separate calibration partition is used only afterward.','uncertainty':'Logistic probability calibration and class-conditional90%prediction sets on the untouched calibration partition; unknown cell groups, ambiguous/empty sets or declared sequence-domain shift return unranked. Report actual final-study coverage/error, not an exchangeability guarantee under shift.','primary_excludes':'All AlphaGenome scores, V4 AVI-taught weights and V3 supervised weights.','final_test':'Entire Kosicki2025 study; labels sealed until final model/config/scorer/uncertainty hashes are frozen.'}
    path=OUT/'model_search_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan
def prepare_runtime():
    d=load()
    for name,h in d['audit']['outputs'].items():assert sha(DATA/name)==h
    d['pair_group']=np.array([r['study']+'|'+r['experiment']+'|'+r['cell'] for r in d['rows']]);del d['rows'],d['bank']
    features=json.loads((OUT/'fixed_features_receipt.json').read_text());assert sha(features['feature_path'])==features['feature_sha256'];bank_x=np.load(features['feature_path'],mmap_mode='r')
    with np.load(DATA/'hard_pairs.npz') as z:pairs=z['pairs']
    assert sha(DATA/'hard_pairs.npz')==json.loads((OUT/'partition_audit.json').read_text())['hard_pairs_sha256']
    return d,features,bank_x,pairs
def standardized_bank(d,receipt,bank_x,recipe,tr,stage):
    cols=columns(receipt,recipe['view']);c,v=branches(cols,receipt['slices'],len(d['cell_names']));mu,sd=statistics(bank_x,d['bank_index'][tr],cols)
    if recipe['architecture']=='contrast':mu[v]=0
    signature=hashlib.sha256((receipt['feature_sha256']+tr.tobytes().hex()+recipe['view']+recipe['architecture']+mu.tobytes().hex()+sd.tobytes().hex()).encode()).hexdigest()
    cache=DATA/'standardized_banks';cache.mkdir(exist_ok=True);p=cache/(stage+'_'+recipe['view']+'_'+recipe['architecture']+'_'+signature[:16]+'.npy');meta=p.with_suffix('.json')
    if p.exists() and meta.exists():assert json.loads(meta.read_text())['signature']==signature
    else:
        x=np.lib.format.open_memmap(p,mode='w+',dtype='float32',shape=(len(bank_x),len(cols)))
        for start in range(0,len(bank_x),1024):x[start:start+1024]=np.clip((bank_x[np.ix_(np.arange(start,min(start+1024,len(bank_x))),cols)]-mu)/sd,-20,20)
        x.flush();del x;write_json(meta,{'signature':signature,'shape':[len(bank_x),len(cols)],'fitted_observations':len(tr),'no_outcome_inputs':True})
    return np.load(p,mmap_mode='r'),cols,mu,sd,c,v
def getx(x,d,ix):return np.c_[x[d['bank_index'][ix]],np.eye(len(d['cell_names']),dtype='float32')[d['cell'][ix]]]
def auxiliary_activity(d,tr):
    target=np.zeros(d['n'],dtype='float32');mask=np.isfinite(d['activity']).astype('float32');metadata={}
    for key in np.unique(d['context'][tr]):
        ids=tr[(d['context'][tr]==key)&np.isfinite(d['activity'][tr])]
        if not len(ids):continue
        mu=float(d['activity'][ids].mean());sd=max(float(d['activity'][ids].std()),1e-4);use=np.flatnonzero((d['context']==key)&np.isfinite(d['activity']));target[use]=np.clip((d['activity'][use]-mu)/sd,-5,5);metadata[str(key)]={'mean':mu,'sd':sd}
    return target,mask,metadata
def fit_neural(d,receipt,bank_x,pairs,recipe,seed,plan,stage='development',fixed_epochs=None):
    architecture=recipe['architecture'];tr=np.flatnonzero(d['split']<=1) if stage=='final' else np.flatnonzero(d['split']==0)
    if recipe['neural_only']:tr=tr[np.array(['neuronal_' in d['cell_names'][i] for i in d['cell'][tr]])]
    va=np.flatnonzero(d['split']==1);x,cols,mu,sd,c,v=standardized_bank(d,receipt,bank_x,recipe,tr,stage)
    rw,pw,eff,weightmeta=weights_and_targets(d,tr);act,actmask,actmeta=auxiliary_activity(d,tr);weighted_positive=float(np.sum(rw[tr]*pw[tr]*d['y'][tr]));weighted_negative=float(np.sum(rw[tr]*(1-d['y'][tr])));prev=np.clip(weighted_positive/(weighted_positive+weighted_negative),1e-5,1-1e-5)
    initial=initialize(x.shape[1]+len(d['cell_names']),architecture,c,v,len(d['study_names']),prev,seed);p={k:mx.array(v) for k,v in initial.items()};c_m=mx.array(c);v_m=mx.array(v);opt=optim.AdamW(.0003,weight_decay=.01);rng=np.random.default_rng(seed+1100);pairrng=np.random.default_rng(seed+3100);trainflag=np.zeros(d['n'],dtype=bool);trainflag[tr]=True;hp=pairs[trainflag[pairs[:,0]]&trainflag[pairs[:,1]]];first_gradient_norm=None
    def huber(a):return mx.where(mx.abs(a)<=1,.5*a*a,mx.abs(a)-.5)
    def loss(pp,xx,ss,yy,ww,positive,ee,aa,am,masks,pairx,pairsid):
        z,e,a=forward(pp,xx,ss,architecture,c_m,v_m,masks);value=mx.mean(ww*((1-yy)*mx.logaddexp(0,z)+positive*yy*mx.logaddexp(0,-z)))+.1*mx.mean(ww*huber(e-ee))+.1*mx.sum(ww*am*huber(a-aa))/mx.maximum(mx.sum(am),1)
        if recipe['pair_loss']:
            hz=forward(pp,pairx,pairsid,architecture,c_m,v_m)[0];half=hz.shape[0]//2;value=value+recipe['pair_loss']*mx.mean(mx.logaddexp(0,-(hz[:half]-hz[half:])))
        return value
    grad=mx.compile(mx.value_and_grad(loss));best=-np.inf;stale=0;best_epoch=0;history=[];started=time.time();models=DATA/'models';models.mkdir(exist_ok=True);path=models/f'{stage}_{recipe["name"]}_{seed}.npz';epochs=fixed_epochs if fixed_epochs is not None else 40
    for epoch in range(1,epochs+1):
        order=rng.permutation(tr);losses=[]
        for start in range(0,len(order),512):
            ids=order[start:start+512];masks=[mx.array(rng.binomial(1,.9,(len(ids),width)).astype('float32')/.9) for width in ([128,32] if architecture=='mlp' else [128,32,128,32])]
            if recipe['pair_loss']:
                assert len(hp);chosen=hp[pairrng.integers(0,len(hp),size=128)];pi=np.r_[chosen[:,0],chosen[:,1]];px=mx.array(getx(x,d,pi));ps=mx.array(d['study'][pi])
            else:px=mx.zeros((2,x.shape[1]+len(d['cell_names'])));ps=mx.zeros((2,),dtype=mx.int32)
            value,g=grad(p,mx.array(getx(x,d,ids)),mx.array(d['study'][ids]),mx.array(d['y'][ids]),mx.array(rw[ids]),mx.array(pw[ids]),mx.array(eff[ids]),mx.array(act[ids]),mx.array(actmask[ids]),masks,px,ps);norm=mx.sqrt(sum(mx.sum(a*a) for a in g.values()))
            if first_gradient_norm is None:first_gradient_norm=float(norm);assert np.isfinite(first_gradient_norm) and first_gradient_norm>0
            scale=mx.minimum(1.,5/(norm+1e-8));p=opt.apply_gradients({k:a*scale for k,a in g.items()},p);mx.eval(p,opt.state);losses.append(float(value));assert np.isfinite(losses[-1])
        if stage=='development':
            pred=[]
            for start in range(0,len(va),1024):
                ids=va[start:start+1024];z=np.array(forward(p,mx.array(getx(x,d,ids)),mx.array(d['study'][ids]),architecture,c_m,v_m)[0]);pred.extend((1/(1+np.exp(-np.clip(z,-30,30)))).tolist())
            pred=np.array(pred);metric=evaluate(d,va,pred,pairs,recipe['neural_only']);score=metric['selection_score'];history.append({'epoch':epoch,'loss':float(np.mean(losses)),'validation':metric})
            improved=score>best+1e-5
        else:metric=None;score=0;history.append({'epoch':epoch,'loss':float(np.mean(losses))});improved=epoch==epochs
        if improved:
            best=score;best_epoch=epoch;stale=0
            np.savez_compressed(path,**{k:np.array(a) for k,a in p.items()},feature_columns=cols,feature_mean=mu,feature_sd=sd,context_columns=c,variant_columns=v,cell_names=np.array(d['cell_names']),study_names=np.array(d['study_names']),recipe_json=np.array(json.dumps(recipe)),selected_epoch=np.array(epoch),fit_row_indices=tr)
            if stage=='development':np.savez_compressed(path.with_name(path.stem+'_validation.npz'),row_index=va,prediction=pred)
        else:stale+=1
        write_json(OUT/f'{stage}_{recipe["name"]}_{seed}_history.json',history)
        print('V5 NEURAL',stage,recipe['name'],seed,'epoch',epoch,'loss',round(float(np.mean(losses)),5),'score',round(score,5),'seconds',round(time.time()-started,1),flush=True)
        if stage=='development' and stale>=6:break
    with np.load(path) as z:changed={k:int(np.count_nonzero(z[k]!=initial[k])) for k in initial}
    assert sum(changed.values())>0
    record={'name':recipe['name'],'seed':seed,'stage':stage,'selected_epoch':best_epoch,'checkpoint':str(path),'checkpoint_sha256':sha(path),'normalization_training_rows':len(tr),'training_hard_pairs':len(hp),'recipe':recipe,'validation':next((r.get('validation') for r in history if r['epoch']==best_epoch),None),'weight_metadata':weightmeta,'activity_metadata':actmeta,'first_gradient_norm':first_gradient_norm,'changed_neural_parameter_scalars':sum(changed.values()),'changed_by_parameter':changed,'seconds':time.time()-started}
    write_json(OUT/f'{stage}_{recipe["name"]}_{seed}_record.json',record);return record
def fit_tree(d,receipt,bank_x,pairs,depth,stage='development'):
    recipe={'name':f'tree_depth{depth}','view':'pattern_motif','architecture':'mlp'};tr=np.flatnonzero(d['split']<=1) if stage=='final' else np.flatnonzero(d['split']==0);va=np.flatnonzero(d['split']==1);x,cols,mu,sd,c,v=standardized_bank(d,receipt,bank_x,recipe,tr,stage);rw,pw,eff,metadata=weights_and_targets(d,tr);model=HistGradientBoostingClassifier(max_iter=250,max_depth=depth,learning_rate=.05,l2_regularization=1.,min_samples_leaf=40,max_bins=64,early_stopping=False,random_state=42);train=getx(x,d,tr);model.fit(train,d['y'][tr],sample_weight=rw[tr]*np.where(d['y'][tr]==1,pw[tr],1.));del train
    path=DATA/'models'/f'{stage}_tree_depth{depth}.pkl'
    with path.open('wb') as f:pickle.dump({'model':model,'columns':cols,'mean':mu,'sd':sd,'cell_names':d['cell_names'],'fit_row_indices':tr},f)
    metric=None
    if stage=='development':
        p=model.predict_proba(getx(x,d,va))[:,1];metric=evaluate(d,va,p,pairs);np.savez_compressed(path.with_name(path.stem+'_validation.npz'),row_index=va,prediction=p)
    record={'name':recipe['name'],'stage':stage,'checkpoint':str(path),'checkpoint_sha256':sha(path),'validation':metric,'depth':depth};write_json(OUT/f'{stage}_{recipe["name"]}_record.json',record);print('V5 TREE',depth,metric,flush=True);return record
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze-plan-only',action='store_true');parser.add_argument('--wait-for-features-pid',type=int);args=parser.parse_args();plan=freeze_plan()
    if args.freeze_plan_only:print('Finite plan frozen before fitting',sha(OUT/'model_search_plan.json'));return
    if args.wait_for_features_pid:
        while True:
            check=subprocess.run(['ps','-p',str(args.wait_for_features_pid),'-o','command='],text=True,capture_output=True)
            if check.returncode or 'sleep_frailty_ml_v6/prepare_features.py' not in check.stdout:break
            time.sleep(10)
    assert not (OUT/'final_test_evaluation_started.json').exists();d,receipt,bank_x,pairs=prepare_runtime();mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2)
    marker=OUT/'model_search_started.json'
    if not marker.exists():write_json(marker,{'plan_sha256':sha(OUT/'model_search_plan.json'),'feature_receipt_sha256':sha(OUT/'fixed_features_receipt.json'),'partition_audit_sha256':sha(OUT/'partition_audit.json'),'source_sha256':sha(Path(__file__)),'final_labels_opened':False})
    else:
        m=json.loads(marker.read_text());assert m['plan_sha256']==sha(OUT/'model_search_plan.json') and m['source_sha256']==sha(Path(__file__))
    records=[]
    for recipe in RECIPES:
        for seed in plan['seeds']:
            path=OUT/f'development_{recipe["name"]}_{seed}_record.json'
            if path.exists():r=json.loads(path.read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256']
            else:r=fit_neural(d,receipt,bank_x,pairs,recipe,seed,plan)
            records.append(r);write_json(OUT/'development_records.json',records);mx.metal.clear_cache()
    for tree in plan['trees']:
        path=OUT/f'development_tree_depth{tree["depth"]}_record.json'
        if path.exists():r=json.loads(path.read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256']
        else:r=fit_tree(d,receipt,bank_x,pairs,tree['depth'])
        records.append(r);write_json(OUT/'development_records.json',records)
    print('FINITE DEVELOPMENT FITS COMPLETE',len(records),flush=True)
if __name__=='__main__':main()
