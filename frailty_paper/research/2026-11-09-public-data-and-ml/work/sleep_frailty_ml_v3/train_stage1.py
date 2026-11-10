"""Finite multi-study search; no v3 test evaluation or test-driven decisions.

Dense measured activities are targets only. AVI teaching uses a fixed prior
permitted-score snapshot while the broader annotation downloads independently.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
from pathlib import Path
import sys,json,time,subprocess,pickle,argparse
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from common import *
assert_no_reserved_evaluation()

p=argparse.ArgumentParser();p.add_argument('--wait-for-features-pid',type=int);args=p.parse_args()
if args.wait_for_features_pid:
    print('Waiting for verified feature job',args.wait_for_features_pid,flush=True)
    while True:
        state=subprocess.run(['ps','-p',str(args.wait_for_features_pid),'-o','command='],capture_output=True,text=True)
        if state.returncode or 'sleep_frailty_ml_v3/prepare_features.py' not in state.stdout:break
        time.sleep(10)
data,features,manifest,tr,va,te=load();y=data['label'];study=data['study'];MODELS.mkdir(exist_ok=True)
teacher,sources=prior_avi(manifest);np.savez_compressed(WORK/'stage1_teacher_snapshot_local.npz',avi_raw=teacher)
known=np.isfinite(teacher);teacher_mu=float(np.nanmean(teacher[tr]));teacher_sd=float(np.nanstd(teacher[tr]));assert teacher_sd>0
teacher_scaled=np.clip((np.where(known,teacher,teacher_mu)-teacher_mu)/teacher_sd,-5,5).astype('float32')
config={'split':'train 6-22, validation4-5, blind test1-3','selection':'Mean per-study normalized validation AP skill; no test evaluation','studies':STUDIES,'classical_configs':{'logistic_C':[.1,1.,10.],'positive_weighting':['none','sqrt_negative_over_positive'],'boost_depth':[3,5],'boost_learning_rate':.05,'boost_trees':300,'boost_l2':1.,'boost_min_leaf':40,'boost_max_bins':64},'neural_configs':{'widths':[128,32],'activation':'GELU','dropout':.1,'seeds':[42,73,101],'lr':.0003,'weight_decay':.01,'batch':512,'max_epochs':60,'patience':8,'task_weights':{'assay_only':[0,0,0],'activity':[.2,0,0],'activity_effect':[.2,.1,0],'activity_avi':[.2,0,.05]},'study_outputs':'Separate study-specific class, activity and effect heads, shared DNA trunk; AVI head shared and has no AVI input','encoder':'Original frozen HyenaDNA; assay-fitted v1/v2 encoders are not used'},'teacher_sources':sources,'teacher_available_training':{name:int(known[tr[study[tr]==i]].sum()) for i,name in enumerate(STUDIES)},'teacher_scope':'Fixed prior permitted AVI score snapshot. Broader AVI retrieval is still independent; no opportunistic changing teacher labels within this search.','input_features':'DNA-pattern counts/frequencies or original pretrained paired DNA features only. Actual activity/effect/AVI target values are not inference inputs.','preparation_audit_sha256':sha(OUT/'data_audit.json'),'plan_sha256':sha(WORK/'EXPERIMENT_PLAN.md')}
(OUT/'stage1_locked_config.json').write_text(json.dumps(config,indent=2)+'\n');print('LOCKED STAGE1',json.dumps(config),flush=True)
results=[];predictions={};history=[]
def record(name,family,inputs,params,pred,artifact,extra=None):
    row={'name':name,'family':family,'features':inputs,'config':params,'validation':study_metrics(y[va],pred,study[va]),'artifact':artifact}
    if extra:row.update(extra)
    results.append(row);predictions[name]=np.asarray(pred,dtype='float64');(OUT/'stage1_search_results.json').write_text(json.dumps(results,indent=2)+'\n');np.savez_compressed(OUT/'stage1_validation_predictions.npz',row_index=va,label=y[va],study=study[va],**predictions)
    print('VALIDATION',name,json.dumps(row['validation']),flush=True)

# Explicit per-study outputs, not a pooled study-mixture ranking.
for weighting in ['none','sqrt']:
    for C in [.1,1.,10.]:
        name=f'pattern_logistic_{weighting}_C{C:g}';models={};pred=np.empty(len(va),dtype='float64')
        for sid in range(2):
            tt=tr[study[tr]==sid];vv=np.flatnonzero(study[va]==sid);pw=float(np.sqrt((len(tt)-y[tt].sum())/y[tt].sum()));cw=None if weighting=='none' else {0:1.,1:pw}
            fitted=make_pipeline(StandardScaler(),LogisticRegression(C=C,max_iter=1500,class_weight=cw,random_state=42));fitted.fit(features['dna_pattern'][tt],y[tt]);models[sid]=fitted;pred[vv]=fitted.predict_proba(features['dna_pattern'][va[vv]])[:,1]
        path=MODELS/(name+'.pkl')
        with path.open('wb') as f:pickle.dump(models,f)
        record(name,'study_logistic','dna_pattern',{'C':C,'weighting':weighting},pred,str(path.relative_to(OUT)))
    for depth in [3,5]:
        name=f'pattern_boost_{weighting}_depth{depth}';models={};pred=np.empty(len(va),dtype='float64')
        for sid in range(2):
            tt=tr[study[tr]==sid];vv=np.flatnonzero(study[va]==sid);pw=float(np.sqrt((len(tt)-y[tt].sum())/y[tt].sum()));weight=None if weighting=='none' else np.where(y[tt]>0,pw,1.)
            fitted=HistGradientBoostingClassifier(max_depth=depth,learning_rate=.05,max_iter=300,min_samples_leaf=40,l2_regularization=1.,max_bins=64,early_stopping=False,random_state=42);fitted.fit(features['dna_pattern'][tt],y[tt],sample_weight=weight);models[sid]=fitted;pred[vv]=fitted.predict_proba(features['dna_pattern'][va[vv]])[:,1]
        path=MODELS/(name+'.pkl')
        with path.open('wb') as f:pickle.dump(models,f)
        record(name,'study_hist_boosting','dna_pattern',{'depth':depth,'weighting':weighting,'learning_rate':.05,'trees':300,'l2':1.,'min_leaf':40,'max_bins':64},pred,str(path.relative_to(OUT)))

sys.path.insert(0,str(OLD/'runtime'));import mlx.core as mx;import mlx.optimizers as optim
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
f=features['hyena_original'];scaler=StandardScaler().fit(f[tr]);xt=np.clip(scaler.transform(f[tr]),-20,20).astype('float32');xv=np.clip(scaler.transform(f[va]),-20,20).astype('float32');st=study[tr];sv=study[va]
pw=np.array([np.sqrt((sum(st==i)-y[tr[st==i]].sum())/y[tr[st==i]].sum()) for i in range(2)],dtype='float32');source_weights=np.array([len(tr)/(2*sum(st==i)) for i in range(2)],dtype='float32')
targets=np.c_[data['activity_mean'],np.log1p(data['absolute_effect'])].astype('float32');target_mean=np.stack([targets[tr[st==i]].mean(0) for i in range(2)]);target_sd=np.stack([targets[tr[st==i]].std(0) for i in range(2)]);target_sd[target_sd<1e-6]=1.;normalized=np.clip((targets-target_mean[study])/target_sd[study],-5,5).astype('float32')
def gelu(x):return .5*x*(1+mx.tanh(.7978845608028654*(x+.044715*x**3)))
def forward(params,x,sid,m1=None,m2=None):
    h=gelu(x@params['shared.w1']+params['shared.b1']);h=h if m1 is None else h*m1
    h=gelu(h@params['shared.w2']+params['shared.b2']);h=h if m2 is None else h*m2
    row=mx.arange(x.shape[0]);logits=(h@params['class.w']+params['class.b'])[row,sid];activity=(h@params['activity.w']+params['activity.b'])[row,sid];effect=(h@params['effect.w']+params['effect.b'])[row,sid];avi=h@params['avi.w']+params['avi.b']
    return logits,activity,effect,avi
def sigmoid(x):return 1/(1+np.exp(-np.clip(np.asarray(x,dtype='float64'),-30,30)))
def huber(x):return mx.where(mx.abs(x)<=1,.5*x*x,mx.abs(x)-.5)
pw_m=mx.array(pw);sw_m=mx.array(source_weights)
for task,task_weights in config['neural_configs']['task_weights'].items():
    members=[]
    for seed in [42,73,101]:
        rng=np.random.default_rng(seed);name=f'multitask_{task}_s{seed}';members.append(name)
        params={'shared.w1':mx.array(rng.normal(0,1/np.sqrt(512),(512,128)).astype('float32')),'shared.b1':mx.zeros(128),'shared.w2':mx.array(rng.normal(0,1/np.sqrt(128),(128,32)).astype('float32')),'shared.b2':mx.zeros(32)}
        for head in ['class','activity','effect']:
            params[head+'.w']=mx.array(rng.normal(0,.01,(32,2)).astype('float32'));params[head+'.b']=mx.array(np.array([np.log(y[tr[st==i]].mean()/(1-y[tr[st==i]].mean())) for i in range(2)],dtype='float32')) if head=='class' else mx.zeros(2)
        params['avi.w']=mx.zeros(32);params['avi.b']=mx.array(0.);initial={k:np.array(v) for k,v in params.items()};opt=optim.AdamW(.0003,weight_decay=.01)
        def loss(params,x,sid,yy,tt,avi_target,mask,drop1,drop2):
            z,a,e,v=forward(params,x,sid,drop1,drop2);w=sw_m[sid];main=mx.mean(w*((1-yy)*mx.logaddexp(0,z)+pw_m[sid]*yy*mx.logaddexp(0,-z)))
            return main+task_weights[0]*mx.mean(w*huber(a-tt[:,0]))+task_weights[1]*mx.mean(w*huber(e-tt[:,1]))+task_weights[2]*mx.sum(w*mask*huber(v-avi_target))/mx.maximum(mx.sum(w*mask),1.)
        grad=mx.compile(mx.value_and_grad(loss));best=-np.inf;stale=0
        for epoch in range(1,61):
            order=rng.permutation(len(tr));losses=[]
            for start in range(0,len(order),512):
                j=order[start:start+512];ix=tr[j];m1=rng.binomial(1,.9,(len(j),128)).astype('float32')/.9;m2=rng.binomial(1,.9,(len(j),32)).astype('float32')/.9
                value,g=grad(params,mx.array(xt[j]),mx.array(st[j]),mx.array(y[ix]),mx.array(normalized[ix]),mx.array(teacher_scaled[ix]),mx.array(known[ix].astype('float32')),mx.array(m1),mx.array(m2));norm=mx.sqrt(sum(mx.sum(v*v) for v in g.values()));scale=mx.minimum(1.,5/(norm+1e-8));g={k:v*scale for k,v in g.items()};params=opt.apply_gradients(g,params);mx.eval(params,opt.state);losses.append(float(value));assert np.isfinite(losses[-1])
            output=forward(params,mx.array(xv),mx.array(sv));pred=sigmoid(np.array(output[0]));m=study_metrics(y[va],pred,sv);row={'name':name,'epoch':epoch,'train_joint_loss':float(np.mean(losses)),'validation':m}
            if task_weights[0]>0:
                pa=np.array(output[1])*target_sd[sv,0]+target_mean[sv,0];row['validation_activity_R2']={STUDIES[i]:float(r2_score(data['activity_mean'][va[sv==i]],pa[sv==i])) for i in range(2)}
            history.append(row);(OUT/'stage1_neural_history.json').write_text(json.dumps(history,indent=2)+'\n')
            if m['macro_AP_skill']>best+1e-6:
                best=m['macro_AP_skill'];stale=0;best_epoch=epoch;best_pred=pred.copy();saved={k:np.array(v) for k,v in params.items()}
            else:stale+=1
            if epoch%5==0:print('Neural',name,epoch,'validation macro skill',m['macro_AP_skill'],'best',best,flush=True)
            if stale>=8:break
        path=MODELS/(name+'.npz');np.savez_compressed(path,**saved,feature_mean=scaler.mean_,feature_sd=scaler.scale_,target_mean=target_mean,target_sd=target_sd,teacher_mean=teacher_mu,teacher_sd=teacher_sd)
        changes={k:float(np.linalg.norm(saved[k]-initial[k])) for k in params};assert changes['shared.w1']>0 and changes['shared.w2']>0 and changes['class.w']>0
        record(name,'frozen_multi_study_head','hyena_original',{'task':task,'task_weights':task_weights,'seed':seed,'best_epoch':best_epoch,'native_encoder_frozen':True},best_pred,str(path.relative_to(OUT)),{'parameter_change_norms':changes,'teacher_available_training':config['teacher_available_training']})
    record('ensemble_'+task,'seed_ensemble','hyena_original',{'members':members,'rule':'Arithmetic mean of the three fixed seeds'},np.mean([predictions[n] for n in members],0),None)
winner=max(results,key=lambda r:r['validation']['macro_AP_skill']);(OUT/'stage1_selection.json').write_text(json.dumps({'selected_record':winner,'candidate_count':len(results),'rule':config['selection'],'stage':'DNA/dense-target/prior-AVI-auxiliary search complete; broad-AVI comparisons and final selection still pending','test_evaluated':False},indent=2)+'\n')
print('STAGE1 COMPLETE',winner['name'],'validation skill',winner['validation']['macro_AP_skill'],flush=True)
