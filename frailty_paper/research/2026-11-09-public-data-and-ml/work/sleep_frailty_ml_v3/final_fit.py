"""Refit locked recipes using all development data and fixed epoch counts.

No early stopping, parameter selection, calibration or evaluation on chr1–3.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
import json,pickle,time,shutil
from sklearn.base import clone
from common import *
from inference_records import StoredFeatureInference
from contrast_network import mx,initialize,forward as contrast_forward,objective as contrast_objective
import mlx.optimizers as optim
def gelu(x):return .5*x*(1+mx.tanh(.7978845608028654*(x+.044715*x**3)))
def huber(x):return mx.where(mx.abs(x)<=1,.5*x*x,mx.abs(x)-.5)
def bio_forward(p,x,sid,masks=None):
    h=gelu(x@p['shared.w1']+p['shared.b1']);h=h if masks is None else h*masks[0];h=gelu(h@p['shared.w2']+p['shared.b2']);h=h if masks is None else h*masks[1];row=mx.arange(len(sid));return [(h@p[k+'.w']+p[k+'.b'])[row,sid] for k in ['class','activity','effect']]+[h@p['avi.w']+p['avi.b']]
def main():
    assert_no_reserved_evaluation()
    plan_path=OUT/'final_training_locked_plan.json';plan=json.loads(plan_path.read_text());assert not plan['blind_test_evaluated'];assert plan['protocol_sha256']==sha(OUT/'final_selection_protocol.json') and plan['development_comparison_sha256']==sha(OUT/'development_comparison.json')
    for s in ['salomon','siraj']:assert plan['source_sequence_sha256'][s]==sha(WORK/(s+'_sequences_local.npz'))
    engine=StoredFeatureInference();data=engine.data;features=engine.features;sid=data['study'];y=data['label'];tr=np.flatnonzero(data['chromosome']>=4);assert set(data['chromosome'][tr])==set(plan['training_chromosomes']);active=features.active;avi=engine.avi_features()[:,0];known=np.isfinite(avi);tmean=float(np.nanmean(avi[tr]));tsd=float(np.nanstd(avi[tr]));tt=np.clip((np.where(known,avi,tmean)-tmean)/tsd,-5,5).astype('float32');records={};history=[]
    assert tsd>0;mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
    for recipe in plan['recipes']:
        family=recipe['family_id'];prefix='final_'+family;marker=OUT/(prefix+'_fit_finished.json')
        if marker.exists():
            previous=json.loads(marker.read_text());assert previous['plan_sha256']==sha(plan_path)
            for r in previous['records']:
                if r['artifact']:assert sha(OUT/r['artifact'])==r['checkpoint_sha256']
                records[r['name']]=r
            print('Reuse verified final fit',family,flush=True);continue
        fitted=[]
        if recipe['kind']=='conditional_stack':
            source=recipe['source'];source_path=OUT/source['artifact'];assert sha(source_path)==source['sha256'];base_name='final_'+recipe['base_family'];base_name=base_name if base_name in records else base_name+'_ensemble';assert base_name in records
            target=MODELS/(prefix+'.pkl');shutil.copy2(source_path,target)
            with target.open('rb') as f:models=pickle.load(f)
            fitted.append({'name':prefix,'family':'conditional_avi_stacker','features':'conditional_expert_avi','config':{'base_model':base_name,'include_DNA_expert':recipe['include_DNA_expert'],'training':'All measured-active development rows using held-fold base predictions; no reserved rows','calibration_qualification':'Prediction-time base ensemble is refitted; output is an uncalibrated research ranking, not a clinical probability'},'artifact':str(target.relative_to(OUT)),'checkpoint_sha256':sha(target),'nonzero_trained_coefficient_counts':{STUDIES[s]:int(np.count_nonzero(models[s]['classifier'].coef_)) for s in range(2)},'training_indices_sha256':source['training_indices_sha256']})
        elif recipe['kind']=='classical':
            source=recipe['source_record'];r={**source,'name':prefix,'development_source_name':source['name']};r.pop('validation',None)
            if source['family']=='direct_avi':r['artifact']=None;r['config']={**r['config'],'imputation_raw_training_median':float(np.nanmedian(avi[tr]))}
            else:
                with (OUT/source['artifact']).open('rb') as f:original=pickle.load(f)
                models={}
                for s in range(2):
                    ix=tr[sid[tr]==s];m=clone(original[s]);x=engine.matrix(source['features'],ix);weight=source['config'].get('weighting','none');pw=float(np.sqrt((len(ix)-y[ix].sum())/y[ix].sum()));cw=None if weight=='none' else {0:1.,1:pw}
                    if source['family']=='study_logistic':step=list(m.named_steps)[-1];m.set_params(**{step+'__class_weight':cw});m.fit(x,y[ix])
                    else:m.fit(x,y[ix],**{('histgradientboostingclassifier__sample_weight' if hasattr(m,'named_steps') else 'sample_weight'):None if cw is None else np.where(y[ix]>0,pw,1.)})
                    models[s]=m
                target=MODELS/(prefix+'.pkl')
                with target.open('wb') as f:pickle.dump(models,f)
                r['artifact']=str(target.relative_to(OUT));r['checkpoint_sha256']=sha(target)
            fitted.append(r)
        else:
            contrast=recipe['network_family']=='context_variant_mixture';task=recipe['task'];D=None
            if contrast:
                features.fit(tr);count=np.zeros((2,2),dtype='int64');pos=count.copy()
                for s in range(2):
                    for g in range(2):ix=tr[(sid[tr]==s)&(active[tr]==g)];count[s,g]=len(ix);pos[s,g]=int(y[ix].sum())
                assert pos[1,0]==0;prevalence=(pos/count).astype('float32');gate_prev=(count[:,1]/count.sum(1)).astype('float32');pw=np.maximum(np.sqrt(np.divide(count-pos,pos,out=np.ones((2,2)),where=pos>0)),1.).astype('float32');gw=np.zeros((2,2),dtype='float32');gw[0]=len(tr)/(4*count[0]);gw[1,1]=len(tr)/(2*count[1,1]);sw=(len(tr)/(2*count.sum(1))).astype('float32');constants={'pw':mx.array(pw),'log_pw':mx.array(np.log(pw)),'group_weights':mx.array(gw),'study_weights':mx.array(sw)}
                am=np.array([data['activity_mean'][tr[sid[tr]==s]].mean() for s in range(2)],dtype='float32');asd=np.array([data['activity_mean'][tr[sid[tr]==s]].std() for s in range(2)],dtype='float32');asd[asd<1e-5]=1.;at=np.clip((data['activity_mean']-am[sid])/asd[sid],-5,5).astype('float32');absvalue=np.log1p(data['absolute_effect']);es=np.array([absvalue[tr[sid[tr]==s]].std() for s in range(2)],dtype='float32');es[es<.01]=.01;et=np.clip(absvalue/es[sid],0,5).astype('float32');avi_weight=.05 if task=='effect_avi' else 0.
                transforms={'context_mean':features.mean,'context_sd':features.sd,'variant_scale':features.vscale,'log_positive_weight':np.log(pw),'activity_mean':am,'activity_sd':asd,'effect_scale':es}
            else:
                def raw(ix):return engine.matrix(recipe['features'],ix,recipe['AVI_as_input'])
                D=raw(tr[:1]).shape[1];sums=np.zeros(D);sq=np.zeros(D);count=np.zeros(D,dtype='int64')
                for start in range(0,len(tr),1024):
                    x=raw(tr[start:start+1024]);ok=np.isfinite(x);safe=np.where(ok,x,0).astype('float64');sums+=safe.sum(0);sq+=(safe*safe).sum(0);count+=ok.sum(0)
                mu=np.divide(sums,count,out=np.zeros(D),where=count>0);sd=np.sqrt(np.maximum(np.divide(sq,count,out=np.zeros(D),where=count>0)-mu*mu,0));sd[sd<1e-5]=1.;mu=mu.astype('float32');sd=sd.astype('float32')
                def standard(ix):x=raw(ix);return np.clip((np.where(np.isfinite(x),x,mu)-mu)/sd,-20,20).astype('float32')
                targets=np.c_[data['activity_mean'],np.log1p(data['absolute_effect'])].astype('float32');tm=np.stack([targets[tr[sid[tr]==s]].mean(0) for s in range(2)]);tstd=np.stack([targets[tr[sid[tr]==s]].std(0) for s in range(2)]);tstd[tstd<1e-5]=1.;targets=np.clip((targets-tm[sid])/tstd[sid],-5,5).astype('float32');pw=np.array([np.sqrt((np.sum(sid[tr]==s)-y[tr[sid[tr]==s]].sum())/y[tr[sid[tr]==s]].sum()) for s in range(2)],dtype='float32');sw=np.array([len(tr)/(2*np.sum(sid[tr]==s)) for s in range(2)],dtype='float32');pwm=mx.array(pw);swm=mx.array(sw);weights={'activity_effect':[.2,.1,0.],'activity_avi':[.2,0.,.05]}[task];transforms={'feature_mean':mu,'feature_sd':sd,'target_mean':tm,'target_sd':tstd}
            member_names=[]
            for seed in recipe['seeds']:
                name=prefix+f'_s{seed}';member_names.append(name)
                if contrast:initial=initialize(features.dc,features.dv,seed,gate_prev,prevalence,pw);rng=np.random.default_rng(seed+500)
                else:
                    rng=np.random.default_rng(seed+100);r1=np.random.default_rng(seed+1);r2=np.random.default_rng(seed+2);r3=np.random.default_rng(seed+3);initial={'shared.w1':r1.normal(0,1/np.sqrt(D),(D,128)).astype('float32'),'shared.b1':np.zeros(128,dtype='float32'),'shared.w2':r2.normal(0,1/np.sqrt(128),(128,32)).astype('float32'),'shared.b2':np.zeros(32,dtype='float32')}
                    for k in ['class','activity','effect']:initial[k+'.w']=r3.normal(0,.01,(32,2)).astype('float32');initial[k+'.b']=np.array([np.log(y[tr[sid[tr]==s]].mean()/(1-y[tr[sid[tr]==s]].mean())) for s in range(2)],dtype='float32') if k=='class' else np.zeros(2,dtype='float32')
                    initial['avi.w']=np.zeros(32,dtype='float32');initial['avi.b']=np.array(0.,dtype='float32')
                params={k:mx.array(v) for k,v in initial.items()};opt=optim.AdamW(.0003,weight_decay=.01)
                if contrast:gradient=mx.compile(mx.value_and_grad(lambda p,c,v,s,yy,aa,at0,ef,t,mask,masks:contrast_objective(p,c,v,s,yy,aa,at0,ef,t,mask,masks,constants,avi_weight)))
                else:
                    def loss(p,x,s,yy,truth,teacher,mask,masks):
                        z,a,e,v=bio_forward(p,x,s,masks);w=swm[s];main=mx.mean(w*((1-yy)*mx.logaddexp(0,z)+pwm[s]*yy*mx.logaddexp(0,-z)));return main+weights[0]*mx.mean(w*huber(a-truth[:,0]))+weights[1]*mx.mean(w*huber(e-truth[:,1]))+weights[2]*mx.sum(w*mask*huber(v-teacher))/mx.maximum(mx.sum(w*mask),1.)
                    gradient=mx.compile(mx.value_and_grad(loss))
                for epoch in range(1,recipe['fixed_epochs_by_seed'][str(seed)]+1):
                    order=rng.permutation(tr);losses=[];start_time=time.time()
                    for start in range(0,len(order),512):
                        ix=order[start:start+512];widths=[128,32,128,32] if contrast else [128,32];masks=[mx.array(rng.binomial(1,.9,(len(ix),d)).astype('float32')/.9) for d in widths]
                        if contrast:
                            c,v=features.transform(ix);value,g=gradient(params,mx.array(c),mx.array(v),mx.array(sid[ix]),mx.array(y[ix]),mx.array(active[ix]),mx.array(at[ix]),mx.array(et[ix]),mx.array(tt[ix]),mx.array(known[ix].astype('float32')),masks)
                        else:value,g=gradient(params,mx.array(standard(ix)),mx.array(sid[ix]),mx.array(y[ix]),mx.array(targets[ix]),mx.array(tt[ix]),mx.array(known[ix].astype('float32')),masks)
                        norm=mx.sqrt(sum(mx.sum(v*v) for v in g.values()));scale=mx.minimum(1.,5/(norm+1e-8));params=opt.apply_gradients({k:v*scale for k,v in g.items()},params);mx.eval(params,opt.state);losses.append(float(value));assert np.isfinite(losses[-1])
                    history.append({'model':name,'epoch':epoch,'joint_training_loss':float(np.mean(losses)),'seconds':time.time()-start_time,'validation_or_test_evaluated':False});(OUT/'final_fit_training_history.json').write_text(json.dumps(history,indent=2)+'\n');print('FINAL FIT',name,'fixed epoch',epoch,'loss',history[-1]['joint_training_loss'],flush=True)
                saved={k:np.array(v) for k,v in params.items()};change={k:int(np.count_nonzero(saved[k]!=initial[k])) for k in initial};assert sum(change.values())>0;target=MODELS/(name+'.npz');np.savez_compressed(target,**saved,**transforms,teacher_mean=tmean,teacher_sd=tsd)
                cfg={'task':task,'seed':seed,'fixed_epochs':recipe['fixed_epochs_by_seed'][str(seed)],'locked_config':'final_'+family+'_inference_config.json'};(OUT/cfg['locked_config']).write_text(json.dumps({'AVI_as_input':recipe['AVI_as_input'],'final_plan_sha256':sha(plan_path)},indent=2)+'\n');fitted.append({'name':name,'family':recipe['network_family'],'features':recipe['features'],'config':cfg,'artifact':str(target.relative_to(OUT)),'checkpoint_sha256':sha(target),'changed_scalar_counts':change,'native_source_encoders_frozen':True})
            fitted.append({'name':prefix+'_ensemble','family':'context_variant_ensemble' if contrast else 'biological_seed_ensemble','features':recipe['features'],'config':{'members':member_names,'rule':'Fixed three-seed arithmetic mean'},'artifact':None})
        for r in fitted:records[r['name']]=r
        marker.write_text(json.dumps({'family_id':family,'records':fitted,'plan_sha256':sha(plan_path),'blind_test_evaluated':False},indent=2)+'\n')
    resolved={}
    for role,family in {**plan['roles'],**plan['controls'],**plan['matched_neural_ablations'],**plan['conditional_stacking_controls']}.items():
        candidates=[r for r in records.values() if r['name'] in ['final_'+family,'final_'+family+'_ensemble']];assert len(candidates)==1;resolved[role]=candidates[0]['name']
    registry={'records':records,'role_model_names':resolved,'plan_sha256':sha(plan_path),'training_chromosomes':plan['training_chromosomes'],'training_rows':len(tr),'training_indices_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'training_counts':{STUDIES[s]:{'n':int(np.sum(sid[tr]==s)),'calls':int(y[tr[sid[tr]==s]].sum())} for s in range(2)},'source_code_sha256':sha(Path(__file__)),'blind_test_evaluated':False,'scope':'Final fixed-epoch research fits. Scores are uncalibrated laboratory readout rankings, not disease probabilities.'};(OUT/'final_model_registry.json').write_text(json.dumps(registry,indent=2)+'\n');print('FINAL REFITS COMPLETE; VERIFICATION AND BLIND EVALUATION REMAIN',flush=True)
if __name__=='__main__':main()
