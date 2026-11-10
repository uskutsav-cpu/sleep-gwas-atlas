"""Train a variant-focused mixture using observed activity only as a target.

Validation only. The reserved chromosomes 1–3 are never evaluated here.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
import sys,json,csv,time,argparse,subprocess
import numpy as np
from common import *
assert_no_reserved_evaluation()
from contrast_features import ContrastFeatures
from contrast_network import mx,initialize,forward,objective
import mlx.optimizers as optim

parser=argparse.ArgumentParser();parser.add_argument('--teacher',choices=['prior','full'],default='prior');parser.add_argument('--validation-chromosomes',default='4,5');parser.add_argument('--seeds',default='42,73,101');parser.add_argument('--tasks',default='auto');parser.add_argument('--wait-for-full-runner',type=int);args=parser.parse_args()
if args.wait_for_full_runner:
    print('Waiting for verified full-AVI runner',args.wait_for_full_runner,flush=True)
    while True:
        state=subprocess.run(['ps','-p',str(args.wait_for_full_runner),'-o','command='],capture_output=True,text=True)
        if state.returncode or 'run_full_avi_searches.py' not in state.stdout:break
        time.sleep(10)
    assert (OUT/'all_scheduled_searches_finished.json').is_file()
assert (OUT/'contrast_implementation_check.json').is_file()
features=ContrastFeatures();data=features.data;y=data['label'];study=data['study'];active=features.active;val_chroms=[int(c) for c in args.validation_chromosomes.split(',')];seeds=[int(s) for s in args.seeds.split(',')];assert all(4<=c<=22 for c in val_chroms)
tr,va=development_split(data,val_chroms)
features.fit(tr);MODELS.mkdir(exist_ok=True);prefix='contrast_'+args.teacher+'_val'+'-'.join(map(str,val_chroms))
if args.teacher=='prior':teacher,sources=prior_avi(features.manifest)
else:
    path=OUT/'avi_exact_allele_scores.csv';audit=json.loads((OUT/'avi_retrieval_audit.json').read_text());assert not audit['failed_groups'] and audit['output_sha256']==sha(path) and audit['manifest_sha256']==sha(WORK/'avi_query_pairs.csv');lookup={r['variant_GRCh38']:r for r in csv.DictReader(path.open())};teacher=[]
    for r in features.manifest:
        row=lookup[r['variant_GRCh38']];assert row['status']!='query_failed';teacher.append(float(row['avi_raw']) if row['status']=='exact_allele_match' else np.nan)
    teacher=np.asarray(teacher,dtype='float32');sources=[{'file':str(path),'sha256':sha(path)}]
known=np.isfinite(teacher);tmean=float(np.nanmean(teacher[tr]));tsd=float(np.nanstd(teacher[tr]));assert tsd>0;tt=np.clip((np.where(known,teacher,tmean)-tmean)/tsd,-5,5).astype('float32')
count=np.zeros((2,2),dtype='int64');pos=count.copy()
for s in range(2):
    for g in range(2):
        ix=tr[(study[tr]==s)&(active[tr]==g)];count[s,g]=len(ix);pos[s,g]=int(y[ix].sum())
assert np.all(count>0) and pos[1,0]==0
prevalence=np.divide(pos,count).astype('float32');gate_prev=(count[:,1]/count.sum(1)).astype('float32');pw=np.sqrt(np.divide(count-pos,pos,out=np.ones((2,2),dtype='float64'),where=pos>0)).astype('float32');pw=np.maximum(pw,1.)
gw=np.zeros((2,2),dtype='float32');gw[0]=len(tr)/(4*count[0]);gw[1,1]=len(tr)/(2*count[1,1]);sw=(len(tr)/(2*count.sum(1))).astype('float32')
activity_mean=np.array([data['activity_mean'][tr[study[tr]==s]].mean() for s in range(2)],dtype='float32');activity_sd=np.array([data['activity_mean'][tr[study[tr]==s]].std() for s in range(2)],dtype='float32');activity_sd[activity_sd<1e-5]=1.;act_target=np.clip((data['activity_mean']-activity_mean[study])/activity_sd[study],-5,5).astype('float32')
abs_value=np.log1p(data['absolute_effect']);effect_scale=np.array([abs_value[tr[study[tr]==s]].std() for s in range(2)],dtype='float32');effect_scale[effect_scale<.01]=.01;effect_target=np.clip(abs_value/effect_scale[study],0,5).astype('float32')
constants={'pw':mx.array(pw),'log_pw':mx.array(np.log(pw)),'group_weights':mx.array(gw),'study_weights':mx.array(sw)}
task_options={'effect':0.,'effect_avi':.05};task_names=(['effect_avi'] if args.teacher=='full' and val_chroms==[4,5] else list(task_options)) if args.tasks=='auto' else args.tasks.split(',');assert all(t in task_options for t in task_names)
config={'architecture':'Context gate + bias-free variant trunk with multiplicative context modulation, study-specific active/inactive experts','context_dim':features.dc,'variant_dim':features.dv,'variant_normalization':'training RMS only, no mean subtraction; zero contrast remains zero','training_chromosomes':sorted(set(data['chromosome'][tr].tolist())),'validation_chromosomes':val_chroms,'test_chromosomes_reserved':[1,2,3],'seeds':seeds,'tasks':{t:task_options[t] for t in task_names},'fixed_loss_weights':{'overall_call':1.,'conditional_call':1.,'activity_gate':.25,'activity_regression':.1,'absolute_effect_regression':.1},'conditional_counts':count.tolist(),'conditional_positive_calls':pos.tolist(),'conditional_positive_weight':pw.tolist(),'conditional_group_weight':gw.tolist(),'teacher':args.teacher,'teacher_sources':sources,'teacher_available_train':{STUDIES[s]:int(known[tr[study[tr]==s]].sum()) for s in range(2)},'lr':.0003,'weight_decay':.01,'batch':512,'max_epochs':80,'patience':10,'dropout':.1,'selection':'Mean of macro all-call normalized AP skill and macro conditional-active normalized AP skill; all components retained','activity_flag_is_inference_input':False,'calibration_qualification':'Positive-loss log-weight correction is approximate; scores are not clinically calibrated probabilities','test_evaluated':False}
(OUT/(prefix+'_locked_config.json')).write_text(json.dumps(config,indent=2)+'\n');print('LOCKED CONTRAST CONFIG',json.dumps(config),flush=True)
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
def predict(p):
    vals={k:[] for k in ['overall','predicted_active','conditional_scores','absolute_effect']}
    for start in range(0,len(va),512):
        ix=va[start:start+512];c,v=features.transform(ix);out=forward(p,mx.array(c),mx.array(v),mx.array(study[ix]),constants['log_pw'])
        for k in vals:vals[k].append(np.array(out[k]))
    return {k:np.concatenate(v) for k,v in vals.items()}
def evaluate(pred):
    all_metrics=study_metrics(y[va],pred['overall'],study[va]);conditional={};mixture_active={}
    for s,name in enumerate(STUDIES):
        mask=(study[va]==s)&(active[va]==1);assert len(np.unique(y[va[mask]]))==2
        conditional[name]=metrics(y[va[mask]],pred['conditional_scores'][mask,1]);mixture_active[name]=metrics(y[va[mask]],pred['overall'][mask])
    cskill=float(np.mean([m['AP_skill_over_prevalence'] for m in conditional.values()]));return {'all_calls':all_metrics,'conditional_active':conditional,'mixture_within_active':mixture_active,'macro_conditional_AP_skill':cskill,'joint_selection_skill':.5*(all_metrics['macro_AP_skill']+cskill)}
results=[];overall_predictions={};conditional_predictions={};history=[]
def record(name,family,cfg,pred,artifact,extra=None):
    evaluation=evaluate(pred);row={'name':name,'family':family,'features':'context_variant_biological','config':cfg,'validation':evaluation['all_calls'],'variant_validation':evaluation,'artifact':artifact}
    if extra:row.update(extra)
    results.append(row);overall_predictions[name]=pred['overall'];conditional_predictions[name]=pred['conditional_scores'][:,1];(OUT/(prefix+'_search_results.json')).write_text(json.dumps(results,indent=2)+'\n');np.savez_compressed(OUT/(prefix+'_validation_predictions.npz'),row_index=va,label=y[va],study=study[va],**overall_predictions);np.savez_compressed(OUT/(prefix+'_conditional_validation_predictions.npz'),row_index=va,label=y[va],study=study[va],active=active[va],**conditional_predictions);print('CONTRAST VALIDATION',name,json.dumps(evaluation),flush=True)
for task,avi_weight in config['tasks'].items():
    members=[];member_predictions=[]
    for seed in seeds:
        name=prefix+'_'+task+f'_s{seed}';members.append(name);initial=initialize(features.dc,features.dv,seed,gate_prev,prevalence,pw);p={k:mx.array(v) for k,v in initial.items()};rng=np.random.default_rng(seed+500);opt=optim.AdamW(.0003,weight_decay=.01);best=-np.inf;stale=0;started=time.time()
        gradient=mx.compile(mx.value_and_grad(lambda pp,c,v,sid,yy,aa,at,ef,teacher,mask,masks:objective(pp,c,v,sid,yy,aa,at,ef,teacher,mask,masks,constants,avi_weight)))
        for epoch in range(1,81):
            order=rng.permutation(tr);losses=[]
            for start in range(0,len(order),512):
                ix=order[start:start+512];c,v=features.transform(ix);masks=[mx.array(rng.binomial(1,.9,(len(ix),d)).astype('float32')/.9) for d in [128,32,128,32]]
                loss,g=gradient(p,mx.array(c),mx.array(v),mx.array(study[ix]),mx.array(y[ix]),mx.array(active[ix]),mx.array(act_target[ix]),mx.array(effect_target[ix]),mx.array(tt[ix]),mx.array(known[ix].astype('float32')),masks);norm=mx.sqrt(sum(mx.sum(x*x) for x in g.values()));scale=mx.minimum(1.,5/(norm+1e-8));p=opt.apply_gradients({k:x*scale for k,x in g.items()},p);mx.eval(p,opt.state);losses.append(float(loss));assert np.isfinite(losses[-1])
            pred=predict(p);evaluation=evaluate(pred);history.append({'name':name,'epoch':epoch,'joint_training_loss':float(np.mean(losses)),'validation':evaluation});(OUT/(prefix+'_training_history.json')).write_text(json.dumps(history,indent=2)+'\n')
            if evaluation['joint_selection_skill']>best+1e-6:best=evaluation['joint_selection_skill'];stale=0;best_epoch=epoch;best_pred={k:v.copy() for k,v in pred.items()};saved={k:np.array(v) for k,v in p.items()}
            else:stale+=1
            if epoch%5==0:print('Contrast',name,epoch,'joint validation skill',evaluation['joint_selection_skill'],'best',best,'seconds',round(time.time()-started,1),flush=True)
            if stale>=10:break
        path=MODELS/(name+'.npz');np.savez_compressed(path,**saved,context_mean=features.mean,context_sd=features.sd,variant_scale=features.vscale,log_positive_weight=np.log(pw),activity_mean=activity_mean,activity_sd=activity_sd,effect_scale=effect_scale,teacher_mean=tmean,teacher_sd=tsd)
        changes={k:float(np.linalg.norm(saved[k]-initial[k])) for k in p};assert changes['variant.w1']>0 and changes['context.w1']>0 and changes['conditional.w']>0
        record(name,'context_variant_mixture',{'task':task,'seed':seed,'best_epoch':best_epoch,'locked_config':prefix+'_locked_config.json'},best_pred,str(path.relative_to(OUT)),{'parameter_change_norms':changes});member_predictions.append(best_pred)
    averaged={k:np.mean([r[k] for r in member_predictions],axis=0) for k in member_predictions[0]};record(prefix+'_ensemble_'+task,'context_variant_ensemble',{'members':members,'rule':'Mean outputs across the fixed seeds'},averaged,None)
winner=max(results,key=lambda r:r['variant_validation']['joint_selection_skill']);(OUT/(prefix+'_selection.json')).write_text(json.dumps({'selected_record':winner,'candidate_count':len(results),'selection_rule':config['selection'],'test_evaluated':False,'scope':'Variant-focused validation search; final model-family selection and blind evaluation pending'},indent=2)+'\n');print('CONTRAST SEARCH COMPLETE',winner['name'],flush=True)
