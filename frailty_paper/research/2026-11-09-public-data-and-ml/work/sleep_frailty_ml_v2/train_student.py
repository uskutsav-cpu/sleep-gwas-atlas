"""Conditional native-weight fine-tuning with an auxiliary permitted AVI target.

Wait for the first validation search. No benchmark or external result is read.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
from pathlib import Path
import sys,json,time,subprocess,argparse,hashlib
import numpy as np
from common import *

p=argparse.ArgumentParser();p.add_argument('--wait-for-search-pid',type=int);p.add_argument('--exploratory-auxiliary',action='store_true');args=p.parse_args()
if args.wait_for_search_pid:
    print('Waiting for verified search PID',args.wait_for_search_pid,flush=True)
    while True:
        ps=subprocess.run(['ps','-p',str(args.wait_for_search_pid),'-o','command='],capture_output=True,text=True)
        if ps.returncode or 'sleep_frailty_ml_v2/train_search.py' not in ps.stdout:break
        time.sleep(10)
selection=json.loads((OUT/'model_selection.json').read_text());records=json.loads((OUT/'search_results.json').read_text())
assert len(records)==selection['candidate_count']
best_avi=max(r['validation']['average_precision'] for r in records if r['features'].endswith('_avi') or r['features']=='avi')
best_dna=max(r['validation']['average_precision'] for r in records if not (r['features'].endswith('_avi') or r['features']=='avi'))
if best_avi<=best_dna and not args.exploratory_auxiliary:
    (OUT/'student_status.json').write_text(json.dumps({'status':'not_triggered','reason':'Permitted AVI did not improve validation AP over the best DNA-only candidate','best_avi_validation_AP':best_avi,'best_dna_validation_AP':best_dna},indent=2)+'\n')
    print('Auxiliary student not triggered: AVI did not improve validation',flush=True);sys.exit(0)
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD))
import mlx.core as mx
import mlx.optimizers as optim
import model as reference
import model_mlx as model
mx.metal.set_memory_limit(3*1024**3);mx.metal.set_cache_limit(256*1024**2)
data,seq,avi,manifest,tr,va,te=load();y=data['label'];f=seq['hyena_original'];mu=f[tr].mean(0);sd=f[tr].std(0);sd[sd<1e-5]=1
xt=np.clip((f[tr]-mu)/sd,-20,20);xv=np.clip((f[va]-mu)/sd,-20,20)
mask=np.isfinite(avi[:,0]).astype('float32');median=float(np.nanmedian(avi[tr,0]));teacher_mean=float(np.nanmean(avi[tr,0]));teacher_sd=float(np.nanstd(avi[tr,0]));assert teacher_sd>0
target=np.clip((np.nan_to_num(avi[:,0],nan=teacher_mean)-teacher_mean)/teacher_sd,-5,5)
rng=np.random.default_rng(1207);mx.random.seed(1207);pw=float(np.sqrt((len(tr)-y[tr].sum())/y[tr].sum()))
base=reference.read_weights();weights={k:mx.array(v) for k,v in base.items()};initial_head=np.load(ROOT/'outputs/sleep_frailty_ai_v1/models/frozen_hyenadna.npz')
head={k:mx.array(initial_head[k]) for k in initial_head.files if k.startswith('head.')};head['aux.w']=mx.zeros(512);head['aux.b']=mx.array(0.)
history=[];best=-1;stale=0;path=MODELS/'avi_auxiliary_student.npz';mu_m=mx.array(mu);sd_m=mx.array(sd)
config={'seed':1207,'trigger':'Best AVI-assisted validation AP exceeds best DNA-only validation AP','best_avi_validation_AP':best_avi,'best_dna_validation_AP':best_dna,'target':'Measured NGN2 call plus auxiliary standardized/clipped raw AVI','auxiliary_loss_weight':.05,'avi_input_to_assay_head':False,'teacher_mean_training':teacher_mean,'teacher_sd_training':teacher_sd,'teacher_clip':[-5,5],'warmup_max_epochs':20,'warmup_patience':5,'finetune_max_epochs':4,'finetune_patience':2,'batch_size':32,'head_lr':.0003,'encoder_lr':.00001,'fine_tune':'last native block and final normalization; original first block and nucleotide embedding fixed','checkpoint_selection':'validation assay average precision only','independence':'Learning to approximate AVI is auxiliary supervision, not experimental or disease validation'}
if args.exploratory_auxiliary:
    config.update({'trigger':'Exploratory additional experiment after round1 failed to show AVI-input benefit','post_round1_extension':True,'evaluation_status':'Reused benchmark and the first external study have now been observed; subsequent checks on them are development comparisons, not fresh independent validation','architecture_and_loss':'Same auxiliary architecture, loss and learning rates specified before round1 performance; the decision to execute despite the failed trigger is new and exploratory'})
(OUT/'student_config.json').write_text(json.dumps(config,indent=2)+'\n')
def snap(p):return {k:mx.array(np.array(v)) for k,v in p.items()}
def sigmoid(x):return 1/(1+np.exp(-np.clip(np.asarray(x,dtype='float64'),-30,30)))
def losses(features,p,yy,tt,mm):
    assay=model.bce(model.head(features,p),yy,pw);prediction=features@p['aux.w']+p['aux.b'];delta=prediction-tt
    # Huber loss limits influence of the deliberately clipped teacher target.
    aux=mx.sum(mm*mx.where(mx.abs(delta)<=1,.5*delta**2,mx.abs(delta)-.5))/mx.maximum(mx.sum(mm),1.)
    return assay+.05*aux
def save(w,h):
    np.savez_compressed(path,**{'encoder::'+k:np.array(v) for k,v in w.items() if k!='lm_head.weight'},**{k:np.array(v) for k,v in h.items()},feature_mean=mu,feature_sd=sd)
def record(phase,epoch,value,p,weights_now,head_now):
    global best,stale,best_p,best_phase,best_epoch
    m=metrics(y[va],p);history.append({'phase':phase,'epoch':epoch,'joint_training_loss':value,'validation':m});(OUT/'student_training_history.json').write_text(json.dumps(history,indent=2)+'\n')
    if m['average_precision']>best+1e-6:
        best=m['average_precision'];stale=0;best_p=p.copy();best_phase=phase;best_epoch=epoch;save(weights_now,head_now)
    else:stale+=1
    print('Student validation',phase,epoch,json.dumps(m),flush=True)
opt=optim.AdamW(.0003,weight_decay=.001);grad=mx.compile(mx.value_and_grad(lambda p,x,yy,tt,mm:losses(x,p,yy,tt,mm)))
for epoch in range(1,21):
    order=rng.permutation(len(tr));ll=[]
    for a in range(0,len(order),512):
        j=order[a:a+512];i=tr[j];value,g=grad(head,mx.array(xt[j]),mx.array(y[i]),mx.array(target[i]),mx.array(mask[i]));head=opt.apply_gradients(g,head);mx.eval(head,opt.state);ll.append(float(value))
    vp=sigmoid(np.array(model.head(mx.array(xv),head)));record('frozen_encoder_auxiliary',epoch,float(np.mean(ll)),vp,weights,head)
    if stale>=5:break
# Start native training from the warmup checkpoint selected on assay validation.
warm=np.load(path);head={k:mx.array(warm[k]) for k in warm.files if k.startswith('head.') or k.startswith('aux.')}
enc={k:v for k,v in weights.items() if (k.startswith('layers.1.') and 'pos_emb' not in k and 'modulation.deltas' not in k) or k.startswith('ln_f.')};params={**enc,**head}
eo=optim.AdamW(.00001,weight_decay=.001);ho=optim.AdamW(.0003,weight_decay=.001)
def loss_fn(p,a,b,pos,yy,tt,mm):
    x=model.strand_features(a,b,pos,{**weights,**p});x=mx.clip((x-mu_m)/sd_m,-20,20)
    return losses(x,p,yy,tt,mm)
grad=mx.compile(mx.value_and_grad(loss_fn));stale=0;changed_last={}
for epoch in range(1,5):
    order=rng.permutation(tr);ll=[];start=time.time()
    for step,a in enumerate(range(0,len(order),32)):
        ix=order[a:a+32];value,g=grad(params,mx.array(data['ref'][ix].astype('int32')),mx.array(data['alt'][ix].astype('int32')),mx.array(data['position'][ix]),mx.array(y[ix]),mx.array(target[ix]),mx.array(mask[ix]));norm=mx.sqrt(sum(mx.sum(v*v) for v in g.values()));scale=mx.minimum(1.,1/(norm+1e-8));g={k:v*scale for k,v in g.items()}
        ep={k:params[k] for k in enc};hp={k:params[k] for k in params if k not in enc};ep=eo.apply_gradients({k:g[k] for k in ep},ep);hp=ho.apply_gradients({k:g[k] for k in hp},hp);params={**ep,**hp};mx.eval(params,eo.state,ho.state);ll.append(float(value));assert np.isfinite(ll[-1])
        if step%100==0:print('Student native update',epoch,a,'/',len(order),'loss',ll[-1],'seconds',round(time.time()-start,1),flush=True)
    current={**weights,**{k:params[k] for k in enc}};hh={k:params[k] for k in params if k not in enc};fn=mx.compile(lambda a,b,p:model.head(mx.clip((model.strand_features(a,b,p,current)-mu_m)/sd_m,-20,20),hh));vp=[]
    for a in range(0,len(va),64):
        ix=va[a:a+64];vp.extend(sigmoid(np.array(fn(mx.array(data['ref'][ix].astype('int32')),mx.array(data['alt'][ix].astype('int32')),mx.array(data['position'][ix])))))
    record('native_auxiliary_finetune',epoch,float(np.mean(ll)),np.array(vp),current,hh)
    changed_last={k:float(np.linalg.norm(np.array(params[k])-base[k])) for k in enc};assert sum(v>0 for v in changed_last.values())>5
    if stale>=2:break
saved=np.load(path);selected_changes={k:float(np.linalg.norm(saved['encoder::'+k]-base[k])) for k in enc}
status={'status':'trained','best_phase':best_phase,'best_epoch':best_epoch,'best_validation_AP':best,'selected_checkpoint_native_parameter_change_norms':selected_changes,'last_attempt_native_parameter_change_norms':changed_last,'checkpoint_sha256':sha(path),'native_parameter_scalars_trained':int(sum(v.size for v in enc.values()))}
(OUT/'student_status.json').write_text(json.dumps(status,indent=2)+'\n')
record_new={'name':'avi_auxiliary_student','family':'distilled_hyenadna','features':'dna_with_avi_auxiliary_training','config':config,'validation':metrics(y[va],best_p),'artifact':str(path.relative_to(OUT))};records.append(record_new)
(OUT/'search_results.json').write_text(json.dumps(records,indent=2)+'\n');old=np.load(OUT/'validation_predictions.npz');np.savez_compressed(OUT/'validation_predictions.npz',**{k:old[k] for k in old.files},avi_auxiliary_student=best_p)
(OUT/'model_selection_before_student.json').write_text(json.dumps(selection,indent=2)+'\n')
winner=max(records,key=lambda r:r['validation']['average_precision']);selection.update({'selected_model':winner['name'],'selected_record':winner,'candidate_count':len(records),'student_experiment_completed':True});(OUT/'model_selection.json').write_text(json.dumps(selection,indent=2)+'\n')
print('STUDENT COMPLETE',json.dumps(status),flush=True)
