"""Teach the tiny DNA model using permitted AVI predictions only.

Teacher agreement is imitation, not experimental/disease validation. Never
reads the reserved teacher regions or modifies the protected V3 model.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
import sys,time,argparse,subprocess
from common import *
from scipy.stats import spearmanr
from sklearn.metrics import r2_score,roc_auc_score,average_precision_score
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import mlx.core as mx;import mlx.optimizers as optim;import model as original;import model_mlx as dna

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-for-preparation-pid',type=int);a=parser.parse_args()
    if a.wait_for_preparation_pid:
        print('Waiting for verified teacher preparation',a.wait_for_preparation_pid,flush=True)
        while True:
            s=subprocess.run(['ps','-p',str(a.wait_for_preparation_pid),'-o','command='],capture_output=True,text=True)
            if s.returncode or 'sleep_frailty_ml_v4/prepare_teacher_data.py' not in s.stdout:break
            time.sleep(10)
    receipt=json.loads((OUT/'teacher_data_receipt.json').read_text());plan=json.loads((OUT/'training_plan.json').read_text());assert receipt['base_weights_sha256']==sha(OLD/'hyenadna/model.safetensors');datasets=[]
    for row in receipt['datasets']:
        assert sha(row['file'])==row['sha256']
        with np.load(row['file']) as z:d={k:z[k] for k in z.files}
        assert not np.any(d['chromosome']<=3);datasets.append(d)
    feat=np.concatenate([d['features'][d['split']==0] for d in datasets]);mu=feat.mean(0);sd=feat.std(0);sd[sd<1e-5]=1.;raw_targets=np.concatenate([d['target'][d['split']==0] for d in datasets]);transformed=np.c_[np.arcsinh(raw_targets[:,0]),raw_targets[:,1]];tm=transformed.mean(0).astype('float32');ts=transformed.std(0).astype('float32');assert np.all(ts>0);high_threshold=float(np.quantile(raw_targets[:,0],.9));del feat
    for d in datasets:
        d['x']=np.clip((d['features']-mu)/sd,-20,20).astype('float32');d['y']=(np.c_[np.arcsinh(d['target'][:,0]),d['target'][:,1]]-tm)/ts;d['y']=d['y'].astype('float32');d['high']=(d['target'][:,0]>=high_threshold).astype('float32')
    train_counts=[int(np.sum(d['split']==0)) for d in datasets];total=sum(train_counts);positive=sum(int(d['high'][d['split']==0].sum()) for d in datasets);pw=float(np.sqrt((total-positive)/positive));rng=np.random.default_rng(41207);h0={'head.w1':rng.normal(0,1/np.sqrt(512),(512,64)).astype('float32'),'head.b1':np.zeros(64,dtype='float32'),'teacher.w':rng.normal(0,.01,(64,2)).astype('float32'),'teacher.b':np.zeros(2,dtype='float32'),'high.w':rng.normal(0,.01,64).astype('float32'),'high.b':np.array(np.log(positive/(total-positive)),dtype='float32')};head={k:mx.array(v) for k,v in h0.items()};base=original.read_weights();weights={k:mx.array(v) for k,v in base.items()};mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2);models=DATA/'models';models.mkdir(exist_ok=True);path=models/'teacher_student_selected.npz';history=[];best=-np.inf;stale=0;mu_m=mx.array(mu);sd_m=mx.array(sd)
    def forward(x,p):
        h=mx.tanh(x@p['head.w1']+p['head.b1']);return h@p['teacher.w']+p['teacher.b'],h@p['high.w']+p['high.b']
    def loss(x,p,y,high):
        z,logit=forward(x,p);v=z-y;huber=mx.mean(mx.where(mx.abs(v)<=1,.5*v*v,mx.abs(v)-.5));return huber+.2*mx.mean((1-high)*mx.logaddexp(0,logit)+pw*high*mx.logaddexp(0,-logit))
    def evaluate(current_weights,head_params,native):
        preds=[];targets=[];anchors=[];high=[];highpred=[]
        for d in datasets:
            ix=np.flatnonzero(d['split']==1)
            for start in range(0,len(ix),64 if native else 512):
                jj=ix[start:start+(64 if native else 512)]
                if native:x=mx.clip((dna.strand_features(mx.array(d['ref'][jj].astype('int32')),mx.array(d['alt'][jj].astype('int32')),mx.array(d['position'][jj]),current_weights)-mu_m)/sd_m,-20,20)
                else:x=mx.array(d['x'][jj])
                z,logit=forward(x,head_params);preds.append(np.array(z));targets.append(d['target'][jj]);anchors.extend(d['anchor_id'][jj].tolist());high.extend(d['high'][jj].tolist());highpred.extend((1/(1+np.exp(-np.clip(np.array(logit),-30,30)))).tolist())
        p=np.concatenate(preds);truth=np.concatenate(targets);decoded=np.c_[np.sinh(np.clip(p[:,0]*ts[0]+tm[0],-20,20)),p[:,1]*ts[1]+tm[1]];anchors=np.asarray(anchors);group=[]
        for anchor in np.unique(anchors):
            use=anchors==anchor
            if len(np.unique(truth[use,0]))>1 and len(np.unique(decoded[use,0]))>1:group.append(float(spearmanr(truth[use,0],decoded[use,0]).statistic))
        assert group and np.isfinite(group).all();return {'new_teacher_validation_rows':len(truth),'novel_anchor_count':len(group),'mean_within_anchor_Spearman':float(np.mean(group)),'pooled_rawAVI_Spearman':float(spearmanr(truth[:,0],decoded[:,0]).statistic),'pooled_rawAVI_R2':float(r2_score(truth[:,0],decoded[:,0])),'pooled_PHRED_R2':float(r2_score(truth[:,1],decoded[:,1])),'high_AVI_base_rate':float(np.mean(high)),'high_AVI_AUROC':float(roc_auc_score(high,highpred)),'high_AVI_AP':float(average_precision_score(high,highpred))}
    def record(phase,epoch,value,w,h):
        nonlocal best,stale
        metrics=evaluate(w,h,phase=='native');history.append({'phase':phase,'epoch':epoch,'training_loss':value,'validation':metrics});write_json(OUT/'teacher_training_history.json',history);score=metrics['mean_within_anchor_Spearman']
        if score>best+1e-5:
            best=score;stale=0;np.savez_compressed(path,**{'encoder::'+k:np.array(v) for k,v in w.items()},**{k:np.array(v) for k,v in h.items()},feature_mean=mu,feature_sd=sd,target_mean=tm,target_sd=ts,high_threshold=np.array(high_threshold),selected_phase=np.array(phase),selected_epoch=np.array(epoch))
        else:stale+=1
        print('TEACHER VALIDATION',phase,epoch,json.dumps(metrics),flush=True)
    config={'training_plan_sha256':sha(OUT/'training_plan.json'),'data_receipt_sha256':sha(OUT/'teacher_data_receipt.json'),'base_sha256':receipt['base_weights_sha256'],'training_rows':total,'new_validation_rows':sum(int(np.sum(d['split']==1)) for d in datasets),'training_high_AVI_threshold':high_threshold,'target_transform':'asinhAVIraw,PHRED; fitted train means/sds','teacher_values_are_inference_inputs':False,'native_trainable':'Last encoder block excluding positional embeddings/modulation constants; final normalization; prediction head','selected_metric':'Mean Spearman within novel validation anchors','max_warmup_epochs':20,'max_native_epochs':2,'seed':41207,'actual_lab_outcomes_used':False,'reserved_teacher_test_used':False};write_json(OUT/'teacher_locked_config.json',config);print('LOCKED TEACHER TRAINING',json.dumps(config),flush=True)
    def orders(batch):
        tasks=[]
        for k,d in enumerate(datasets):
            ids=rng.permutation(np.flatnonzero(d['split']==0));tasks.extend((k,ids[i:i+batch]) for i in range(0,len(ids),batch))
        rng.shuffle(tasks);return tasks
    opt=optim.AdamW(.0003,weight_decay=.001);grad=mx.compile(mx.value_and_grad(lambda p,x,y,h:loss(x,p,y,h)))
    for epoch in range(1,21):
        values=[]
        for k,ix in orders(512):
            d=datasets[k];value,g=grad(head,mx.array(d['x'][ix]),mx.array(d['y'][ix]),mx.array(d['high'][ix]));head=opt.apply_gradients(g,head);mx.eval(head,opt.state);values.append(float(value));assert np.isfinite(values[-1])
        record('frozen',epoch,float(np.mean(values)),weights,head)
        if stale>=4:break
    with np.load(path) as z:head={k:mx.array(z[k]) for k in z.files if k.startswith(('head.','teacher.','high.'))};warmup_record=next(r for r in history if r['phase']==str(z['selected_phase']) and r['epoch']==int(z['selected_epoch']))
    native_keys=[k for k in base if (k.startswith('layers.1.') and 'pos_emb' not in k and 'modulation.deltas' not in k) or k.startswith('ln_f.')];params={**{k:weights[k] for k in native_keys},**head};eo=optim.AdamW(.00001,weight_decay=.001);ho=optim.AdamW(.0003,weight_decay=.001)
    def native_loss(p,a,b,pos,y,high):
        f=dna.strand_features(a,b,pos,{**weights,**p});x=mx.clip((f-mu_m)/sd_m,-20,20);return loss(x,p,y,high)
    grad=mx.compile(mx.value_and_grad(native_loss));stale=0
    for epoch in range(1,3):
        values=[];started=time.time()
        for step,(k,ix) in enumerate(orders(32)):
            d=datasets[k];value,g=grad(params,mx.array(d['ref'][ix].astype('int32')),mx.array(d['alt'][ix].astype('int32')),mx.array(d['position'][ix]),mx.array(d['y'][ix]),mx.array(d['high'][ix]));norm=mx.sqrt(sum(mx.sum(v*v) for v in g.values()));scale=mx.minimum(1.,1/(norm+1e-8));ep={q:params[q] for q in native_keys};hp={q:params[q] for q in head};ep=eo.apply_gradients({q:g[q]*scale for q in ep},ep);hp=ho.apply_gradients({q:g[q]*scale for q in hp},hp);params={**ep,**hp};mx.eval(params,eo.state,ho.state);values.append(float(value));assert np.isfinite(values[-1])
            if step%200==0:print('NATIVE TEACHER UPDATE',epoch,step,'loss',values[-1],'seconds',round(time.time()-started,1),flush=True)
        current={**weights,**{k:params[k] for k in native_keys}};hh={k:params[k] for k in head};record('native',epoch,float(np.mean(values)),current,hh)
        if stale>=1:break
    changed_last={k:int(np.count_nonzero(np.array(params[k])!=base[k])) for k in native_keys};assert sum(changed_last.values())>0
    with np.load(path) as z:selected_changed={k:int(np.count_nonzero(z['encoder::'+k]!=base[k])) for k in native_keys};phase=str(z['selected_phase']);epoch=int(z['selected_epoch'])
    write_json(OUT/'teacher_training_status.json',{'status':'Teacher-only pilot trained','selected_phase':phase,'selected_epoch':epoch,'best_novel_anchor_Spearman':best,'warmup_validation':warmup_record['validation'],'last_attempt_changed_native_scalars':sum(changed_last.values()),'selected_changed_native_scalars':sum(selected_changed.values()),'checkpoint':str(path),'checkpoint_sha256':sha(path),'training_rows':total,'new_exact_teacher_pairs':receipt['new_exact_teacher_rows'],'scope':'Predicting permitted AlphaGenome AVI outputs. No measured biological improvement or disease prediction established.','next':'Run matched actual-assay development probes and find a new independent measured source before claiming stronger biology. Reserved teacher regions remain untouched.','old_V3_registry_unchanged':sha(V3/'final_model_registry.json')==json.loads((OUT/'teacher_region_plan.json').read_text())['old_V3_registry_sha256']});print('TEACHER PILOT TRAINING COMPLETE; BIOLOGICAL VALIDATION STILL REQUIRED',flush=True)
if __name__=='__main__':main()
