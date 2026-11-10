"""Train frozen and last-block-finetuned HyenaDNA with a locked chromosome split."""
from pathlib import Path
import sys,json,csv,time,hashlib,copy
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np
import mlx.core as mx
import mlx.optimizers as optim
import model as reference
import model_mlx as model

ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent
OUT=ROOT/'outputs/sleep_frailty_ai_v1';MODELS=OUT/'models';MODELS.mkdir(exist_ok=True)
SEED=42;BATCH=32
rng=np.random.default_rng(SEED);mx.random.seed(SEED)
mx.metal.set_memory_limit(3*1024**3)
mx.metal.set_cache_limit(256*1024**2)

def metrics(y,p):
    y=np.asarray(y);p=np.asarray(p);order=np.argsort(-p,kind='stable');yy=y[order];ss=p[order]
    ends=np.r_[np.flatnonzero(np.diff(ss)!=0),len(ss)-1];tp=np.cumsum(yy)[ends];fp=ends+1-tp
    pos=float(y.sum());neg=len(y)-pos
    rec=tp/pos;prec=tp/(tp+fp)
    return {'n':len(y),'detected_calls':int(pos),'prevalence':pos/len(y),'auroc':float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/neg])),'average_precision':float(np.sum(np.diff(np.r_[0,rec])*prec)),'brier_uncalibrated':float(np.mean((p-y)**2))}

def threshold_from_validation(y,p):
    order=np.argsort(-p);yy=y[order];ss=p[order];tp=np.cumsum(yy);prec=tp/np.arange(1,len(y)+1);rec=tp/y.sum();f1=2*prec*rec/np.maximum(prec+rec,1e-12)
    ends=np.r_[np.flatnonzero(np.diff(ss)!=0),len(ss)-1];j=ends[np.argmax(f1[ends])]
    return float(ss[j])

def fixed_threshold(y,p,t):
    calls=p>=t;tp=int(np.sum(calls&(y==1)));fp=int(np.sum(calls&(y==0)));fn=int(np.sum(~calls&(y==1)))
    return {'validation_selected_threshold':t,'test_precision':tp/max(tp+fp,1),'test_recall':tp/max(tp+fn,1),'test_true_positives':tp,'test_false_positives':fp}

def sigmoid(z):return 1/(1+np.exp(-np.clip(z,-30,30)))
def snap(p):return {k:mx.array(np.array(v)) for k,v in p.items()}
def save(name,w,h,mu,sd):
    arrays={**{f'encoder::{k}':np.array(v) for k,v in w.items() if k!='lm_head.weight'},**{k:np.array(v) for k,v in h.items()},'feature_mean':mu,'feature_sd':sd}
    np.savez_compressed(MODELS/name,**arrays)

data=np.load(WORK/'prepared.npz');Xr=data['ref'].astype('int32');Xa=data['alt'].astype('int32');positions=data['position'];y=data['label'];chrom=data['chromosome']
tr=np.flatnonzero(chrom<=16);va=np.flatnonzero((chrom>=17)&(chrom<=18));te=np.flatnonzero(chrom>=19)
assert set(tr).isdisjoint(va) and set(tr).isdisjoint(te) and set(va).isdisjoint(te)
npw=reference.read_weights();w={k:mx.array(v) for k,v in npw.items()}
feature_fn=mx.compile(lambda ref,alt,pos:model.strand_features(ref,alt,pos,w))
positive_weight=float(np.sqrt((len(tr)-y[tr].sum())/y[tr].sum()))
config={'seed':SEED,'batch_size':BATCH,'head_width':32,'feature_dim':512,'positive_loss_weight':positive_weight,'head_warmup_lr':.001,'head_warmup_max_epochs':80,'head_warmup_patience':10,'finetune_head_lr':.0003,'finetune_encoder_lr':.00001,'finetune_max_epochs':8,'finetune_patience':2,'fine_tuned_parameters':'Last encoder block and final layer normalization; first block and embedding frozen','checkpoint_selection':'validation average precision only','test_chromosomes':[19,20,21,22],'checkpoint_sha256':hashlib.sha256((WORK/'hyenadna/model.safetensors').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((ROOT/'work/salomon_neuronal_variant_readout.tsv.gz').read_bytes()).hexdigest(),'model_revision':'e8c1effa8673814e257e627d2e1eda9ea5a373f6','backend':'MLX 0.22.1 GPU with complex FFT; NumPy/Autograd numerical reference'}
config['prepared_input_sha256']=hashlib.sha256((WORK/'prepared.npz').read_bytes()).hexdigest()
(OUT/'training_config.json').write_text(json.dumps(config,indent=2)+'\n')
print('Locked configuration',json.dumps(config),flush=True)

def features(indices,fn=feature_fn,tag='features'):
    cache=WORK/(tag.replace(' ','_')+'_features.npy')
    metadata=cache.with_suffix('.cache.json')
    key={'checkpoint_sha256':config['checkpoint_sha256'],'prepared_input_sha256':config['prepared_input_sha256'],'model_code_sha256':hashlib.sha256((WORK/'model_mlx.py').read_bytes()).hexdigest(),'index_sha256':hashlib.sha256(np.ascontiguousarray(indices).tobytes()).hexdigest()}
    if cache.exists() and metadata.exists():
        try:
            saved=json.loads(metadata.read_text())
            if saved['key']==key and saved['cache_sha256']==hashlib.sha256(cache.read_bytes()).hexdigest():
                existing=np.load(cache)
                if existing.shape==(len(indices),512) and np.isfinite(existing).all():return existing
        except (ValueError,KeyError,OSError):pass
    result=[];start=time.time()
    for step,a in enumerate(range(0,len(indices),64)):
        ix=indices[a:a+64];f=fn(mx.array(Xr[ix]),mx.array(Xa[ix]),mx.array(positions[ix]));result.append(np.array(f))
        if step%50==0:print(tag,a,len(indices),'elapsed',round(time.time()-start,1),flush=True)
    output=np.concatenate(result).astype('float32')
    np.save(cache,output)
    metadata.write_text(json.dumps({'key':key,'cache_sha256':hashlib.sha256(cache.read_bytes()).hexdigest()},indent=2)+'\n')
    mx.metal.clear_cache()
    return output

# Test sequences are not embedded/evaluated until checkpoints have been selected.
train_f=features(tr,tag='pretrained train');val_f=features(va,tag='pretrained validation')
mu=train_f.mean(0);sd=train_f.std(0);sd[sd<1e-5]=1
train_s=np.clip((train_f-mu)/sd,-20,20);val_s=np.clip((val_f-mu)/sd,-20,20)
mu_m=mx.array(mu);sd_m=mx.array(sd)

# Fair baseline: both local allele differences and wider sequence-context counts.
def kmer_features(ix):
    out=[]
    for j in ix:
        vectors=[]
        for arr in [Xr[j]-7,Xa[j]-7]:
            c=[]
            for k in [1,2,3]:
                ids=sum(arr[o:len(arr)-k+1+o].astype('int32')*(4**(k-o-1)) for o in range(k));c.extend(np.bincount(ids,minlength=4**k))
            vectors.append(np.array(c,dtype='float32'))
        out.append(np.r_[(vectors[0]+vectors[1])/2,np.abs(vectors[0]-vectors[1])])
    return np.stack(out)

kt=kmer_features(tr);kv=kmer_features(va);km=kt.mean(0);ks=kt.std(0);ks[ks<1e-5]=1
kt=np.clip((kt-km)/ks,-20,20);kv=np.clip((kv-km)/ks,-20,20)
baseline={'weight':mx.zeros(kt.shape[1]),'bias':mx.array(float(np.log(y[tr].mean()/(1-y[tr].mean()))))};bo=optim.Adam(.01)
bf=lambda p,x,yy:model.bce(x@p['weight']+p['bias'],yy,positive_weight)+.01*mx.sum(p['weight']**2)
bg=mx.compile(mx.value_and_grad(bf));best=-1;baseline_history=[]
for epoch in range(1,201):
    loss,g=bg(baseline,mx.array(kt),mx.array(y[tr]));baseline=bo.apply_gradients(g,baseline);mx.eval(baseline,bo.state)
    if epoch%5==0:
        vp=sigmoid(np.array(mx.array(kv)@baseline['weight']+baseline['bias']));m=metrics(y[va],vp);baseline_history.append({'epoch':epoch,'validation_average_precision':m['average_precision']})
        if m['average_precision']>best:best=m['average_precision'];best_baseline=snap(baseline);best_baseline_val=vp.copy()
np.savez_compressed(MODELS/'kmer_baseline.npz',weight=np.array(best_baseline['weight']),bias=np.array(best_baseline['bias']),mean=km,sd=ks)
(OUT/'baseline_training_history.json').write_text(json.dumps(baseline_history,indent=2)+'\n')
print('Kmer baseline validation',best,flush=True)

head={k:mx.array(v) for k,v in reference.init_head(512,float(y[tr].mean())).items()};ho=optim.AdamW(.001,weight_decay=.001)
hf=lambda p,x,yy:model.bce(model.head(x,p),yy,positive_weight)
hg=mx.compile(mx.value_and_grad(hf));history=[];best=-1;stale=0
for epoch in range(1,81):
    losses=[]
    for start in range(0,len(tr),512):
        if start==0:order=rng.permutation(len(tr))
        ix=order[start:start+512];loss,g=hg(head,mx.array(train_s[ix]),mx.array(y[tr[ix]]));head=ho.apply_gradients(g,head);mx.eval(head,ho.state);losses.append(float(loss))
    vp=sigmoid(np.array(model.head(mx.array(val_s),head)));m=metrics(y[va],vp)
    history.append({'phase':'frozen_head','epoch':epoch,'train_weighted_bce':float(np.mean(losses)),'validation_average_precision':m['average_precision'],'validation_auroc':m['auroc']})
    (OUT/'training_history.json').write_text(json.dumps(history,indent=2)+'\n')
    if m['average_precision']>best+1e-6:best=m['average_precision'];best_head=snap(head);best_frozen_val=vp.copy();best_frozen_epoch=epoch;stale=0
    else:stale+=1
    if epoch%5==0:print('Frozen head',epoch,'validation AP',round(m['average_precision'],4),'best',round(best,4),flush=True)
    if stale>=10:break
save('frozen_hyenadna.npz',w,best_head,mu,sd)
print('Frozen checkpoint selected',best_frozen_epoch,best,flush=True)

enc={k:v for k,v in w.items() if (k.startswith('layers.1.') and 'pos_emb' not in k and 'modulation.deltas' not in k) or k.startswith('ln_f.')}
params={**enc,**best_head};enc_opt=optim.AdamW(.00001,weight_decay=.001);head_opt=optim.AdamW(.0003,weight_decay=.001)
def loss_fn(p,ref,alt,pos,yy):
    feature=model.strand_features(ref,alt,pos,{**w,**p})
    feature=mx.clip((feature-mu_m)/sd_m,-20,20)
    return model.bce(model.head(feature,p),yy,positive_weight)
fg=mx.compile(mx.value_and_grad(loss_fn));best=-1;stale=0;best_finetuned=None
for epoch in range(1,9):
    order=rng.permutation(tr);start_time=time.time();losses=[]
    for step,a in enumerate(range(0,len(order),BATCH)):
        ix=order[a:a+BATCH];loss,g=fg(params,mx.array(Xr[ix]),mx.array(Xa[ix]),mx.array(positions[ix]),mx.array(y[ix]))
        norm=mx.sqrt(sum(mx.sum(v*v) for v in g.values()));factor=mx.minimum(1.,1./(norm+1e-8));g={k:v*factor for k,v in g.items()}
        ep={k:params[k] for k in enc};hp={k:v for k,v in params.items() if k.startswith('head.')}
        ep=enc_opt.apply_gradients({k:g[k] for k in ep},ep);hp=head_opt.apply_gradients({k:g[k] for k in hp},hp);params={**ep,**hp};mx.eval(params,enc_opt.state,head_opt.state)
        losses.append(float(loss))
        if not np.isfinite(losses[-1]):raise ValueError('nonfinite loss')
        if step%100==0:print('Finetune',epoch,a,len(order),'loss',round(losses[-1],4),'elapsed',round(time.time()-start_time,1),flush=True)
    current={**w,**{k:params[k] for k in enc}}
    predict_f=mx.compile(lambda ref,alt,pos:model.head(mx.clip((model.strand_features(ref,alt,pos,current)-mu_m)/sd_m,-20,20),params))
    vp=[]
    for a in range(0,len(va),64):
        ix=va[a:a+64];vp.extend(sigmoid(np.array(predict_f(mx.array(Xr[ix]),mx.array(Xa[ix]),mx.array(positions[ix])))))
    vp=np.array(vp);m=metrics(y[va],vp)
    history.append({'phase':'finetune','epoch':epoch,'train_weighted_bce':float(np.mean(losses)),'validation_average_precision':m['average_precision'],'validation_auroc':m['auroc'],'seconds':time.time()-start_time})
    (OUT/'training_history.json').write_text(json.dumps(history,indent=2)+'\n')
    print('Finetune validation',epoch,m,flush=True)
    if m['average_precision']>best+1e-6:
        best=m['average_precision'];best_finetuned=snap(params);best_finetuned_val=vp.copy();best_finetuned_epoch=epoch;stale=0
        save('finetuned_hyenadna.npz',{**w,**{k:params[k] for k in enc}},{k:v for k,v in params.items() if k.startswith('head.')},mu,sd)
    else:stale+=1
    if stale>=2:break

# Test stage starts only here, after all checkpoint/model choices above.
print('FINAL TEST STAGE: no subsequent training or model redesign',flush=True)
test_f=features(te,tag='pretrained final test');test_s=np.clip((test_f-mu)/sd,-20,20)
k_test=np.clip((kmer_features(te)-km)/ks,-20,20)
predictions={'constant':np.full(len(te),y[tr].mean()),'kmer_logistic':sigmoid(np.array(mx.array(k_test)@best_baseline['weight']+best_baseline['bias'])),'frozen_hyenadna':sigmoid(np.array(model.head(mx.array(test_s),best_head)))}
final_w={**w,**{k:best_finetuned[k] for k in enc}};final_h={k:v for k,v in best_finetuned.items() if k.startswith('head.')}
final_fn=mx.compile(lambda ref,alt,pos:model.head(mx.clip((model.strand_features(ref,alt,pos,final_w)-mu_m)/sd_m,-20,20),final_h))
pp=[]
for a in range(0,len(te),64):
    ix=te[a:a+64];pp.extend(sigmoid(np.array(final_fn(mx.array(Xr[ix]),mx.array(Xa[ix]),mx.array(positions[ix])))))
predictions['finetuned_hyenadna']=np.array(pp)
validation={'constant':np.full(len(va),y[tr].mean()),'kmer_logistic':best_baseline_val,'frozen_hyenadna':best_frozen_val,'finetuned_hyenadna':best_finetuned_val}
results={}
for name,p in predictions.items():
    m=metrics(y[te],p);t=threshold_from_validation(y[va],validation[name]);m.update(fixed_threshold(y[te],p,t));results[name]=m
changes={k:float(np.linalg.norm(np.array(best_finetuned[k])-npw[k])) for k in enc}
results_meta={'models':results,'frozen_epoch':best_frozen_epoch,'finetuned_epoch':best_finetuned_epoch,'fine_tuned_scalar_count':int(sum(v.size for v in enc.values())),'head_scalar_count':int(sum(v.size for v in best_head.values())),'native_pretrained_parameter_change_norms':changes,'test_evaluated_after_selection':True,'limitations':['One selected reporter study, four supervised holdout chromosomes, no separate disease outcome cohort','Nonsignificant calls are imperfect labels; probabilities uncalibrated due class-weighted loss','Reference genome may be seen in unsupervised pretraining; assay sampling used prior sequence predictors','Custom numerical port verified against source-formula reference and gradients, not official PyTorch executable','No independent performance proof for sleep/frailty clinical outcomes or the project variants']}
(OUT/'test_results.json').write_text(json.dumps(results_meta,indent=2)+'\n')
with (OUT/'test_predictions.csv').open('w',newline='') as f:
    fields=['row_index','chromosome','source_detected_call',*predictions];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    for j,i in enumerate(te):writer.writerow({'row_index':int(i),'chromosome':int(chrom[i]),'source_detected_call':int(y[i]),**{k:float(v[j]) for k,v in predictions.items()}})
print(json.dumps(results_meta,indent=2),flush=True)
