"""Stream fixed biological encoders into small, truly trained multi-study heads.

No v3 test evaluation. Reference chromatin is a physical covariate, not an MPRA
target. The auxiliary AVI head never receives AVI as an input.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
from pathlib import Path
import sys,csv,json,time,argparse
import numpy as np
from common import *
assert_no_reserved_evaluation()
parser=argparse.ArgumentParser();parser.add_argument('--representation',choices=['motif','motif_chromatin','hyena_chromatin'],required=True);parser.add_argument('--teacher',choices=['prior','full'],default='prior');parser.add_argument('--avi-input',action='store_true');parser.add_argument('--tasks',default='activity_effect,activity_avi');parser.add_argument('--validation-chromosomes',default='4,5');parser.add_argument('--seeds',default='42,73,101');args=parser.parse_args()
tasks={'assay_only':[0,0,0],'activity':[.2,0,0],'activity_effect':[.2,.1,0],'activity_avi':[.2,0,.05]};selected_tasks=args.tasks.split(',');assert all(t in tasks for t in selected_tasks)
if args.avi_input:assert args.teacher=='full' and all(tasks[t][2]==0 for t in selected_tasks),'AVI input cannot be echoed into an auxiliary AVI target'
data,original_features,manifest,tr,va,te=load();y=data['label'];study=data['study'];nsal=int(np.sum(study==0));MODELS.mkdir(exist_ok=True)
val_chroms=[int(c) for c in args.validation_chromosomes.split(',')];seeds=[int(s) for s in args.seeds.split(',')];tr,va=development_split(data,val_chroms)
chrom_receipt=json.loads((OUT/'reference_chromatin_feature_receipt.json').read_text());cp=WORK/'reference_chromatin_features_local.npz';assert sha(cp)==chrom_receipt['feature_sha256']
with np.load(cp) as z:chromatin=z['features']
assert len(chromatin)==len(y)
motifs=[]
if args.representation.startswith('motif'):
    receipt=json.loads((OUT/'motif_feature_receipt.json').read_text())
    for stem in ['salomon','siraj']:
        r=next(r for r in receipt['datasets'] if r['study']==stem);assert r['input_sha256']==sha(WORK/(stem+'_sequences_local.npz')) and r['features_sha256']==sha(r['path']);motifs.append(np.load(r['path'],mmap_mode='r'))
    D=motifs[0].shape[1];assert motifs[1].shape[1]==D
else:D=512
if args.representation.endswith('_chromatin'):D+=chromatin.shape[1]
if args.teacher=='prior':teacher,teacher_sources=prior_avi(manifest);avi_matrix=None
else:
    path=OUT/'avi_exact_allele_scores.csv';audit=json.loads((OUT/'avi_retrieval_audit.json').read_text());assert not audit['failed_groups'] and audit['output_sha256']==sha(path) and audit['manifest_sha256']==sha(WORK/'avi_query_pairs.csv')
    aa={r['variant_GRCh38']:r for r in csv.DictReader(path.open())};values=[]
    for row in manifest:
        r=aa[row['variant_GRCh38']];assert r['status']!='query_failed';values.append([float(r['avi_raw']),float(r['avi_phred']),1.] if r['status']=='exact_allele_match' else [np.nan,np.nan,0.])
    avi_matrix=np.asarray(values,dtype='float32');teacher=avi_matrix[:,0];teacher_sources=[{'file':str(path),'sha256':sha(path)}]
if args.avi_input:D+=3
def raw(indices):
    indices=np.asarray(indices);parts=[]
    if motifs:
        x=np.empty((len(indices),motifs[0].shape[1]),dtype='float32');lo=indices<nsal;x[lo]=motifs[0][indices[lo]];x[~lo]=motifs[1][indices[~lo]-nsal];parts.append(x)
    else:parts.append(original_features['hyena_original'][indices])
    if args.representation.endswith('_chromatin'):parts.append(chromatin[indices])
    if args.avi_input:parts.append(avi_matrix[indices])
    return np.concatenate(parts,1).astype('float32')
# Streaming moments avoid materializing a >1 GB feature matrix or using test data.
sums=np.zeros(D,dtype='float64');squares=np.zeros(D,dtype='float64');counts=np.zeros(D,dtype='int64')
for start in range(0,len(tr),1024):
    x=raw(tr[start:start+1024]);ok=np.isfinite(x);safe=np.where(ok,x,0).astype('float64');sums+=safe.sum(0);squares+=(safe*safe).sum(0);counts+=ok.sum(0)
mu=np.divide(sums,counts,out=np.zeros(D),where=counts>0);variance=np.divide(squares,counts,out=np.zeros(D),where=counts>0)-mu*mu;sd=np.sqrt(np.maximum(variance,0));sd[sd<1e-5]=1.;mu=mu.astype('float32');sd=sd.astype('float32')
def standardized(indices):
    x=raw(indices);return np.clip((np.where(np.isfinite(x),x,mu)-mu)/sd,-20,20).astype('float32')
known=np.isfinite(teacher);tmean=float(np.nanmean(teacher[tr]));tsd=float(np.nanstd(teacher[tr]));assert tsd>0;tt=np.clip((np.where(known,teacher,tmean)-tmean)/tsd,-5,5).astype('float32')
targets=np.c_[data['activity_mean'],np.log1p(data['absolute_effect'])].astype('float32');tm=np.stack([targets[tr[study[tr]==i]].mean(0) for i in range(2)]);tstd=np.stack([targets[tr[study[tr]==i]].std(0) for i in range(2)]);tstd[tstd<1e-5]=1.;targets=np.clip((targets-tm[study])/tstd[study],-5,5).astype('float32')
prefix=args.representation+'_'+args.teacher+('_avi_input' if args.avi_input else '')+('' if val_chroms==[4,5] else '_val'+'-'.join(map(str,val_chroms)));config={'representation':args.representation,'feature_dim':D,'input_availability_counts_training':counts.tolist(),'teacher':args.teacher,'AVI_as_input':args.avi_input,'teacher_sources':teacher_sources,'tasks':{t:tasks[t] for t in selected_tasks},'seeds':seeds,'training_chromosomes':sorted(set(data['chromosome'][tr].tolist())),'validation_chromosomes':val_chroms,'test_chromosomes_reserved':[1,2,3],'widths':[128,32],'dropout':.1,'max_epochs':60,'patience':8,'batch_size':512,'lr':.0003,'weight_decay':.01,'selection':'Macro per-study normalized validation AP skill; test not evaluated','input_statistics':'Streaming training-only means and standard deviations; unavailable values imputed to training mean, all-unavailable columns to zero. Peak availability indicators distinguish absent peaks and unmapped coordinates.','encoder_frozen':True,'motif_feature_precision':'float16 storage, float32 training batches','test_evaluated':False}
(OUT/(prefix+'_locked_config.json')).write_text(json.dumps(config,indent=2)+'\n');print('BIOLOGICAL HEAD CONFIG',prefix,'features',D,flush=True)
sys.path.insert(0,str(OLD/'runtime'));import mlx.core as mx;import mlx.optimizers as optim
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
pw=np.array([np.sqrt((sum(study[tr]==i)-y[tr[study[tr]==i]].sum())/y[tr[study[tr]==i]].sum()) for i in range(2)],dtype='float32');sw=np.array([len(tr)/(2*sum(study[tr]==i)) for i in range(2)],dtype='float32');pw_m=mx.array(pw);sw_m=mx.array(sw)
def gelu(x):return .5*x*(1+mx.tanh(.7978845608028654*(x+.044715*x**3)))
def huber(x):return mx.where(mx.abs(x)<=1,.5*x*x,mx.abs(x)-.5)
def forward(p,x,sid,m1=None,m2=None):
    h=gelu(x@p['shared.w1']+p['shared.b1']);h=h if m1 is None else h*m1;h=gelu(h@p['shared.w2']+p['shared.b2']);h=h if m2 is None else h*m2;row=mx.arange(x.shape[0]);return [(h@p[k+'.w']+p[k+'.b'])[row,sid] for k in ['class','activity','effect']]+[h@p['avi.w']+p['avi.b']]
def sigmoid(x):return 1/(1+np.exp(-np.clip(np.asarray(x,dtype='float64'),-30,30)))
def predict(p):
    values=[]
    for a in range(0,len(va),512):
        ix=va[a:a+512];values.extend(sigmoid(np.array(forward(p,mx.array(standardized(ix)),mx.array(study[ix]))[0])))
    return np.array(values)
results=[];predictions={};history=[]
def record(name,kind,cfg,pred,artifact,extra=None):
    row={'name':name,'family':kind,'features':args.representation,'config':cfg,'validation':study_metrics(y[va],pred,study[va]),'artifact':artifact}
    if extra:row.update(extra)
    results.append(row);predictions[name]=pred;(OUT/(prefix+'_search_results.json')).write_text(json.dumps(results,indent=2)+'\n');np.savez_compressed(OUT/(prefix+'_validation_predictions.npz'),row_index=va,label=y[va],study=study[va],**predictions);print('BIO VALIDATION',name,json.dumps(row['validation']),flush=True)
for task in selected_tasks:
    weights=tasks[task];members=[]
    for seed in seeds:
        name=prefix+'_'+task+f'_s{seed}';members.append(name);rng=np.random.default_rng(seed+100);r1=np.random.default_rng(seed+1);r2=np.random.default_rng(seed+2);r3=np.random.default_rng(seed+3)
        params={'shared.w1':mx.array(r1.normal(0,1/np.sqrt(D),(D,128)).astype('float32')),'shared.b1':mx.zeros(128),'shared.w2':mx.array(r2.normal(0,1/np.sqrt(128),(128,32)).astype('float32')),'shared.b2':mx.zeros(32)}
        for k in ['class','activity','effect']:
            params[k+'.w']=mx.array(r3.normal(0,.01,(32,2)).astype('float32'));params[k+'.b']=mx.array(np.array([np.log(y[tr[study[tr]==i]].mean()/(1-y[tr[study[tr]==i]].mean())) for i in range(2)],dtype='float32')) if k=='class' else mx.zeros(2)
        params['avi.w']=mx.zeros(32);params['avi.b']=mx.array(0.);initial={k:np.array(v) for k,v in params.items()};opt=optim.AdamW(.0003,weight_decay=.01)
        def loss(p,x,sid,yy,truth,at,mask,m1,m2):
            z,a,e,v=forward(p,x,sid,m1,m2);w=sw_m[sid];main=mx.mean(w*((1-yy)*mx.logaddexp(0,z)+pw_m[sid]*yy*mx.logaddexp(0,-z)))
            return main+weights[0]*mx.mean(w*huber(a-truth[:,0]))+weights[1]*mx.mean(w*huber(e-truth[:,1]))+weights[2]*mx.sum(w*mask*huber(v-at))/mx.maximum(mx.sum(w*mask),1.)
        grad=mx.compile(mx.value_and_grad(loss));best=-np.inf;stale=0;started=time.time()
        for epoch in range(1,61):
            order=rng.permutation(tr);ll=[]
            for start in range(0,len(order),512):
                ix=order[start:start+512];m1=rng.binomial(1,.9,(len(ix),128)).astype('float32')/.9;m2=rng.binomial(1,.9,(len(ix),32)).astype('float32')/.9
                value,g=grad(params,mx.array(standardized(ix)),mx.array(study[ix]),mx.array(y[ix]),mx.array(targets[ix]),mx.array(tt[ix]),mx.array(known[ix].astype('float32')),mx.array(m1),mx.array(m2));norm=mx.sqrt(sum(mx.sum(v*v) for v in g.values()));scale=mx.minimum(1.,5/(norm+1e-8));params=opt.apply_gradients({k:v*scale for k,v in g.items()},params);mx.eval(params,opt.state);ll.append(float(value));assert np.isfinite(ll[-1])
            pred=predict(params);m=study_metrics(y[va],pred,study[va]);history.append({'name':name,'epoch':epoch,'joint_training_loss':float(np.mean(ll)),'validation':m});(OUT/(prefix+'_training_history.json')).write_text(json.dumps(history,indent=2)+'\n')
            if m['macro_AP_skill']>best+1e-6:best=m['macro_AP_skill'];stale=0;best_pred=pred.copy();best_epoch=epoch;saved={k:np.array(v) for k,v in params.items()}
            else:stale+=1
            if epoch%5==0:print('Bio head',name,epoch,'macro validation skill',m['macro_AP_skill'],'best',best,'seconds',round(time.time()-started,1),flush=True)
            if stale>=8:break
        path=MODELS/(name+'.npz');np.savez_compressed(path,**saved,feature_mean=mu,feature_sd=sd,target_mean=tm,target_sd=tstd,teacher_mean=tmean,teacher_sd=tsd);changes={k:float(np.linalg.norm(saved[k]-initial[k])) for k in params};assert changes['shared.w1']>0 and changes['class.w']>0
        record(name,'biological_feature_head',{'task':task,'seed':seed,'best_epoch':best_epoch,'locked_config':prefix+'_locked_config.json'},best_pred,str(path.relative_to(OUT)),{'parameter_change_norms':changes})
    record(prefix+'_ensemble_'+task,'biological_seed_ensemble',{'members':members,'rule':'Arithmetic mean of three fixed seeds'},np.mean([predictions[n] for n in members],0),None)
winner=max(results,key=lambda r:r['validation']['macro_AP_skill']);(OUT/(prefix+'_selection.json')).write_text(json.dumps({'selected_record':winner,'candidate_count':len(results),'test_evaluated':False,'stage':'Additional biological-feature validation search, final combined selection pending'},indent=2)+'\n');print('BIOLOGICAL SEARCH COMPLETE',prefix,winner['name'],flush=True)
