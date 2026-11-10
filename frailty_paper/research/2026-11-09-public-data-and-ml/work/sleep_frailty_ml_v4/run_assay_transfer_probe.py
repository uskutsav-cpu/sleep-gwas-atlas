"""Matched measured-laboratory development probe after AVI-only teaching.

Both variants receive the same permitted AVI input, original/taught DNA
features and a64d latent projection with matched dimension. The control
projection is untrained. These are reused development data, not fresh tests.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
import sys,argparse,subprocess,time
from common import *
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import mlx.core as mx;import mlx.optimizers as optim;import model,model_mlx
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-for-student-pid',type=int);a=parser.parse_args()
    if a.wait_for_student_pid:
        print('Waiting for verified AVI student',a.wait_for_student_pid,flush=True)
        while True:
            s=subprocess.run(['ps','-p',str(a.wait_for_student_pid),'-o','command='],capture_output=True,text=True)
            if s.returncode or 'sleep_frailty_ml_v4/train_teacher_student.py' not in s.stdout:break
            time.sleep(10)
    status=json.loads((OUT/'teacher_training_status.json').read_text());path=Path(status['checkpoint']);assert sha(path)==status['checkpoint_sha256'];mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2)
    with np.load(path) as z:w={k.removeprefix('encoder::'):mx.array(z[k]) for k in z.files if k.startswith('encoder::')};fm=z['feature_mean'];fs=z['feature_sd'];head_w=z['head.w1'];head_b=z['head.b1']
    random_w=np.random.default_rng(41207).normal(0,1/np.sqrt(512),(512,64)).astype('float32');data=[];old_features=[];new_features=[];known_AVI=[];manifest=[];avip=V3/'avi_exact_allele_scores.csv';audit=json.loads((V3/'avi_retrieval_audit.json').read_text());assert sha(avip)==audit['output_sha256'];lookup={r['variant_GRCh38']:r for r in csv.DictReader(avip.open())}
    for stem in ['salomon','siraj']:
        sp=ROOT/'work/sleep_frailty_ml_v3'/(stem+'_sequences_local.npz');fp=sp.with_name(stem+'_features_local.npz');receipt=next(r for r in json.loads((V3/'feature_receipt.json').read_text())['datasets'] if r['study']==stem);assert sha(sp)==receipt['input_sha256'] and sha(fp)==receipt['features_sha256']
        with np.load(sp) as z:d={k:z[k] for k in ['ref','alt','position','chromosome','label','activity_mean','absolute_effect']}
        with np.load(fp) as z:original=z['hyena_original']
        rows=list(csv.DictReader(sp.with_name(stem+'_manifest_local.csv').open()));ix=np.flatnonzero(d['chromosome']>=4);data.append({k:v[ix] for k,v in d.items() if k not in ['ref','alt','position']});original=original[ix];old_features.append(np.c_[original,np.tanh(np.clip((original-fm)/fs,-20,20)@random_w)])
        target=DATA/(stem+'_taught_assay_features.npz');signature=sha(path)+'_'+sha(sp)
        if target.exists():
            with np.load(target) as z:assert str(z['signature'])==signature and np.array_equal(z['row_index'],ix);new=z['features']
        elif status['selected_changed_native_scalars']==0:new=original.copy();np.savez_compressed(target,features=new,row_index=ix,signature=np.array(signature))
        else:
            compiled=mx.compile(lambda aa,bb,pp:model_mlx.strand_features(aa,bb,pp,w));new=np.empty_like(original);started=time.time()
            for start in range(0,len(ix),32):
                jj=ix[start:start+32];new[start:start+len(jj)]=np.array(compiled(mx.array(d['ref'][jj].astype('int32')),mx.array(d['alt'][jj].astype('int32')),mx.array(d['position'][jj])))
                if start%8192==0:print('Taught laboratory embeddings',stem,start,'/',len(ix),'seconds',round(time.time()-started,1),flush=True)
            np.savez_compressed(target,features=new,row_index=ix,signature=np.array(signature))
        new_features.append(np.c_[new,np.tanh(np.clip((new-fm)/fs,-20,20)@head_w+head_b)]);rr=[rows[i] for i in ix];manifest.extend(rr)
        for r in rr:
            q=lookup[r['variant_GRCh38']];known_AVI.append([float(q['avi_raw']),float(q['avi_phred']),1.] if q['status']=='exact_allele_match' else [np.nan,np.nan,0.])
    y=np.concatenate([r['label'] for r in data]);activity=np.concatenate([r['activity_mean'] for r in data]);effect=np.log1p(np.concatenate([r['absolute_effect'] for r in data]));chrom=np.concatenate([r['chromosome'] for r in data]);sid=np.concatenate([np.full(len(r['label']),i,dtype='int32') for i,r in enumerate(data)]);active=np.array([int(r['active_screen']) for r in manifest]);tr=np.flatnonzero(chrom>=6);va=np.flatnonzero((chrom>=4)&(chrom<=5));assert not np.any(chrom<=3);avi=np.array(known_AVI,dtype='float32');inputs={'original_control':np.c_[np.concatenate(old_features),avi],'AVI_taught':np.c_[np.concatenate(new_features),avi]};assert inputs['original_control'].shape==inputs['AVI_taught'].shape
    targets=np.c_[activity,effect];tm=np.stack([targets[tr[sid[tr]==s]].mean(0) for s in range(2)]);ts=np.stack([targets[tr[sid[tr]==s]].std(0) for s in range(2)]);ts[ts<1e-5]=1.;targets=np.clip((targets-tm[sid])/ts[sid],-5,5).astype('float32');pw=np.array([np.sqrt((np.sum(sid[tr]==s)-y[tr[sid[tr]==s]].sum())/y[tr[sid[tr]==s]].sum()) for s in range(2)],dtype='float32');sw=np.array([len(tr)/(2*np.sum(sid[tr]==s)) for s in range(2)],dtype='float32');pwm=mx.array(pw);swm=mx.array(sw);models=DATA/'models';history=[];results=[];predictions={}
    def metric(truth,p):
        ap=float(average_precision_score(truth,p));prev=float(np.mean(truth));return {'n':len(truth),'calls':int(np.sum(truth)),'AUROC':float(roc_auc_score(truth,p)),'AP':ap,'prevalence':prev,'APskill':(ap-prev)/(1-prev)}
    def evaluate(p):
        all_metrics={str(s):metric(y[va[sid[va]==s]],p[sid[va]==s]) for s in range(2)};within={str(s):metric(y[va[(sid[va]==s)&(active[va]==1)]],p[(sid[va]==s)&(active[va]==1)]) for s in range(2)};return {'studies':all_metrics,'within_active':within,'macro_APskill':float(np.mean([m['APskill'] for m in all_metrics.values()]))}
    def gelu(x):return .5*x*(1+mx.tanh(.7978845608028654*(x+.044715*x**3)))
    def huber(x):return mx.where(mx.abs(x)<=1,.5*x*x,mx.abs(x)-.5)
    def forward(p,x,s,masks=None):
        h=gelu(x@p['w1']+p['b1']);h=h if masks is None else h*masks[0];h=gelu(h@p['w2']+p['b2']);h=h if masks is None else h*masks[1];row=mx.arange(len(s));return [(h@p[k+'.w']+p[k+'.b'])[row,s] for k in ['class','activity','effect']]
    for name,x in inputs.items():
        mean=np.nanmean(x[tr],axis=0);sd=np.nanstd(x[tr],axis=0);sd[sd<1e-5]=1.;x=np.clip((np.where(np.isfinite(x),x,mean)-mean)/sd,-20,20).astype('float32');member_predictions=[]
        for seed in [42,73,101]:
            rng=np.random.default_rng(seed+100);r1=np.random.default_rng(seed+1);r2=np.random.default_rng(seed+2);r3=np.random.default_rng(seed+3);D=x.shape[1];p={'w1':mx.array(r1.normal(0,1/np.sqrt(D),(D,128)).astype('float32')),'b1':mx.zeros(128),'w2':mx.array(r2.normal(0,1/np.sqrt(128),(128,32)).astype('float32')),'b2':mx.zeros(32)}
            for k in ['class','activity','effect']:p[k+'.w']=mx.array(r3.normal(0,.01,(32,2)).astype('float32'));p[k+'.b']=mx.array(np.array([np.log(y[tr[sid[tr]==s]].mean()/(1-y[tr[sid[tr]==s]].mean())) for s in range(2)],dtype='float32')) if k=='class' else mx.zeros(2)
            def loss(pp,xx,s,yy,tt,masks):
                z,a,e=forward(pp,xx,s,masks);weight=swm[s];return mx.mean(weight*((1-yy)*mx.logaddexp(0,z)+pwm[s]*yy*mx.logaddexp(0,-z)))+.2*mx.mean(weight*huber(a-tt[:,0]))+.1*mx.mean(weight*huber(e-tt[:,1]))
            grad=mx.compile(mx.value_and_grad(loss));opt=optim.AdamW(.0003,weight_decay=.01);best=-np.inf;stale=0
            for epoch in range(1,31):
                values=[]
                for start in range(0,len(tr),512):
                    if start==0:order=rng.permutation(tr)
                    ix=order[start:start+512];masks=[mx.array(rng.binomial(1,.9,(len(ix),d)).astype('float32')/.9) for d in [128,32]];value,g=grad(p,mx.array(x[ix]),mx.array(sid[ix]),mx.array(y[ix]),mx.array(targets[ix]),masks);p=opt.apply_gradients(g,p);mx.eval(p,opt.state);values.append(float(value));assert np.isfinite(values[-1])
                pred=[]
                for start in range(0,len(va),512):
                    ix=va[start:start+512];z=np.array(forward(p,mx.array(x[ix]),mx.array(sid[ix]))[0]);pred.extend((1/(1+np.exp(-np.clip(z,-30,30)))).tolist())
                pred=np.array(pred);m=evaluate(pred);history.append({'variant':name,'seed':seed,'epoch':epoch,'loss':float(np.mean(values)),'validation':m});write_json(OUT/'assay_transfer_history.json',history)
                if m['macro_APskill']>best+1e-6:best=m['macro_APskill'];stale=0;bp=pred.copy();best_epoch=epoch;params={k:np.array(v) for k,v in p.items()}
                else:stale+=1
                if stale>=5:break
            target=models/(f'assay_probe_{name}_{seed}.npz');np.savez_compressed(target,**params,feature_mean=mean,feature_sd=sd);results.append({'variant':name,'seed':seed,'selected_epoch':best_epoch,'metrics':evaluate(bp),'checkpoint':str(target),'checkpoint_sha256':sha(target)});member_predictions.append(bp);print('MEASURED ASSAY DEVELOPMENT PROBE',name,seed,json.dumps(evaluate(bp)),flush=True)
        pred=np.mean(member_predictions,axis=0);results.append({'variant':name,'seed':'ensemble','metrics':evaluate(pred)});predictions[name]=pred
    ensembles={r['variant']:r for r in results if r['seed']=='ensemble'};delta=ensembles['AVI_taught']['metrics']['macro_APskill']-ensembles['original_control']['metrics']['macro_APskill'];write_json(OUT/'assay_transfer_results.json',{'results':results,'paired_macro_APskill_difference':delta,'teacher_checkpoint_sha256':status['checkpoint_sha256'],'fair_comparison':'Same579input dimensions, three seeds, optimizer, batches, targets and stopping criterion. Both receive actual permittedAVI input. Control has original encoder/untrained64d projection; adapted has selected taught encoder/trained projection.','scope':'Measured laboratory outcomes, but chromosomes4–5are reused development data. This is NOT a new independent validation study or patient genetic-score check.','next':'Expand only if the real-readout comparisons justify it; retain V3 unchanged and seek a new independent measured source.'});np.savez_compressed(OUT/'assay_transfer_validation_predictions.npz',row_index=va,label=y[va],study=sid[va],active=active[va],**predictions);print('ASSAY TRANSFER DEVELOPMENT COMPARISON COMPLETE',delta,flush=True)
if __name__=='__main__':main()
