"""Validation-only model search using permitted AVI and measured neuron labels.

No benchmark evaluation is performed by this script. Saved predictions are
validation-only, so later fit choices cannot accidentally use benchmark results.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
from pathlib import Path
import sys,time,json,pickle,argparse,subprocess,copy
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from common import *

parser=argparse.ArgumentParser();parser.add_argument('--wait-for-avi-pid',type=int);args=parser.parse_args()
if args.wait_for_avi_pid:
    print('Waiting for verified AVI retrieval PID',args.wait_for_avi_pid,flush=True)
    while True:
        state=subprocess.run(['ps','-p',str(args.wait_for_avi_pid),'-o','command='],capture_output=True,text=True)
        if state.returncode or 'sleep_frailty_ml_v2/fetch_avi.py' not in state.stdout:break
        time.sleep(10)
data,seq,avi,manifest,tr,va,te=load();y=data['label'];MODELS.mkdir(exist_ok=True)
features={'avi':avi,'kmer':seq['kmer'],'kmer_avi':np.c_[seq['kmer'],avi],'hyena':seq['hyena_original'],'hyena_avi':np.c_[seq['hyena_original'],avi]}
results=[];predictions={};history=[]
def record(name,family,feature,config,p,artifact):
    row={'name':name,'family':family,'features':feature,'config':config,'validation':metrics(y[va],p),'artifact':artifact}
    results.append(row);predictions[name]=np.asarray(p,dtype='float64')
    (OUT/'search_results.json').write_text(json.dumps(results,indent=2)+'\n')
    np.savez_compressed(OUT/'validation_predictions.npz',row_index=va,label=y[va],**predictions)
    print('Validation',name,json.dumps(row['validation']),flush=True)
start=time.time()
# Baseline from v1 is eligible for the recommendation and must be replayed.
sys.path.insert(0,str(OLD));from predict import predict as old_predict
alphabet=np.array(list('ACGT'))
vr=[{'variant':str(i),'reference_sequence':''.join(alphabet[data['ref'][i]-7]),'alternate_sequence':''.join(alphabet[data['alt'][i]-7])} for i in va]
v1=old_predict(ROOT/'outputs/sleep_frailty_ai_v1/models/kmer_baseline.npz',vr,'numpy')
record('v1_simple_baseline','v1_logistic','v1_kmer',{},v1,'../sleep_frailty_ai_v1/models/kmer_baseline.npz')
raw=np.where(np.isfinite(avi[va,0]),avi[va,0],np.nanmedian(avi[tr,0]))
# Direct score ranking is a comparator, not a calibrated probability.
direct=1/(1+np.exp(-np.clip(raw,-30,30)))
record('avi_direct_rank','direct_avi','avi',{'missing_imputation':'training median','score_transform':'monotone sigmoid of raw AVI, not probability'},direct,None)

for feature in ['avi','kmer','kmer_avi']:
    for C in [.01,.1,1.,10.]:
        name=f'logistic_{feature}_C{C:g}';path=MODELS/(name+'.pkl')
        model=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),LogisticRegression(C=C,max_iter=1500,class_weight=None,random_state=42))
        model.fit(features[feature][tr],y[tr]);p=model.predict_proba(features[feature][va])[:,1]
        with path.open('wb') as f:pickle.dump(model,f)
        record(name,'logistic',feature,{'C':C,'class_weight':None},p,str(path.relative_to(OUT)))
for feature in ['kmer','kmer_avi']:
    for depth in [2,3,5]:
        for lr in [.03,.08]:
            name=f'boost_{feature}_d{depth}_lr{lr:g}';path=MODELS/(name+'.pkl')
            config={'max_depth':depth,'learning_rate':lr,'max_iter':200,'min_samples_leaf':40,'l2_regularization':10.,'max_bins':64,'early_stopping':False,'random_state':42}
            model=HistGradientBoostingClassifier(**config);model.fit(features[feature][tr],y[tr]);p=model.predict_proba(features[feature][va])[:,1]
            with path.open('wb') as f:pickle.dump(model,f)
            record(name,'hist_gradient_boosting',feature,config,p,str(path.relative_to(OUT)))

# Neural weights really optimize by reverse-mode differentiation and AdamW.
sys.path.insert(0,str(OLD/'runtime'));import mlx.core as mx;import mlx.optimizers as optim
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
def sigmoid(z):return 1/(1+np.exp(-np.clip(np.asarray(z,dtype='float64'),-30,30)))
positive_weight=float(np.sqrt((len(tr)-y[tr].sum())/y[tr].sum()))
def forward(x,p):return mx.tanh(x@p['w1']+p['b1'])@p['w2']+p['b2']
for feature in ['hyena','hyena_avi']:
    imp=SimpleImputer(strategy='median',keep_empty_features=True);x=imp.fit_transform(features[feature][tr]);scaler=StandardScaler().fit(x)
    xt=np.clip(scaler.transform(x),-20,20).astype('float32');xv=np.clip(scaler.transform(imp.transform(features[feature][va])),-20,20).astype('float32')
    for width in [32,128]:
        names=[]
        for seed in [42,73,101]:
            rng=np.random.default_rng(seed);mx.random.seed(seed);name=f'neural_{feature}_w{width}_s{seed}';names.append(name);path=MODELS/(name+'.npz')
            params={'w1':mx.array(rng.normal(0,1/np.sqrt(xt.shape[1]),(xt.shape[1],width)).astype('float32')),'b1':mx.zeros(width),'w2':mx.array(rng.normal(0,.01,width).astype('float32')),'b2':mx.array(float(np.log(y[tr].mean()/(1-y[tr].mean()))))}
            initial={k:np.array(v) for k,v in params.items()};opt=optim.AdamW(.0003,weight_decay=.01)
            def loss(p,x,yb):
                z=forward(x,p);return mx.mean((1-yb)*mx.logaddexp(0,z)+positive_weight*yb*mx.logaddexp(0,-z))
            grad=mx.compile(mx.value_and_grad(loss));best=-1;stale=0
            for epoch in range(1,51):
                order=rng.permutation(len(tr));losses=[]
                for j in range(0,len(order),512):
                    ix=order[j:j+512];value,g=grad(params,mx.array(xt[ix]),mx.array(y[tr[ix]]));params=opt.apply_gradients(g,params);mx.eval(params,opt.state);losses.append(float(value))
                p=sigmoid(np.array(forward(mx.array(xv),params)));m=metrics(y[va],p)
                history.append({'name':name,'epoch':epoch,'training_loss':float(np.mean(losses)),'validation_average_precision':m['average_precision'],'validation_auroc':m['auroc']})
                (OUT/'neural_training_history.json').write_text(json.dumps(history,indent=2)+'\n')
                if m['average_precision']>best+1e-6:
                    best=m['average_precision'];stale=0;best_p=p.copy();best_epoch=epoch;saved={k:np.array(v) for k,v in params.items()}
                else:stale+=1
                if stale>=7:break
            changes={k:float(np.linalg.norm(saved[k]-initial[k])) for k in initial};assert changes['w1']>0 and changes['w2']>0
            np.savez_compressed(path,**saved,imputation=imp.statistics_,mean=scaler.mean_,sd=scaler.scale_)
            record(name,'neural_head',feature,{'width':width,'seed':seed,'best_epoch':best_epoch,'maximum_epochs':50,'patience':7,'lr':.0003,'weight_decay':.01,'positive_weight':positive_weight,'parameter_change_norms':changes,'frozen_encoder':True},best_p,str(path.relative_to(OUT)))
        record(f'ensemble_{feature}_w{width}','neural_seed_ensemble',feature,{'members':names,'combination':'arithmetic_mean'},np.mean([predictions[n] for n in names],axis=0),None)

best=max(results,key=lambda row:row['validation']['average_precision'])
selection={'selected_model':best['name'],'rule':'Largest validation average precision across documented candidates; benchmark not evaluated by this selector','selected_record':best,'candidate_count':len(results),'training_seconds':time.time()-start,'development_limitations':['Validation reused for a finite multi-model search with only 78 positive calls; optimism likely','Benchmark chromosomes 19-22 were already observed in v1; external assay validation remains required','Only permissible downloaded AVI scores used; their training/benchmark overlap is unresolved','Not genetic-score validation, disease causality or patient risk prediction'],'input_hashes':{'prepared':sha(OLD/'prepared.npz'),'sequence_features':sha(WORK/'sequence_features.npz'),'avi_table':sha(OUT/'avi_exact_allele_scores.csv'),'plan':sha(WORK/'EXPERIMENT_PLAN.md')}}
(OUT/'model_selection.json').write_text(json.dumps(selection,indent=2)+'\n')
print('SEARCH COMPLETE',json.dumps(selection,indent=2),flush=True)
