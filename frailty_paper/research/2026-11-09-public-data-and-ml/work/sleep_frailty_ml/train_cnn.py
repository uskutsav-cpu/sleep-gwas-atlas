from pathlib import Path
import sys,json,time
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np
import mlx.core as mx
import mlx.optimizers as optim
import model_cnn as model

WORK=Path(__file__).resolve().parent;ROOT=WORK.parents[1];OUT=ROOT/'outputs/sleep_frailty_ai_v1'
mx.metal.set_memory_limit(1024**3);mx.metal.set_cache_limit(128*1024**2);mx.random.seed(73);rng=np.random.default_rng(73)
d=np.load(WORK/'prepared.npz');ref=d['ref'];alt=d['alt'];y=d['label'];chrom=d['chromosome'];tr=np.flatnonzero(chrom<=16);va=np.flatnonzero((chrom>=17)&(chrom<=18));te=np.flatnonzero(chrom>=19)
def metrics(y,p):
    ix=np.argsort(-p,kind='stable');ss=p[ix];yy=y[ix];ends=np.r_[np.flatnonzero(np.diff(ss)!=0),len(ss)-1];tp=np.cumsum(yy)[ends];fp=ends+1-tp;rec=tp/y.sum();precision=tp/(ends+1)
    return {'n':len(y),'detected_calls':int(y.sum()),'prevalence':float(y.mean()),'auroc':float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/(len(y)-y.sum())])),'average_precision':float(np.sum(np.diff(np.r_[0,rec])*precision)),'brier_uncalibrated':float(np.mean((p-y)**2))}
def sigmoid(x):return 1/(1+np.exp(-np.clip(x,-30,30)))
positive_weight=float(np.sqrt((len(tr)-y[tr].sum())/y[tr].sum()));p=model.initialize();optimizer=optim.AdamW(.001,weight_decay=.001)
def loss(p,a,b,t):
    z=model.logits(a,b,p);return mx.mean((1-t)*mx.logaddexp(0,z)+positive_weight*t*mx.logaddexp(0,-z))
vg=mx.compile(mx.value_and_grad(loss));best=-1;stale=0;history=[]
def predict(ix,pp):
    fn=mx.compile(lambda a,b:model.logits(a,b,pp));out=[]
    for a in range(0,len(ix),128):
        ii=ix[a:a+128];out.extend(sigmoid(np.array(fn(mx.array(ref[ii]),mx.array(alt[ii])))))
    return np.array(out)

for epoch in range(1,17):
    order=rng.permutation(tr);losses=[];start=time.time()
    for a in range(0,len(order),128):
        ii=order[a:a+128];ll,g=vg(p,mx.array(ref[ii]),mx.array(alt[ii]),mx.array(y[ii]));norm=mx.sqrt(sum(mx.sum(x*x) for x in g.values()));g={k:v*mx.minimum(1,1/(norm+1e-8)) for k,v in g.items()};p=optimizer.apply_gradients(g,p);mx.eval(p,optimizer.state);losses.append(float(ll))
    pred=predict(va,p);m=metrics(y[va],pred);history.append({'epoch':epoch,'train_weighted_bce':float(np.mean(losses)),'validation_average_precision':m['average_precision'],'validation_auroc':m['auroc'],'seconds':time.time()-start});(OUT/'cnn_training_history.json').write_text(json.dumps(history,indent=2)+'\n');print('CNN',epoch,history[-1],flush=True)
    if m['average_precision']>best+1e-6:
        best=m['average_precision'];best_p={k:mx.array(np.array(v)) for k,v in p.items()};best_val=pred.copy();best_epoch=epoch;stale=0;np.savez_compressed(OUT/'models/cnn.npz',**{k:np.array(v) for k,v in best_p.items()})
    else:stale+=1
    if stale>=3:break
(OUT/'cnn_selection.json').write_text(json.dumps({'selected_epoch':best_epoch,'validation_average_precision':best,'parameter_count':sum(v.size for v in best_p.values()),'architecture':'Shared 3-layer local DNA CNN, symmetric allele and strand pooling; trained from scratch','seed':73,'plan':'CNN_COMPARATOR_PLAN.md'},indent=2)+'\n')
print('CNN FITTING FINISHED; final test now evaluated without further fitting',flush=True)
pred=predict(te,best_p);m=metrics(y[te],pred)
order=np.argsort(-best_val);yy=y[va][order];ss=best_val[order];tp=np.cumsum(yy);precision=tp/np.arange(1,len(va)+1);recall=tp/y[va].sum();f1=2*precision*recall/np.maximum(precision+recall,1e-12);ends=np.r_[np.flatnonzero(np.diff(ss)!=0),len(ss)-1];threshold=float(ss[ends[np.argmax(f1[ends])]])
call=pred>=threshold;tp=int(np.sum(call&(y[te]==1)));fp=int(np.sum(call&(y[te]==0)));fn=int(np.sum(~call&(y[te]==1)))
m.update({'validation_selected_threshold':threshold,'test_precision':tp/max(tp+fp,1),'test_recall':tp/max(tp+fn,1),'test_true_positives':tp,'test_false_positives':fp})
(OUT/'cnn_test_results.json').write_text(json.dumps(m,indent=2)+'\n');np.savez_compressed(OUT/'cnn_test_predictions.npz',row_index=te,predictions=pred)
print('CNN final test results saved; no further fitting',flush=True)
