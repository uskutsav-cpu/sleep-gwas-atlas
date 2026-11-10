"""Fixed simpler quantitative baselines; all scalers fit training loci only."""
import gc,joblib
from train import *
from sklearn.linear_model import Ridge
from sklearn.ensemble import HistGradientBoostingRegressor
COLS=np.r_[np.arange(32),np.arange(512,680),np.arange(835,1145)].astype('int32')
def compact(d,ids,t,avi=None):
    a=np.asarray(d['x'][d['bank_index'][ids]][:,COLS],dtype='float32');a=np.clip((a-t['mean'][COLS])/t['sd'][COLS],-12,12);parts=[a,flags(d,ids,t)]
    if avi is not None:
        v=np.array(avi[d['bank_index'][ids]],copy=True);v[:,:2]=np.where(np.isfinite(v[:,:2]),v[:,:2],t['avi_median']);v[:,:2]=(v[:,:2]-t['avi_mean'])/t['avi_sd'];parts.append(np.clip(v,-12,12))
    return np.concatenate(parts,1).astype('float32')
def predict_classical(d,ids,t,model,avi=None):
    parts=[]
    for start in range(0,len(ids),2048):parts.append(np.expm1(np.clip(model.predict(compact(d,ids[start:start+2048],t,avi)),0,12)))
    return np.concatenate(parts)
def fit_classical(d,fit,va,family,stage='development',avi=None,omit=None):
    name=f'{stage}_{family}'+('_avi' if avi is not None else '')+(f'_omit{omit}' if omit is not None else '');record=OUT/(name+'_record.json')
    if record.exists():
        r=json.loads(record.read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256'];return r
    t=fit_transform(d,fit)
    if avi is not None:t=add_avi_transform(d,fit,t,avi)
    ids=fit[d['unit'][fit]==0];dim=len(COLS)+d['n_context_flags']+(3 if avi is not None else 0);y=d['log_effect'][ids];w=t['weights'][ids];started=time.time()
    if family.startswith('ridge_'):
        # Exact weighted least-squares sufficient statistics avoid retaining
        # a second complete observation-by-feature array in eight-GB memory.
        signature=hashlib.sha256(ids.tobytes()+np.asarray([avi is not None],dtype='uint8').tobytes()).hexdigest();cache=DATA/'runtime_tmp'/('ridge_gram_'+signature+'.npz');cache.parent.mkdir(exist_ok=True)
        if cache.exists():
            with np.load(cache) as z:gram=z['gram'];xy=z['xy'];sx=z['sx'];sy=float(z['sy']);sw=float(z['sw'])
        else:
            gram=np.zeros((dim,dim),dtype='float64');xy=np.zeros(dim);sx=np.zeros(dim);sy=0.;sw=0.
            for start in range(0,len(ids),2048):
                stop=min(start+2048,len(ids));a=compact(d,ids[start:stop],t,avi).astype('float64');ww=w[start:stop].astype('float64');yy=y[start:stop].astype('float64');weighted=a*ww[:,None];gram+=a.T@weighted;xy+=weighted.T@yy;sx+=weighted.sum(0);sy+=float(ww@yy);sw+=float(ww.sum())
            np.savez_compressed(cache,gram=gram,xy=xy,sx=sx,sy=sy,sw=sw)
        penalty=float(family.split('_')[-1]);centered=gram-np.outer(sx,sx)/sw;centered.flat[::dim+1]+=penalty;coef=np.linalg.solve(centered,xy-sx*sy/sw);model=Ridge(alpha=penalty);model.coef_=coef;model.intercept_=sy/sw-(sx/sw)@coef;model.n_features_in_=dim
    else:
        scratch=DATA/'runtime_tmp'/(name+'_input.npy');a=np.lib.format.open_memmap(scratch,mode='w+',shape=(len(ids),dim),dtype='float32')
        for start in range(0,len(ids),2048):a[start:start+2048]=compact(d,ids[start:start+2048],t,avi)
        a.flush();model=HistGradientBoostingRegressor(max_depth=int(family.split('_')[-1]),max_iter=160,learning_rate=.08,max_leaf_nodes=31,min_samples_leaf=50,l2_regularization=10.,max_bins=63,early_stopping=False,random_state=42);model.fit(a,y,sample_weight=w);del a;gc.collect();scratch.unlink()
    gc.collect();path=MODELS/(name+'.joblib');joblib.dump(model,path,compress=3);tx=MODELS/(name+'_transform.npz');transform_save(tx,t)
    if stage=='development':
        p=predict_classical(d,va,t,model,avi);np.savez_compressed(MODELS/(name+'_validation.npz'),row_index=va,prediction=p);m=evaluate(d,fit,va,p);write_json(OUT/(name+'_metrics.json'),m)
    r={'name':name,'recipe':family,'stage':stage,'selected_epoch':None,'fit_rows':len(ids),'AVI_added_input':avi is not None,'omitted_study':omit,'checkpoint':str(path),'checkpoint_sha256':sha(path),'transform':str(tx),'transform_sha256':sha(tx),'seconds':time.time()-started,'final_study_loaded':False};write_json(record,r);print('BASELINE FIT COMPLETE',name,r['seconds'],flush=True);return r
def main():
    freeze_plan();d=load();fit=np.flatnonzero(d['split']==0);va=np.flatnonzero(d['split']==1)
    write_json(OUT/'development_cell_median_metrics.json',evaluate(d,fit,va,baseline(d,fit,va)))
    for family in ['ridge_100','ridge_1000','ridge_10000','tree_3','tree_5']:fit_classical(d,fit,va,family)
if __name__=='__main__':main()
