"""Reuse identical encodings/control fits, then serially train the blood arm."""
import argparse,subprocess,gzip
from train_development import *
from predict_checkpoints import predict_member
PRIOR=DATA.parent/'sleep_frailty_ml_v5';PRIOR_OUT=PRIOR/'outputs'
def reuse(d,receipt,bank_x,pairs):
    old_file=PRIOR/'measured_observations.csv';new_file=DATA/'measured_observations.csv';h=hashlib.sha256()
    with new_file.open('rb') as f:
        remaining=old_file.stat().st_size
        while remaining:
            b=f.read(min(1024*1024,remaining));assert b;h.update(b);remaining-=len(b)
    assert h.hexdigest()==sha(old_file)
    old_n=json.loads((PRIOR_OUT/'measured_data_audit.json').read_text())['observations'];oldtr=np.flatnonzero((d['split']==0)&(np.arange(d['n'])<old_n));newtr=np.flatnonzero(d['split']==0);assert np.array_equal(np.unique(d['bank_index'][oldtr]),np.unique(d['bank_index'][newtr]));cache=DATA/'standardized_banks';cache.mkdir(exist_ok=True);cached=[]
    records=json.loads((PRIOR_OUT/'development_records.json').read_text());assert len(records)==17
    for name in ['local_mlp','biological_mlp','contrast_standard','contrast_neural_only','tree_depth3']:
        r=next(r for r in records if r['name']==name);path=Path(r['checkpoint'])
        if path.suffix=='.npz':
            with np.load(path) as z:cols=z['feature_columns'];mu=z['feature_mean'];sd=z['feature_sd'];oldfit=z['fit_row_indices'];recipe=json.loads(str(z['recipe_json']))
        else:
            import pickle
            with path.open('rb') as f:p=pickle.load(f)
            cols=p['columns'];mu=p['mean'];sd=p['sd'];oldfit=p['fit_row_indices'];recipe={'view':'pattern_motif','architecture':'mlp','neural_only':False}
        fit=newtr
        if recipe['neural_only']:fit=fit[np.array(['neuronal_' in d['cell_names'][i] for i in d['cell'][fit]])]
        assert np.array_equal(np.unique(d['bank_index'][oldfit]),np.unique(d['bank_index'][fit]))
        def signature(ids):return hashlib.sha256((receipt['feature_sha256']+ids.tobytes().hex()+recipe['view']+recipe['architecture']+mu.tobytes().hex()+sd.tobytes().hex()).encode()).hexdigest()
        oldsig=signature(oldfit);newsig=signature(fit);stem='development_'+recipe['view']+'_'+recipe['architecture']+'_';source=PRIOR/'standardized_banks'/(stem+oldsig[:16]+'.npy');dest=cache/(stem+newsig[:16]+'.npy');meta=json.loads(source.with_suffix('.json').read_text());assert meta['signature']==oldsig
        xx=np.load(source,mmap_mode='r');maximum=0.
        for start in range(0,len(bank_x),2048):
            ids=np.arange(start,min(start+2048,len(bank_x)));expected=np.clip((bank_x[np.ix_(ids,cols)]-mu)/sd,-20,20);maximum=max(maximum,float(np.max(np.abs(expected-xx[start:start+len(ids)]))))
        assert maximum==0.,(name,maximum)
        if not dest.exists():dest.symlink_to(source)
        write_json(dest.with_suffix('.json'),{'signature':newsig,'shape':list(xx.shape),'fitted_observations':len(fit),'no_outcome_inputs':True,'reused_from':str(source),'sha256':sha(source),'all_rows_recomputed_and_compared':len(bank_x),'maximum_difference':maximum,'same_unique_sequence_training_indices':True})
        cached.append({'name':name,'file':str(dest),'source':str(source),'new_signature':newsig,'verified_rows':len(bank_x),'maximum_error':maximum});print('Reused fully checked fixed normalization',name,flush=True)
    va=np.flatnonzero(d['split']==1);cells=np.array([d['cell_names'][i] for i in d['cell'][va]]);reused=[]
    for r in records:
        if r['name']!='contrast_neural_only':continue
        with np.load(r['checkpoint']) as z:fit=z['fit_row_indices'];expected=np.flatnonzero((d['split']==0)&np.array(['neuronal_' in d['cell_names'][i] for i in d['cell']]));assert np.array_equal(fit,expected);assert z['cell_names'].tolist()==d['cell_names'] and z['study_names'].tolist()==d['study_names']
        p,_=predict_member(r,bank_x,d['bank_index'][va],cells);source=Path(r['checkpoint']);target=DATA/'models'/source.name;target.parent.mkdir(exist_ok=True)
        if not target.exists():target.symlink_to(source)
        np.savez_compressed(target.with_name(target.stem+'_validation.npz'),row_index=va,prediction=p)
        updated={**r,'checkpoint':str(target),'validation':evaluate(d,va,p,pairs,True),'reused_from_prior_version':{'record':str(PRIOR_OUT/f'development_contrast_neural_only_{r["seed"]}_record.json'),'reason':'New K562 rows never enter this neural-only fit. All actual fitting row indices, sequence values, targets, cells, source indices, normalization, initialization, minibatches and neural-only stopping target are unchanged. Validation predictions recomputed for the expanded full validation universe.','original_checkpoint_sha256':r['checkpoint_sha256']}}
        write_json(OUT/f'development_contrast_neural_only_{r["seed"]}_record.json',updated);reused.append({'seed':r['seed'],'checkpoint_sha256':r['checkpoint_sha256'],'fitting_rows':len(fit),'expanded_validation_rows':len(va)})
    write_json(OUT/'immutable_reuse_audit.json',{'passed':True,'old_observation_prefix_sha256':h.hexdigest(),'prior_observations':old_n,'new_observations':d['n'],'fixed_feature_bank_sha256':receipt['feature_sha256'],'normalizations':cached,'neural_only_fit_reuse':reused,'original_V5_outcomes_or_final_test_used':False})
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prior-training-pid',type=int,default=41440);a=parser.parse_args()
    while True:
        r=subprocess.run(['ps','-p',str(a.prior_training_pid),'-o','command='],capture_output=True,text=True)
        if r.returncode or 'sleep_frailty_ml_v5/train_development.py' not in r.stdout:break
        time.sleep(10)
    assert not (OUT/'final_test_evaluation_started.json').exists();assert not (PRIOR_OUT/'final_test_evaluation_started.json').exists()
    subprocess.run([sys.executable,str(ROOT/'work/sleep_frailty_ml_v5/select_and_refit.py'),'--select-only'],check=True)
    subprocess.run([sys.executable,str(WORK/'audit_partitions.py')],check=True)
    freeze_plan();d,receipt,bank_x,pairs=prepare_runtime();reuse(d,receipt,bank_x,pairs);del d,bank_x,pairs
    subprocess.run([sys.executable,'-u',str(WORK/'train_development.py')],check=True)
    subprocess.run([sys.executable,'-u',str(WORK/'select_and_refit.py')],check=True)
    print('EXPANDED BLOOD-ARM TRAINING AND FIXED REFITS COMPLETE; CALIBRATION/FINAL STUDY STILL REQUIRED',flush=True)
if __name__=='__main__':main()
