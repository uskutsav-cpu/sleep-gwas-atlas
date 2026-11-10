"""Numerical and raw-input checks on development rows before blind evaluation."""
from common import *
from inference_records import StoredFeatureInference
from raw_features import RawDNAInference
from final_fit import bio_forward,mx
from contrast_network import forward as contrast_forward
def main():
    path=OUT/'final_model_registry.json';registry=json.loads(path.read_text());assert registry['plan_sha256']==sha(OUT/'final_training_locked_plan.json');records=registry['records'];engine=StoredFeatureInference();tr=np.flatnonzero(engine.data['chromosome']>=4);assert registry['training_indices_sha256']==hashlib.sha256(tr.tobytes()).hexdigest() and registry['training_rows']==len(tr)
    raw_records=json.loads((OUT/'raw_inference_validation_example.json').read_text());raw=RawDNAInference(raw_records);receipt=json.loads((OUT/'raw_inference_check.json').read_text());ix=np.array(receipt['rows']);assert len(set(ix))==len(ix) and np.all(engine.data['chromosome'][ix]>=4);jj=np.arange(len(ix));checks={}
    for name,r in records.items():
        if r['artifact']:assert sha(OUT/r['artifact'])==r['checkpoint_sha256']
        pred=engine.predict(r,ix,records);rp=raw.predict(r,jj,records);errors={k:float(np.max(np.abs(pred[k]-rp[k]))) for k in pred};assert max(errors.values())<1e-4,(name,errors);gpu_error=None
        if r['family']=='biological_feature_head':
            with np.load(OUT/r['artifact']) as z:
                cfg=json.loads((OUT/r['config']['locked_config']).read_text());x=engine.matrix(r['features'],ix,cfg['AVI_as_input']);x=np.clip((np.where(np.isfinite(x),x,z['feature_mean'])-z['feature_mean'])/z['feature_sd'],-20,20).astype('float32');params={k:mx.array(z[k]) for k in z.files if k.startswith(('shared.','class.','activity.','effect.','avi.'))};logit=np.array(bio_forward(params,mx.array(x),mx.array(engine.data['study'][ix]))[0]);prob=1/(1+np.exp(-logit.astype('float64')));gpu_error=float(np.max(np.abs(prob-pred['overall'])));assert gpu_error<3e-5
        elif r['family']=='context_variant_mixture':
            with np.load(OUT/r['artifact']) as z:
                engine.features.load_transform(z);c,v=engine.features.transform(ix);params={k:mx.array(z[k]) for k in z.files if k.startswith(('context.','variant.','modulation.','gate.','conditional.','activity.','effect.','avi.'))};output=contrast_forward(params,mx.array(c),mx.array(v),mx.array(engine.data['study'][ix]),mx.array(z['log_positive_weight']));gpu_error=max(float(np.max(np.abs(np.array(output[k])-pred[k]))) for k in pred);assert gpu_error<3e-5
        elif r['family']=='conditional_avi_stacker':
            import pickle
            from conditional_stacking import inputs,predict_models
            base=engine.predict(records[r['config']['base_model']],ix,records);expert=base['conditional_scores'][:,1] if 'conditional_scores' in base else base['overall'];avi=engine.avi_features()[ix];sid=engine.data['study'][ix];include=r['config']['include_DNA_expert']
            with (OUT/r['artifact']).open('rb') as f:models=pickle.load(f)
            x,flag=inputs(expert,avi,include);manual=np.empty(len(ix))
            for s,m in models.items():
                use=np.flatnonzero(sid==s);q=m['quantile'].transform(m['imputer'].transform(x[use]));logit=(np.c_[q,flag[use]]@m['classifier'].coef_.T+m['classifier'].intercept_).ravel();manual[use]=1/(1+np.exp(-logit))
            error=float(np.max(np.abs(manual-pred['conditional_scores'][:,1])));assert error<1e-12;single=np.array([predict_models(models,expert[j:j+1],avi[j:j+1],sid[j:j+1],include)[0] for j in range(len(ix))]);assert np.max(np.abs(single-manual))<1e-12
        if 'changed_scalar_counts' in r:assert sum(r['changed_scalar_counts'].values())>0
        checks[name]={'raw_input_vs_stored_feature_errors':errors,'independent_CPU_vs_GPU_prediction_max_error':gpu_error,'checkpoint_hash_checked':bool(r['artifact'])}
    output={'models_checked':len(checks),'checks':checks,'registry_sha256':sha(path),'final_plan_sha256':registry['plan_sha256'],'development_rows_used':ix.tolist(),'no_reserved_test_scores_evaluated':True,'scope':'Implementation and checkpoint integrity checks on development rows, not generalization or biological validation.'};(OUT/'final_fit_verification.json').write_text(json.dumps(output,indent=2)+'\n');print('FINAL CHECKPOINT AND RAW INPUT VERIFICATION PASSED',len(checks),flush=True)
if __name__=='__main__':main()
