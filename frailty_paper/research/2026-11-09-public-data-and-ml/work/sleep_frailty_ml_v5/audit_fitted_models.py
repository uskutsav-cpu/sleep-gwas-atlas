"""Replay every development checkpoint and check every final fitted row."""
import subprocess
from train_development import prepare_runtime
from predict_checkpoints import predict_member,load_member
from metrics import point
from common import *
def main():
    assert not (OUT/'final_test_evaluation_started.json').exists();d,receipt,bank_x,pairs=prepare_runtime();va=np.flatnonzero(d['split']==1);cells=np.array([d['cell_names'][i] for i in d['cell'][va]]);records=json.loads((OUT/'development_records.json').read_text());assert len(records)==17;replay=[]
    for r in records:
        p,dist=predict_member(r,bank_x,d['bank_index'][va],cells);path=Path(r['checkpoint'])
        with np.load(path.with_name(path.stem+'_validation.npz')) as z:assert np.array_equal(z['row_index'],va);expected=z['prediction']
        error=float(np.max(np.abs(p-expected)));assert error<2e-5,(r['name'],r.get('seed'),error);replay.append({'family':r['name'],'seed':r.get('seed','tree'),'rows':len(p),'maximum_probability_error':error})
        print('Replayed development checkpoint',r['name'],r.get('seed','tree'),error,flush=True)
    registry=json.loads((OUT/'final_fit_registry.json').read_text());fit_checks=[]
    for r in registry['records']:
        member=load_member(r)
        if member['kind']=='neural':params=member['parameters'];ix=params['fit_row_indices'];assert all(np.isfinite(v).all() for k,v in params.items() if np.issubdtype(v.dtype,np.number));assert r['first_gradient_norm']>0 and r['changed_neural_parameter_scalars']>0
        else:ix=member['fit_row_indices']
        assert len(np.unique(ix))==len(ix) and np.all(d['split'][ix]<=1)
        if r['name']=='contrast_neural_only':assert all('neuronal_' in d['cell_names'][c] for c in d['cell'][ix])
        fit_checks.append({'family':r['name'],'seed':r.get('seed','tree'),'fitted_rows':len(ix),'calibration_rows_fitted':0,'checkpoint_sha256':r['checkpoint_sha256']})
    with np.load(OUT/'development_ensemble_predictions.npz') as z:assert np.array_equal(z['row_index'],va);pred={k:z[k] for k in z.files if k!='row_index'}
    rchecks=[];rroot=DATA/'independent_R_development';rroot.mkdir(exist_ok=True)
    for key in np.unique(d['context'][va]):
        for active_only in [False,True]:
            use=(d['context'][va]==key)&((d['active'][va]==1) if active_only else True);ids=va[use];truth=d['y'][ids]
            if len(np.unique(truth))<2:continue
            name=str(key).replace('|','_')+('_active' if active_only else '_all');src=rroot/(name+'.csv');dst=rroot/(name+'_metrics.csv')
            with src.open('w',newline='') as f:w=csv.writer(f);w.writerow(['row_id','label',*pred]);w.writerows(zip(ids,truth,*(p[use] for p in pred.values())))
            subprocess.run(['Rscript',str(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R'),str(src),str(dst)],check=True,capture_output=True,text=True)
            for row in csv.DictReader(dst.open()):
                m=point(truth,pred[row['model']][use]);error=max(abs(float(row['auroc'])-m['AUROC']),abs(float(row['average_precision'])-m['AP']));assert error<1e-12;rchecks.append({'context':str(key),'active_only':active_only,'model':row['model'],'maximum_metric_error':error})
    write_json(OUT/'final_checkpoint_audit.json',{'passed':True,'development_checkpoint_replays':replay,'final_fit_partition_and_weight_checks':fit_checks,'independent_R_point_metric_checks':rchecks,'independent_R_metric_count':2*len(rchecks),'data_source_plan_sha256':sha(OUT/'model_search_plan.json'),'reserved_outcomes_opened':False,'V3_registry_unchanged':sha(ROOT/'outputs/sleep_frailty_ai_v3/final_model_registry.json')==d['audit']['V3_registry_sha256']});print('FINAL CHECKPOINT AND INDEPENDENT METRIC AUDIT PASSED',len(replay),'checkpoints;',2*len(rchecks),'Rmetrics',flush=True)
if __name__=='__main__':main()
