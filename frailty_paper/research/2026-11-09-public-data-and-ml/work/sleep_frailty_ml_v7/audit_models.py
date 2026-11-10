"""Replay actual final checkpoints and cross-study arithmetic independently."""
import subprocess
from finalize import *
def main():
    lock=json.loads((OUT/'final_model_and_prediction_lock.json').read_text());results=json.loads((OUT/'reserved_study_results.json').read_text());d=load();ids=np.flatnonzero((d['split']==1)&(d['unit']==0))[:32];checks=[]
    audits={arm:lock[arm]['records'] for arm in ['sequence_only','AVI_added']}
    audits['development_controls']=[{'record':json.loads(p.read_text())} for p in OUT.glob('development_*_record.json') if json.loads(p.read_text())['recipe'] in freeze_plan()['recipes']]
    for arm,members in audits.items():
        avi=np.load(DATA/'avi_bank.npy',mmap_mode='r') if arm=='AVI_added' else None
        for member in members:
            r=member['record'];assert sha(r['checkpoint'])==r['checkpoint_sha256'];t=transform_load(Path(r['transform']));assert sha(r['transform'])==r['transform_sha256']
            if r['recipe'] not in freeze_plan()['recipes']:checks.append({'checkpoint':r['checkpoint'],'hash_verified':True,'kind':'classical comparator'});continue
            with np.load(r['checkpoint']) as z:p={k:z[k] for k in z.files}
            config=freeze_plan()['recipes'][r['recipe']];c,v=features(d,ids,t,config['no_cell'],avi);cpu=forward(p,c,v,backend='numpy');gpu=np.array(forward({k:mx.array(a) for k,a in p.items()},mx.array(c),mx.array(v)));err=float(np.max(abs(cpu-gpu)));assert err<1e-4,(r['name'],err);assert np.max(abs(forward(p,c,np.zeros_like(v),backend='numpy')))==0
            # Gradient direction check on a real fitted continuous head.
            pm={k:mx.array(a) for k,a in p.items()};u=d['unit'][ids].astype('int32');y=d['log_effect'][ids]/t['effect_scales'][u];w=t['weights'][ids];masks=[mx.ones((len(ids),128)),mx.ones((len(ids),64)),mx.ones((len(ids),128)),mx.ones((len(ids),64))]
            def objective(pp):return loss(pp,mx.array(c),mx.array(v),mx.array(u),mx.array(y),mx.array(w),masks,len(ids),0,0,0.,config.get('raw_mae',False),mx.array(t['effect_scales']))
            value,g=mx.value_and_grad(objective)(pm);rng=np.random.default_rng(942);direction=rng.normal(size=p['effect_w'].shape).astype('float32');direction/=np.linalg.norm(direction);eps=.0001;plus=dict(pm);minus=dict(pm);plus['effect_w']=pm['effect_w']+eps*mx.array(direction);minus['effect_w']=pm['effect_w']-eps*mx.array(direction);numeric=(float(objective(plus))-float(objective(minus)))/(2*eps);analytic=float(np.sum(np.array(g['effect_w'])*direction));absolute=abs(numeric-analytic);assert absolute<max(.003,.1*abs(analytic)),(numeric,analytic);checks.append({'checkpoint':r['checkpoint'],'MLX_NumPy_max_error':err,'fitted_head_numeric_gradient':numeric,'fitted_head_backprop_gradient':analytic,'gradient_absolute_error':absolute,'zero_contrast_magnitude_exactly_zero':True,'actual_backprop_updates':r['actual_backprop_updates']})
    csvout=OUT/'independent_R_results.csv';subprocess.run(['/opt/homebrew/bin/Rscript',str(WORK/'check_results.R'),str(OUT/'reserved_all_observations.csv'),str(csvout)],check=True,env=environment());independent=[]
    for row in csv.DictReader(csvout.open()):
        original=results['arms'][row['arm']]['cells'][row['cell']];errors={k:abs(float(row[k])-original[k]) for k in ['n','MAE','RMSE','baseline_MAE','local_mean_Spearman','local_groups','interval_coverage','median_interval_width','ranked_fraction','ranked_observations']};errors['global_Spearman']=abs(float(row['global_Spearman'])-original['global_Spearman']);assert max(errors.values())<1e-8,errors;independent.append({'arm':row['arm'],'cell':row['cell'],'maximum_error':max(errors.values())})
    additional=json.loads((OUT/'additional_neural_results.json').read_text());extra_csv=OUT/'independent_R_additional_results.csv';subprocess.run(['/opt/homebrew/bin/Rscript',str(WORK/'check_additional.R'),str(OUT/'additional_neural_observations.csv'),str(extra_csv)],check=True,env=environment());extra_checks=[]
    for row in csv.DictReader(extra_csv.open()):
        original=additional['conditions'][row['condition']][row['arm']];errors={k:abs(float(row[k])-original[k]) for k in ['n','MAE','RMSE','baseline_MAE','global_Spearman','local_mean_Spearman','local_groups','interval_coverage','median_interval_width','ranked_fraction','ranked_observations']};assert max(errors.values())<1e-8,errors;extra_checks.append({'arm':row['arm'],'condition':row['condition'],'maximum_error':max(errors.values())})
    write_json(OUT/'final_checkpoint_and_R_audit.json',{'passed':True,'fitted_models':checks,'independent_R':independent,'independent_R_additional_source':extra_checks,'R_results_sha256':sha(csvout),'additional_R_results_sha256':sha(extra_csv),'previous_V6_lock':assert_previous_intact(),'scope':'Implementation and arithmetic verification, not independent biology or participant-level validation.'});print('FINAL CHECKPOINT AND INDEPENDENT R AUDIT PASSED',len(checks),flush=True)
if __name__=='__main__':main()
