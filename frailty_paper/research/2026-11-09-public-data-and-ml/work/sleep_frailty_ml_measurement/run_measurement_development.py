"""Finite, outcome-protected comparison; no added-input or fresh-data fishing."""
import subprocess,time
from measurement_common import *
from finite_measurement_plan import get_plan

def main():
    assert_prior_versions();plan=get_plan();assert json.loads((DEST/'folded_median_pretraining_audit.json').read_text())['passed'];out=DEST/'development';out.mkdir(exist_ok=True);start=time.monotonic();names=[]
    for recipe in plan['recipes']:
        for source in plan['held_studies']:
            name=f'{recipe}_{source}_{plan["seed"]}';names.append(name)
            if (out/(name+'_record.json')).exists():continue
            print('START MEASUREMENT FIT',name,flush=True)
            with (out/(name+'.log')).open('w') as handle:subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('train_measurement_development.py')),'--recipe',recipe,'--held-study',source],stdout=handle,stderr=subprocess.STDOUT,check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
            print('DONE MEASUREMENT FIT',name,flush=True)
    for source in plan['held_studies']:
        a,b=[json.loads((out/f'{recipe}_{source}_42_record.json').read_text()) for recipe in plan['recipes']]
        for key in ['training_index_sha256','prediction_index_sha256','initial_parameter_sha256','normalization','actual_batch_schedule_sha256','completed_epochs','actual_backprop_updates']:assert a[key]==b[key],(source,key)
    choices=[]
    for recipe in plan['recipes']:
        curves={s:json.loads((out/f'{recipe}_{s}_42_curve.json').read_text()) for s in plan['held_studies'][1:]}
        for epoch in range(1,plan['epochs']+1):
            per={s:r[epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s,r in curves.items()};choices.append({'recipe':recipe,'common_epoch':epoch,'mean_source_line_median_error_skill':float(np.mean(list(per.values()))),'per_source':per})
    best=max(choices,key=lambda r:r['mean_source_line_median_error_skill']);by={recipe:max([r for r in choices if r['recipe']==recipe],key=lambda r:r['mean_source_line_median_error_skill']) for recipe in plan['recipes']};held={recipe:json.loads((out/f'{recipe}_none_42_curve.json').read_text())[r['common_epoch']-1] for recipe,r in by.items()}
    write_json(DEST/'finite_measurement_assessment.json',{'all_16_fits_complete':True,'matched_parameters_rows_batches_normalization_verified':True,'selected':best,'best_by_recipe':by,'all_common_epochs':choices,'held_region_metrics':held,
               'primary_selection_is_measured_magnitude_MAE':True,'no_new_Ulirsch_effects_used':True,'new_independent_validation_established':False,'latent_component_or_noise_prediction_is_not_primary_validation':True,
               'plan_sha256':sha(DEST/'finite_measurement_plan.json'),'scientific_objective_complete':False,'elapsed_seconds':time.monotonic()-start});assert_prior_versions();print('FINITE MEASUREMENT COMPARISON COMPLETE',json.dumps(best),flush=True)

if __name__=='__main__':main()
