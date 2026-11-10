"""One finite projection comparison, without source or architecture search loops."""
import subprocess,time
from channel_common import *

def main():
    assert_frozen_previous();plan=get_plan();assert json.loads((DEST/'channel_feature_receipt.json').read_text())['passed'];directory=DEST/'development';directory.mkdir(exist_ok=True);start=time.monotonic();names=[]
    for recipe in plan['recipes']:
        for source in plan['held_studies']:
            name=f'{recipe}_{source}_{plan["seed"]}';names.append(name)
            if (directory/(name+'_record.json')).exists():continue
            print('START CHANNEL FIT',name,flush=True)
            with (directory/(name+'.log')).open('w') as handle:subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('train_channel_development.py')),'--recipe',recipe,'--held-study',source],check=True,stdout=handle,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
            print('DONE CHANNEL FIT',name,flush=True)
    for source in plan['held_studies']:
        a,b=[json.loads((directory/f'{recipe}_{source}_42_record.json').read_text()) for recipe in plan['recipes']]
        for key in ['training_index_sha256','prediction_index_sha256','initial_parameter_sha256','actual_batch_schedule_sha256','completed_epochs','actual_backprop_updates','scalar_feature_mean','scalar_feature_sd','channel_feature_mean','channel_feature_sd']:assert a[key]==b[key],(source,key)
    choices=[]
    for recipe in plan['recipes']:
        curves={s:json.loads((directory/f'{recipe}_{s}_42_curve.json').read_text()) for s in plan['held_studies'][1:]}
        for epoch in range(1,plan['epochs']+1):
            per={s:curve[epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s,curve in curves.items()};choices.append({'recipe':recipe,'common_epoch':epoch,'mean_source_line_median_error_skill':float(np.mean(list(per.values()))),'per_source':per})
    best=max(choices,key=lambda r:r['mean_source_line_median_error_skill']);best_by_recipe={recipe:max([r for r in choices if r['recipe']==recipe],key=lambda r:r['mean_source_line_median_error_skill']) for recipe in plan['recipes']};held={}
    for recipe,selected in best_by_recipe.items():held[recipe]=json.loads((directory/f'{recipe}_none_42_curve.json').read_text())[selected['common_epoch']-1]
    write_json(DEST/'finite_channel_assessment.json',{'all_16_fits_complete':True,'matched_parameters_rows_batches_and_normalization_pass':True,'selected':best,'best_by_recipe':best_by_recipe,'all_common_epochs':choices,'held_region_metrics':held,
            'plan_sha256':sha(DEST/'finite_channel_plan.json'),'previous_adverse_tests_preserved':True,'new_untouched_validation_performed':False,'scientific_objective_complete':False,'elapsed_seconds':time.monotonic()-start})
    assert_frozen_previous();print('FINITE CHANNEL COMPARISON COMPLETE',json.dumps(best),flush=True)

if __name__=='__main__':main()
