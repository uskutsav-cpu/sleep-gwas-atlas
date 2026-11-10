"""Finite original-direction ablation, preserving the new outside test."""
import subprocess,time
from direction_common import *
from finite_direction_plan import get_plan


def main():
    assert_versions();plan=get_plan();out=DEST/'development';out.mkdir(exist_ok=True);start=time.monotonic()
    for source in plan['held_studies']:
        if source=='Siraj': continue
        name=f'source_signed_{source}_42'
        if (out/(name+'_record.json')).exists(): continue
        print('START ORIGINAL DIRECTION FIT',source,flush=True)
        with (out/(name+'.log')).open('w') as handle:
            subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('train_direction_development.py')),'--recipe','source_signed','--held-study',source],stdout=handle,stderr=subprocess.STDOUT,check=True,
               env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
        print('DONE ORIGINAL DIRECTION FIT',source,flush=True)
    curves={};matched=[]
    for source in plan['held_studies']:
        control=MEASUREMENT/'development'/f'observation_only_{source}_42'
        candidate=out/f'source_signed_{source}_42' if source!='Siraj' else control
        a,b=[json.loads(Path(str(path)+'_record.json').read_text()) for path in [control,candidate]]
        for key in ['training_index_sha256','prediction_index_sha256','initial_parameter_sha256','normalization','actual_batch_schedule_sha256','completed_epochs','actual_backprop_updates']:
            assert a[key]==b[key],(source,key)
        matched.append(source)
        curves[source]={arm:json.loads(Path(str(path)+'_curve.json').read_text()) for arm,path in [('control',control),('source_signed',candidate)]}
    choices=[]
    for epoch in range(1,plan['epochs']+1):
        skills={s:curves[s]['source_signed'][epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in plan['held_studies'][1:]}
        controls={s:curves[s]['control'][epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in plan['held_studies'][1:]}
        choices.append({'common_epoch':epoch,'candidate_mean_source_skill':float(np.mean(list(skills.values()))),
           'control_mean_source_skill':float(np.mean(list(controls.values()))),'per_source_candidate':skills,'per_source_control':controls})
    best=max(choices,key=lambda x:x['candidate_mean_source_skill'])
    write_json(DEST/'finite_direction_assessment.json',{'all_seven_new_fits_complete':True,'Siraj_omission_control_reused_by_exact_training_equivalence':True,
       'matched_initialization_rows_batches_and_normalization_verified':matched,'selected':best,'all_common_epochs':choices,
       'new_outside_test_used':False,'scientific_goal_complete':False,'plan_sha256':sha(DEST/'finite_direction_plan.json'),
       'elapsed_seconds':time.monotonic()-start})
    assert_versions();print('FINITE ORIGINAL-DIRECTION COMPARISON COMPLETE',json.dumps(best),flush=True)


if __name__=='__main__': main()
