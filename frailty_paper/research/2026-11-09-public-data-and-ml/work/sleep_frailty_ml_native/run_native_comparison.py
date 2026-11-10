"""Run the fixed gated comparison; no expansion following a failed neural gate."""
import subprocess,time
from native_common import *


def fit(recipe,source,plan):
    directory=DEST/'development';directory.mkdir(exist_ok=True);name=f'{recipe}_{source}_42'
    if (directory/(name+'_record.json')).exists():return
    assert not (directory/(name+'_started.json')).exists(),('Incomplete existing fit requires process/checkpoint inspection; do not silently restart',name)
    print('START NATIVE COMPARISON FIT',name,flush=True)
    with (directory/(name+'.log')).open('w') as handle:
        subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('train_native_development.py')),
           '--recipe',recipe,'--held-study',source],stdout=handle,stderr=subprocess.STDOUT,check=True,
           env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
    print('DONE NATIVE COMPARISON FIT',name,flush=True)


def matched_roles(roles):
    for source in roles:
        records=[json.loads((DEST/'development'/f'{recipe}_{source}_42_record.json').read_text()) for recipe in ['frozen_reader','adapt_last_stage']]
        for key in ['training_index_sha256','prediction_index_sha256','initial_parameter_sha256','normalization','context_transformations',
              'actual_batch_schedule_sha256','actual_backprop_updates','completed_epochs','new_original_Siraj_signed_training_rows','new_original_neuronal_signed_training_rows']:
            assert records[0][key]==records[1][key],(source,key)


def curve(recipe,source):return json.loads((DEST/'development'/f'{recipe}_{source}_42_curve.json').read_text())


def phase_one_gate(plan):
    sources=plan['phase_one_roles'][1:];choices=[]
    for epoch in range(1,plan['epochs']+1):
        a={s:curve('adapt_last_stage',s)[epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in sources}
        b={s:curve('frozen_reader',s)[epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in sources}
        choices.append({'common_epoch':epoch,'adapted_mean_neural_source_skill':float(np.mean(list(a.values()))),
           'frozen_mean_neural_source_skill':float(np.mean(list(b.values()))),'adapted_per_source':a,'frozen_per_source':b})
    selected=max(choices,key=lambda r:r['adapted_mean_neural_source_skill']);e=selected['common_epoch']-1
    a=curve('adapt_last_stage','none')[e]['comparisons'];b=curve('frozen_reader','none')[e]['comparisons']
    criteria={'mean_actual_magnitude_skill_at_least_one_percent':selected['adapted_mean_neural_source_skill']>=.01,
       'improvement_vs_matched_frozen_at_least_half_percentage_point':selected['adapted_mean_neural_source_skill']-selected['frozen_mean_neural_source_skill']>=.005,
       'no_neural_source_error_worsens_over_five_percent':all(1-selected['adapted_per_source'][s]<=1.05*(1-selected['frozen_per_source'][s]) for s in sources),
       'nearby_ordering_drop_at_most_point_zero_one':a['local']['ordering_accuracy_nonzero_measured_differences']>=b['local']['ordering_accuracy_nonzero_measured_differences']-.01,
       'within_fixed_cell_pair_rho_drop_at_most_point_zero_one':a['cross']['macro_within_fixed_study_cell_pair_Spearman']>=b['cross']['macro_within_fixed_study_cell_pair_Spearman']-.01,
       'nearby_difference_error_beats_zero_difference':a['local']['difference_MAE']<a['local']['zero_difference_MAE']}
    return {'passed':all(criteria.values()),'criteria':criteria,'selected':selected,'all_prespecified_epochs':choices,
       'held_region_pair_diagnostics':{'adapted':a,'frozen':b},'new_outside_test_used':False,
       'phase_two_permitted':all(criteria.values()),'scientific_goal_complete':False}


def main():
    assert_versions_native();plan=get_plan();assert json.loads((DEST/'native_runtime_and_gradient_audit.json').read_text())['passed'];start=time.monotonic()
    for source in plan['phase_one_roles']:
        for recipe in plan['recipes']:fit(recipe,source,plan)
    matched_roles(plan['phase_one_roles']);gate=phase_one_gate(plan);write_json(DEST/'phase_one_neural_gate.json',gate)
    if not gate['passed']:
        write_json(DEST/'native_comparison_assessment.json',{'phase_one_complete':True,'phase_two_run':False,'decision':'Stop this finite comparison; the prespecified neural generalization gate failed.',
           'gate':gate,'matched_rows_batches_initial_parameters_verified':True,'new_test_used':False,'scientific_goal_complete':False,'elapsed_seconds':time.monotonic()-start})
        assert_versions_native();print('NATIVE COMPARISON STOPPED AT PRESPECIFIED NEURAL GATE',json.dumps(gate['selected']),flush=True);return
    for source in plan['phase_two_roles_if_gate_passes']:
        for recipe in plan['recipes']:fit(recipe,source,plan)
    roles=plan['phase_one_roles']+plan['phase_two_roles_if_gate_passes'];matched_roles(roles);epoch=gate['selected']['common_epoch'];allsource={}
    for recipe in plan['recipes']:
        per={s:curve(recipe,s)[epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in roles if s!='none'}
        allsource[recipe]={'mean_seven_source_actual_magnitude_skill':float(np.mean(list(per.values()))),'per_source':per}
    write_json(DEST/'native_comparison_assessment.json',{'phase_one_complete':True,'phase_two_run':True,'common_epoch':epoch,
       'gate':gate,'all_source_comparison':allsource,'matched_rows_batches_initial_parameters_verified':True,
       'new_test_used':False,'new_confidence_or_independent_prediction_success_established':False,
       'scientific_goal_complete':False,'elapsed_seconds':time.monotonic()-start})
    assert_versions_native();print('NATIVE COMPARISON COMPLETE',json.dumps(allsource),flush=True)


if __name__=='__main__':main()
