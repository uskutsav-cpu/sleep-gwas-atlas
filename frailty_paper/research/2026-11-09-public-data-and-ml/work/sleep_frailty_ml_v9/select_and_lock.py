"""One common development choice; preserve every adverse source comparison."""
from v9_common import *

FINAL = ROOT/'outputs/sleep_frailty_ai_v9_validation'


def main():
    assert_previous_sources()
    execution = json.loads((OUT/'finite_development_execution_receipt.json').read_text())
    assert len(execution['fits']) == 16 and not execution['fresh_model_performance_read']
    plan = json.loads((OUT/'finite_development_plan.json').read_text())
    studies = [s for s in plan['held_studies'] if s != 'none']
    comparisons = []
    for recipe in plan['recipes']:
        curves = {s: json.loads((OUT/f'{recipe}_{s}_42_curve.json').read_text()) for s in studies}
        assert all(len(c) == plan['epochs'] for c in curves.values())
        for epoch in range(1, plan['epochs']+1):
            per_source = {s: c[epoch-1]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s,c in curves.items()}
            comparisons.append({'recipe':recipe, 'epoch':epoch, 'mean_source_skill':float(np.mean(list(per_source.values()))), 'per_source_skill':per_source})
    selected = max(comparisons, key=lambda x:x['mean_source_skill'])
    FINAL.mkdir(parents=True, exist_ok=True)
    result = {'selection':selected, 'comparisons':comparisons, 'development_not_replication':True,
              'interpretation':'Only a small development improvement. No claim of useful independent prediction; untouched-study and uncertainty tests remain.',
              'baseline_skill':0., 'baseline_wins_development':selected['mean_source_skill']<=0,
              'finite_plan_sha256':sha(OUT/'finite_development_plan.json'),
              'curve_sha256':{f'{r}_{s}_42_curve.json':sha(OUT/f'{r}_{s}_42_curve.json') for r in plan['recipes'] for s in studies}}
    path = FINAL/'common_epoch_selection.json'
    assert not path.exists()
    write_json(path,result)
    lock = {'recipe':selected['recipe'], 'epochs':selected['epoch'], 'seeds':[42,73,101],
            'held_studies':['none']+studies, 'source_data_sha256':sha(DATA/'training_arrays.npz'),
            'selection_sha256':sha(path), 'development_plan_sha256':sha(OUT/'finite_development_plan.json'),
            'code_sha256':{p.name:sha(p) for p in [WORK/'line_aware_network.py',WORK/'training_utilities.py']},
            'training_roles':[0,1], 'calibration_role':2,
            'fresh_model_performance_read':False,
            'stopping_rule':'No further architecture/source/epoch changes before fresh evaluation. Any later change makes this test revealed development.',
            'matched_comparison':'Sequence and AVI arms use identical rows, seeds, epochs, batches and shared initial weights; appended AVI weights start at zero.',
            'external_protocol_sha256':sha(V8/'outputs/V8_external_validation_protocol.json'),
            'baseline_amendment':'Also require comparison with the training-only culture-line median; unknown line falls back to cell median. No unknown-line release.',
            'output_location':'Internal disk to avoid EXFAT allocation overhead; prior external outputs remain unchanged.'}
    write_json(FINAL/'final_settings_lock.json',lock)
    print(json.dumps({'selected':selected,'final_output':str(FINAL)}))


if __name__ == '__main__':main()
