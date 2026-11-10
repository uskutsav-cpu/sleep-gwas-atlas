"""Select one common neural candidate, retaining the simple baseline result."""
from v8_common import *


def main():
    assert (OUT/'activity_round_execution_receipt.json').exists()
    studies = ['Kosicki2025', 'IGVF240K', 'Salomon', 'Siraj']
    candidates = []
    for family, recipes in [('initial', ['frozen', 'pooled_tune']), ('activity', ['activity', 'activity_precision', 'local_activity_precision'])]:
        for recipe in recipes:
            curves = []
            for study in studies:
                name = f'{recipe}_{study}_42' if family == 'initial' else f'round2_{recipe}_{study}_development_sequence_42'
                curves.append(json.loads((OUT/(name+'_curve.json')).read_text()))
            assert len({len(c) for c in curves}) == 1
            values = np.array([[r['macro_MAE_skill_supported_cells'] for r in c] for c in curves])
            for e in range(values.shape[1]):
                candidates.append({'family': family, 'recipe': recipe, 'epoch': e+1, 'mean_whole_study_MAE_skill': float(values[:, e].mean()), 'per_study_MAE_skill': dict(zip(studies, values[:, e].tolist()))})
    best = max(candidates, key=lambda r: r['mean_whole_study_MAE_skill'])
    result = {
        'selected_family': best['family'], 'selected_recipe': best['recipe'], 'selected_epoch': best['epoch'],
        'selected_neural_candidate': best, 'all_common_epoch_candidates': candidates,
        'simple_baseline_MAE_skill': 0.,
        'neural_candidate_improves_mean_development_transfer': best['mean_whole_study_MAE_skill'] > 0,
        'currently_recommended_by_development_error': 'neural_candidate' if best['mean_whole_study_MAE_skill'] > 0 else 'training_cell_median_baseline',
        'interpretation': 'This chooses the best tested NEURAL research candidate for a fully declared external experiment, not an adopted useful model. The simple median remains the performance leader if all neural candidates have negative transfer skill. No post hoc favorable-study or epoch selection is presented as independent validation.',
        'untried_end_to_end_source_omission_results_used': False,
        'fresh_study_performance_used': False,
    }
    write_json(OUT/'activity_round_model_selection.json', result)
    lock = {
        'selected_family': best['family'], 'selected_recipe': best['recipe'], 'selected_epoch': best['epoch'],
        'selection_sha256': sha(OUT/'activity_round_model_selection.json'),
        'external_protocol_sha256': sha(OUT/'V8_external_validation_protocol.json'),
        'seeds': [42, 73, 101], 'matched_input_arms': ['sequence', 'AVI'],
        'source_omission_studies': ['Salomon', 'Kosicki2025', 'Siraj', 'IGVF240K', 'Koesterich2023'],
        'final_training_roles': [0, 1], 'calibration_role': 2,
        'fresh_performance_calculated': False,
    }
    path = OUT/'activity_round_final_settings_lock.json'
    if path.exists():
        assert json.loads(path.read_text()) == lock
    else:
        write_json(path, lock)
    print('NEURAL CANDIDATE LOCKED', best, 'median_baseline_still_leads', best['mean_whole_study_MAE_skill'] <= 0, flush=True)


if __name__ == '__main__':
    main()
