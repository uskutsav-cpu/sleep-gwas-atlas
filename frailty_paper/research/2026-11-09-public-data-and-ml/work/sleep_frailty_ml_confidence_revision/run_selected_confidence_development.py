"""Finite matched calibration-method comparison on revealed development only."""
import datetime
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'work/sleep_frailty_ml_next'))
from calibrate_coverage_confidence import fit_error_package, uncertainty
from coverage_validation_common import load_ensemble, confidence_inputs
from coverage_inputs import load_inputs
from v9_common import np, json, hashlib, sha
sys.path.insert(0, str(Path(__file__).resolve().parent))
from selected_calibration import calibrate_selected, apply_selected

OUT = ROOT / 'outputs/sleep_frailty_ml_confidence_revision'
OLD = ROOT / 'outputs/sleep_frailty_ml_next/validation'


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def metrics(y, p, baseline, width, accepted, groups):
    result = {'rows': len(y), 'accepted_rows': int(accepted.sum()),
              'accepted_fraction': float(accepted.mean()),
              'accepted_regions': len(np.unique(groups[accepted]))}
    if not accepted.any():
        result['useful_answers'] = False
        return result
    actual, pred, b, g = y[accepted], p[accepted], baseline[accepted], groups[accepted]
    radius = width[accepted]
    _, inverse, counts = np.unique(g, return_inverse=True, return_counts=True)
    weights = 1. / counts[inverse]
    error, be = abs(actual - pred), abs(actual - b)
    covered = abs(actual - pred) <= radius
    mae, baseline_mae = float(error.mean()), float(be.mean())
    result.update(raw_observed_MAE=mae, baseline_MAE=baseline_mae,
                  raw_MAE_skill=1 - mae / baseline_mae if baseline_mae > 0 else None,
                  row_coverage=float(covered.mean()), region_balanced_coverage=float(np.average(covered, weights=weights)),
                  median_halfwidth=float(np.median(radius)),
                  median_baseline_magnitude=float(np.median(b)),
                  fraction_radius_at_most_original_width_cap=float(np.mean(radius <= np.maximum(.2, b))),
                  useful_answers=False)
    return result


def main():
    OUT.mkdir(exist_ok=True)
    previous = json.loads((OUT / 'archived_confidence_bottleneck_diagnostic.json').read_text())
    assert previous['passed_identity_and_policy_reconstruction']
    policy_path = OLD / 'confidence_policy_lock.json'
    policy = json.loads(policy_path.read_text())
    settings = json.loads((OLD / 'final_settings_lock.json').read_text())
    old_assessment = json.loads((OLD / 'sequence_confidence_assessment.json').read_text())
    plan = {
        'purpose': 'One finite development experiment isolating interval calibration/answer selection. No change to the point predictor.',
        'archived_base_model': str((OLD / 'final_settings_lock.json').relative_to(ROOT)),
        'base_model_settings_sha256': sha(OLD / 'final_settings_lock.json'),
        'control_confidence_policy_sha256': sha(policy_path),
        'outer_studies': settings['held_studies'][1:],
        'shared_error_critics': policy['critic_parameters'],
        'matched_error_inputs_and_fit_scale_partition': 'Unchanged original outer-source-omitted base predictions. Other sources original calibration groups only; hash even error fit, odd scale. Outer source absent from all fitting stages.',
        'candidate_selection': 'Known exact cell and line, original support/OOD and estimated 5% plus 0.005 error-gain criteria. Compare raw learned 90th-error radius with the original max(0.2, typical effect) cap BEFORE looking at scale outcomes.',
        'candidate_scale': 'Calibrate selected residual/raw-radius ratios within the same cell; max pooled and adequately represented source 90th region-balanced quantile. No unconditional absolute-error floor. Require original 100 rows/20 regions; sparse cell yields no answer rather than different-cell fallback.',
        'after_calibration': 'Do not filter by the resulting interval width. Report oversized ranges as failure of usefulness; do not silently change the calibrated selected population.',
        'control': 'Exactly reproduce the archived policy and nested outer-source summaries using the same newly fitted critics before interpreting candidate results.',
        'scope': 'Previously revealed development sources. This is not new independent validation and does not make the point predictor stronger.',
        'readiness_only': {'minimum_qualifying_sources': 2, 'minimum_accepted_rows': 30,
             'minimum_accepted_regions': 20, 'minimum_accepted_fraction': .01,
             'minimum_row_and_region_balanced_coverage': .9, 'minimum_selected_MAE_skill': .05,
             'minimum_fraction_accepted_intervals_within_original_cap': .9},
        'coverage_theory': 'Empirical only. Selected/group-dependent observations and changing laboratories are not asserted exchangeable; no general conformal guarantee.',
        'stopping': 'Exactly these two policies. No threshold, critic, source, epoch or quantile search after outcomes.',
        'references': ['https://academic.oup.com/jrsssb/article/87/4/1239/8113856',
                       'https://academic.oup.com/jrsssb/article/87/4/1100/8058684'],
        'new_Ulirsch_outcomes_read': False,
        'code_sha256': {p.name: sha(p) for p in [Path(__file__), Path(__file__).with_name('selected_calibration.py')]}}
    lock = OUT / 'selected_confidence_development_plan.json'
    if lock.exists():
        assert json.loads(lock.read_text()) == plan
    else:
        write(lock, plan)
    start = time.monotonic()
    d, pairs, bank, x = load_inputs()
    results = {}
    for source in plan['outer_studies']:
        target = OUT / f'selected_confidence_{source}_development.json'
        if target.exists():
            existing = json.loads(target.read_text())
            assert existing['plan_sha256'] == sha(lock)
            results[source] = existing
            continue
        meta, models, pred = load_ensemble(source, 'sequence')
        del models
        ids = pred['observation_index']
        features, baseline, _, cs, ls, counts = confidence_inputs(pred, x[d['bank_index'][ids]], d['cell'][ids], d['line'][ids], meta)
        other_meta, models, other = load_ensemble(source, 'sequence', True)
        del models
        oi = other['observation_index']
        assert np.all(d['split'][oi] == 2)
        sid = d['study_names'].tolist().index(source)
        assert not np.any(d['study'][oi] == sid)
        assert not np.intersect1d(d['group'][oi], d['group'][d['study'] == sid]).size
        of, ob, _, ocs, ols, oc = confidence_inputs(other, x[d['bank_index'][oi]], d['cell'][oi], d['line'][oi], other_meta)
        oi, of, ob, op, ols, oc = oi[ocs], of[ocs], ob[ocs], other['magnitude'][ocs], ols[ocs], oc[ocs]
        fit = np.array([int(hashlib.sha256(('error-partition-20261108|' + str(g)).encode()).hexdigest()[:16], 16) % 2 == 0 for g in d['group'][oi]])
        scale = ~fit
        package = fit_error_package(of, d['effect'][oi], op, ob, d['cell'][oi], d['study'][oi], d['group'][oi], fit, scale, policy)
        control = uncertainty(package, features, d['cell'][ids], baseline, cs, ls, counts, policy)
        cm = metrics(d['effect'][ids], pred['magnitude'], baseline, control['interval_halfwidth'], control['pretest_policy_accepted'], d['group'][ids])
        reported = old_assessment['nested_outer_source_development'][source]
        assert cm['accepted_rows'] == reported['accepted_rows'], (source, cm, reported)
        assert abs(float(np.mean(abs(d['effect'][ids] - pred['magnitude']) <= control['interval_halfwidth'])) - reported['coverage']) < 1e-8
        if cm['accepted_rows']:
            assert abs(cm['raw_observed_MAE'] - reported['accepted']['MAE']) < 1e-6
            assert abs(cm['row_coverage'] - reported['accepted']['coverage']) < 1e-8
        calibrated = calibrate_selected(package, of[scale], d['effect'][oi[scale]], op[scale], ob[scale], d['cell'][oi[scale]],
              d['study'][oi[scale]], d['group'][oi[scale]], np.ones(int(scale.sum()), dtype=bool), ols[scale], oc[scale], policy)
        candidate = apply_selected(package, calibrated, features, baseline, d['cell'][ids], cs, ls, counts, policy)
        sm = metrics(d['effect'][ids], pred['magnitude'], baseline, candidate['interval_halfwidth'], candidate['accepted'], d['group'][ids])
        r = plan['readiness_only']
        criteria = {'enough_answers': sm['accepted_rows'] >= r['minimum_accepted_rows'],
            'enough_regions': sm['accepted_regions'] >= r['minimum_accepted_regions'],
            'enough_answer_fraction': sm['accepted_fraction'] >= r['minimum_accepted_fraction'],
            'row_coverage': sm.get('row_coverage', 0) >= r['minimum_row_and_region_balanced_coverage'],
            'region_coverage': sm.get('region_balanced_coverage', 0) >= r['minimum_row_and_region_balanced_coverage'],
            'selected_prediction_better_than_baseline': (sm.get('raw_MAE_skill') or 0) >= r['minimum_selected_MAE_skill'],
            'useful_widths': sm.get('fraction_radius_at_most_original_width_cap', 0) >= r['minimum_fraction_accepted_intervals_within_original_cap']}
        result = {'source': source, 'plan_sha256': sha(lock), 'archived_control_reproduced': True,
            'control': cm, 'selected_calibration_candidate': sm, 'candidate_readiness_criteria': criteria,
            'candidate_source_qualifies': all(criteria.values()), 'conditional_scale': calibrated,
            'outer_source_absent_from_all_fitting_stages': True,
            'error_fit_index_sha256': hashlib.sha256(oi[fit].tobytes()).hexdigest(),
            'scale_index_sha256': hashlib.sha256(oi[scale].tobytes()).hexdigest(),
            'new_study_outcomes_read': False, 'independent_validation_success_established': False}
        write(target, result)
        np.savez_compressed(OUT / f'selected_confidence_{source}_decisions.npz', observation_index=ids,
            control_accepted=control['pretest_policy_accepted'], candidate_eligible=candidate['precalibration_eligible'],
            candidate_accepted=candidate['accepted'], candidate_interval_halfwidth=candidate['interval_halfwidth'])
        results[source] = result
        print('SELECTED CONFIDENCE DEVELOPMENT SOURCE', source, json.dumps(sm), flush=True)
    qualified = [source for source, r in results.items() if r['candidate_source_qualifies']]
    write(OUT / 'selected_confidence_development_assessment.json', {'comparison_complete': True,
        'all_archived_control_results_reproduced': True, 'qualifying_sources': qualified,
        'candidate_passes_development_readiness_only': len(qualified) >= plan['readiness_only']['minimum_qualifying_sources'],
        'sources': results, 'plan_sha256': sha(lock), 'new_independent_study_used': False,
        'current_native_training_changed': False, 'new_point_prediction_or_confidence_success_established': False,
        'scientific_goal_complete': False, 'elapsed_seconds': time.monotonic() - start})
    print('FINITE SELECTED CONFIDENCE COMPARISON COMPLETE', qualified, flush=True)


if __name__ == '__main__':
    main()
