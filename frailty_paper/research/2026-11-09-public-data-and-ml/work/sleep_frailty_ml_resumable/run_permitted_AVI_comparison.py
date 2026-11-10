"""Run one independently gated, matched permissible-score development test.

This new driver has not been used for scientific training. Its prerequisites
include the completed original native decision and an actual new-trainer
execution/continuation/resource audit. No fresh-study outcomes are read here.
"""
import argparse
import contextlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from permitted_AVI_experiment import OUT, PLAN_PATH, get_AVI_plan, readiness, sha
from process_control import command

WORK = Path(__file__).parent
DEVELOPMENT = OUT / 'AVI_development'


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def live_fit(script, recipe, source, arm):
    rows = subprocess.run(['ps', '-axo', 'pid,command'], check=True,
                          capture_output=True, text=True).stdout.splitlines()
    matches = []
    for row in rows:
        fields = row.strip().split(None, 1)
        if len(fields) != 2:
            continue
        pid, cmd = fields
        if str(script) in cmd and all(token in ' ' + cmd + ' ' for token in
                [f' --recipe {recipe} ', f' --held-study {source} ', f' --arm {arm} ']):
            matches.append(int(pid))
    assert len(matches) <= 1, ('Duplicate matched-AVI workers', matches)
    return matches[0] if matches else None


@contextlib.contextmanager
def lock():
    path = OUT / 'permitted_AVI_driver_lock.json'
    script = str(Path(__file__).resolve())
    if path.exists():
        previous = json.loads(path.read_text())
        if script in command(previous['PID']):
            raise SystemExit(f'Existing matched-AVI driver confirmed live: {previous["PID"]}.')
        path.unlink()
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, 'w') as handle:
        json.dump({'PID': os.getpid(), 'script': script}, handle)
    try:
        yield
    finally:
        if path.exists() and json.loads(path.read_text())['PID'] == os.getpid():
            path.unlink()


def stem(plan, source, arm):
    return f'{plan["recipes"][0]}_{source}_{arm}_{plan["seed"]}'


def fit(plan, source, arm):
    name = stem(plan, source, arm)
    record = DEVELOPMENT / (name + '_record.json')
    if record.exists():
        saved = json.loads(record.read_text())
        assert saved['plan_sha256'] == sha(PLAN_PATH) and saved['completed_epochs'] == plan['epochs']
        return
    script = (WORK / 'train_permitted_AVI_development.py').resolve()
    existing = live_fit(script, plan['recipes'][0], source, arm)
    if existing:
        print('FOLLOW EXISTING MATCHED AVI FIT', name, existing, flush=True)
        while str(script) in command(existing):
            time.sleep(10)
        if record.exists():
            saved = json.loads(record.read_text())
            assert saved['plan_sha256'] == sha(PLAN_PATH) and saved['completed_epochs'] == plan['epochs']
            return
    started = (DEVELOPMENT / (name + '_started.json')).exists()
    if started:
        assert (DEVELOPMENT / (name + '_recovery.pt')).exists(), 'Full state is required for explicit resumption.'
    print('START MATCHED AVI INPUT FIT', name, flush=True)
    with (DEVELOPMENT / (name + '.log')).open('a' if started else 'w') as handle:
        subprocess.run([sys.executable, '-u', str(script), '--recipe', plan['recipes'][0],
                        '--held-study', source, '--arm', arm] + (['--resume'] if started else []),
                       check=True, stdout=handle, stderr=subprocess.STDOUT,
                       env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'OMP_NUM_THREADS': '2',
                            'OPENBLAS_NUM_THREADS': '2', 'VECLIB_MAXIMUM_THREADS': '2',
                            'PYTORCH_MPS_HIGH_WATERMARK_RATIO': '0.4', 'PYTORCH_MPS_LOW_WATERMARK_RATIO': '0.3'})
    assert record.exists()
    print('DONE MATCHED AVI INPUT FIT', name, flush=True)


def verify_matched_roles(plan, roles):
    common = ['recipe', 'seed', 'epochs', 'training_index_sha256', 'prediction_index_sha256',
              'initial_parameter_sha256', 'normalization', 'cell_medians', 'line_medians',
              'line_parent', 'supported_cells', 'supported_lines', 'context_transformations',
              'AVI_scaling', 'AVI_bank_sha256', 'readonly_input_receipt', 'plan_sha256',
              'code_sha256', 'actual_batch_schedule_sha256', 'actual_backprop_updates',
              'completed_epochs', 'new_original_Siraj_signed_training_rows',
              'new_original_neuronal_signed_training_rows', 'total_known_signed_training_rows',
              'trainable_head_parameters', 'trainable_reader_parameters']
    for source in roles:
        records = [json.loads((DEVELOPMENT / (stem(plan, source, arm) + '_record.json')).read_text())
                   for arm in plan['arms']]
        assert [record['arm'] for record in records] == ['scoreless', 'scores']
        for key in common:
            assert records[0][key] == records[1][key], (source, key)


def final_metrics(plan, source, arm):
    curve = json.loads((DEVELOPMENT / (stem(plan, source, arm) + '_curve.json')).read_text())
    assert len(curve) == plan['epochs']
    assert [point['epoch'] for point in curve] == list(range(1, plan['epochs'] + 1))
    return curve[-1]  # Fixed before score fitting; no best-AVI-epoch selection.


def development_gate(plan):
    gate = plan['development_gate']
    sources = plan['phase_one_roles'][1:]
    assert len(sources) == 3
    values = {arm: {source: final_metrics(plan, source, arm)['macro_MAE_skill_vs_line_baseline_supported_cells']
                    for source in sources} for arm in plan['arms']}
    assert all(value is not None and np.isfinite(value) for arm in values.values() for value in arm.values())
    means = {arm: float(np.mean(list(per.values()))) for arm, per in values.items()}
    pairs = {arm: final_metrics(plan, 'none', arm)['comparisons'] for arm in plan['arms']}
    a, b = pairs['scores'], pairs['scoreless']
    criteria = {
        'score_arm_mean_neural_source_magnitude_skill': means['scores'] >= gate['mean_three_neural_source_score_arm_skill_at_least'],
        'gain_vs_exact_matched_scoreless_control': means['scores'] - means['scoreless'] >= gate['mean_neural_source_skill_gain_vs_matched_scoreless_at_least'],
        'individual_neural_source_error_not_excessively_worse': all(
            1 - values['scores'][s] <= gate['individual_neural_source_error_ratio_vs_matched_scoreless_at_most'] *
            (1 - values['scoreless'][s]) for s in sources),
        'nearby_ordering_preserved': a['local']['ordering_accuracy_nonzero_measured_differences'] >=
            b['local']['ordering_accuracy_nonzero_measured_differences'] - gate['held_region_nearby_ordering_drop_at_most'],
        'fixed_cell_pair_ranking_preserved': a['cross']['macro_within_fixed_study_cell_pair_Spearman'] >=
            b['cross']['macro_within_fixed_study_cell_pair_Spearman'] - gate['held_region_fixed_cell_pair_rho_drop_at_most'],
        'nearby_difference_error_beats_zero_difference': a['local']['difference_MAE'] < a['local']['zero_difference_MAE']}
    return {'passed': all(criteria.values()), 'criteria': criteria, 'fixed_epoch': plan['epochs'],
            'per_source_skill': values, 'mean_three_neural_source_skill': means,
            'held_region_pair_diagnostics': pairs, 'best_AVI_epoch_selected': False,
            'outside_study_outcomes_read': False, 'keep_AVI_for_project_predictions': False,
            'interpretation': 'Development screening only. Even a pass requires all-source checks, useful confidence and the locked untouched-study benefit test.'}


def all_source_candidate_gate(plan, per_source, phase_one):
    sources = plan['phase_one_roles'][1:] + plan['phase_two_roles_if_gate_passes']
    assert len(sources) == len(set(sources)) == 7
    assert set(per_source) == set(plan['arms']) == {'scoreless','scores'}
    for arm in plan['arms']:
        assert set(per_source[arm]) == set(sources), ('Incomplete whole-source comparison',arm)
        assert all(value is not None and math.isfinite(value) and value <= 1+1e-12
                   for value in per_source[arm].values()),arm
    thresholds = plan['development_gate']['all_source_candidate_gate']
    means = {arm:math.fsum(per_source[arm][source] for source in sources)/len(sources)
             for arm in plan['arms']}
    criteria = {
        'earlier_neural_and_pair_gate_passed': bool(phase_one['passed']),
        'mean_seven_source_magnitude_skill':means['scores'] >= thresholds['mean_seven_source_score_arm_skill_at_least'],
        'mean_seven_source_gain_vs_matched_control':means['scores']-means['scoreless'] >= thresholds['mean_seven_source_gain_vs_matched_control_at_least'],
        'no_source_error_worsens_over_declared_limit':all(
            1-per_source['scores'][s] <= thresholds['individual_source_error_ratio_vs_matched_control_at_most']*(1-per_source['scoreless'][s]) for s in sources)}
    passed = all(criteria.values())
    return {'passed':passed,'criteria':criteria,'mean_seven_source_skill':means,
            'per_source_skill':per_source,'fresh_study_prediction_lock_permitted_by_development':passed,
            'outside_study_outcomes_read':False,'keep_AVI_for_project_predictions':False,
            'interpretation':'Only a full seven-source candidate pass permits preparation of the separate fresh-study test. Useful confidence and external benefit remain unproven.'}


def main():
    plan = get_AVI_plan()
    DEVELOPMENT.mkdir(exist_ok=True)
    start = time.monotonic()
    for source in plan['phase_one_roles']:
        for arm in plan['arms']:
            fit(plan, source, arm)
    verify_matched_roles(plan, plan['phase_one_roles'])
    gate = development_gate(plan)
    write(DEVELOPMENT / 'matched_AVI_phase_one_gate.json', gate)
    if not gate['passed']:
        result = {'phase_one_complete': True, 'phase_two_run': False, 'gate': gate,
                  'fresh_study_candidate_qualified': False,
                  'decision': 'Stop this finite added-input comparison at its declared development gate.',
                  'matched_training_verified': True, 'keep_AVI_for_project_predictions': False,
                  'outside_study_outcomes_read': False, 'scientific_goal_complete': False,
                  'elapsed_seconds': time.monotonic() - start}
        write(DEVELOPMENT / 'matched_AVI_development_assessment.json', result)
        print('MATCHED AVI DEVELOPMENT STOPPED AT FIXED GATE', json.dumps(gate['criteria']), flush=True)
        return
    for source in plan['phase_two_roles_if_gate_passes']:
        for arm in plan['arms']:
            fit(plan, source, arm)
    roles = plan['phase_one_roles'] + plan['phase_two_roles_if_gate_passes']
    verify_matched_roles(plan, roles)
    per = {arm: {source: final_metrics(plan, source, arm)['macro_MAE_skill_vs_line_baseline_supported_cells']
                 for source in roles if source != 'none'} for arm in plan['arms']}
    candidate = all_source_candidate_gate(plan, per, gate)
    result = {'phase_one_complete': True, 'phase_two_run': True, 'gate': gate, 'fixed_epoch': plan['epochs'],
              'final_all_source_candidate_gate': candidate,
              'fresh_study_candidate_qualified': candidate['passed'],
              'per_source_skill': per, 'mean_seven_source_skill': {arm: float(np.mean(list(v.values()))) for arm, v in per.items()},
              'matched_training_verified': True, 'outside_study_outcomes_read': False,
              'keep_AVI_for_project_predictions': False, 'scientific_goal_complete': False,
              'remaining': ('Independently audit results, fix confidence or explicit no-confidence release, and freeze both arms and final-study predictions before any untouched-study outcome comparison.' if candidate['passed'] else 'Stop this candidate: the final seven-source gate failed. Preserve all results and leave the fresh study unopened.'),
              'elapsed_seconds': time.monotonic() - start}
    write(DEVELOPMENT / 'matched_AVI_development_assessment.json', result)
    print('MATCHED AVI DEVELOPMENT COMPLETE; NO INPUT PROMOTION', 'seven_source_gate', candidate['passed'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-readiness', action='store_true')
    args = parser.parse_args()
    if args.check_readiness:
        print(json.dumps(readiness()), flush=True)
    else:
        with lock():
            main()
