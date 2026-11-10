"""Freeze one matched score-input development experiment after the native test.

The candidate rule is declared before any new AVI training or reserved-study
evaluation. This file neither resumes the interrupted native run nor opens the
reserved measurements.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'
AVI = OUT / 'AVI_preparation'
PLAN_PATH = AVI / 'matched_AVI_development_plan.json'
SYSTEM_RESERVE = 768 * 1024**2


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def check_audit(path):
    record = json.loads(path.read_text())
    assert record['passed'], str(path)
    hashes = record.get('code_sha256', {})
    if isinstance(hashes, str):
        hashes = {'audit_corrected_permitted_AVI.py': hashes}
    for filename, digest in hashes.items():
        candidates = [ROOT / 'work/sleep_frailty_ml_resumable' / filename,
                      ROOT / 'work/sleep_frailty_ml_next' / filename]
        source = next((p for p in candidates if p.exists()), None)
        assert source and sha(source) == digest, (path.name, filename)
    return record


def readiness():
    free = shutil.disk_usage(ROOT).free
    assessment_path = OUT / 'native_comparison_assessment.json'
    assessment = json.loads(assessment_path.read_text()) if assessment_path.exists() else None
    blockers = []
    if free < SYSTEM_RESERVE:
        blockers.append('The system volume has less than the recorded 768 MiB launch reserve; preserve partial work and recover storage before admission.')
    if not assessment or not assessment.get('phase_one_complete'):
        blockers.append('The original native comparison has not reached its prespecified scientific decision.')
    audits = {}
    for name in ['corrected_permitted_AVI_independent_audit.json',
                 'permitted_AVI_train_only_scaling_audit.json',
                 'permitted_AVI_expanded_head_audit.json',
                 'permitted_AVI_native_runtime_audit.json',
                 'readonly_native_input_equivalence_audit.json']:
        path = AVI / name
        if not path.exists():
            blockers.append('Missing numerical audit: ' + name)
            continue
        check_audit(path)
        audits[name] = sha(path)
    external_path = AVI / 'fresh_study_benefit_protocol.json'
    external_hash = None
    if not external_path.exists():
        blockers.append('The fresh-study input-benefit criteria have not been locked.')
    else:
        external = json.loads(external_path.read_text())
        assert external['reserved_outcome_values_read_by_this_script'] is False
        assert external['new_model_predictions_or_fits_made'] is False
        assert sha(Path(__file__).with_name('freeze_permitted_AVI_external_protocol.py')) == external['code_sha256']
        for relative_path, digest in external['evidence_sha256'].items():
            assert sha(ROOT / relative_path) == digest, relative_path
        external_hash = sha(external_path)
    statistics_path = Path(__file__).with_name('external_magnitude_statistics.py')
    statistics_proof_path = AVI / 'external_magnitude_statistics_synthetic_checks.json'
    statistics_hash = statistics_proof_hash = None
    if not statistics_path.exists() or not statistics_proof_path.exists():
        blockers.append('The external paired-effect statistical routine still needs its synthetic checks.')
    else:
        statistics_proof = json.loads(statistics_proof_path.read_text())
        assert statistics_proof['passed'] and statistics_proof['synthetic_only']
        assert statistics_proof['real_experimental_outcomes_read'] is False
        assert statistics_proof['code_sha256'] == sha(statistics_path)
        if external_hash is not None:
            assert statistics_proof['external_protocol_sha256'] == external_hash
        statistics_hash, statistics_proof_hash = sha(statistics_path), sha(statistics_proof_path)
    gate_proof_path = AVI / 'seven_source_candidate_gate_synthetic_checks.json'
    gate_source_path = Path(__file__).with_name('run_permitted_AVI_comparison.py')
    gate_proof_hash = None
    if not gate_proof_path.exists():
        blockers.append('The final seven-source development decision still needs its synthetic checks.')
    else:
        gate_proof = json.loads(gate_proof_path.read_text())
        assert gate_proof['passed'] and gate_proof['synthetic_only']
        assert gate_proof['real_model_or_reserved_outcomes_read'] is False
        assert gate_proof['source_sha256'] == sha(gate_source_path)
        gate_proof_hash = sha(gate_proof_path)
    trainer_audit = AVI / 'matched_AVI_trainer_execution_audit.json'
    if not trainer_audit.exists():
        blockers.append('The new matched trainer still needs its execution and full-state continuation audit; head-only numerical checks do not substitute for that.')
    else:
        record = check_audit(trainer_audit)
        for requirement in ['both_arms_executed', 'common_initialization_and_batch_schedule_verified',
                            'actual_saved_state_continuation_verified', 'fitting_device_resource_preflight_passed']:
            assert record.get(requirement) is True, requirement
        audits[trainer_audit.name] = sha(trainer_audit)
    resource_path = AVI / 'matched_AVI_resource_admission_revision.json'
    resource_hash = None
    if not resource_path.exists():
        blockers.append('The revised resource admission needs its source derivation and measured native-runtime evidence.')
    else:
        resource_proof = check_audit(resource_path)
        assert resource_proof['launch_reserve_bytes'] == SYSTEM_RESERVE
        assert resource_proof['live_saved_state_reserve_bytes'] == 300 * 1024**2
        for relative_path, digest in resource_proof['evidence_sha256'].items():
            assert sha(ROOT / relative_path) == digest, relative_path
        resource_hash = sha(resource_path)
    return {'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'scientific_training_ready': not blockers, 'blockers': blockers,
            'free_system_volume_bytes': free, 'launch_reserve_bytes': SYSTEM_RESERVE,
            'resource_admission_revision_sha256': resource_hash,
            'current_numerical_audits': audits,
            'external_benefit_protocol_sha256': external_hash,
            'external_statistics_source_sha256': statistics_hash,
            'external_statistics_synthetic_check_sha256': statistics_proof_hash,
            'seven_source_candidate_gate_synthetic_check_sha256': gate_proof_hash,
            'candidate_rule': 'If the original native gate passes, use its declared adapted candidate and common epoch. Otherwise use the original frozen reader at exactly three epochs for a separate input-benefit question. Do not extend or rescue the failed adaptation comparison.',
            'reserved_outcomes_read': False, 'scientific_training_launched': False}


def get_AVI_plan(*, for_execution_audit=False):
    status = readiness()
    blockers = status['blockers']
    if for_execution_audit:
        # Only bounded diagnostic runs in their separate audit directory may
        # execute before their own execution audit exists. The storage/native
        # decision/data/model audits still apply. The trainer requires an
        # explicit positive update limit for every such run.
        blockers = [message for message in blockers if not message.startswith('The new matched trainer still needs')]
    if blockers:
        raise RuntimeError('New AVI scientific fit not started: ' + ' '.join(blockers))
    assessment_path = OUT / 'native_comparison_assessment.json'
    assessment = json.loads(assessment_path.read_text())
    original_path = OUT / 'finite_stable_native_plan.json'
    plan = json.loads(original_path.read_text())
    adapted = bool(assessment['gate']['passed'])
    recipe = 'adapt_last_stage' if adapted else 'frozen_reader'
    epochs = assessment['gate']['selected']['common_epoch'] if adapted else 3
    plan.update({
        'experiment': 'One separately frozen matched permissible composite-score input comparison',
        'original_native_plan_sha256': sha(original_path),
        'original_native_assessment_sha256': sha(assessment_path),
        'recipes': [recipe], 'arms': ['scoreless', 'scores'], 'epochs': epochs,
        'candidate_rule_declared_before_new_score_training': status['candidate_rule'],
        'candidate_interpretation': 'Development candidate for the added-input question; the frozen fallback is not claimed to be a useful or independently validated sequence predictor.',
        'permitted_score_bank_sha256': sha(AVI / 'corrected_permitted_AVI_physical_bank.npy'),
        'score_inputs': 'Train-only scaled raw/PHRED composite AVI, plus identical availability information in both arms. Unknown inputs are neutral model features with a flag, not zero experimental effects.',
        'score_scope': 'Only the designated downloadable permissive composite AVI SNV artifact. No ordinary API/Atlas tracks, embeddings, attributions, pseudolabels, or protein predictions enter training. Composite benefits cannot be assigned uniquely to AlphaGenome.',
        'matched_model': 'Same expanded native head in both arms: 239 molecular/physical/score features. Preserve original common initial weights; zero-initialize all three added columns in both arms. The scoreless arm zeros only raw/PHRED values; availability is identical.',
        'data_loader': 'Separately audited read-only batch indexing of the same immutable banks; every used row and the full training normalization matched exactly. This AVI experiment has separate identities and does not modify the original native source or its preserved checkpoints; native allocation-only continuation is recorded in its own execution revision.',
        'optimization_and_supervision': 'Unchanged actual observed magnitude absolute error, original signed/activity auxiliaries, local/cross-cell pair losses, seed, batches, update schedule, optimizer settings, native lengths and sequence policy from the original plan. Calibration rows stay excluded.',
        'epoch_selection': 'Fix epochs from the candidate rule before any AVI fit. No score-guided epoch, architecture, loss, source, threshold or learning-rate search.',
        'development_gate': {
            'phase_one_roles': plan['phase_one_roles'],
            'mean_three_neural_source_score_arm_skill_at_least': .01,
            'mean_neural_source_skill_gain_vs_matched_scoreless_at_least': .005,
            'individual_neural_source_error_ratio_vs_matched_scoreless_at_most': 1.05,
            'held_region_nearby_ordering_drop_at_most': .01,
            'held_region_fixed_cell_pair_rho_drop_at_most': .01,
            'nearby_difference_error_must_beat_zero_difference': True,
            'other_four_source_roles_only_if_gate_passes': plan['phase_two_roles_if_gate_passes'],
            'all_source_candidate_gate': {
                'mean_seven_source_score_arm_skill_at_least': .01,
                'mean_seven_source_gain_vs_matched_control_at_least': .005,
                'individual_source_error_ratio_vs_matched_control_at_most': 1.05,
                'phase_one_gate_must_also_pass': True,
                'failure_action': 'Stop this added-input candidate before reserved-study prediction or outcome joining; preserve all seven source results.'},
            'stop_rule': 'Fail any phase-one criterion: preserve the negative result, stop this finite AVI development comparison, and do not inspect the untouched study to rescue it.'},
        'external_validation': 'Passing development alone does not retain the score. Fix both fitted arms, confidence and acceptance rules, matched eligible inputs and predictions before joining the reserved study outcomes. Retain only after a separately locked real benefit test on that study; otherwise reject. Do not claim current numerical audits establish such benefit.',
        'external_benefit_protocol_sha256': status['external_benefit_protocol_sha256'],
        'external_statistics_source_sha256': status['external_statistics_source_sha256'],
        'external_statistics_synthetic_check_sha256': status['external_statistics_synthetic_check_sha256'],
        'finite_comparison_driver_sha256': sha(Path(__file__).with_name('run_permitted_AVI_comparison.py')),
        'seven_source_candidate_gate_synthetic_check_sha256': status['seven_source_candidate_gate_synthetic_check_sha256'],
        'engineering_requirements_before_scientific_launch': 'Audit the derived trainer and matched initialization, updates and full-state continuation with real training examples; this input/model audit does not by itself prove the new trainer is ready.',
        'numerical_audits_sha256': {key: value for key, value in status['current_numerical_audits'].items()
                                  if key != 'matched_AVI_trainer_execution_audit.json'},
        'launch_system_reserve_bytes': SYSTEM_RESERVE,
        'resource_admission_revision_sha256': status['resource_admission_revision_sha256'],
        'recipe_and_epochs_not_a_rescue_of_original_adaptation_experiment': True,
        'new_neuron_activity_package_inserted': False,
        'reserved_Ulirsch_outcomes_used': False, 'scientific_goal_complete': False,
        'plan_builder_sha256': sha(Path(__file__))})
    if PLAN_PATH.exists():
        assert json.loads(PLAN_PATH.read_text()) == plan, 'A frozen AVI plan cannot be silently replaced.'
    else:
        PLAN_PATH.write_text(json.dumps(plan, indent=2) + '\n')
    return plan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze-plan', action='store_true')
    args = parser.parse_args()
    if args.freeze_plan:
        plan = get_AVI_plan()
        print(json.dumps({'plan_frozen': True, 'recipe': plan['recipes'][0], 'epochs': plan['epochs']}))
    else:
        result = readiness()
        (AVI / 'matched_AVI_preparation_status.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result))


if __name__ == '__main__':
    main()
