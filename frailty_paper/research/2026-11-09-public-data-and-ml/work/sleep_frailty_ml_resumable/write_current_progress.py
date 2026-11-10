"""Snapshot actual processes, audited partial evidence and remaining science."""
import datetime
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def read(path):
    return json.loads(path.read_text()) if path.exists() else None


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    raw = subprocess.run(['ps', '-axo', 'pid,etime,rss,command'], check=True, capture_output=True, text=True).stdout
    scripts = [str(ROOT / 'work/sleep_frailty_ml_resumable' / name)
               for name in ['run_resumable_comparison.py', 'train_resumable_development.py',
                            'readonly/run_resumable_comparison.py', 'readonly_native_runtime.py']]
    scripts += [str(ROOT / 'work/sleep_frailty_ml_resumable' / name)
                for name in ['readonly_768/run_resumable_comparison.py', 'readonly_native_runtime_768.py']]
    live = []
    engineering_live = []
    engineering_script = 'bench_readonly_native_training.py'
    for line in raw.splitlines():
        values = line.strip().split(None, 3)
        if len(values) == 4 and any(s in values[3] for s in scripts) and 'python' in values[3] and '-c ' not in values[3]:
            target = engineering_live if ' --audit-tag ' in values[3] else live
            target.append({'PID': int(values[0]), 'elapsed': values[1], 'RSS_KiB': int(values[2]), 'command': values[3]})
        if len(values) == 4:
            arguments = shlex.split(values[3])
            if arguments and Path(arguments[0]).name.startswith('python') and any(Path(arg).name == engineering_script for arg in arguments[1:]):
                engineering_live.append({'PID': int(values[0]), 'elapsed': values[1], 'RSS_KiB': int(values[2]), 'command': values[3]})
    development = OUT / 'development'
    live_pids = {row['PID'] for row in live}
    progress = []
    for path in development.glob('*_live_progress.json'):
        row = read(path)
        if row['PID'] in live_pids:
            progress.append(row)
    epochs = []
    for path in sorted(development.glob('*_epoch*_receipt.json')):
        stem = path.name.removesuffix('_receipt.json')
        receipt = read(path)
        checkpoint = path.with_name(stem + '.safetensors')
        predictions = path.with_name(stem + '_predictions.npz')
        assert sha(checkpoint) == receipt['checkpoint_sha256']
        assert sha(predictions) == receipt['prediction_sha256']
        metric = read(OUT / (stem + '_independent_saved_metric_audit.json'))
        network = read(OUT / (stem + '_independent_audit.json'))
        paired = read(OUT / (stem + '_independent_pair_audit.json'))
        if metric:
            assert metric['checkpoint_sha256'] == receipt['checkpoint_sha256']
        if network:
            assert network['checkpoint_sha256'] == receipt['checkpoint_sha256']
        if paired:
            assert paired['prediction_sha256'] == receipt['prediction_sha256']
        epochs.append({'checkpoint': stem, 'saved_artifact_hashes_verified': True,
                       'independent_metric_audit_passed': bool(metric and metric['passed']),
                       'independent_forward_audit_passed': bool(network and network['passed']),
                       'independent_measured_pair_audit_passed': bool(paired and paired['passed']),
                       'independent_pair_diagnostics': paired['comparisons'] if paired else None,
                       'development_macro_raw_MAE_skill': metric['independent_supported_cell_macro_raw_MAE_skill'] if metric else None,
                       'independent_study_success_established': False})
    neuronal = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'
    spatial = read(neuronal / 'additional_activity_spatial_and_reference_receipt.json')
    core = read(neuronal / 'additional_activity_closed_core_receipt.json')
    independent_core = read(neuronal / 'additional_activity_independent_no_prefix_audit.json')
    held_region = read(OUT / 'completed_held_region_matched_comparison.json')
    avi = OUT / 'AVI_preparation'
    avi_prelaunch = read(avi / 'matched_AVI_comparison_prelaunch_check.json')
    avi_prelaunch_changed_sources = [name for name, digest in (avi_prelaunch or {}).get('code_sha256', {}).items()
                                    if sha(Path(__file__).with_name(name)) != digest]
    interruption_history = [read(path) for path in OUT.glob('verified_storage_interruption*.json')]
    interruption_history.sort(key=lambda row: row['time_UTC'])
    interruption = interruption_history[-1] if interruption_history else None
    assessment = read(OUT / 'native_comparison_assessment.json')
    storage_stopped = bool(interruption and not live and not assessment)
    resource_proof = read(OUT / 'readonly_native_sustained_resource_audit.json')
    if assessment:
        next_action = 'The finite native comparison has reached its recorded decision. Follow the separate prespecified permitted-score input comparison only after its own readiness checks; do not reopen the native epoch or adaptation search.'
    elif storage_stopped and interruption.get('next_unstarted_fit'):
        next_action = 'The next fit has not started because disk space was below its recorded launch reserve. Preserve every completed fit and audit. Recover sufficient storage or establish an explicitly reviewed resource-only execution policy before continuing the unstarted role; do not retrain completed fits or change scientific gates.'
    elif storage_stopped and resource_proof and resource_proof['passed']:
        next_action = 'The read-only allocation passed the sustained resource check. Complete the attested wrapper recovery audit, then continue the pinned interrupted checkpoint under the explicitly recorded execution revision. Keep the original model, targets, batches and gates.'
    elif storage_stopped:
        next_action = 'Free several GiB on the system volume before resuming the interrupted fit, or establish an independently audited allocation-only execution revision. Do not change its model or insert new inputs.'
    else:
        next_action = 'Finish the fixed entire-source comparison without restarting it or adding inputs. Apply its prespecified gate. Keep independently audited auxiliary data separate until a new comparison is justified and fixed.'
    result = {'time_UTC': timestamp, 'actual_live_processes': live, 'verified_live_worker_progress': progress,
              'scientific_fits_completed': len(list(development.glob('*_record.json'))),
              'completed_saved_epochs': epochs, 'comparison_assessment': assessment,
              'verified_storage_interruption': interruption,
              'verified_storage_interruption_history': interruption_history,
              'scientific_comparison_currently_stopped_for_storage': storage_stopped,
              'storage_recovery': {
                  'cache_cleanup': read(OUT / 'npm_download_cache_cleanup_20261109.json'),
                  'repository_index_metadata_cleanup': read(OUT / 'conda_index_metadata_cache_cleanup_20261109.json'),
                  'remaining_download_cache_cleanup': read(OUT / 'remaining_conda_cache_cleanup_result_20261109.json'),
                  'lossless_scientific_asset_compression': read(OUT / 'lossless_scientific_asset_filesystem_compression_20261109.json'),
                  'lossless_text_asset_compression': read(OUT / 'lossless_text_assets_filesystem_compression_20261109.json'),
                  'lossless_additional_asset_compression': read(OUT / 'lossless_additional_assets_filesystem_compression_20261109.json'),
                  'lossless_remaining_binary_asset_compression': read(OUT / 'lossless_remaining_binary_assets_filesystem_compression_20261109.json'),
                  'verified_redundant_download_chunks_recovery': read(OUT / 'verified_redundant_download_chunks_recovery_20261109.json'),
                  'public_source_offload_prepared_not_approved_or_executed': read(OUT / 'proposed_public_source_archive_offload_20261109.json'),
                  'adapted_Retallick_launch': read(OUT / 'readonly_768_adapt_Retallick_launch_20261109.json'),
                  'adapted_Retallick_verified_progress': read(OUT / 'readonly_768_adapt_Retallick_launch_verification_20261109.json'),
                  'next_unstarted_fit_relaunch': read(OUT / 'readonly_scientific_next_fit_launch_20261109.json'),
                  'resource_only_768_MiB_admission': read(OUT / 'readonly_runtime_768_admission_audit.json'),
                  'resource_only_768_MiB_launch': read(OUT / 'readonly_768_scientific_launch_20261109.json'),
                  'relaunch': read(OUT / 'storage_recovery_relaunch_20261109.json'),
                  'restoration_and_advancing_updates': read(OUT / 'storage_recovery_resume_verification_20261109.json')},
              'completed_held_region_matched_comparison': held_region,
              'completed_Salomon_matched_source_comparison': read(OUT / 'completed_Salomon_matched_comparison.json'),
              'completed_Kosicki_frozen_arm_not_full_comparison': read(OUT / 'completed_Kosicki2025_frozen_arm.json'),
              'completed_Kosicki_matched_comparison': read(OUT / 'completed_Kosicki2025_matched_comparison.json'),
              'completed_Retallick_frozen_arm_not_full_comparison': read(OUT / 'completed_Retallick2026_frozen_arm.json'),
              'completed_Retallick_matched_comparison': read(OUT / 'completed_Retallick2026_matched_comparison.json'),
              'independent_native_comparison_decision_replay': read(OUT / 'native_comparison_independent_decision_audit.json'),
              'Retallick_initial_metadata_matching_not_completed_comparison': read(OUT / 'Retallick2026_initial_matched_arm_audit.json'),
              'development_cell_and_nominal_unit_coverage_audit': read(OUT / 'development_context_coverage_audit.json'),
              'nearby_experimental_replicate_repeatability_diagnostic': read(OUT / 'nearby_assay_replicate_repeatability.json'),
              'readonly_native_engineering': {
                  'actual_live_benchmark_processes': engineering_live,
                  'four_update_reference_comparison_and_recovery': read(OUT / 'readonly_native_recovery_audit.json'),
                  'resource_profile': read(OUT / 'recovery_audit_runs/readonly_native_adapt_resource1000/resource_profile.json'),
                  'resource_completion': read(OUT / 'recovery_audit_runs/readonly_native_adapt_resource1000/initial_segment_result.json'),
                  'sustained_resource_audit': resource_proof,
                  'attested_execution_revision': read(OUT / 'readonly_native_execution_revision.json'),
                  'attested_resource_only_768_MiB_execution_revision': read(OUT / 'readonly_native_execution_revision_768.json'),
                  'actual_wrapper_recovery_audit': read(OUT / 'readonly_runtime_recovery_audit.json'),
                  'scientific_migrations': [read(path) for path in development.glob('*_readonly_migration.json')],
                  'completed_recovery_payload_audits': [read(path) for path in OUT.glob('*_completed_recovery_audit.json')],
                  'used_in_scientific_comparison': bool(list(development.glob('*_readonly_migration.json'))
                                                      or any('readonly_native_runtime' in row['command'] for row in live)),
                  'prediction_improvement_established': False},
              'corrected_permitted_AVI_inputs': {
                  'preparation': read(avi / 'corrected_permitted_AVI_preparation_receipt.json'),
                  'independent_input_audit': read(avi / 'corrected_permitted_AVI_independent_audit.json'),
                  'training_only_scaling_audit': read(avi / 'permitted_AVI_train_only_scaling_audit.json'),
                  'expanded_head_forward_and_gradient_audit': read(avi / 'permitted_AVI_expanded_head_audit.json'),
                  'full_native_runtime_audit': read(avi / 'permitted_AVI_native_runtime_audit.json'),
                  'lower_memory_input_equivalence_audit': read(avi / 'readonly_native_input_equivalence_audit.json'),
                  'separate_trainer_source_derivation': read(avi / 'matched_AVI_trainer_derivation_receipt.json'),
                  'subsequent_resource_only_trainer_guard_derivation': read(avi / 'matched_AVI_trainer_resource_guard_derivation.json'),
                  'resource_admission_revision_not_AVI_runtime_qualification': read(avi / 'matched_AVI_resource_admission_revision.json'),
                  'separate_training_readiness': read(avi / 'matched_AVI_preparation_status.json'),
                  'separate_driver_prelaunch_checks_historical_snapshot': avi_prelaunch,
                  'prelaunch_snapshot_matches_current_sources': bool(avi_prelaunch and not avi_prelaunch_changed_sources),
                  'sources_changed_since_historical_prelaunch_checks': avi_prelaunch_changed_sources,
                  'seven_source_candidate_gate_synthetic_checks_only': read(avi / 'seven_source_candidate_gate_synthetic_checks.json'),
                  'execution_auditor_preparation_not_execution_success': read(avi / 'execution_auditor_preparation_and_live_refusal.json'),
                  'lower_memory_coordinator_preparation_not_execution_success': read(avi / 'AVI_coordinator_Torch_isolation_preparation.json'),
                  'locked_fresh_study_benefit_protocol': read(avi / 'fresh_study_benefit_protocol.json'),
                  'external_statistics_synthetic_checks_only': read(avi / 'external_magnitude_statistics_synthetic_checks.json'),
                  'used_in_current_comparison': False,
                  'prediction_benefit_established': False},
              'recovery_state_engineering_audit_passed': read(OUT / 'recovery_state_audit.json')['passed'],
              'prior_interrupted_partial_files_preserved': True,
              'confidence_method_experiment_on_archived_development': read(ROOT / 'outputs/sleep_frailty_ml_confidence_revision/confidence_method_decision.json'),
              'additional_neuron_activity_data': {'initial_spatial_reference_screen': spatial,
                                                  'shared_core_screen': core, 'independent_no_prefix_audit': independent_core,
                                                  'prepared_package': read(neuronal / 'audited_additional_neuron_activity_package_receipt.json'),
                                                  'target_description_correction': read(neuronal / 'activity_target_interpretation_correction.json'),
                                                  'used_in_current_comparison': False},
              'additional_neural_validation_candidate': read(ROOT / 'outputs/sleep_frailty_ml_neural_validation_candidates/Lee2025/candidate_availability_and_reservation.json'),
              'additional_neural_design_header_followup': read(ROOT / 'outputs/sleep_frailty_ml_neural_validation_candidates/Lee2025/design_availability_followup.json'),
              'requirements': {
                  'actual_measured_effect_size': ('Incomplete: all eight fixed phase-one fits finished. The independently replayed accuracy gate failed; the additional four source roles were not run. No useful cross-study effect-size candidate was established.' if assessment and not assessment['gate']['passed'] else ('Incomplete: matched held-region fits are complete; updating the reader did not improve their final-pass error. Entire-source comparisons remain.' if held_region and held_region['held_region_matched_fits_complete'] else 'Incomplete: matched held-region and entire-source comparisons remain.')),
                  'nearby_and_cell_specific_differences': 'Incomplete: trained auxiliary objectives; useful performance has not been established.',
                  'useful_confidence_across_experiments': 'Incomplete: the finite selected-calibration prototype on the archived model failed all seven source readiness checks; no externally useful calibrated confidence result.',
                  'permitted_AlphaGenome_added_input': 'Incomplete: repaired permissible inputs, scaling, expanded head, and full native runtime passed numerical audits. No scientific fit has used them; matched training and new-study benefit testing remain.'},
              'reserved_outcomes_read_for_model_selection': False, 'candidate_promoted': False, 'scientific_goal_complete': False,
              'next': next_action}
    (OUT / 'current_progress_and_readiness.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = [
        'Sleep–Frailty molecular-effect model: evidence at ' + timestamp,
        '',
        'Actual live scientific processes: ' + (', '.join(str(r['PID']) for r in live) or 'none'),
        'Actual live read-only engineering benchmarks: ' + (', '.join(str(r['PID']) for r in engineering_live) or 'none'),
        'Completed full fits: ' + str(result['scientific_fits_completed']),
        'Completed saved epochs: ' + str(len(epochs)),
        'Current storage stop: ' + str(storage_stopped),
        'A saved epoch is a partial result, not a completed comparison or independent validation.',
        '',
        'The comparison tests whether updating one part of an existing biology model helps more than leaving it fixed.',
        'Both versions receive the same measured targets, examples, training schedule and evaluation rules.',
        'Saved predictions have independent numerical checks; any small development gain does not prove new-study performance.',
        '',
        'The additional neuron data and corrected permissible score inputs were audited separately. They have not been used in the scientific comparison.',
        'The original interrupted work and negative findings remain preserved.',
        '',
        *result['requirements'].values(),
        '',
        'No model is promoted for project ranking or individual sleep/frailty risk prediction.',
        '',
        result['next'],
    ]
    (OUT / 'PLAIN_LANGUAGE_PROGRESS.txt').write_text('\n'.join(lines) + '\n')
    print(json.dumps({k: result[k] for k in ['time_UTC', 'actual_live_processes', 'verified_live_worker_progress', 'scientific_fits_completed', 'completed_saved_epochs', 'scientific_goal_complete']}), flush=True)


if __name__ == '__main__':
    main()
