#!/usr/bin/env python3
"""Prepare (default) or later execute only 16 independently admitted controls.

Preparation imports no numerical runtime, opens no GWAS outcome columns, and
launches no estimator. An explicit reviewed admission bound to the plan and all
190 historical completion receipts is required for every future worker.
"""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from canonical_calibration_common_v4_3 import ADMISSION, PACKAGE, PLAN, ROOT, SHARED_LOCK, TERMINATION_REQUEST, admit, assert_no_termination, check_hashes, deferred_termination_signals, exclusive_heavy_lock, safe_diagnostic, sha, utc, write_new

HISTORICAL = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
BLOCKS = PACKAGE / 'manifests/canonical_200_interval_preparation_v4.json'
CONTROLS = PACKAGE / 'statistical_validation/genomicsem_native_calibration_controls_v4.tsv'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/statistical_validation/native_canonical_calibration_v4_3')
ALLOWED = ['CORE_INSOMNIA_BMI_DUPLICATE', 'FIRST_FROZEN_VALIDATION_PAIR_DUPLICATE']


def prepare():
    if PLAN.exists() or ADMISSION.exists():
        raise RuntimeError('EXISTING_FROZEN_PREPARATION_REQUIRES_REVIEW_NO_OVERWRITE')
    native = json.loads(HISTORICAL.read_text())
    blocks = json.loads(BLOCKS.read_text())
    if blocks['native_plan_sha256'] != sha(HISTORICAL) or blocks['build'] != 'GRCh37/hg19' or blocks['block_count'] != 200 or not blocks['all_source_reference_hashes_unchanged_after']:
        raise RuntimeError('CANONICAL_PARTITION_NOT_ADMITTED_FOR_PREPARATION')
    if len(blocks['reference_sha256']) != 22 or any(native['dependencies_sha256'].get(p) != h for p, h in blocks['reference_sha256'].items()):
        raise RuntimeError('CONSUMED_REFERENCE_IDENTITIES_DO_NOT_MATCH')
    if sha(blocks['interval_tsv_path']) != blocks['interval_tsv_sha256'] or sha(blocks['coordinate_map_path']) != blocks['coordinate_map_sha256']:
        raise RuntimeError('FROZEN_COORDINATES_OR_INTERVALS_CHANGED')
    verified = {r['path']: r['actual_sha256'] for r in native['inputs_verified'] if r['match']}
    with CONTROLS.open(newline='') as f:
        controls = list(csv.DictReader(f, delimiter='\t'))
    if [r['control_id'] for r in controls] != ALLOWED:
        raise RuntimeError('PRESPECIFIED_CONTROL_MEMBERSHIP_CHANGED')
    dependencies = dict(native['dependencies_sha256'])
    small_paths = [Path(__file__), PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py',
                   PACKAGE / 'scripts/canonical_calibration_common_v4_3.py', PACKAGE / 'scripts/canonical_native_adapter.py',
                   PACKAGE / 'scripts/canonical_calibration_common.py',
                   PACKAGE / 'scripts/canonical_calibration_capture_v4_3.py', PACKAGE / 'scripts/canonical_calibration_arithmetic_v4_3.py',
                   PACKAGE / 'scripts/canonical_calibration_fixture_v4_3.py',
                   PACKAGE / 'scripts/canonical_calibration_operational_fixture_v4_3.py',
                   PACKAGE / 'scripts/canonical_calibration_cleanup_fixture_v4_3.py',
                   PACKAGE / 'scripts/canonical_calibration_signal_diagnostic_fixture_v4_3.py',
                   PACKAGE / 'scripts/52_acquire_extension_raw_sources_v4.py',
                   PACKAGE / 'reviews/independent_canonical_executor_correction_review_v4_2.md',
                   PACKAGE / 'reviews/independent_canonical_executor_correction_review_v4_2.sha256',
                   PACKAGE / 'scripts/41_prepare_canonical_blocks.py', HISTORICAL, BLOCKS, CONTROLS,
                   PACKAGE / 'statistical_validation/common_boundary_adapter_prototype_v4.py',
                   PACKAGE / 'statistical_validation/canonical_native_calibration_fixture_receipt_v4_4.json',
                   PACKAGE / 'statistical_validation/canonical_operational_correction_fixture_receipt_v4_3.json',
                   PACKAGE / 'statistical_validation/canonical_cleanup_and_final_guard_fixture_receipt_v4_3.json',
                   PACKAGE / 'statistical_validation/canonical_signal_diagnostic_post_hash_fixture_receipt_v4_3.json',
                   Path(blocks['interval_tsv_path']), Path(blocks['coordinate_map_path'])]
    dependencies.update({str(p): sha(p) for p in small_paths})
    python = next(p for p in native['dependencies_sha256'] if p.endswith('/.ldsc-env/bin/python'))
    ldsc_dir = str(ROOT.parent / 'ldsc-code')
    refs = str(Path(next(iter(blocks['reference_sha256']))).parent) + '/'
    jobs = []
    source_hashes = {}
    for control in controls:
        inputs = [control['input1'], control['input2']]
        identities = {control['input1']: control['input1_sha256'], control['input2']: control['input2_sha256']}
        if any(verified.get(p) != h for p, h in identities.items()):
            raise RuntimeError('CONTROL_SOURCES_NOT_EXACTLY_HISTORICALLY_ADMITTED')
        source_hashes.update(identities)
        for method in ('one_step', 'two_step30'):
            for partition in ('stock_default200', 'canonical200'):
                for duplicate in ('A', 'B'):
                    job_id = control['control_id'] + '__' + method + '__' + partition + '__' + duplicate
                    out = SSD / job_id
                    ldsc_args = ['--rg', ','.join(inputs), '--ref-ld-chr', refs, '--w-ld-chr', refs,
                                 '--n-blocks', '200', '--print-delete-vals', '--out', str(out / 'stock_ldsc')]
                    if method == 'two_step30':
                        ldsc_args += ['--two-step', '30']
                    jobs.append({'job_id': job_id, 'control_id': control['control_id'], 'method': method, 'partition': partition,
                                 'duplicate': duplicate, 'inputs': inputs, 'input_sha256': identities, 'output_dir': str(out), 'ldsc_args': ldsc_args})
    # Compressed input identity reads only; no GWAS rows/phenotype results read.
    current_inputs = check_hashes(source_hashes)
    if len(jobs) != 16 or len(current_inputs) != 3:
        raise RuntimeError('EXPECTED_16_CONTROLS_AND_THREE_UNIQUE_SOURCES')
    record = {'schema': 'native_canonical_calibration_operational_plan_v1', 'prepared_utc': utc(), 'self_path': str(PLAN),
              'operational_revision': 'v4_3',
              'superseded_preparation_plan_sha256': sha(PACKAGE / 'manifests/native_canonical_calibration_plan_v4_2.json'),
              'superseded_preparation_seal_sha256': sha(PACKAGE / 'reviews/genomicsem_native_calibration_preparation_receipt_v4_2.json'),
              'operational_corrections': ['Exact relocated historical capture-wrapper identity mapping', 'Historical output and stage-plan hash verification',
                                          'Failure receipts and lock retention through owned-group cleanup', 'Residual-child and final-resource guards',
                                          'Cross-arm SNP/mask/Nbar/M and weighted-array comparisons', 'Effective method and free-intercept assertion before regression',
                                          'Deferred SIGTERM/SIGHUP/SIGINT through ownership and cleanup', 'Non-throwing diagnostics', 'Final guards after all output hashing'],
              'scope': 'METHOD_CALIBRATION_ONLY', 'allowed_control_ids': ALLOWED, 'jobs': jobs,
              'historical_plan_path': str(HISTORICAL), 'historical_plan_sha256': sha(HISTORICAL),
              'block_receipt_path': str(BLOCKS), 'block_receipt_sha256': sha(BLOCKS),
              'build': blocks['build'], 'interval_tsv_sha256': blocks['interval_tsv_sha256'],
              'ordered_all_reference_records_sha256': blocks['ordered_all_reference_records_sha256'],
              'ordered_construction_records_sha256': blocks['ordered_construction_records_sha256'],
              'coordinate_map_path': blocks['coordinate_map_path'], 'coordinate_map_sha256': blocks['coordinate_map_sha256'],
              'partition_amendment': blocks['partition_amendment'], 'static_construction_exclusions': blocks['static_construction_exclusions'],
              'reference_sha256': blocks['reference_sha256'], 'dependencies_sha256': dependencies,
              'inputs_current_sha256': current_inputs, 'python': python, 'ldsc_dir': ldsc_dir, 'ssd_output_root': str(SSD),
              'environment': native['environment'], 'worker_count': 1, 'blas_threads': 1,
              'shared_lock_path': str(SHARED_LOCK), 'shared_lock_policy': 'LOCK_EX|LOCK_NB for entire executor lifetime, inherited descriptor required by every native/arithmetic worker.',
              'internal_floor_bytes': 3 * 1024**3, 'ssd_floor_bytes': 5 * 1024**3, 'worker_rss_limit_bytes': 2 * 1024**3,
              'output_limit_bytes': 8 * 1024**3, 'per_worker_seconds_limit': 2 * 3600, 'runtime_poll_seconds': 2,
              'resource_storage_assumption': '16 controls retain exact weighted arrays; 8 GiB total cap, plus 5 GiB SSD free floor; admission review must assess capacity before launch.',
              'operational_resource_amendment': 'Capture output cap raised from prototype 4 GiB to 8 GiB to retain all 16 prespecified control designs/deletes; this separate operational amendment requires root review before execution. Scientific thresholds/inputs/partition remain unchanged.',
              'arithmetic': {'point_delete_atol': 1e-10, 'point_delete_rtol': 1e-8, 'covariance_atol': 1e-14, 'covariance_rtol': 1e-6,
                             'duplicate_zero_delta_se_atol': 1e-12, 'independent_row_solves': 'replicate A of all 8 control/arm/partition combinations',
                             'full_fit_weights_and_Nbar_frozen_for_deletion': True, 'no_calibrated_difference_p_values': True},
              'cross_arm_identity_rule': 'Exact final SNP order, mask indices, Nbar/M and source arrays. Record initial/final weighted-array hashes and compare arrays under the existing point/delete tolerance; native sequential Gencov and two-step response may inherit binary64 block-sum roundoff. Independent root review required before admission.',
              'native_jobs_launched': 0, 'GWAS_outcome_columns_read': False, 'allow_41_covariance_outcomes': False,
              'allow_calibrated_biological_p_values': False, 'realistic_LD_sampling_calibration_pass': False,
              'remaining_method_gate': 'Verified source-matched signed LD/genotype asset and frozen realistic-LD simulations; synthetic/duplicate arithmetic cannot establish real genomic correlation length, selection, or reference validity.'}
    write_new(PLAN, record)
    write_new(ADMISSION, {'schema': 'native_canonical_calibration_execution_admission_v1', 'prepared_utc': utc(),
                         'plan_sha256': sha(PLAN), 'scope': 'METHOD_CALIBRATION_ONLY', 'execution_admitted': False,
                         'independent_binding_review_pass': False, 'resource_plan_review_pass': False,
                         'independent_review_path': None, 'independent_review_sha256': None,
                         'historical_campaign_receipt_sha256': {}, 'allow_41_covariance_outcomes': False,
                         'allow_calibrated_biological_p_values': False,
                         'blocked_reason': 'Historical 190-job campaign active; independent executor binding/resource review and exact completion receipts required before any fit.'})
    print(json.dumps({'status': 'PREPARED_EXECUTION_BLOCKED', 'jobs': len(jobs), 'plan_sha256': sha(PLAN), 'native_jobs_launched': 0}))


def monitor_helpers():
    path = PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py'
    spec = importlib.util.spec_from_file_location('canonical_calibration_owned_process_helpers', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def snapshot(output_root):
    return {'recorded_utc': utc(), 'internal_free_bytes': shutil.disk_usage('/System/Volumes/Data').free,
            'ssd_free_bytes': shutil.disk_usage(output_root.parent).free}


def output_bytes(output_root):
    return sum(p.stat().st_size for p in output_root.rglob('*') if p.is_file())


def owned_group_exists(pid):
    try:
        os.killpg(pid, 0)
        return True
    except ProcessLookupError:
        return False


def await_owned_cleanup(ownership, plan):
    """Never return, including after interruption, while an owned group survives."""
    helper = monitor_helpers()
    while ownership:
        entry = ownership[0]
        proc = entry['proc']
        try:
            result = helper.terminate_owned(proc)
            if owned_group_exists(proc.pid):
                raise RuntimeError('OWNED_GROUP_EXISTS_AFTER_HELPER_TEARDOWN')
        except BaseException as error:
            entry['cleanup_errors'].append(repr(error))
            try:
                # Only the session/group created by this executor is targeted.
                # A signal-0 probe verifies disappearance independently of ps.
                if owned_group_exists(proc.pid):
                    os.killpg(proc.pid, signal.SIGKILL)
                proc.wait(timeout=10)
                if owned_group_exists(proc.pid):
                    raise RuntimeError('OWNED_GROUP_REMAINS_AFTER_INDEPENDENT_KILL')
                result = {'helper_error': repr(error), 'independent_signal': 'SIGKILL_IF_GROUP_PRESENT',
                          'remaining_group_members': [], 'signal_zero_group_absent': True}
            except BaseException as fallback_error:
                entry['cleanup_errors'].append(repr(fallback_error))
                safe_diagnostic(lambda: 'SHARED_HEAVY_WORKER_LOCK_RETAINED: owned group ' + str(proc.pid) + ': ' + repr(fallback_error))
                try:
                    time.sleep(2)
                except BaseException:
                    pass
                continue
        entry['verified_recovery'] = result
        ownership.pop(0)
        recovery = {'completed_utc': utc(), 'owned_process_group_id': proc.pid, 'cleanup_errors': entry['cleanup_errors'],
                    'process_group_teardown': result, 'group_disappearance_verified_before_shared_lock_release': True}
        path = entry['output_dir'] / 'owned_group_cleanup_recovery.json'
        try:
            write_new(path, recovery)
        except BaseException as proof_error:
            entry['cleanup_errors'].append('RECOVERY_PROOF_WRITE: ' + repr(proof_error))
            safe_diagnostic(lambda: 'OWNED_GROUP_GONE_RECOVERY_PROOF_WRITE_FAILED: ' + repr(proof_error))


def final_violations(plan, state, size, elapsed, admission_sha, plan_sha):
    violations = []
    if state['internal_free_bytes'] < plan['internal_floor_bytes']: violations.append('INTERNAL_RESOURCE_FLOOR')
    if state['ssd_free_bytes'] < plan['ssd_floor_bytes']: violations.append('SSD_RESOURCE_FLOOR')
    if size > plan['output_limit_bytes']: violations.append('OUTPUT_LIMIT')
    if elapsed > plan['per_worker_seconds_limit']: violations.append('WORKER_DEADLINE')
    if sha(ADMISSION) != admission_sha: violations.append('ADMISSION_CHANGED_DURING_RUN')
    if sha(PLAN) != plan_sha: violations.append('PLAN_CHANGED_DURING_RUN')
    if TERMINATION_REQUEST: violations.append('DEFERRED_TERMINATION_SIGNAL: ' + str(TERMINATION_REQUEST))
    return violations


def monitored_worker(plan, admission_sha, job_id, command, output_dir, lock_fd, ownership):
    helper = monitor_helpers()
    root = Path(plan['ssd_output_root'])
    plan_sha = sha(PLAN)
    env = {**os.environ, 'TMPDIR': str(root / 'tmp'), 'XDG_CACHE_HOME': str(root / 'cache'),
           'MPLCONFIGDIR': str(root / 'cache/matplotlib'), 'PYTHONDONTWRITEBYTECODE': '1',
           'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
    for path in ('tmp', 'cache/matplotlib'):
        (root / path).mkdir(parents=True, exist_ok=True)
    logfile = output_dir / 'worker_stdout.log'
    proc = None
    started = time.monotonic()
    peak = 0
    reason = None
    teardown = None
    returncode = None
    samples = []
    run_error = None
    entry = None
    first_teardown = None
    final_state = None
    violations = []
    output_hashes = {}
    output_hash_error = None
    try:
        assert_no_termination()
        state = snapshot(root)
        if state['internal_free_bytes'] < plan['internal_floor_bytes'] or state['ssd_free_bytes'] < plan['ssd_floor_bytes']:
            reason = 'PRELAUNCH_RESOURCE_FLOOR_NO_WORKER_LAUNCHED'
            raise RuntimeError(reason)
        with logfile.open('xb') as stream:
            proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True, pass_fds=(lock_fd,))
            entry = {'proc': proc, 'output_dir': output_dir, 'cleanup_errors': []}
            ownership.append(entry)
            while True:
                returncode = proc.poll()
                rss, pids = helper.owned_rss(proc.pid)
                peak = max(peak, rss)
                state = snapshot(root)
                size = output_bytes(root)
                elapsed = time.monotonic() - started
                current = final_violations(plan, state, size, elapsed, admission_sha, plan_sha)
                if rss > plan['worker_rss_limit_bytes']: current.insert(0, 'OWNED_WORKER_RSS_LIMIT')
                if current: reason = current[0]
                if len(samples) < 100 or reason or returncode is not None:
                    samples.append({**state, 'elapsed_seconds': elapsed, 'rss_bytes': rss, 'owned_pids': pids, 'output_bytes': size})
                if reason or returncode is not None: break
                time.sleep(plan['runtime_poll_seconds'])
    except BaseException as error:
        run_error = {'type': type(error).__name__, 'repr': repr(error)}
        if reason is None: reason = 'MONITOR_OR_WORKER_LAUNCH_EXCEPTION'
    finally:
        if proc is not None:
            if entry is None:
                entry = {'proc': proc, 'output_dir': output_dir, 'cleanup_errors': []}
                ownership.append(entry)
            natural_success = returncode == 0 and reason is None and run_error is None
            try:
                teardown = helper.terminate_owned(proc)
                first_teardown = {'helper_return': teardown, 'signal_zero_group_absent': not owned_group_exists(proc.pid)}
                if not first_teardown['signal_zero_group_absent']:
                    raise RuntimeError('OWNED_GROUP_EXISTS_AFTER_FIRST_TEARDOWN')
                ownership.remove(entry)
                if natural_success and teardown.get('initial_group_members'):
                    reason = 'RESIDUAL_OWNED_GROUP_AFTER_PARENT_EXIT'
            except BaseException as error:
                entry['cleanup_errors'].append(repr(error))
                if first_teardown is None: first_teardown = {}
                first_teardown['exception'] = {'type': type(error).__name__, 'repr': repr(error)}
                if reason is None: reason = 'OWNED_GROUP_TEARDOWN_EXCEPTION'
            # Preserve first failure/teardown evidence BEFORE retries. A failing
            # helper may already have killed/reaped some of the original group.
            try:
                write_new(output_dir / 'first_teardown_attempt.json', {'recorded_utc': utc(), 'job_id': job_id,
                          'first_teardown_attempt': first_teardown, 'run_error': run_error, 'stop_reason': reason})
                if reason:
                    write_new(output_dir / 'monitor_failure_receipt.json', {'recorded_utc': utc(), 'job_id': job_id,
                              'plan_sha256': plan_sha, 'returncode_before_cleanup': returncode, 'stop_reason': reason,
                              'first_teardown_attempt': first_teardown, 'run_error': run_error,
                              'shared_lock_retained_until_verified_group_cleanup': True})
            except BaseException as proof_error:
                entry['cleanup_errors'].append('FIRST_PROOF_WRITE: ' + repr(proof_error))
                if reason is None: reason = 'MONITOR_PROOF_WRITE_EXCEPTION'
            finally:
                if ownership:
                    await_owned_cleanup(ownership, plan)
                if entry.get('verified_recovery') is not None:
                    teardown = entry['verified_recovery']
            returncode = proc.returncode
        try:
            # Bulk output audit belongs BEFORE the terminal floors/deadline and
            # identity observation. Signals during hashing remain deferred.
            output_hashes = {str(p): sha(p) for p in output_dir.rglob('*') if p.is_file()}
        except BaseException as error:
            output_hash_error = {'type': type(error).__name__, 'repr': repr(error)}
            if reason is None: reason = 'POST_TEARDOWN_OUTPUT_HASH_EXCEPTION'
        try:
            final_state = snapshot(root)
            final_size = output_bytes(root)
            violations = final_violations(plan, final_state, final_size, time.monotonic() - started, admission_sha, plan_sha)
            if violations and reason is None: reason = violations[0]
        except BaseException as error:
            violations.append('POST_TEARDOWN_GUARD_EXCEPTION: ' + repr(error))
            if reason is None: reason = 'POST_TEARDOWN_GUARD_EXCEPTION'
        if returncode != 0 and reason is None: reason = 'NATIVE_WORKER_NONZERO_EXIT'
        receipt = {'job_id': job_id, 'command': command, 'completed_utc': utc(), 'plan_sha256': plan_sha,
                   'admission_sha256': admission_sha, 'returncode': returncode, 'stop_reason': reason,
                   'elapsed_seconds': time.monotonic() - started, 'peak_owned_rss_bytes': peak,
                   'resource_samples': samples, 'process_group_teardown': teardown,
                   'first_teardown_attempt': first_teardown, 'run_error': run_error,
                   'cleanup_errors': [] if entry is None else entry['cleanup_errors'],
                   'final_resource_state_after_teardown': final_state, 'final_guard_violations': violations,
                   'bulk_output_hashes_completed_before_final_guards': output_hash_error is None,
                   'output_hash_error': output_hash_error, 'deferred_termination_requests': list(TERMINATION_REQUEST),
                   'owned_group_disappearance_verified': proc is None or not owned_group_exists(proc.pid),
                   'all_output_sha256': output_hashes}
        write_new(output_dir / 'monitor_receipt.json', receipt)
    if returncode != 0 or reason:
        raise RuntimeError('CALIBRATION_WORKER_FAILED_PRESERVED_PARTIAL_OUTPUT: ' + job_id)


def execute(lock_fd, ownership):
    assert_no_termination()
    plan, _ = admit()
    admission_sha = sha(ADMISSION)
    root = Path(plan['ssd_output_root'])
    if root.exists():
        raise RuntimeError('EXISTING_CALIBRATION_NAMESPACE_REQUIRES_REVIEW_NO_OVERWRITE')
    root.mkdir(parents=True, exist_ok=False)
    for job in plan['jobs']:
        assert_no_termination()
        # Recheck admission/history before EACH worker, never overlap native work.
        admit()
        out = Path(job['output_dir'])
        out.mkdir(exist_ok=False)
        command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_capture_v4_3.py'),
                   '--plan', str(PLAN), '--admission', str(ADMISSION), '--plan-sha', sha(PLAN), '--job-id', job['job_id'], '--shared-lock-fd', str(lock_fd)]
        monitored_worker(plan, admission_sha, job['job_id'], command, out, lock_fd, ownership)
        receipt = json.loads((out / 'native_capture_receipt.json').read_text())
        if receipt['job'] != job or receipt['plan_sha256'] != sha(PLAN):
            raise RuntimeError('NATIVE_CAPTURE_RECEIPT_NOT_BOUND')
        print(json.dumps({'completed_method_control': job['job_id'], 'biological_outcomes_computed': 0}), flush=True)
    verify(lock_fd, ownership)
    assert_no_termination()


def verify(lock_fd, ownership):
    assert_no_termination()
    plan, _ = admit()
    for job in plan['jobs']:
        out = Path(job['output_dir'])
        r = json.loads((out / 'monitor_receipt.json').read_text())
        if r['returncode'] != 0 or r['stop_reason'] is not None or r['process_group_teardown']['remaining_group_members'] != []:
            raise RuntimeError('ALL_16_NATIVE_CONTROLS_MUST_COMPLETE_BEFORE_ARITHMETIC')
        if r['job_id'] != job['job_id'] or r['plan_sha256'] != sha(PLAN) or r['admission_sha256'] != sha(ADMISSION) or r['final_guard_violations'] or r['owned_group_disappearance_verified'] is not True or r['bulk_output_hashes_completed_before_final_guards'] is not True or r['output_hash_error'] is not None or r['deferred_termination_requests']:
            raise RuntimeError('CONTROL_MONITOR_IDENTITY_OR_FINAL_GUARD_MISMATCH')
        check_hashes(r['all_output_sha256'])
        if not (out / 'native_capture_receipt.json').exists():
            raise RuntimeError('MISSING_NATIVE_CAPTURE')
    out = Path(plan['ssd_output_root']) / 'independent_arithmetic_worker'
    out.mkdir(exist_ok=False)
    command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_arithmetic_v4_3.py'), '--plan', str(PLAN), '--admission', str(ADMISSION), '--shared-lock-fd', str(lock_fd)]
    monitored_worker(plan, sha(ADMISSION), 'INDEPENDENT_CAPTURED_ROW_ARITHMETIC', command, out, lock_fd, ownership)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=('prepare', 'execute', 'verify'), nargs='?', default='prepare')
    args = p.parse_args()
    if args.mode == 'prepare':
        prepare()
    else:
        ownership = []
        def before_release():
            if ownership:
                # Teardown must remain possible even if the plan becomes
                # unreadable. Cleanup uses only recorded process ownership.
                await_owned_cleanup(ownership, None)
        # Restore signal dispositions only AFTER every registered owned group
        # is gone and the shared-lock cleanup callback has completed.
        with deferred_termination_signals():
            with exclusive_heavy_lock(before_release=before_release) as lock_fd:
                {'execute': execute, 'verify': verify}[args.mode](lock_fd, ownership)


if __name__ == '__main__':
    main()
