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

from canonical_calibration_common_v4_4 import ADMISSION, PACKAGE, PLAN, ROOT, SHARED_LOCK, TERMINATION_REQUEST, admit, assert_no_termination, check_hashes, deferred_termination_signals, exclusive_heavy_lock, safe_diagnostic, sha, utc, write_new

from terminal_commit_common_v2 import TerminalCommit

HISTORICAL = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
BLOCKS = PACKAGE / 'manifests/canonical_200_interval_preparation_v4.json'
CONTROLS = PACKAGE / 'statistical_validation/genomicsem_native_calibration_controls_v4.tsv'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/statistical_validation/native_canonical_calibration_v4_4')
ALLOWED = ['CORE_INSOMNIA_BMI_DUPLICATE', 'FIRST_FROZEN_VALIDATION_PAIR_DUPLICATE']


def prepare():
    if PLAN.exists() or ADMISSION.exists():raise RuntimeError('EXISTING_FROZEN_PREPARATION_REQUIRES_REVIEW_NO_OVERWRITE')
    prior=PACKAGE/'manifests/native_canonical_calibration_plan_v4_3.json'
    prior_sha='cfe5ef24f3455cc0508785e69b291d023f524b2e4872754598b864365f884c21'
    if sha(prior)!=prior_sha:raise RuntimeError('PRESERVED_CANONICAL_PLAN_CHANGED')
    plan=json.loads(prior.read_text());check_hashes(plan['dependencies_sha256'])
    old_root=plan['ssd_output_root']
    for job in plan['jobs']:
        job['output_dir']=job['output_dir'].replace(old_root+'/',str(SSD)+'/',1)
        job['ldsc_args']=[arg.replace(old_root+'/',str(SSD)+'/',1) if arg.startswith(old_root+'/') else arg for arg in job['ldsc_args']]
    files=[prior,Path(__file__),PACKAGE/'manifests/global_SSD_resource_reservation_v4_3.json']
    files.extend(PACKAGE/'scripts'/n for n in ['canonical_calibration_common_v4_4.py','canonical_calibration_capture_v4_4.py','canonical_calibration_arithmetic_v4_4.py','terminal_commit_common_v2.py','native_stage_completion_v4_3.py'])
    files.extend(PACKAGE/'reviews'/n for n in ['terminal_commit_common_review_seal_v2.json','genomicsem_native_stage_completion_review_v4_4.md','genomicsem_native_stage_completion_review_v4_4.json','independent_core_numerical_adjudication_v4.sha256','independent_whole_extension_v4.sha256','independent_whole_validation_v4.sha256','genomicsem_native_calibration_preparation_receipt_v4_3.json'])
    for file in files:plan['dependencies_sha256'][str(file)]=sha(file)
    plan.update(self_path=str(PLAN),operational_revision='v4_4',prepared_utc=utc(),ssd_output_root=str(SSD),preserved_v4_3_plan=str(prior),preserved_v4_3_plan_sha256=prior_sha,global_resource_ledger_path=str(files[2]),global_resource_ledger_sha256=sha(files[2]),global_reservation_bytes=300*(1<<30),terminal_protocol='DURABLE_PRIVATE_PENDING_BEFORE_WORK;ALL_CAPTURE_AND_MONITOR_OUTPUT_IDENTITIES;POST_PERSISTENCE_ADMISSION_HISTORY_RESOURCE_TERMINATION_CHECKS;ABSENT_PENDING_AND_FAILURE_ADDENDA_REQUIRED')
    if len(plan['jobs'])!=16:raise RuntimeError('EXACT16_CANONICAL_CONTROLS_REQUIRED')
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'ssd_free_bytes':shutil.disk_usage('/Volumes/Extreme SSD').free}
    if plan['resource_preflight']['internal_free_bytes']<plan['internal_floor_bytes'] or plan['resource_preflight']['ssd_free_bytes']<plan['ssd_floor_bytes']:raise RuntimeError('UNCHANGED_CANONICAL_FLOORS_NOT_MET')
    check_hashes(plan['dependencies_sha256']);write_new(PLAN,plan)
    print(json.dumps({'plan':str(PLAN),'sha256':sha(PLAN),'jobs':16,'workers_launched':0,'scope':'METHOD_CALIBRATION_ONLY'}))


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
        command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_capture_v4_4.py'),
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
    command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_arithmetic_v4_4.py'), '--plan', str(PLAN), '--admission', str(ADMISSION), '--shared-lock-fd', str(lock_fd)]
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
                if args.mode!='execute':raise RuntimeError('NEW_TERMINAL_EPOCH_REQUIRES_SINGLE_COMPLETE_EXECUTE_ROUTE')
                pending=SSD.parent/'native_canonical_calibration_pending_v4_4.json'
                seal=SSD.parent/'native_canonical_calibration_terminal_seal_v4_4.json'
                terminal=TerminalCommit(pending,seal,{'plan_sha256':sha(PLAN),'admission_sha256':sha(ADMISSION),'executor_sha256':sha(Path(__file__))})
                execute(lock_fd,ownership)
                plan,_=admit();root=Path(plan['ssd_output_root'])
                output_hashes={str(f):sha(f) for f in root.rglob('*') if f.is_file() and not f.is_symlink()}
                historical=__import__('canonical_calibration_common_v4_4').historical_completion(plan)
                result=SSD.parent/'native_canonical_calibration_execution_receipt_v4_4.json'
                record={'status':'ALL16_CANONICAL_SOFTWARE_CONTROLS_AND_CAPTURED_ROW_ARITHMETIC_COMPLETE','plan_sha256':sha(PLAN),'admission_sha256':sha(ADMISSION),'all_output_sha256':output_hashes,'historical190_receipt_sha256':historical,'completed_utc':utc(),'biological_outcomes_computed':0,'calibrated_difference_p_values':False}
                write_new(result,record)
                def final_identity():
                    if ownership:raise RuntimeError('CANONICAL_OWNED_GROUPS_NOT_EMPTY')
                    current,_=admit()
                    if current!=plan or __import__('canonical_calibration_common_v4_4').historical_completion(plan)!=historical:raise RuntimeError('CANONICAL_TERMINAL_HISTORY_OR_PLAN_CHANGED')
                    check_hashes(output_hashes)
                def final_resource():
                    assert_no_termination()
                    state=snapshot(root);size=output_bytes(root)+sum(f.stat().st_size for f in [pending,seal,result] if f.is_file())
                    reasons=final_violations(plan,state,size,0,sha(ADMISSION),sha(PLAN))
                    total=sum(f.stat().st_size for f in SSD.parents[1].rglob('*') if f.is_file() and not f.is_symlink())
                    if reasons or total>plan['global_reservation_bytes']:raise RuntimeError('CANONICAL_POST_PERSISTENCE_RESOURCE_GUARD: '+str(reasons))
                    assert_no_termination()
                if not terminal.commit({str(result):sha(result)},identity_gate=final_identity,resource_gate=final_resource,termination_gate=assert_no_termination):raise RuntimeError('CANONICAL_TERMINAL_COMMIT_FAILED_PENDING_REMAINS')


if __name__ == '__main__':
    main()
