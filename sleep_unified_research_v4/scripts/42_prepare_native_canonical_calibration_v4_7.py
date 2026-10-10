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

from canonical_calibration_common_v4_7 import ADMISSION, PACKAGE, PLAN, ROOT, SHARED_LOCK, TERMINATION_REQUEST, admit, assert_no_termination, check_hashes, deferred_termination_signals, exclusive_heavy_lock, safe_diagnostic, sha, utc, write_new, regular_sha, intended_sha, write_intended_new, terminal_intended_sha, freeze_companions, regular_tree, canonical_epoch_bytes, global_campaign_bytes, preserved_failed_attempt

from terminal_commit_common_v2 import TerminalCommit, require_committed

HISTORICAL = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
BLOCKS = PACKAGE / 'manifests/canonical_200_interval_preparation_v4.json'
CONTROLS = PACKAGE / 'statistical_validation/genomicsem_native_calibration_controls_v4.tsv'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/statistical_validation/native_canonical_calibration_v4_7')
ALLOWED = ['CORE_INSOMNIA_BMI_DUPLICATE', 'FIRST_FROZEN_VALIDATION_PAIR_DUPLICATE']


def prepare():
    if PLAN.exists() or ADMISSION.exists():raise RuntimeError('EXISTING_FROZEN_PREPARATION_REQUIRES_REVIEW_NO_OVERWRITE')
    prior=PACKAGE/'manifests/native_canonical_calibration_plan_v4_6.json'
    prior_sha='3bbbaa75b5520e569c3386eb32cc28967cd5595ce1beb7dbea0cb7202c23b6ce'
    if regular_sha(prior)!=prior_sha:raise RuntimeError('PRESERVED_CANONICAL_V4_6_PLAN_CHANGED')
    plan=json.loads(prior.read_text());check_hashes(plan['dependencies_sha256'])
    failed=PACKAGE/'logs/canonical_v4_6_failed_preserved_stop_inventory_v4_7.json'
    failed_sha='c80990a119349b0deda5b4d1c6ee15ceb1d16f3c79d0f9f457894b8998674215'
    if regular_sha(failed)!=failed_sha:raise RuntimeError('EXACT_PRESERVED_V4_6_FAILURE_REQUIRED')
    failure=json.loads(failed.read_text())
    if failure['status']!='FAILED_PRESERVED' or failure['preserved_plan_sha256']!=prior_sha or failure['completed_native_controls_credited']!=0:raise RuntimeError('PRESERVED_FAILURE_SCOPE_CHANGED')
    preserved={path:item['sha256'] for path,item in failure['all_preserved_file_identity'].items()}
    regular_hashes(preserved)
    old_root=plan['ssd_output_root']
    for job in plan['jobs']:
        job['output_dir']=job['output_dir'].replace(old_root+'/',str(SSD)+'/',1)
        job['ldsc_args']=[arg.replace(old_root+'/',str(SSD)+'/',1) if arg.startswith(old_root+'/') else arg for arg in job['ldsc_args']]
    files=[Path(__file__),prior,failed]
    files.extend(PACKAGE/'scripts'/n for n in ['canonical_calibration_common_v4_7.py','canonical_calibration_capture_v4_7.py','canonical_calibration_arithmetic_v4_7.py','terminal_commit_common_v2.py'])
    files.append(PACKAGE/'reviews/independent_canonical_executor_prelaunch_review_v4_6_seal.json')
    for file in files:plan['dependencies_sha256'][str(file)]=regular_sha(file)
    plan['dependencies_sha256'].update(preserved)
    plan.update(self_path=str(PLAN),operational_revision='v4_7',prepared_utc=utc(),ssd_output_root=str(SSD),
                preserved_v4_6_plan=str(prior),preserved_v4_6_plan_sha256=prior_sha,preserved_failed_attempt_receipt=str(failed),preserved_failed_attempt_receipt_sha256=failed_sha,
                preserved_v4_6_file_sha256=preserved,preserved_v4_6_review_artifact_sha256=failure['old_review_artifact_sha256'],
                global_reservation_bytes=300*(1<<30),all_retained_canonical_epochs_output_cap_bytes=plan['output_limit_bytes'],
                output_meter_policy='ALL_RETAINED_CANONICAL_EPOCHS_AND_CANONICAL_PENDING_TERMINAL_EXECUTION_METADATA_AND_APPLEDOUBLE;UNCHANGED8GIB',
                transport_companion_policy='EXACT_REGISTERED_TARGET_APPLEDOUBLE;REGULAR_NO_SYMLINK;VALID_HEADER_ENTRY_BOUNDS;FIRST_SHA_IMMUTABLE;NO_IGNORED_FILES_OR_DELETION',
                arithmetic_report_route=str(SSD/'independent_control_arithmetic.json'),
                terminal_protocol='PRIVATE_TERMINAL2_PENDING;INTENDED_RESULT_AND_SEAL_SHA_BEFORE_WRITES;MONITOR_SHA_BEFORE_PARSE;ALL_REGISTERED_WORKER_AND_COMPANION_SHA;POST_PERSISTENCE_HISTORY_IDENTITY_RESOURCE_SIGNAL_GATES',
                preparation_only=True,execution_admitted=False,old_attempt_lock_probe_qualification='BUSY with separately root-admitted corelarge35 stage; old canonical group97558 absent; no lock acquired by preparation')
    if len(plan['jobs'])!=16 or plan['output_limit_bytes']!=8*(1<<30) or plan['worker_rss_limit_bytes']!=2*(1<<30) or plan['internal_floor_bytes']!=3*(1<<30) or plan['ssd_floor_bytes']!=5*(1<<30) or plan['per_worker_seconds_limit']!=7200:raise RuntimeError('CANONICAL_CONTROL_OR_RESOURCE_RULE_CHANGED')
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'ssd_free_bytes':shutil.disk_usage('/Volumes/Extreme SSD').free,
                               'all_retained_canonical_epoch_bytes':canonical_epoch_bytes(SSD),'global_campaign_bytes':global_campaign_bytes(SSD)}
    if plan['resource_preflight']['internal_free_bytes']<plan['internal_floor_bytes'] or plan['resource_preflight']['ssd_free_bytes']<plan['ssd_floor_bytes'] or plan['resource_preflight']['all_retained_canonical_epoch_bytes']>plan['output_limit_bytes'] or plan['resource_preflight']['global_campaign_bytes']>plan['global_reservation_bytes']:raise RuntimeError('UNCHANGED_CANONICAL_RESOURCE_GATES_NOT_MET')
    check_hashes(plan['dependencies_sha256']);regular_hashes(preserved);preserved_failed_attempt(plan)
    write_intended_new(PLAN,plan,intended_sha(plan))
    print(json.dumps({'plan':str(PLAN),'sha256':sha(PLAN),'jobs':16,'arithmetic_workers_prespecified':1,'workers_launched':0,'scope':'METHOD_CALIBRATION_ONLY','execution_admitted':False}))


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
    return canonical_epoch_bytes(output_root)


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
    campaign_bytes=global_campaign_bytes(SSD)
    if campaign_bytes>plan['global_reservation_bytes']:violations.append('300GIB_CANONICAL_GLOBAL_CAMPAIGN_GUARD')
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
            audit_paths=[p for p in output_dir.rglob('*') if p.is_file() or p.is_symlink()]
            if job_id=='INDEPENDENT_CAPTURED_ROW_ARITHMETIC':
                report=root/'independent_control_arithmetic.json'
                audit_paths.append(report)
                companion=report.with_name('._'+report.name)
                if companion.exists() or companion.is_symlink():audit_paths.append(companion)
            output_hashes = {str(p): regular_sha(p) for p in audit_paths}
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
        monitor_path=output_dir/'monitor_receipt.json'
        monitor_intended_sha=intended_sha(receipt)
        write_intended_new(monitor_path,receipt,monitor_intended_sha)
    if returncode != 0 or reason:
        raise RuntimeError('CALIBRATION_WORKER_FAILED_PRESERVED_PARTIAL_OUTPUT: ' + job_id)
    return {'path':str(monitor_path),'sha256':monitor_intended_sha}


def regular_hashes(mapping):
    for path,expected in mapping.items():
        file=Path(path)
        if file.is_symlink() or not file.is_file() or sha(file)!=expected:
            raise RuntimeError('FROZEN_REGULAR_CANONICAL_EVIDENCE_CHANGED: '+str(file))
        if any(parent.is_symlink() for parent in file.parents):
            raise RuntimeError('CANONICAL_EVIDENCE_PARENT_IS_SYMLINK')


def worker_directories(out,native_job):
    return [out,out/'weighted_capture'] if native_job is not None else [out]


def global_directories(plan):
    root=Path(plan['ssd_output_root'])
    return [root,root/'tmp',root/'cache',root/'cache/matplotlib',root/'independent_arithmetic_worker']+[Path(j['output_dir']) for j in plan['jobs']]+[Path(j['output_dir'])/'weighted_capture' for j in plan['jobs']]


def freeze_worker_evidence(plan,admission_sha,job_id,command,out,evidence,native_job=None,monitor_identity=None):
    monitor_path=out/'monitor_receipt.json'
    if monitor_identity is None or monitor_identity.get('path')!=str(monitor_path):raise RuntimeError('INTENDED_MONITOR_IDENTITY_REQUIRED_BEFORE_PARSE')
    monitor_map={str(monitor_path):monitor_identity['sha256']};regular_hashes(monitor_map)
    r=json.loads(monitor_path.read_text());regular_hashes(monitor_map)
    if r['returncode']!=0 or r['stop_reason'] is not None or r['process_group_teardown']['remaining_group_members']!=[]:raise RuntimeError('CANONICAL_WORKER_NOT_COMPLETE_AND_REAPED')
    if r['job_id']!=job_id or r['plan_sha256']!=sha(PLAN) or r['admission_sha256']!=admission_sha or r['command']!=command or r['final_guard_violations'] or r['owned_group_disappearance_verified'] is not True or r['bulk_output_hashes_completed_before_final_guards'] is not True or r['output_hash_error'] is not None or r['deferred_termination_requests']:raise RuntimeError('CANONICAL_MONITOR_FIXED_IDENTITY_OR_FINAL_GUARD_DIFFERS')
    mapping=dict(r['all_output_sha256']);regular_hashes(mapping)
    if str(monitor_path) in mapping:raise RuntimeError('CANONICAL_MONITOR_NOT_A_PREEXISTING_WORKER_OUTPUT')
    external=[];deleted=[]
    if native_job is None:
        if job_id!='INDEPENDENT_CAPTURED_ROW_ARITHMETIC':raise RuntimeError('UNREGISTERED_CANONICAL_WORKER')
        external=[Path(plan['ssd_output_root'])/'independent_control_arithmetic.json']
        deleted=[out/'arithmetic_pending.json']
    allowed=worker_directories(out,native_job)
    observed=regular_tree(out,allowed)
    for path in external:
        if not path.is_file():raise RuntimeError('REQUIRED_EXTERNAL_ARITHMETIC_REPORT_MISSING')
        observed.add(str(path))
        companion=path.with_name('._'+path.name)
        if companion.exists() or companion.is_symlink():observed.add(str(companion))
    ordinary={path for path in mapping if not Path(path).name.startswith('._')}
    if any(not Path(path).is_relative_to(out) and Path(path) not in external for path in ordinary):raise RuntimeError('CANONICAL_WORKER_OUTPUT_SCOPE_DIFFERS')
    if {path for path in observed if not Path(path).name.startswith('._')}!=ordinary|{str(monitor_path)}:raise RuntimeError('CANONICAL_WORKER_EXACT_ORDINARY_OUTPUT_INVENTORY_DIFFERS')
    companions={path:digest for path,digest in mapping.items() if Path(path).name.startswith('._')}
    targets=[Path(path) for path in ordinary]+[monitor_path]+[d for d in allowed if d!=out]+deleted
    freeze_companions(targets,companions,allowed_absent_targets=deleted)
    if observed!=ordinary|{str(monitor_path)}|set(companions):raise RuntimeError('CANONICAL_WORKER_EXACT_OUTPUT_INVENTORY_DIFFERS')
    mapping.update(companions);mapping.update(monitor_map);regular_hashes(mapping)
    if native_job is not None:
        capture_path=out/'native_capture_receipt.json'
        if str(capture_path) not in mapping:raise RuntimeError('REQUIRED_NATIVE_CAPTURE_NOT_IN_FROZEN_WORKER_OUTPUTS')
        cap=json.loads(capture_path.read_text());regular_hashes(mapping)
        if cap['job']!=native_job or cap['plan_sha256']!=sha(PLAN):raise RuntimeError('NATIVE_CAPTURE_FIXED_JOB_OR_PLAN_DIFFERS')
    else:
        intent_path=out/'arithmetic_report_intent.json'
        if str(intent_path) not in mapping or str(external[0]) not in mapping:raise RuntimeError('ARITHMETIC_REPORT_AND_INTENT_NOT_FROZEN')
        intent=json.loads(intent_path.read_text());regular_hashes(mapping)
        if intent['schema']!='canonical_arithmetic_report_intent_v4_7' or intent['executor_sha256']!=sha(PACKAGE/'scripts/canonical_calibration_arithmetic_v4_7.py') or intent['plan_sha256']!=sha(PLAN) or intent['admission_sha256']!=admission_sha or intent['report_path']!=str(external[0]) or intent['intended_report_sha256']!=mapping[str(external[0])]:raise RuntimeError('ARITHMETIC_INTENDED_EXTERNAL_REPORT_DIFFERS')
        native_sources={str(Path(j['output_dir'])/'native_capture_receipt.json'):evidence['all_output_sha256'][str(Path(j['output_dir'])/'native_capture_receipt.json')] for j in plan['jobs']}
        if intent['capture_receipt_sha256']!=native_sources:raise RuntimeError('ARITHMETIC_CONSUMED_CAPTURE_RECEIPTS_DIFFERS')
        binding={'plan_sha256':sha(PLAN),'admission_sha256':admission_sha,'executor_sha256':sha(PACKAGE/'scripts/canonical_calibration_arithmetic_v4_7.py')}
        require_committed(out/'arithmetic_pending.json',out/'arithmetic_terminal_seal.json',binding,{str(external[0]):intent['intended_report_sha256']})
        regular_hashes(mapping)
    if job_id in evidence['workers'] or set(mapping)&set(evidence['all_output_sha256']):raise RuntimeError('CANONICAL_DUPLICATE_WORKER_OR_OUTPUT_EVIDENCE')
    evidence['workers'][job_id]={'output_dir':str(out),'sha256':mapping,'registered_directories':list(map(str,allowed)),'external_output_paths':list(map(str,external)),'registered_deleted_markers':list(map(str,deleted))}
    evidence['all_output_sha256'].update(mapping)


def freeze_global_companions(plan,evidence):
    root=Path(plan['ssd_output_root'])
    targets=[root,root/'tmp',root/'cache',root/'cache/matplotlib',root/'independent_arithmetic_worker']+[Path(j['output_dir']) for j in plan['jobs']]
    frozen=evidence.setdefault('global_transport_sha256',{})
    freeze_companions(targets,frozen)
    for path,digest in frozen.items():
        if path in evidence['all_output_sha256'] and evidence['all_output_sha256'][path]!=digest:raise RuntimeError('CANONICAL_GLOBAL_COMPANION_REBOUND')
        evidence['all_output_sha256'].setdefault(path,digest)


def check_frozen_worker_evidence(evidence):
    regular_hashes(evidence['all_output_sha256'])
    for item in evidence['workers'].values():
        root=Path(item['output_dir']);observed=regular_tree(root,item['registered_directories'])
        for path in map(Path,item['external_output_paths']):
            observed.add(str(path));companion=path.with_name('._'+path.name)
            if companion.exists() or companion.is_symlink():observed.add(str(companion))
        if observed!=set(item['sha256']):raise RuntimeError('CANONICAL_FROZEN_WORKER_INVENTORY_CHANGED')


def final_evidence_inventory(plan,evidence):
    root=Path(plan['ssd_output_root']);check_frozen_worker_evidence(evidence)
    expected_jobs={job['job_id'] for job in plan['jobs']}|{'INDEPENDENT_CAPTURED_ROW_ARITHMETIC'}
    captures={str(Path(job['output_dir'])/'native_capture_receipt.json') for job in plan['jobs']}
    monitors={str(Path(job['output_dir'])/'monitor_receipt.json') for job in plan['jobs']}|{str(root/'independent_arithmetic_worker/monitor_receipt.json')}
    if len(plan['jobs'])!=16 or len(captures)!=16 or len(monitors)!=17 or set(evidence['workers'])!=expected_jobs:raise RuntimeError('EXACT16_NATIVE_AND1_ARITHMETIC_WORKER_CARDINALITY_REQUIRED')
    observed=regular_tree(root,global_directories(plan))
    companion=root.with_name('._'+root.name)
    if companion.exists() or companion.is_symlink():observed.add(str(companion))
    if {path for path in observed if Path(path).name=='native_capture_receipt.json'}!=captures or {path for path in observed if Path(path).name=='monitor_receipt.json'}!=monitors:raise RuntimeError('CANONICAL_REQUIRED_CAPTURE_OR_MONITOR_INVENTORY_DIFFERS')
    if not (captures|monitors)<=set(evidence['all_output_sha256']):raise RuntimeError('CANONICAL_REQUIRED_EVIDENCE_NOT_FROZEN')
    report=root/'independent_control_arithmetic.json'
    if str(report) not in evidence['all_output_sha256']:raise RuntimeError('CANONICAL_ARITHMETIC_REPORT_NOT_FROZEN')
    if observed!=set(evidence['all_output_sha256']):raise RuntimeError('CANONICAL_GLOBAL_FROZEN_OUTPUT_INVENTORY_CHANGED')

def execute(lock_fd, ownership, evidence):
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
        command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_capture_v4_7.py'),
                   '--plan', str(PLAN), '--admission', str(ADMISSION), '--plan-sha', sha(PLAN), '--job-id', job['job_id'], '--shared-lock-fd', str(lock_fd)]
        monitor_identity=monitored_worker(plan, admission_sha, job['job_id'], command, out, lock_fd, ownership)
        freeze_worker_evidence(plan,admission_sha,job['job_id'],command,out,evidence,job,monitor_identity)
        freeze_global_companions(plan,evidence)
        print(json.dumps({'completed_method_control': job['job_id'], 'biological_outcomes_computed': 0}), flush=True)
    verify(lock_fd, ownership, evidence)
    assert_no_termination()


def verify(lock_fd, ownership, evidence):
    assert_no_termination()
    plan, _ = admit()
    check_frozen_worker_evidence(evidence)
    for job in plan['jobs']:
        out = Path(job['output_dir'])
        regular_hashes(evidence['workers'][job['job_id']]['sha256'])
        r = json.loads((out / 'monitor_receipt.json').read_text())
        regular_hashes(evidence['workers'][job['job_id']]['sha256'])
        if r['returncode'] != 0 or r['stop_reason'] is not None or r['process_group_teardown']['remaining_group_members'] != []:
            raise RuntimeError('ALL_16_NATIVE_CONTROLS_MUST_COMPLETE_BEFORE_ARITHMETIC')
        if r['job_id'] != job['job_id'] or r['plan_sha256'] != sha(PLAN) or r['admission_sha256'] != sha(ADMISSION) or r['final_guard_violations'] or r['owned_group_disappearance_verified'] is not True or r['bulk_output_hashes_completed_before_final_guards'] is not True or r['output_hash_error'] is not None or r['deferred_termination_requests']:
            raise RuntimeError('CONTROL_MONITOR_IDENTITY_OR_FINAL_GUARD_MISMATCH')
        check_hashes(r['all_output_sha256'])
        if not (out / 'native_capture_receipt.json').exists():
            raise RuntimeError('MISSING_NATIVE_CAPTURE')
    out = Path(plan['ssd_output_root']) / 'independent_arithmetic_worker'
    out.mkdir(exist_ok=False)
    command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_arithmetic_v4_7.py'), '--plan', str(PLAN), '--admission', str(ADMISSION), '--shared-lock-fd', str(lock_fd)]
    monitor_identity=monitored_worker(plan, sha(ADMISSION), 'INDEPENDENT_CAPTURED_ROW_ARITHMETIC', command, out, lock_fd, ownership)
    freeze_worker_evidence(plan,sha(ADMISSION),'INDEPENDENT_CAPTURED_ROW_ARITHMETIC',command,out,evidence,monitor_identity=monitor_identity)
    freeze_global_companions(plan,evidence)
    report=Path(plan['ssd_output_root'])/'independent_control_arithmetic.json'
    report_map={str(report):evidence['all_output_sha256'][str(report)]};regular_hashes(report_map)
    arithmetic=json.loads(report.read_text());regular_hashes(report_map)
    if arithmetic['all_checks_pass'] is not True or arithmetic['scope']!='METHOD_CALIBRATION_IMPLEMENTATION_CHECK_ONLY' or arithmetic['allow_41_covariance_outcomes'] or arithmetic['calibrated_biological_p_values_computed']:raise RuntimeError('CANONICAL_ARITHMETIC_REPORT_NOT_PASS_OR_SCOPE_CHANGED')
    check_frozen_worker_evidence(evidence)


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
                pending=SSD.parent/'native_canonical_calibration_pending_v4_7.json'
                seal=SSD.parent/'native_canonical_calibration_terminal_seal_v4_7.json'
                initial_identity={'plan_sha256':sha(PLAN),'admission_sha256':sha(ADMISSION),'executor_sha256':sha(Path(__file__))}
                terminal=TerminalCommit(pending,seal,initial_identity)
                evidence={'workers':{},'all_output_sha256':{}}
                execute(lock_fd,ownership,evidence)
                plan,_=admit();root=Path(plan['ssd_output_root'])
                check_frozen_worker_evidence(evidence)
                freeze_global_companions(plan,evidence)
                final_evidence_inventory(plan,evidence)
                output_hashes=dict(evidence['all_output_sha256'])
                historical=__import__('canonical_calibration_common_v4_7').historical_completion(plan)
                result=SSD.parent/'native_canonical_calibration_execution_receipt_v4_7.json'
                record={'status':'ALL16_CANONICAL_SOFTWARE_CONTROLS_AND_CAPTURED_ROW_ARITHMETIC_COMPLETE','plan_sha256':sha(PLAN),'admission_sha256':sha(ADMISSION),'all_output_sha256':output_hashes,'frozen_worker_evidence':evidence['workers'],'historical190_receipt_sha256':historical,'completed_utc':utc(),'biological_outcomes_computed':0,'calibrated_difference_p_values':False}
                result_sha=intended_sha(record)
                expected_seal_sha=terminal_intended_sha(terminal,{str(result):result_sha})
                terminal_transport={}
                write_intended_new(result,record,result_sha)
                def terminal_transport_gate():
                    if seal.exists() or seal.is_symlink():regular_hashes({str(seal):expected_seal_sha})
                    freeze_companions([pending,result,seal],terminal_transport)
                def final_identity():
                    if sha(PLAN)!=initial_identity['plan_sha256'] or sha(ADMISSION)!=initial_identity['admission_sha256'] or sha(Path(__file__))!=initial_identity['executor_sha256']:raise RuntimeError('CANONICAL_TERMINAL_INITIAL_PLAN_OR_ADMISSION_OR_EXECUTOR_CHANGED')
                    if ownership:raise RuntimeError('CANONICAL_OWNED_GROUPS_NOT_EMPTY')
                    current,_=admit()
                    if current!=plan or __import__('canonical_calibration_common_v4_7').historical_completion(plan)!=historical:raise RuntimeError('CANONICAL_TERMINAL_HISTORY_OR_PLAN_CHANGED')
                    regular_hashes(output_hashes);regular_hashes({str(result):result_sha})
                    terminal_transport_gate()
                    freeze_global_companions(plan,evidence)
                    if evidence['all_output_sha256']!=output_hashes:raise RuntimeError('CANONICAL_NEW_GLOBAL_COMPANION_AFTER_RESULT_FREEZE')
                    final_evidence_inventory(plan,evidence)
                def final_resource():
                    assert_no_termination()
                    terminal_transport_gate()
                    state=snapshot(root);size=output_bytes(root)
                    reasons=final_violations(plan,state,size,0,initial_identity['admission_sha256'],initial_identity['plan_sha256'])
                    total=global_campaign_bytes(SSD)
                    if reasons or total>plan['global_reservation_bytes']:raise RuntimeError('CANONICAL_POST_PERSISTENCE_RESOURCE_GUARD: '+str(reasons))
                    assert_no_termination()
                if not terminal.commit({str(result):result_sha},identity_gate=final_identity,resource_gate=final_resource,termination_gate=assert_no_termination):raise RuntimeError('CANONICAL_TERMINAL_COMMIT_FAILED_PENDING_REMAINS')


if __name__ == '__main__':
    main()
