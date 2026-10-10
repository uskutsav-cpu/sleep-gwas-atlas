#!/usr/bin/env python3
"""Preserve a stopped new attempt and resume only unfinished frozen jobs.

The original runner, invoker, 190-job plan, estimator arguments and floors
remain unchanged. The original runner verifies and skips successful receipts.
No archival research input or historical result is moved or overwritten.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import time

import canonical_calibration_common_v4_3 as common

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
SSD = common.SHARED_LOCK.parent
HISTORICAL = P / 'manifests/ssd_native_execution_plan_v4_3.json'
PLAN = P / 'manifests/native_extension_checkpoint_continuation_v4_3.json'
OUT = SSD / 'native/extension_checkpoint_continuation_v4_3'
PRIOR_MONITOR = P / 'logs/extension_native_monitor_receipt_v4.json'
NEW_MONITOR = P / 'logs/extension_native_monitor_receipt_v4_3.json'
REVIEW_NAMES = (
    'independent_partial_extension_checkpoint_v4_2.tsv',
    'independent_partial_extension_checkpoint_receipt_v4_2.json',
    'independent_partial_extension_checkpoint_review_v4_2.md',
    'independent_partial_extension_checkpoint_checker_v4_2.py',
    'independent_partial_extension_checkpoint_v4_2.sha256',
)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def source_bindings(plan):
    if common.sha(HISTORICAL) != plan['historical_plan_sha256']:
        raise RuntimeError('ORIGINAL_190_PLAN_CHANGED')
    common.check_hashes(plan['dependencies_sha256'])
    load(P/'scripts/extension_replay_common_v2.py', '_continuation_physical_ssd').physical_mount()
    h = json.loads(HISTORICAL.read_text())
    runner = load(ROOT/'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py', '_continuation_original04')
    runner.PACKAGE = Path(h['ssd_support_package'])
    if runner.jobs() != h['jobs'] or runner.execution_dependencies() != common.expected_historical_dependencies(h):
        raise RuntimeError('ORIGINAL_JOBS_OR_DEPENDENCIES_CHANGED')
    return h


def verify_completed(plan, include_new=False):
    h = source_bindings(plan)
    support = Path(h['ssd_support_package'])
    expected_deps = common.expected_historical_dependencies(h)
    identities = {r['path']:r['actual_sha256'] for r in h['inputs_verified'] if r['match']}
    original = load(ROOT/'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py', '_continuation_command04')
    verifier = load(P/'scripts/sensitivity_executor_v4_3_1.py', '_continuation_output_verifier')
    checks = {}
    jobs = plan['completed_jobs'] + (plan['remaining_jobs'] if include_new else [])
    for job in jobs:
        prefix = support/'native/extension_reproduction_v1'/job['job_id']
        rp = Path(str(prefix)+'.execution_receipt.json')
        if str(rp) in plan['completed_receipt_sha256'] and common.sha(rp) != plan['completed_receipt_sha256'][str(rp)]:
            raise RuntimeError('SUCCESSFUL_CHECKPOINT_RECEIPT_CHANGED')
        r = json.loads(rp.read_text())
        command = [str(original.PYTHON),
                   '-u',str(support/'scripts/native_ldsc_capture.py'),'--ldsc-dir',str(ROOT.parent/'ldsc-code'),
                   '--'+job['kind'],','.join(job['inputs']),'--ref-ld-chr',str(original.REF)+'/',
                   '--w-ld-chr',str(original.REF)+'/',
                   '--print-delete-vals','--out',str(prefix)]+job['options']
        inputs = {s:identities[s] for s in job['inputs']}
        if r['job']!=job or r['command']!=command or r['returncode']!=0 or not r['scientific_cardinality_gate_pass'] or not r['execution_identity_gate_pass']:
            raise RuntimeError('CHECKPOINT_JOB_COMMAND_OR_SUCCESS_DIFFERS')
        if r['input_sha256']!=inputs or r['input_sha256_after']!=inputs or r['dependency_sha256_before']!=expected_deps or r['dependency_sha256_after']!=expected_deps:
            raise RuntimeError('CHECKPOINT_INPUT_OR_DEPENDENCY_IDENTITY_DIFFERS')
        verifier.baseline_output_binding_gate(prefix,r)
        if r['all_output_sha256'].get(str(prefix)+'.full_precision.json')!=r['output_sha256']:
            raise RuntimeError('CHECKPOINT_FULL_PRECISION_BINDING_MISSING')
        checks[str(rp)] = common.sha(rp)
    if len(checks)!=(112 if include_new else 107) or sum(j['estimates'] for j in jobs if j['kind']=='rg')!=(1200 if include_new else 700):
        raise RuntimeError('SUCCESSFUL_CHECKPOINT_CARDINALITY_DIFFERS')
    return checks


def prepare():
    if PLAN.exists(): raise RuntimeError('PREPARATION_ALREADY_FROZEN')
    prior = json.loads(PRIOR_MONITOR.read_text())
    if prior['stage']!='extension' or prior['returncode']!=-15 or prior['stop_reason']!='INTERNAL_FULL_NATIVE_FLOOR_REACHED' or prior['process_group_teardown']['remaining_group_members']:
        raise RuntimeError('PRIOR_ACTUAL_RESOURCE_STOP_NOT_CLEARED')
    h = json.loads(HISTORICAL.read_text())
    support = Path(h['ssd_support_package'])
    jobs = [j for j in h['jobs'] if j['stage']=='extension']
    completed, remaining, receipt_hashes = [], [], {}
    for j in jobs:
        rp = support/'native/extension_reproduction_v1'/(j['job_id']+'.execution_receipt.json')
        if rp.exists(): completed.append(j); receipt_hashes[str(rp)]=common.sha(rp)
        else: remaining.append(j)
    if len(completed)!=107 or [j['job_id'] for j in remaining]!=['extension_rg_snoring','extension_rg_sleep_apnea','extension_rg_sleep_efficiency','extension_rg_accel_sleep_duration','extension_rg_sleep_timing']:
        raise RuntimeError('EXACT_REMAINING_FROZEN_QUEUE_DIFFERS')
    prefix = support/'native/extension_reproduction_v1/extension_rg_snoring'
    partials = sorted(p for p in prefix.parent.glob(prefix.name+'*') if not p.name.startswith('._'))
    if len(partials)!=284 or any(p.is_symlink() or not p.is_file() for p in partials):
        raise RuntimeError('UNSEALED_SNORING_PREFIX_SCOPE_DIFFERS')
    sidecars = [q.with_name('._'+q.name) for q in partials]
    if len(sidecars)!=284 or any(not q.is_file() or q.is_symlink() for q in sidecars):
        raise RuntimeError('EXACT_TRANSPORT_SIDECAR_SCOPE_DIFFERS')
    review_files = [P/'reviews'/n for n in REVIEW_NAMES]
    review_files += [P/'scripts/61_native_checkpoint_continuation.py',P/'manifests/native_extension_checkpoint_continuation_v4_2.json',P/'logs/native_sidecar_rename_control_v2.json']
    dependencies = {str(p):common.sha(p) for p in [HISTORICAL,PRIOR_MONITOR,Path(__file__),P/'scripts/31_invoke_original_native_runner.py',P/'scripts/30_prepare_and_run_ssd_native_campaign.py',P/'scripts/sensitivity_executor_v4_3_1.py',P/'scripts/canonical_calibration_common_v4_3.py',P/'scripts/58_run_finngen_feasibility_stage_v2.py',P/'scripts/extension_replay_common_v2.py',*review_files]}
    dependencies[prior['stdout_path']] = prior['stdout_sha256']
    plan = dict(schema='exact_original_extension_checkpoint_continuation_v4_3',stage='extension',historical_plan=str(HISTORICAL),historical_plan_sha256=common.sha(HISTORICAL),dependencies_sha256=dependencies,
                prior_failed_monitor=str(PRIOR_MONITOR),prior_failed_monitor_sha256=common.sha(PRIOR_MONITOR),original190_jobs=h['jobs'],completed_jobs=completed,remaining_jobs=remaining,completed_receipt_sha256=receipt_hashes,
                interrupted_artifacts=[dict(original_path=str(p),preserved_path=str(OUT/'interrupted_snoring'/p.name),bytes=p.stat().st_size,sha256=common.sha(p)) for p in partials],
                interrupted_sidecars=[dict(original_path=str(q),preserved_path=str(OUT/'interrupted_snoring'/q.name),byte_backup_path=str(OUT/'interrupted_snoring_sidecar_original_bytes'/('sidecar_'+str(i).zfill(4)+'.bin')),bytes=q.stat().st_size,sha256=common.sha(q),biological_estimate=False) for i,q in enumerate(sidecars,1)],
                namespace=str(OUT),new_monitor=str(NEW_MONITOR),shared_heavy_worker_lock=str(common.SHARED_LOCK),
                guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,observed_aggregate_worker_RSS_limit_bytes=2<<30,new_output_limit_bytes=4<<30,deadline_seconds=36*3600,poll_seconds=2),
                original_runner_unchanged=True,scientific_arguments_unchanged=True,successful_fits_repeated=False,partial_new_attempt_only_preserved=True,execution_requires_separate_independent_review=True,prepared_utc=common.utc())
    verify_completed(plan)
    write_new(PLAN,plan)
    common.safe_diagnostic(lambda:json.dumps(dict(plan_sha256=common.sha(PLAN),completed_jobs=107,unfinished_jobs=5,workers_launched=0)))


def preserve_partial(plan):
    """Bind both files in each pair; back up literal sidecar bytes before rename."""
    journal = OUT/'interrupted_snoring_preservation.jsonl'
    receipt = OUT/'interrupted_snoring_preservation_receipt.json'
    if journal.exists() or receipt.exists(): raise RuntimeError('PRIOR_PRESERVATION_ATTEMPT_REQUIRES_REVIEW')
    (OUT/'interrupted_snoring').mkdir(parents=True,exist_ok=False)
    (OUT/'interrupted_snoring_sidecar_original_bytes').mkdir(exist_ok=False)
    all_files = plan['interrupted_artifacts'] + plan['interrupted_sidecars']
    if len(plan['interrupted_artifacts'])!=284 or len(plan['interrupted_sidecars'])!=284:
        raise RuntimeError('EXACT_BIOLOGICAL_AND_TRANSPORT_COUNTS_REQUIRED')
    for r in all_files:
        source,target=Path(r['original_path']),Path(r['preserved_path'])
        if source.is_symlink() or not source.is_file() or source.stat().st_size!=r['bytes'] or common.sha(source)!=r['sha256'] or target.exists():
            raise RuntimeError('UNSEALED_PARTIAL_IDENTITY_OR_DESTINATION_DIFFERS')
        if source.stat().st_dev!=target.parent.stat().st_dev: raise RuntimeError('SAME_VOLUME_RENAME_REQUIRED')
    def event(f,name,r):
        f.write(json.dumps(dict(event=name,**r))+'\n');f.flush();os.fsync(f.fileno())
    companions = {r['original_path']:r for r in plan['interrupted_sidecars']}
    with journal.open('x') as f:
        for r in plan['interrupted_sidecars']:
            event(f,'LITERAL_SIDECAR_BACKUP_INTENT',r)
            with Path(r['original_path']).open('rb') as source,Path(r['byte_backup_path']).open('xb') as target:
                for block in iter(lambda:source.read(65536),b''):target.write(block)
                target.flush();os.fsync(target.fileno())
            if Path(r['byte_backup_path']).stat().st_size!=r['bytes'] or common.sha(r['byte_backup_path'])!=r['sha256']:
                raise RuntimeError('SIDECAR_ORIGINAL_BYTE_BACKUP_DIFFERS')
            event(f,'LITERAL_SIDECAR_BACKUP_READBACK_VERIFIED',r)
        for r in plan['interrupted_artifacts']:
            source,target=Path(r['original_path']),Path(r['preserved_path'])
            side=companions[str(source.with_name('._'+source.name))]
            pair=dict(biological=r,transport_sidecar=side)
            event(f,'PAIRED_SAME_VOLUME_PRESERVATION_INTENT',pair)
            if common.sha(source)!=r['sha256'] or common.sha(side['original_path'])!=side['sha256']:
                raise RuntimeError('PAIR_CHANGED_BEFORE_TRANSPORT')
            source.rename(target)
            side_source,side_target=Path(side['original_path']),Path(side['preserved_path'])
            if side_source.exists():
                if side_target.exists():raise RuntimeError('AMBIGUOUS_COMPANION_RENAME_PRESERVED_NO_LAUNCH')
                side_source.rename(side_target)
                transport='EXPLICIT_COMPANION_RENAME'
            else:transport='OBSERVED_IMPLICIT_COMPANION_RENAME'
            if common.sha(target)!=r['sha256'] or not side_target.is_file() or common.sha(side_target)!=side['sha256']:
                raise RuntimeError('PAIRED_TRANSPORT_READBACK_DIFFERS_ORIGINAL_BACKUP_RETAINED')
            event(f,'PAIR_PRESERVED_READBACK_VERIFIED',dict(pair,transport=transport))
    for r in all_files:
        if Path(r['original_path']).exists() or Path(r['preserved_path']).stat().st_size!=r['bytes'] or common.sha(r['preserved_path'])!=r['sha256']:
            raise RuntimeError('FINAL_ALL568_TRANSPORT_IDENTITIES_NOT_EXACT')
    for r in plan['interrupted_sidecars']:
        if common.sha(r['byte_backup_path'])!=r['sha256']:raise RuntimeError('FINAL_SIDECAR_BACKUP_CHANGED')
    write_new(receipt,dict(status='ALL284_BIOLOGICAL_AND284_SIDECAR_ARTIFACTS_PRESERVED_WITHOUT_DELETION',artifacts=plan['interrupted_artifacts'],transport_sidecars=plan['interrupted_sidecars'],sidecars_count_as_biological_estimates=False,journal_sha256=common.sha(journal),rollback='Before any continuation fit exists, restore each biological and paired sidecar to its original absent path; original sidecar literal backups are retained. Require exact size/SHA before and after. Never replace a new continuation output.',completed_utc=common.utc()))
    return {str(receipt):common.sha(receipt),str(journal):common.sha(journal)}


def execute(expected_sha, admission_path, admission_sha):
    if common.sha(PLAN)!=expected_sha or common.sha(admission_path)!=admission_sha: raise RuntimeError('PLAN_OR_ROOT_ADMISSION_CHANGED')
    plan=json.loads(PLAN.read_text());admission=json.loads(admission_path.read_text())
    if admission.get('execution_admitted') is not True or admission.get('plan_sha256')!=expected_sha or admission.get('executor_sha256')!=common.sha(Path(__file__)):
        raise RuntimeError('EXACT_INDEPENDENTLY_REVIEWED_ROOT_ADMISSION_REQUIRED')
    reviews=admission.get('independent_review_sha256',{})
    if not isinstance(reviews,dict) or not any(p.endswith('.json') for p in reviews) or not any(p.endswith('.md') for p in reviews): raise RuntimeError('NONEMPTY_INDEPENDENT_CONTINUATION_REVIEW_REQUIRED')
    common.check_hashes(reviews);verify_completed(plan)
    if OUT.exists() or NEW_MONITOR.exists():raise RuntimeError('PRIOR_CONTINUATION_ATTEMPT_PRESERVED_NO_RETRY')
    monitor=load(P/'scripts/30_prepare_and_run_ssd_native_campaign.py','_continuation_monitor')
    protected=load(P/'scripts/sensitivity_executor_v4_3_1.py','_continuation_supervisor')
    protected.SSD=OUT;protected.TERMINATION_REQUEST=common.TERMINATION_REQUEST
    started=time.monotonic();ownership=[True,[]];peak=[0]
    historical=json.loads(HISTORICAL.read_text());support=Path(historical['ssd_support_package'])
    def limits(_plan,when,rss):
        state=monitor.snapshot();peak[0]=max(peak[0],rss)
        output=sum(p.stat().st_size for root in [support/'native',OUT] for p in root.rglob('*') if p.is_file())
        state.update(observed_worker_RSS_bytes=rss,native_output_bytes=output,elapsed_stage_seconds=time.monotonic()-when)
        reason=monitor.final_limits(state,rss,output,state['elapsed_stage_seconds'])
        return state,reason
    protected.limits=limits
    record=dict(stage='extension',plan_sha256=plan['historical_plan_sha256'],continuation_plan_sha256=expected_sha,admission_sha256=admission_sha,started_utc=common.utc(),returncode=None,stop_reason=None,prior_failed_monitor_sha256=plan['prior_failed_monitor_sha256'],unchanged_successful_jobs_reused=107,remaining_frozen_jobs=5,old_failed_attempt_reclassified=False)
    def identity_gate():
        if common.sha(PLAN)!=expected_sha or common.sha(admission_path)!=admission_sha:
            raise RuntimeError('PLAN_OR_ROOT_ADMISSION_CHANGED_BEFORE_ACTION')
        common.check_hashes(reviews);source_bindings(plan);common.assert_no_termination()
    def failure(message):
        if record['stop_reason'] is None:record['stop_reason']=message
        else:record.setdefault('additional_failures',[]).append(message)
    def failure_addendum():
        try:write_new(Path(str(NEW_MONITOR)+'.failure.json'),record)
        except BaseException as e:
            common.safe_diagnostic(lambda:json.dumps(dict(event='CONTINUATION_FAILURE_RECEIPT_STORAGE_FAILED',record=record,secondary_error=repr(e))))
    with common.deferred_termination_signals():
        with common.exclusive_heavy_lock(before_release=lambda:protected.await_owned_cleanup(monitor,ownership,expected_sha)) as fd:
            proxy=load(P/'scripts/58_run_finngen_feasibility_stage_v2.py','_continuation_fd_proxy')
            protected.subprocess=proxy.InheritedMutexSubprocess(fd)
            try:
                identity_gate()
                protected.stage_resource_gate(plan,started,'BEFORE_CONTINUATION_PRESERVATION')
                common.check_hashes(reviews);verify_completed(plan)
                OUT.mkdir(parents=True,exist_ok=False)
                for name in ['logs_v4','receipts_v4']:(OUT/name).mkdir()
                record['interrupted_preservation_sha256']=preserve_partial(plan)
                identity_gate()
                protected.stage_resource_gate(plan,started,'AFTER_PARTIAL_PRESERVATION_IDENTITY_IO')
                command=[__import__('sys').executable,'-B',str(P/'scripts/31_invoke_original_native_runner.py'),'extension']
                worker_receipt=OUT/'receipts_v4/extension_original_runner_continuation.worker.json'
                protected.worker(command,OUT/'invoker_output',worker_receipt,plan,PLAN,expected_sha,started,monitor,ownership)
                r=json.loads(worker_receipt.read_text())
                if r['status']!='WORKER_COMPLETE_VERIFIED':raise RuntimeError('ORIGINAL_RUNNER_CONTINUATION_FAILED')
                verify_completed(plan,include_new=True)
                paths=[support/'native/extension_reproduction_v1'/(j['job_id']+'.execution_receipt.json') for j in plan['completed_jobs']+plan['remaining_jobs']]
                if len(paths)!=112 or any(not p.is_file() for p in paths):raise RuntimeError('FULL112_EXTENSION_RECEIPTS_MISSING')
                record['native_receipt_sha256']={str(p):common.sha(p) for p in paths}
                record.update(returncode=0,stop_reason=None,worker_receipt=str(worker_receipt),worker_receipt_sha256=common.sha(worker_receipt),process_group_teardown=r['process_group_teardown'])
                source_bindings(plan);common.check_hashes(reviews)
                if common.sha(PLAN)!=expected_sha or common.sha(admission_path)!=admission_sha:raise RuntimeError('TERMINAL_PLAN_OR_ADMISSION_CHANGED')
                record['final_resource_snapshot']=protected.stage_resource_gate(plan,started,'AFTER_ALL_NATIVE_RECEIPT_AND_DEPENDENCY_HASHES')
            except BaseException as e:
                failure(type(e).__name__+': '+str(e))
            finally:
                errors=protected.await_owned_cleanup(monitor,ownership,expected_sha)
                if errors:failure('OWNED_CLEANUP_RECEIPT_ERROR');record['cleanup_errors']=errors
                record.update(completed_utc=common.utc(),elapsed_seconds=time.monotonic()-started,owned_cleanup_verified=ownership[0],peak_observed_owned_rss_bytes=peak[0],termination_requests=list(common.TERMINATION_REQUEST))
                if common.TERMINATION_REQUEST:failure('DEFERRED_TERMINATION_SIGNAL')
                try:write_new(NEW_MONITOR,record)
                except BaseException as e:
                    failure('PRIMARY_MONITOR_SAVE_FAILURE: '+repr(e));failure_addendum()
                try:
                    identity_gate()
                    protected.stage_resource_gate(plan,started,'AFTER_CONTINUATION_MONITOR_PERSISTENCE')
                    common.assert_no_termination()
                except BaseException as e:
                    failure('POST_PERSISTENCE_GATE: '+repr(e));failure_addendum()
    if record['returncode']!=0 or record['stop_reason'] is not None:raise SystemExit('CONTINUATION_STOP_PRESERVED_REQUIRES_REVIEW')
    common.safe_diagnostic(lambda:json.dumps(dict(stage='extension',completed_native_commands=112,unchanged_reused=107,new_completed=5,numerical_review_still_required=True)))


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',action='store_true');g.add_argument('--execute',action='store_true')
    p.add_argument('--plan-sha');p.add_argument('--admission',type=Path);p.add_argument('--admission-sha');a=p.parse_args()
    if a.prepare:prepare()
    elif a.plan_sha and a.admission and a.admission_sha:execute(a.plan_sha,a.admission,a.admission_sha)
    else:p.error('Exact reviewed plan and root admission bindings are required')


if __name__=='__main__':main()
