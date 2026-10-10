#!/usr/bin/env python3
"""Route one admitted diagnostic through the independently reviewed worker.

The adapter changes only its private output root and passes the already-held
mutex descriptor to its child. It does not invoke sensitivity fits or change
the pinned worker, cleanup, signals, guards or estimator implementation.
"""
import argparse
import importlib.util
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import time
from terminal_commit_common_v2 import TerminalCommit
from extension_replay_common_v4 import physical_mount

import canonical_calibration_common_v4_5 as common


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class InheritedMutexSubprocess:
    """Local module proxy; the original subprocess module is never modified."""
    def __init__(self, fd):
        self.fd = fd
    def Popen(self, *args, **kwargs):
        if 'pass_fds' in kwargs:
            raise RuntimeError('UNEXPECTED_WORKER_DESCRIPTOR_OVERRIDE')
        return subprocess.Popen(*args, **kwargs, pass_fds=(self.fd,))
    def __getattr__(self, name):
        return getattr(subprocess, name)


def require_admission(plan_path, expected_hash):
    if common.sha(plan_path) != expected_hash:
        raise RuntimeError('FROZEN_FINNGEN_PLAN_CHANGED')
    plan = json.loads(plan_path.read_text())
    if plan['scope'] != 'FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY' or plan['new_rg_commands'] != 0:
        raise RuntimeError('FINNGEN_SCOPE_CHANGED')
    if len(plan['jobs']) != 1 or plan['stage'] != 'observed_h2':
        raise RuntimeError('ONE_FROZEN_DIAGNOSTIC_STAGE_REQUIRED')
    namespace = Path(plan['namespace'])
    approved = common.SHARED_LOCK.parent/'new_source_feasibility/finngen_R13_F5_INSOMNIA'
    if namespace != approved/'observed_h2_v3_3' or namespace.is_symlink() or Path(plan['source_namespace'])!=approved or plan['shared_heavy_worker_lock'] != str(common.SHARED_LOCK):
        raise RuntimeError('PRIVATE_FEASIBILITY_NAMESPACE_OR_SHARED_LOCK_CHANGED')
    guard = plan['guard']
    if (guard['worker_count'], guard['BLAS_threads'], guard['internal_floor_bytes'], guard['SSD_floor_bytes'],
        guard['observed_aggregate_worker_RSS_limit_bytes'], guard['new_output_limit_bytes'], guard['deadline_seconds']) != (1,1,3<<30,5<<30,2<<30,2<<30,7200):
        raise RuntimeError('FROZEN_FEASIBILITY_RESOURCE_LIMITS_CHANGED')
    admission_path = Path(plan['admission'])
    a = json.loads(admission_path.read_text())
    if a['execution_admitted'] is not True or a['plan_sha256'] != expected_hash or a['scope'] != plan['scope']:
        raise RuntimeError('EXACT_FINNGEN_EXECUTION_ADMISSION_REQUIRED')
    if a['independent_binding_review_pass'] is not True or a['resource_plan_review_pass'] is not True:
        raise RuntimeError('INDEPENDENT_FINNGEN_REVIEW_NOT_CLEARED')
    reviews = a.get('independent_review_sha256')
    required = plan['required_independent_review_paths']
    if len(required) != 7 or len(set(required)) != 7 or not isinstance(reviews, dict) or set(reviews) != set(required):
        raise RuntimeError('EXACT_NONEMPTY_SCIENTIFIC_AND_OPERATIONAL_REVIEW_BINDINGS_REQUIRED')
    if any(not isinstance(value,str) or len(value)!=64 for value in reviews.values()):
        raise RuntimeError('INVALID_INDEPENDENT_REVIEW_DIGEST')
    if a.get('executor_sha256')!=common.sha(Path(__file__)):raise RuntimeError('FINNGEN_EXACT_EXECUTOR_ADMISSION_REQUIRED')
    common.check_hashes(reviews)
    common.check_hashes(plan['dependencies_sha256'])
    prior_preprocessing_gate(plan)
    gate_path = Path(plan['baseline_gate_plan'])
    if common.sha(gate_path) != plan['baseline_gate_plan_sha256']:
        raise RuntimeError('BASELINE_GATE_PLAN_CHANGED')
    return plan, a, common.sha(admission_path)


def prior_preprocessing_gate(plan):
    from terminal_commit_common_v2 import require_committed
    proof = plan['prior_preprocessing_terminal']
    regular_hashes(proof['metadata_sha256'])
    prior_plan = json.loads(Path(proof['plan_path']).read_text())
    prior = json.loads(Path(proof['stage_receipt']).read_text())
    result = json.loads(Path(proof['preprocessing_receipt']).read_text())
    if prior['status'] != 'QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED' or prior['plan_sha256'] != proof['plan_sha256'] or result['plan_sha256'] != proof['plan_sha256']:
        raise RuntimeError('EXACT_FROZEN_PREPROCESSING_PRODUCER_REQUIRED')
    binding = {'plan_sha256':proof['plan_sha256'], 'admission_sha256':prior['admission_sha256'], 'executor_sha256':proof['executor_sha256']}
    require_committed(prior_plan['terminal_pending'], prior_plan['terminal_seal'], binding,
                      {proof['stage_receipt']:proof['stage_receipt_sha256']})
    if common.sha(prior_plan['terminal_seal']) != proof['terminal_seal_sha256'] or prior['result_sha256'] != proof['result_sha256']:
        raise RuntimeError('FROZEN_PREPROCESSING_TERMINAL_OR_OUTPUT_MAP_CHANGED')
    regular_hashes(proof['result_sha256'])
    if result['status'] != 'QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS' or result['derivative'] != plan['derivative'] or result['derivative_sha256'] != plan['input_sha256'][plan['derivative']]:
        raise RuntimeError('H2_INPUT_NOT_EXACT_COMPLETED_PREPROCESSING_DERIVATIVE')
    if not result['full_source_gzip_CRC_and_EOF_verified'] or not result['derivative_full_gzip_CRC_and_EOF_verified'] or result['scientific_source_admitted'] or result['independent_replication_established']:
        raise RuntimeError('QUALIFIED_PREPROCESSING_SCOPE_OR_FULLSTREAM_PROOF_CHANGED')


def result_gate(plan, plan_hash):
    prior_preprocessing_gate(plan)
    job = plan['jobs'][0]
    path = Path(job['result_receipt'])
    r = json.loads(path.read_text())
    if r['plan_sha256'] != plan_hash:
        raise RuntimeError('FINNGEN_RESULT_PLAN_BINDING_CHANGED')
    if plan['stage'] == 'preprocessing':
        if r['status'] != 'QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS' or not r['retained_minimum_diagnostic_pass']:
            raise RuntimeError('FINNGEN_PREPROCESSING_FAILED_NO_FIT')
        if r['source_before'] != plan['source_identity'] or r['source_after'] != plan['source_identity']:
            raise RuntimeError('FINNGEN_SOURCE_IDENTITY_CHANGED')
        if r['derivative'] != plan['derivative'] or common.sha(plan['derivative']) != r['derivative_sha256']:
            raise RuntimeError('FINNGEN_DERIVATIVE_CHANGED')
        if not r['full_source_gzip_CRC_and_EOF_verified'] or not r['derivative_full_gzip_CRC_and_EOF_verified']:
            raise RuntimeError('FINNGEN_FULL_STREAM_VERIFICATION_ABSENT')
        if r['scientific_source_admitted'] or r['independent_replication_established'] or r['verified_per_variant_N'] or r['verified_per_variant_INFO']:
            raise RuntimeError('FINNGEN_QUALIFICATION_OVERRIDDEN')
    else:
        if len(r['estimates']) != 1 or len(r['final_intersections']) != 1:
            raise RuntimeError('EXACTLY_ONE_STANDALONE_H2_REQUIRED')
        expected = plan['jobs'][0]['ldsc_args']
        if '--h2' not in expected or '--rg' in expected or '--pop-prev' in expected or '--samp-prev' in expected:
            raise RuntimeError('ONLY_OBSERVED_SCALE_H2_DIAGNOSTIC_ADMITTED')
        if r['input_sha256_before'] != plan['input_sha256'] or r['input_sha256_after'] != plan['input_sha256']:
            raise RuntimeError('H2_INPUT_IDENTITY_CHANGED')
    return {str(path):common.sha(path)}


def regular_hashes(mapping):
    for path,expected in mapping.items():
        file=Path(path)
        if file.is_symlink() or not file.is_file() or common.sha(file)!=expected:raise RuntimeError('FINNGEN_FROZEN_REGULAR_OUTPUT_CHANGED: '+path)
        if any(parent.is_symlink() for parent in file.parents):raise RuntimeError('FINNGEN_OUTPUT_PARENT_SYMLINK')


def namespace_bytes(path):
    return sum(p.stat().st_size for p in Path(path).rglob('*') if p.is_file() and not p.is_symlink())


def execute(plan_path, plan_hash):
    plan, admission, admission_hash = require_admission(plan_path, plan_hash)
    fixed_executor=common.sha(Path(__file__))
    supervisor=load('_reviewed_finngen_worker',plan['worker_code'])
    monitor=load('_frozen_finngen_monitor',plan['monitor_code'])
    supervisor.SSD=Path(plan['namespace']);supervisor.TERMINATION_REQUEST=common.TERMINATION_REQUEST
    baseline=json.loads(Path(plan['baseline_gate_plan']).read_text());ownership=[True,[]]
    started=time.monotonic();target=Path(plan['stage_receipt']);job=plan['jobs'][0]
    if target.exists() or target.is_symlink() or Path(job['result_receipt']).exists() or Path(job['result_receipt']).is_symlink():raise RuntimeError('PRIOR_FEASIBILITY_STAGE_PRESERVED_NO_RETRY')
    record=dict(schema='supervised_finngen_source_feasibility_stage_v3',stage=plan['stage'],scope=plan['scope'],plan_sha256=plan_hash,started_utc=common.utc(),status='FAILED_PRESERVED',new_rg_commands=0)
    proofs={};baseline_receipts={};terminal=None
    def limits(_plan,start,rss):
        g=plan['guard'];state={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(common.SHARED_LOCK.parent).free,'observed_worker_RSS_bytes':rss,'Finn_namespace_bytes':namespace_bytes(plan['source_namespace']),'new_campaign_namespace_bytes':namespace_bytes(common.SHARED_LOCK.parent),'elapsed_stage_seconds':time.monotonic()-start}
        for bad,reason in [(state['internal_free_bytes']<g['internal_floor_bytes'],'INTERNAL_SPACE_GUARD'),(state['SSD_free_bytes']<g['SSD_floor_bytes'],'SSD_SPACE_GUARD'),(rss>g['observed_aggregate_worker_RSS_limit_bytes'],'AGGREGATE_FINN_RSS_GUARD'),(state['Finn_namespace_bytes']>g['new_output_limit_bytes'],'2GIB_COMPLETE_FINN_NAMESPACE_GUARD'),(state['new_campaign_namespace_bytes']>plan['global_reservation_bytes'],'300GIB_CAMPAIGN_GUARD'),(state['elapsed_stage_seconds']>g['deadline_seconds'],'2H_FINN_STAGE_DEADLINE')]:
            if bad:return state,reason
        return state,None
    supervisor.limits=limits
    def identity():
        physical_mount()
        source=Path(plan['source'])
        if source.is_symlink() or not source.is_file() or source.stat().st_size!=plan['source_identity']['bytes']:raise RuntimeError('FINNGEN_RAW_SOURCE_NOT_REGULAR_EXACT_SIZE')
        current,_,current_admit=require_admission(plan_path,plan_hash)
        if current!=plan or current_admit!=admission_hash or common.sha(Path(__file__))!=fixed_executor:raise RuntimeError('FINNGEN_FIXED_PLAN_ADMISSION_OR_EXECUTOR_CHANGED')
        if not ownership[0] or ownership[1]:raise RuntimeError('FINNGEN_OWNED_GROUPS_NOT_EMPTY')
        if supervisor.baseline_gate(baseline)!=baseline_receipts:raise RuntimeError('FINNGEN_FULL190_BASELINE_CHANGED')
        regular_hashes(proofs)
        if result_gate(plan,plan_hash)!=result_identity:raise RuntimeError('FINNGEN_RESULT_CONSUMPTION_CHANGED')
        worker_path=Path(job['worker_receipt']);worker=json.loads(worker_path.read_text())
        if worker['status']!='WORKER_COMPLETE_VERIFIED' or worker['plan_sha256']!=plan_hash or worker['process_group_teardown']['remaining_group_members'] or not worker['owned_cleanup_verified']:raise RuntimeError('FINNGEN_WORKER_NOT_COMPLETE_OR_REAPED')
        regular_hashes(worker['output_sha256'])
        if common.sha(plan['source'])!=plan['source_identity']['sha256']:raise RuntimeError('FINNGEN_RAW_CHANGED_AFTER_CONSUMPTION')
    with common.deferred_termination_signals():
        with common.exclusive_heavy_lock(before_release=lambda:supervisor.await_owned_cleanup(monitor,ownership,plan_hash)) as fd:
            supervisor.subprocess=InheritedMutexSubprocess(fd)
            try:
                common.assert_no_termination();physical_mount()
                terminal=TerminalCommit(Path(plan['terminal_pending']),Path(plan['terminal_seal']),{'plan_sha256':plan_hash,'admission_sha256':admission_hash,'executor_sha256':fixed_executor})
                baseline_receipts=supervisor.baseline_gate(baseline)
                supervisor.stage_resource_gate(plan,started,'AFTER_COMPLETE190_PREFLIGHT')
                substitutions={'{HEAVY_LOCK_FD}':str(fd),'{PLAN_SHA256}':plan_hash};command=[substitutions.get(x,x) for x in job['command_template']]
                supervisor.worker(command,Path(job['output_prefix']),Path(job['worker_receipt']),plan,plan_path,plan_hash,started,monitor,ownership)
                result_identity=result_gate(plan,plan_hash);proofs.update(result_identity)
                worker_path=Path(job['worker_receipt']);worker=json.loads(worker_path.read_text());proofs[str(worker_path)]=common.sha(worker_path)
                proofs.update(worker['output_sha256'])
                if plan['stage']=='preprocessing':proofs[plan['derivative']]=common.sha(plan['derivative'])
                regular_hashes(proofs);supervisor.await_owned_cleanup(monitor,ownership,plan_hash);identity()
                supervisor.stage_resource_gate(plan,started,'AFTER_FROZEN_RESULT_HASHES')
                record.update(status='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED',result_sha256=proofs,historical_baseline_receipt_sha256=baseline_receipts,admission_sha256=admission_hash,biological_or_replication_claim_admitted=False)
            except BaseException as error:record['failure']=type(error).__name__+': '+str(error)
            finally:
                supervisor.await_owned_cleanup(monitor,ownership,plan_hash);record.update(completed_utc=common.utc(),owned_cleanup_verified=ownership[0],elapsed_seconds=time.monotonic()-started)
                common.write_new(target,record)
                if record['status']=='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED':
                    try:
                        if terminal is None or not terminal.commit({str(target):common.sha(target)},identity_gate=identity,resource_gate=lambda:supervisor.stage_resource_gate(plan,started,'POST_PERSISTENCE_FINN_TERMINAL'),termination_gate=common.assert_no_termination):raise RuntimeError('FINNGEN_TERMINAL_COMMIT_FAILED')
                    except BaseException as error:
                        record['status']='FAILED_PRESERVED';record['terminal_error']=type(error).__name__+': '+str(error)
                        try:common.write_new(Path(str(target)+'.failure.json'),record)
                        except BaseException:pass
    if record['status']!='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED':raise SystemExit('FINNGEN_STAGE_STOP_PRESERVED_REQUIRES_REVIEW')
    common.safe_diagnostic(lambda:json.dumps(dict(stage=plan['stage'],status=record['status'])))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--expected-plan-sha256', required=True)
    p.add_argument('--execute', action='store_true')
    args = p.parse_args()
    if not args.execute:
        p.error('Explicit --execute and independently bound admission required')
    execute(args.plan, args.expected_plan_sha256)


if __name__ == '__main__':
    main()
