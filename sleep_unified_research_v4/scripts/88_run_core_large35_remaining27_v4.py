#!/usr/bin/env python3
"""Admitted 27-member continuation of the preserved original-core preprocessing."""
import argparse
import datetime
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import time
from extension_replay_common_v3 import sha,write_new,check_bindings,physical_mount
from terminal_commit_common_v2 import TerminalCommit

P=Path(__file__).resolve().parents[1];ROOT=P.parent
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'core_pipeline/large35_bounded_replay_v4_remaining27'


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def admission_gate(plan,a):
    if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('CORE_PLAN_OR_ROOT_ADMISSION_CHANGED')
    r=json.loads(a.admission.read_text())
    if r.get('execution_admitted') is not True or r.get('plan_sha256')!=a.plan_sha or r.get('executor_sha256')!=sha(Path(__file__)):
        raise RuntimeError('EXACT_CORE_ROOT_ADMISSION_REQUIRED')
    if not r.get('independent_review_sha256'):raise RuntimeError('CORE_INDEPENDENT_REVIEW_REQUIRED')
    for path,digest in r['independent_review_sha256'].items():
        if sha(path)!=digest:raise RuntimeError('CORE_INDEPENDENT_REVIEW_CHANGED')
    check_bindings(plan)


def runtime_gate(plan):
    p=Path(plan['qualified_runtime_receipt'])
    if sha(p)!=plan['qualified_runtime_receipt_sha256']:raise RuntimeError('QUALIFIED_CORE_RUNTIME_RECEIPT_CHANGED')
    r=json.loads(p.read_text())
    for path,m in r['regular_files'].items():
        f=Path(path)
        if f.is_symlink() or not f.is_file() or f.stat().st_size!=m['bytes'] or sha(f)!=m['sha256']:raise RuntimeError('READ_ONLY_CORE_RUNTIME_FILE_CHANGED: '+path)
    for path,m in r['symlinks'].items():
        f=Path(path)
        if not f.is_symlink() or os.readlink(f)!=m['literal_link'] or str(f.resolve())!=m['resolved_target']:raise RuntimeError('CORE_RUNTIME_LITERAL_SYMLINK_CHANGED')
        if m['target_is_file'] and sha(f.resolve())!=m['resolved_sha256']:raise RuntimeError('CORE_RUNTIME_SYMLINK_TARGET_CHANGED')
    return {'runtime_receipt_sha256':sha(p),'regular_files_verified':len(r['regular_files']),'literal_symlinks_verified':len(r['symlinks']),'historical_per_trait_binary_attestation':False}


def namespace_bytes(path):return sum(p.stat().st_size for p in Path(path).rglob('*') if p.is_file() and not p.is_symlink())


def predecessor_exclusion_gate(plan):
    for path,digest in plan['excluded_v3_evidence_sha256'].items():
        p=Path(path)
        if p.is_symlink() or not p.is_file() or any(parent.is_symlink() for parent in p.parents) or sha(p)!=digest:
            raise RuntimeError('PRESERVED_V3_PREFIX_OR_SNORING_EVIDENCE_CHANGED: '+path)


def acquire_mutex(plan,protected,started):
    while True:
        protected.stage_resource_gate(plan,started,'WAITING_FOR_CORE_HEAVY_MUTEX')
        fd=os.open(plan['guard']['shared_heavy_worker_lock'],os.O_RDWR|os.O_CREAT,0o600)
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(fd);protected.check_termination('BEFORE_CORE_MUTEX_WAIT');time.sleep(plan['guard']['checkpoint_poll_seconds']);continue
        return fd


def run(a):
    if sha(a.plan)!=a.plan_sha:raise RuntimeError('FROZEN_CORE_PLAN_CHANGED')
    plan=json.loads(a.plan.read_text());admission_gate(plan,a)
    if plan['member_count']!=27 or len(plan['members'])!=27 or plan['private_namespace']!=str(OUT) or any(m['prefilter_command'] is not None or m['original_design']['prefilter'] is not None or m['trait_id']=='bmi' for m in plan['members']):raise RuntimeError('EXACT_REMAINING27_PLAN_SCOPE_REQUIRED')
    predecessor_exclusion_gate(plan)
    if plan['guard']['shared_heavy_worker_lock']!=str(SSD/'native_heavy_worker.lock'):raise RuntimeError('EXACT_SHARED_CORE_MUTEX_REQUIRED')
    protected=module('_core_owned_worker',P/'scripts/sensitivity_executor_v4_4.py');protected.SSD=OUT
    monitor=module('_core_native_monitor',P/'scripts/30_prepare_and_run_ssd_native_campaign.py')
    prior={getattr(signal,n):signal.signal(getattr(signal,n),protected.catchable_termination) for n in ['SIGINT','SIGTERM','SIGHUP']}
    ownership=[True,[]];fd=None;terminal=None;started=time.monotonic();active=[None]
    master={'status':'FAILED_PRESERVED_NO_AUTOMATIC_RETRY','plan_sha256':a.plan_sha,'admission_sha256':a.admission_sha,'completed_members':[],'worker_receipt_sha256':{},'expected_worker_commands':{},'generated_metadata_sha256':{},'estimator_calls':0}
    def limits(_plan,start,rss):
        g=plan['guard'];state={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free,
            'observed_worker_RSS_bytes':rss,'core_namespace_bytes':namespace_bytes(SSD/'core_pipeline'),'new_campaign_namespace_bytes':namespace_bytes(SSD),
            'elapsed_stage_seconds':time.monotonic()-start,'elapsed_worker_seconds':0 if active[0] is None else time.monotonic()-active[0]}
        reason=None
        for bad,label in [(state['internal_free_bytes']<g['internal_floor_bytes'],'INTERNAL_SPACE_GUARD'),(state['SSD_free_bytes']<g['SSD_floor_bytes'],'SSD_SPACE_GUARD'),(rss>g['observed_aggregate_worker_RSS_limit_bytes'],'AGGREGATE_CORE_RSS_GUARD'),(state['core_namespace_bytes']>g['new_output_limit_bytes'],'16GIB_ALL_CORE_NAMESPACE_GUARD'),(state['new_campaign_namespace_bytes']>g['global_reservation_bytes'],'300GIB_CAMPAIGN_GUARD'),(state['elapsed_stage_seconds']>g['deadline_seconds'],'96H_CORE_STAGE_DEADLINE'),(state['elapsed_worker_seconds']>g['per_worker_deadline_seconds'],'2H_CORE_WORKER_DEADLINE')]:
            if bad:reason=label;break
        return state,reason
    protected.limits=limits
    def final_identity():
        if not ownership[0] or ownership[1]:raise RuntimeError('CORE_OWNED_GROUPS_NOT_EMPTY')
        admission_gate(plan,a);predecessor_exclusion_gate(plan);physical_mount()
        if runtime_gate(plan)!=master['runtime_before']:raise RuntimeError('READ_ONLY_CORE_RUNTIME_CHANGED_AFTER_EXECUTION')
        if protected.baseline_gate(plan)!=master['historical190_gate_receipts']:raise RuntimeError('CORE_FINAL190_BASELINE_CHANGED')
        if [m['trait'] for m in master['completed_members']]!=[m['trait_id'] for m in plan['members']] or len(master['completed_members'])!=27:raise RuntimeError('ALL27_EXACT_ORDERED_CORE_CONTENT_PROOFS_REQUIRED')
        if len(master['worker_receipt_sha256'])!=135:raise RuntimeError('ALL135_EXACT_CORE_WORKER_RECEIPTS_REQUIRED')
        for item in master['completed_members']:
            member=next(m for m in plan['members'] if m['trait_id']==item['trait'])
            for key in ['harmonized','munged','source_gate_receipt','comparison_receipt','harmonization_qc','bounded_receipt','column_manifest','spool_inventory','spool_cleanup_receipt']:
                if not Path(member[key]).is_file() or Path(member[key]).is_symlink() or any(p.is_symlink() for p in Path(member[key]).parents) or sha(member[key])!=item['output_sha256'][member[key]]:raise RuntimeError('CORE_FINAL_REPLAY_OUTPUT_CHANGED')
            cleanup=json.loads(Path(member['spool_cleanup_receipt']).read_text())
            generated=master['generated_metadata_sha256'][member['trait_id']]
            if item['output_sha256'][member['bounded_receipt']]!=generated[str(Path(member['ephemeral_spool'])/'global_harmonization_candidate_receipt.json')] or item['output_sha256'][member['column_manifest']]!=generated[str(Path(member['ephemeral_spool'])/'columns/column_preparation_manifest.json')]:raise RuntimeError('CORE_FROZEN_GENERATED_METADATA_COPY_CHANGED')
            if cleanup['status']!='OWN_GENERATED_SCRATCH_DISPOSED_AFTER_EXACT_FULL_CONTENT_QC_PASS' or cleanup['plan_sha256']!=a.plan_sha or cleanup['trait']!=member['trait_id'] or not cleanup['spool_absent'] or Path(member['ephemeral_spool']).exists() or Path(member['ephemeral_spool']).is_symlink():raise RuntimeError('CORE_GENERATED_SCRATCH_DISPOSAL_NOT_VERIFIED')
            if cleanup['inventory_sha256']!=item['output_sha256'][member['spool_inventory']] or cleanup['copied_candidate_receipt_sha256']!=item['output_sha256'][member['bounded_receipt']] or cleanup['copied_column_manifest_sha256']!=item['output_sha256'][member['column_manifest']]:raise RuntimeError('CORE_DISPOSAL_METADATA_CHANGED')
            for proof,digest in cleanup['scientific_success_proof_sha256'].items():
                if item['output_sha256'].get(proof)!=digest:raise RuntimeError('CORE_DISPOSAL_SUCCESS_PROOF_DIFFERS')
            for path,digest in item['consumed_prefilter_sha256'].items():
                if sha(path)!=digest:raise RuntimeError('CORE_FINAL_PREFILTER_DERIVATIVE_CHANGED')
            raw=member['original_design']['raw']
            if not Path(raw['resolved_path']).is_file() or Path(raw['resolved_path']).is_symlink() or sha(raw['resolved_path'])!=raw['sealed_verified_sha256']:raise RuntimeError('CORE_FINAL_RAW_SOURCE_CHANGED')
        for path,digest in master['worker_receipt_sha256'].items():
            if Path(path).is_symlink() or not Path(path).is_file() or sha(path)!=digest:raise RuntimeError('CORE_FINAL_OWNED_WORKER_RECEIPT_CHANGED')
            record=json.loads(Path(path).read_text())
            if record.get('command')!=master['expected_worker_commands'][path] or record.get('status')!='WORKER_COMPLETE_VERIFIED' or record.get('plan_sha256')!=a.plan_sha or not record.get('owned_cleanup_verified') or record['process_group_teardown']['remaining_group_members']:
                raise RuntimeError('CORE_FINAL_OWNED_WORKER_NOT_VERIFIED')
            for output,expected in record['output_sha256'].items():
                if Path(output).is_symlink() or not Path(output).is_file() or any(p.is_symlink() for p in Path(output).parents) or sha(output)!=expected:raise RuntimeError('CORE_FINAL_WORKER_OUTPUT_OR_JOURNAL_CHANGED')
    try:
        physical_mount();master['historical190_gate_receipts']=protected.baseline_gate(plan)
        for member in plan['members']:
            if any(p.exists() or p.is_symlink() for p in [Path(member[k]) for k in ['harmonized','munged','source_gate_receipt','comparison_receipt','harmonization_qc','bounded_receipt','column_manifest','spool_inventory','spool_cleanup_receipt']]):raise RuntimeError('PRIOR_CORE_OUTPUT_PRESERVED_NO_RETRY')
        terminal=TerminalCommit(OUT/'core_large35_remaining27_pending_v4.json',OUT/'core_large35_remaining27_terminal_seal_v4.json',{'plan_sha256':a.plan_sha,'admission_sha256':a.admission_sha,'executor_sha256':sha(Path(__file__))})
        master['runtime_before']=runtime_gate(plan)
        protected.stage_resource_gate(plan,started,'AFTER_READ_ONLY_CORE_RUNTIME_PREFLIGHT')
        for member in plan['members']:
            trait=member['trait_id'];fd=acquire_mutex(plan,protected,started)
            try:
                def validation(mode,key):return [plan['harmonization_python'],'-B',str(P/'scripts/85_validate_core_large35_bounded_replay_v3.py'),'--plan',str(a.plan),'--plan-sha',a.plan_sha,'--trait',trait,'--mode',mode,'--out',member[key]]
                jobs=[('source',validation('source','source_gate_receipt'),Path(member['source_gate_receipt']))]
                if member['prefilter_command'] is not None:
                    command=member['prefilter_command'];jobs.append(('prefilter',command,Path(command[command.index('--output')+1])))
                jobs.extend([('harmonize',member['harmonize_command'],Path(member['harmonized'])),('munge',member['munge_command'],Path(member['munged_prefix'])),('compare',validation('compare','comparison_receipt'),Path(member['comparison_receipt']))])
                for label,command,prefix in jobs:
                    admission_gate(plan,a);protected.stage_resource_gate(plan,started,'BEFORE_CORE_'+label.upper())
                    active[0]=time.monotonic()
                    record_path=OUT/'receipts_v4'/(trait+'__'+label+'.worker.json')
                    master['expected_worker_commands'][str(record_path)]=command
                    protected.worker(command,prefix,record_path,plan,a.plan,a.plan_sha,started,monitor,ownership)
                    master['worker_receipt_sha256'][str(record_path)]=sha(record_path)
                    if label=='harmonize':
                        generated=[Path(member['ephemeral_spool'])/'global_harmonization_candidate_receipt.json',Path(member['ephemeral_spool'])/'columns/column_preparation_manifest.json']
                        if any(not f.is_file() or f.is_symlink() or any(parent.is_symlink() for parent in f.parents) for f in generated):raise RuntimeError('GENERATED_HARMONIZER_METADATA_NOT_REGULAR')
                        master['generated_metadata_sha256'][trait]={str(f):sha(f) for f in generated}
                    active[0]=None;admission_gate(plan,a);protected.stage_resource_gate(plan,started,'AFTER_CORE_'+label.upper())
                cleanup_command=[plan['harmonization_python'],'-B',str(P/'scripts/86_cleanup_core_bounded_spool_v3.py'),'--plan',str(a.plan),'--plan-sha',a.plan_sha,'--trait',trait]
                comparison_worker=OUT/'receipts_v4'/(trait+'__compare.worker.json')
                generated=master['generated_metadata_sha256'][trait]
                cleanup_command+=['--expected-candidate-sha256',generated[str(Path(member['ephemeral_spool'])/'global_harmonization_candidate_receipt.json')],
                    '--expected-columns-sha256',generated[str(Path(member['ephemeral_spool'])/'columns/column_preparation_manifest.json')],
                    '--expected-comparison-worker-sha256',master['worker_receipt_sha256'][str(comparison_worker)]]
                admission_gate(plan,a);protected.stage_resource_gate(plan,started,'BEFORE_OWNED_SCRATCH_DISPOSAL')
                active[0]=time.monotonic();cleanup_worker=OUT/'receipts_v4'/(trait+'__scratch_cleanup.worker.json')
                master['expected_worker_commands'][str(cleanup_worker)]=cleanup_command
                protected.worker(cleanup_command,Path(member['spool_cleanup_receipt']),cleanup_worker,plan,a.plan,a.plan_sha,started,monitor,ownership)
                master['worker_receipt_sha256'][str(cleanup_worker)]=sha(cleanup_worker)
                active[0]=None;admission_gate(plan,a);protected.stage_resource_gate(plan,started,'AFTER_OWNED_SCRATCH_DISPOSAL')
                result=json.loads(Path(member['comparison_receipt']).read_text())
                if result.get('status')!='EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS':raise RuntimeError('CORE_FULL_CONTENT_QC_REPLAY_NOT_EXACT')
                outputs={member[k]:sha(member[k]) for k in ['harmonized','munged','source_gate_receipt','comparison_receipt','harmonization_qc','bounded_receipt','column_manifest','spool_inventory','spool_cleanup_receipt']}
                prefilter={}
                if member['original_design']['prefilter'] is not None:
                    command=member['harmonize_command']
                    for flag in ['--infile','--prefilter-provenance']:
                        path=command[command.index(flag)+1];prefilter[path]=sha(path)
                master['completed_members'].append({'trait':trait,'output_sha256':outputs,'consumed_prefilter_sha256':prefilter})
                protected.stage_resource_gate(plan,started,'AFTER_CORE_OUTPUT_HASHES')
                protected.safe_print({'completed_original_core_traits':len(master['completed_members']),'trait':trait},flush=True)
            finally:
                if not ownership[0]:
                    errors=protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
                    if errors:raise RuntimeError('CORE_SOURCE_CLEANUP_RECEIPT_ERRORS: '+str(errors))
                if fd is not None:os.close(fd);fd=None
                active[0]=None
        final_identity();protected.stage_resource_gate(plan,started,'AFTER_ALL_CORE_FINAL_IDENTITY')
        master['status']='ALL27_REMAINING_ORIGINAL_CORE_BOUNDED_HARMONIZED_MUNGED_CONTENT_QC_REPLAY_PASS'
    except BaseException as e:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY';master['error']=type(e).__name__+': '+str(e)
    finally:
        if not ownership[0]:
            master['recovery_receipt_errors']=protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
            if master['recovery_receipt_errors']:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        master['owned_cleanup_verified']=ownership[0];master['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();master['termination_requests']=list(protected.TERMINATION_REQUEST)
        if protected.TERMINATION_REQUEST:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        path=OUT/'core_large35_remaining27_execution_receipt_v4.json'
        try:
            write_new(path,master)
            if master['status']=='ALL27_REMAINING_ORIGINAL_CORE_BOUNDED_HARMONIZED_MUNGED_CONTENT_QC_REPLAY_PASS':
                if terminal is None or not terminal.commit({str(path):sha(path)},identity_gate=final_identity,resource_gate=lambda:protected.stage_resource_gate(plan,started,'POST_PERSISTENCE_CORE_TERMINAL'),termination_gate=lambda:protected.check_termination('POST_PERSISTENCE_CORE_TERMINAL')):
                    raise RuntimeError('CORE_TERMINAL_COMMIT_FAILED_PENDING_REMAINS')
        except BaseException as e:
            master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY';master['terminal_error']=type(e).__name__+': '+str(e)
            try:write_new(Path(str(path)+'.failure.json'),master)
            except BaseException:pass
        finally:
            if fd is not None:os.close(fd)
            for signum,handler in prior.items():signal.signal(signum,handler)
    if master['status']!='ALL27_REMAINING_ORIGINAL_CORE_BOUNDED_HARMONIZED_MUNGED_CONTENT_QC_REPLAY_PASS':raise SystemExit('CORE_REPLAY_STOP_PRESERVED_REQUIRES_REVIEW')


def main():
    p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha',required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--admission-sha',required=True);a=p.parse_args()
    if not a.execute:raise SystemExit('Prepared only; no workers launched.')
    run(a)


if __name__=='__main__':main()
