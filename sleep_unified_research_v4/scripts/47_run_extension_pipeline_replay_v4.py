#!/usr/bin/env python3
"""Root-admitted, exclusive, deferred-signal historical preprocessing replay."""
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
from extension_replay_common_v4 import sha,write_new,check_bindings,acquisition_receipt_gate,physical_mount,acquisition_operational_gate,checkpoint_binding_gate,acquisition_family_gate
from terminal_commit_common_v2 import TerminalCommit

ROOT=Path(__file__).resolve().parents[2];P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v4'


def load_protected():
    spec=importlib.util.spec_from_file_location('_pinned_deferred_owned_worker',P/'scripts/sensitivity_executor_v4_4.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def namespace_bytes(root):
    return sum(p.stat().st_size for p in Path(root).rglob('*') if p.is_file() and not p.is_symlink())


def no_preserved_outputs_gate(plan):
    for member in plan['members']:
        for key in ['harmonized','harmonization_qc','harmonization_receipt','munged','source_gate_receipt','comparison_receipt']:
            path=Path(member[key])
            try:path.resolve().relative_to(OUT.resolve())
            except ValueError:raise RuntimeError('OUTPUT_ESCAPES_PRIVATE_PIPELINE_NAMESPACE')
            if path.exists() or path.is_symlink():raise RuntimeError('PRIOR_PIPELINE_OUTPUT_PRESERVED_NO_OVERWRITE: '+str(path))
        for key in ['harmonized','munged_prefix']:
            path=Path(member[key])
            if list(path.parent.glob(path.name+'*')):raise RuntimeError('UNSEALED_PIPELINE_PREFIX_PRESERVED')


def admission_gate(plan,admission,plan_sha,executor_sha):
    if admission.get('execution_admitted') is not True or admission.get('plan_sha256')!=plan_sha or admission.get('executor_sha256')!=executor_sha:
        raise RuntimeError('ROOT_OPERATIONAL_ADMISSION_NOT_EXACT')
    expected={m['acquisition_receipt'] for m in plan['members']}
    initial=admission.get('acquisition_receipt_sha256',{})
    if not isinstance(initial,dict) or not set(initial).issubset(expected):raise RuntimeError('INITIAL_ACQUISITION_RECEIPT_MAP_NOT_FROZEN_PANEL_SUBSET')
    if admission.get('checkpoint_receipt_binding_policy')!=plan['checkpoint_receipt_binding_policy']:
        raise RuntimeError('ROOT_CHECKPOINT_POLICY_NOT_EXPLICITLY_ADMITTED')
    binding=plan['acquisition_execution_binding']
    for key,expected in [('acquisition_execution_plan_sha256',plan['acquisition_plan_sha256']),('acquisition_executor_sha256',binding['executor_sha256']),('acquisition_root_admission_sha256',binding['sha256'][binding['root_admission']])]:
        if admission.get(key)!=expected:raise RuntimeError('ROOT_ACQUISITION_OPERATIONAL_BINDING_DIFFERS: '+key)
    if not admission.get('independent_review_sha256'):raise RuntimeError('INDEPENDENT_ROOT_REVIEW_BINDING_REQUIRED')
    for path,h in admission['independent_review_sha256'].items():
        if sha(path)!=h:raise RuntimeError('ROOT_OPERATIONAL_REVIEW_CHANGED')
    for path,h in initial.items():
        if sha(path)!=h:raise RuntimeError('INITIAL_SUCCESSFUL_ACQUISITION_RECEIPT_CHANGED')
    acquisition_operational_gate(plan)


def seal_master(protected,master,receipt_path,plan,started,ownership,terminal,identity_gate):
    """All success receipts remain provisional until durable terminal commit."""
    try:
        if master['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS':
            master['resource_final']=protected.stage_resource_gate(plan,started,'BEFORE_PIPELINE_PROVISIONAL_RECEIPT')
        write_new(receipt_path,master)
        if master['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS':
            if terminal is None or not ownership[0]:raise RuntimeError('TERMINAL_PENDING_OR_OWNED_CLEANUP_UNCLEARED')
            if not terminal.commit({str(receipt_path):sha(receipt_path)},identity_gate=identity_gate,
                    resource_gate=lambda:protected.stage_resource_gate(plan,started,'POST_PERSISTENCE_PIPELINE_TERMINAL'),
                    termination_gate=lambda:protected.check_termination('POST_PERSISTENCE_PIPELINE_TERMINAL')):
                raise RuntimeError('PIPELINE_TERMINAL_COMMIT_FAILED_PENDING_REMAINS')
    except BaseException as e:
        master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        master['receipt_or_postseal_error']=type(e).__name__+': '+str(e)
        master['termination_requests']=list(protected.TERMINATION_REQUEST)
        try:write_new(Path(str(receipt_path)+'.failure.json'),master)
        except BaseException as secondary:
            protected.safe_print(dict(event='PIPELINE_FAILURE_RECEIPT_WRITE_FAILED',primary_error=str(e),fallback_error=str(secondary),owned_cleanup_verified=ownership[0],pending_authoritative_veto=True),flush=True)


def wait_family_commit(plan,protected,started,expected_receipts):
    """Wait unlocked for the producer's final durable two-copy commitment."""
    wait_started=time.monotonic()
    while True:
        protected.stage_resource_gate(plan,started,'WAITING_FOR_FULL100_SOURCE_FAMILY_COMMIT')
        result=acquisition_family_gate(plan,expected_receipts,require_complete=False)
        if result is not None:return result
        if time.monotonic()-wait_started>plan['guard']['family_terminal_assembly_seconds']:
            raise RuntimeError('SOURCE_FAMILY_FINAL_TERMINAL_SEAL_NOT_COMPLETE')
        protected.check_termination('BEFORE_SOURCE_FAMILY_TERMINAL_WAIT_SLEEP')
        time.sleep(plan['guard']['checkpoint_poll_seconds'])


def wait_checkpoint(plan,member,a,admission,protected,started):
    """Wait for the next frozen source without owning the heavy-worker mutex."""
    receipt=Path(member['acquisition_receipt']);polls=0;assembly_start=None
    while True:
        protected.stage_resource_gate(plan,started,'WAITING_FOR_FROZEN_SOURCE_CHECKPOINT')
        if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('PLAN_OR_ROOT_ADMISSION_CHANGED_WHILE_WAITING')
        if polls%12==0:
            check_bindings(plan);admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
            protected.safe_print(dict(waiting_for_acquisition_source=member['index'],trait=member['extension_trait_id'],heavy_mutex_held=False),flush=True)
        acquisition_family_gate(plan,require_complete=False)
        if receipt.exists():
            try:r=json.loads(receipt.read_text())
            except json.JSONDecodeError:
                if assembly_start is None:assembly_start=time.monotonic()
                if time.monotonic()-assembly_start>plan['guard']['receipt_assembly_seconds']:raise RuntimeError('SOURCE_RECEIPT_INCOMPLETE_PRESERVED_NO_RETRY')
            else:
                if r.get('status')!='EXACT_IMMUTABLE_SOURCE_ACQUIRED':raise RuntimeError('SOURCE_ACQUISITION_CHECKPOINT_FAILED_PRESERVED_NO_RETRY')
                digest=sha(receipt)
                fixed=admission.get('acquisition_receipt_sha256',{}).get(str(receipt))
                if fixed is not None and fixed!=digest:raise RuntimeError('ROOT_PREBOUND_ACQUISITION_RECEIPT_CHANGED')
                acquisition_operational_gate(plan);acquisition_receipt_gate(plan,member,digest)
                if sha(receipt)!=digest:raise RuntimeError('DISCOVERED_SOURCE_RECEIPT_NOT_IMMUTABLE')
                return dict(plan_sha256=a.plan_sha,pipeline_root_admission_sha256=a.admission_sha,extension_trait_id=member['extension_trait_id'],acquisition_receipt=str(receipt),acquisition_receipt_sha256=digest,acquisition_execution_binding_sha256=plan['acquisition_execution_binding']['sha256'],recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        polls+=1;protected.check_termination('BEFORE_WAITING_FOR_CHECKPOINT_SLEEP');time.sleep(plan['guard']['checkpoint_poll_seconds'])


def acquire_source_mutex(plan,protected,started):
    """Busy workers cause bounded unlocked waiting, never a second worker."""
    while True:
        protected.stage_resource_gate(plan,started,'WAITING_FOR_SHARED_HEAVY_MUTEX')
        fd=os.open(plan['guard']['shared_heavy_worker_lock'],os.O_RDWR|os.O_CREAT,0o600)
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(fd);protected.check_termination('BEFORE_HEAVY_MUTEX_WAIT_SLEEP');time.sleep(plan['guard']['checkpoint_poll_seconds']);continue
        return fd


def run(a):
    if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('PLAN_OR_ROOT_ADMISSION_CHANGED')
    plan=json.loads(a.plan.read_text());admission=json.loads(a.admission.read_text())
    admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
    if plan['guard']['shared_heavy_worker_lock']!=str(SSD/'native_heavy_worker.lock'):raise RuntimeError('SHARED_HEAVY_MUTEX_PATH_DIFFERS')
    protected=load_protected();protected.SSD=OUT
    prior={getattr(signal,n):signal.signal(getattr(signal,n),protected.catchable_termination) for n in ['SIGINT','SIGTERM','SIGHUP']}
    ownership=[True,[]];fd=None;master=dict(status='FAILED_PRESERVED_NO_AUTOMATIC_RETRY',plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,completed_members=[],discovered_acquisition_receipt_sha256={},checkpoint_binding_sha256={},estimator_calls=0)
    started=time.monotonic();active_worker=[None];terminal=None
    def limits(_plan,family_started,rss):
        g=plan['guard'];internal=shutil.disk_usage('/System/Volumes/Data').free;free=shutil.disk_usage(OUT).free
        total=namespace_bytes(OUT);global_bytes=namespace_bytes(SSD)
        state=dict(internal_free_bytes=internal,SSD_free_bytes=free,observed_worker_RSS_bytes=rss,pipeline_namespace_bytes=total,new_campaign_namespace_bytes=global_bytes,elapsed_stage_seconds=time.monotonic()-family_started,elapsed_worker_seconds=0 if active_worker[0] is None else time.monotonic()-active_worker[0])
        reason=None
        for condition,label in [(internal<g['internal_floor_bytes'],'INTERNAL_SPACE_GUARD'),(free<g['SSD_floor_bytes'],'SSD_SPACE_GUARD'),(rss>g['observed_aggregate_worker_RSS_limit_bytes'],'AGGREGATE_WORKER_RSS_GUARD'),(total>g['new_output_limit_bytes'],'32GIB_PIPELINE_OUTPUT_GUARD'),(global_bytes>g['global_reservation_bytes'],'300GIB_NEW_CAMPAIGN_RESERVATION_GUARD'),(state['elapsed_stage_seconds']>g['deadline_seconds'],'96H_PIPELINE_DEADLINE'),(state['elapsed_worker_seconds']>g['per_worker_deadline_seconds'],'2H_WORKER_DEADLINE')]:
            if condition:reason=label;break
        return state,reason
    def final_identity_gate():
        if not ownership[0] or ownership[1]:raise RuntimeError('PIPELINE_OWNED_WORKERS_NOT_CLEANED')
        physical_mount();check_bindings(plan,include_archived=True)
        if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('FINAL_PLAN_OR_ADMISSION_CHANGED')
        admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
        evidence=acquisition_family_gate(plan,master['discovered_acquisition_receipt_sha256'],require_complete=True)
        if evidence!=master['acquisition_family_terminal_evidence']:raise RuntimeError('FULL100_FAMILY_TERMINAL_EVIDENCE_CHANGED')
        if protected.baseline_gate(plan)!=master['historical190_gate_receipts']:raise RuntimeError('FINAL190_BASELINE_CHANGED')
        for member in plan['members']:
            bp=OUT/'receipts'/(member['extension_trait_id']+'.checkpoint_binding.json')
            if sha(bp)!=master['checkpoint_binding_sha256'][str(bp)]:raise RuntimeError('FINAL_SOURCE_CHECKPOINT_BINDING_CHANGED')
            checkpoint_binding_gate(plan,member,json.loads(bp.read_text()),a.plan_sha,a.admission_sha)
        for completed in master['completed_members']:
            member=next(m for m in plan['members'] if m['extension_trait_id']==completed['trait'])
            if sha(member['comparison_receipt'])!=completed['comparison_sha256'] or sha(member['munged'])!=completed['new_munged_sha256']:raise RuntimeError('FINAL_REPLAY_OUTPUT_IDENTITY_CHANGED')
    protected.limits=limits
    spec=importlib.util.spec_from_file_location('_unchanged_owned_monitor',P/'scripts/30_prepare_and_run_ssd_native_campaign.py')
    monitor=importlib.util.module_from_spec(spec);spec.loader.exec_module(monitor)
    try:
        protected.check_termination('BEFORE_PIPELINE_PREFLIGHT')
        master['physical_SSD_preflight']=physical_mount();check_bindings(plan,include_archived=True);no_preserved_outputs_gate(plan)
        master['historical190_gate_receipts']=protected.baseline_gate(plan)
        protected.stage_resource_gate(plan,started,'AFTER_BASELINE_AND_PIPELINE_PREFLIGHT')
        attempt=OUT/'pipeline_attempt_v4.json'
        if attempt.exists():raise RuntimeError('PRIOR_PIPELINE_ATTEMPT_PRESERVED_NO_AUTOMATIC_RETRY')
        write_new(attempt,dict(owner_pid=os.getpid(),plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        terminal=TerminalCommit(OUT/'pipeline_pending_v4.json',OUT/'pipeline_terminal_seal_v4.json',
                {'plan_sha256':a.plan_sha,'admission_sha256':a.admission_sha,'executor_sha256':sha(Path(__file__))})
        for member in plan['members']:
            t=member['extension_trait_id'];binding=wait_checkpoint(plan,member,a,admission,protected,started)
            binding_path=OUT/'receipts'/(t+'.checkpoint_binding.json');write_new(binding_path,binding);binding_sha=sha(binding_path)
            master['discovered_acquisition_receipt_sha256'][member['acquisition_receipt']]=binding['acquisition_receipt_sha256']
            master['checkpoint_binding_sha256'][str(binding_path)]=binding_sha
            protected.stage_resource_gate(plan,started,'AFTER_IMMUTABLE_CHECKPOINT_BINDING_BEFORE_MUTEX')
            fd=acquire_source_mutex(plan,protected,started)
            try:
                checkpoint_binding_gate(plan,member,binding,a.plan_sha,a.admission_sha);check_bindings(plan)
                protected.stage_resource_gate(plan,started,'AFTER_SOURCE_CHECKPOINT_IDENTITY_IO_BEFORE_COMMANDS')
                def validation(mode,out):
                    return [plan['python'],'-u','-B',str(P/'scripts/extension_replay_validate_v4.py'),'--plan',str(a.plan),'--plan-sha',a.plan_sha,'--admission',str(a.admission),'--admission-sha',a.admission_sha,'--checkpoint-binding',str(binding_path),'--checkpoint-binding-sha',binding_sha,'--trait',t,'--mode',mode,'--out',out]
                jobs=[('source',validation('source',member['source_gate_receipt']),Path(member['source_gate_receipt'])),('harmonize',member['harmonize_command'],Path(member['harmonized'])),('munge',member['munge_command'],Path(member['munged_prefix'])),('compare',validation('compare',member['comparison_receipt']),Path(member['comparison_receipt']))]
                for label,command,prefix in jobs:
                    protected.check_termination('BEFORE_'+label.upper())
                    if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha or sha(binding_path)!=binding_sha:raise RuntimeError('PLAN_ADMISSION_OR_CHECKPOINT_BINDING_CHANGED_BEFORE_COMMAND')
                    checkpoint_binding_gate(plan,member,binding,a.plan_sha,a.admission_sha)
                    active_worker[0]=time.monotonic()
                    protected.worker(command,prefix,OUT/'receipts_v4'/(t+'__'+label+'.worker.json'),plan,a.plan,a.plan_sha,started,monitor,ownership)
                    check_bindings(plan);admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
                    if sha(binding_path)!=binding_sha:raise RuntimeError('CHECKPOINT_BINDING_CHANGED_AFTER_COMMAND')
                    checkpoint_binding_gate(plan,member,binding,a.plan_sha,a.admission_sha)
                    protected.stage_resource_gate(plan,started,'AFTER_'+label.upper()+'_FINAL_IDENTITY_IO');active_worker[0]=None
                comparison=json.loads(Path(member['comparison_receipt']).read_text())
                if comparison['status']!='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH':raise RuntimeError('FULL_MUNGED_STREAM_COMPARISON_NOT_EXACT')
                master['completed_members'].append(dict(trait=t,acquisition_receipt_sha256=binding['acquisition_receipt_sha256'],checkpoint_binding_sha256=binding_sha,comparison_sha256=sha(member['comparison_receipt']),new_munged_sha256=sha(member['munged']),archived_munged_sha256=member['archived_munged_sha256']))
                protected.stage_resource_gate(plan,started,'AFTER_SOURCE_COMPLETION_HASHES')
                protected.safe_print(dict(completed_sources=len(master['completed_members']),trait=t),flush=True)
            finally:
                # Hold this source's flock through quarantine, even if final
                # metadata queries or source/result hashing have failed.
                recovery=[]
                if not ownership[0]:recovery=protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
                if fd is not None:os.close(fd);fd=None
                active_worker[0]=None
                if recovery:master['recovery_receipt_errors']=recovery;raise RuntimeError('SOURCE_CLEANUP_RECEIPT_FAILURE_PRESERVED')
        master['physical_SSD_final']=physical_mount();check_bindings(plan,include_archived=True)
        if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('FINAL_PLAN_OR_ADMISSION_CHANGED')
        admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
        for member in plan['members']:
            bp=OUT/'receipts'/(member['extension_trait_id']+'.checkpoint_binding.json')
            if sha(bp)!=master['checkpoint_binding_sha256'][str(bp)]:raise RuntimeError('FINAL_SOURCE_CHECKPOINT_BINDING_CHANGED')
            checkpoint_binding_gate(plan,member,json.loads(bp.read_text()),a.plan_sha,a.admission_sha)
        for completed in master['completed_members']:
            member=next(m for m in plan['members'] if m['extension_trait_id']==completed['trait'])
            if sha(member['comparison_receipt'])!=completed['comparison_sha256'] or sha(member['munged'])!=completed['new_munged_sha256']:raise RuntimeError('FINAL_REPLAY_OUTPUT_IDENTITY_CHANGED')
        protected.stage_resource_gate(plan,started,'AFTER_ALL100_FINAL_IDENTITY_IO')
        if len(master['completed_members'])!=100 or len(master['discovered_acquisition_receipt_sha256'])!=100:raise RuntimeError('FINAL_FULL100_CARDINALITY_DIFFERS')
        master['acquisition_family_terminal_evidence']=wait_family_commit(plan,protected,started,master['discovered_acquisition_receipt_sha256'])
        master['status']='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS';master['raw_source_chain_scope']='Original immutable source and unchanged preprocessing replay; exact serialized full-template comparison; no h2/rg fits.'
    except BaseException as e:
        master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY';master['error']=type(e).__name__+': '+str(e)
    finally:
        if not ownership[0]:
            master['recovery_receipt_errors']=protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
            if master['recovery_receipt_errors']:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        master['owned_cleanup_verified']=ownership[0];master['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();master['termination_requests']=list(protected.TERMINATION_REQUEST)
        if protected.TERMINATION_REQUEST:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        try:seal_master(protected,master,OUT/'pipeline_execution_receipt_v4.json',plan,started,ownership,terminal,final_identity_gate)
        finally:
            if fd is not None:os.close(fd)
            for signum,handler in prior.items():signal.signal(signum,handler)
    if master['status']!='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS':raise SystemExit('PIPELINE_STOP_PRESERVED_REQUIRES_REVIEW')


def main():
    p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--plan',type=Path,default=OUT/'extension_pipeline_replay_plan_v4.json');p.add_argument('--plan-sha');p.add_argument('--admission',type=Path);p.add_argument('--admission-sha');a=p.parse_args()
    if not a.execute:raise SystemExit('Prepared only; no worker launched. Root review and explicit execution required.')
    if not a.plan_sha or a.admission is None or not a.admission_sha:raise SystemExit('Exact reviewed plan and root admission hashes required; no worker launched.')
    run(a)


if __name__=='__main__':main()
