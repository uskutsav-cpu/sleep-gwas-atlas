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
from extension_replay_common import sha,write_new,check_bindings,acquisition_receipt_gate,physical_mount

ROOT=Path(__file__).resolve().parents[2];P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v1'


def load_protected():
    spec=importlib.util.spec_from_file_location('_pinned_deferred_owned_worker',P/'scripts/sensitivity_executor_v4_3_1.py')
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
    if not admission.get('execution_admitted') or admission.get('plan_sha256')!=plan_sha or admission.get('executor_sha256')!=executor_sha:
        raise RuntimeError('ROOT_OPERATIONAL_ADMISSION_NOT_EXACT')
    expected={m['acquisition_receipt'] for m in plan['members']}
    if set(admission.get('acquisition_receipt_sha256',{}))!=expected:
        raise RuntimeError('ALL100_EXACT_ACQUISITION_RECEIPT_SHA_BINDINGS_REQUIRED')
    if not admission.get('independent_review_sha256'):
        raise RuntimeError('INDEPENDENT_ROOT_REVIEW_BINDING_REQUIRED')
    for path,h in admission['independent_review_sha256'].items():
        if sha(path)!=h:raise RuntimeError('ROOT_OPERATIONAL_REVIEW_CHANGED')


def seal_master(protected,master,receipt_path,plan,started,ownership):
    """Receipt/storage/late deferred-signal failures cannot become success."""
    try:
        if master['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS':
            master['resource_final']=protected.stage_resource_gate(plan,started,'AFTER_ALL_PIPELINE_IDENTITY_IO_BEFORE_SUCCESS_SEAL')
        write_new(receipt_path,master)
        protected.check_termination('AFTER_PIPELINE_RECEIPT_SEAL')
    except BaseException as e:
        master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        master['receipt_or_postseal_error']=type(e).__name__+': '+str(e)
        master['termination_requests']=list(protected.TERMINATION_REQUEST)
        try:write_new(Path(str(receipt_path)+'.failure.json'),master)
        except BaseException as secondary:
            protected.safe_print(dict(event='PIPELINE_FAILURE_RECEIPT_WRITE_FAILED',primary_error=str(e),fallback_error=str(secondary),owned_cleanup_verified=ownership[0],record=master),flush=True)
        # The caller finishes only after all owned groups are verified gone,
        # then reports failure; an earlier provisional success is not admitted.


def run(a):
    if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('PLAN_OR_ROOT_ADMISSION_CHANGED')
    plan=json.loads(a.plan.read_text());admission=json.loads(a.admission.read_text())
    admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
    protected=load_protected();protected.SSD=OUT
    prior={getattr(signal,n):signal.signal(getattr(signal,n),protected.catchable_termination) for n in ['SIGINT','SIGTERM','SIGHUP']}
    ownership=[True,[]];fd=None;master=dict(status='FAILED_PRESERVED_NO_AUTOMATIC_RETRY',plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,completed_members=[],estimator_calls=0)
    started=time.monotonic();active_worker=[started]
    def limits(_plan,family_started,rss):
        g=plan['guard'];internal=shutil.disk_usage('/System/Volumes/Data').free;free=shutil.disk_usage(OUT).free
        total=namespace_bytes(OUT);global_bytes=namespace_bytes(SSD)
        state=dict(internal_free_bytes=internal,SSD_free_bytes=free,observed_worker_RSS_bytes=rss,pipeline_namespace_bytes=total,new_campaign_namespace_bytes=global_bytes,elapsed_stage_seconds=time.monotonic()-family_started,elapsed_worker_seconds=time.monotonic()-active_worker[0])
        reason=None
        for condition,label in [(internal<g['internal_floor_bytes'],'INTERNAL_SPACE_GUARD'),(free<g['SSD_floor_bytes'],'SSD_SPACE_GUARD'),(rss>g['observed_aggregate_worker_RSS_limit_bytes'],'AGGREGATE_WORKER_RSS_GUARD'),(total>g['new_output_limit_bytes'],'32GIB_PIPELINE_OUTPUT_GUARD'),(global_bytes>g['global_reservation_bytes'],'300GIB_NEW_CAMPAIGN_RESERVATION_GUARD'),(state['elapsed_stage_seconds']>g['deadline_seconds'],'96H_PIPELINE_DEADLINE'),(state['elapsed_worker_seconds']>g['per_worker_deadline_seconds'],'2H_WORKER_DEADLINE')]:
            if condition:reason=label;break
        return state,reason
    protected.limits=limits
    monitor_spec=importlib.util.spec_from_file_location('_unchanged_owned_monitor',P/'scripts/30_prepare_and_run_ssd_native_campaign.py')
    monitor=importlib.util.module_from_spec(monitor_spec);monitor_spec.loader.exec_module(monitor)
    try:
        fd=os.open(plan['guard']['shared_heavy_worker_lock'],os.O_RDWR|os.O_CREAT,0o600);fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        protected.check_termination('BEFORE_PIPELINE_PREFLIGHT')
        master['physical_SSD_preflight']=physical_mount()
        check_bindings(plan,include_archived=True)
        no_preserved_outputs_gate(plan)
        master['historical190_gate_receipts']=protected.baseline_gate(plan)
        for member in plan['members']:acquisition_receipt_gate(plan,member,admission['acquisition_receipt_sha256'][member['acquisition_receipt']])
        protected.check_termination('AFTER_ALL100_CHECKPOINT_AND_BASELINE_PREFLIGHT')
        protected.stage_resource_gate(plan,started,'AFTER_ALL100_CHECKPOINT_AND_BASELINE_IDENTITY_IO')
        attempt=OUT/'pipeline_attempt_v1.json'
        if attempt.exists():raise RuntimeError('PRIOR_PIPELINE_ATTEMPT_PRESERVED_NO_AUTOMATIC_RETRY')
        write_new(attempt,dict(owner_pid=os.getpid(),plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        for member in plan['members']:
            protected.check_termination('BEFORE_SOURCE_COMMANDS')
            t=member['extension_trait_id'];check_bindings(plan)
            def validation(mode,out):
                return [plan['python'],'-u','-B',str(P/'scripts/extension_replay_validate.py'),'--plan',str(a.plan),'--plan-sha',a.plan_sha,'--admission',str(a.admission),'--admission-sha',a.admission_sha,'--trait',t,'--mode',mode,'--out',out]
            jobs=[('source',validation('source',member['source_gate_receipt']),Path(member['source_gate_receipt'])),
                  ('harmonize',member['harmonize_command'],Path(member['harmonized'])),
                  ('munge',member['munge_command'],Path(member['munged_prefix'])),
                  ('compare',validation('compare',member['comparison_receipt']),Path(member['comparison_receipt']))]
            for label,command,prefix in jobs:
                protected.check_termination('BEFORE_'+label.upper())
                if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('PLAN_OR_ADMISSION_CHANGED_BEFORE_COMMAND')
                active_worker[0]=time.monotonic()
                protected.worker(command,prefix,OUT/'receipts_v4'/(t+'__'+label+'.worker.json'),plan,a.plan,a.plan_sha,started,monitor,ownership)
                protected.check_termination('AFTER_'+label.upper())
                check_bindings(plan)
                admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
                protected.stage_resource_gate(plan,started,'AFTER_'+label.upper()+'_FINAL_IDENTITY_IO')
            comparison=json.loads(Path(member['comparison_receipt']).read_text())
            if comparison['status']!='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH':raise RuntimeError('FULL_MUNGED_STREAM_COMPARISON_NOT_EXACT')
            master['completed_members'].append(dict(trait=t,comparison_sha256=sha(member['comparison_receipt']),new_munged_sha256=sha(member['munged']),archived_munged_sha256=member['archived_munged_sha256']))
            protected.safe_print(dict(completed_sources=len(master['completed_members']),trait=t),flush=True)
        master['physical_SSD_final']=physical_mount()
        check_bindings(plan,include_archived=True)
        if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('FINAL_PLAN_OR_ADMISSION_CHANGED')
        admission_gate(plan,admission,a.plan_sha,sha(Path(__file__)))
        protected.stage_resource_gate(plan,started,'AFTER_FINAL_PIPELINE_HASHES')
        if len(master['completed_members'])!=100:raise RuntimeError('FINAL_PIPELINE_CARDINALITY_OR_RESOURCE_GATE_FAILS')
        master['status']='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS';master['raw_source_chain_scope']='Exact source MD5/SHA/bytes, original harmonic code/local EOF verification, unchanged stock munge and complete decompressed template identity; no h2/rg fits.'
    except BaseException as e:
        master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY';master['error']=type(e).__name__+': '+str(e)
    finally:
        # No resource/hash/persistence/signal failure bypasses owned quarantine.
        if not ownership[0]:
            master['recovery_receipt_errors']=protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
            if master['recovery_receipt_errors']:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        master['owned_cleanup_verified']=ownership[0]
        master['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();master['termination_requests']=list(protected.TERMINATION_REQUEST)
        if protected.TERMINATION_REQUEST:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        try:seal_master(protected,master,OUT/'pipeline_execution_receipt_v1.json',plan,started,ownership)
        finally:
            if fd is not None:os.close(fd)
            for signum,handler in prior.items():signal.signal(signum,handler)
    if master['status']!='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS':raise SystemExit('PIPELINE_STOP_PRESERVED_REQUIRES_REVIEW')


def main():
    p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--plan',type=Path,default=OUT/'extension_pipeline_replay_plan_v1.json');p.add_argument('--plan-sha');p.add_argument('--admission',type=Path);p.add_argument('--admission-sha');a=p.parse_args()
    if not a.execute:raise SystemExit('Prepared only; no worker launched. Root review and explicit execution required.')
    if not a.plan_sha or a.admission is None or not a.admission_sha:raise SystemExit('Exact reviewed plan and root admission hashes required; no worker launched.')
    run(a)


if __name__=='__main__':main()
