#!/usr/bin/env python3
"""Supervise thirteen original local filter replays and full ordered comparisons."""
import argparse
import datetime
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import time
from terminal_commit_common_v2 import TerminalCommit
from extension_replay_common_v4 import physical_mount,write_new,sha

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'validation_pipeline_replay_v2'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(v.is_symlink() for v in p.parents):raise RuntimeError('REGULAR_PIPELINE_EVIDENCE_REQUIRED: '+str(p))
    return p

def frozen_hashes(mapping):
    if not mapping:raise RuntimeError('NONEMPTY_PIPELINE_HASH_MAP_REQUIRED')
    for path,digest in mapping.items():
        if sha(regular(path))!=digest:raise RuntimeError('FROZEN_PIPELINE_EVIDENCE_CHANGED: '+path)

def admission_gate(plan,a):
    if sha(regular(a.plan))!=a.plan_sha or sha(regular(a.admission))!=a.admission_sha:raise RuntimeError('EXACT_PIPELINE_PLAN_AND_ADMISSION_REQUIRED')
    r=json.loads(a.admission.read_text())
    if r['execution_admitted'] is not True or r['plan_sha256']!=a.plan_sha or r['executor_sha256']!=sha(Path(__file__)):raise RuntimeError('EXACT_ROOT_PIPELINE_ADMISSION_REQUIRED')
    required={str(P/'reviews'/name) for name in plan['required_review_filenames']}
    if len(required)!=3 or not required.issubset(r['independent_review_sha256']):raise RuntimeError('INDEPENDENT_WHOLE_PIPELINE_REVIEW_REQUIRED')
    frozen_hashes(r['independent_review_sha256']);frozen_hashes(plan['dependencies_sha256'])
    if plan['member_count']!=13 or len(plan['members'])!=13 or plan['command_count']!=39 or plan['new_estimator_fits']!=0 or plan['jobs'] or plan['new_network_transfers']!=0:
        raise RuntimeError('EXACT_ORIGINAL13_PREPROCESSING_ONLY_FAMILY_REQUIRED')
    if [m['index'] for m in plan['members']]!=list(range(1,14)) or len({m['source_id'] for m in plan['members']})!=13 or plan['adapter_path']!=str(P/'scripts/96_replay_original_validation_collector_v2.py'):
        raise RuntimeError('EXACT_UNIQUE_ORDERED13_MEMBERS_AND_CORRECTED_ADAPTER_REQUIRED')
    for m in plan['members']:
        sid=m['source_id'];row=m['frozen_source_metadata']
        for key,relative in [('new_munged',row['replication_munged_path']),('new_receipt',row['replication_receipt_path'])]:
            q=Path(relative)
            if q.is_absolute() or '..' in q.parts or m[key]!=str(OUT/'workspace'/q):raise RuntimeError('EXACT_ORIGINAL_PRIVATE_MEMBER_ROUTE_REQUIRED')
        targets={'new_qc':OUT/'workspace'/('discovery_extension/results/replication/qc/'+sid+'.tsv'),
            'adapter_receipt':OUT/'receipts_v4'/(sid+'.adapter.json'),'source_gate_receipt':OUT/'receipts_v4'/(sid+'.source.json'),
            'comparison_receipt':OUT/'receipts_v4'/(sid+'.comparison.json')}
        if any(m[k]!=str(v) or any(parent.is_symlink() for parent in v.parents) or v.is_symlink() for k,v in targets.items()):
            raise RuntimeError('EXACT_PRIVATE_REGULAR_RECEIPT_QC_ROUTES_REQUIRED')
    if plan['workspace']!=str(OUT/'workspace') or plan['output_namespace']!=str(OUT) or plan['pending_path']!=str(OUT/'validation_pipeline_pending_v2.json') or plan['terminal_seal_path']!=str(OUT/'validation_pipeline_terminal_v2.json') or plan['primary_receipt_path']!=str(OUT/'validation_pipeline_execution_receipt_v2.json'):
        raise RuntimeError('EXACT_PRIVATE_PIPELINE_ROUTE_REQUIRED')
    g=plan['guard']
    if (g['worker_count'],g['BLAS_threads'],g['internal_floor_bytes'],g['SSD_floor_bytes'],g['observed_aggregate_worker_RSS_limit_bytes'],g['new_output_limit_bytes'],g['global_reservation_bytes'],g['deadline_seconds'],g['per_worker_deadline_seconds'],g['shared_heavy_worker_lock'])!=(1,1,3<<30,5<<30,2<<30,2<<30,300<<30,96*3600,7200,str(SSD/'native_heavy_worker.lock')):
        raise RuntimeError('UNCHANGED_FULL_PIPELINE_RESOURCE_GUARDS_REQUIRED')
    return r

def namespace_bytes(path):return sum(p.stat().st_size for p in Path(path).rglob('*') if p.is_file() and not p.is_symlink())

class InheritedMutexSubprocess:
    def __init__(self,fd):self.fd=fd
    def Popen(self,*args,**kwargs):
        import subprocess
        if 'pass_fds' in kwargs:raise RuntimeError('NO_WORKER_DESCRIPTOR_OVERRIDE')
        return subprocess.Popen(*args,**kwargs,pass_fds=(self.fd,))
    def __getattr__(self,name):
        import subprocess
        return getattr(subprocess,name)

def commands(plan,a,member):
    sid=member['source_id'];validator=P/'scripts/97_validate_original_validation_replay_v2.py'
    def validation(mode,key):return [plan['python'],'-B',str(validator),'--plan',str(a.plan),'--plan-sha256',a.plan_sha,'--source-id',sid,'--mode',mode,'--out',member[key]]
    return [('source',validation('source','source_gate_receipt'),Path(member['source_gate_receipt'])),
            ('collector',[plan['python'],'-B',plan['adapter_path'],'--plan',str(a.plan),'--plan-sha256',a.plan_sha,'--source-id',sid],Path(member['new_munged'])),
            ('compare',validation('compare','comparison_receipt'),Path(member['comparison_receipt']))]

def run(a):
    plan=json.loads(regular(a.plan).read_text());admission_gate(plan,a)
    adapter=module('_validation_local_adapter',plan['adapter_path'])
    validator=module('_validation_full_compare',P/'scripts/97_validate_original_validation_replay_v2.py')
    protected=module('_validation_owned_worker',P/'scripts/sensitivity_executor_v4_4.py');protected.SSD=OUT
    monitor=module('_validation_native_monitor',P/'scripts/30_prepare_and_run_ssd_native_campaign.py')
    handlers={getattr(signal,n):signal.signal(getattr(signal,n),protected.catchable_termination) for n in ['SIGINT','SIGTERM','SIGHUP']}
    ownership=[True,[]];fd=None;active=[None];started=time.monotonic();terminal=None;fixed_primary=None;fixed_seal=None
    proof={};completed=[];workers={};frozen_collector={};baseline=None
    master=dict(schema='original13_validation_pipeline_execution_v2',status='FAILED_PRESERVED_NO_AUTOMATIC_RETRY',
        plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,new_estimator_fits=0,completed_members=completed,worker_receipt_sha256=workers)
    def limits(_plan,start,rss):
        g=plan['guard'];state=dict(internal_free_bytes=shutil.disk_usage('/System/Volumes/Data').free,
            SSD_free_bytes=shutil.disk_usage(SSD).free,observed_worker_RSS_bytes=rss,pipeline_namespace_bytes=sum(namespace_bytes(q) for q in SSD.glob('validation_pipeline_replay_v*') if q.is_dir() and not q.is_symlink()),
            new_campaign_namespace_bytes=namespace_bytes(SSD),elapsed_stage_seconds=time.monotonic()-start,
            elapsed_worker_seconds=0 if active[0] is None else time.monotonic()-active[0])
        for bad,label in [(state['internal_free_bytes']<g['internal_floor_bytes'],'INTERNAL_SPACE_GUARD'),(state['SSD_free_bytes']<g['SSD_floor_bytes'],'SSD_SPACE_GUARD'),(rss>g['observed_aggregate_worker_RSS_limit_bytes'],'2GIB_OWNED_VALIDATION_RSS_GUARD'),(state['pipeline_namespace_bytes']>g['new_output_limit_bytes'],'2GIB_ALL_VALIDATION_PIPELINE_GUARD'),(state['new_campaign_namespace_bytes']>g['global_reservation_bytes'],'300GIB_GLOBAL_NAMESPACE_GUARD'),(state['elapsed_stage_seconds']>g['deadline_seconds'],'96H_FAMILY_DEADLINE'),(state['elapsed_worker_seconds']>g['per_worker_deadline_seconds'],'2H_WORKER_DEADLINE')]:
            if bad:return state,label
        return state,None
    protected.limits=limits
    def identity():
        admission_gate(plan,a);physical_mount();adapter.runtime_gate(plan);adapter.completed_acquisition(plan)
        if ownership[0] is not True or ownership[1]:raise RuntimeError('ALL_OWNED_GROUPS_MUST_BE_REAPED')
        if protected.baseline_gate(plan)!=baseline:raise RuntimeError('ALL190_NATIVE_BASELINES_CHANGED')
        if len(completed)!=13 or len(workers)!=39 or [x['source_id'] for x in completed]!=[m['source_id'] for m in plan['members']]:raise RuntimeError('ALL13_ORDERED_SOURCE_AND39_WORKER_PROOFS_REQUIRED')
        expected_workers={str(OUT/'receipts_v4'/(m['source_id']+'__'+label+'.worker.json')):command
            for m in plan['members'] for label,command,prefix in commands(plan,a,m)}
        if set(workers)!=set(expected_workers):raise RuntimeError('EXACT_REGISTERED39_WORKER_RECEIPTS_REQUIRED')
        frozen_hashes(proof)
        for member in plan['members']:
            validator.current_source(plan,member,adapter)
            record=json.loads(Path(member['comparison_receipt']).read_text())
            if record['status']!='EXACT_ORIGINAL_VALIDATION_ORDERED_CONTENT_QC_REPLAY_PASS' or record['plan_sha256']!=a.plan_sha or record['source_id']!=member['source_id']:
                raise RuntimeError('EXACT_CURRENT_FULL13_COMPARISON_FAMILY_REQUIRED')
            frozen_hashes(frozen_collector[member['source_id']])
            if any(record['consumed_output_sha256'].get(path)!=digest for path,digest in frozen_collector[member['source_id']].items()):
                raise RuntimeError('COMPARISON_MUST_CONSUME_CONTROLLER_FROZEN_COLLECTOR_OUTPUTS')
        for path,digest in workers.items():
            r=json.loads(regular(path).read_text())
            if r['status']!='WORKER_COMPLETE_VERIFIED' or r['command']!=expected_workers[path] or r['plan_sha256']!=a.plan_sha or r['returncode']!=0 or r['stop_reason'] is not None or r['owned_cleanup_verified'] is not True or r['process_group_teardown']['remaining_group_members']:
                raise RuntimeError('EXACT_COMPLETE_REAPED39_WORKERS_REQUIRED')
            required_worker_outputs={str(OUT/'logs_v4'/(Path(path).stem+'.stdout.log')),str(path)+'.failure_journal.jsonl'}
            if not required_worker_outputs.issubset(r['output_sha256']) or r['failure_journal']!=str(path)+'.failure_journal.jsonl':
                raise RuntimeError('EXACT_CURRENT_WORKER_STDOUT_AND_DURABLE_JOURNAL_REQUIRED')
            frozen_hashes(r['output_sha256'])
            if any(proof.get(output)!=h for output,h in r['output_sha256'].items()):raise RuntimeError('COMPLETE_WORKER_OUTPUT_AND_JOURNAL_MAP_REQUIRED')
            if Path(str(path)+'.failure.json').exists() or Path(str(path)+'.failure.json').is_symlink():raise RuntimeError('WORKER_FAILURE_VETOES_FAMILY')
        if fixed_primary is not None:frozen_hashes({plan['primary_receipt_path']:fixed_primary})
        if terminal is not None:
            frozen_hashes({str(terminal.pending):terminal.pending_sha})
            if terminal.seal.exists() or terminal.seal.is_symlink():
                if fixed_seal is None:raise RuntimeError('UNEXPECTED_PRIVATE_TERMINAL_SEAL')
                frozen_hashes({str(terminal.seal):fixed_seal})
    try:
        physical_mount();adapter.runtime_gate(plan);adapter.completed_acquisition(plan);baseline=protected.baseline_gate(plan)
        for member in plan['members']:
            for key in ['new_munged','new_qc','new_receipt','adapter_receipt','source_gate_receipt','comparison_receipt']:
                if Path(member[key]).exists() or Path(member[key]).is_symlink():raise RuntimeError('PRIOR_PIPELINE_ATTEMPT_PRESERVED_NO_RETRY')
        terminal=TerminalCommit(plan['pending_path'],plan['terminal_seal_path'],dict(plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,executor_sha256=sha(Path(__file__))))
        for member in plan['members']:
            while True:
                protected.stage_resource_gate(plan,started,'WAITING_FOR_VALIDATION_HEAVY_MUTEX');protected.check_termination('BEFORE_VALIDATION_MUTEX')
                fd=os.open(plan['guard']['shared_heavy_worker_lock'],os.O_RDWR|os.O_CREAT,0o600)
                try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);break
                except BlockingIOError:os.close(fd);fd=None;time.sleep(plan['guard']['checkpoint_poll_seconds'])
            try:
                protected.subprocess=InheritedMutexSubprocess(fd)
                for label,command,prefix in commands(plan,a,member):
                    admission_gate(plan,a);adapter.runtime_gate(plan);adapter.completed_acquisition(plan);protected.stage_resource_gate(plan,started,'BEFORE_'+label.upper())
                    active[0]=time.monotonic();worker=OUT/'receipts_v4'/(member['source_id']+'__'+label+'.worker.json')
                    protected.worker(command,prefix,worker,plan,a.plan,a.plan_sha,started,monitor,ownership)
                    active[0]=None;workers[str(worker)]=sha(regular(worker));r=json.loads(worker.read_text());proof[str(worker)]=workers[str(worker)]
                    frozen_hashes(r['output_sha256']);proof.update(r['output_sha256'])
                    if label=='collector':
                        frozen_collector[member['source_id']]={member[key]:sha(regular(member[key])) for key in ['new_munged','new_qc','new_receipt','adapter_receipt','source_gate_receipt']}
                        proof.update(frozen_collector[member['source_id']])
                    if label=='compare':frozen_hashes(frozen_collector[member['source_id']])
                    protected.stage_resource_gate(plan,started,'AFTER_'+label.upper())
                r=json.loads(Path(member['comparison_receipt']).read_text())
                if r['status']!='EXACT_ORIGINAL_VALIDATION_ORDERED_CONTENT_QC_REPLAY_PASS':raise RuntimeError('FULL_VALIDATION_CONTENT_QC_COMPARISON_FAILED')
                completed.append(dict(source_id=member['source_id'],comparison_receipt_sha256=sha(member['comparison_receipt'])))
                protected.safe_print(json.dumps(dict(completed_original_validation_sources=len(completed),total=13)),flush=True)
            finally:
                protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
                if ownership[0] is not True or ownership[1]:raise RuntimeError('VALIDATION_CLEANUP_BEFORE_MUTEX_RELEASE_REQUIRED')
                if fd is not None:os.close(fd);fd=None
                active[0]=None
        identity();master.update(status='ALL13_ORIGINAL_VALIDATION_RAW_CONTENT_QC_REPLAY_PASS',all_result_sha256=proof,
            frozen_collector_output_sha256=frozen_collector,historical190_gate_receipts=baseline,original217_classes_unchanged=True,
            original41_native_estimates_unchanged=True,fully_independent_two_trait_replications=0)
    except BaseException as error:master['error']=type(error).__name__+': '+str(error)
    finally:
        protected.await_owned_cleanup(monitor,ownership,a.plan_sha)
        master.update(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),owned_cleanup_verified=ownership[0],termination_requests=list(protected.TERMINATION_REQUEST))
        if ownership[0] is not True or ownership[1] or protected.TERMINATION_REQUEST:master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        try:
            fixed_primary=hashlib.sha256((json.dumps(master,indent=2,allow_nan=False)+'\n').encode()).hexdigest()
            write_new(plan['primary_receipt_path'],master);frozen_hashes({plan['primary_receipt_path']:fixed_primary})
            if master['status']=='ALL13_ORIGINAL_VALIDATION_RAW_CONTENT_QC_REPLAY_PASS':
                intended=dict(status='REVIEWED_STAGE_TERMINAL_SEAL',binding=terminal.binding,result_receipt_sha256={plan['primary_receipt_path']:fixed_primary},pending_path=str(terminal.pending),pending_sha256=terminal.pending_sha,success_requires_absent_PENDING_and_all_failure_addenda=True)
                fixed_seal=hashlib.sha256((json.dumps(intended,indent=2,allow_nan=False)+'\n').encode()).hexdigest()
                if not terminal.commit({plan['primary_receipt_path']:fixed_primary},identity_gate=identity,resource_gate=lambda:protected.stage_resource_gate(plan,started,'POST_PERSISTENCE_VALIDATION_TERMINAL'),termination_gate=lambda:protected.check_termination('FINAL_VALIDATION_TERMINAL')):
                    raise RuntimeError('VALIDATION_PIPELINE_TERMINAL_FAILED_PENDING_REMAINS')
        except BaseException as error:
            master['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
            try:write_new(Path(plan['primary_receipt_path']+'.failure.json'),dict(error=repr(error),status=master['status']))
            except BaseException:pass
        finally:
            if fd is not None:os.close(fd)
            for signum,handler in handlers.items():signal.signal(signum,handler)
    if master['status']!='ALL13_ORIGINAL_VALIDATION_RAW_CONTENT_QC_REPLAY_PASS':raise SystemExit('VALIDATION_PIPELINE_STOP_PRESERVED_REQUIRES_REVIEW')

def main():
    p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha',required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--admission-sha',required=True)
    a=p.parse_args()
    if not a.execute:raise SystemExit('Explicit independently reviewed execution required')
    run(a)

if __name__=='__main__':main()
