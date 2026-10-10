#!/usr/bin/env python3
"""Actual controller with metadata-only fake worker/monitor/lock and resources.

No Popen, process group, actual shared lock, physical mount query or source
parser is used. Pure controller/Terminal2 control flow remains the candidate.
"""
import contextlib
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
sys.path.insert(0,str(S))
import terminal_commit_common_v2 as terminal_helper

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def run():
    started=time.monotonic();out=R/'independent_finngen_controller_contradiction_controls_v4_3_8';out.mkdir(exist_ok=False)
    cases=[];candidate=S/'58_run_finngen_feasibility_stage_v3_8.py';candidate_hash=sha(candidate)
    labels=['healthy', 'result_mutation_during_parse', 'result_changed_at_derivative_hash', 'current_journal_after_stage_write', 'current_worker_after_stage_write', 'source_after_stage_write', 'stage_payload_on_write', 'terminal_seal_payload_on_write', 'postpersist_internal_guard', 'postpersist_deadline_guard', 'postpersist_signal', 'failed_supplemental_writes', 'full_stream_false', 'INFO_claim_true', 'no_retry', 'metadata_errors_nonempty', 'plan_unchanged_false', 'teardown_cleanup_error', 'wrong_command', 'wrong_returncode', 'stop_reason_nonempty', 'worker_mutation_after_parse', 'worker_failure_addendum']
    for label in labels:
        d=out/label;d.mkdir();C=load('_finngen_controller_'+label,candidate)
        originals={n:getattr(C.common,n) for n in ['SHARED_LOCK','exclusive_heavy_lock','write_new','sha']}
        original_save=terminal_helper.save_new
        C.common.TERMINATION_REQUEST.clear();flags={};events=[];cleanup=[0];worker_calls=[0]
        campaign=d/'campaign';approved=campaign/'new_source_feasibility/finngen_R13_F5_INSOMNIA'
        ns=approved/'pipeline_replay_v3_8'
        for folder in ['derived','receipts_v4','logs_v4','tmp','cache']: (ns/folder).mkdir(parents=True,exist_ok=True)
        C.common.SHARED_LOCK=campaign/'native_heavy_worker.lock'
        source=approved/'source_fixture.txt';source.write_text('PRIVATE METADATA FIXTURE ONLY\n')
        source_id=dict(bytes=source.stat().st_size,sha256=sha(source),md5=hashlib.md5(source.read_bytes()).hexdigest())
        derivative=ns/'derived/insomnia.sumstats.gz';result=ns/'derived/insomnia.preprocessing.json'
        stage=d/'stage.json';plan_path=d/'plan.json';admit=d/'admission.json'
        worker_path=ns/'receipts_v4/worker.json';journal=Path(str(worker_path)+'.failure_journal.jsonl')
        stdout=ns/'logs_v4/stdout.log';pending=ns/'pending.json';seal=ns/'terminal.json'
        baseline=d/'baseline.json';save(baseline,dict(fixture_only=True))
        dep=d/'dependency.txt';dep.write_text('PRIVATE BOUND METADATA\n')
        reviews=[]
        for i in range(7):
            path=d/('review_'+str(i)+'.txt');path.write_text('PRIVATE REVIEW FIXTURE\n');reviews.append(str(path))
        plan=dict(stage='preprocessing',scope='FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY',new_rg_commands=0,
            namespace=str(ns),source_namespace=str(approved),source=str(source),source_identity=source_id,
            derivative=str(derivative),preprocessing_receipt=str(result),terminal_pending=str(pending),terminal_seal=str(seal),
            stage_receipt=str(stage),admission=str(admit),executor_path=str(candidate),worker_code='MOCK_WORKER',monitor_code='MOCK_MONITOR',
            global_reservation_bytes=300<<30,baseline_gate_plan=str(baseline),baseline_gate_plan_sha256=sha(baseline),
            shared_heavy_worker_lock=str(C.common.SHARED_LOCK),dependencies_sha256={str(dep):sha(dep),str(candidate):candidate_hash},
            required_independent_review_paths=reviews,guard=dict(worker_count=1,BLAS_threads=1,
                internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,observed_aggregate_worker_RSS_limit_bytes=2<<30,
                new_output_limit_bytes=2<<30,deadline_seconds=7200,poll_seconds=2),
            jobs=[dict(job_id='finngen_insomnia_preprocessing',output_prefix=str(ns/'derived/insomnia'),worker_receipt=str(worker_path),
                result_receipt=str(result),command_template=['MOCK_METADATA_ONLY','--plan',str(plan_path),'--plan-sha256','{PLAN_SHA256}','--inherited-heavy-lock-fd','{HEAVY_LOCK_FD}'])])
        save(plan_path,plan);ph=sha(plan_path)
        save(admit,dict(execution_admitted=True,plan_sha256=ph,scope=plan['scope'],independent_binding_review_pass=True,
            resource_plan_review_pass=True,executor_sha256=candidate_hash,independent_review_sha256={p:sha(p) for p in reviews}))
        def await_cleanup(monitor,ownership,plan_hash):
            cleanup[0]+=1;assert ownership==[True,[]];events.append('MOCK_OWNED_EMPTY_CLEANUP')
        supervisor=SimpleNamespace(baseline_gate=lambda b:{'MOCK_COMPLETE190_METADATA': 'f'*64},await_owned_cleanup=await_cleanup)
        def resource_gate(plan,start,where):
            if C.common.TERMINATION_REQUEST:raise RuntimeError('MOCK_DEFERRED_TERMINATION_AT_'+where)
            state,reason=supervisor.limits(plan,start,0)
            if reason:raise RuntimeError('ACTUAL_CANDIDATE_RESOURCE_LIMIT_'+where+': '+reason)
            return state
        supervisor.stage_resource_gate=resource_gate
        def worker(command,prefix,receipt,plan,path,plan_hash,start,monitor,ownership):
            worker_calls[0]+=1;events.append('MOCK_METADATA_WORKER');assert supervisor.subprocess.fd==917
            assert command[-1]=='917' and command[command.index('--plan-sha256')+1]==ph
            derivative.write_text('PRIVATE DERIVATIVE MARKER; NO GWAS CONTENT\n');stdout.write_text('PRIVATE MOCK STDOUT\n')
            journal.write_text('{"private_mock_worker_only":true}\n')
            record=dict(status='QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS',plan_sha256=ph,retained_minimum_diagnostic_pass=True,
                source_before=source_id,source_after=source_id,derivative=str(derivative),derivative_sha256=sha(derivative),
                full_source_gzip_CRC_and_EOF_verified=label!='full_stream_false',derivative_full_gzip_CRC_and_EOF_verified=True,
                scientific_source_admitted=False,independent_replication_established=False,verified_per_variant_N=False,
                verified_per_variant_INFO=label=='INFO_claim_true')
            save(result,record);outputs={str(p):sha(p) for p in [derivative,result,stdout,journal]}
            save(receipt,dict(status='WORKER_COMPLETE_VERIFIED',plan_sha256=ph,
                command=command+['UNFROZEN_ARGUMENT'] if label=='wrong_command' else command,
                returncode=1 if label=='wrong_returncode' else 0,
                stop_reason='MOCK_NONEMPTY_STOP' if label=='stop_reason_nonempty' else None,
                plan_unchanged=label!='plan_unchanged_false',
                metadata_errors=['MOCK_METADATA_ERROR'] if label=='metadata_errors_nonempty' else [],
                process_group_teardown=dict(remaining_group_members=[],cleanup_error='MOCK_CLEANUP_ERROR' if label=='teardown_cleanup_error' else None),
                owned_cleanup_verified=True,output_sha256=outputs))
            flags['worker_returned']=True
            if label=='worker_failure_addendum':save(Path(str(receipt)+'.failure.json'),dict(status='MOCK_FAILURE_VETO'))
        supervisor.worker=worker
        C.load=lambda n,p:supervisor if p=='MOCK_WORKER' else SimpleNamespace()
        C.physical_mount=lambda:None
        C.shutil=SimpleNamespace(disk_usage=lambda p:SimpleNamespace(free=(2<<30) if flags.get('internal_bad') and str(p)=='/System/Volumes/Data' else 900<<30))
        C.time=SimpleNamespace(monotonic=lambda:time.monotonic()+(7201 if flags.get('deadline_bad') else 0))
        @contextlib.contextmanager
        def mocklock(before_release=None):
            events.append('MOCK_LOCK_ENTER')
            try:yield 917
            finally:
                if before_release:before_release()
                events.append('MOCK_LOCK_EXIT')
        C.common.exclusive_heavy_lock=mocklock
        def write_new(path,data):
            if label=='failed_supplemental_writes' and str(path).endswith('.failure.json'):raise OSError('MOCK_SUPPLEMENTAL_STORAGE_FAILURE')
            originals['write_new'](path,data)
            if Path(path)==stage:
                if label=='stage_payload_on_write':stage.write_text(stage.read_text()+' ')
                if label=='current_journal_after_stage_write':journal.write_text(journal.read_text()+' ')
                if label=='current_worker_after_stage_write':worker_path.write_text(worker_path.read_text()+' ')
                if label=='source_after_stage_write':source.write_text(source.read_text()+' ')
                if label=='postpersist_internal_guard':flags['internal_bad']=True
                if label=='postpersist_deadline_guard':flags['deadline_bad']=True
                if label=='postpersist_signal':C.common.TERMINATION_REQUEST.append(15)
        C.common.write_new=write_new
        def metadata_sha(path):
            h=originals['sha'](path)
            if label=='result_changed_at_derivative_hash' and Path(path)==derivative and flags.get('worker_returned') and not flags.get('mutation'):
                result.write_text(result.read_text()+' ');flags['mutation']=True
            return h
        C.common.sha=metadata_sha
        if label in ['result_mutation_during_parse','worker_mutation_after_parse']:
            def parsed(payload):
                value=json.loads(payload)
                if isinstance(value,dict) and not flags.get('mutation'):
                    if label=='result_mutation_during_parse' and value.get('status')=='QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS':
                        result.write_text(result.read_text()+' ');flags['mutation']=True
                    if label=='worker_mutation_after_parse' and value.get('status')=='WORKER_COMPLETE_VERIFIED':
                        worker_path.write_text(worker_path.read_text()+' ');flags['mutation']=True
                return value
            C.json=SimpleNamespace(loads=parsed,dumps=json.dumps)
        def seal_save(path,data):
            if label=='failed_supplemental_writes' and str(path).endswith('.failure.json'):raise OSError('MOCK_SUPPLEMENTAL_STORAGE_FAILURE')
            original_save(path,data)
            if Path(path)==seal and label in ['terminal_seal_payload_on_write','failed_supplemental_writes']:
                seal.write_text(seal.read_text()+' ')
        terminal_helper.save_new=seal_save
        accepted=False;error=None;consumer=False;retry_rejected=None
        try:
            try:
                C.execute(plan_path,ph);accepted=True
            except BaseException as e:error=type(e).__name__+': '+str(e)
            if stage.exists():
                binding=dict(plan_sha256=ph,admission_sha256=sha(admit),executor_sha256=candidate_hash)
                try:terminal_helper.require_committed(pending,seal,binding,{str(stage):sha(stage)});consumer=True
                except BaseException:pass
            if label=='no_retry':
                try:C.execute(plan_path,ph);retry_rejected=False
                except BaseException:retry_rejected=True
            expected=label in ['healthy','no_retry']
            cases.append(dict(label=label,execute_returned_success=accepted,consumer_accepts_stage=consumer,
                expected_success=expected,passed=accepted==expected and consumer==expected and (retry_rejected is not False),
                exception=error,PENDING_remains=pending.exists(),retry_rejected=retry_rejected,mock_worker_calls=worker_calls[0],
                cleanup_calls=cleanup[0],events=events,actual_workers=0,actual_locks=0,fixture_plan_sha256=ph))
        finally:
            for n,v in originals.items():setattr(C.common,n,v)
            terminal_helper.save_new=original_save;C.common.TERMINATION_REQUEST.clear()
    receipt=dict(schema='independent_actual_finngen_controller_contradiction_mock_controls_v4_3_8',candidate_sha256=candidate_hash,
        candidate_unchanged=sha(candidate)==candidate_hash,controls=cases,all23_expected_controls_pass=all(c['passed'] for c in cases),
        private_regular_file_sha256={str(p):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and not p.is_symlink()},
        no_actual_worker_body_reference_network_mutex=True,scope='Actual58 execute/Terminal2 control flow with fake worker, monitor, lock, physical mount, baseline and disk/time probes only.',
        elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    save(R/'independent_finngen_controller_contradiction_controls_receipt_v4_3_8.json',receipt)
    print(json.dumps(dict(controls=len(cases),all_pass=receipt['all23_expected_controls_pass'],actual_workers=0,actual_locks=0)))

if __name__=='__main__':run()
