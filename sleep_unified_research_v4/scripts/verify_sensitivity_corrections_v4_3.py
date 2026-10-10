#!/usr/bin/env python3
"""Read-only bindings and isolated stdlib failure/cleanup controls; no estimators."""
import argparse
import ast
import contextlib
import copy
import datetime
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import resource
import signal
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan-sha',required=True);a=p.parse_args()
    checks=[]
    def check(value,name):
        if not value:raise RuntimeError('CONTROL_FAILED: '+name)
        checks.append(name)
    plan_path=SSD/'sensitivity_operational_plan_v4_3.json';plan=json.loads(plan_path.read_text())
    check(sha(plan_path)==a.plan_sha,'exact_v4_plan_SHA')
    for path,digest in plan['dependencies_sha256'].items():check(sha(path)==digest,'fresh_bound_dependency '+path)
    for path,digest in plan['input_sha256'].items():check(sha(path)==digest,'fresh_target_input '+path)
    old=json.loads(Path(plan['superseded_operational_plan']).read_text())
    check(sha(plan['superseded_operational_plan'])==plan['superseded_operational_plan_sha256'],'sealed_v1_plan_unchanged')
    for key in ['scientific_member_manifest_sha256','input_sha256','guard','baseline_dependency_sha256','original_190_jobs','native_command_count','new_fitted_estimates','stock_merge_only_command_count','stock_merge_only_intersection_identities']:
        check(plan[key]==old[key],'scientific_inputs_membership_guard_unchanged '+key)
    for kind in ['jobs','audit_jobs']:
        for new,prior in zip(plan[kind],old[kind]):
            for key,value in prior.items():
                if key in ['out_prefix','ldsc_args']:continue
                check(new[key]==value,'no_scientific_job_change '+new.get('job_id',new.get('audit_id'))+' '+key)
            i=new['ldsc_args'].index('--out');x=list(new['ldsc_args']);x[i+1]=prior['out_prefix']
            check(x==prior['ldsc_args'],'sole_argv_namespace_change '+new.get('job_id',new.get('audit_id')))
    executor=P/'scripts/sensitivity_executor_v4_3.py'
    spec=importlib.util.spec_from_file_location('_isolated_sensitivity_executor_v4_controls',executor)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    normalized=mod.relocated_baseline_dependencies(plan)
    relocation=plan['baseline_dependency_relocation'];prior=plan['baseline_dependency_sha256']
    check(len(normalized)==len(prior)==56,'56_baseline_dependencies_after_single_relocation')
    check(set(prior)-set(normalized)=={relocation['original_path']} and set(normalized)-set(prior)=={relocation['relocated_path']},'exactly_one_path_key_relocation')
    check(all(normalized[p]==h for p,h in prior.items() if p!=relocation['original_path']),'all_other_55_dependencies_unchanged')
    count=0
    mod.baseline_monitor_binding_gate(json.loads((P/'logs/core_native_monitor_receipt_v4.json').read_text()),'core',plan)
    check(True,'actual_core_monitor_exact_stage_plan_binding')
    for key,value in [('stage','extension'),('plan_sha256','0'*64)]:
        bad=json.loads((P/'logs/core_native_monitor_receipt_v4.json').read_text());bad[key]=value
        try:mod.baseline_monitor_binding_gate(bad,'core',plan)
        except RuntimeError:check(True,'reject_wrong_monitor_'+key)
        else:raise RuntimeError('wrong monitor admitted')
    for j in plan['original_190_jobs']:
        if j['stage']!='core':continue
        prefix=Path(plan['baseline_support_package'])/'native/core_reproduction_v1'/j['job_id']
        r=json.loads(Path(str(prefix)+'.execution_receipt.json').read_text())
        mod.baseline_output_binding_gate(prefix,r);check(True,'actual_core_full_precision_output_bound '+j['job_id'])
        if count==0:
            for label,mutate in [('missing_full',lambda q:q['all_output_sha256'].pop(str(prefix)+'.full_precision.json')),('wrong_full_hash',lambda q:q.update(output_sha256='0'*64)),('wrong_namespace',lambda q:q['all_output_sha256'].update({'/tmp/foreign.full_precision.json':'0'*64}))]:
                bad=copy.deepcopy(r);mutate(bad)
                try:mod.baseline_output_binding_gate(prefix,bad)
                except RuntimeError:check(True,'reject_'+label)
                else:raise RuntimeError('bad output identity admitted')
        check(r['dependency_sha256_before']==r['dependency_sha256_after']==normalized,'actual_core_receipt_relocation '+j['job_id']);count+=1
    check(count==57,'all_57_actual_core_receipts_checked_without_future_completion_claim')
    c=json.loads(Path(plan['core_precision_adjudication']['comparison']).read_text())
    mod.admit_exact_core_adjudication(plan,c);check(True,'exact_preserved_core_exception_admitted_with_all_source_hashes')
    for label,mutate in [
        ('other_pair_not_admitted',lambda q:q['arithmetic_failures'][0].update(identity='insomnia__bmi')),
        ('extra_failure_not_admitted',lambda q:q['arithmetic_failures'].append(copy.deepcopy(q['arithmetic_failures'][0]))),
        ('non_P_failure_not_admitted',lambda q:q['arithmetic_failures'][0].update(rg_se=False)),
        ('different_tolerance_not_admitted',lambda q:q['arithmetic_tolerance'].update(relative=2e-12)),
        ('changed_original_source_not_admitted',lambda q:q['sources_before'][next(iter(q['sources_before']))].update(sha256='0'*64)),
    ]:
        q=copy.deepcopy(c);mutate(q)
        try:mod.admit_exact_core_adjudication(plan,q)
        except RuntimeError:check(True,label)
        else:raise RuntimeError('unadjudicated failure or changed evidence was admitted')
    names=['sensitivity_executor_v4_3.py','prepare_sensitivity_final_gates_v4_3.py','verify_sensitivity_corrections_v4_3.py']
    for name in names:ast.parse((P/'scripts'/name).read_text(),filename=name);check(True,'syntax '+name)
    # Everything below uses isolated software fixtures and mocked Popen. No
    # subprocess, stock merge, GWAS calculation or biological data is generated.
    original=dict(check_dependencies=mod.check_dependencies,check_inputs=mod.check_inputs,baseline_gate=mod.baseline_gate,compare=mod.compare_intersections,SSD=mod.SSD,lock=mod.SHARED_HEAVY_WORKER_LOCK,limits=mod.limits,sha=mod.sha,save=mod.save,
                  popen=mod.subprocess.Popen,getpgid=mod.os.getpgid,sleep=mod.time.sleep,argv=sys.argv,execute=mod.execute_under_shared_lock)
    try:
        with tempfile.TemporaryDirectory(prefix='sensitivity_v4_controls_',dir=SSD/'tmp') as directory:
            tmp=Path(directory);mod.SSD=tmp
            for name in ['logs_v4','receipts_v4','fits_v4','proofs']:(tmp/name).mkdir()
            lock=tmp/'isolated_shared.lock';mod.SHARED_HEAVY_WORKER_LOCK=lock
            invocations=[]
            def protected(_args):
                fd=os.open(str(lock),os.O_RDWR)
                try:
                    try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    except BlockingIOError:invocations.append(True);return
                    raise RuntimeError('mutex not held')
                finally:os.close(fd)
            mod.execute_under_shared_lock=protected
            sys.argv=['software_executor','--plan-sha','unused']
            try:mod.main()
            except SystemExit as e:check('Explicit --execute' in str(e),'no_execute_flag_stops')
            check(not lock.exists() and not invocations,'no_lock_or_dispatch_without_execute_flag')
            sys.argv=['software_executor','--plan-sha','unused','--execute'];mod.main()
            check(len(invocations)==1,'v3_mutex_held_through_execution_callback')
            fd=os.open(str(lock),os.O_RDWR);fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
            try:
                try:mod.main()
                except RuntimeError as e:check(str(e)=='SHARED_HEAVY_WORKER_BUSY_NO_WORKER_LAUNCHED','busy_mutex_stops_before_dispatch')
                else:raise RuntimeError('busy mutex admitted')
                check(len(invocations)==1,'occupied_mutex_created_no_worker')
                fixture=tmp/'software_fixture_plan.json';fixture.write_text('{}\n');fixture_hash=sha(fixture)
                launched=[]
                class Proc:
                    pid=123456789;returncode=0
                    def poll(self):return 0
                inject_signal=[]
                def fake_popen(*_a,**_k):
                    launched.append(True)
                    if inject_signal:mod.catchable_termination(inject_signal[0],None)
                    instance=Proc();instance.pid=123456789+len(launched);return instance
                mod.subprocess.Popen=fake_popen;mod.os.getpgid=lambda pid:pid;mod.time.sleep=lambda _seconds:None
                def held():
                    other=os.open(str(lock),os.O_RDWR)
                    try:
                        try:fcntl.flock(other,fcntl.LOCK_EX|fcntl.LOCK_NB)
                        except BlockingIOError:return True
                        return False
                    finally:os.close(other)
                def scenario(label,final_limit_fail=False,plan_hash_fail=False,save_fail=False,teardown_fail=False,all_saves_fail=False,preflight_fail=False,deferred_signal=None,pending_before_launch=False,deadline_during_hash=False):
                    mod.TERMINATION_REQUEST.clear();inject_signal[:]=[deferred_signal] if deferred_signal else []
                    if pending_before_launch:mod.catchable_termination(signal.SIGTERM,None)
                    rec=tmp/'receipts_v4'/(label+'.worker.json');prefix=tmp/'fits_v4'/label
                    queried=[];cleaned=[];expired=[]
                    def controlled_limits(*_a):
                        queried.append(True)
                        if preflight_fail or (final_limit_fail and len(queried)>1):raise OSError('CONTROLLED_RESOURCE_QUERY_ERROR')
                        return dict(software_control=True),'STAGE_DEADLINE' if expired else None
                    def controlled_sha(path):
                        if deadline_during_hash and str(path).endswith('.stdout.log'):expired.append(True)
                        if plan_hash_fail and Path(path)==fixture:raise OSError('CONTROLLED_PLAN_HASH_ERROR')
                        return original['sha'](path)
                    def controlled_save(path,data):
                        if all_saves_fail or (save_fail and Path(path)==rec):raise OSError('CONTROLLED_RECEIPT_WRITE_ERROR')
                        return original['save'](path,data)
                    class Monitor:
                        @staticmethod
                        def terminate_owned(proc):
                            check(held(),label+'_shared_lock_held_during_cleanup');cleaned.append(proc.pid)
                            if teardown_fail and len(cleaned)==1:raise OSError('CONTROLLED_TEARDOWN_ERROR')
                            return dict(initial_group_members=[],signals=[],remaining_group_members=[])
                    mod.limits=controlled_limits;mod.sha=controlled_sha;mod.save=controlled_save
                    ownership=[True,[]];before=len(launched);captured=io.StringIO()
                    with contextlib.redirect_stderr(captured):
                        try:mod.worker(['SOFTWARE_FIXTURE_NO_PROCESS'],prefix,rec,{'guard':{'poll_seconds':2}},fixture,fixture_hash,0,Monitor,ownership)
                        except (RuntimeError,OSError):pass
                        else:raise RuntimeError('controlled worker failure did not stop')
                        if not ownership[0]:
                            errors=mod.await_owned_cleanup(Monitor,ownership,fixture_hash)
                            check(bool(errors)==all_saves_fail,label+'_cleanup_receipt_failure_preserved')
                    check(ownership==[True,[]],label+'_all_owned_groups_verified_before_release')
                    no_launch=preflight_fail or pending_before_launch
                    check(len(launched)-before==(0 if no_launch else 1),label+'_mock_only_expected_launch_count')
                    if not no_launch:check(bool(cleaned),label+'_cleanup_not_skipped_by_final_metadata_failure')
                    journal=Path(str(rec)+'.failure_journal.jsonl')
                    events=[json.loads(line)['event'] for line in journal.read_text().splitlines()]
                    check(events[0]=='PRELAUNCH' and 'WORKER_TEARDOWN_AND_FINAL_QUERIES' in events,label+'_durable_failure_journal')
                    if not no_launch:check('OWNED_WORKER_LAUNCHED' in events,label+'_durable_owned_PID_witness')
                    if all_saves_fail:check('FAILURE_RECEIPT_STORAGE_FAILED' in captured.getvalue(),label+'_explicit_storage_failure_emitted')
                    else:
                        target=Path(str(rec)+'.failure.json') if save_fail else rec
                        record=json.loads(target.read_text())
                        check(record['status']=='WORKER_FAILED_PRESERVED',label+'_failure_receipt_not_success')
                        if deadline_during_hash:check(record['stop_reason']=='STAGE_DEADLINE' and len(queried)==3,label+'_post_output_and_journal_hash_gate_blocks_success')
                        if final_limit_fail:check(any('FINAL_RESOURCE_QUERY_FAILED' in x for x in record['metadata_errors']),label+'_resource_query_error_recorded')
                        if plan_hash_fail:check(any('FINAL_PLAN_HASH_FAILED' in x for x in record['metadata_errors']),label+'_hash_error_recorded')
                        if save_fail:check('primary_receipt_save_error' in record,label+'_primary_write_failure_in_fallback')
                scenario('deadline_after_output_hash',deadline_during_hash=True)
                scenario('final_resource_query',final_limit_fail=True)
                scenario('final_hash',plan_hash_fail=True)
                scenario('primary_receipt_write',save_fail=True)
                scenario('combined_failed_teardown',final_limit_fail=True,plan_hash_fail=True,save_fail=True,teardown_fail=True)
                scenario('all_receipt_storage_failure',final_limit_fail=True,teardown_fail=True,all_saves_fail=True)
                scenario('preflight_query_failure',preflight_fail=True)
                for signum in [signal.SIGINT,signal.SIGTERM,signal.SIGHUP]:scenario('deferred_signal_'+str(signum),deferred_signal=signum,teardown_fail=True)
                scenario('pending_signal_before_launch',pending_before_launch=True)
                mod.TERMINATION_REQUEST.clear()
                # Family witness expires resources while the final proof is hashed.
                # Admission/dependencies/intersections are explicitly mocked; no
                # source bodies, result-free audits or fits are executed.
                from types import SimpleNamespace
                stage_plan=tmp/'stage_fixture.json'
                fixture_plan={'guard':{'shared_heavy_worker_lock':str(lock),'deadline_seconds':129600},'input_sha256':{},'audit_jobs':[],'jobs':[]}
                stage_plan.write_text(json.dumps(fixture_plan));stage_sha=original['sha'](stage_plan)
                mod.check_dependencies=lambda _p:None;mod.check_inputs=lambda *_a:None;mod.baseline_gate=lambda _p:{}
                mod.compare_intersections=lambda *_a:{'status':'MOCK_SOFTWARE_ONLY'}
                expired=[];stage_queries=[]
                def stage_hash(path):
                    if Path(path).name=='final_stock_intersection_comparison_v4_3.json':expired.append(True)
                    return original['sha'](path)
                mod.sha=stage_hash
                def stage_limits(*_a):
                    stage_queries.append(bool(expired));return {'software_control':True},'STAGE_DEADLINE' if expired else None
                mod.limits=stage_limits;mod.save=original['save']
                mod.execute_under_shared_lock=original['execute']
                try:mod.execute_under_shared_lock(SimpleNamespace(plan=stage_plan,plan_sha=stage_sha,merge_only=True))
                except RuntimeError as e:check('FINAL_RESOURCE_GATE' in str(e),'final_family_identity_hash_deadline_stops')
                else:raise RuntimeError('family final hash deadline admitted')
                check(not (tmp/'receipts_v4/merge_only_stage_receipt_v4_3.json').exists(),'no_family_success_seal_after_final_hash_deadline')
                check((tmp/'receipts_v4/stage_failure_v4_3.json').exists(),'family_final_hash_failure_receipt_preserved')
                check(stage_queries==[False,True],'family_resource_gate_follows_final_proof_clock_hashes')
                check(held(),'shared_lock_never_released_across_failure_controls')
            finally:os.close(fd)
    finally:
        mod.check_dependencies=original['check_dependencies'];mod.check_inputs=original['check_inputs'];mod.baseline_gate=original['baseline_gate'];mod.compare_intersections=original['compare'];mod.SSD=original['SSD'];mod.SHARED_HEAVY_WORKER_LOCK=original['lock'];mod.limits=original['limits'];mod.sha=original['sha'];mod.save=original['save']
        mod.subprocess.Popen=original['popen'];mod.os.getpgid=original['getpgid'];mod.time.sleep=original['sleep'];sys.argv=original['argv'];mod.execute_under_shared_lock=original['execute']
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    check(rss<=128*(1<<20),'bounded_stdlib_peak_RSS_128MiB')
    check(sha(plan_path)==a.plan_sha,'v3_plan_unchanged_after_controls')
    check(sha(plan['superseded_operational_plan'])==plan['superseded_operational_plan_sha256'],'v1_plan_unchanged_after_controls')
    out=SSD/'proofs/sensitivity_correction_controls_v4_3.json'
    receipt=dict(schema='sensitivity_operational_correction_controls_v4',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 status='CORRECTION_BINDINGS_AND_FAILURE_CONTROLS_PASS',plan_sha256=a.plan_sha,verifier_sha256=sha(Path(__file__)),check_count=len(checks),checks=checks,
                 self_peak_RSS_bytes=rss,real_worker_subprocesses_launched=0,estimator_calls=0,stock_merge_audits_launched=0,
                 scope='Actual historical binding checks plus isolated mocked failure/cleanup software controls; no future baseline completion or scientific validity claim')
    with out.open('x') as f:json.dump(receipt,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=receipt['status'],check_count=len(checks),sha256=sha(out),self_peak_RSS_bytes=rss,real_worker_subprocesses_launched=0),indent=2))


if __name__=='__main__':main()
