#!/usr/bin/env python3
"""Bounded stdlib plan/software controls; never load LDSC, merge GWAS or fit.

Software fixtures below test file-lock and fail-closed execution behavior only.
They are not GWAS data, scientific calibration or biological findings.
"""
import argparse
import ast
import contextlib
import datetime
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import resource
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
MEMBER_SHA='3ce7b50bb4e1289ef5a69d8e46dcef759963703ef7b317e156bdee0d4c40e5e1'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,default=SSD/'sensitivity_operational_plan_v1.json');p.add_argument('--plan-sha',required=True);a=p.parse_args()
    checks=[]
    def check(condition,name):
        if not condition:raise RuntimeError('CONTROL_FAILED: '+name)
        checks.append(name)
    check(sha(a.plan)==a.plan_sha,'operational_plan_exact_SHA')
    plan=json.loads(a.plan.read_text())
    member=Path(plan['scientific_member_manifest']);m=json.loads(member.read_text())
    check(sha(member)==MEMBER_SHA==plan['scientific_member_manifest_sha256'],'admitted_scientific_manifest_exact_SHA')
    for path,digest in plan['dependencies_sha256'].items():check(sha(path)==digest,'fresh_dependency_SHA: '+path)
    for path,digest in plan['input_sha256'].items():check(sha(path)==digest,'fresh_target_input_SHA: '+path)
    check(len(plan['original_190_jobs'])==190,'original_campaign_190_retained')
    check(len(plan['jobs'])==26 and sum(j['estimates'] for j in plan['jobs'])==62,'26_commands_62_estimates')
    check(len(plan['audit_jobs'])==52 and sum(j['expected_identities'] for j in plan['audit_jobs'])==124,'52_merge_commands_124_identity_instances')
    check(plan['estimator_calls_in_preparation']==0 and plan['pinned_code_edits']==0,'preparation_only_no_estimator_or_pinned_edit')
    check(plan['no_biological_testing_family_or_QC_threshold_changes'] is True,'unchanged_biological_families_and_QC')
    check(plan['guard']['internal_floor_bytes']==3*(1<<30) and plan['guard']['SSD_floor_bytes']==5*(1<<30),'unchanged_3GiB_5GiB_floors')
    check(plan['guard']['worker_count']==1 and plan['guard']['BLAS_threads']==1,'one_worker_one_BLAS_thread')
    check(plan['guard']['observed_aggregate_worker_RSS_limit_bytes']==2*(1<<30),'aggregate_RSS_2GiB')
    check(plan['guard']['new_output_limit_bytes']==4*(1<<30) and plan['guard']['deadline_seconds']==36*3600,'4GiB_outputs_and_persistent_36h_deadline')
    check(plan['guard']['shared_heavy_worker_lock']==str(SSD.parent/'native_heavy_worker.lock'),'exact_shared_heavy_worker_mutex_path')
    for j,expected in zip(plan['jobs'],m['jobs']):
        for key,value in expected.items():check(j[key]==value,'scientific_job_unchanged '+j['job_id']+' '+key)
        check(j['out_prefix'].startswith(str(SSD/'fits')+'/'),'new_fit_namespace '+j['job_id'])
        cli=j['ldsc_args']
        check(cli[:2]==['--'+j['kind'],','.join(j['inputs'])],'exact_ordered_input_argv '+j['job_id'])
        check(cli[2:10]==['--ref-ld-chr',plan['reference_prefix'],'--w-ld-chr',plan['reference_prefix'],'--n-blocks','200','--print-delete-vals','--out'], 'fixed_reference_blocks_delete_argv '+j['job_id'])
        check(cli[10:]==[j['out_prefix']]+j['options'],'only_declared_option_changes '+j['job_id'])
        if j['job_id'].startswith('lipid_two_step_'):check(j['kind']=='rg' and j['estimates']==3 and j['options']==['--two-step','30'],'explicit_lipid_two_step30 '+j['job_id'])
        elif j['kind']=='rg':check(j['estimates']==2 and j['options']==[],'binary_default_estimator '+j['job_id'])
        else:
            basename=Path(j['inputs'][0]).name
            case,total,K=(2182,376169,.002) if basename=='ms.sumstats.gz' else (2993,290130,.02)
            check(j['options']==['--samp-prev',str(case/total),'--pop-prev',str(K)],'exact_configured_case_fraction_K '+j['job_id'])
    for j in plan['audit_jobs']:
        check(j['out_prefix'].startswith(str(SSD/'proofs/intersections')+'/'),'new_audit_namespace '+j['audit_id'])
        check(j['historical_fit_instrumented'] is False,'baseline_replay_scope '+j['audit_id'])
    proof_path=Path(plan['derivative_proof']);proof=json.loads(proof_path.read_text())
    check(sha(proof_path)==plan['derivative_proof_sha256'],'sealed_full_stream_proof_SHA')
    check(sha(P/'scripts/36_materialize_total_n.py')==proof['worker_sha256'],'sealed_materializer_worker_binding')
    for r in proof['results']:
        check(r['source_sha256_before']==r['source_sha256_after']==sha(r['source']),'original_stream_before_after_SHA '+Path(r['source']).name)
        check(r['counts']['template_rows']==r['full_stream_verified_rows']==1217311,'entire_template_compared '+Path(r['source']).name)
        check(r['original_and_derivative_invariant_match'] and r['literal_missingness_preserved'] and r['gzip']['independent_second_compression_matches'],'all_non_N_bytes_missingness_order_gzip '+Path(r['source']).name)
    names=['36_materialize_total_n.py','37_stock_ldsc_intersection_audit.py','38_sensitivity_ldsc_capture.py','39_prepare_sensitivity_execution.py','40_run_estimator_sensitivities.py','sensitivity_capture_common.py','verify_sensitivity_operational_plan.py']
    for name in names:ast.parse((P/'scripts'/name).read_text(),filename=name);check(True,'syntax_parse '+name)
    # Import only our stdlib executor; execute no worker and do not import LDSC.
    executor_path=P/'scripts/40_run_estimator_sensitivities.py'
    spec=importlib.util.spec_from_file_location('_software_control_executor',executor_path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    with tempfile.TemporaryDirectory(prefix='sensitivity_controls_',dir=SSD/'tmp') as directory:
        tmp=Path(directory);lock=tmp/'isolated_test.lock';mod.SHARED_HEAVY_WORKER_LOCK=lock
        original_argv=sys.argv;invocations=[]
        def protected(_args):
            fd=os.open(str(lock),os.O_RDWR)
            try:
                try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:invocations.append('LOCK_HELD_INSIDE_EXECUTION');return
                raise RuntimeError('second heavy worker could enter')
            finally:os.close(fd)
        try:
            mod.execute_under_shared_lock=protected
            sys.argv=['executor','--plan-sha','unused']
            try:mod.main()
            except SystemExit as e:check('Explicit --execute' in str(e),'explicit_execute_required')
            else:raise RuntimeError('missing execute flag did not stop')
            check(not lock.exists() and not invocations,'no_execution_or_lock_without_explicit_flag')
            sys.argv=['executor','--plan-sha','unused','--execute'];mod.main()
            check(invocations==['LOCK_HELD_INSIDE_EXECUTION'],'shared_mutex_held_across_execution_callback')
            first=os.open(str(lock),os.O_RDWR);fcntl.flock(first,fcntl.LOCK_EX|fcntl.LOCK_NB)
            try:
                try:mod.main()
                except RuntimeError as e:check(str(e)=='SHARED_HEAVY_WORKER_BUSY_NO_WORKER_LAUNCHED','busy_mutex_fails_closed')
                else:raise RuntimeError('occupied heavy worker mutex was bypassed')
                check(len(invocations)==1,'no_callback_or_worker_when_shared_mutex_busy')
            finally:os.close(first)
            def fail(_args):raise RuntimeError('CONTROLLED_PREWORKER_FAILURE')
            mod.execute_under_shared_lock=fail
            try:mod.main()
            except RuntimeError as e:check(str(e)=='CONTROLLED_PREWORKER_FAILURE','preworker_failure_preserved')
            fd=os.open(str(lock),os.O_RDWR);fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);os.close(fd)
            check(True,'shared_mutex_released_after_no_worker_failure')
            # A mocked teardown failure must retain the shared lock through retry;
            # no process is created and no estimator is imported by this control.
            fd=os.open(str(lock),os.O_RDWR);fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
            retry=[]
            class Proc:pid=123456789
            class Monitor:
                @staticmethod
                def terminate_owned(proc):
                    other=os.open(str(lock),os.O_RDWR)
                    try:
                        try:fcntl.flock(other,fcntl.LOCK_EX|fcntl.LOCK_NB)
                        except BlockingIOError:pass
                        else:raise RuntimeError('mutex released before owned cleanup')
                    finally:os.close(other)
                    retry.append(proc.pid)
                    if len(retry)==1:raise RuntimeError('CONTROLLED_TEARDOWN_FAILURE')
                    return dict(initial_group_members=[],signals=[],remaining_group_members=[])
            mod.SSD=tmp;(tmp/'receipts').mkdir();mod.time.sleep=lambda _seconds:None
            ownership=[False,[Proc()]]
            try:
                with contextlib.redirect_stderr(io.StringIO()):mod.await_owned_cleanup(Monitor,ownership,'isolated_software_control')
                check(ownership==[True,[]] and len(retry)==2,'shared_mutex_retained_through_failed_owned_cleanup_retry')
            finally:os.close(fd)
        finally:sys.argv=original_argv
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    check(rss<=128*(1<<20),'bounded_stdlib_self_peak_RSS_128MiB')
    check(sha(a.plan)==a.plan_sha,'plan_unchanged_after_controls')
    receipt=dict(schema='result_free_sensitivity_operational_controls_v1',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 plan_sha256=a.plan_sha,verifier_sha256=sha(Path(__file__)),status='PLAN_AND_SOFTWARE_CONTROLS_PASS',
                 check_count=len(checks),checks=checks,self_peak_RSS_bytes=rss,worker_subprocesses_launched=0,LDSC_imported=False,
                 stock_merge_audits_launched=0,estimator_calls=0,scope='structural plan binding, fresh hashes and mocked mutex controls; not native or scientific validation')
    out=SSD/'proofs/sensitivity_operational_controls_v1.json'
    with out.open('x') as f:json.dump(receipt,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=receipt['status'],check_count=len(checks),receipt=str(out),sha256=sha(out),self_peak_RSS_bytes=rss,estimator_calls=0),indent=2))


if __name__=='__main__':main()
