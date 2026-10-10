#!/usr/bin/env python3
"""One-worker guarded execution of the separately frozen 26-command diagnostic.

Preparation never calls this executor. An explicit --execute is required. Every
baseline and result-free intersection proof must pass before the first fit.
"""
import argparse
import datetime
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
SHARED_HEAVY_WORKER_LOCK=SSD.parent/'native_heavy_worker.lock'


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path,data):
    with Path(path).open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')


def check_inputs(plan,inputs):
    for path in inputs:
        if sha(path)!=plan['input_sha256'][path]:raise RuntimeError('SOURCE_MUTATION: '+path)


def check_dependencies(plan):
    for path,expected in plan['dependencies_sha256'].items():
        if sha(path)!=expected:raise RuntimeError('DEPENDENCY_MUTATION: '+path)


def baseline_gate(plan):
    support=Path(plan['baseline_support_package'])
    receipts={}
    for stage in ['core','extension','validation']:
        monitor=P/'logs'/(stage+'_native_monitor_receipt_v4.json')
        comparison=P/'logs'/(stage+'_native_comparison_receipt_v4.json')
        if not monitor.exists() or not comparison.exists():raise RuntimeError('HISTORICAL_190_OR_NUMERICAL_REVIEW_NOT_COMPLETE: '+stage)
        m=json.loads(monitor.read_text());c=json.loads(comparison.read_text())
        if m.get('returncode')!=0 or m.get('stop_reason') is not None:raise RuntimeError('BASELINE_STAGE_UNCLEARED: '+stage)
        if not c.get('all_consumed_hashes_unchanged_after') or c.get('arithmetic_failures') or c.get('meaningful_printed_discrepancies') or c.get('fdr_boundary_status_changes'):
            raise RuntimeError('BASELINE_NUMERICAL_DISCREPANCY_REQUIRES_REVIEW: '+stage)
        if c['plan_sha256']!=plan['baseline_execution_plan_sha256'] or c['monitor_sha256']!=sha(monitor):raise RuntimeError('BASELINE_REVIEW_BINDING_DIFFERS')
        receipts[str(monitor)]=sha(monitor);receipts[str(comparison)]=sha(comparison)
    for job in plan['original_190_jobs']:
        prefix=support/'native'/(job['stage']+'_reproduction_v1')/job['job_id']
        path=Path(str(prefix)+'.execution_receipt.json')
        if not path.exists():raise RuntimeError('MISSING_BASELINE_JOB: '+job['job_id'])
        r=json.loads(path.read_text())
        if r['job']!=job or r['returncode']!=0 or not r['scientific_cardinality_gate_pass'] or not r['execution_identity_gate_pass']:
            raise RuntimeError('UNVERIFIED_BASELINE_JOB: '+job['job_id'])
        if r['input_sha256']!=r['input_sha256_after'] or r['dependency_sha256_before']!=r['dependency_sha256_after']:
            raise RuntimeError('BASELINE_IDENTITY_DIFFERS')
        expected_sources={p:plan['baseline_input_sha256'][p] for p in job['inputs']}
        expected_command=[plan['python'],'-u',str(support/'scripts/native_ldsc_capture.py'),'--ldsc-dir',plan['ldsc_dir'],
                          '--'+job['kind'],','.join(job['inputs']),'--ref-ld-chr',plan['reference_prefix'],
                          '--w-ld-chr',plan['reference_prefix'],'--print-delete-vals','--out',str(prefix)]+job['options']
        if r['input_sha256']!=expected_sources or r['dependency_sha256_before']!=plan['baseline_dependency_sha256'] or r['command']!=expected_command:
            raise RuntimeError('BASELINE_FROZEN_SOURCE_DEPENDENCY_OR_COMMAND_DIFFERS')
        if not r['all_output_sha256'] or any(sha(p)!=h for p,h in r['all_output_sha256'].items()):raise RuntimeError('BASELINE_OUTPUT_MUTATED')
        receipts[str(path)]=sha(path)
    if len(plan['original_190_jobs'])!=190:raise RuntimeError('HISTORICAL_CARDINALITY_DIFFERS')
    return receipts


def limits(plan,started,rss):
    g=plan['guard'];internal=shutil.disk_usage('/System/Volumes/Data').free;out=shutil.disk_usage(SSD).free
    total=sum(p.stat().st_size for p in SSD.rglob('*') if p.is_file())
    state=dict(recorded_utc=utc(),internal_free_bytes=internal,SSD_free_bytes=out,observed_worker_RSS_bytes=rss,
               sensitivity_namespace_bytes=total,elapsed_stage_seconds=time.monotonic()-started)
    reason=None
    if internal<g['internal_floor_bytes']:reason='INTERNAL_SPACE_GUARD'
    elif out<g['SSD_floor_bytes']:reason='SSD_SPACE_GUARD'
    elif rss>g['observed_aggregate_worker_RSS_limit_bytes']:reason='AGGREGATE_WORKER_RSS_GUARD'
    elif total>g['new_output_limit_bytes']:reason='SENSITIVITY_OUTPUT_SIZE_GUARD'
    elif state['elapsed_stage_seconds']>g['deadline_seconds']:reason='STAGE_DEADLINE'
    return state,reason


def worker(command,output_prefix,record_path,plan,plan_path,plan_hash,started,monitor,ownership):
    if record_path.exists():
        old=json.loads(record_path.read_text())
        if old.get('status')!='WORKER_COMPLETE_VERIFIED' or old.get('plan_sha256')!=plan_hash or old.get('command')!=command:
            raise RuntimeError('PRIOR_WORKER_REQUIRES_REVIEW; no overwrite')
        if any(sha(p)!=h for p,h in old['output_sha256'].items()):raise RuntimeError('CHECKPOINT_OUTPUT_CHANGED')
        print('CHECKPOINT_VERIFIED '+record_path.name,flush=True)
        return
    if list(output_prefix.parent.glob(output_prefix.name+'*')):raise RuntimeError('UNSEALED_WORKER_OUTPUT_PRESERVED')
    stdout=SSD/'logs'/(record_path.stem+'.stdout.log')
    if stdout.exists():raise RuntimeError('UNSEALED_STDOUT_PRESERVED')
    snapshot,reason=limits(plan,started,0)
    if reason:raise RuntimeError(reason)
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','TMPDIR':str(SSD/'tmp'),'XDG_CACHE_HOME':str(SSD/'cache'),
         'MPLCONFIGDIR':str(SSD/'cache/matplotlib'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1',
         'MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'}
    proc=None;peak=0;stop=None;teardown=None;samples=[];exc=None
    try:
        with stdout.open('x') as log:
            proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            if os.getpgid(proc.pid)!=proc.pid:raise RuntimeError('OWNED_SESSION_IDENTITY_FAILED')
            while proc.poll() is None:
                rss,pids=monitor.owned_rss(proc.pid);peak=max(peak,rss)
                state,stop=limits(plan,started,rss);samples.append({**state,'owned_pids':pids})
                if stop:break
                if len(samples)%30==0:print(json.dumps(dict(worker=record_path.stem,pid=proc.pid,RSS_bytes=rss,
                                                           internal_free_bytes=state['internal_free_bytes'],elapsed_stage_seconds=state['elapsed_stage_seconds'])),flush=True)
                time.sleep(plan['guard']['poll_seconds'])
    except BaseException as error:
        exc=error;stop=stop or type(error).__name__+': '+str(error)
    finally:
        if proc is not None:
            try:
                teardown=monitor.terminate_owned(proc)
                if teardown['initial_group_members'] and stop is None:stop='UNEXPECTED_DESCENDANTS_AFTER_WORKER_EXIT'
            except BaseException as error:
                stop='TEARDOWN_FAILURE: '+str(error);teardown=dict(error=str(error),verified=False)
                ownership[0]=False
                ownership[1].append(proc)
        final,final_reason=limits(plan,started,peak);stop=stop or final_reason
        identity=sha(plan_path)==plan_hash
        outputs={str(p):sha(p) for p in sorted(output_prefix.parent.glob(output_prefix.name+'*')) if p.is_file()}
        outputs[str(stdout)]=sha(stdout) if stdout.exists() else None
        record=dict(status='WORKER_COMPLETE_VERIFIED' if proc and proc.returncode==0 and stop is None and identity else 'WORKER_FAILED_PRESERVED',
                    command=command,plan_sha256=plan_hash,completed_utc=utc(),returncode=proc.returncode if proc else None,
                    owned_process_group_id=proc.pid if proc else None,
                    stop_reason=stop,plan_unchanged=identity,peak_observed_owned_RSS_bytes=peak,resource_preflight=snapshot,
                    resource_final=final,process_group_teardown=teardown,output_sha256=outputs,samples=samples)
        save(record_path,record)
    if exc:raise exc
    if record['status']!='WORKER_COMPLETE_VERIFIED':raise RuntimeError('WORKER_STOPPED_PRESERVED: '+str(stop))


def compare_intersections(plan,plan_hash):
    pairs=[];files={}
    for job in plan['jobs']:
        arms={}
        for arm in ['baseline','sensitivity']:
            audit=next(a for a in plan['audit_jobs'] if a['arm']==arm and a['sensitivity_job_id']==job['job_id'])
            path=Path(audit['out_prefix']+'.intersection.json');r=json.loads(path.read_text());files[str(path)]=sha(path)
            if r['plan_sha256']!=plan_hash or r['estimator_calls']!=0 or len(r['final_intersections'])!=job['estimates']:
                raise RuntimeError('STOCK_AUDIT_BINDING_OR_COUNT_DIFFERS')
            arms[arm]=r['final_intersections']
        for a,b in zip(arms['baseline'],arms['sensitivity']):
            fields=['final_ordered_SNP_count','final_ordered_SNP_sha256']
            if job['kind']=='rg':fields+=['final_ordered_SNP_allele_compatibility_sha256','final_ordered_aligned_Z1_Z2_float64_big_endian_sha256']
            else:fields+=['final_ordered_Z_float64_big_endian_sha256']
            if any(a[k]!=b[k] for k in fields):raise RuntimeError('FINAL_INTERSECTION_OR_ALLELE_OR_Z_DIFFERS: '+job['job_id'])
            nkey='final_ordered_N1_N2_float64_big_endian_sha256' if job['kind']=='rg' else 'final_ordered_N_float64_big_endian_sha256'
            sameN=a[nkey]==b[nkey]
            expected_same=job['job_id'].startswith('lipid_two_step_')
            if sameN!=expected_same:raise RuntimeError('N_CONVENTION_IDENTITY_UNEXPECTED')
            if not expected_same and a.get('two_step_hsq_mask')!=b.get('two_step_hsq_mask'):raise RuntimeError('BINARY_H2_MASK_DIFFERS')
            pairs.append(dict(job_id=job['job_id'],outcome=Path(b.get('p2',b.get('input'))).name,
                              final_ordered_identity_fields_equal=fields,N_identity_same=sameN,
                              baseline=a,sensitivity=b))
    if len(pairs)!=62:raise RuntimeError('PAIRED_INTERSECTION_CARDINALITY_DIFFERS')
    return dict(schema='frozen_62_paired_stock_intersection_proof_v1',completed_utc=utc(),plan_sha256=plan_hash,
                paired_estimates=62,identity_instances=124,proof_sha256=files,pairs=pairs,
                estimator_calls=0,baseline_scope='stock merge replay, not instrumentation of completed original fit',status='ALL_FINAL_ORDERED_SNP_ALLELE_Z_IDENTITIES_MATCH')


def native_result_gate(job,plan,plan_hash):
    path=Path(job['out_prefix']+'.full_precision.json');r=json.loads(path.read_text())
    audit=next(a for a in plan['audit_jobs'] if a['arm']=='sensitivity' and a['sensitivity_job_id']==job['job_id'])
    expected=json.loads(Path(audit['out_prefix']+'.intersection.json').read_text())['final_intersections']
    if r['plan_sha256']!=plan_hash or len(r['estimates'])!=job['estimates'] or r['final_intersections']!=expected:
        raise RuntimeError('NATIVE_CAPTURE_COUNT_OR_INTERSECTION_DIFFERS')
    arithmetic=[]
    def scalar_deletes(values):
        flat=[]
        for value in values:
            if isinstance(value,list):
                if len(value)!=1:raise RuntimeError('UNEXPECTED_DELETE_ARRAY_SHAPE')
                value=value[0]
            if not isinstance(value,(int,float)) or not math.isfinite(value):raise RuntimeError('NONFINITE_DELETE_VALUE')
            flat.append(value)
        if len(flat)!=200:raise RuntimeError('UNEXPECTED_DELETE_ARRAY_COUNT')
        return flat
    def centered_se(values):
        center=math.fsum(values)/len(values)
        return math.sqrt((len(values)-1)/len(values)*math.fsum((x-center)**2 for x in values))
    def close(actual,expected,p=False):
        # Same independent within-fit arithmetic tolerance as the historical v4 collator;
        # never applied to differences between estimators or N conventions.
        return math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-300 if p else 1e-15)
    for row in r['estimates']:
        if job['kind']=='rg':
            if row.get('status')!='NATIVE_ESTIMATE_RETURNED':raise RuntimeError('PAIR_ESTIMATION_FAILURE_PRESERVED')
            objects=[row['hsq1'],row['hsq2'],row['gencov']]
            fields=[row.get(x) for x in ['rg_ratio','rg_se','z','p']]
        else:objects=[row];fields=[]
        fields += [o.get(k) for o in objects for k in ['tot','tot_se','intercept','intercept_se']]
        if any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in fields):raise RuntimeError('EMPTY_OR_NONFINITE_ESTIMATE_PRESERVED')
        if any(o['tot_se']<=0 for o in objects) or (job['kind']=='rg' and row['rg_se']<=0):raise RuntimeError('INVALID_SE_PRESERVED')
        deletes=[];checks={}
        for i,o in enumerate(objects):
            values=scalar_deletes(o['tot_delete_values']);deletes.append(values)
            checks['object_'+str(i)+'_tot_se']=close(centered_se(values),o['tot_se'])
            intercept_values=scalar_deletes(o['intercept_delete_values'])
            checks['object_'+str(i)+'_intercept_se']=close(centered_se(intercept_values),o['intercept_se'])
        if job['kind']=='rg':
            a,b,c=objects
            ratio=c['tot']/math.sqrt(a['tot']*b['tot'])
            checks['rg_ratio']=close(ratio,row['rg_ratio'])
            if any(x<=0 or y<=0 for x,y in zip(deletes[0],deletes[1])):raise RuntimeError('NONPOSITIVE_DELETE_RATIO_DENOMINATOR_REQUIRES_REVIEW')
            ratios=[z/math.sqrt(x*y) for x,y,z in zip(*deletes)]
            se=centered_se(ratios)
            jack=math.fsum(200*ratio-199*x for x in ratios)/200
            checks['rg_se']=close(se,row['rg_se']);checks['rg_jknife']=close(jack,row['rg_jknife'])
            checks['z']=close(ratio/se,row['z']);checks['p']=close(math.erfc(abs(ratio/se)/math.sqrt(2)),row['p'],True)
        if not all(checks.values()):raise RuntimeError('WITHIN_FIT_NUMERICAL_ARITHMETIC_CONTROL_FAILED: '+str(checks))
        arithmetic.append(dict(outcome=row.get('p2',row.get('input')),checks=checks,
                               cross_fit_boundary_alignment_established=False))
    validation=SSD/'receipts'/(job['job_id']+'.numerical_validation.json')
    record=dict(plan_sha256=plan_hash,full_precision_sha256=sha(path),final_intersection_matches_stock_audit=True,
                status='WITHIN_FIT_NUMERICAL_AND_IDENTITY_CONTROLS_PASS',checks=arithmetic,
                control_tolerance_scope='historical within-fit arithmetic only; no tolerance for actual sensitivity differences')
    if validation.exists():
        if json.loads(validation.read_text())!=record:raise RuntimeError('PRIOR_NUMERICAL_VALIDATION_DIFFERS')
    else:save(validation,record)
    return sha(path)


def await_owned_cleanup(monitor,ownership,plan_hash):
    """Keep the shared flock held while any owned worker teardown is unresolved.

    Catch interruption during cleanup too: an interrupted parent may not release
    the shared heavy-worker mutex while a launched descendant remains active.
    An external SIGKILL cannot be handled by Python and needs host-level review.
    """
    while ownership[1]:
        proc=ownership[1][0]
        try:
            teardown=monitor.terminate_owned(proc)
        except BaseException as error:
            print('SHARED_HEAVY_WORKER_LOCK_HELD: owned group '+str(proc.pid)+' remains unresolved: '+str(error),file=sys.stderr,flush=True)
            try:time.sleep(30)
            except BaseException:pass
            continue
        ownership[1].pop(0)
        save(SSD/'receipts'/('owned_group_'+str(proc.pid)+'.recovery.json'),
             dict(completed_utc=utc(),plan_sha256=plan_hash,owned_process_group_id=proc.pid,
                  status='OWNED_GROUP_CLEANUP_VERIFIED_BEFORE_SHARED_LOCK_RELEASE',process_group_teardown=teardown))
    ownership[0]=True


def execute_under_shared_lock(args):
    if sha(args.plan)!=args.plan_sha:raise RuntimeError('PLAN_SHA_DIFFERS')
    plan=json.loads(args.plan.read_text());check_dependencies(plan);check_inputs(plan,plan['input_sha256'])
    if plan['guard']['shared_heavy_worker_lock']!=str(SHARED_HEAVY_WORKER_LOCK):raise RuntimeError('SHARED_LOCK_PATH_DIFFERS')
    baseline=baseline_gate(plan)
    monitor_path=P/'scripts/30_prepare_and_run_ssd_native_campaign.py'
    spec=importlib.util.spec_from_file_location('_frozen_owned_process_monitor',monitor_path);monitor=importlib.util.module_from_spec(spec);spec.loader.exec_module(monitor)
    # An exclusive lock prevents two sensitivity queues from running concurrently.
    lock=SSD/'sensitivity_executor.lock';fd=os.open(str(lock),os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    os.write(fd,(str(os.getpid())+'\n').encode());os.close(fd)
    ownership=[True,[]];completed=[]
    try:
        clock_path=SSD/'receipts/stage_clock_v1.json'
        if clock_path.exists():
            clock=json.loads(clock_path.read_text())
            if clock['plan_sha256']!=args.plan_sha:raise RuntimeError('STAGE_CLOCK_BINDING_DIFFERS')
        else:
            clock=dict(plan_sha256=args.plan_sha,initial_started_utc=utc(),initial_started_epoch=time.time(),deadline_seconds=plan['guard']['deadline_seconds'])
            save(clock_path,clock)
        elapsed=time.time()-clock['initial_started_epoch']
        if elapsed<0:raise RuntimeError('SYSTEM_CLOCK_PRECEDES_PRESERVED_STAGE_START')
        started=time.monotonic()-elapsed
        for audit in plan['audit_jobs']:
            check_dependencies(plan);check_inputs(plan,audit['inputs'])
            cmd=[plan['python'],'-B',str(P/'scripts/37_stock_ldsc_intersection_audit.py'),'--plan',str(args.plan),'--plan-sha',args.plan_sha,'--audit-id',audit['audit_id']]
            worker(cmd,Path(audit['out_prefix']),SSD/'receipts'/(audit['audit_id']+'.worker.json'),plan,args.plan,args.plan_sha,started,monitor,ownership)
            check_dependencies(plan);check_inputs(plan,audit['inputs'])
        proof_path=SSD/'proofs/final_stock_intersection_comparison_v1.json'
        proof=compare_intersections(plan,args.plan_sha)
        if proof_path.exists():
            old=json.loads(proof_path.read_text())
            if any(old[k]!=proof[k] for k in ['plan_sha256','paired_estimates','identity_instances','proof_sha256','pairs','status']):raise RuntimeError('PRIOR_FINAL_INTERSECTION_PROOF_DIFFERS')
        else:save(proof_path,proof)
        if not args.merge_only:
            for job in plan['jobs']:
                check_dependencies(plan);check_inputs(plan,job['inputs'])
                cmd=[plan['python'],'-u','-B',str(P/'scripts/38_sensitivity_ldsc_capture.py'),'--plan',str(args.plan),'--plan-sha',args.plan_sha,'--job-id',job['job_id']]+job['ldsc_args']
                worker(cmd,Path(job['out_prefix']),SSD/'receipts'/(job['job_id']+'.worker.json'),plan,args.plan,args.plan_sha,started,monitor,ownership)
                digest=native_result_gate(job,plan,args.plan_sha)
                check_dependencies(plan);check_inputs(plan,job['inputs'])
                completed.append(dict(job_id=job['job_id'],estimates=job['estimates'],full_precision_sha256=digest))
                print(json.dumps(dict(native_job_complete=job['job_id'],completed_commands=len(completed))),flush=True)
        check_dependencies(plan);check_inputs(plan,plan['input_sha256'])
        name='merge_only_stage_receipt_v1.json' if args.merge_only else 'sensitivity_execution_receipt_v1.json'
        save(SSD/'receipts'/name,dict(completed_utc=utc(),plan_sha256=args.plan_sha,baseline_gate_receipt_sha256=baseline,
                                    final_intersection_proof_sha256=sha(proof_path),completed_jobs=completed,
                                    stage_clock_sha256=sha(clock_path),
                                    new_fit_count=sum(x['estimates'] for x in completed),merge_only=args.merge_only,
                                    historical_protocol_and_threshold_changes=0,elapsed_seconds=time.monotonic()-started))
    finally:
        if not ownership[0]:await_owned_cleanup(monitor,ownership,args.plan_sha)
        lock.unlink()


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,default=SSD/'sensitivity_operational_plan_v1.json');p.add_argument('--plan-sha',required=True);p.add_argument('--execute',action='store_true');p.add_argument('--merge-only',action='store_true');args=p.parse_args()
    if not args.execute:raise SystemExit('Prepared only. Explicit --execute is required; no worker launched.')
    # Shared with all native and raw-pipeline heavy workers; keep held for the
    # entire execution, including stock merge audits and owned-group teardown.
    # Opening does not truncate or change another worker's lock metadata.
    fd=os.open(str(SHARED_HEAVY_WORKER_LOCK),os.O_RDWR|os.O_CREAT,0o600)
    try:
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('SHARED_HEAVY_WORKER_BUSY_NO_WORKER_LAUNCHED')
        execute_under_shared_lock(args)
    finally:
        os.close(fd)


if __name__=='__main__':main()
