#!/usr/bin/env python3
"""Additive v4_3 executor: exact baseline relocation/adjudication, robust cleanup.

Preparation never calls this executor. An explicit --execute is required. Every
baseline and result-free intersection proof must pass before the first fit.
"""
import argparse
import builtins
import datetime
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import native_stage_completion_v4_3 as stage_completion
from terminal_commit_common_v2 import TerminalCommit

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
SHARED_HEAVY_WORKER_LOCK=SSD.parent/'native_heavy_worker.lock'
TERMINATION_REQUEST=[]


def catchable_termination(signum,frame):
    # Defer requests through child registration, teardown and receipt writes.
    # The handler never raises asynchronously or performs diagnostic I/O.
    TERMINATION_REQUEST.append(signum)


def check_termination(where):
    if TERMINATION_REQUEST:raise RuntimeError('DEFERRED_TERMINATION_'+where+': '+str(TERMINATION_REQUEST))


def safe_print(*args,**kwargs):
    try:builtins.print(*args,**kwargs)
    except BaseException:pass


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path,data):
    # Serialize before file creation and exclusively create without overwriting.
    # The SSD does not support hard links. A failed partial receipt is preserved
    # and cannot pass JSON/checkpoint parsing; a separate failure receipt follows.
    payload=json.dumps(data,indent=2,allow_nan=False)+'\n'
    with Path(path).open('x') as f:f.write(payload);f.flush();os.fsync(f.fileno())


def admit_exact_core_adjudication(plan,comparison):
    """Admit only the unchanged, individually reviewed snoring/BMI P failure."""
    binding=plan['core_precision_adjudication']
    for path,digest in binding['evidence_sha256'].items():
        if sha(path)!=digest:raise RuntimeError('CORE_ADJUDICATION_EVIDENCE_CHANGED: '+path)
    if comparison['arithmetic_failures']!=[binding['exact_preserved_failure']]:raise RuntimeError('CORE_FAILURE_IS_NOT_EXACT_ADJUDICATED_RECORD')
    a=json.loads(Path(binding['adjudication']).read_text())
    whole=json.loads(Path(binding['whole_core']).read_text())
    command=json.loads(Path(binding['command_binding']).read_text())
    tolerance={'relative':1e-12,'absolute_estimates':1e-15,'absolute_p':1e-300}
    fixed={'rtol':1e-12,'atol_statistics':1e-15,'atol_p':1e-300}
    if comparison['arithmetic_tolerance']!=tolerance or whole['fixed_tolerances']!=fixed or not a['fixed_tolerance_unchanged']:
        raise RuntimeError('CORE_ADJUDICATION_TOLERANCE_CHANGED')
    if a['whole_core_receipt_sha256']!=sha(binding['whole_core']) or a['report_sha256']!=sha(binding['report']) or a['command_binding_receipt_sha256']!=sha(binding['command_binding']):
        raise RuntimeError('CORE_ADJUDICATION_CROSS_BINDING_DIFFERS')
    if a['verdict']!='PROCESSED_INPUT_CORE_NUMERICAL_PASS_WITH_PRECISION_ADJUDICATION' or not a['initial_collator_failure_retained'] or a['native_meaningful_numerical_discrepancy_identified']:
        raise RuntimeError('CORE_ADJUDICATION_NOT_UNCHANGED_PRECISION_PASS')
    if (a['job_count'],a['h2_count'],a['rg_count'],a['primary_fdr_count'],a['sensitivity_fdr_count'],a['fdr_boundary_changes'])!=(57,45,396,153,8,0):
        raise RuntimeError('CORE_ADJUDICATION_CARDINALITY_DIFFERS')
    if not whole['all_independent_checks_pass'] or not command['all_pass'] or not command['consumed_hashes_unchanged']:
        raise RuntimeError('CORE_INDEPENDENT_NUMERICAL_OR_COMMAND_PROOF_UNCLEARED')
    if (whole['job_count'],whole['h2_count'],whole['rg_count'],whole['native_primary_bh_positives'],whole['native_sensitivity_bh_positives'],whole['bh_boundary_changes'])!=(57,45,396,153,8,0):
        raise RuntimeError('WHOLE_CORE_PROOF_CARDINALITY_DIFFERS')
    detail=a['snoring_bmi']
    if detail!=whole['snoring_bmi_adjudication'] or detail['pair_id']!='snoring__bmi' or detail['pseudo_p_fixed_tolerance_pass'] or not detail['centered_p_fixed_tolerance_pass'] or not detail['captured_z_erfc_fixed_tolerance_pass'] or detail['tolerance_relaxed']:
        raise RuntimeError('EXACT_SNORING_BMI_PRECISION_ADJUDICATION_DIFFERS')
    if len(whole['pairs'])!=396 or len(whole['h2'])!=45 or len(whole['jobs'])!=57 or any(not r['all_pass'] for r in whole['pairs']+whole['h2']):
        raise RuntimeError('WHOLE_CORE_ROW_PROOFS_UNCLEARED')
    # These maps cover original logs/tables, native outputs, commands and fixed
    # dependencies. They do not pretend that the reviewer reread raw GWAS bodies.
    for recorded in [comparison['sources_before'],whole['inputs'],command['inputs']]:
        for path,meta in recorded.items():
            expected=meta['sha256'] if isinstance(meta,dict) else meta
            if sha(path)!=expected:raise RuntimeError('CORE_ADJUDICATED_SOURCE_CHANGED: '+path)
    return binding['evidence_sha256']


def check_inputs(plan,inputs):
    for path in inputs:
        if sha(path)!=plan['input_sha256'][path]:raise RuntimeError('SOURCE_MUTATION: '+path)


def check_dependencies(plan):
    for path,expected in plan['dependencies_sha256'].items():
        if sha(path)!=expected:raise RuntimeError('DEPENDENCY_MUTATION: '+path)


def relocated_baseline_dependencies(plan):
    original=str(ROOT/'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    relocated=str(Path(plan['baseline_support_package'])/'scripts/native_ldsc_capture.py')
    expected=dict(plan['baseline_dependency_sha256'])
    pinned='bd1b87e403fca2bad4c30b6a662a64f7c3d1a47f6aa9884f5157e370cb95707d'
    if expected.get(original)!=pinned or relocated in expected or sha(original)!=pinned or sha(relocated)!=pinned:
        raise RuntimeError('EXACT_BASELINE_CAPTURE_RELOCATION_DIFFERS')
    expected[relocated]=expected.pop(original)
    if len(expected)!=len(plan['baseline_dependency_sha256']) or expected!=plan['baseline_relocated_dependency_sha256']:
        raise RuntimeError('BASELINE_RELOCATION_CHANGED_OTHER_DEPENDENCIES')
    return expected


def baseline_monitor_binding_gate(m,stage,plan):
    if m.get('stage')!=stage or m.get('plan_sha256')!=plan['baseline_execution_plan_sha256']:
        raise RuntimeError('BASELINE_MONITOR_STAGE_OR_PLAN_DIFFERS: '+stage)


def baseline_output_binding_gate(prefix,r):
    outputs=r.get('all_output_sha256',{})
    full=str(prefix)+'.full_precision.json'
    if full not in outputs or outputs[full]!=r.get('output_sha256'):
        raise RuntimeError('BASELINE_EXACT_FULL_PRECISION_OUTPUT_BINDING_REQUIRED')
    allowed={str(prefix)+suffix for suffix in ['.full_precision.json','.log','.stdout.log']}
    job=r['job']
    if job['kind']=='h2':allowed.update(str(prefix)+suffix for suffix in ['.delete','.part_delete'])
    elif job['kind']=='rg':
        for other in job['inputs'][1:]:
            stem=str(prefix)+Path(job['inputs'][0]).name+'_'+Path(other).name
            allowed.update(stem+suffix for suffix in ['.hsq1.delete','.hsq2.delete','.gencov.delete'])
    else:raise RuntimeError('UNKNOWN_BASELINE_JOB_KIND')
    for path,digest in outputs.items():
        candidate=Path(path)
        if str(candidate)!=path or candidate.parent!=prefix.parent or path not in allowed:
            raise RuntimeError('BASELINE_OUTPUT_OUTSIDE_DECLARED_JOB_STAGE_NAMESPACE: '+path)
        if not candidate.is_file() or candidate.is_symlink() or sha(candidate)!=digest:
            raise RuntimeError('BASELINE_OUTPUT_MUTATED: '+path)


def baseline_gate(plan):
    support=Path(plan['baseline_support_package'])
    relocated=relocated_baseline_dependencies(plan)
    receipts={}
    for stage in ['core','extension','validation']:
        monitor=stage_completion.stage_monitor_path(stage)
        comparison=P/'logs'/(stage+'_native_comparison_receipt_v4.json')
        if not monitor.exists() or not comparison.exists():raise RuntimeError('HISTORICAL_190_OR_NUMERICAL_REVIEW_NOT_COMPLETE: '+stage)
        m=stage_completion.stage_monitor(stage);c=json.loads(comparison.read_text())
        baseline_monitor_binding_gate(m,stage,plan)
        if m.get('returncode')!=0 or m.get('stop_reason') is not None or m['process_group_teardown']['remaining_group_members']:
            raise RuntimeError('BASELINE_STAGE_UNCLEARED: '+stage)
        if not c.get('all_consumed_hashes_unchanged_after') or c.get('meaningful_printed_discrepancies') or c.get('fdr_boundary_status_changes'):
            raise RuntimeError('BASELINE_NUMERICAL_DISCREPANCY_REQUIRES_REVIEW: '+stage)
        if c.get('arithmetic_failures'):
            if stage!='core':raise RuntimeError('BASELINE_NUMERICAL_DISCREPANCY_REQUIRES_REVIEW: '+stage)
            receipts.update(admit_exact_core_adjudication(plan,c))
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
        if r['input_sha256']!=expected_sources or r['dependency_sha256_before']!=relocated or r['command']!=expected_command:
            raise RuntimeError('BASELINE_FROZEN_SOURCE_DEPENDENCY_OR_COMMAND_DIFFERS')
        baseline_output_binding_gate(prefix,r)
        receipts[str(path)]=sha(path)
    if len(plan['original_190_jobs'])!=190:raise RuntimeError('HISTORICAL_CARDINALITY_DIFFERS')
    return receipts


def limits(plan,started,rss):
    g=plan['guard'];internal=shutil.disk_usage('/System/Volumes/Data').free;out=shutil.disk_usage(SSD).free
    total=sum(p.stat().st_size for p in SSD.rglob('*') if p.is_file())
    state=dict(recorded_utc=utc(),internal_free_bytes=internal,SSD_free_bytes=out,observed_worker_RSS_bytes=rss,
               sensitivity_namespace_bytes=total,elapsed_stage_seconds=time.monotonic()-started)
    global_bytes=sum(p.stat().st_size for p in SSD.parent.rglob('*') if p.is_file() and not p.is_symlink())
    state['new_campaign_namespace_bytes']=global_bytes
    reason=None
    if internal<g['internal_floor_bytes']:reason='INTERNAL_SPACE_GUARD'
    elif out<g['SSD_floor_bytes']:reason='SSD_SPACE_GUARD'
    elif rss>g['observed_aggregate_worker_RSS_limit_bytes']:reason='AGGREGATE_WORKER_RSS_GUARD'
    elif total>g['new_output_limit_bytes']:reason='SENSITIVITY_OUTPUT_SIZE_GUARD'
    elif global_bytes>g['global_reservation_bytes']:reason='300GIB_CAMPAIGN_GUARD'
    elif state['elapsed_stage_seconds']>g['deadline_seconds']:reason='STAGE_DEADLINE'
    return state,reason


def stage_resource_gate(plan,started,where):
    check_termination(where)
    state,reason=limits(plan,started,0)
    if reason:raise RuntimeError('FINAL_RESOURCE_GATE_'+where+': '+reason)
    check_termination(where+'_AFTER_QUERY')
    return state


def worker(command,output_prefix,record_path,plan,plan_path,plan_hash,started,monitor,ownership):
    if record_path.exists():
        old=json.loads(record_path.read_text())
        if old.get('status')!='WORKER_COMPLETE_VERIFIED' or old.get('plan_sha256')!=plan_hash or old.get('command')!=command:
            raise RuntimeError('PRIOR_WORKER_REQUIRES_REVIEW; no overwrite')
        if any(sha(p)!=h for p,h in old['output_sha256'].items()):raise RuntimeError('CHECKPOINT_OUTPUT_CHANGED')
        stage_resource_gate(plan,started,'AFTER_CHECKPOINT_OUTPUT_HASHES')
        safe_print('CHECKPOINT_VERIFIED '+record_path.name,flush=True)
        return
    if list(output_prefix.parent.glob(output_prefix.name+'*')):raise RuntimeError('UNSEALED_WORKER_OUTPUT_PRESERVED')
    stdout=SSD/'logs_v4'/(record_path.stem+'.stdout.log')
    if stdout.exists():raise RuntimeError('UNSEALED_STDOUT_PRESERVED')
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','TMPDIR':str(SSD/'tmp'),'XDG_CACHE_HOME':str(SSD/'cache'),
         'MPLCONFIGDIR':str(SSD/'cache/matplotlib'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1',
         'MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'}
    proc=None;peak=0;stop=None;teardown=None;samples=[];exc=None
    snapshot=None;final=None;identity=False;outputs={};errors=[];journal=None
    journal_path=Path(str(record_path)+'.failure_journal.jsonl')
    def record_error(label,problem):
        nonlocal stop
        message=label+': '+type(problem).__name__+': '+str(problem)
        errors.append(message);stop=stop or message
    def journal_event(event,**values):
        journal.write(json.dumps(dict(event=event,recorded_utc=utc(),plan_sha256=plan_hash,**values),allow_nan=False)+'\n')
        journal.flush();os.fsync(journal.fileno())
    try:
        # A durable prelaunch record and immediately recorded PID remain even if
        # later resource queries, hashing or primary receipt writing fail.
        journal=journal_path.open('x')
        journal_event('PRELAUNCH',command=command)
        check_termination('BEFORE_WORKER')
        snapshot,reason=limits(plan,started,0)
        if reason:raise RuntimeError(reason)
        with stdout.open('x') as log:
            proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            ownership[0]=False;ownership[1].append(proc)
            journal_event('OWNED_WORKER_LAUNCHED',owned_process_group_id=proc.pid)
            check_termination('AFTER_CHILD_REGISTRATION')
            if os.getpgid(proc.pid)!=proc.pid:raise RuntimeError('OWNED_SESSION_IDENTITY_FAILED')
            while proc.poll() is None:
                rss,pids=monitor.owned_rss(proc.pid);peak=max(peak,rss)
                state,stop=limits(plan,started,rss);samples.append({**state,'owned_pids':pids})
                if TERMINATION_REQUEST:stop='DEFERRED_TERMINATION_DURING_WORKER: '+str(TERMINATION_REQUEST)
                if stop:break
                if len(samples)%30==0:safe_print(json.dumps(dict(worker=record_path.stem,pid=proc.pid,RSS_bytes=rss,
                                                           internal_free_bytes=state['internal_free_bytes'],elapsed_stage_seconds=state['elapsed_stage_seconds'])),flush=True)
                time.sleep(plan['guard']['poll_seconds'])
            check_termination('AFTER_WORKER')
    except BaseException as error:
        exc=error;stop=stop or type(error).__name__+': '+str(error)
    finally:
        # Cleanup is independent of, and precedes, every final metadata query.
        if proc is not None:
            try:
                teardown=monitor.terminate_owned(proc)
                if teardown['initial_group_members'] and stop is None:stop='UNEXPECTED_DESCENDANTS_AFTER_WORKER_EXIT'
                ownership[1].remove(proc);ownership[0]=not ownership[1]
            except BaseException as problem:
                record_error('TEARDOWN_FAILURE',problem);teardown=dict(error=str(problem),verified=False)
        try:
            final,final_reason=limits(plan,started,peak);stop=stop or final_reason
        except BaseException as problem:record_error('FINAL_RESOURCE_QUERY_FAILED',problem)
        try:
            identity=sha(plan_path)==plan_hash
            if not identity:stop=stop or 'FINAL_PLAN_HASH_DIFFERS'
        except BaseException as problem:record_error('FINAL_PLAN_HASH_FAILED',problem)
        try:
            paths=sorted(output_prefix.parent.glob(output_prefix.name+'*'))
            for path in paths+[stdout]:
                try:
                    if path.is_file():outputs[str(path)]=sha(path)
                    elif path==stdout and proc is not None:raise RuntimeError('worker stdout missing')
                except BaseException as problem:record_error('FINAL_OUTPUT_HASH_FAILED '+str(path),problem)
        except BaseException as problem:record_error('FINAL_OUTPUT_INVENTORY_FAILED',problem)
        if journal is not None:
            try:journal_event('WORKER_TEARDOWN_AND_FINAL_QUERIES',owned_process_group_id=proc.pid if proc else None,
                              returncode=proc.returncode if proc else None,stop_reason=stop,metadata_errors=errors,
                              owned_cleanup_verified=ownership[0],process_group_teardown=teardown)
            except BaseException as problem:record_error('FAILURE_JOURNAL_FINAL_APPEND_FAILED',problem)
            try:journal.close()
            except BaseException as problem:record_error('FAILURE_JOURNAL_CLOSE_FAILED',problem)
            try:outputs[str(journal_path)]=sha(journal_path)
            except BaseException as problem:record_error('FAILURE_JOURNAL_HASH_FAILED',problem)
        try:
            final,final_reason=limits(plan,started,peak);stop=stop or final_reason
        except BaseException as problem:record_error('POST_HASH_FINAL_RESOURCE_QUERY_FAILED',problem)
        if TERMINATION_REQUEST:stop=stop or 'DEFERRED_TERMINATION_BEFORE_WORKER_SEAL: '+str(TERMINATION_REQUEST)
        record=dict(status='WORKER_COMPLETE_VERIFIED' if proc and proc.returncode==0 and stop is None and identity else 'WORKER_FAILED_PRESERVED',
                    command=command,plan_sha256=plan_hash,completed_utc=utc(),returncode=proc.returncode if proc else None,
                    owned_process_group_id=proc.pid if proc else None,
                    stop_reason=stop,plan_unchanged=identity,peak_observed_owned_RSS_bytes=peak,resource_preflight=snapshot,
                    resource_final=final,process_group_teardown=teardown,output_sha256=outputs,samples=samples,
                    metadata_errors=errors,failure_journal=str(journal_path),owned_cleanup_verified=ownership[0])
        try:save(record_path,record)
        except BaseException as problem:
            record['status']='WORKER_FAILED_PRESERVED';record['primary_receipt_save_error']=type(problem).__name__+': '+str(problem)
            fallback=Path(str(record_path)+'.failure.json')
            try:save(fallback,record)
            except BaseException as secondary:
                # Physical failure of all storage cannot be repaired by Python.
                # Keep any durable prelaunch/PID journal and emit the essential
                # receipt to stderr; owned cleanup still runs before lock release.
                safe_print(json.dumps(dict(event='FAILURE_RECEIPT_STORAGE_FAILED',record_path=str(record_path),
                                      failure_journal=str(journal_path),owned_process_group_id=proc.pid if proc else None,
                                      primary_error=str(problem),fallback_error=str(secondary),owned_cleanup_verified=ownership[0])),file=sys.stderr,flush=True)
            raise RuntimeError('PRIMARY_WORKER_RECEIPT_SAVE_FAILED_REVIEW_REQUIRED') from problem
    if exc:raise exc
    check_termination('AFTER_WORKER_RECEIPT')
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


def native_result_gate(job,plan,plan_hash,started):
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
    validation=SSD/'receipts_v4'/(job['job_id']+'.numerical_validation.json')
    record=dict(plan_sha256=plan_hash,full_precision_sha256=sha(path),final_intersection_matches_stock_audit=True,
                status='WITHIN_FIT_NUMERICAL_AND_IDENTITY_CONTROLS_PASS',checks=arithmetic,
                control_tolerance_scope='historical within-fit arithmetic only; no tolerance for actual sensitivity differences')
    stage_resource_gate(plan,started,'AFTER_NATIVE_RESULT_HASHES_BEFORE_VALIDATION_SEAL')
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
    recording_errors=[]
    while ownership[1]:
        proc=ownership[1][0]
        try:
            teardown=monitor.terminate_owned(proc)
        except BaseException as error:
            try:safe_print('SHARED_HEAVY_WORKER_LOCK_HELD: owned group '+str(proc.pid)+' remains unresolved: '+str(error),file=sys.stderr,flush=True)
            except BaseException:pass
            try:time.sleep(30)
            except BaseException:pass
            continue
        ownership[1].pop(0)
        try:
            save(SSD/'receipts_v4'/('owned_group_'+str(proc.pid)+'.recovery.json'),
                 dict(completed_utc=utc(),plan_sha256=plan_hash,owned_process_group_id=proc.pid,
                      status='OWNED_GROUP_CLEANUP_VERIFIED_BEFORE_SHARED_LOCK_RELEASE',process_group_teardown=teardown))
        except BaseException as error:
            # A receipt write failure must not skip cleanup of another owned group.
            recording_errors.append(dict(owned_process_group_id=proc.pid,error=type(error).__name__+': '+str(error)))
    ownership[0]=True
    return recording_errors


def execute_under_shared_lock(args):
    check_termination('BEFORE_EXECUTION')
    if sha(args.plan)!=args.plan_sha:raise RuntimeError('PLAN_SHA_DIFFERS')
    plan=json.loads(args.plan.read_text());check_dependencies(plan);check_inputs(plan,plan['input_sha256'])
    if plan['guard']['shared_heavy_worker_lock']!=str(SHARED_HEAVY_WORKER_LOCK):raise RuntimeError('SHARED_LOCK_PATH_DIFFERS')
    baseline=baseline_gate(plan)
    check_termination('AFTER_BASELINE_GATE')
    monitor_path=P/'scripts/30_prepare_and_run_ssd_native_campaign.py'
    spec=importlib.util.spec_from_file_location('_frozen_owned_process_monitor',monitor_path);monitor=importlib.util.module_from_spec(spec);spec.loader.exec_module(monitor)
    # An exclusive lock prevents two sensitivity queues from running concurrently.
    lock=SSD/'sensitivity_executor_v4_5.lock';fd=os.open(str(lock),os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    os.write(fd,(str(os.getpid())+'\n').encode());os.close(fd)
    ownership=[True,[]];completed=[];consumed_output_sha256={}
    try:
        clock_path=SSD/'receipts_v4/stage_clock_v4_5.json'
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
            check_termination('BEFORE_AUDIT_HASHES')
            check_dependencies(plan);check_inputs(plan,audit['inputs'])
            check_termination('AFTER_AUDIT_HASHES')
            cmd=[plan['python'],'-B',str(P/'scripts/37_stock_ldsc_intersection_audit.py'),'--plan',str(args.plan),'--plan-sha',args.plan_sha,'--audit-id',audit['audit_id']]
            worker_path=SSD/'receipts_v4'/(audit['audit_id']+'.worker.json')
            worker(cmd,Path(audit['out_prefix']),worker_path,plan,args.plan,args.plan_sha,started,monitor,ownership)
            consumed_output_sha256[str(worker_path)]=sha(worker_path)
            for path,digest in json.loads(worker_path.read_text())['output_sha256'].items():
                if path in consumed_output_sha256 and consumed_output_sha256[path]!=digest:raise RuntimeError('AUDIT_OUTPUT_REOBSERVATION_CHANGED')
                consumed_output_sha256[path]=digest
            check_dependencies(plan);check_inputs(plan,audit['inputs'])
            stage_resource_gate(plan,started,'AFTER_AUDIT_FINAL_HASHES')
        proof_path=SSD/'proofs/final_stock_intersection_comparison_v4_5.json'
        proof=compare_intersections(plan,args.plan_sha)
        stage_resource_gate(plan,started,'AFTER_INTERSECTION_PROOF_HASHES_BEFORE_SEAL')
        if proof_path.exists():
            old=json.loads(proof_path.read_text())
            if any(old[k]!=proof[k] for k in ['plan_sha256','paired_estimates','identity_instances','proof_sha256','pairs','status']):raise RuntimeError('PRIOR_FINAL_INTERSECTION_PROOF_DIFFERS')
        else:save(proof_path,proof)
        consumed_output_sha256[str(proof_path)]=sha(proof_path)
        for path,digest in proof['proof_sha256'].items():
            if path in consumed_output_sha256 and consumed_output_sha256[path]!=digest:raise RuntimeError('PAIRED_AUDIT_PROOF_CHANGED')
            consumed_output_sha256[path]=digest
        if not args.merge_only:
            for job in plan['jobs']:
                check_termination('BEFORE_FIT_HASHES')
                check_dependencies(plan);check_inputs(plan,job['inputs'])
                check_termination('AFTER_FIT_HASHES')
                cmd=[plan['python'],'-u','-B',str(P/'scripts/38_sensitivity_ldsc_capture.py'),'--plan',str(args.plan),'--plan-sha',args.plan_sha,'--job-id',job['job_id']]+job['ldsc_args']
                worker_path=SSD/'receipts_v4'/(job['job_id']+'.worker.json')
                worker(cmd,Path(job['out_prefix']),worker_path,plan,args.plan,args.plan_sha,started,monitor,ownership)
                consumed_output_sha256[str(worker_path)]=sha(worker_path)
                for path,digest in json.loads(worker_path.read_text())['output_sha256'].items():
                    if path in consumed_output_sha256 and consumed_output_sha256[path]!=digest:raise RuntimeError('FIT_OUTPUT_REOBSERVATION_CHANGED')
                    consumed_output_sha256[path]=digest
                digest=native_result_gate(job,plan,args.plan_sha,started)
                full_precision=job['out_prefix']+'.full_precision.json'
                if consumed_output_sha256.get(full_precision)!=digest:raise RuntimeError('FIT_FULL_PRECISION_CHANGED_AFTER_WORKER_SEAL')
                validation=SSD/'receipts_v4'/(job['job_id']+'.numerical_validation.json')
                consumed_output_sha256[str(validation)]=sha(validation)
                check_dependencies(plan);check_inputs(plan,job['inputs'])
                stage_resource_gate(plan,started,'AFTER_FIT_FINAL_HASHES')
                completed.append(dict(job_id=job['job_id'],estimates=job['estimates'],full_precision_sha256=digest))
                safe_print(json.dumps(dict(native_job_complete=job['job_id'],completed_commands=len(completed))),flush=True)
        check_dependencies(plan);check_inputs(plan,plan['input_sha256'])
        check_termination('AFTER_FINAL_STAGE_HASHES')
        name='merge_only_stage_receipt_v4_5.json' if args.merge_only else 'sensitivity_execution_receipt_v4_5.json'
        stage_record=dict(completed_utc=utc(),plan_sha256=args.plan_sha,baseline_gate_receipt_sha256=baseline,
                          final_intersection_proof_sha256=sha(proof_path),completed_jobs=completed,
                          stage_clock_sha256=sha(clock_path),new_fit_count=sum(x['estimates'] for x in completed),
                          merge_only=args.merge_only,historical_protocol_and_threshold_changes=0,
                          consumed_output_sha256=consumed_output_sha256,
                          elapsed_seconds=time.monotonic()-started)
        stage_record['resource_final']=stage_resource_gate(plan,started,'AFTER_ALL_STAGE_IDENTITY_IO_BEFORE_SUCCESS_SEAL')
        save(SSD/'receipts_v4'/name,stage_record)
        check_termination('AFTER_FINAL_STAGE_SEAL')
        return SSD/'receipts_v4'/name,started
    except BaseException as problem:
        failure=dict(status='STAGE_FAILURE_PRESERVED_NO_AUTOMATIC_RETRY',recorded_utc=utc(),plan_sha256=args.plan_sha,
                     error=type(problem).__name__+': '+str(problem),termination_requests=list(TERMINATION_REQUEST),
                     owned_process_group_ids=[p.pid for p in ownership[1]],completed_jobs=completed)
        try:save(SSD/'receipts_v4/stage_failure_v4_5.json',failure)
        except BaseException as write_error:safe_print(dict(failure,failure_receipt_error=str(write_error)),file=sys.stderr,flush=True)
        raise
    finally:
        if not ownership[0]:
            recording_errors=await_owned_cleanup(monitor,ownership,args.plan_sha)
            if recording_errors:raise RuntimeError('OWNED_CLEANUP_VERIFIED_BUT_RECOVERY_RECEIPTS_FAILED; queue lock retained: '+str(recording_errors))
        lock.unlink()


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,default=SSD/'sensitivity_operational_plan_v4_5.json');p.add_argument('--plan-sha',required=True);p.add_argument('--execute',action='store_true');p.add_argument('--merge-only',action='store_true');args=p.parse_args()
    if not args.execute:raise SystemExit('Prepared only. Explicit --execute is required; no worker launched.')
    prior={getattr(signal,n):signal.signal(getattr(signal,n),catchable_termination) for n in ['SIGINT','SIGTERM','SIGHUP']}
    # Shared with all native and raw-pipeline heavy workers; keep held for the
    # entire execution, including stock merge audits and owned-group teardown.
    # Opening does not truncate or change another worker's lock metadata.
    fd=os.open(str(SHARED_HEAVY_WORKER_LOCK),os.O_RDWR|os.O_CREAT,0o600)
    try:
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('SHARED_HEAVY_WORKER_BUSY_NO_WORKER_LAUNCHED')
        terminal=TerminalCommit(SSD/'receipts_v4/stage_pending_v4_5.json',SSD/'receipts_v4/stage_terminal_seal_v4_5.json',dict(plan_sha256=args.plan_sha,executor_sha256=sha(Path(__file__)),merge_only=args.merge_only))
        target,started=execute_under_shared_lock(args)
        plan=json.loads(args.plan.read_text())
        def final_identity():
            if sha(args.plan)!=args.plan_sha:raise RuntimeError('TERMINAL_PLAN_CHANGED')
            check_dependencies(plan);check_inputs(plan,plan['input_sha256'])
            stage=json.loads(target.read_text());prior=stage['baseline_gate_receipt_sha256']
            if baseline_gate(plan)!=prior:raise RuntimeError('TERMINAL_BASELINE_EVIDENCE_CHANGED')
            for path,digest in stage['consumed_output_sha256'].items():
                f=Path(path)
                if f.is_symlink() or not f.is_file() or sha(f)!=digest:raise RuntimeError('TERMINAL_SENSITIVITY_OUTPUT_CHANGED: '+path)
            expected_commands=plan['stock_merge_only_command_count']+(0 if args.merge_only else plan['native_command_count'])
            workers=[path for path in stage['consumed_output_sha256'] if path.endswith('.worker.json')]
            if len(workers)!=expected_commands:raise RuntimeError('TERMINAL_SENSITIVITY_WORKER_COUNT_DIFFERS')
            for path in workers:
                r=json.loads(Path(path).read_text())
                if r.get('status')!='WORKER_COMPLETE_VERIFIED' or r.get('plan_sha256')!=args.plan_sha or not r.get('owned_cleanup_verified') or r['process_group_teardown']['remaining_group_members']:raise RuntimeError('TERMINAL_SENSITIVITY_WORKER_UNCLEARED')
            if len(stage['completed_jobs'])!=(0 if args.merge_only else 26):raise RuntimeError('TERMINAL_SENSITIVITY_JOB_COUNT_DIFFERS')
        if not terminal.commit({str(target):sha(target)},identity_gate=final_identity,
                               resource_gate=lambda:stage_resource_gate(plan,started,'AFTER_TERMINAL_RECEIPT_PERSISTENCE'),
                               termination_gate=lambda:check_termination('TERMINAL_COMMIT')):
            raise RuntimeError('STAGE_TERMINAL_COMMIT_FAILED_PRESERVED_NO_RETRY')
    finally:
        os.close(fd)
        for signum,handler in prior.items():signal.signal(signum,handler)


if __name__=='__main__':main()
