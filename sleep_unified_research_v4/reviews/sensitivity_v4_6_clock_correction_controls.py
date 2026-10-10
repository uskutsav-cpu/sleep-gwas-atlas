#!/usr/bin/env python3
"""Narrow v4_6 diff and actual controller clock gates; tiny metadata stubs only."""
import copy
import fcntl
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode=True
P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4')
R=P/'reviews'
S=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
CODE=P/'scripts/sensitivity_executor_v4_6.py'
CODEPIN='118cc27c4e882c66be29714b1786d9171c1d204215ce4feaa00cf95c878783df'
PLAN=S/'sensitivity_operational_plan_v4_6.json'
PLANPIN='b8f9c4c0dfecf4f5f4559df3f468b1f2128ef932caece54691124a8788432efb'
PREP=P/'scripts/73_prepare_estimator_sensitivity_v4_6.py'
PREPPIN='8549321f27e46f559999cbd784b3b93da3c0ed21192e14bb997cfb8945ecf48a'
BASE=R/'sensitivity_v4_5_prelaunch_controls.py'
BASEPIN='1a9b55687e44a18757e28857338813f7b5cdce0b2a0f8a66948e0391a25822ad'
FIX=R/'sensitivity_v4_6_clock_correction_controls'
RECORDS=[]

def sha(q):return hashlib.sha256(Path(q).read_bytes()).hexdigest()
def new(q,x):
    with Path(q).open('x') as f:
        json.dump(x,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def record(case,good,**detail):
    RECORDS.append(dict(case=case,control_pass=bool(good),**detail));assert good,RECORDS[-1]
def route(x):
    if isinstance(x,str):return x.replace('/fits_v4_5/','/fits_v4_6/').replace('/intersections_v4_5/','/intersections_v4_6/')
    if isinstance(x,list):return [route(i) for i in x]
    if isinstance(x,dict):return {k:route(v) for k,v in x.items()}
    return x

def main():
    assert sha(CODE)==CODEPIN and sha(PLAN)==PLANPIN and sha(PREP)==PREPPIN and sha(BASE)==BASEPIN
    before=P/'scripts/sensitivity_executor_v4_5.py'
    old=before.read_text();current=CODE.read_text()
    added='        consumed_output_sha256[str(clock_path)]=sha(clock_path)\n'
    record('only_clock_insertion_and_operational_epoch_diff',current.count(added)==1 and current.replace('_v4_6','_v4_5').replace(added,'')==old)
    record('clock_bound_before_elapsed_and_audits',current.index(added)<current.index('        elapsed=time.time()-clock')<current.index("        for audit in plan['audit_jobs']:"))
    a=json.loads((S/'sensitivity_operational_plan_v4_5.json').read_text());b=json.loads(PLAN.read_text())
    record('exact_science_argv_and_matrix_inherited',b['jobs']==route(a['jobs']) and b['audit_jobs']==route(a['audit_jobs']) and len(b['jobs'])==26 and sum(j['estimates'] for j in b['jobs'])==62 and len(b['audit_jobs'])==52 and sum(j['expected_identities'] for j in b['audit_jobs'])==124)
    allowed={'jobs','audit_jobs','dependencies_sha256','executor','schema','prepared_utc','resource_preflight','preserved_unexecuted_v4_5_plan','preserved_unexecuted_v4_5_plan_sha256','clock_binding_correction'}
    record('all_remaining_scientific_runtime_guard_bindings_equal',all(a.get(k)==b.get(k) for k in set(a)|set(b) if k not in allowed))
    removed={str(P/'scripts/sensitivity_executor_v4_5.py'),str(P/'scripts/73_prepare_estimator_sensitivity_v4_5.py')}
    newdeps={str(CODE):CODEPIN,str(PREP):PREPPIN,str(S/'sensitivity_operational_plan_v4_5.json'):'152cddd729cef3f9d4a148e7c527a59251cd73528fbb2867bac1fcdca8d1a4de',str(R/'sensitivity_v4_5_stage_clock_witness_receipt_v2.json'):'033116a587eaa577bdfa5fc3e7400ca663896d3132445d17a18fc1b3f7a7e964'}
    expected={k:v for k,v in a['dependencies_sha256'].items() if k not in removed};expected.update(newdeps)
    record('exact_dependency_relocation_and_preserved_plan_witness',b['dependencies_sha256']==expected and all(sha(k)==v for k,v in newdeps.items()) and b['preserved_unexecuted_v4_5_plan_sha256']==newdeps[str(S/'sensitivity_operational_plan_v4_5.json')])
    priorseal=R/'sensitivity_v4_5_independent_prelaunch_review_seal.json'
    assert sha(priorseal)=='1af2859afca77c4e1c388dc3905a69a8d44a076fff55da6400248e74c90d10df'
    frozen=json.loads(priorseal.read_text())
    record('preserved_v4_5_eleven_file_seal_rechecked',all(sha(k)==v for k,v in frozen['file_sha256'].items()),files_verified=len(frozen['file_sha256']))
    receipt=json.loads((R/'sensitivity_v4_5_prelaunch_controls_receipt.json').read_text())
    record('old33_control_evidence_inherited_not_rerun',receipt['control_count']==33 and receipt['all_controls_pass'] and all(x['control_pass'] for x in receipt['controls']) and receipt['actual_workers']==receipt['actual_fits']==receipt['actual_stock_merge_audits']==0)
    assert sha(P/'scripts/terminal_commit_common_v2.py')=='9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd'
    spec=importlib.util.spec_from_file_location('_narrow_clock_base',BASE);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    base.CODE=CODE;base.CODEPIN=CODEPIN;base.PLAN=PLAN;base.PIN=PLANPIN;base.FIX=FIX
    FIX.mkdir(exist_ok=False)
    setup_source=inspect.getsource(base.setup).replace('_v4_5','_v4_6')
    exec(compile(setup_source,str(__file__)+':private_setup','exec'),base.__dict__)
    factory_source=inspect.getsource(base.full_route_control)
    token='    args=SimpleNamespace(plan=pp,plan_sha=ph,merge_only=False)\n'
    assert factory_source.count(token)==1
    factory_source=factory_source.split(token)[0].replace('def full_route_control(case):','def route_factory(case):')+'    return d,m,p,pp,ph,calls\n'
    exec(compile(factory_source,str(__file__)+':metadata_stub_factory','exec'),base.__dict__)
    for case in ['success','mutate_clock_after_persistence']:
        d,m,p,pp,ph,calls=base.route_factory(case)
        clock=m.SSD/'receipts_v4/stage_clock_v4_6.json'
        target=m.SSD/'receipts_v4/sensitivity_execution_receipt_v4_6.json'
        pending=m.SSD/'receipts_v4/stage_pending_v4_6.json';seal=m.SSD/'receipts_v4/stage_terminal_seal_v4_6.json'
        baseline_calls=[0];locked=[True];mutation=[None];failed_writes=[]
        def lockheld():
            fd=os.open(m.SHARED_HEAVY_WORKER_LOCK,os.O_RDWR)
            try:
                try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);return False
                except BlockingIOError:return True
            finally:os.close(fd)
        def baseline(plan):
            baseline_calls[0]+=1;locked[0]=locked[0] and lockheld()
            if baseline_calls[0]==3 and case!='success':
                assert seal.is_file() and pending.is_file()
                state=json.loads(target.read_text());before_hash=sha(clock)
                assert state['consumed_output_sha256'][str(clock)]==state['stage_clock_sha256']==before_hash
                clock.write_text('POST_PERSISTENCE_METADATA_CLOCK_MUTATION\n')
                mutation[0]=dict(seal_persisted_before_mutation=True,pending_owned_before_mutation=True,recorded_sha256=before_hash,mutated_sha256=sha(clock))
            return {'MOCK190_METADATA_ONLY':'not scientific evidence'}
        m.baseline_gate=baseline
        original_worker=m.worker
        def metadata_worker(*args,**kw):
            locked[0]=locked[0] and lockheld()
            return original_worker(*args,**kw)
        m.worker=metadata_worker
        terminal=base.module(P/'scripts/terminal_commit_common_v2.py');original_save=terminal.save_new
        def save(path,value):
            if case!='success' and str(path).endswith('.failure.json'):
                failed_writes.append(str(path));raise OSError('FIXTURE_SUPPLEMENTAL_WRITE_FAILURE')
            return original_save(path,value)
        terminal.save_new=save;m.TerminalCommit=terminal.TerminalCommit
        original_argv=sys.argv;sys.argv=[str(CODE),'--execute','--plan',str(pp),'--plan-sha',ph]
        try:
            try:m.main();success=True;error=None
            except BaseException as e:success=False;error=type(e).__name__+': '+str(e)
        finally:sys.argv=original_argv
        consumer=base.module(P/'scripts/terminal_commit_common_v2.py')
        try:consumer.require_committed(pending,seal,dict(plan_sha256=ph,executor_sha256=CODEPIN,merge_only=False),{str(target):sha(target)});accepted=True
        except BaseException as e:accepted=False;consumer_error=type(e).__name__+': '+str(e)
        state=json.loads(target.read_text())
        record(case,success==(case=='success') and accepted==(case=='success') and len(calls)==78 and [x[0] for x in calls]==['audit']*52+['fit']*26 and len(state['completed_jobs'])==26 and state['new_fit_count']==62 and str(clock) in state['consumed_output_sha256'] and locked[0] and not lockheld() and (case=='success' or (pending.is_file() and failed_writes and mutation[0] is not None and mutation[0]['recorded_sha256']!=mutation[0]['mutated_sha256'])),controller_success=success,consumer_accepts=accepted,controller_error=error,consumer_error=None if accepted else consumer_error,baseline_metadata_callbacks=baseline_calls[0],metadata_stub_calls=len(calls),actual_worker_calls=0,actual_estimator_calls=0,pending_retained=pending.is_file(),provisional_terminal_seal_preserved=seal.is_file(),failed_supplemental_writes=failed_writes,mutation=mutation[0],private_shared_mutex_held_through_callbacks=locked[0],private_mutex_released_after_return=not lockheld())
    for q,h in [(CODE,CODEPIN),(PLAN,PLANPIN),(PREP,PREPPIN),(BASE,BASEPIN),(priorseal,'1af2859afca77c4e1c388dc3905a69a8d44a076fff55da6400248e74c90d10df')]:assert sha(q)==h
    out=R/'sensitivity_v4_6_clock_correction_controls_receipt.json'
    new(out,dict(schema='independent_sensitivity_v4_6_narrow_clock_controls',executor_sha256=CODEPIN,plan_sha256=PLANPIN,prepare_sha256=PREPPIN,script_sha256=sha(__file__),inherited_v4_5_controls_receipt_sha256='085334ec51c2d71be9771a81bc0524be8a984b198142b38aa1178e17dd003b87',inherited_controls=33,narrow_control_count=len(RECORDS),all_controls_pass=all(x['control_pass'] for x in RECORDS),controls=RECORDS,private_setup_source_sha256=hashlib.sha256(setup_source.encode()).hexdigest(),metadata_worker_factory_source_sha256=hashlib.sha256(factory_source.encode()).hexdigest(),fixture_regular_file_sha256={str(q):sha(q) for q in FIX.rglob('*') if q.is_file() and not q.is_symlink()},GWAS_body_reads=0,reference_body_reads=0,real_data_decompressions=0,actual_workers=0,actual_fits=0,actual_stock_merge_audits=0,execution_admission_granted=False))
    print(json.dumps(dict(narrow_controls=len(RECORDS),receipt_sha256=sha(out),actual_workers=0,actual_fits=0),indent=2))

if __name__=='__main__':main()
