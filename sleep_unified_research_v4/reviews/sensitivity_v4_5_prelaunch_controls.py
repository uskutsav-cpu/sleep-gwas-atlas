#!/usr/bin/env python3
"""Narrow independent metadata-only sensitivity v4_5 controls; no actual workers."""
import ast
import collections
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode=True
P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4')
PLAN=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities/sensitivity_operational_plan_v4_5.json')
PIN='152cddd729cef3f9d4a148e7c527a59251cd73528fbb2867bac1fcdca8d1a4de'
CODE=P/'scripts/sensitivity_executor_v4_5.py'
CODEPIN='50d8c449b2b3d781bf2fed2cc163b9ca3ae8b375444aa9309b628be6a1b0454e'
FIX=P/'reviews/sensitivity_v4_5_prelaunch_controls'
RECORDS=[]
sys.path.insert(0,str(P/'scripts'))

def digest(q):return hashlib.sha256(Path(q).read_bytes()).hexdigest()
def new(q,x):
    with Path(q).open('x') as f:json.dump(x,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def module(q):
    s=importlib.util.spec_from_file_location('_sensitivity_private_'+str(len(RECORDS)),q)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def record(case,good,**detail):
    RECORDS.append(dict(case=case,control_pass=bool(good),**detail));assert good,RECORDS[-1]
def route(v):
    if isinstance(v,str):return v.replace('/fits_v4_3/','/fits_v4_5/').replace('/intersections_v4_3/','/intersections_v4_5/')
    if isinstance(v,list):return [route(x) for x in v]
    if isinstance(v,dict):return {k:route(x) for k,x in v.items()}
    return v
def defs(q):
    return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(q).read_text()).body if isinstance(n,ast.FunctionDef)}

def static_checks():
    assert digest(PLAN)==PIN and digest(CODE)==CODEPIN
    n=json.loads(PLAN.read_text());oldpath=Path(n['preserved_v4_3_1_plan']);o=json.loads(oldpath.read_text())
    assert digest(oldpath)=='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
    record('exact_scientific_argv_and_audit_route',n['jobs']==route(o['jobs']) and n['audit_jobs']==route(o['audit_jobs']))
    fields=['input_sha256','derivative_proof','derivative_proof_sha256','environment','python','ldsc_dir','reference_prefix',
            'original_190_jobs','baseline_input_sha256','baseline_dependency_sha256','baseline_relocated_dependency_sha256',
            'scientific_member_manifest_sha256','new_fitted_estimates','native_command_count','stock_merge_only_command_count',
            'stock_merge_only_intersection_identities','core_precision_adjudication']
    record('old_science_and_baseline_bindings_unchanged',all(n[k]==o[k] for k in fields))
    record('old_dependencies_retained',all(n['dependencies_sha256'].get(k)==v for k,v in o['dependencies_sha256'].items()))
    g=copy.deepcopy(n['guard']);g.pop('global_reservation_bytes')
    record('old_guards_unchanged_plus_global300',g==o['guard'] and n['guard']['global_reservation_bytes']==300<<30)
    sm=json.loads(Path(n['scientific_member_manifest']).read_text())
    assert digest(Path(n['scientific_member_manifest']))==n['scientific_member_manifest_sha256']=='3ce7b50bb4e1289ef5a69d8e46dcef759963703ef7b317e156bdee0d4c40e5e1'
    simplified=[{k:j[k] for k in ['job_id','kind','inputs','options','estimates']} for j in n['jobs']]
    record('frozen_26_62_matrix',simplified==sm['jobs'] and len(simplified)==26 and sum(j['estimates'] for j in simplified)==62)
    record('frozen52_124_audits',len(n['audit_jobs'])==52 and sum(a['expected_identities'] for a in n['audit_jobs'])==124
           and len({a['audit_id'] for a in n['audit_jobs']})==52)
    dp=json.loads(Path(n['derivative_proof']).read_text());assert digest(n['derivative_proof'])==n['derivative_proof_sha256']
    record('unchanged_finite_N_derivative_proof',dp['status']=='RESULT_FREE_ONLY_FINITE_N_CHANGED'
           and [r['new_finite_N'] for r in dp['results']]==[376169,290130]
           and all(r['source_sha256_before']==r['source_sha256_after'] and r['original_and_derivative_invariant_match']
             and r['literal_missingness_preserved'] and r['counts']['finite_N_changed_rows']+r['counts']['nonfinite_N_preserved_rows']==1217311
             and r['full_stream_verified_rows']==1217311 for r in dp['results']))
    a=defs(P/'scripts/sensitivity_executor_v4_3_1.py');b=defs(CODE)
    inherited=['compare_intersections','native_result_gate','worker','await_owned_cleanup','catchable_termination','check_termination']
    record('inherited_science_worker_controls_AST',all(a[k]==b[k] for k in inherited),functions=inherited)
    helper=module(P/'scripts/native_stage_completion_v4_3.py')
    record('exact_baseline_stage_selection',helper.stage_monitor_path('extension').name=='extension_native_monitor_receipt_v4_3.json'
           and helper.stage_monitor_path('core').name=='core_native_monitor_receipt_v4.json'
           and helper.stage_monitor_path('validation').name=='validation_native_monitor_receipt_v4.json')
    ledger=json.loads(Path(n['global_resource_ledger_path']).read_text())
    record('ledger3_arithmetic',ledger==n['reservation_arithmetic'] and digest(n['global_resource_ledger_path'])==n['global_resource_ledger_sha256']
           and sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes']
           and ledger['reserved_total_bytes']+ledger['unallocated_margin_bytes']==300<<30,
           reserved_total_bytes=ledger['reserved_total_bytes'],margin_bytes=ledger['unallocated_margin_bytes'])
    checked={};excluded={}
    for path,wanted in n['dependencies_sha256'].items():
        if '/ref/' in path and not path.endswith('.provenance.json'):
            excluded[path]='reference body not read';continue
        q=Path(path);assert q.stat().st_size<=12<<20,path
        assert digest(q)==wanted,path;checked[path]=wanted
    record('metadata_dependency_recheck',bool(checked),checked_count=len(checked),reference_bodies_excluded=len(excluded))
    return dict(verified_metadata_code_sha256=checked,excluded_reference_sha256=excluded,
        scientific_input_sha256_not_rehashed=n['input_sha256'],new_estimates=62,native_commands=26,audit_commands=52,
        derivative_proof_sha256=n['derivative_proof_sha256'],old_plan_sha256=digest(oldpath),
        historical190_job_counts=dict(collections.Counter(j['stage'] for j in n['original_190_jobs'])))

def setup(case):
    d=FIX/case;d.mkdir();m=module(CODE);m.SSD=d/'SSD/sensitivities';m.SSD.mkdir(parents=True)
    m.SHARED_HEAVY_WORKER_LOCK=m.SSD.parent/'native_heavy_worker.lock'
    for sub in ['receipts_v4','logs_v4','proofs','fits_v4_5','tmp','cache']:(m.SSD/sub).mkdir()
    p=copy.deepcopy(json.loads(PLAN.read_text()));p['input_sha256']={};p['dependencies_sha256']={}
    p['guard']['shared_heavy_worker_lock']=str(m.SHARED_HEAVY_WORKER_LOCK)
    for group in ['jobs','audit_jobs']:
        for j in p[group]:j['out_prefix']=str(m.SSD/('fits_v4_5' if group=='jobs' else 'proofs')/(j.get('job_id') if group=='jobs' else j['audit_id']))
    pp=d/'fixture_plan.json';new(pp,p);ph=digest(pp)
    m.check_dependencies=lambda p:None;m.check_inputs=lambda p,i:None
    m.baseline_gate=lambda p:{'MOCK190_METADATA_ONLY':'not scientific evidence'}
    m.safe_print=lambda *a,**k:None
    class Disk:
        def disk_usage(self,path):return SimpleNamespace(free=1<<50)
    m.shutil=Disk()
    return d,m,p,pp,ph

def full_route_control(case):
    d,m,p,pp,ph=setup(case);calls=[]
    def worker(cmd,prefix,rp,plan,plan_path,plan_hash,started,monitor,ownership):
        kind='audit' if '--audit-id' in cmd else 'fit';ident=cmd[cmd.index('--audit-id' if kind=='audit' else '--job-id')+1];calls.append((kind,ident))
        if kind=='fit':assert len([x for x in calls if x[0]=='audit'])==52
        output=Path(str(prefix)+('.intersection.json' if kind=='audit' else '.full_precision.json'))
        if kind=='audit':
            audit=next(a for a in plan['audit_jobs'] if a['audit_id']==ident);job=next(j for j in plan['jobs'] if j['job_id']==audit['sensitivity_job_id'])
            records=[]
            for i in range(job['estimates']):
                same=job['job_id'].startswith('lipid_two_step_');arm=audit['arm']
                records.append({'final_ordered_SNP_count':1,'final_ordered_SNP_sha256':'a'*64,
                    'final_ordered_SNP_allele_compatibility_sha256':'b'*64,'final_ordered_aligned_Z1_Z2_float64_big_endian_sha256':'c'*64,
                    'final_ordered_Z_float64_big_endian_sha256':'c'*64,'final_ordered_N1_N2_float64_big_endian_sha256':('d' if same or arm=='baseline' else 'e')*64,
                    'final_ordered_N_float64_big_endian_sha256':('d' if same or arm=='baseline' else 'e')*64,
                    'two_step_hsq_mask':{'fixture_only':True},'p2':'METADATA_FIXTURE_'+str(i),'input':'METADATA_FIXTURE_'+str(i),
                    'fixture_only_not_scientific_evidence':True})
            new(output,{'plan_sha256':ph,'estimator_calls':0,'final_intersections':records,'fixture_only_not_scientific_evidence':True})
        else:new(output,{'fixture_only_not_scientific_evidence':True,'job_id':ident})
        journal=Path(str(rp)+'.failure_journal.jsonl');journal.write_text('{"fixture_only":true}\n')
        new(rp,{'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':ph,'owned_cleanup_verified':True,
                'process_group_teardown':{'remaining_group_members':[]},'command':cmd,
                'output_sha256':{str(output):digest(output),str(journal):digest(journal)}})
    m.worker=worker
    def native(job,plan,planhash,started):
        q=Path(job['out_prefix']+'.full_precision.json')
        if case=='capture_changed_after_worker':q.write_text('CHANGED_METADATA_CAPTURE\n')
        new(m.SSD/'receipts_v4'/(job['job_id']+'.numerical_validation.json'),{'fixture_only_not_scientific_evidence':True})
        return digest(q)
    m.native_result_gate=native
    args=SimpleNamespace(plan=pp,plan_sha=ph,merge_only=False)
    try:target,started=m.execute_under_shared_lock(args);success=True
    except RuntimeError:success=False
    if case=='full_route':
        s=json.loads(target.read_text());record('all52_audits_precede26_stub_fits',success and [k for k,i in calls]==['audit']*52+['fit']*26
            and s['new_fit_count']==62 and len(s['completed_jobs'])==26
            and len([q for q in s['consumed_output_sha256'] if q.endswith('.worker.json')])==78,
            metadata_stub_calls=len(calls),actual_workers=0,actual_estimators=0)
        proof=json.loads((m.SSD/'proofs/final_stock_intersection_comparison_v4_5.json').read_text())
        record('actual_compare_function_metadata_only124',proof['identity_instances']==124 and proof['paired_estimates']==62)
    else:record(case,not success and len(calls)==53,metadata_stub_calls=len(calls),actual_workers=0,actual_estimators=0)

def terminal_control(case):
    d,m,p,pp,ph=setup('terminal_'+case);call_count=[0];baseline_calls=[0];roles={};cleanup_locked=[False]
    target=m.SSD/'receipts_v4/sensitivity_execution_receipt_v4_5.json'
    def lockheld():
        fd=os.open(m.SHARED_HEAVY_WORKER_LOCK,os.O_RDWR)
        try:
            try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);return False
            except BlockingIOError:return True
        finally:os.close(fd)
    def execute(args):
        call_count[0]+=1;cleanup_locked[0]=lockheld();consumed={};completed=[]
        for a in p['audit_jobs']:
            out=Path(a['out_prefix']+'.intersection.json');new(out,{'fixture_only_not_scientific_evidence':True});consumed[str(out)]=digest(out);roles.setdefault('audit',out)
        proof=m.SSD/'proofs/final_stock_intersection_comparison_v4_5.json';new(proof,{'fixture_only_not_scientific_evidence':True});consumed[str(proof)]=digest(proof);roles['proof']=proof
        for i in range(78):
            w=m.SSD/'receipts_v4'/f'fixture_{i:02d}.worker.json';j=Path(str(w)+'.failure_journal.jsonl');j.write_text('{"fixture_only":true}\n')
            r={'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':ph,'owned_cleanup_verified':True,
                'process_group_teardown':{'remaining_group_members':[]},'output_sha256':{str(j):digest(j)}}
            if i==0:
                if case=='uncleared_cleanup':r['owned_cleanup_verified']=False
                if case=='remaining_group':r['process_group_teardown']['remaining_group_members']=[{'metadata_placeholder_no_process':True}]
                if case=='worker_bad_status':r['status']='WORKER_FAILED_PRESERVED'
                if case=='worker_bad_plan':r['plan_sha256']='f'*64
                roles['worker']=w;roles['journal']=j
            new(w,r);consumed[str(w)]=digest(w);consumed[str(j)]=digest(j)
        for j in p['jobs']:
            q=Path(j['out_prefix']+'.full_precision.json');new(q,{'fixture_only_not_scientific_evidence':True});consumed[str(q)]=digest(q);roles.setdefault('capture',q)
            completed.append({'job_id':j['job_id'],'estimates':j['estimates'],'full_precision_sha256':digest(q)})
        if case=='missing_worker':consumed.pop(str(roles['worker']))
        if case=='extra_worker':
            q=m.SSD/'receipts_v4/extra.worker.json';new(q,{'fixture_only_not_scientific_evidence':True});consumed[str(q)]=digest(q)
        if case=='missing_job':completed.pop()
        if case=='extra_job':completed.append(completed[-1])
        new(target,{'baseline_gate_receipt_sha256':{'MOCK190_METADATA_ONLY':'not scientific evidence'},'completed_jobs':completed,
            'consumed_output_sha256':consumed,'fixture_only_not_scientific_evidence':True})
        return target,time.monotonic()
    m.execute_under_shared_lock=execute
    def baseline(plan):
        baseline_calls[0]+=1
        assert lockheld()
        if baseline_calls[0]==2:
            if case.startswith('mutate_'):roles[case[len('mutate_'):]].write_text('POST_PERSISTENCE_METADATA_MUTATION\n')
            if case=='symlink_capture':
                q=roles['capture'];clone=d/'capture_clone';clone.write_bytes(q.read_bytes());q.unlink();q.symlink_to(clone)
            if case=='deferred_signal':m.catchable_termination(signal.SIGTERM,None)
        return {'MOCK190_METADATA_ONLY':'not scientific evidence'}
    m.baseline_gate=baseline
    gate=m.stage_resource_gate
    def resource(plan,start,where):
        if case=='postpersist_resource' and baseline_calls[0]>=2:raise RuntimeError('FIXTURE_POST_PERSIST_RESOURCE')
        return gate(plan,start,where)
    m.stage_resource_gate=resource
    tm=module(P/'scripts/terminal_commit_common_v2.py');save=tm.save_new
    def terminal_write(path,value):
        if case!='success' and str(path).endswith('.failure.json'):raise OSError('FIXTURE_SUPPLEMENTAL_WRITE_FAILURE')
        return save(path,value)
    tm.save_new=terminal_write;m.TerminalCommit=tm.TerminalCommit
    args=[str(CODE),'--execute','--plan',str(pp),'--plan-sha',ph]
    oldargv=sys.argv;sys.argv=args
    try:
        try:m.main();success=True;error=None
        except BaseException as e:success=False;error=type(e).__name__+': '+str(e)
    finally:sys.argv=oldargv
    consumer=module(P/'scripts/terminal_commit_common_v2.py')
    pending=m.SSD/'receipts_v4/stage_pending_v4_5.json';seal=m.SSD/'receipts_v4/stage_terminal_seal_v4_5.json'
    binding={'plan_sha256':ph,'executor_sha256':CODEPIN,'merge_only':False}
    try:consumer.require_committed(pending,seal,binding,{str(target):digest(target)});accepted=True
    except BaseException:accepted=False
    fd=os.open(m.SHARED_HEAVY_WORKER_LOCK,os.O_RDWR)
    try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);released=True
    except BlockingIOError:released=False
    finally:os.close(fd)
    expect=case=='success'
    record('terminal_'+case,success==expect and accepted==expect and released and cleanup_locked[0],
        controller_success=success,consumer_accepts=accepted,pending_exists_or_symlink=pending.exists() or pending.is_symlink(),
        private_shared_mutex_held_during_mock_execution=cleanup_locked[0],private_mutex_released_after_return=released,
        supplemental_write_failure_injected=case!='success',error=error,actual_worker_calls=0)
    if expect:
        sys.argv=args+['--merge-only']
        try:
            try:m.main();rejected=False
            except RuntimeError:rejected=True
        finally:sys.argv=oldargv
        record('second_split_invocation_rejected',rejected and call_count[0]==1,actual_worker_calls=0)

def resource_control():
    d,m,p,pp,ph=setup('global_guard')
    original=m.limits;g=p['guard'];m.SSD.parent.joinpath('global_marker_metadata').write_text('x')
    gp=copy.deepcopy(p);gp['guard']['global_reservation_bytes']=0
    state,reason=original(gp,time.monotonic(),0)
    record('new_global_meter_and_guard',reason=='300GIB_CAMPAIGN_GUARD' and state['new_campaign_namespace_bytes']>0)

def main():
    FIX.mkdir(exist_ok=False);metadata=static_checks()
    full_route_control('full_route');full_route_control('capture_changed_after_worker')
    for case in ['success','mutate_capture','mutate_worker','mutate_journal','mutate_audit','mutate_proof','symlink_capture',
                 'missing_worker','extra_worker','missing_job','extra_job','uncleared_cleanup','remaining_group',
                 'worker_bad_status','worker_bad_plan','postpersist_resource','deferred_signal']:
        terminal_control(case)
    resource_control()
    assert digest(PLAN)==PIN and digest(CODE)==CODEPIN
    for q,h in metadata['verified_metadata_code_sha256'].items():assert digest(q)==h,q
    files={str(q):digest(q) for q in FIX.rglob('*') if q.is_file() and not q.is_symlink()}
    links={str(q):os.readlink(q) for q in FIX.rglob('*') if q.is_symlink()}
    out=P/'reviews/sensitivity_v4_5_prelaunch_controls_receipt.json'
    new(out,{'schema':'independent_sensitivity_v4_5_metadata_prelaunch_controls','plan_sha256':PIN,'executor_sha256':CODEPIN,
        'control_script_sha256':digest(__file__),'control_count':len(RECORDS),'all_controls_pass':all(x['control_pass'] for x in RECORDS),
        'controls':RECORDS,'metadata':metadata,'fixture_regular_file_sha256':files,'fixture_symlink_targets':links,
        'mock_scientific_proofs_are_not_outcomes':True,'GWAS_body_reads':0,'reference_body_reads':0,'real_data_decompressions':0,
        'actual_workers':0,'actual_fits':0,'actual_stock_merge_audits':0,'execution_admission_granted':False})
    print(json.dumps({'control_count':len(RECORDS),'receipt_sha256':digest(out),'actual_workers':0,'actual_fits':0},indent=2))

if __name__=='__main__':main()
