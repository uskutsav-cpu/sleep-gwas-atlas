#!/usr/bin/env python3
"""Exact controller preflight with metadata stubs: no GWAS/reference read or worker.

The real controller, admission/result gates, TerminalCommit, resource logic,
deferred-signal helper and lock lifetime run against owned invented metadata.
Only lock/flock, baseline gate, source/statistical outputs, worker and monitor
are simulated. No subprocess or actual mutex is constructed.
"""
import contextlib
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1]
FIX=P/'reviews/independent_finngen_operational_controls_v4_3_1'
PLAN=P/'manifests/finngen_preprocessing_plan_v4_3_1.json'
PIN='8e6ffda56823fdbf4f1a7b64d854072383f2732d69c9a02dd465fcb6aee6ed06'
OUT=P/'reviews/independent_finngen_operational_prelaunch_v4_3_1.json'
RECORDS=[]
BINDINGS={}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def bind(p,expected=None):
    h=sha(p);assert expected is None or h==expected,(str(p),h,expected)
    BINDINGS[str(Path(p))]=h
    return h

def new(p,v):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def record(name,good,**d):
    RECORDS.append(dict(name=name,pass_control=bool(good),**d));assert good,RECORDS[-1]

def module(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def metadata():
    bind(PLAN,PIN);plan=json.loads(PLAN.read_text())
    old=P/'manifests/finngen_preprocessing_plan_v4_2.json';bind(old);prior=json.loads(old.read_text())
    record('one_unchanged_preprocessor_zero_estimators',plan['stage']=='preprocessing' and len(plan['jobs'])==1
           and plan['jobs'][0]['command_template'][2]==str(P/'scripts/56_preprocess_finngen_insomnia_feasibility_v2.py')
           and plan['jobs'][0]['command_template'][0]==prior['jobs'][0]['command_template'][0]
           and plan['new_rg_commands']==0)
    record('unchanged_source_identity_N_and_scientific_qualifiers',plan['source_identity']==prior['source_identity']
           and plan['assumed_effective_N']==prior['assumed_effective_N']==185146.70377332723
           and plan['source_generation']==prior['source_generation']=='1777989563097164'
           and plan['scientific_source_admitted'] is False and plan['independent_replication_established'] is False)
    record('unchanged_guards',plan['guard']==prior['guard'] and plan['global_reservation_bytes']==300*2**30)
    deps=plan['dependencies_sha256'];unread=[]
    for path,digest in deps.items():
        f=Path(path)
        metadata_file=f.suffix in ['.py','.json','.md','.sha256'] or (f.suffix=='.tsv' and '/config/' in path)
        executable=f.name=='python' and '/.ldsc-env/bin/' in path
        if metadata_file or executable:bind(f,digest)
        else:unread.append(dict(path=path,expected_sha256=digest))
    for name in ['canonical_calibration_common_v4_5.py','extension_replay_common_v4.py',
                 'native_stage_completion_v4_3.py','terminal_commit_common_v2.py',
                 'sensitivity_executor_v4_4.py','58_run_finngen_feasibility_stage_v3_1.py',
                 '59_prepare_finngen_feasibility_stage_v3_1.py','canonical_calibration_common_v4_3.py']:
        path=str(P/'scripts'/name)
        record('required_direct_or_transitive_code_bound_'+name,path in deps and sha(path)==deps[path])
    for key in ['reference','chain','allele_list','liftover_code']:
        record('scientific_asset_identity_preserved_'+key,plan[key]==prior[key]
               and deps[plan[key]]==prior['dependencies_sha256'][prior[key]])
    required=plan['required_independent_review_paths']
    record('exact_seven_science_ops_successor_review_paths',len(required)==len(set(required))==7
           and required[:4]==prior['required_independent_review_paths']
           and required[4:]==[str(P/'reviews'/n) for n in [
             'independent_finngen_operational_prelaunch_v4_3_1.md','independent_finngen_operational_prelaunch_v4_3_1.json',
             'independent_finngen_operational_prelaunch_seal_v4_3_1.json']])
    proof=json.loads((P/'reviews/independent_finngen_row_science_receipt_v4.json').read_text())
    record('independent_31_science_controls_preserved',proof.get('check_count',proof.get('control_count'))==31
           or len(proof.get('checks',proof.get('controls',[])))==31)
    baseline=json.loads(Path(plan['baseline_gate_plan']).read_text());bind(plan['baseline_gate_plan'],plan['baseline_gate_plan_sha256'])
    counts={s:sum(j['stage']==s for j in baseline['original_190_jobs']) for s in ['core','extension','validation']}
    record('full190_baseline_membership',counts=={'core':57,'extension':112,'validation':21},counts=counts)
    addendum=P/'logs/native190_root_independent_completion_addendum_v4.json';j=json.loads(addendum.read_text())
    for path,digest in j['artifact_sha256'].items():bind(path,digest)
    record('existing_full190_independent_evidence_preserved',j['command_counts']==dict(counts,total=190)
           and j['fully_independent_two_trait_validation_count']==0 and not j['raw_chain_complete'])
    expected={str(P/'scripts/'+Path(path).name):digest for path,digest in []} if False else None
    record('old_unexecuted_v3_and_v2_plan_code_preserved',all(str(P/q) in deps for q in [
           'scripts/58_run_finngen_feasibility_stage_v3.py','scripts/59_prepare_finngen_feasibility_stage_v3.py',
           'manifests/finngen_preprocessing_plan_v4_3.json','scripts/58_run_finngen_feasibility_stage_v2.py']))
    record('actual_admission_and_stage_receipt_absent',not Path(plan['admission']).exists()
           and not Path(plan['stage_receipt']).exists() and not Path(plan['terminal_pending']).exists()
           and not Path(plan['terminal_seal']).exists())
    return plan,unread

class FakeHandle:
    def __init__(self,state):self.state=state
    def __enter__(self):self.state['handle_open']=True;return self
    def fileno(self):return 1933
    def __exit__(self,*a):self.state['handle_open']=False;self.state['released']=True

class FakeLockPath:
    def __init__(self,path,state):self.path=path;self.parent=path.parent;self.state=state
    def __str__(self):return str(self.path)
    def open(self,*a,**k):return FakeHandle(self.state)

def setup(name,fault):
    root=FIX/name;root.mkdir()
    scripts=P/'scripts'
    controller=module(scripts/'58_run_finngen_feasibility_stage_v3_1.py','_finn_controller_'+name)
    common=module(scripts/'canonical_calibration_common_v4_5.py','_finn_common_'+name)
    worker=module(scripts/'sensitivity_executor_v4_4.py','_finn_worker_'+name)
    controller.common=common
    state=dict(handle_open=False,released=False,flock_calls=0,worker_calls=0,cleanup_calls=0,
               persisted=False,post_identity_calls=0,post_resource_calls=0,baseline_calls=0,
               source_reads=0,mutex_open_calls=0,dispatch=None,active_cleanup_fault=False)
    source_root=root/'new_source_feasibility/finngen_R13_F5_INSOMNIA';ns=source_root/'pipeline_replay_v3_1'
    for n in ['tmp','cache','derived','receipts_v4','logs_v4']:(ns/n).mkdir(parents=True)
    source=source_root/'invented_source_metadata.txt';source.write_text('INVENTED CONTROL METADATA, NO GWAS ROWS\n')
    reviewpaths=[]
    for i in range(7):
        r=root/('review_'+str(i)+'.json');new(r,{'invented_review':i});reviewpaths.append(str(r))
    baseline=root/'baseline.json';new(baseline,{'invented_no_fit_baseline':True})
    derivative=ns/'derived/insomnia.sumstats.gz';result=ns/'derived/insomnia.preprocessing.json'
    wr=ns/'receipts_v4/worker.json';stage=root/'stage.json';pp=root/'plan.json';ad=root/'admission.json'
    plan=dict(scope='FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY',new_rg_commands=0,stage='preprocessing',
      namespace=str(ns),source_namespace=str(source_root),shared_heavy_worker_lock=str(root/'native_heavy_worker.lock'),
      source=str(source),source_identity=dict(bytes=source.stat().st_size,sha256=sha(source),md5='invented-control'),
      derivative=str(derivative),admission=str(ad),stage_receipt=str(stage),baseline_gate_plan=str(baseline),
      baseline_gate_plan_sha256=sha(baseline),worker_code='STUB_WORKER',monitor_code='STUB_MONITOR',
      global_reservation_bytes=300*2**30,required_independent_review_paths=reviewpaths,
      terminal_pending=str(ns/'pending.json'),terminal_seal=str(ns/'seal.json'),dependencies_sha256={},
      guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,
                 observed_aggregate_worker_RSS_limit_bytes=2<<30,new_output_limit_bytes=2<<30,deadline_seconds=7200,poll_seconds=2),
      jobs=[dict(command_template=['INVENTED_STUB','{PLAN_SHA256}','{HEAVY_LOCK_FD}'],
                 output_prefix=str(ns/'derived/insomnia'),worker_receipt=str(wr),result_receipt=str(result))])
    new(pp,plan);ph=sha(pp)
    admission=dict(execution_admitted=True,plan_sha256=ph,scope=plan['scope'],independent_binding_review_pass=True,
                   resource_plan_review_pass=True,executor_sha256=sha(scripts/'58_run_finngen_feasibility_stage_v3_1.py'),
                   independent_review_sha256={r:sha(r) for r in reviewpaths})
    if fault=='empty_review':admission['independent_review_sha256']={}
    if fault=='wrong_executor':admission['executor_sha256']='0'*64
    new(ad,admission)
    common.SHARED_LOCK=FakeLockPath(root/'native_heavy_worker.lock',state)
    def fake_flock(fd,flags):
        assert state['handle_open'] and fd==1933;state['flock_calls']+=1
        if fault=='lock_busy':raise BlockingIOError('invented lock busy')
    common.fcntl=SimpleNamespace(flock=fake_flock,LOCK_EX=2,LOCK_NB=4)
    common.time=SimpleNamespace(sleep=lambda _:None)
    # Real deferred handler/context, no real signal is sent to another process.
    controller.physical_mount=lambda:{'invented_mount_control':True}
    real_ns=controller.namespace_bytes
    def nsbytes(path):
        if state['persisted'] and fault=='Finn_cap':return 2**31+1 if str(path)==str(source_root) else real_ns(path)
        if state['persisted'] and fault=='global_cap':return 300*2**30+1 if str(path)==str(root) else real_ns(path)
        return real_ns(path)
    controller.namespace_bytes=nsbytes
    controller.shutil=SimpleNamespace(disk_usage=lambda path:SimpleNamespace(
        free=(3*2**30-1 if state['persisted'] and fault=='internal_floor' and str(path)=='/System/Volumes/Data'
              else 5*2**30-1 if state['persisted'] and fault=='SSD_floor' and str(path)!='//System/Volumes/Data' and str(path)!='/System/Volumes/Data'
              else 100*2**30)))
    clock=[0.0]
    controller.time=SimpleNamespace(monotonic=lambda:clock[0])
    original_write=common.write_new
    def write_hook(path,value):
        if Path(path)==stage:
            if fault=='receipt_write_failure':raise OSError('invented primary receipt write failure')
            original_write(path,value);state['persisted']=True
            if fault=='plan_after_persistence':
                with pp.open('a') as f:f.write(' ')
            if fault=='admission_after_persistence':
                with ad.open('a') as f:f.write(' ')
            if fault=='output_after_persistence':derivative.write_text('CHANGED INVENTED METADATA')
            if fault=='worker_after_persistence':
                with wr.open('a') as f:f.write(' ')
            if fault=='result_after_persistence':
                with result.open('a') as f:f.write(' ')
            if fault=='source_after_persistence':source.write_text('CHANGED INVENTED SOURCE METADATA')
            if fault=='SIGTERM':common.catchable_termination(signal.SIGTERM,None)
            if fault=='SIGHUP':common.catchable_termination(signal.SIGHUP,None)
            if fault=='SIGINT':common.catchable_termination(signal.SIGINT,None)
            if fault=='deadline':clock[0]=7201.0
        else:original_write(path,value)
    common.write_new=write_hook
    def baseline_gate(_):
        state['baseline_calls']+=1
        if state['persisted']:state['post_identity_calls']+=1
        if fault=='baseline_incomplete':raise RuntimeError('invented 189 baseline gate failure')
        return {'invented_FULL190_PROOF':'changed' if state['persisted'] and fault=='baseline_after_persistence' else 'fixed'}
    worker.baseline_gate=baseline_gate
    original_cleanup=worker.await_owned_cleanup
    def cleanup(monitor,ownership,plan_hash):
        state['cleanup_calls']+=1
        if state['active_cleanup_fault']:
            assert state['handle_open'], 'logical lock released around unresolved cleanup'
            if state['cleanup_calls']<4:raise KeyboardInterrupt('invented cleanup interruption')
            state['active_cleanup_fault']=False;ownership[:]=[True,[]]
        return original_cleanup(monitor,ownership,plan_hash)
    worker.await_owned_cleanup=cleanup
    def fake_worker(command,prefix,record_path,_,plan_path,plan_hash,started,monitor,ownership):
        state['worker_calls']+=1;state['dispatch']=command
        assert state['handle_open'] and worker.subprocess.fd==1933 and command[-1]=='1933'
        if fault=='worker_failure':raise RuntimeError('invented worker failed')
        if fault=='cleanup_exception':
            state['active_cleanup_fault']=True;ownership[:]=[False,[SimpleNamespace(pid=1001)]]
            raise RuntimeError('invented owned worker requires cleanup')
        derivative.write_text('INVENTED DERIVATIVE METADATA, NO SNP ROWS\n')
        source_id=plan['source_identity']
        r=dict(plan_sha256=ph,status='QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS',retained_minimum_diagnostic_pass=True,
               source_before=source_id,source_after=source_id,derivative=str(derivative),derivative_sha256=sha(derivative),
               full_source_gzip_CRC_and_EOF_verified=True,derivative_full_gzip_CRC_and_EOF_verified=True,
               scientific_source_admitted=False,independent_replication_established=False,
               verified_per_variant_N=False,verified_per_variant_INFO=False)
        if fault=='false_CRC':r['full_source_gzip_CRC_and_EOF_verified']=False
        if fault=='source_admitted':r['scientific_source_admitted']=True
        if fault=='replication_claim':r['independent_replication_established']=True
        if fault=='fabricated_INFO':r['verified_per_variant_INFO']=True
        if fault=='fabricated_N':r['verified_per_variant_N']=True
        new(result,r)
        new(wr,dict(status='WORKER_COMPLETE_VERIFIED',plan_sha256=ph,owned_cleanup_verified=True,
                    process_group_teardown={'remaining_group_members':[]},output_sha256={str(derivative):sha(derivative),str(result):sha(result)}))
    worker.worker=fake_worker
    controller.load=lambda name,path:worker if path=='STUB_WORKER' else SimpleNamespace()
    outcome=None
    try:controller.execute(pp,ph)
    except BaseException as e:outcome=type(e).__name__+': '+str(e)
    pending=Path(plan['terminal_pending']);seal=Path(plan['terminal_seal'])
    committed=not pending.exists() and seal.exists() and not Path(str(stage)+'.failure.json').exists()
    if fault=='healthy':
        record(name,committed and outcome is None and state['worker_calls']==1 and state['released']
               and state['post_identity_calls']>=1,virtual_state=state)
    elif fault in ['empty_review','wrong_executor','lock_busy']:
        record(name,not committed and state['worker_calls']==0,virtual_state=state,outcome=outcome)
    else:
        record(name,not committed and pending.exists() and state['released'] and outcome is not None,
               virtual_state=state,outcome=outcome,PENDING_veto_preserved=True)
    assert not state['handle_open']

def proxy_control():
    c=module(P/'scripts/58_run_finngen_feasibility_stage_v3_1.py','_finn_proxy')
    real=c.subprocess;calls=[]
    fake=SimpleNamespace(Popen=lambda *a,**k:calls.append((a,k)),STDOUT='invented')
    c.subprocess=fake;proxy=c.InheritedMutexSubprocess(1933)
    proxy.Popen(['INVENTED_NO_PROCESS'],start_new_session=True)
    record('exact_local_descriptor_proxy_no_worker',calls[0][1]['pass_fds']==(1933,)
           and real.Popen is not fake.Popen and proxy.STDOUT=='invented')
    rejected=False
    try:proxy.Popen([],pass_fds=(2,))
    except RuntimeError:rejected=True
    record('descriptor_override_rejected_no_worker',rejected and len(calls)==1)

def main():
    t=time.monotonic();assert not FIX.exists() and not OUT.exists();FIX.mkdir()
    sys.path.insert(0,str(P/'scripts'))
    _,unread=metadata();proxy_control()
    for fault in ['healthy','empty_review','wrong_executor','lock_busy','baseline_incomplete','worker_failure',
                  'false_CRC','source_admitted','replication_claim','fabricated_INFO','fabricated_N',
                  'receipt_write_failure','cleanup_exception','plan_after_persistence','admission_after_persistence',
                  'worker_after_persistence','result_after_persistence','output_after_persistence','source_after_persistence',
                  'baseline_after_persistence','SIGTERM','SIGHUP','SIGINT','internal_floor','SSD_floor','Finn_cap','global_cap','deadline']:
        setup('controller_'+fault,fault)
    for path,h in BINDINGS.items():assert sha(path)==h,('metadata changed',path)
    fixtures={str(p):sha(p) for p in FIX.rglob('*') if p.is_file()}
    assert sum(p.stat().st_size for p in FIX.rglob('*') if p.is_file())<8*2**20 and time.monotonic()-t<120
    result=dict(schema='independent_finngen_controller_prelaunch_v4_3_1',verdict='QUALIFIED_OPERATIONAL_PRELAUNCH_PASS',
                exact_plan_sha256=PIN,controls=RECORDS,control_count=len(RECORDS),all_control_expectations_pass=True,
                consumed_code_metadata_before_after_sha256=BINDINGS,excluded_real_body_dependencies=unread,
                fixture_metadata_sha256=fixtures,actual_workers=0,actual_mutex_operations=0,production_body_reads=0,
                real_source_decompressions=0,fits=0,execution_admission_granted=False,
                full190_gate_rerun=False,simulated_controller_controls_are_not_actual_preprocessing_or_science=True,
                qualifications=['absent INFO and measured per-SNP N','assumed constant N_eff','clinical rights and construct not admitted',
                                'cohort independence not established','not replication','private immutable PENDING/deferred-signal contract',
                                'periodic guards can overshoot; no kernel quota','host SIGKILL outside Python cleanup'],
                elapsed_seconds=time.monotonic()-t)
    new(OUT,result);print(json.dumps(dict(verdict=result['verdict'],controls=len(RECORDS),receipt_sha256=sha(OUT),elapsed=result['elapsed_seconds'])))

if __name__=='__main__':main()
