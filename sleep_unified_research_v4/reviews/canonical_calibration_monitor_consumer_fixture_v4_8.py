"""Author tiny actual-consumer checks only; no worker, body, mutex or admission."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import struct
import sys

R=Path(__file__).resolve().parent;P=R.parent;S=P/'scripts'
F=R/'canonical_calibration_monitor_consumer_fixtures_v4_8'
OUT=R/'canonical_calibration_monitor_consumer_fixture_receipt_v4_8.json'
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_8 as c
spec=importlib.util.spec_from_file_location('author_canonical8',S/'42_prepare_native_canonical_calibration_v4_8.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
checks=[];cases=[]


def require(value,label):
    checks.append({'label':label,'pass':bool(value)})
    if not value:raise AssertionError(label)


def save(path,record):
    path.parent.mkdir(parents=True,exist_ok=True)
    c.write_intended_new(path,record,c.intended_sha(record))


def companion(path):
    with path.open('xb') as f:
        f.write(struct.pack('>II16sH',0x00051607,0x00020000,b'\0'*16,1)+struct.pack('>III',9,38,4)+b'TEST')


def base(name,apple=False):
    root=F/name/'native_canonical_calibration_v4_8';out=root/'native'
    (out/'weighted_capture').mkdir(parents=True)
    m.PLAN=root.parent/'plan.json';m.ADMISSION=root.parent/'admission.json'
    save(m.PLAN,{'private_fixture':name});save(m.ADMISSION,{'private_fixture':name})
    job={'job_id':'AUTHOR_NO_WORKER','output_dir':str(out)}
    plan={'jobs':[job],'ssd_output_root':str(root)};command=['NO_PROCESS','METADATA_ONLY']
    cap=out/'native_capture_receipt.json';save(cap,{'job':job,'plan_sha256':c.sha(m.PLAN)})
    monitor=out/'monitor_receipt.json'
    record={'job_id':job['job_id'],'plan_sha256':c.sha(m.PLAN),'admission_sha256':c.sha(m.ADMISSION),
      'command':command,'returncode':0,'stop_reason':None,'process_group_teardown':{'remaining_group_members':[]},
      'first_teardown_attempt':{'signal_zero_group_absent':True},'run_error':None,'cleanup_errors':[],
      'final_guard_violations':[],'owned_group_disappearance_verified':True,
      'bulk_output_hashes_completed_before_final_guards':True,'output_hash_error':None,
      'deferred_termination_requests':[],'all_output_sha256':{str(cap):c.sha(cap)}}
    if apple:companion(out/'._monitor_receipt.json')
    return {'root':root,'out':out,'job':job,'plan':plan,'command':command,'cap':cap,'monitor':monitor,'record':record,
            'evidence':{'workers':{},'all_output_sha256':{}}}


def mutate(record,name):
    if name=='cleanup_errors':record['cleanup_errors']=['EXPLICIT_PROOF_FAILURE']
    elif name=='run_error':record['run_error']={'type':'RuntimeError','repr':'EXPLICIT_FAILURE'}
    elif name=='teardown_cleanup_error':record['process_group_teardown']['cleanup_error']='EXPLICIT_FAILURE'
    elif name=='first_teardown_exception':record['first_teardown_attempt']['exception']={'type':'RuntimeError','repr':'EXPLICIT_FAILURE'}
    elif name=='first_teardown_exception_null_present':record['first_teardown_attempt']['exception']=None


def freeze(f):
    save(f['monitor'],f['record'])
    m.freeze_worker_evidence(f['plan'],c.sha(m.ADMISSION),f['job']['job_id'],f['command'],f['out'],f['evidence'],f['job'],
                            {'path':str(f['monitor']),'sha256':c.sha(f['monitor'])})


def final_consume(f):
    # Supply an already-frozen contradictory record to isolate the second
    # consumer. No post-freeze mutation or ordinary producer path is claimed.
    save(f['monitor'],f['record'])
    mapping={str(f['cap']):c.sha(f['cap']),str(f['monitor']):c.sha(f['monitor'])}
    side=f['out']/'._monitor_receipt.json'
    if side.exists():mapping[str(side)]=c.sha(side)
    f['evidence']={'workers':{f['job']['job_id']:{'output_dir':str(f['out']),'sha256':dict(mapping),
                    'registered_directories':[str(f['out']),str(f['out']/'weighted_capture')],
                    'external_output_paths':[],'registered_deleted_markers':[]}},'all_output_sha256':dict(mapping)}
    original_admit=m.admit;original_worker=m.monitored_worker;calls=[]
    class HealthyReachedArithmeticWithoutLaunch(Exception):pass
    def no_worker(*args,**kwargs):
        calls.append('arithmetic_boundary');raise HealthyReachedArithmeticWithoutLaunch('NO_PROCESS_LAUNCHED')
    m.admit=lambda:(f['plan'],{})
    m.monitored_worker=no_worker
    try:
        try:m.verify(-1,[],f['evidence']);return None,calls
        except BaseException as e:return e,calls
    finally:m.admit=original_admit;m.monitored_worker=original_worker


def main():
    require(not F.exists() and not OUT.exists(),'fresh_owned_fixture_namespace');F.mkdir()
    names=['42_prepare_native_canonical_calibration','canonical_calibration_common','canonical_calibration_capture','canonical_calibration_arithmetic']
    sources=[S/(n+'_v4_'+str(v)+'.py') for v in [7,8] for n in names]
    before={str(p):c.sha(p) for p in sources}
    for path in sources:compile(path.read_text(),str(path),'exec')
    fields=['cleanup_errors','run_error','teardown_cleanup_error','first_teardown_exception']
    for name in ['healthy','registered_AppleDouble']+fields+['first_teardown_exception_null_present']:
        f=base('initial_'+name,name=='registered_AppleDouble');mutate(f['record'],name)
        try:freeze(f);error=None
        except BaseException as e:error=repr(e)
        healthy=name in ['healthy','registered_AppleDouble']
        require((error is None)==healthy,'initial_actual_consumer_'+name)
        if not healthy:
            require('CANONICAL_MONITOR_RECORDED_FAILURE' in error,'initial_failure_oracle_'+name)
            require(not f['evidence']['workers'] and not f['evidence']['all_output_sha256'],'initial_no_evidence_credited_'+name)
        else:m.check_frozen_worker_evidence(f['evidence'])
        cases.append({'path':'freeze_worker_evidence','case':name,'error':error,'pass':True,'ordinary_producer_failure_claim':False})
    for name in ['healthy','registered_AppleDouble']+fields+['first_teardown_exception_null_present']:
        f=base('final_'+name,name=='registered_AppleDouble');mutate(f['record'],name)
        error,calls=final_consume(f)
        healthy=name in ['healthy','registered_AppleDouble']
        require(bool(calls)==healthy,'final_actual_consumer_arithmetic_boundary_'+name)
        if healthy:require('HealthyReachedArithmeticWithoutLaunch' in type(error).__name__,'final_healthy_no_worker_'+name)
        else:require('CANONICAL_MONITOR_RECORDED_FAILURE' in repr(error),'final_failure_oracle_'+name)
        cases.append({'path':'verify','case':name,'error':repr(error),'stub_calls':calls,'pass':True,
                      'isolated_already_frozen_contradictory_metadata':not healthy,'actual_workers':0})
    def function(path,name):
        return next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
    def normalized(node):
        return ast.dump(node,include_attributes=False)
    require((S/'canonical_calibration_capture_v4_7.py').read_text().replace('v4_7','v4_8')==(S/'canonical_calibration_capture_v4_8.py').read_text(),
            'capture_import_namespace_only')
    require(normalized(function(S/'canonical_calibration_arithmetic_v4_7.py','validate'))==normalized(function(S/'canonical_calibration_arithmetic_v4_8.py','validate')),
            'independent_arithmetic_validate_exact_AST')
    for name in ['freeze_worker_evidence','verify']:
        old=function(S/'42_prepare_native_canonical_calibration_v4_7.py',name)
        new=function(S/'42_prepare_native_canonical_calibration_v4_8.py',name)
        class RemoveVeto(ast.NodeTransformer):
            count=0
            def visit_Expr(self,node):
                if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='assert_monitor_error_fields_clear':
                    self.count+=1;return None
                return self.generic_visit(node)
        remover=RemoveVeto();new=remover.visit(new)
        require(remover.count==1,'one_shared_failure_oracle_each_consumer_'+name)
        # Only declared namespace string constants differ outside the insertion.
        class OldNamespace(ast.NodeTransformer):
            def visit_Constant(self,node):
                if isinstance(node.value,str):node.value=node.value.replace('v4_8','v4_7')
                return node
        new=OldNamespace().visit(new)
        require(normalized(old)==normalized(new),'consumer_AST_only_shared_veto_and_namespace_'+name)
    for name in ['monitored_worker','monitor_helpers','snapshot','output_bytes','owned_group_exists','await_owned_cleanup','final_violations','check_frozen_worker_evidence','final_evidence_inventory']:
        old=function(S/'42_prepare_native_canonical_calibration_v4_7.py',name)
        new=function(S/'42_prepare_native_canonical_calibration_v4_8.py',name)
        require(normalized(old)==normalized(new),'producer_ownership_resource_inventory_AST_unchanged_'+name)
    for name in ['catchable_termination','deferred_termination_signals','safe_diagnostic','exclusive_heavy_lock','validate_inherited_lock','historical_completion','admit','intended_sha','write_intended_new','terminal_intended_sha','freeze_companions','canonical_epoch_bytes','preserved_failed_attempt']:
        require(normalized(function(S/'canonical_calibration_common_v4_7.py',name))==normalized(function(S/'canonical_calibration_common_v4_8.py',name)),
                'common_terminal_history_transport_guard_AST_unchanged_'+name)
    # Preparation-only map checking is tested without any dependency body IO.
    d={str(i):str(i) for i in range(50)};fake={'dependencies_sha256':{**d,'fresh_code':'HASH'}};seen=[];orig=m.check_hashes
    m.check_hashes=lambda x:seen.append(x)
    try:
        m.preparation_dependency_hashes(fake,d)
        require(seen==[{'fresh_code':'HASH'}],'prep_exact50_deferred_not_fresh_hashed')
        bad=dict(d);bad['0']='MUTATED'
        try:m.preparation_dependency_hashes(fake,bad);error=None
        except BaseException as e:error=repr(e)
        require(error is not None and len(seen)==1,'prep_deferred_map_mismatch_fails_before_any_hash')
        cases.extend([{'case':'preparation_exact50_deferred','pass':True},{'case':'preparation_map_mismatch','error':error,'pass':True}])
    finally:m.check_hashes=orig
    require(before=={str(p):c.sha(p) for p in sources},'old_and_new_sources_unchanged_by_checks')
    receipt={'schema':'author_canonical_v4_8_monitor_consumer_fixture','status':'AUTHOR_NARROW_SANITY_PASS_NOT_INDEPENDENT_ADMISSION',
      'cases':cases,'case_count':len(cases),'checks':checks,'check_count':len(checks),'all_checks_pass':all(x['pass'] for x in checks),
      'source_sha256':before,'workers_launched':0,'mutex_acquisitions':0,'GWAS_reference_runtime_body_reads':0,
      'ordinary_producer_false_success_demonstrated':False,'qualification':'Both real consumer functions run on tiny metadata. Healthy final consumer stops at an explicit nonlaunch stub. Future full16+1 execution and independent review remain required.'}
    save(OUT,receipt)
    print(json.dumps({'receipt':str(OUT),'sha256':c.sha(OUT),'cases':len(cases),'checks':len(checks),'all_checks_pass':True}))


if __name__=='__main__':main()
