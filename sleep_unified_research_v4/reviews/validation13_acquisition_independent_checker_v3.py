#!/usr/bin/env python3
"""Metadata and tiny synthetic controls only: no curl, flock or GWAS/ref bodies."""
import ast
import base64
import contextlib
import copy
import datetime
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import types
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / 'sleep_unified_research_v4'
SCRIPT = PKG / 'scripts/90_acquire_validation_raw_sources_v3.py'
PLAN = PKG / 'manifests/validation_raw_acquisition_plan_v4_3.json'
DESIGN = PKG / 'statistical_validation/validation_raw_replay_design_v1.json'
SEAL = PKG / 'statistical_validation/validation_raw_replay_design_seal_v1.json'
OUT = Path(__file__).with_name('validation13_acquisition_independent_controls_v3')
RECEIPT = Path(__file__).with_name('validation13_acquisition_independent_receipt_v3.json')
EXPECTED_PLAN = '8e93d548163bda284b6292b9780962060bbc3ff4afde44e5a916c80f3b69b708'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(65536), b''):
            h.update(block)
    return h.hexdigest()

def save(p,v):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f:
        json.dump(v,f,indent=2); f.write('\n')

def metadata():
    assert sha(PLAN) == EXPECTED_PLAN
    plan=json.loads(PLAN.read_text()); design=json.loads(DESIGN.read_text()); seal=json.loads(SEAL.read_text())
    assert len(design['metadata_dependencies_sha256']) == 53
    deps={**design['metadata_dependencies_sha256'], **seal['file_sha256'], **plan['bound_sources']}
    observed={}
    for name,digest in deps.items():
        p=Path(name)
        # Only explicitly listed metadata/code/runtime identities are read.
        assert p == Path('/usr/bin/curl') or p.suffix in {'.json','.tsv','.py','.txt','.md'}, name
        assert p.is_file() and not p.is_symlink() and not any(x.is_symlink() for x in p.parents), name
        assert p.stat().st_size <= 20*1024**2, name
        observed[name]=sha(p); assert observed[name] == digest, name
    assert [x['index'] for x in plan['members']] == list(range(1,14))
    assert [x['source_id'] for x in plan['members']] == [x['source_id'] for x in design['members']]
    assert len(set(x['source_id'] for x in plan['members'])) == 13
    assert sum(x['expected_bytes'] for x in plan['members']) == plan['compressed_network_bytes'] == 10058648185
    finn=mvp=0
    for m,d in zip(plan['members'],design['members']):
        h=d['historical_source_body']; q=d['HEAD_probe']; fields=q['headers']; t=m['expected_transport_headers']
        assert (m['url'],m['expected_bytes'],m['expected_md5'],m['expected_sha256']) == (h['source'],h['observed_size_bytes'],h['observed_md5'],h['observed_sha256'])
        old=json.loads(Path(d['original_receipt_path']).read_text())
        assert h==old['source_verification']
        # Design carries the immutable historical body proof, itself bound to original receipt bytes.
        assert d['original_receipt_sha256'] == sha(d['original_receipt_path'])
        assert q['status']==200 and int(fields['content-length'])==m['expected_bytes']
        assert t['etag']==fields['etag'].strip('"') and t['last-modified']==fields['last-modified']
        assert m['body_path']==d['proposed_body_path'].replace('/validation_raw_replay_v1/','/validation_raw_replay_v3/')
        if m['source_id'].startswith('finngen_'):
            finn+=1
            assert '?generation='+t['x-goog-generation'] in m['url']
            assert q['identity_checks']=={'etag':True,'generation':True,'length':True}
            assert int(t['x-goog-stored-content-length'])==m['expected_bytes']
            assert base64.b64decode(t['x-goog-hash'].split('md5=',1)[1]).hex()==m['expected_md5']
        else:
            mvp+=1
            assert d['current_official_yaml_MD5_if_MVP']==m['expected_md5']
            y=PKG/'source_provenance'/('sleep_clinician_'+m['source_id'].removeprefix('gwas_catalog_')+'_yaml_primary_metadata_v4.txt')
            assert re.search(r'^data_file_md5sum:\s*'+m['expected_md5']+r'\s*$',y.read_text(),re.M)
            assert q['identity_checks']['etag'] is False and 'differs' in m['transport_qualification']
    assert (finn,mvp)==(11,2)
    raw100=json.loads((PKG/'manifests/extension_raw_acquisition_plan_v4_4.json').read_text())
    assert plan['exclusive_family_lock_path']==raw100['exclusive_family_lock_path']
    assert (plan['internal_floor_bytes'],plan['ssd_floor_bytes'],plan['maximum_owned_transfer_rss_bytes'])==(3<<30,5<<30,2<<30)
    assert (plan['per_body_seconds_limit'],plan['family_seconds_limit'],plan['SSD_reservation_bytes'])==(7200,96*3600,300<<30)
    for key in ['automatic_retry','source_substitution','scientific_membership_changes','independent_two_trait_replication_claimed','human_or_restricted_access_requested']:
        assert plan[key] is False
    assert plan['source_bodies_read_in_preparation']==plan['workers_launched']==0
    fns=lambda p:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
    a,b=fns(SCRIPT),fns(PKG/'scripts/52_acquire_extension_raw_sources_v8.py')
    immutable=['utc','safe_print','catchable_termination','hashes','write_new','family_lock','load_monitor','physical_mount','global_namespace_gate','assert_bindings','header_fields','safe_state','cleanup_proof','retain_lock_until_gone']
    assert all(a[k]==b[k] for k in immutable)
    prior_maps={}; prior_seals={}
    for version,sealdigest,count in [(1,'d40c0524acfe924c3ae86cbdfec1fe8d7df30afee3c5aa32cb377e6361223de9',7),(2,'0f599f9da1f8af6a595b5b57667d306cd5241cb1a8ee88f85067496998079197',6)]:
        oldseal=PKG/f'reviews/validation13_acquisition_independent_review_seal_v{version}.json'
        assert sha(oldseal)==sealdigest
        previous=json.loads(oldseal.read_text()); actual_map=previous['review_artifact_sha256']
        assert len(actual_map)==count
        for name,digest in actual_map.items():
            assert sha(name)==digest and plan['bound_sources'][name]==digest
        assert plan['bound_sources'][str(oldseal)]==sealdigest
        prior_maps[str(version)]=actual_map;prior_seals[str(version)]=sealdigest
    oldplan=json.loads((PKG/'manifests/validation_raw_acquisition_plan_v4_2.json').read_text())
    invariant_keys=['member_count','compressed_network_bytes','internal_floor_bytes','ssd_floor_bytes','SSD_reservation_bytes','maximum_owned_transfer_rss_bytes','monitor_output_limit_bytes','per_body_seconds_limit','family_seconds_limit','runtime_poll_seconds','exclusive_family_lock_path','curl_version','global_resource_ledger_path','global_resource_ledger_sha256','automatic_retry','source_substitution','scientific_membership_changes','independent_two_trait_replication_claimed','shared_single_transfer_family','human_or_restricted_access_requested']
    assert all(plan[k]==oldplan[k] for k in invariant_keys)
    for old,new in zip(oldplan['members'],plan['members']):
        assert {k:v for k,v in old.items() if k!='body_path'}=={k:v for k,v in new.items() if k!='body_path'}
        assert new['body_path']==old['body_path'].replace('/validation_raw_replay_v2/','/validation_raw_replay_v3/')
    oldfns=fns(PKG/'scripts/82_acquire_validation_raw_sources_v2.py')
    assert a['assert_headers']==oldfns['assert_headers']
    normalized_acquire=a['acquire'].replace('validation_raw_acquisition_v4_3','validation_raw_acquisition_v4_2')
    assert normalized_acquire==oldfns['acquire']
    assert sha(SCRIPT)=='8d0765930a346d1097d15a7b5fc3bea294f782b4eb8f52a98f894bbf693152cc'
    return plan,{'dependencies_checked':len(observed),'design_metadata_dependencies':53,'dependency_sha256':observed,'Finn_generation_and_MD5_bindings':finn,'MVP_original_body_and_official_yaml_MD5_bindings':mvp,'immutable_helper_AST_matches':immutable,'prior_review_seal_sha256':prior_seals,'prior_review_artifact_sha256':prior_maps,'all13_prior_review_artifacts_fresh_hash_match_and_in_plan':True,'source_member_fields_unchanged_except_private_body_route':True,'unchanged_transport_function_AST_after_route_normalization':True,'unchanged_guard_keys':invariant_keys}

def prepare_reader_controls(plan):
    # Execute the actual preparation reader only: suppress all mkdir/write/mount/runtime actions.
    sys.path.insert(0,str(PKG/'scripts'))
    ns={'__name__':'metadata_only_prepare_reader','__file__':str(SCRIPT)}
    exec(compile(ast.parse(SCRIPT.read_text()),str(SCRIPT),'exec'),ns)
    ns.update(physical_mount=lambda:None,load_monitor=lambda:types.SimpleNamespace(snapshot=lambda:{'internal_free_bytes':10<<30,'ssd_free_bytes':500<<30}))
    originals={v:(PKG/f'reviews/validation13_acquisition_independent_review_seal_v{v}.json').read_text() for v in (1,2)}
    original_loads=json.loads
    results=[]
    for case in ['healthy','empty_v1_map','wrong_v1_artifact_digest','empty_v2_map','wrong_v2_artifact_digest']:
        writes=[];mkdirs=[]
        def load_hook(value,*args,**kwargs):
            parsed=original_loads(value,*args,**kwargs)
            for v,text in originals.items():
                if value==text and case.endswith(str(v)+'_map'):
                    parsed['review_artifact_sha256']={}
                if value==text and case==f'wrong_v{v}_artifact_digest':
                    first=next(iter(parsed['review_artifact_sha256']));parsed['review_artifact_sha256'][first]='0'*64
            return parsed
        ns['write_new']=lambda path,value:writes.append((str(path),copy.deepcopy(value)))
        error=None
        with patch('json.loads',load_hook),patch.object(Path,'mkdir',lambda self,*a,**k:mkdirs.append(str(self))),patch('subprocess.run',lambda *a,**k:types.SimpleNamespace(stdout=plan['curl_version'])),contextlib.redirect_stdout(io.StringIO()):
            try:ns['prepare']()
            except BaseException as exc:error=type(exc).__name__+': '+str(exc)
        if case=='healthy':
            assert error is None and len(writes)==2
            assert writes[0][1]['bound_sources']==plan['bound_sources'] and writes[1][1]==writes[0][1]
            assert writes[0][1]['members']==plan['members']
        else:assert error is not None and not writes
        results.append({'case':case,'error':error,'captured_plan_writes':len(writes),'suppressed_mkdir_calls':len(mkdirs),'real_writes_or_network_or_locks':0,'actual_prepare_reader_executed':True})
    return results

class ScaleSyntheticByteBudget(ast.NodeTransformer):
    def __init__(self,total):self.total=total;self.changed=0
    def visit_Constant(self,node):
        if node.value==10058648185:
            self.changed+=1; return ast.copy_location(ast.Constant(self.total),node)
        return node

def control(plan,case):
    base=OUT/case; base.mkdir(parents=True)
    pkg=base/'PKG'; ssd=base/'SSD'; folder=ssd/'validation_raw_replay_v3'
    for p in [pkg/'logs',pkg/'source_provenance/validation_raw_acquisition_v4_3',folder/'raw',folder/'logs',folder/'receipts',ssd/'tmp']:
        p.mkdir(parents=True)
    payloads=[('SYNTHETIC_FIXTURE_'+str(i)+'\n').encode() for i in range(1,14)]
    scaled=sum(map(len,payloads)); tree=ast.parse(SCRIPT.read_text()); transform=ScaleSyntheticByteBudget(scaled);tree=ast.fix_missing_locations(transform.visit(tree))
    # Original function bodies, import/definition path only; __main__ is never executed.
    sys.path.insert(0,str(PKG/'scripts'))
    ns={'__name__':'metadata_only_synthetic_fixture','__file__':str(SCRIPT)}
    exec(compile(tree,str(SCRIPT),'exec'),ns)
    assert transform.changed==6, transform.changed
    p=copy.deepcopy(plan)
    p.update(private_namespace=str(folder),compressed_network_bytes=scaled,bound_sources={},
             executor_path=str(SCRIPT),pending_path=str(folder/'pending.json'),terminal_seal_path=str(folder/'seal.json'),
             exclusive_family_lock_path=str(ssd/'DO_NOT_OPEN_OR_LOCK'),curl_version='SYNTHETIC_CURL_IDENTITY_ONLY')
    for m,data in zip(p['members'],payloads):
        m['body_path']=str(folder/'raw'/m['filename']);m['expected_bytes']=len(data)
        m['expected_md5']=hashlib.md5(data).hexdigest();m['expected_sha256']=hashlib.sha256(data).hexdigest()
        if 'x-goog-stored-content-length' in m['expected_transport_headers']:
            m['expected_transport_headers']['x-goog-stored-content-length']=str(len(data))
            m['expected_transport_headers']['x-goog-hash']='md5='+base64.b64encode(hashlib.md5(data).digest()).decode()
    planpath=base/'fixture_plan.json';save(planpath,p);planhash=sha(planpath)
    reviewmd=base/'review.md'; reviewmd.write_text('SYNTHETIC_NOT_OPERATIONAL_ADMISSION\n')
    reviewjson=base/'review.json';save(reviewjson,{'synthetic':True})
    admission=base/'fixture_admission.json';save(admission,{'execution_admitted':True,'plan_sha256':planhash,'executor_sha256':sha(SCRIPT),'independent_review_artifact_sha256':{str(reviewmd):sha(reviewmd),str(reviewjson):sha(reviewjson)}})
    state={'internal_free_bytes':10<<30,'ssd_free_bytes':500<<30}
    monitor=types.SimpleNamespace(INTERNAL_FLOOR=3<<30,SSD_FLOOR=5<<30,RSS_LIMIT=2<<30,OUTPUT_LIMIT=4<<30,
        snapshot=lambda:dict(state),final_limits=lambda *args:None,group_members=lambda pid:[],
        terminate_owned=lambda proc:{'worker_launched':True,'remaining_group_members':[],'teardown_verified':True})
    spies={'fake_processes':0,'actual_subprocesses':0,'actual_flocks':0,'quarantine_spy_calls':0,'mutation_applied':False}
    ns.update(PACKAGE=pkg,SSD=ssd,FOLDER=folder,PLAN=planpath,ADMISSION=admission,
              load_monitor=lambda:monitor,physical_mount=lambda:None,
              family_lock=lambda path:contextlib.nullcontext(),TERMINATION_REQUEST=[])
    class FakeProc:
        pid=19000001;returncode=0
        def poll(self):return self.returncode
        def wait(self,timeout=None):return self.returncode
    def fake_popen(command,**kwargs):
        spies['fake_processes']+=1
        partial=Path(command[command.index('--output')+1]);header=Path(command[command.index('--dump-header')+1])
        m=next(m for m in p['members'] if str(partial)==m['body_path']+'.partial'); data=payloads[m['index']-1]
        partial.write_bytes(data)
        h={'content-length':str(len(data)),**m['expected_transport_headers']}
        if case=='wrong_transport_header' and m['index']==1:h['etag']='WRONG'
        header.write_text('HTTP/1.1 200 OK\n'+''.join(k+': '+v+'\n' for k,v in h.items()))
        return FakeProc()
    original_write=ns['write_new']; family_paths=[pkg/'logs/validation_raw_acquisition_family_receipt_v4_3.json',folder/'validation_raw_acquisition_family_receipt_v4_3.json']
    intended={};seal_intended={}
    def write_hook(path,value):
        path=Path(path)
        if path in family_paths:intended[str(path)]=hashlib.sha256((json.dumps(value,indent=2)+'\n').encode()).hexdigest()
        if case=='cleanup_and_receipt_persistence_failure' and path.parent==folder/'receipts':raise OSError('SYNTHETIC_RECEIPT_FAILURE')
        original_write(path,value)
        if path==family_paths[1] and case in ('mutate_both_before_hash_capture','mutate_one_before_hash_capture'):
            for item in family_paths if case=='mutate_both_before_hash_capture' else family_paths[:1]:
                altered=json.loads(item.read_text()); altered.update(status='MUTATED_BEFORE_FIRST_HASH_CAPTURE',completed_source_count=12)
                item.write_text(json.dumps(altered,indent=2)+'\n')
            spies['mutation_applied']=True
        if path==family_paths[1] and case=='missing_primary_before_hash_capture':
            family_paths[0].unlink();spies['mutation_applied']=True
    ns['write_new']=write_hook
    if case=='cleanup_and_receipt_persistence_failure':
        ns['cleanup_proof']=lambda *a:{'teardown_verified':False,'remaining_group_members':[19000001],'cleanup_error':'SYNTHETIC'}
        ns['retain_lock_until_gone']=lambda *a:spies.__setitem__('quarantine_spy_calls',spies['quarantine_spy_calls']+1)
    terminal_module=sys.modules['terminal_commit_common_v2']; original_save=terminal_module.save_new
    def terminal_save(path,value):
        if Path(path)==Path(p['terminal_seal_path']):
            seal_intended['sha256']=hashlib.sha256((json.dumps(value,indent=2,allow_nan=False)+'\n').encode()).hexdigest()
            seal_intended['pending_sha256']=value['pending_sha256']
        original_save(path,value)
        if Path(path)!=Path(p['terminal_seal_path']):return
        if case=='samebytes_repo_logs_parent_after_seal':
            old=pkg/'logs'; retained=pkg/'logs_retained';old.rename(retained);old.symlink_to(retained, target_is_directory=True)
        elif case=='samebytes_primary_leaf_after_seal':
            old=family_paths[0];retained=old.with_suffix('.retained.json');old.rename(retained);old.symlink_to(retained)
        elif case=='source_body_after_seal':Path(p['members'][0]['body_path']).write_bytes(b'ALTERED_SYNTHETIC')
        elif case=='source_mirror_after_seal':
            q=pkg/'source_provenance/validation_raw_acquisition_v4_3'/(p['members'][0]['source_id']+'.json');q.write_text(q.read_text()+' ')
        elif case=='source_failure_after_seal':
            q=folder/'receipts'/(p['members'][0]['source_id']+'.json.failure.json');save(q,{'synthetic_failure':True})
        elif case=='family_failure_after_seal':save(Path(str(family_paths[0])+'.failure.json'),{'synthetic_failure':True})
        elif case=='deferred_signal_after_seal':ns['TERMINATION_REQUEST'].append(15)
        elif case=='resource_failure_after_seal':state['internal_free_bytes']=0
        elif case=='primary_corruption_after_seal':family_paths[0].write_text(family_paths[0].read_text()+' ')
        elif case=='missing_source_body_after_seal':Path(p['members'][0]['body_path']).unlink()
        elif case=='missing_primary_after_seal':family_paths[0].unlink()
        elif case=='pending_leaf_symlink_after_seal':
            old=Path(p['pending_path']);retained=old.with_suffix('.retained.json');old.rename(retained);old.symlink_to(retained)
        elif case=='seal_leaf_symlink_after_seal':
            old=Path(p['terminal_seal_path']);retained=old.with_suffix('.retained.json');old.rename(retained);old.symlink_to(retained)
        elif case=='pending_corruption_after_seal':
            old=Path(p['pending_path']);old.write_text(old.read_text()+' ')
        elif case in ['seal_content_before_first_hash_capture','seal_binding_before_first_hash_capture','seal_primary_map_before_first_hash_capture','seal_qualification_before_first_hash_capture','seal_whitespace_before_first_hash_capture']:
            old=Path(p['terminal_seal_path']);altered=json.loads(old.read_text())
            if case=='seal_content_before_first_hash_capture':altered['pending_sha256']='0'*64
            elif case=='seal_binding_before_first_hash_capture':altered['binding']['executor_sha256']='0'*64
            elif case=='seal_primary_map_before_first_hash_capture':altered['result_receipt_sha256'][str(family_paths[0])]='0'*64
            elif case=='seal_qualification_before_first_hash_capture':altered['success_requires_absent_PENDING_and_all_failure_addenda']=False
            old.write_text(json.dumps(altered,indent=2)+'\n'+(' ' if case=='seal_whitespace_before_first_hash_capture' else ''))
        if case!='success':spies['mutation_applied']=True
    error=None
    # Fake process/runtime identity and logical ownership context: no actual curl or flock.
    with patch('subprocess.Popen',fake_popen),patch('subprocess.run',lambda *a,**k:types.SimpleNamespace(stdout=p['curl_version'])),patch('os.getpgid',lambda pid:pid),patch.object(terminal_module,'save_new',terminal_save),contextlib.redirect_stdout(io.StringIO()):
        try:ns['execute'](planhash)
        except BaseException as exc:error=type(exc).__name__+': '+str(exc)
    pending=Path(p['pending_path']);sealpath=Path(p['terminal_seal_path'])
    result={'case':case,'synthetic_total_body_bytes':scaled,'byte_budget_constant_only_AST_scaling_count':transform.changed,
        'pending_exists':pending.exists(),'seal_exists':sealpath.exists(),'seal_failure_exists':Path(str(sealpath)+'.failure.json').exists(),
        'error':error,'terminal_success':not pending.exists() and sealpath.is_file() and not Path(str(sealpath)+'.failure.json').exists(),
        'family_primary_intended_sha256':intended,'family_primary_actual_sha256':{str(q):sha(q) for q in family_paths if q.is_file()},
        'family_primary_statuses':[json.loads(q.read_text())['status'] for q in family_paths if q.is_file()],**spies}
    result['private_seal_intended_before_write']=seal_intended
    result['private_seal_actual_sha256']=sha(sealpath) if sealpath.is_file() else None
    consumer_error=None;consumer_sha=None
    if sealpath.is_file():
        binding=json.loads(admission.read_text())
        expected_binding={'plan_sha256':planhash,'admission_sha256':sha(admission),'executor_sha256':sha(SCRIPT)}
        try:consumer_sha=terminal_module.require_committed(pending,sealpath,expected_binding,intended)
        except BaseException as exc:consumer_error=type(exc).__name__+': '+str(exc)
    result['read_only_consumer_sha256']=consumer_sha;result['read_only_consumer_error']=consumer_error
    if case=='success':
        assert consumer_sha==seal_intended['sha256']==result['private_seal_actual_sha256']
    else:assert pending.exists() and consumer_sha is None
    save(base/'control_receipt.json',result); return result

def main():
    assert not OUT.exists() and not RECEIPT.exists(), 'Preserve previous independent review; no rerun/overwrite'
    plan,proof=metadata()
    prepare_controls=prepare_reader_controls(plan)
    cases=['success','mutate_both_before_hash_capture','mutate_one_before_hash_capture','samebytes_repo_logs_parent_after_seal','samebytes_primary_leaf_after_seal','source_body_after_seal','source_mirror_after_seal','source_failure_after_seal','family_failure_after_seal','deferred_signal_after_seal','resource_failure_after_seal','wrong_transport_header','cleanup_and_receipt_persistence_failure','missing_primary_before_hash_capture','primary_corruption_after_seal','missing_source_body_after_seal','missing_primary_after_seal','pending_leaf_symlink_after_seal','seal_leaf_symlink_after_seal','pending_corruption_after_seal','seal_content_before_first_hash_capture','seal_binding_before_first_hash_capture','seal_primary_map_before_first_hash_capture','seal_qualification_before_first_hash_capture','seal_whitespace_before_first_hash_capture']
    controls=[control(plan,c) for c in cases]
    vulnerable=set()
    for r in controls:
        assert r['terminal_success']==(r['case']=='success' or r['case'] in vulnerable), r
        if r['case'] in vulnerable:assert r['mutation_applied']
        if r['case']=='cleanup_and_receipt_persistence_failure':assert r['quarantine_spy_calls']>=1
    _,after_proof=metadata()
    assert proof['dependency_sha256']==after_proof['dependency_sha256']
    assert sha(PLAN)==EXPECTED_PLAN
    assert sha(SCRIPT)==plan['bound_sources'][str(SCRIPT)]
    result={'schema':'independent_original13_acquisition_metadata_and_tiny_controls_v3','recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'verdict':'QUALIFIED_PASS_90_V3_PRELAUNCH_TERMINAL_DELTA_ONLY','plan_sha256':EXPECTED_PLAN,'executor_sha256':sha(SCRIPT),'checker_sha256':sha(__file__),
        'metadata_proof':proof,'actual_prepare_reader_controls':prepare_controls,'dependencies_fresh_rechecked_after_controls':True,'controls':controls,'measured_gap_cases':sorted(vulnerable),'GWAS_or_reference_body_bytes_read':0,'actual_network_transfers':0,'actual_flock_acquisitions':0,
        'synthetic_fixture_curl_processes_only':True,'production_code_or_plan_mutated':False,'scientific_estimates_or_replications':0}
    save(RECEIPT,result)
    print(json.dumps({'verdict':result['verdict'],'metadata_dependencies':proof['dependencies_checked'],'controls':len(controls),'measured_gaps':sorted(vulnerable),'receipt':str(RECEIPT)},indent=2))

if __name__=='__main__':main()
