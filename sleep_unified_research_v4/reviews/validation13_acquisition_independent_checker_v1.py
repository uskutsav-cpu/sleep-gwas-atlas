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
SCRIPT = PKG / 'scripts/81_acquire_validation_raw_sources_v1.py'
PLAN = PKG / 'manifests/validation_raw_acquisition_plan_v4_1.json'
DESIGN = PKG / 'statistical_validation/validation_raw_replay_design_v1.json'
SEAL = PKG / 'statistical_validation/validation_raw_replay_design_seal_v1.json'
OUT = Path(__file__).with_name('validation13_acquisition_independent_controls_v1')
RECEIPT = Path(__file__).with_name('validation13_acquisition_independent_receipt_v1.json')
EXPECTED_PLAN = 'fd753352967f01b2cfda3747044158540ce29103ee7725071789efb4855965e5'

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
        assert m['body_path']==d['proposed_body_path']
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
    return plan,{'dependencies_checked':len(observed),'design_metadata_dependencies':53,'dependency_sha256':observed,'Finn_generation_and_MD5_bindings':finn,'MVP_original_body_and_official_yaml_MD5_bindings':mvp,'immutable_helper_AST_matches':immutable}

class ScaleSyntheticByteBudget(ast.NodeTransformer):
    def __init__(self,total):self.total=total;self.changed=0
    def visit_Constant(self,node):
        if node.value==10058648185:
            self.changed+=1; return ast.copy_location(ast.Constant(self.total),node)
        return node

def control(plan,case):
    base=OUT/case; base.mkdir(parents=True)
    pkg=base/'PKG'; ssd=base/'SSD'; folder=ssd/'validation_raw_replay_v1'
    for p in [pkg/'logs',pkg/'source_provenance/validation_raw_acquisition_v4_1',folder/'raw',folder/'logs',folder/'receipts',ssd/'tmp']:
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
    original_write=ns['write_new']; family_paths=[pkg/'logs/validation_raw_acquisition_family_receipt_v4_1.json',folder/'validation_raw_acquisition_family_receipt_v4_1.json']
    intended={}
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
    ns['write_new']=write_hook
    if case=='cleanup_and_receipt_persistence_failure':
        ns['cleanup_proof']=lambda *a:{'teardown_verified':False,'remaining_group_members':[19000001],'cleanup_error':'SYNTHETIC'}
        ns['retain_lock_until_gone']=lambda *a:spies.__setitem__('quarantine_spy_calls',spies['quarantine_spy_calls']+1)
    terminal_module=sys.modules['terminal_commit_common_v2']; original_save=terminal_module.save_new
    def terminal_save(path,value):
        original_save(path,value)
        if Path(path)!=Path(p['terminal_seal_path']):return
        if case=='samebytes_repo_logs_parent_after_seal':
            old=pkg/'logs'; retained=pkg/'logs_retained';old.rename(retained);old.symlink_to(retained, target_is_directory=True)
        elif case=='samebytes_primary_leaf_after_seal':
            old=family_paths[0];retained=old.with_suffix('.retained.json');old.rename(retained);old.symlink_to(retained)
        elif case=='source_body_after_seal':Path(p['members'][0]['body_path']).write_bytes(b'ALTERED_SYNTHETIC')
        elif case=='source_mirror_after_seal':
            q=pkg/'source_provenance/validation_raw_acquisition_v4_1'/(p['members'][0]['source_id']+'.json');q.write_text(q.read_text()+' ')
        elif case=='source_failure_after_seal':
            q=folder/'receipts'/(p['members'][0]['source_id']+'.json.failure.json');save(q,{'synthetic_failure':True})
        elif case=='family_failure_after_seal':save(Path(str(family_paths[0])+'.failure.json'),{'synthetic_failure':True})
        elif case=='deferred_signal_after_seal':ns['TERMINATION_REQUEST'].append(15)
        elif case=='resource_failure_after_seal':state['internal_free_bytes']=0
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
    save(base/'control_receipt.json',result); return result

def main():
    assert not OUT.exists() and not RECEIPT.exists(), 'Preserve previous independent review; no rerun/overwrite'
    plan,proof=metadata()
    cases=['success','mutate_both_before_hash_capture','mutate_one_before_hash_capture','samebytes_repo_logs_parent_after_seal','samebytes_primary_leaf_after_seal','source_body_after_seal','source_mirror_after_seal','source_failure_after_seal','family_failure_after_seal','deferred_signal_after_seal','resource_failure_after_seal','wrong_transport_header','cleanup_and_receipt_persistence_failure']
    controls=[control(plan,c) for c in cases]
    vulnerable={'mutate_both_before_hash_capture','samebytes_repo_logs_parent_after_seal'}
    for r in controls:
        assert r['terminal_success']==(r['case']=='success' or r['case'] in vulnerable), r
        if r['case'] in vulnerable:assert r['mutation_applied']
        if r['case']=='cleanup_and_receipt_persistence_failure':assert r['quarantine_spy_calls']>=1
    assert sha(PLAN)==EXPECTED_PLAN
    assert sha(SCRIPT)==plan['bound_sources'][str(SCRIPT)]
    result={'schema':'independent_original13_acquisition_metadata_and_tiny_controls_v1','recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'verdict':'MEASURED_TERMINAL_GAPS_DO_NOT_LAUNCH_81_V1','plan_sha256':EXPECTED_PLAN,'executor_sha256':sha(SCRIPT),'checker_sha256':sha(__file__),
        'metadata_proof':proof,'controls':controls,'measured_gap_cases':sorted(vulnerable),'GWAS_or_reference_body_bytes_read':0,'actual_network_transfers':0,'actual_flock_acquisitions':0,
        'synthetic_fixture_curl_processes_only':True,'production_code_or_plan_mutated':False,'scientific_estimates_or_replications':0}
    save(RECEIPT,result)
    print(json.dumps({'verdict':result['verdict'],'metadata_dependencies':proof['dependencies_checked'],'controls':len(controls),'measured_gaps':sorted(vulnerable),'receipt':str(RECEIPT)},indent=2))

if __name__=='__main__':main()
