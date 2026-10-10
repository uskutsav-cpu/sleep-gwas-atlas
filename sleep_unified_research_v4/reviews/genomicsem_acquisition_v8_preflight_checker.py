"""Independent v8 successor bindings, no production GWAS reads or transfers."""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import tempfile
import time
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN=P/'manifests/extension_raw_acquisition_plan_v4_8.json';OLD=P/'manifests/extension_raw_acquisition_plan_v4_7.json'
EXECUTOR=S/'52_acquire_extension_raw_sources_v8.py'
PINS={str(PLAN):'77bf2da98b851ac1aa9d414b2b695fb04a65b0c4ca987f7d85db2097e9e7ef8b',str(EXECUTOR):'5d4fd07d1e91458fdfb818b501d738a37c0e7575a07871f9e992d525d31aedd4',str(OLD):'0958a8383ada38c25eceb6451a7d074527753fae5dffce278d02e8bb2ba1345c',str(P/'manifests/global_SSD_resource_reservation_v4_4.json'):'04e2ea61d9a06281468674c1a8326367ffb7e507857814365e63108ee38ebbc7',str(R/'genomicsem_acquisition_v7_preflight_seal_v2.json'):'05c890036b61660314fdbff5c9dce9ff59385882c4305b82de7ae2575936ceb0'}
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
checks=[]
def check(value,label,detail=None):
    checks.append({'control':label,'pass':bool(value),'detail':detail});assert value,label
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
for path,h in PINS.items():check(sha(path)==h,'pin_'+Path(path).name)
for sealpath in [R/'genomicsem_acquisition_v7_preflight_seal_v2.json',R/'genomicsem_acquisition_v7_preflight_seal.json',R/'genomicsem_extension_pipeline_v5_preflight_seal.json']:
    seal=json.loads(sealpath.read_text())
    for path,identity in seal['artifacts'].items():check(sha(path)==identity['sha256'],'inherited_sealed_'+Path(path).name)
plan,old=[json.loads(x.read_text()) for x in [PLAN,OLD]]
for x in [SSD/'manifests'/PLAN.name]:check(sha(x)==PINS[str(PLAN)],'exact_SSD_plan_copy')
forbidden={m['body_path'] for q in [plan,old] for m in q['members']}
forbidden|={x+'.partial' for x in forbidden};forbidden|={q['explicit_resume']['original_partial_path'] for q in [plan,old]}
before={}
for path,h in plan['bound_sources'].items():
    check(path not in forbidden and not Path(path).name.endswith(('.gz','.bgz','.bed','.bim','.fam')),'non_body_binding_'+Path(path).name)
    check(Path(path).is_file() and not Path(path).is_symlink() and sha(path)==h,'bound_metadata_'+Path(path).name)
    before[path]={'sha256':h,'bytes':Path(path).stat().st_size}
check(len(plan['members'])==len(old['members'])==100,'exact100_ordered_members')
for a,b in zip(old['members'],plan['members']):
    expected=dict(a)
    if a['index']>4:expected['body_path']=str(SSD/'extension_raw_replay_v8/raw'/a['filename'])
    check(expected==b,'original_source_tuple_and_private_path_'+str(b['index']))
check(set(plan['reused_checkpoints'])=={'1','2','3','4'} and plan['network_sources_to_transfer']==96,'exact4_reused96_new')
check(all(plan['bound_sources'].get(k)==v for k,v in old['bound_sources'].items()),'all152_prior_bound_declarations_retained')
for k in ['compressed_network_bytes','SSD_reservation_bytes','largest_body_bytes','internal_floor_bytes','ssd_floor_bytes','maximum_owned_transfer_rss_bytes','per_body_seconds_limit','family_seconds_limit','transfer_worker_count','automatic_retry','retained_raw_bodies','protected_GWAS_bodies_committed','scientific_membership_or_threshold_changes','scope','curl_path','curl_version','monitor_path','monitor_output_limit_bytes','runtime_poll_seconds','exclusive_family_lock_path','terminal_protocol','all_failed_provisional_primary_copies_must_be_rejected']:
    check(plan[k]==old[k],'unchanged_guard_or_science_'+k)
check(plan['internal_floor_bytes']==3*1024**3 and plan['ssd_floor_bytes']==5*1024**3 and plan['maximum_owned_transfer_rss_bytes']==2*1024**3 and plan['SSD_reservation_bytes']==300*1024**3,'unchanged3_5floors_2GiBRSS_300GiBcap')
check(plan['per_body_seconds_limit']==7200 and plan['family_seconds_limit']==345600 and plan['automatic_retry'] is False and plan['transfer_worker_count']==1,'unchanged2h_96h_one_curl_manual_only')
r=plan['explicit_resume'];member=plan['members'][4]
check(r['source_index']==5 and r['source_trait']==member['extension_trait_id'],'exact_source5_trait')
check(r['prefix_bytes']==559470675 and r['prefix_sha256']=='a578bc085bdde5ade891948632ab195c635a1a0bd9fd00545eedf30cbfef4802','exact_source5_preserved_prefix_binding')
check(member['expected_bytes']==2299805399 and member['expected_sha256']=='3336d342199b26d0b60fb7a4ea3a39c63e2184b64c7fe09bde9816eff20e8136','same_original_full_source5_identity')
check(r['original_partial_path']==old['members'][4]['body_path']+'.partial','old_source5_prefix_namespace_preserved')
prefix=Path(r['original_partial_path']);check(not prefix.is_symlink() and prefix.stat().st_size==r['prefix_bytes'],'old_source5_prefix_size_stat_only_no_body_hash')
check(plan['new_network_bytes']==sum(m['expected_bytes'] for m in plan['members'][4:])-r['prefix_bytes']==217937033696,'exact_remaining_network_byte_arithmetic')
check(plan['additional_preserved_prefix_bytes']==old['additional_preserved_prefix_bytes']+r['prefix_bytes']==4091803839,'all_preserved_prefix_byte_arithmetic')
failedpath=P/'logs/extension_raw_acquisition_family_receipt_v4_7.json';failed=json.loads(failedpath.read_text());sourcefail=json.loads(Path(r['prior_failed_receipt']).read_text())
check(sha(failedpath)==plan['prior_v7_failed_family_sha256'] and failed['status']=='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT' and failed['completed_source_count']==4,'v7_failed_family_preserved_not_reclassified')
check(sourcefail['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and sourcefail['returncode']==56 and sourcefail['stop_reason']=='AssertionError: CURL_NONZERO_EXIT: 56','failed_source5_curl56_preserved')
check(sourcefail['teardown']['teardown_verified'] and sourcefail['teardown']['remaining_group_members']==[],'failed_source5_recorded_cleanup_empty')
check(sourcefail['retained_partial_sha256']==r['prefix_sha256'] and sourcefail['retained_partial_bytes']==r['prefix_bytes'],'prefix_bound_to_actual_failed_source5_receipt')
check(Path(old['pending_path']).is_file() and not Path(old['terminal_seal_path']).exists(),'v7_PENDING_preserved_no_family_terminal_credit')
check(plan['bound_sources'][str(P/'manifests/extension_raw_acquisition_admission_v4_7.json')]==sha(P/'manifests/extension_raw_acquisition_admission_v4_7.json'),'prior_v7_root_admission_exact_binding')
successful={x['path']:x['sha256'] for x in failed['source_receipts']};gateplan=json.loads(Path(plan['reused_checkpoint_gate_plan']).read_text());common=module(S/'extension_replay_common_v4.py','_v8_individual_origins')
accepted=[]
for m in plan['members'][:4]:
    fixed=plan['reused_checkpoints'][str(m['index'])];receipt=json.loads(Path(fixed['receipt_path']).read_text())
    check(sha(fixed['receipt_path'])==fixed['receipt_sha256']==successful[fixed['receipt_path']],'four_reused_exact_receipt_and_prior_failure_map_'+str(m['index']))
    check(receipt['member']==m and receipt['plan_sha256']==fixed['original_acquisition_plan_sha256'],'four_exact_members_and_origin_'+str(m['index']))
    prior=old['reused_checkpoints'].get(str(m['index']))
    check(fixed==prior if prior else fixed['original_acquisition_plan_sha256']==PINS[str(OLD)],'unchanged_three_old_origins_and_one_exact_v7_origin_'+str(m['index']))
    gm=next(x for x in gateplan['members'] if x['index']==m['index'])
    check(common.acquisition_receipt_gate(gateplan,gm,fixed['receipt_sha256'])==receipt,'actual_receipt_header_stat_gate_no_body_hash_'+str(m['index']))
    accepted.append({'index':m['index'],'receipt_sha256':fixed['receipt_sha256'],'origin_plan_sha256':fixed['original_acquisition_plan_sha256'],'body_whole_hash_independently_reread':False})
ledger=json.loads(Path(plan['global_resource_ledger_path']).read_text());priorledger=json.loads(Path(old['global_resource_ledger_path']).read_text())
check(sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes']==320204645117,'all_current_ledger_components_sum')
check(ledger['reserved_total_bytes']<ledger['ceiling_bytes']==322122547200 and ledger['unallocated_margin_bytes']==1917902083,'unchanged_ceiling_with_exact_margin')
check(ledger['reserved_total_bytes']-priorledger['reserved_total_bytes']==r['prefix_bytes'],'ledger_increase_only_preserved_source5_prefix')
check({k:v for k,v in ledger['component_bytes'].items() if k!='preserved_prefix_bytes'}=={k:v for k,v in priorledger['component_bytes'].items() if k!='preserved_prefix_bytes'},'all_other_component_caps_unchanged')
check(ledger['component_bytes']['preserved_prefix_bytes']==plan['additional_preserved_prefix_bytes'] and ledger['source5_failed_receipt_sha256']==sha(r['prior_failed_receipt']),'prefix_and_failed_receipt_frozen_in_global_ledger')

functions=lambda path:{n.name:n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef)}
a,b=functions(S/'52_acquire_extension_raw_sources_v7.py'),functions(EXECUTOR)
changed={'prepare','reused_checkpoint_gate','acquire','finalize_family','execute'}
for name in a:
    if name not in changed:check(ast.dump(a[name],include_attributes=False)==ast.dump(b[name],include_attributes=False),'exact_unchanged_method_AST_'+name)
class EpochOnly(ast.NodeTransformer):
    def visit_Constant(self,n):
        if isinstance(n.value,str):n.value=n.value.replace('extension_raw_acquisition_v4_8','extension_raw_acquisition_v4_7').replace('extension_raw_acquisition_family_receipt_v4_8','extension_raw_acquisition_family_receipt_v4_7')
        return n
for name in ['acquire','finalize_family']:
    check(ast.dump(a[name],include_attributes=False)==ast.dump(EpochOnly().visit(copy.deepcopy(b[name])),include_attributes=False),name+'_AST_only_epoch_paths')
class OwnershipOnly(ast.NodeTransformer):
    def visit_Dict(self,n):
        for k,v in zip(n.keys,n.values):
            if isinstance(k,ast.Constant) and k.value=='immutable_attempt' and isinstance(v,ast.Constant):
                check(v.value==8,'corrected_ownership_attempt8');v.value=6
        return self.generic_visit(n)
check(ast.dump(a['execute'],include_attributes=False)==ast.dump(OwnershipOnly().visit(copy.deepcopy(b['execute'])),include_attributes=False),'execute_AST_only_ownership_epoch_literal')
class OriginImportOnly(ast.NodeTransformer):
    def visit_Constant(self,n):
        if isinstance(n.value,str):n.value=n.value.replace('extension_replay_common_v4.py','extension_replay_common_v3.py')
        return n
check(ast.dump(a['reused_checkpoint_gate'],include_attributes=False)==ast.dump(OriginImportOnly().visit(copy.deepcopy(b['reused_checkpoint_gate'])),include_attributes=False),'origin_gate_AST_only_reviewed_common4_import')

mod=module(EXECUTOR,'_v8_bounded_meta_fixtures');observations=[]
with tempfile.TemporaryDirectory(prefix='raw-v8-meta-review-') as td:
    folder=Path(td)
    # Execute exact source HTTP assertion AST with actual source5 constants and
    # invented header bytes. No curl, body or acquire() invocation occurs.
    try_node=next(n for n in b['acquire'].body if isinstance(n,ast.Try))
    header_nodes=[n for n in try_node.body if (isinstance(n,ast.Assert) and 'observed' in ast.unparse(n.test)) or (isinstance(n,ast.If) and ast.unparse(n.test)=='offset')]
    code=compile(ast.fix_missing_locations(ast.Module(body=header_nodes,type_ignores=[])),str(EXECUTOR),'exec')
    base={'status':206,'x-amz-version-id':member['expected_s3_version_id'],'etag':member['expected_etag'],'content-length':str(member['expected_bytes']-r['prefix_bytes']),'content-range':f"bytes {r['prefix_bytes']}-{member['expected_bytes']-1}/{member['expected_bytes']}"}
    check(base['content-range']=='bytes 559470675-2299805398/2299805399' and base['content-length']=='1740334724','exact_source5_HTTP206_range_and_suffix_length')
    h=folder/'headers.txt';h.write_text('HTTP/2 206\n'+''.join(k+': '+v+'\n' for k,v in base.items() if k!='status'))
    check(mod.header_fields(h)==base,'real_tiny_header_parser_exact206')
    exec(code,{'observed':base,'offset':r['prefix_bytes'],'member':member});check(True,'exact_source5_HTTP206_assertions_accept')
    for k,v in [('status',200),('x-amz-version-id','wrong'),('etag','wrong'),('content-length','2299805399'),('content-range','bytes 559470676-2299805398/2299805399')]:
        bad={**base,k:v}
        try:exec(code,{'observed':bad,'offset':r['prefix_bytes'],'member':member})
        except AssertionError:check(True,'source5_HTTP206_wrong_'+k+'_rejects')
        else:raise AssertionError('source5 malformed HTTP accepted '+k)
    for good in [True,False]:
        home=folder/('good' if good else 'wrong_counts');home.mkdir();package=home/'package';(package/'logs').mkdir(parents=True)
        pending=home/'pending.json';terminal=home/'terminal.json';mod.write_new(pending,{'only_fixture':True});items=[]
        for j in range(100):
            path=home/('source_%03d.json'%j);mod.write_new(path,{'synthetic_source':j});items.append({'path':str(path),'sha256':sha(path)})
        fixture={'reused_checkpoints':{str(j):{} for j in range(1,5)},'terminal_protocol':'PRIVATE_SYNTHETIC_PROTOCOL','terminal_seal_path':str(terminal),'pending_path':str(pending),'internal_floor_bytes':3*1024**3,'ssd_floor_bytes':5*1024**3,'family_seconds_limit':1000}
        master={'status':'ALL100_EXACT_SOURCE_BODIES_ACQUIRED','source_receipts':items,'reused_exact_source_count':4 if good else 3,'new_exact_source_count':96 if good else 97}
        monitor=SimpleNamespace(snapshot=lambda:{'internal_free_bytes':4*1024**3,'ssd_free_bytes':6*1024**3})
        with patch.object(mod,'PACKAGE',package),patch.object(mod,'FOLDER',home),patch.object(mod,'assert_bindings',lambda *_:None),patch.object(mod,'physical_mount',lambda:None),patch.object(mod,'global_namespace_gate',lambda *_:0),patch('subprocess.Popen',side_effect=AssertionError('NO_TRANSFER_PERMITTED')):
            committed=mod.finalize_family(fixture,'SYNTHETIC_PLAN',master,monitor,time.monotonic(),pending,sha(pending))
        copies=[package/'logs/extension_raw_acquisition_family_receipt_v4_8.json',home/'extension_raw_acquisition_family_receipt_v4_8.json']
        if good:
            seal=json.loads(terminal.read_text());check(committed and not pending.exists() and seal['primary_receipt_sha256']=={str(p):sha(p) for p in copies} and len({sha(p) for p in copies})==1,'healthy4_reused96_new_exact_two_copies')
        else:check(not committed and pending.exists() and all(Path(str(p)+'.failure.json').exists() for p in [*copies,terminal]),'old3_reused97_new_rejects_pending_and_failure_addenda')
        observations.append({'healthy4_96':good,'committed':committed,'pending_retained':pending.exists(),'real_workers':0,'identity_mount_namespace_callbacks':'EXPLICIT_FIXTURE_MOCKS'})
for path,meta in before.items():check(sha(path)==meta['sha256'] and Path(path).stat().st_size==meta['bytes'],'bound_metadata_unchanged_after_'+Path(path).name)
for path,h in PINS.items():check(sha(path)==h,'pin_unchanged_after_'+Path(path).name)
result={'status':'PASS_NARROW_V8_PRELAUNCH_METADATA_ORIGIN_RESUME_AND_CARDINALITY','checks':checks,'check_count':len(checks),'plan_sha256':PINS[str(PLAN)],'executor_sha256':PINS[str(EXECUTOR)],'fresh_metadata_before_after':before,'fresh_metadata_count':len(before),'metadata_bytes_per_pass':sum(x['bytes'] for x in before.values()),'four_preserved_individual_origin_receipt_header_stat_proofs':accepted,'actual_prefix_stat_only':True,'actual_body_prefix_SHA_are_root_preparation_receipt_claims_not_review_body_rereads':True,'new_count_controls':observations,'inherited634_v7_controls_and_v6_12_terminal8_cleanup_controls':'REUSED_VIA_EXACT_UNCHANGED_METHOD_AST_NOT_REPEATED','no_live_HEAD_needed':'Original immutable version source tuple/header is frozen; new actual206range remains execution gate, only tiny assertions exercised here.','v7_failed_family_preserved':True,'execution_admitted':False,'actual_full_family_acquired':False,'full_gzip_EOF_or_pipeline_certified':False,'production_raw_body_reads':0,'reference_body_reads':0,'real_transfer_or_fit_operations':0,'production_mutex_operations':0,'peak_reviewer_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
with (R/'genomicsem_acquisition_v8_controls_receipt.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':result['status'],'checks':len(checks),'fresh_metadata':len(before),'bytes_per_pass':result['metadata_bytes_per_pass'],'peak_RSS':result['peak_reviewer_RSS_bytes']}))
