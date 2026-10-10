#!/usr/bin/env python3
"""Bounded frozen v4 blueprint and metadata fixtures, no production body IO."""
import ast
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import tempfile
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN=SSD/'extension_pipeline_replay_v4/extension_pipeline_replay_plan_v3.json'
OLD=SSD/'extension_pipeline_replay_v3/extension_pipeline_replay_plan_v3.json'
RAW=P/'manifests/extension_raw_acquisition_plan_v4_7.json'
PINNED={str(PLAN):'7c941a230287e04b42526b215b1f13373b850885fd9c8abacd2f0ab83731c690',
 str(OLD):'565b34e997b2cab841f3a101cb4f8649c4c4660e56f5b1421ce4f68f29afe4ac',
 str(RAW):'0958a8383ada38c25eceb6451a7d074527753fae5dffce278d02e8bb2ba1345c',
 str(S/'46_prepare_extension_pipeline_replay_v4.py'):'9de2ccafe1b51076f89482bacb8d2a2ed9f46f66ec0f0eb2b06d4e2679a4f6ba',
 str(S/'47_run_extension_pipeline_replay_v4.py'):'016e35eca161740cbba3c811f5a62574543260f9bea65ed58c1e28a104c8e31c',
 str(S/'extension_replay_common_v4.py'):'319fd5da9dbac1ef197aeb107c566394d8a43bff25241af0c63e6053ef19c2e7',
 str(S/'extension_replay_validate_v4.py'):'ec6033bc749eca5fa7d0911d0a6ad9de6690fe256be7f2f97beef4dfb5496500',
 str(S/'sensitivity_executor_v4_4.py'):'7ac802596f29adfacb017f406676e6a30ee00bf890904976afda7f88a3ffa9c5',
 str(S/'terminal_commit_common_v2.py'):'9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd',
 str(R/'genomicsem_extension_checkpoint_review_seal_v3.json'):'adb876ce4f45a85178b7e446ed46b87b1577655d32898dfbdece7c0c2532839a'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(path,name):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m


class PathOnly(ast.NodeTransformer):
 def visit_Constant(self,n):
  if isinstance(n.value,str):
   for a,b in [('extension_pipeline_replay_v4','extension_pipeline_replay_v3'),('pipeline_attempt_v4.json','pipeline_attempt_v3.json'),('pipeline_pending_v4.json','pipeline_pending_v3.json'),('pipeline_terminal_seal_v4.json','pipeline_terminal_seal_v3.json'),('pipeline_execution_receipt_v4.json','pipeline_execution_receipt_v3.json'),('extension_replay_validate_v4.py','extension_replay_validate_v3.py'),('extension_pipeline_replay_plan_v4.json','extension_pipeline_replay_plan_v2.json')]:n.value=n.value.replace(a,b)
  return n
 def visit_ImportFrom(self,n):
  if n.module=='extension_replay_common_v4':n.module='extension_replay_common_v3'
  return n


def run():
 checks=[]
 def check(value,label):
  if not value:raise AssertionError(label)
  checks.append({'case':label,'pass':True})
 for p,h in PINNED.items():check(sha(p)==h,'pin '+p)
 old,plan,raw=[json.loads(p.read_text()) for p in [OLD,PLAN,RAW]]
 oldseal=json.loads((R/'genomicsem_extension_checkpoint_review_seal_v3.json').read_text())
 for p,m in oldseal['artifacts'].items():check(sha(p)==m['sha256'],'sealed_v3_artifact '+p)
 check(len(old['members'])==len(plan['members'])==len(raw['members'])==100,'exact100_ordered_members')
 for a,b,r in zip(old['members'],plan['members'],raw['members']):
  expected=copy.deepcopy(a)
  for k,v in list(expected.items()):
   if isinstance(v,str) and v.startswith(str(OLD.parent)+'/'):expected[k]=v.replace(str(OLD.parent)+'/',str(PLAN.parent)+'/',1)
   elif k in ['harmonize_command','munge_command']:expected[k]=[x.replace(str(OLD.parent)+'/',str(PLAN.parent)+'/',1) if x.startswith(str(OLD.parent)+'/') else r['body_path'] if x==a['raw'] else x for x in v]
  expected['raw'],expected['acquisition_member']=r['body_path'],r
  reuse=raw['reused_checkpoints'].get(str(r['index']))
  expected['acquisition_receipt']=reuse['receipt_path'] if reuse else str(Path(raw['pending_path']).parent/'receipts'/(r['extension_trait_id']+'.json'))
  expected['acquisition_origin_plan_sha256']=reuse['original_acquisition_plan_sha256'] if reuse else PINNED[str(RAW)]
  check(expected==b,'path_origin_only_member_and_science_argv '+b['extension_trait_id'])
  check({k:v for k,v in a['acquisition_member'].items() if k!='body_path'}=={k:v for k,v in r.items() if k!='body_path'},'original_source_tuple '+b['extension_trait_id'])
 for k in ['reference','reference_sha256','w_hm3','template_rows','environment','archived_input_sha256','original_190_jobs','baseline_execution_plan_sha256','baseline_input_sha256','baseline_dependency_sha256','baseline_relocated_dependency_sha256','guard','pipeline_terminal_protocol']:
  check(old[k]==plan[k],'unchanged_blueprint_field '+k)
 check(all(plan['dependencies_sha256'].get(p)==h for p,h in old['dependencies_sha256'].items()),'all_preserved_v3_dependency_declarations_retained')
 check(all(plan['dependencies_sha256'].get(p)==h for p,h in raw['bound_sources'].items()),'all_v7_bound_source_declarations_retained')
 for origin,identity in old['acquisition_operational_identity_by_origin'].items():check(plan['acquisition_operational_identity_by_origin'][origin]==identity,'old_origin_identity_retained '+origin)
 check(plan['acquisition_operational_identity_by_origin'][PINNED[str(RAW)]]=={'source_folder':str(Path(raw['pending_path']).parent),'curl':'/usr/bin/curl','per_body_seconds_limit':7200,'resume_offset_by_index':{'4':1126899820}},'exact_new_v7_source4_offset_origin')
 contract=plan['acquisition_family_terminal_contract'];primary=P/'logs/extension_raw_acquisition_family_receipt_v4_7.json'
 check(contract==dict(pending_path=raw['pending_path'],terminal_seal_path=raw['terminal_seal_path'],primary_receipt_paths=[str(primary),str(Path(raw['pending_path']).parent/primary.name)],reused_exact_source_count=3,new_exact_source_count=97,checkpoint_credit_while_pending='INDIVIDUAL_EXACT_SUCCESSFUL_SOURCE_ONLY',full_family_credit='BOTH_HASH_BOUND_PRIMARY_COPIES_AND_FULL100_TERMINAL_SEAL_AND_ABSENT_PENDING_AND_FAILURE_ADDENDA'),'exact_v7_full100_3plus97_contract')
 check(plan['global_resource_ledger_path']==raw['global_resource_ledger_path'] and plan['global_resource_ledger_sha256']==raw['global_resource_ledger_sha256'],'exact_v7_ledger_declaration')
 metadata=dict(PINNED);metadata.update(plan['acquisition_execution_binding']['sha256'])
 metadata.update({p:h for p,h in plan['dependencies_sha256'].items() if p not in old['dependencies_sha256']})
 for name in ['independent_core_numerical_adjudication_v4.sha256','independent_whole_extension_v4.sha256','independent_whole_validation_v4.sha256']:
  p=str(R/name);check(p in metadata and metadata[p]==plan['dependencies_sha256'][p],'full190_independent_review_seal '+name)
 forbidden=set(plan['archived_input_sha256'])|{m['raw'] for m in plan['members']}|{plan['reference'],plan['w_hm3']}
 for p,h in metadata.items():check(p not in forbidden and Path(p).suffix not in ['.gz','.bgz','.bed','.bim','.fam'] and sha(p)==h,'static_non_body_metadata '+p)
 a=ast.parse((S/'47_run_extension_pipeline_replay_v3.py').read_text());b=PathOnly().visit(ast.parse((S/'47_run_extension_pipeline_replay_v4.py').read_text()))
 check(ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False),'entire_runner_AST_only_version_path_import_and_default_deltas')
 a=ast.parse((S/'extension_replay_validate_v3.py').read_text());b=PathOnly().visit(ast.parse((S/'extension_replay_validate_v4.py').read_text()))
 check(ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False),'entire_validator_AST_import_delta_only')
 funcs=lambda p:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
 a,b=funcs(S/'extension_replay_common_v3.py'),funcs(S/'extension_replay_common_v4.py')
 for name in a:
  if name!='acquisition_family_gate':check(a[name]==b[name],'exact_common_science_origin_function '+name)
 sys.path.insert(0,str(S));common=module(S/'extension_replay_common_v4.py','_pipeline4_common');runner=module(S/'47_run_extension_pipeline_replay_v4.py','_pipeline4_runner')
 accepted=[]
 for member in plan['members'][:3]:
  fixed=raw['reused_checkpoints'][str(member['index'])]
  r=common.acquisition_receipt_gate(plan,member,fixed['receipt_sha256']);check(r['member']==member['acquisition_member'] and r['plan_sha256']==member['acquisition_origin_plan_sha256'],'actual_old3_receipt_header_stat_origin '+str(member['index']))
  accepted.append({'index':member['index'],'receipt_sha256':fixed['receipt_sha256'],'origin':r['plan_sha256'],'body_reread':False})
 common.acquisition_operational_gate(plan);check(True,'exact_v7_root_admission_and_review_metadata_gate')
 fixture=module(P/'statistical_validation/genomicsem_extension_checkpoint_identity_controls_v3.py','_family_fixture4')
 controls=[]
 with tempfile.TemporaryDirectory(prefix='pipeline4-review-') as td:
  def newfamily(label):
   f=fixture.family(Path(td)/label);f['plan']['acquisition_family_terminal_contract'].update(reused_exact_source_count=3,new_exact_source_count=97)
   f['primary'].update(reused_exact_source_count=3,new_exact_source_count=97)
   f['s']['primary_receipt_sha256']={str(p):fixture.save(p,f['primary']) for p in f['copies']};fixture.save(f['seal'],f['s']);return f
  def reject(f,code,label):
   try:common.acquisition_family_gate(f['plan'],f['map'],require_complete=True)
   except RuntimeError as e:check(code in str(e),label);controls.append({'case':label,'rejected':code})
   else:raise AssertionError('REJECTION_MISSING '+label)
  with patch.object(common,'acquisition_operational_gate',lambda *_:None):
   f=newfamily('healthy');check(common.acquisition_family_gate(f['plan'],f['map'],require_complete=True) is not None,'healthy_exact100_3plus97_family')
   f=newfamily('wrong_contract');f['plan']['acquisition_family_terminal_contract'].update(reused_exact_source_count=2,new_exact_source_count=98);reject(f,'SOURCE_FAMILY_FROZEN3_PLUS97_COUNTS_REQUIRED','old2plus98_contract_veto')
   f=newfamily('wrong_primary');f['primary'].update(reused_exact_source_count=2,new_exact_source_count=98);f['s']['primary_receipt_sha256']={str(p):fixture.save(p,f['primary']) for p in f['copies']};fixture.save(f['seal'],f['s']);reject(f,'SOURCE_FAMILY_PRIMARY_FULL100_BINDING_DIFFERS','old2plus98_primary_veto')
   wrong_primary=f
   f=newfamily('pending');fixture.save(f['pending'],{});check(common.acquisition_family_gate(f['plan'],require_complete=False) is None,'individual_checkpoint_credit_while_pending');reject(f,'SOURCE_FAMILY_TERMINAL_PENDING_OR_UNSEALED','pending_full_family_veto')
   f=newfamily('missing_one');f['s']['source_receipts'].pop();fixture.save(f['seal'],f['s']);reject(f,'SOURCE_FAMILY_FULL100_RECEIPT_SET_DIFFERS','missing_source_full100_map_veto')
   f=newfamily('changed_one');fixture.save(f['paths'][50],{'changed':True});reject(f,'SOURCE_FAMILY_FINAL_CHECKPOINT_CHANGED','changed_source_checkpoint_map_veto')
   f=newfamily('failed_family');fixture.save(f['copies'][0],dict(f['primary'],status='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'));reject(f,'SOURCE_ACQUISITION_FAMILY_FAILED_PRESERVED','failed_source_family_veto')
   f=newfamily('addendum');fixture.save(str(f['seal'])+'.failure.json',{});reject(f,'SOURCE_ACQUISITION_FAMILY_FAILURE_ADDENDUM','terminal_failure_addendum_veto')
   # Narrow integration of the new3/97 veto with the unchanged terminal path.
   protected=module(S/'sensitivity_executor_v4_4.py','_pipeline4_protected');protected.TERMINATION_REQUEST.clear()
   folder=Path(td)/'pipeline_commit';folder.mkdir();target=folder/'receipt.json';pending=folder/'pending.json';seal=folder/'seal.json'
   terminal=runner.TerminalCommit(pending,seal,{'fixture_only':True});master={'status':'ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS'}
   with patch.object(protected,'limits',lambda *_:({'fixture_only':True},None)):
    runner.seal_master(protected,master,target,{},0,[True,[]],terminal,lambda:common.acquisition_family_gate(wrong_primary['plan'],wrong_primary['map'],require_complete=True))
   check(master['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and pending.exists() and Path(str(target)+'.failure.json').exists(),'new3plus97_veto_after_pipeline_receipt_retains_pending')
 default=PLAN.parent/'extension_pipeline_replay_plan_v4.json';check(not default.exists() and PLAN.is_file(),'explicit_v3_filename_required_default_v4_unavailable')
 for p,h in metadata.items():check(sha(p)==h,'static_metadata_after '+p)
 return {'recorded_utc':datetime.now(timezone.utc).isoformat(),'status':'QUALIFIED_PIPELINE_V4_EXPLICIT_PLAN_PRELAUNCH_PASS','prepared_plan':str(PLAN),'prepared_plan_sha256':PINNED[str(PLAN)],'static_metadata_sha256_before_and_after':metadata,'static_metadata_file_count':len(metadata),'declared_dependencies':plan['dependencies_sha256'],'declarations_not_body_rereads':True,'check_count':len(checks),'checks':checks,'new_family_veto_controls':controls,'actual_three_receipt_origin_metadata':accepted,'prior304_controls_and_two_v2_rejection_witnesses_reused':True,'former_resource_and_late_family_paths_AST_unchanged':True,'helper_sha256':sha(__file__),'interface_qualification':'Mandatory explicit --plan to reviewed folder_v4/filename_v3 and exactSHA; omitted default_v4 is unavailable and fails before workers.','execution_admitted':False,'body_or_reference_rereads':0,'real_workers':0,'production_mutex_acquired':False,'fits':0,'actual_full100_family_claim':False,'live_data_or_worker_operations':False,'peak_review_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


if __name__=='__main__':
 out=R/'genomicsem_extension_pipeline_v4_controls_receipt.json'
 if out.exists():raise RuntimeError('DISTINCT_REVIEW_OUTPUT_REQUIRED')
 result=run()
 with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ['status','check_count','static_metadata_file_count','peak_review_RSS_bytes']}))
