#!/usr/bin/env python3
"""Independent exact metadata/AST binding review, not production admission.

64KiB hashing permitted generated preserved evidence; GWAS/reference bodies and
binary/runtime dependencies are recorded without fresh reads. Historical190
receipt bytes/commands/source maps are compared to the already admitted v6
map, without reopening estimates or giant source hashes.
"""
import ast
from collections import Counter
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
PLAN=P/'manifests/native_canonical_calibration_plan_v4_7.json'
PIN='a72501b3daf01fd4102a0f002fd78ca232ff2429b6fe84fb5903109503ab19b9'
FAILED=P/'logs/canonical_v4_6_failed_preserved_stop_inventory_v4_7.json'
FPIN='c80990a119349b0deda5b4d1c6ee15ceb1d16f3c79d0f9f457894b8998674215'
OUT=R/'independent_canonical_v4_7_binding_receipt.json'
START=time.monotonic();CHECKS=[];HASHED={};DEFERRED={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def check(x,label):
 CHECKS.append({'label':label,'passed':bool(x)})
 if not x:raise AssertionError(label)
def regular(p):
 p=Path(p);check(p.is_file() and not p.is_symlink() and not any(d.is_symlink() for d in p.parents),'regular:'+str(p));return p
def hashed(p,expected=None):
 p=regular(p);h=sha(p)
 if expected is not None:check(h==expected,'SHA:'+str(p))
 HASHED[str(p)]=h;return h
def load(p,expected=None):
 h=hashed(p,expected);r=json.loads(Path(p).read_text());check(sha(p)==h,'postparse:'+str(p));return r
def function(path,name):
 t=ast.parse(path.read_text());return ast.dump(next(n for n in t.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name==name),include_attributes=False)
def main():
 check(not OUT.exists(),'immutable_new_receipt');plan=load(PLAN,PIN);failed=load(FAILED,FPIN)
 old=load(Path(plan['preserved_v4_6_plan']),plan['preserved_v4_6_plan_sha256'])
 check(len(plan['jobs'])==len(old['jobs'])==16,'exact16_control_jobs')
 previous_root=old['ssd_output_root'];next_root=plan['ssd_output_root']
 converted=copy.deepcopy(old['jobs'])
 for j in converted:
  j['output_dir']=j['output_dir'].replace(previous_root,next_root);j['ldsc_args']=[x.replace(previous_root,next_root) for x in j['ldsc_args']]
 check(converted==plan['jobs'],'all16_exact_scientific_jobs_only_private_epoch_reroute')
 for field in ['scope','allowed_control_ids','block_receipt_path','block_receipt_sha256','build','interval_tsv_sha256','ordered_all_reference_records_sha256','ordered_construction_records_sha256','coordinate_map_path','coordinate_map_sha256','partition_amendment','static_construction_exclusions','reference_sha256','inputs_current_sha256','python','ldsc_dir','environment','worker_count','blas_threads','shared_lock_path','internal_floor_bytes','ssd_floor_bytes','worker_rss_limit_bytes','output_limit_bytes','per_worker_seconds_limit','arithmetic','cross_arm_identity_rule','global_resource_ledger_path','global_resource_ledger_sha256','global_reservation_bytes']:
  check(plan[field]==old[field],'unchanged_scientific_operational_field:'+field)
 check(plan['internal_floor_bytes']==3*2**30 and plan['ssd_floor_bytes']==5*2**30 and plan['worker_rss_limit_bytes']==2*2**30 and plan['output_limit_bytes']==8*2**30 and plan['global_reservation_bytes']==300*2**30 and plan['per_worker_seconds_limit']==7200,'exact_resource_caps')
 check(plan['all_retained_canonical_epochs_output_cap_bytes']==plan['output_limit_bytes'],'all_retained_epochs8GiB')
 check(plan['scope']=='METHOD_CALIBRATION_ONLY' and plan['preparation_only'] is True and plan['execution_admitted'] is False and plan['native_jobs_launched']==0 and plan['allow_41_covariance_outcomes'] is False and plan['allow_calibrated_biological_p_values'] is False and plan['realistic_LD_sampling_calibration_pass'] is False,'no_outcome_or_calibration_or_launch_admission')
 codes={
 '42_prepare_native_canonical_calibration_v4_7.py':'d93dab3e1e392003751b643ee58d3d092ab5ebf2f1927a5c14dc44a146cc1ca4',
 'canonical_calibration_common_v4_7.py':'bfa7612f11c258d7ffb9c4478b2a9d565716b273396cbc5fc8d663e34461d8af',
 'canonical_calibration_capture_v4_7.py':'9f85bed6e9d647967c31715fecdcbd1b0fdb4509de69bfe488736a92e4cea9f7',
 'canonical_calibration_arithmetic_v4_7.py':'a3148d8f5de1c7666f6ebeaa28b236152158a220c5d5cefec29a5d985d707ab5'}
 for n,h in codes.items():hashed(S/n,h);check(plan['dependencies_sha256'][str(S/n)]==h,'plan_binds_current_code:'+n)
 cap6=S/'canonical_calibration_capture_v4_6.py';cap7=S/'canonical_calibration_capture_v4_7.py'
 check(cap6.read_text().replace('v4_6','v4_7')==cap7.read_text(),'capture_only_helper_import_namespace_changed')
 check(function(S/'canonical_calibration_arithmetic_v4_6.py','validate')==function(S/'canonical_calibration_arithmetic_v4_7.py','validate'),'science_arithmetic_validate_exact_AST')
 for n in ['monitor_helpers','snapshot','owned_group_exists','await_owned_cleanup']:
  check(function(S/'42_prepare_native_canonical_calibration_v4_6.py',n)==function(S/'42_prepare_native_canonical_calibration_v4_7.py',n),'inherited_executor_function_AST:'+n)
 old_guard=ast.parse((S/'42_prepare_native_canonical_calibration_v4_6.py').read_text());new_guard=ast.parse((S/'42_prepare_native_canonical_calibration_v4_7.py').read_text())
 og=next(n for n in old_guard.body if isinstance(n,ast.FunctionDef) and n.name=='final_violations');ng=next(n for n in new_guard.body if isinstance(n,ast.FunctionDef) and n.name=='final_violations')
 for node in ast.walk(og):
  if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='campaign_bytes' for t in node.targets):node.value=ast.parse('global_campaign_bytes(SSD)',mode='eval').body
 check(ast.dump(og,include_attributes=False)==ast.dump(ng,include_attributes=False),'resource_guard_only_strengthened_global_regular_meter_helper')
 for n in ['catchable_termination','deferred_termination_signals','assert_no_termination','safe_diagnostic','exclusive_heavy_lock','validate_inherited_lock','expected_historical_dependencies','historical_completion']:
  check(function(S/'canonical_calibration_common_v4_6.py',n)==function(S/'canonical_calibration_common_v4_7.py',n),'inherited_common_function_AST:'+n)
 check(plan['preserved_failed_attempt_receipt']==str(FAILED) and plan['preserved_failed_attempt_receipt_sha256']==FPIN,'failed_attempt_receipt_bound')
 preserved=failed['all_preserved_file_identity'];check(len(preserved)==failed['file_count']==92,'exact92_preserved_file_count')
 check({p:r['sha256'] for p,r in preserved.items()}==plan['preserved_v4_6_file_sha256'],'complete_preserved_map_exact_plan')
 total=0
 for p,r in preserved.items():
  hashed(p,r['sha256']);check(Path(p).stat().st_size==r['bytes'],'preserved_byte_size:'+p);total+=r['bytes']
 check(total==failed['preserved_bytes']==127954369,'preserved127954369_bytes')
 check(failed['completed_native_controls_credited']==0 and failed['method_calibration_pass'] is False and failed['PENDING_preserved'] is True and Path(failed['pending_path']).is_file(),'old_failure_pending_and_no_controls_credit')
 check(failed['monitor_returncode']==0 and failed['monitor_stop_reason'] is None and failed['monitor_owned_group_disappearance_verified'] is True and failed['failure_scope']=='EXACTLY_MONITOR_APPLEDOUBLE_CREATED_AFTER_MONITOR_OUTPUT_MAP_FREEZE','preserved_transport_only_failure_qualified')
 check(len(failed['old_review_artifact_sha256'])==21 and failed['old_review_artifact_sha256']==plan['preserved_v4_6_review_artifact_sha256'],'old21_review_map_exact')
 for p,h in failed['old_review_artifact_sha256'].items():hashed(p,h)
 seal6=load(R/'independent_canonical_executor_prelaunch_review_v4_6_seal.json',failed['old_review_seal_sha256'])
 seal7=load(R/'genomicsem_canonical_operational_preparation_v4_7_seal.json','85b1957b4a63c554bbd1ad5a6605f4f87d559e5ac11b70d12f1845cc43cc93c8')
 check(seal7['file_count']==21 and len(seal7['file_sha256'])==21,'author21_artifact_map')
 for p,h in seal7['file_sha256'].items():hashed(p,h)
 check(seal7['root_execution_admission'] is False and seal7['workers_launched']==0 and seal7['scientific_sampling_calibration_pass'] is False,'author_seal_not_independent_or_science')
 # Current190 metadata exactness; reuse already sealed native science/output
 # evidence, rather than rereading 190 delete arrays or original GWAS bodies.
 historical=load(plan['historical_plan_path'],plan['historical_plan_sha256'])
 check(len(historical['jobs'])==190,'current_native190_jobs')
 admitted6=load(failed['preserved_admission_path'],failed['preserved_admission_sha256'])
 prior193=admitted6['historical_campaign_receipt_sha256'];check(len(prior193)==193,'preserved_v6_admitted190plus3_receipt_map')
 dep=dict(historical['dependencies_sha256']);original=str(P.parent/'sleep_unified_research_v1/scripts/native_ldsc_capture.py');relocated=str(Path(historical['ssd_support_package'])/'scripts/native_ldsc_capture.py')
 value=dep.pop(original);check(relocated not in dep and historical['support_file_sha256'][relocated]==value,'exact_one_key_original_capture_relocation');dep[relocated]=value;hashed(original,value);hashed(relocated,value)
 identities={r['path']:r['actual_sha256'] for r in historical['inputs_verified'] if r['match']};counts=Counter();current={}
 for j in historical['jobs']:
  output_dir=Path(historical['ssd_support_package'])/'native'/(j['stage']+'_reproduction_v1');path=output_dir/(j['job_id']+'.execution_receipt.json');r=load(path,prior193[str(path)]);current[str(path)]=sha(path);counts[j['stage']]+=1
  check(r['job']==j and r['returncode']==0 and r['scientific_cardinality_gate_pass'] is True and r['execution_identity_gate_pass'] is True,'completed_exact_job:'+j['job_id'])
  check(r['dependency_sha256_before']==dep==r['dependency_sha256_after'],'exact_dependencies:'+j['job_id'])
  inputs={p:identities[p] for p in j['inputs']};check(r['input_sha256']==inputs==r['input_sha256_after'],'exact_source_input_map:'+j['job_id'])
  mapping=r['all_output_sha256'];fp=str(output_dir/(j['job_id']+'.full_precision.json'))
  check(bool(mapping) and mapping[fp]==r['output_sha256'] and all(Path(p).parent==output_dir and Path(p).name.startswith(j['job_id']) for p in mapping),'output_fullprecision_scope:'+j['job_id'])
 check(counts=={'core':57,'extension':112,'validation':21},'exact57plus112plus21')
 sys.path.insert(0,str(S));import native_stage_completion_v4_3 as stage
 for name in counts:
  path=stage.stage_monitor_path(name);r=load(path,prior193[str(path)]);current[str(path)]=sha(path)
  check(r['returncode']==0 and r['stop_reason'] is None and r['process_group_teardown']['remaining_group_members']==[] and r['stage']==name and r['plan_sha256']==plan['historical_plan_sha256'],'current_completed_reaped_monitor:'+name)
 check(current==prior193,'all_current193_metadata_matches_preserved_admitted_complete_campaign')
 # Hash small actual code/metadata deps; record source/reference/binary maps
 # untouched. The candidate's live admit() still requires their full hashes.
 for p,h in plan['dependencies_sha256'].items():
  q=Path(p)
  if p in HASHED:check(HASHED[p]==h,'already_hashed_plan_dependency:'+p)
  elif q.suffix in ['.py','.json','.md','.tsv','.sh','.diff','.txt','.log'] and q.is_relative_to(P):hashed(q,h)
  else:DEFERRED[p]=h
 check(len(plan['dependencies_sha256'])==187,'exact187_frozen_dependency_entries')
 check(not (P/'manifests/native_canonical_calibration_admission_v4_7.json').exists(),'actual_v7_admission_absent')
 check(not Path(plan['ssd_output_root']).exists(),'actual_v7_execution_namespace_absent')
 check(sha(PLAN)==PIN and sha(FAILED)==FPIN,'final_plan_failed_inventory_unchanged')
 elapsed=time.monotonic()-START;rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 check(elapsed<180 and rss<128*1024*1024,'bounded_review_resource')
 record={'schema':'independent_canonical_v4_7_binding_review','status':'BINDING_AND_UNCHANGED_METHOD_PASS_SEPARATE_MONITOR_ORACLE_REJECTION','check_count':len(CHECKS),'checks':CHECKS,'frozen_plan_sha256':PIN,'failed_inventory_sha256':FPIN,'current190_metadata_stage_counts':dict(counts),'current193_receipt_sha256':current,'hashed_code_metadata_preserved_output_sha256':HASHED,'deferred_source_reference_runtime_dependencies_sha256':DEFERRED,'deferred_dependencies_count':len(DEFERRED),'preserved_generated_bytes_verified':total,'peak_rss_bytes':rss,'elapsed_seconds':elapsed,'checker_sha256':sha(__file__),'actual_workers':0,'actual_fits':0,'actual_mutex':0,'GWAS_or_reference_body_reads':0,'fresh190_numerical_audit':False,'execution_admission_granted':False,'qualification':'Current193 receipt bytes/source/command/dependency/output maps match the prior admitted sealed campaign. Native numeric/output-content evidence is inherited unchanged, not rerun. This does not credit any new canonical native controls or realistic-LD scientific calibration.'}
 with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 print(json.dumps({'receipt_sha256':sha(OUT),'checks':len(CHECKS),'metadata193':len(current),'deferred_dependencies':len(DEFERRED),'peak_rss_bytes':rss,'elapsed_seconds':elapsed}))
if __name__=='__main__':main()
