#!/usr/bin/env python3
"""Narrow independent actual8 binding review; no source/reference/runtime reads."""
import ast
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S))
PLAN=P/'manifests/native_canonical_calibration_plan_v4_8.json';PIN='c9ad3842526a61ab24ff8fd39222c408fb611ba2bbf26e7ef0959f1e86d8de51'
OUT=R/'independent_canonical_v4_8_binding_receipt.json';HASHED={};CHECKS=[];START=time.monotonic()
def check(v,label):
 CHECKS.append({'label':label,'passed':bool(v)})
 if not v:raise AssertionError(label)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def hashed(p,h):
 p=Path(p);check(p.is_file() and not p.is_symlink() and not any(q.is_symlink() for q in p.parents),'regular:'+str(p));got=sha(p);check(got==h,'SHA:'+str(p));HASHED[str(p)]=got;return got
def load(p,h):
 hashed(p,h);r=json.loads(Path(p).read_text());check(sha(p)==h,'postparse:'+str(p));return r
def functions(path,normalize=False):
 text=path.read_text()
 if normalize:text=text.replace('v4_7','v4_8')
 return {n.name:n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
def same(a,b):return ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
def main():
 check(not OUT.exists(),'fresh_receipt');plan=load(PLAN,PIN)
 old=load(P/'manifests/native_canonical_calibration_plan_v4_7.json','a72501b3daf01fd4102a0f002fd78ca232ff2429b6fe84fb5903109503ab19b9')
 expected=copy.deepcopy(old['jobs'])
 for j in expected:
  j['output_dir']=j['output_dir'].replace(old['ssd_output_root'],plan['ssd_output_root']);j['ldsc_args']=[s.replace(old['ssd_output_root'],plan['ssd_output_root']) for s in j['ldsc_args']]
 check(expected==plan['jobs'] and len(expected)==16,'exact16_same_science_only_private_epoch_paths')
 for field in ['scope','allowed_control_ids','historical_plan_path','historical_plan_sha256','block_receipt_path','block_receipt_sha256','build','interval_tsv_sha256','ordered_all_reference_records_sha256','ordered_construction_records_sha256','coordinate_map_path','coordinate_map_sha256','partition_amendment','static_construction_exclusions','reference_sha256','inputs_current_sha256','python','ldsc_dir','environment','worker_count','blas_threads','shared_lock_path','shared_lock_policy','internal_floor_bytes','ssd_floor_bytes','worker_rss_limit_bytes','output_limit_bytes','per_worker_seconds_limit','runtime_poll_seconds','arithmetic','cross_arm_identity_rule','global_resource_ledger_path','global_resource_ledger_sha256','global_reservation_bytes','all_retained_canonical_epochs_output_cap_bytes','output_meter_policy','transport_companion_policy','terminal_protocol']:
  check(plan[field]==old[field],'unchanged_science_operation:'+field)
 check(plan['scope']=='METHOD_CALIBRATION_ONLY' and plan['preparation_only'] is True and plan['execution_admitted'] is False and plan['native_jobs_launched']==0 and plan['GWAS_outcome_columns_read'] is False and plan['allow_41_covariance_outcomes'] is False and plan['allow_calibrated_biological_p_values'] is False and plan['realistic_LD_sampling_calibration_pass'] is False,'no_native_or_biological_release')
 source={
 '42_prepare_native_canonical_calibration_v4_8.py':'c90c088b5557dbcc995914eea928453ec12662f80175fa52ee3b1c23a98676d0',
 'canonical_calibration_common_v4_8.py':'778ff2004b2dd06acdbc6eae05f960957187b37624c3bfc6dbf7bfb0e3c79379',
 'canonical_calibration_capture_v4_8.py':'9a4b1823df56af6195c4990ac2193973bb4183d24e3442899cd89720f41784f0',
 'canonical_calibration_arithmetic_v4_8.py':'50dd6f7986184a7f88415a781837ebb831dc2186805865b83801e018ad7b98bd'}
 for n,h in source.items():hashed(S/n,h);check(plan['dependencies_sha256'][str(S/n)]==h,'actual_plan_binds_source:'+n)
 check((S/'canonical_calibration_capture_v4_7.py').read_text().replace('v4_7','v4_8')==(S/'canonical_calibration_capture_v4_8.py').read_text(),'capture_exact_import_namespace_only')
 # Every runtime function remains exactly the prior AST except the two new
 # helper calls after fixed parse/hash. Preparation is separately reviewed.
 prior=functions(S/'42_prepare_native_canonical_calibration_v4_7.py',True);new=functions(S/'42_prepare_native_canonical_calibration_v4_8.py')
 check(set(new)==set(prior)|{'preparation_dependency_hashes'},'only_new_preparation_helper_function')
 for name,node in prior.items():
  if name=='prepare':continue
  candidate=copy.deepcopy(new[name])
  if name in ['freeze_worker_evidence','verify']:
   removed=[]
   for parent in ast.walk(candidate):
    for key,value in ast.iter_fields(parent):
     if not isinstance(value,list):continue
     for pos,n in enumerate(list(value)):
      if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='assert_monitor_error_fields_clear':
       check(pos>=2 and isinstance(value[pos-1],ast.Expr) and isinstance(value[pos-1].value,ast.Call) and isinstance(value[pos-1].value.func,ast.Name) and value[pos-1].value.func.id=='regular_hashes','oracle_after_postparse_fixed_hash:'+name)
       removed.append(n);value.remove(n)
   check(len(removed)==1,'exact_one_oracle_call:'+name)
  check(same(node,candidate),'unchanged_runtime_AST_except_two_oracle_calls:'+name)
 priorc=functions(S/'canonical_calibration_common_v4_7.py',True);newc=functions(S/'canonical_calibration_common_v4_8.py')
 check(set(newc)==set(priorc)|{'assert_monitor_error_fields_clear'},'only_new_common_monitor_oracle')
 for name,node in priorc.items():check(same(node,newc[name]),'unchanged_common_AST_including_admit_full_hash190:'+name)
 priora=functions(S/'canonical_calibration_arithmetic_v4_7.py',True);newa=functions(S/'canonical_calibration_arithmetic_v4_8.py')
 check(set(priora)==set(newa),'arithmetic_same_all_functions')
 for name,node in priora.items():check(same(node,newa[name]),'unchanged_entire_arithmetic_AST:'+name)
 rejection=load(plan['preserved_v4_7_complete_rejection_seal'],'dc310205bd6754708456fb7cdca56449f4995c5e8a6e3d38a7f0db777134c845')
 check(len(rejection['file_sha256'])==rejection['artifact_count']==25 and rejection['status']=='REJECT_PRESERVED_NO_EXECUTION_ADMISSION','whole25_prior_rejection_closure')
 check(rejection['file_sha256']==plan['preserved_v4_7_rejection_artifact_sha256'],'plan_exact_prior_rejection_map')
 for p,h in rejection['file_sha256'].items():hashed(p,h);check(plan['dependencies_sha256'][p]==h,'plan_binds_prior_rejection_leaf:'+p)
 original=load(R/'independent_canonical_executor_prelaunch_review_v4_7_seal.json','272919b1ea6440528c467b346ac4aac9033105b9eb42736a1f09d92f2d07d33a')
 check(len(original['file_sha256'])==original['artifact_count']==23 and all(rejection['file_sha256'].get(p)==h for p,h in original['file_sha256'].items()),'original23_artifacts_retained_in_complete25')
 author7=load(plan['preserved_v4_7_author_seal'],'85b1957b4a63c554bbd1ad5a6605f4f87d559e5ac11b70d12f1845cc43cc93c8')
 check(len(author7['file_sha256'])==21 and author7['file_sha256']==plan['preserved_v4_7_author_artifact_sha256'],'author7_exact21_closure')
 for p,h in author7['file_sha256'].items():hashed(p,h);check(plan['dependencies_sha256'][p]==h,'plan_binds_author7_leaf:'+p)
 author8=load(R/'genomicsem_canonical_operational_preparation_v4_8_seal.json','f46bf64d0a3ce14846dd8becdd6694c7cd7cd534d356e302f695a18c320171ab')
 check(len(author8['file_sha256'])==author8['file_count']==24 and author8['root_execution_admission'] is False and author8['scientific_sampling_calibration_pass'] is False,'author8_exact24_not_admission_or_science')
 for p,h in author8['file_sha256'].items():hashed(p,h)
 failed=load(plan['preserved_failed_attempt_receipt'],'c80990a119349b0deda5b4d1c6ee15ceb1d16f3c79d0f9f457894b8998674215')
 check(len(failed['all_preserved_file_identity'])==92 and failed['completed_native_controls_credited']==0 and failed['method_calibration_pass'] is False and failed['PENDING_preserved'] is True,'all92_failed6_and_zero_credits')
 total=0
 for p,r in failed['all_preserved_file_identity'].items():
  hashed(p,r['sha256']);check(Path(p).stat().st_size==r['bytes'] and plan['preserved_v4_6_file_sha256'][p]==r['sha256'],'failed6_size_and_planmap:'+p);total+=r['bytes']
 check(total==127954369==failed['preserved_bytes'],'failed6_preserved_byte_count')
 check(len(failed['old_review_artifact_sha256'])==21 and failed['old_review_artifact_sha256']==plan['preserved_v4_6_review_artifact_sha256'],'failed6_old21_review_closure')
 for p,h in failed['old_review_artifact_sha256'].items():hashed(p,h)
 add=load(plan['preparation_deferred_identity_addendum'],'fa02743dbba82390c7bb5fc15d9c479c5a9f0e7c16797bd603b8ed5dbfa4c0af');deferred=add['remaining50_source_reference_runtime_dependency_sha256']
 check(len(deferred)==50==plan['preparation_deferred_dependency_count'] and deferred==plan['preparation_deferred_dependency_sha256'] and all(old['dependencies_sha256'][p]==h==plan['dependencies_sha256'][p] for p,h in deferred.items()),'authenticated_exact50_deferred_prior_identities')
 check('NOT_FRESHLY_HASHED' in plan['preparation_dependency_qualification'] and 'FUTURE_ADMIT_MANDATORY_FULL_CURRENT' in plan['preparation_dependency_qualification'],'preparation_deferral_explicit_not_live_admission_bypass')
 check(len(plan['dependencies_sha256'])==232,'all232_plan_dependencies')
 current={p:h for p,h in plan['dependencies_sha256'].items() if p not in deferred}
 for p,h in current.items():hashed(p,h)
 check(len(current)==182,'fresh182_current_nonbody_dependencies')
 # Private helper contracts only; never invoke prepare(), admit(), a body hash
 # gate, runtime Python, or the actual shared mutex.
 sp=importlib.util.spec_from_file_location('_independent8_binding_helper',S/'42_prepare_native_canonical_calibration_v4_8.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
 called=[];m.check_hashes=lambda mapping:called.append(dict(mapping))
 fake={f'NO_BODY_{i}.gz':str(i).zfill(64) for i in range(50)};fp={'dependencies_sha256':dict(fake,OWN_META='a'*64)}
 m.preparation_dependency_hashes(fp,fake);check(called==[{'OWN_META':'a'*64}],'prep_helper_only_hashes_nondeferred_mock_map')
 for label,bad in [('49_entries',dict(list(fake.items())[:49])),('wrong_identity',dict(fake,**{'NO_BODY_0.gz':'f'*64}))]:
  try:m.preparation_dependency_hashes(fp,bad);rejected=False
  except RuntimeError:rejected=True
  check(rejected,'prep_deferral_reject:'+label)
 check(not (P/'manifests/native_canonical_calibration_admission_v4_8.json').exists() and not Path(plan['ssd_output_root']).exists(),'no_actual8_admission_or_worker_namespace')
 for p,h in HASHED.items():check(sha(p)==h,'final_preserved_readback:'+p)
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;elapsed=time.monotonic()-START;check(rss<128*2**20 and elapsed<180,'bounded_review_limits')
 record={'schema':'independent_canonical_v4_8_binding_review','verdict':'ACTUAL_FROZEN_PLAN_BINDING_AND_UNCHANGED_RUNTIME_METHOD_PASS','plan_sha256':PIN,'check_count':len(CHECKS),'checks':CHECKS,'checker_sha256':sha(__file__),'current_nonbody_dependency_count':182,'all_dependency_count':232,'fresh_hashed_code_metadata_generated_output_sha256':HASHED,'exact50_inherited_unread_source_reference_runtime_identity_sha256':deferred,'all92_preserved_bytes_verified':total,'v7_runtime_functions_identical_except_two_oracle_calls':True,'full_common_admit_hash_and_original190_functions_AST_identical':True,'whole_arithmetic_AST_namespace_only':True,'capture_import_namespace_only':True,'prior_complete25_and_original23_rejection_authenticated':True,'prior_author21_and_current_author24_authenticated':True,'fresh_original190_numerical_audit':False,'future_live_full_hash190_admission_still_required':True,'peak_rss_bytes':rss,'elapsed_seconds':elapsed,'body_reads':0,'actual_workers':0,'actual_mutex':0,'actual_preparer_called':False,'actual_admit_called':False,'execution_admission':False}
 with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 print(json.dumps({'receipt_sha256':sha(OUT),'assertions':len(CHECKS),'hashed_unique':len(HASHED),'current182':182,'deferred50':50,'elapsed_seconds':elapsed,'peak_rss_bytes':rss}))
if __name__=='__main__':main()
