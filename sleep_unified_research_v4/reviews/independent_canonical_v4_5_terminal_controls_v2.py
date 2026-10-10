#!/usr/bin/env python3
"""Private metadata-only canonical controller faults; no worker, fit or mutex."""
import copy
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
from types import SimpleNamespace
sys.dont_write_bytecode=True
P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4');R=P/'reviews'
PLAN=P/'manifests/native_canonical_calibration_plan_v4_5.json';PIN='eef6ae0855b8e798b0fbe3a2956b565cd9fdc87ea81dc8007104df51fe8a39c6'
CODE=P/'scripts/42_prepare_native_canonical_calibration_v4_5.py';CODEPIN='3aac60b37a9304f15105c56eb389b938663fd77c4e775ec1be54968544b9b16b'
FIX=R/'independent_canonical_v4_5_terminal_controls_v2';RECORDS=[]
sys.path.insert(0,str(P/'scripts'))
def sha(q):return hashlib.sha256(Path(q).read_bytes()).hexdigest()
def save(q,x):
 q=Path(q);q.parent.mkdir(parents=True,exist_ok=True)
 with q.open('x') as f:json.dump(x,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def module(q,name):
 s=importlib.util.spec_from_file_location(name,q);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def case(name):
 d=FIX/name;d.mkdir();private=d/'exact_executor_copy.py';private.write_bytes(CODE.read_bytes());assert sha(private)==CODEPIN
 m=module(private,'_canonical_'+name);common=sys.modules['canonical_calibration_common_v4_5'];common.TERMINATION_REQUEST.clear()
 m.SSD=d/'campaign/statistical_validation/native_canonical_calibration_v4_5';m.SSD.parent.mkdir(parents=True)
 p=copy.deepcopy(json.loads(PLAN.read_text()));old=p['ssd_output_root'];p['ssd_output_root']=str(m.SSD);pp=d/'plan.json';ap=d/'admission.json'
 p['self_path']=str(pp);p['dependencies_sha256']={}
 for j in p['jobs']:
  j['output_dir']=j['output_dir'].replace(old,str(m.SSD));j['ldsc_args']=[x.replace(old,str(m.SSD)) for x in j['ldsc_args']]
 save(pp,p);save(ap,{'fixture_only_not_admission':True});m.PLAN=pp;m.ADMISSION=ap
 initial={'plan':sha(pp),'admission':sha(ap),'executor':sha(private)}
 m.admit=lambda:(json.loads(pp.read_text()),json.loads(ap.read_text()))
 original_history=common.historical_completion;common.historical_completion=lambda p:{'MOCK190_METADATA_ONLY':'not scientific evidence'}
 held=[False];lifecycle=[]
 @contextmanager
 def no_mutex(*,before_release=None):
  held[0]=True;lifecycle.append('enter_metadata_ownership_context')
  try:yield -777
  finally:
   if before_release:before_release()
   held[0]=False;lifecycle.append('exit_metadata_ownership_context')
 m.exclusive_heavy_lock=no_mutex
 m.monitor_helpers=lambda:(_ for _ in ()).throw(AssertionError('REAL_MONITOR_FORBIDDEN'))
 roles={};worker_calls=[];postpersist=[False];supplemental_failures=[]
 def worker(plan,admission_sha,job_id,command,out,fd,ownership):
  assert held[0] and fd==-777 and ownership==[]
  worker_calls.append(job_id)
  stdout=out/'worker_stdout.log';stdout.write_text('METADATA ONLY; NO WORKER\n');roles.setdefault('output',stdout)
  if job_id=='INDEPENDENT_CAPTURED_ROW_ARITHMETIC':
   ar=m.SSD/'independent_control_arithmetic.json';save(ar,{'fixture_only_not_scientific_evidence':True,'all_checks_pass':True});roles['arithmetic']=ar
  else:
   j=next(j for j in plan['jobs'] if j['job_id']==job_id);capture=out/'native_capture_receipt.json';save(capture,{'job':j,'plan_sha256':sha(pp),'fixture_only_not_scientific_evidence':True});roles.setdefault('capture',capture)
  files={str(q):sha(q) for q in out.rglob('*') if q.is_file()}
  receipt=out/'monitor_receipt.json'
  r={'returncode':0,'stop_reason':None,'process_group_teardown':{'remaining_group_members':[]},'job_id':job_id,'plan_sha256':sha(pp),'admission_sha256':sha(ap),'final_guard_violations':[],'owned_group_disappearance_verified':True,'bulk_output_hashes_completed_before_final_guards':True,'output_hash_error':None,'deferred_termination_requests':[],'all_output_sha256':files,'fixture_only_not_scientific_evidence':True}
  if name=='uncleared_native_monitor' and len(worker_calls)==1:r['process_group_teardown']['remaining_group_members']=[{'fixture_only_no_process':True}]
  save(receipt,r);roles.setdefault('monitor',receipt)
  if job_id=='INDEPENDENT_CAPTURED_ROW_ARITHMETIC':
   roles['arithmetic_monitor']=receipt
   if name=='required_capture_removed_after_completed_workers':roles['capture'].unlink()
  elif name=='required_capture_missing_after_worker' and len(worker_calls)==1:roles['capture'].unlink()
 m.monitored_worker=worker
 m.snapshot=lambda root:{'internal_free_bytes':0 if postpersist[0] and name=='resource_after_persistence' else 1<<50,'ssd_free_bytes':1<<50}
 tm=module(P/'scripts/terminal_commit_common_v2.py','_canonical_terminal_'+name);original_save=tm.save_new
 pending=m.SSD.parent/'native_canonical_calibration_pending_v4_5.json';seal=m.SSD.parent/'native_canonical_calibration_terminal_seal_v4_5.json';result=m.SSD.parent/'native_canonical_calibration_execution_receipt_v4_5.json'
 def terminal_save(q,value):
  if str(q).endswith('.failure.json') and name!='success':supplemental_failures.append(str(q));raise OSError('METADATA_FAULT_SUPPLEMENTAL_WRITE_FAILURE')
  original_save(q,value)
  if Path(q)==seal:
   assert pending.is_file() and held[0];postpersist[0]=True
   if name.startswith('mutate_'):
    key=name.removeprefix('mutate_');target={'plan':pp,'admission':ap,'executor':private}.get(key,roles.get(key))
    with target.open('a') as f:f.write('\nPOST_PERSISTENCE_METADATA_MUTATION\n')
   if name=='symlink_output_after_persistence':
    target=roles['output'];clone=d/'same_output_bytes';clone.write_bytes(target.read_bytes());target.unlink();target.symlink_to(clone)
   if name=='deferred_signal_after_persistence':common.catchable_termination(signal.SIGTERM,None)
 tm.save_new=terminal_save;m.TerminalCommit=tm.TerminalCommit
 m.print=lambda *args,**kwargs:None
 argv=sys.argv;sys.argv=[str(private),'execute'];before_handlers={s:signal.getsignal(s) for s in [signal.SIGTERM,signal.SIGHUP,signal.SIGINT]}
 try:
  try:m.main();succeeded=True;error=None
  except BaseException as e:succeeded=False;error=type(e).__name__+': '+str(e)
 finally:sys.argv=argv;common.historical_completion=original_history
 restored=all(signal.getsignal(s)==h for s,h in before_handlers.items());common.TERMINATION_REQUEST.clear()
 consumer=module(P/'scripts/terminal_commit_common_v2.py','_canonical_consumer_'+name)
 binding={'plan_sha256':initial['plan'],'admission_sha256':initial['admission'],'executor_sha256':initial['executor']}
 try:consumer.require_committed(pending,seal,binding,{str(result):sha(result)});accepted=True
 except BaseException as e:accepted=False
 expected=name in ['success','symlink_output_after_persistence','required_capture_removed_after_completed_workers']
 assert succeeded==accepted==expected,(name,succeeded,accepted,error)
 assert not held[0] and restored
 finding=name in ['symlink_output_after_persistence','required_capture_removed_after_completed_workers']
 record=dict(case=name,fixture_expectation_verified=True,controller_success=succeeded,consumer_accepts=accepted,pending_retained=pending.exists() or pending.is_symlink(),seal_persisted=seal.is_file(),result_exists=result.is_file(),controller_error=error,metadata_worker_stub_calls=len(worker_calls),actual_workers=0,actual_fits=0,actual_mutex_acquisitions=0,supplemental_write_failures=supplemental_failures,signal_handlers_restored=restored,metadata_ownership_context=lifecycle,material_gap_witness=finding)
 if result.is_file():
  r=json.loads(result.read_text());record['final_output_count']=len(r['all_output_sha256']);record['first_required_capture_in_final_inventory']=str(roles.get('capture')) in r['all_output_sha256'];record['output_is_symlink']=roles['output'].is_symlink()
 RECORDS.append(record)

def main():
 assert sha(PLAN)==PIN and sha(CODE)==CODEPIN;FIX.mkdir(exist_ok=False)
 for name in ['success','mutate_plan','mutate_admission','mutate_executor','mutate_output','mutate_capture','mutate_monitor','mutate_arithmetic','mutate_arithmetic_monitor','resource_after_persistence','deferred_signal_after_persistence','required_capture_missing_after_worker','uncleared_native_monitor','symlink_output_after_persistence','required_capture_removed_after_completed_workers']:case(name)
 assert sha(PLAN)==PIN and sha(CODE)==CODEPIN
 out=R/'independent_canonical_v4_5_terminal_controls_receipt_v2.json';save(out,{'schema':'independent_canonical_v4_5_metadata_terminal_controls','executor_sha256':CODEPIN,'plan_sha256':PIN,'checker_sha256':sha(__file__),'bounded_controls':len(RECORDS),'all_fixture_expectations_verified':all(r['fixture_expectation_verified'] for r in RECORDS),'material_gap_witness_count':sum(r['material_gap_witness'] for r in RECORDS),'controls':RECORDS,'fixture_regular_file_sha256':{str(q):sha(q) for q in FIX.rglob('*') if q.is_file() and not q.is_symlink()},'fixture_symlink_targets':{str(q):os.readlink(q) for q in FIX.rglob('*') if q.is_symlink()},'GWAS_body_reads':0,'reference_body_reads':0,'real_data_decompressions':0,'actual_workers':0,'actual_fits':0,'actual_mutex_acquisitions':0,'execution_admission_granted':False});print(json.dumps({'controls':len(RECORDS),'gap_witnesses':sum(r['material_gap_witness'] for r in RECORDS),'receipt_sha256':sha(out)},indent=2))
if __name__=='__main__':main()
