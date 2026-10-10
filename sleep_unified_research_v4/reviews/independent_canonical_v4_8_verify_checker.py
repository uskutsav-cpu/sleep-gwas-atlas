#!/usr/bin/env python3
"""Independently isolate actual final verify consumer with fixed own hashes."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import struct
import sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S))
F=P.parents[1]/'independent_canonical_v4_8_verify_fixtures';OUT=R/'independent_canonical_v4_8_verify_receipt.json';CASES=[];CHECKS=[]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 return sha(p)
def check(v,label):
 CHECKS.append({'label':label,'passed':bool(v)})
 if not v:raise AssertionError(label)
class OwnNoLaunchBoundary(Exception):pass
def main():
 check(not F.exists() and not OUT.exists(),'fresh_private_verify_namespace');F.mkdir()
 source=S/'42_prepare_native_canonical_calibration_v4_8.py';pin='c90c088b5557dbcc995914eea928453ec12662f80175fa52ee3b1c23a98676d0';check(sha(source)==pin,'actual_fixed_v8_source')
 for name in ['healthy','healthy_monitor_companion','cleanup_errors','run_error','teardown_cleanup_error','first_exception','first_exception_null']:
  d=F/name;d.mkdir();root=d/'native_canonical_calibration_v4_8';root.mkdir();out=root/'native';out.mkdir();(out/'weighted_capture').mkdir()
  sp=importlib.util.spec_from_file_location('_independent8_verify_'+name,source);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
  m.PLAN=d/'plan.json';m.ADMISSION=d/'admission.json';save(m.PLAN,{'fixture_only':name});save(m.ADMISSION,{'fixture_only_not_admission':name})
  j={'job_id':'OWN_VERIFY','output_dir':str(out)};plan={'jobs':[j],'ssd_output_root':str(root),'python':'NO_RUNTIME_OR_WORKER'};m.admit=lambda:(plan,{})
  stdout=out/'worker_stdout.log';stdout.write_text('PRIVATE METADATA ONLY\n');cap=out/'native_capture_receipt.json';save(cap,{'job':j,'plan_sha256':sha(m.PLAN)})
  outputs={str(stdout):sha(stdout),str(cap):sha(cap)}
  monitor=out/'monitor_receipt.json';record={'job_id':j['job_id'],'command':['OWN_METADATA_EXPECTED_COMMAND'],'plan_sha256':sha(m.PLAN),'admission_sha256':sha(m.ADMISSION),'returncode':0,'stop_reason':None,'process_group_teardown':{'remaining_group_members':[]},'final_guard_violations':[],'owned_group_disappearance_verified':True,'bulk_output_hashes_completed_before_final_guards':True,'output_hash_error':None,'deferred_termination_requests':[],'cleanup_errors':[],'run_error':None,'first_teardown_attempt':{'signal_zero_group_absent':True},'all_output_sha256':outputs}
  if name=='cleanup_errors':record['cleanup_errors']=['EXPLICIT_ERROR']
  if name=='run_error':record['run_error']={'type':'RuntimeError','repr':'EXPLICIT_ERROR'}
  if name=='teardown_cleanup_error':record['process_group_teardown']['cleanup_error']='EXPLICIT_ERROR'
  if name in ['first_exception','first_exception_null']:record['first_teardown_attempt']['exception']=None if name.endswith('_null') else {'type':'RuntimeError','repr':'EXPLICIT_ERROR'}
  save(monitor,record);frozen=dict(outputs,**{str(monitor):sha(monitor)})
  if name=='healthy_monitor_companion':
   companion=out/'._monitor_receipt.json';companion.write_bytes(struct.pack('>II16sH',0x00051607,0x00020000,b'\0'*16,1)+struct.pack('>III',9,38,4)+b'TEST');frozen[str(companion)]=sha(companion)
  # Intentionally synthesize the fixed already-frozen consumer input, without
  # routing contradictions through the now-correct initial freeze. This tests
  # the independent final oracle, not a claim of actual producer reachability.
  evidence={'workers':{j['job_id']:{'output_dir':str(out),'sha256':frozen,'registered_directories':[str(out),str(out/'weighted_capture')],'external_output_paths':[],'registered_deleted_markers':[]}},'all_output_sha256':dict(frozen)}
  calls=[]
  def no_arithmetic(*a,**k):calls.append('own_no_launch_boundary');raise OwnNoLaunchBoundary('ACTUAL_VERIFY_PASSED_MONITOR_NO_SCIENCE_LAUNCH')
  m.monitored_worker=no_arithmetic
  try:m.verify(-777,[],evidence);error=None;boundary=False
  except BaseException as e:error=type(e).__name__+': '+str(e);boundary=isinstance(e,OwnNoLaunchBoundary)
  expected=name.startswith('healthy')
  check(boundary==expected,name+':actual_final_consumer_outcome')
  if expected:check(calls==['own_no_launch_boundary'],name+':healthy_reaches_only_explicit_no_launch_boundary')
  else:check(not calls and 'CANONICAL_MONITOR_RECORDED_FAILURE' in error,name+':error_veto_before_arithmetic')
  check(all(sha(p)==h for p,h in frozen.items()),name+':fixed_preexisting_receipt_and_output_hashes_unchanged')
  CASES.append({'case':name,'healthy_reaches_no_launch_boundary':boundary,'error':error,'arithmetic_boundary_stub_calls':len(calls),'actual_workers':0,'exact_monitor_sha256':sha(monitor),'frozen_current_regular_sha256':frozen})
 check(sha(source)==pin,'source_still_unchanged')
 result={'schema':'independent_canonical_v4_8_final_verify_controls','verdict':'ALL_FIVE_CONTRADICTIONS_REJECT_BOTH_HEALTHY_FINAL_MONITORS_ACCEPT','cases':CASES,'case_count':len(CASES),'checks':CHECKS,'check_count':len(CHECKS),'checker_sha256':sha(__file__),'executor_sha256':pin,'fixture_sha256':{str(p):sha(p) for p in F.rglob('*') if p.is_file()},'qualification':'Actual final verify function with fixed synthesized already-frozen records; one native metadata job isolates the loop, both healthy variants stop at explicit no-launch arithmetic boundary. Whole16plus1 healthy Terminal2 closure is tested separately. Contradictory records do not originate from ordinary producer or corrected initial freeze.','actual_workers':0,'actual_mutex':0,'body_reads':0,'execution_admission':False}
 save(OUT,result);print(json.dumps({'receipt_sha256':sha(OUT),'cases':len(CASES),'assertions':len(CHECKS)}))
if __name__=='__main__':main()
