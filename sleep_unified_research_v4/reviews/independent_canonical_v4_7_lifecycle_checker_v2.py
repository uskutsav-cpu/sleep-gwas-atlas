#!/usr/bin/env python3
"""Tiny own signal/cleanup/meter checks; fake process ownership, no mutex."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
from types import SimpleNamespace
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S))
import canonical_calibration_common_v4_7 as c
F=P.parents[1]/'independent_canonical_v4_7_lifecycle_fixtures_v2';OUT=R/'independent_canonical_v4_7_lifecycle_receipt.json';CHECKS=[];CASES=[]
def check(value,label):
 CHECKS.append({'label':label,'passed':bool(value)})
 if not value:raise AssertionError(label)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def module(name):
 sp=importlib.util.spec_from_file_location(name,S/'42_prepare_native_canonical_calibration_v4_7.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def cleanup_case(name,fail_first_fallback=False,proof_failure=False):
 d=F/name;d.mkdir();m=module('_independent_lifecycle_'+name);held=[True];group=[True];calls=[];waits=[0];safe_calls=[]
 class Proc:
  pid=777777777
  def wait(self,timeout):
   check(held[0] and timeout==10,name+':wait_inside_ownership');waits[0]+=1
   if fail_first_fallback and waits[0]==1:raise TimeoutError('OWN_FAKE_GROUP_STILL_LIVE')
   group[0]=False;return -9
 def helper(proc):
  check(held[0] and proc.pid==777777777,name+':helper_exact_owned_fake_pid');raise RuntimeError('OWN_HELPER_ERROR')
 m.monitor_helpers=lambda:SimpleNamespace(terminate_owned=helper)
 m.owned_group_exists=lambda pid:group[0]
 # Replace the controller's os binding, never actual os.killpg.
 m.os=SimpleNamespace(killpg=lambda pid,sig:(check(held[0] and pid==777777777 and sig==signal.SIGKILL,name+':fallback_exact_owned_group'),calls.append((pid,int(sig)))))
 m.safe_diagnostic=lambda value:(check(held[0],name+':diagnostic_inside_retained_ownership'),safe_calls.append(value()))
 m.time=SimpleNamespace(sleep=lambda seconds:check(held[0] and seconds==2 and group[0],name+':quarantine_retains_ownership_while_fake_group_live'))
 if proof_failure:m.write_new=lambda *a:(_ for _ in ()).throw(OSError('OWN_RECOVERY_PROOF_ERROR'))
 entry={'proc':Proc(),'output_dir':d,'cleanup_errors':[]};owned=[entry];m.await_owned_cleanup(owned,None)
 check(not group[0] and owned==[],name+':no_owned_group_before_release');held[0]=False
 check(entry['verified_recovery']['signal_zero_group_absent'] is True and len(entry['cleanup_errors'])>=1,name+':recovery_failure_recorded')
 if proof_failure:check(any('RECOVERY_PROOF_WRITE' in e for e in entry['cleanup_errors']),name+':proof_failure_preserved_after_group_gone')
 if fail_first_fallback:check(waits[0]==2 and len(calls)==2 and safe_calls,name+':failed_fallback_retry_remains_protected')
 CASES.append({'case':name,'fake_owned_pid':777777777,'real_process':False,'real_mutex':False,'helper_failure_retained':entry['cleanup_errors'],'fallback_wait_count':waits[0],'ownership_empty':owned==[],'group_gone':not group[0],'diagnostics':safe_calls})
def main():
 check(not F.exists() and not OUT.exists(),'fresh_private_lifecycle_namespace');F.mkdir()
 for sig in [signal.SIGTERM,signal.SIGHUP,signal.SIGINT]:
  c.TERMINATION_REQUEST.clear();before=signal.getsignal(sig)
  with c.deferred_termination_signals():
   os.kill(os.getpid(),sig)
   check(c.TERMINATION_REQUEST and c.TERMINATION_REQUEST[-1]==int(sig),'real_own_signal_deferred:'+str(sig))
   try:c.assert_no_termination();rejected=False
   except RuntimeError:rejected=True
   check(rejected,'protected_gate_rejects_signal:'+str(sig))
  check(signal.getsignal(sig)==before,'handler_restored:'+str(sig));c.TERMINATION_REQUEST.clear()
  CASES.append({'case':'real_own_deferred_signal_'+str(sig),'real_signal_only_own_reviewer_process':True,'actual_worker':False})
 for typ in [BrokenPipeError,KeyboardInterrupt]:
  def fail():raise typ('OWN_DIAGNOSTIC_FAILURE')
  c.safe_diagnostic(fail);check(True,'nonthrowing_diagnostic:'+typ.__name__)
  CASES.append({'case':'safe_diagnostic_'+typ.__name__,'actual_worker':False})
 cleanup_case('helper_exception_fallback');cleanup_case('failed_fallback_retains_ownership',True);cleanup_case('cleanup_recovery_proof_failure',proof_failure=True)
 root=F/'meter/campaign/statistical_validation/native_canonical_calibration_v4_7';old=root.parent/'native_canonical_calibration_v4_6';old.mkdir(parents=True);root.mkdir();(old/'old').write_bytes(b'A'*41);(root/'new').write_bytes(b'B'*19)
 (root.parent/'native_canonical_calibration_pending_v4_6.json').write_bytes(b'C'*23)
 (root.parent/'._native_canonical_calibration_pending_v4_6.json').write_bytes(b'D'*31)
 value=c.canonical_epoch_bytes(root);check(value==41+19+23+31,'all_epochs_and_retained_metadata_byte_meter')
 check(c.global_campaign_bytes(root)==value,'global_meter_complete_private_campaign')
 link=root/'bad_symlink';link.symlink_to(old/'old')
 try:c.canonical_epoch_bytes(root);reject=False
 except RuntimeError:reject=True
 check(reject,'epoch_meter_rejects_symlink');CASES.append({'case':'all_epochs_meter_and_symlink','expected_regular_bytes':value,'symlink_rejected':reject})
 record={'schema':'independent_canonical_v4_7_lifecycle_controls','case_count':len(CASES),'check_count':len(CHECKS),'checks':CHECKS,'cases':CASES,'checker_sha256':sha(__file__),'source_sha256':{str(S/n):sha(S/n) for n in ['42_prepare_native_canonical_calibration_v4_7.py','canonical_calibration_common_v4_7.py']},'fixture_regular_sha256':{str(p):sha(p) for p in F.rglob('*') if p.is_file() and not p.is_symlink()},'actual_workers':0,'actual_mutex':0,'actual_fits':0,'source_reference_reads':0,'qualification':'Fallback cleanup exercised against explicitly fake owned process/group and held metadata context. Real self-signals only target the reviewer. Identical v6 cleanup/signal AST and previously sealed real controls supply inherited real-child evidence; no new production or hostile-writer guarantee.'}
 with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'receipt_sha256':sha(OUT),'cases':len(CASES),'checks':len(CHECKS)}))
if __name__=='__main__':main()
