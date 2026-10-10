#!/usr/bin/env python3
"""Independent private/POSIX metadata controller controls; never start a process.

The science and historical190 producers are explicit stubs. Every production
worker, mutex, GWAS/reference body and execution admission is forbidden here.
"""
import copy
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import struct
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_8 as c
CODE=S/'42_prepare_native_canonical_calibration_v4_8.py'
PLAN=P/'manifests/native_canonical_calibration_plan_v4_8.json'
FIX=P.parents[1]/'independent_canonical_v4_8_terminal_fixtures'
OUT=R/'independent_canonical_v4_8_terminal_receipt.json'
CODE_PIN='c90c088b5557dbcc995914eea928453ec12662f80175fa52ee3b1c23a98676d0'
PLAN_PIN='c9ad3842526a61ab24ff8fd39222c408fb611ba2bbf26e7ef0959f1e86d8de51'
RECORDS=[];CHECKS=[];START=time.monotonic()

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(65536),b''):h.update(block)
 return h.hexdigest()
def check(value,label):
 CHECKS.append({'label':label,'passed':bool(value)})
 if not value:raise AssertionError(label)
def save(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 payload=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()
 with path.open('xb') as f:f.write(payload);f.flush();os.fsync(f.fileno())
 return sha(path)
def module(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def apple(path):
 with Path(path).open('xb') as f:f.write(struct.pack('>II16sH',0x00051607,0x00020000,b'\0'*16,1)+struct.pack('>III',9,38,4)+b'TEST')
def drift(path):
 with Path(path).open('a') as f:f.write(' ')
def run_case(name):
 d=FIX/name;d.mkdir();private=d/'exact_executor_copy.py';private.write_bytes(CODE.read_bytes())
 check(sha(private)==CODE_PIN,name+':exact_executor')
 m=module(private,'_independent_canonical7_'+name);c.TERMINATION_REQUEST.clear()
 m.SSD=d/'campaign/statistical_validation/native_canonical_calibration_v4_8';m.SSD.parent.mkdir(parents=True)
 plan=copy.deepcopy(json.loads(PLAN.read_text()));old=plan['ssd_output_root'];plan['ssd_output_root']=str(m.SSD)
 pp=d/'plan.json';ap=d/'admission.json';plan['self_path']=str(pp);plan['dependencies_sha256']={}
 for j in plan['jobs']:
  j['output_dir']=j['output_dir'].replace(old,str(m.SSD));j['ldsc_args']=[s.replace(old,str(m.SSD)) for s in j['ldsc_args']]
 save(pp,plan);save(ap,{'metadata_fixture_only_not_admission':True});m.PLAN=pp;m.ADMISSION=ap
 initial={'plan_sha256':sha(pp),'admission_sha256':sha(ap),'executor_sha256':sha(private)}
 m.admit=lambda:(json.loads(pp.read_text()),json.loads(ap.read_text()))
 old_history=c.historical_completion;c.historical_completion=lambda p:{'EXPLICIT_MOCK190':'NO_NATIVE_EVIDENCE'}
 held=[False];lifecycle=[];calls=[];roles={};post=[False];failed_writes=[]
 @contextmanager
 def no_mutex(*,before_release=None):
  held[0]=True;lifecycle.append('metadata_ownership_begin')
  try:yield -777
  finally:
   if before_release:before_release()
   held[0]=False;lifecycle.append('metadata_ownership_end')
 m.exclusive_heavy_lock=no_mutex
 m.monitor_helpers=lambda:(_ for _ in ()).throw(AssertionError('REAL_PROCESS_HELPER_FORBIDDEN'))
 def stub(plan,admission_sha,jid,command,out,fd,ownership):
  check(held[0] and fd==-777 and ownership==[],name+':owned_stub_'+str(len(calls)))
  calls.append(jid);native=jid!='INDEPENDENT_CAPTURED_ROW_ARITHMETIC'
  stdout=out/'worker_stdout.log';stdout.write_text('PRIVATE METADATA STUB; ZERO WORKERS\n');roles.setdefault('output',stdout)
  files=[stdout]
  if native:
   j=next(x for x in plan['jobs'] if x['job_id']==jid)
   weighted=out/'weighted_capture';weighted.mkdir();tiny=weighted/'tiny_fixture';tiny.write_bytes(b'METADATA')
   cap=out/'native_capture_receipt.json';save(cap,{'job':j,'plan_sha256':sha(pp),'metadata_fixture_only':True});files.extend([tiny,cap]);roles.setdefault('capture',cap)
  else:
   for sub in ['tmp','cache','cache/matplotlib']:(m.SSD/sub).mkdir(parents=True,exist_ok=True)
   report=m.SSD/'independent_control_arithmetic.json'
   record={'all_checks_pass':True,'scope':'METHOD_CALIBRATION_IMPLEMENTATION_CHECK_ONLY','allow_41_covariance_outcomes':False,'calibrated_biological_p_values_computed':False,'metadata_fixture_only':True}
   if name=='bad_arithmetic_scope':record['scope']='INVALID_SCOPE'
   report_sha=save(report,record);roles['report']=report
   sources={str(Path(j['output_dir'])/'native_capture_receipt.json'):sha(Path(j['output_dir'])/'native_capture_receipt.json') for j in plan['jobs']}
   executor_sha=sha(S/'canonical_calibration_arithmetic_v4_8.py')
   intent={'schema':'canonical_arithmetic_report_intent_v4_8','executor_sha256':executor_sha,'plan_sha256':sha(pp),'admission_sha256':sha(ap),'report_path':str(report),'intended_report_sha256':report_sha,'capture_receipt_sha256':sources}
   proof=out/'arithmetic_report_intent.json';save(proof,intent)
   terminal=module(S/'terminal_commit_common_v2.py','_arithterminal_'+name)
   binding={'plan_sha256':sha(pp),'admission_sha256':sha(ap),'executor_sha256':executor_sha}
   tc=terminal.TerminalCommit(out/'arithmetic_pending.json',out/'arithmetic_terminal_seal.json',binding)
   intended=c.terminal_intended_sha(tc,{str(report):report_sha})
   check(tc.commit({str(report):report_sha},identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None),name+':arithmetic_private_commit')
   check(sha(tc.seal)==intended,name+':arithmetic_seal_intent')
   files.extend([proof,tc.seal,report]);roles['arithmetic_proof']=proof
   if name=='healthy_companions':apple(out/'._arithmetic_pending.json');files.append(out/'._arithmetic_pending.json')
  mapping={str(q):sha(q) for q in files}
  mon=out/'monitor_receipt.json'
  rec={'job_id':jid,'plan_sha256':sha(pp),'admission_sha256':sha(ap),'command':command,'returncode':0,'stop_reason':None,'process_group_teardown':{'remaining_group_members':[]},'final_guard_violations':[],'owned_group_disappearance_verified':True,'bulk_output_hashes_completed_before_final_guards':True,'output_hash_error':None,'deferred_termination_requests':[],'cleanup_errors':[],'run_error':None,'first_teardown_attempt':{'signal_zero_group_absent':True},'all_output_sha256':mapping}
  if len(calls)==1:
   if name=='contradictory_cleanup_errors':rec['cleanup_errors']=['EXPLICIT_ERROR']
   if name=='contradictory_run_error':rec['run_error']={'type':'RuntimeError','repr':'EXPLICIT_ERROR'}
   if name=='contradictory_teardown_error':rec['process_group_teardown']['cleanup_error']='EXPLICIT_ERROR'
   if name=='contradictory_first_exception':rec['first_teardown_attempt']['exception']={'type':'RuntimeError','repr':'EXPLICIT_ERROR'}
   if name=='contradictory_null_exception':rec['first_teardown_attempt']['exception']=None
   if name=='wrong_command':rec['command']=['WRONG']
   if name=='remaining_group':rec['process_group_teardown']['remaining_group_members']=[{'stub_no_real_process':True}]
   if name=='late_failure_flags':rec['final_guard_violations']=['EXPLICIT_FAILURE']
  intended=save(mon,rec);roles.setdefault('monitor',mon)
  if name=='healthy_companions':
   apple(out/'._monitor_receipt.json')
   if native:apple(out/'._weighted_capture')
   else:apple(m.SSD/'._independent_control_arithmetic.json')
  if len(calls)==1:
   if name=='monitor_before_parse':drift(mon)
   if name=='missing_capture':cap.unlink()
   if name=='extra_file':(out/'unregistered.json').write_text('{}')
   if name=='extra_directory':(out/'unknown').mkdir()
   if name=='directory_symlink':
    weighted.rename(d/'same_weighted');weighted.symlink_to(d/'same_weighted',target_is_directory=True)
   if name=='orphan_companion':apple(out/'._orphan')
   if name=='malformed_companion':(out/'._monitor_receipt.json').write_bytes(b'BAD')
   if name=='companion_symlink':(out/'._monitor_receipt.json').symlink_to(mon)
  if not native:
   if name=='missing_report':report.unlink()
   if name=='report_before_parse':drift(report)
   if name=='proof_before_parse':drift(proof)
  return {'path':str(mon),'sha256':intended}
 m.monitored_worker=stub
 # JSON hook mutates the real consumed private file after its read/parse, while
 # preserving the independently fixed intended digest. Never patch global json.
 parse_fired=[False]
 def loads(text,*args,**kwargs):
  value=json.loads(text,*args,**kwargs)
  if not parse_fired[0] and isinstance(value,dict):
   target=None
   if name=='monitor_after_parse' and value.get('job_id') and 'all_output_sha256' in value:target=roles['monitor']
   if name=='capture_after_parse' and value.get('metadata_fixture_only') and 'job' in value:target=roles['capture']
   if name=='report_after_parse' and value.get('all_checks_pass') is True:target=roles['report']
   if target:drift(target);parse_fired[0]=True
  return value
 m.json=SimpleNamespace(loads=loads,dumps=json.dumps)
 m.snapshot=lambda root:{'internal_free_bytes':0 if post[0] and name=='post_resource' else 1<<50,'ssd_free_bytes':1<<50}
 terminal=module(S/'terminal_commit_common_v2.py','_mainterminal_'+name);original_save=terminal.save_new
 pending=m.SSD.parent/'native_canonical_calibration_pending_v4_8.json';seal=m.SSD.parent/'native_canonical_calibration_terminal_seal_v4_8.json';master=m.SSD.parent/'native_canonical_calibration_execution_receipt_v4_8.json'
 def t_save(path,value):
  if str(path).endswith('.failure.json') and name=='supplemental_failure':failed_writes.append(str(path));raise OSError('OWN_FAILURE_WRITE_INJECTION')
  original_save(path,value)
  if Path(path)==seal:
   check(held[0] and pending.is_file(),name+':pending_and_ownership_at_seal');post[0]=True
   if name.startswith('post_mutate_'):
    role=name.removeprefix('post_mutate_');drift({'plan':pp,'admission':ap,'executor':private,'master':master,'seal':seal}.get(role,roles.get(role)))
   if name in ['post_signal','supplemental_failure']:c.catchable_termination(signal.SIGTERM,None)
   if name=='post_extra_file':(roles['output'].parent/'late.json').write_text('{}')
   if name=='post_sidecar':apple(roles['output'].parent/'._worker_stdout.log')
 terminal.save_new=t_save;m.TerminalCommit=terminal.TerminalCommit;m.print=lambda *a,**k:None
 argv=sys.argv;sys.argv=[str(private),'execute'];handlers={s:signal.getsignal(s) for s in [signal.SIGTERM,signal.SIGINT,signal.SIGHUP]}
 try:
  try:m.main();success=True;error=None
  except BaseException as e:success=False;error=type(e).__name__+': '+str(e)
 finally:sys.argv=argv;c.historical_completion=old_history
 restored=all(signal.getsignal(s)==old for s,old in handlers.items());c.TERMINATION_REQUEST.clear()
 consumer=module(S/'terminal_commit_common_v2.py','_consumer_'+name)
 try:consumer.require_committed(pending,seal,initial,{str(master):sha(master)});accepted=True
 except BaseException:accepted=False
 expected=name in ['healthy','healthy_companions']
 check(success==accepted==expected,name+':controller_and_terminal_expected')
 check(not held[0] and restored,name+':ownership_and_handlers_restored')
 rec={'case':name,'expected_and_observed_success':expected,'controller_success':success,'consumer_accepts':accepted,'error':error,'pending_retained':pending.exists(),'stub_calls':len(calls),'actual_workers':0,'actual_mutex_acquisitions':0,'actual_fits':0,'handlers_restored':restored,'metadata_ownership':lifecycle,'supplemental_failed_writes':failed_writes}
 if master.is_file():
  r=json.loads(master.read_text());rec.update(frozen_worker_count=len(r['frozen_worker_evidence']),capture_count=sum(Path(p).name=='native_capture_receipt.json' for p in r['all_output_sha256']),monitor_count=sum(Path(p).name=='monitor_receipt.json' for p in r['all_output_sha256']),external_arithmetic_report_bound=str(m.SSD/'independent_control_arithmetic.json') in r['all_output_sha256'])
  if expected:check(rec['frozen_worker_count']==17 and rec['capture_count']==16 and rec['monitor_count']==17 and rec['external_arithmetic_report_bound'],name+':exact16plus17th_closure')
 RECORDS.append(rec)

def main():
 check(sha(CODE)==CODE_PIN and sha(PLAN)==PLAN_PIN,'initial_frozen_source_plan')
 check(not FIX.exists() and not OUT.exists(),'fresh_private_namespace');FIX.mkdir()
 names=['healthy','healthy_companions','contradictory_cleanup_errors','contradictory_run_error','contradictory_teardown_error','contradictory_first_exception','contradictory_null_exception']
 for name in names:
  run_case(name)
  check(time.monotonic()-START<180,'bounded_deadline:'+name)
  check(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<128*1024*1024,'bounded_RSS:'+name)
 check(sha(CODE)==CODE_PIN and sha(PLAN)==PLAN_PIN,'final_frozen_source_plan')
 files={str(p):sha(p) for p in FIX.rglob('*') if p.is_file() and not p.is_symlink()}
 sizes=sum(p.stat().st_size for p in FIX.rglob('*') if p.is_file() and not p.is_symlink())
 check(sizes<32*1024*1024,'private32MiB_output_cap')
 out={'schema':'independent_canonical_v4_8_terminal_controls','status':'PASS_CORRECTED_CONSUMER_AND_HEALTHY_FULL_PRIVATE_TERMINAL','cases':RECORDS,'case_count':len(RECORDS),'checks':CHECKS,'check_count':len(CHECKS),'checker_sha256':sha(__file__),'frozen_executor_sha256':CODE_PIN,'frozen_plan_sha256':PLAN_PIN,'elapsed_seconds':time.monotonic()-START,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'private_regular_bytes':sizes,'fixture_regular_file_sha256':files,'fixture_symlinks':{str(p):os.readlink(p) for p in FIX.rglob('*') if p.is_symlink()},'actual_workers':0,'actual_fits':0,'actual_mutex':0,'GWAS_or_reference_body_reads':0,'execution_admission_granted':False,'qualification':'Historical190 and scientific workers are explicit stubs. All five contradictory error/null-key fixtures reject before first native metadata stub is credited; both healthy fixtures complete exact16plus1 private terminal. Actual producer stop guarantees unchanged. Private/POSIX cooperative contract, not hostile writer or abrupt-kill durability.'}
 save(OUT,out);print(json.dumps({'receipt_sha256':sha(OUT),'cases':len(RECORDS),'checks':len(CHECKS),'peak_rss_bytes':out['peak_rss_bytes'],'elapsed_seconds':out['elapsed_seconds']}))
if __name__=='__main__':main()
