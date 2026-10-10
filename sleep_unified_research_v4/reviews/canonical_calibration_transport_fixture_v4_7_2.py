#!/usr/bin/env python3
"""Author metadata-only sanity controls; independent review remains mandatory."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import sys
import types

R=Path(__file__).resolve().parent;P=R.parent;S=P/'scripts';F=R/'canonical_calibration_transport_fixtures_v4_7_2'
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_7 as c
from terminal_commit_common_v2 import TerminalCommit
spec=importlib.util.spec_from_file_location('candidate_canonical7',S/'42_prepare_native_canonical_calibration_v4_7.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
checks=[];cases=[]
def check(value,label):
 checks.append(dict(label=label,pass_=bool(value)))
 if not value:raise AssertionError(label)
def save(p,v):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n');return p
def blob(p,v=b'TINY'):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(v);return p
def apple(p):return blob(p,struct.pack('>II16sH',0x00051607,0x00020000,b'\0'*16,1)+struct.pack('>III',9,38,4)+b'TEST')
def base(name,arithmetic=False):
 root=F/name/'native_canonical_calibration_v4_7';root.mkdir(parents=True)
 out=root/('independent_arithmetic_worker' if arithmetic else 'native');out.mkdir()
 planpath=save(root.parent/'plan.json',dict(fixture=name));admission=save(root.parent/'admission.json',dict(fixture=name))
 m.PLAN=planpath;m.ADMISSION=admission
 job=dict(job_id='NATIVE_TINY',output_dir=str(out));plan=dict(jobs=[job],ssd_output_root=str(root))
 command=['NO_WORKER','METADATA_ONLY'];jobid='INDEPENDENT_CAPTURED_ROW_ARITHMETIC' if arithmetic else job['job_id']
 stdout=blob(out/'worker_stdout.log');paths=[stdout]
 evidence=dict(workers={},all_output_sha256={})
 if not arithmetic:
  cap=save(out/'native_capture_receipt.json',dict(job=job,plan_sha256=c.sha(planpath)));paths.append(cap)
  paths.append(blob(out/'weighted_capture/tiny_array'))
 else:
  plan['jobs']=[];sources={}
  for i in range(16):
   source=save(root/('native'+str(i))/'native_capture_receipt.json',dict(synthetic_only=True,index=i));plan['jobs'].append(dict(job_id='N'+str(i),output_dir=str(source.parent)));sources[str(source)]=c.sha(source)
  evidence['all_output_sha256'].update(sources)
  report=root/'independent_control_arithmetic.json';record=dict(all_checks_pass=True,scope='METHOD_CALIBRATION_IMPLEMENTATION_CHECK_ONLY',allow_41_covariance_outcomes=False,calibrated_biological_p_values_computed=False)
  report_sha=c.intended_sha(record)
  proof=out/'arithmetic_report_intent.json';intent=dict(schema='canonical_arithmetic_report_intent_v4_7',executor_sha256=c.sha(S/'canonical_calibration_arithmetic_v4_7.py'),plan_sha256=c.sha(planpath),admission_sha256=c.sha(admission),report_path=str(report),intended_report_sha256=report_sha,capture_receipt_sha256=sources)
  c.write_intended_new(proof,intent,c.intended_sha(intent));c.write_intended_new(report,record,report_sha)
  binding=dict(plan_sha256=c.sha(planpath),admission_sha256=c.sha(admission),executor_sha256=c.sha(S/'canonical_calibration_arithmetic_v4_7.py'))
  tc=TerminalCommit(out/'arithmetic_pending.json',out/'arithmetic_terminal_seal.json',binding)
  expected_seal=c.terminal_intended_sha(tc,{str(report):report_sha})
  check(tc.commit({str(report):report_sha},identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None),name+'_private_terminal_success')
  check(c.sha(tc.seal)==expected_seal,name+'_seal_intent_before_write_matches')
  paths.extend([report,proof,tc.seal])
 mapping={str(p):c.sha(p) for p in paths}
 monitor=out/'monitor_receipt.json'
 receipt=dict(job_id=jobid,plan_sha256=c.sha(planpath),admission_sha256=c.sha(admission),command=command,
  returncode=0,stop_reason=None,process_group_teardown=dict(remaining_group_members=[]),final_guard_violations=[],owned_group_disappearance_verified=True,
  bulk_output_hashes_completed_before_final_guards=True,output_hash_error=None,deferred_termination_requests=[],all_output_sha256=mapping)
 c.write_intended_new(monitor,receipt,c.intended_sha(receipt))
 return dict(root=root,out=out,job=job,plan=plan,command=command,jobid=jobid,evidence=evidence,monitor=monitor,monitor_identity=dict(path=str(monitor),sha256=c.sha(monitor)),native_job=None if arithmetic else job)
def freeze(f):m.freeze_worker_evidence(f['plan'],c.sha(m.ADMISSION),f['jobid'],f['command'],f['out'],f['evidence'],f['native_job'],f['monitor_identity'])
def reject(name,change):
 f=base(name);change(f)
 try:freeze(f);error=None
 except BaseException as e:error=type(e).__name__+': '+str(e)
 check(error is not None,name+'_reject');cases.append(dict(name=name,error=error,pass_=True))

def main():
 check(not F.exists(),'fresh_metadata_fixture_namespace');F.mkdir()
 sources=[S/(n+'_v4_'+str(v)+'.py') for v in [6,7] for n in ['42_prepare_native_canonical_calibration','canonical_calibration_common','canonical_calibration_capture','canonical_calibration_arithmetic']]
 before={str(p):c.sha(p) for p in sources}
 for p in sources:ast.parse(p.read_text())
 a=ast.parse((S/'canonical_calibration_arithmetic_v4_6.py').read_text());b=ast.parse((S/'canonical_calibration_arithmetic_v4_7.py').read_text())
 f=lambda t:next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='validate')
 check(ast.dump(f(a),include_attributes=False)==ast.dump(f(b),include_attributes=False),'scientific_arithmetic_validate_exact_AST')
 check((S/'canonical_calibration_capture_v4_6.py').read_text().replace('v4_6','v4_7')==(S/'canonical_calibration_capture_v4_7.py').read_text(),'capture_only_import_namespace_changed')
 f=base('valid_new_monitor_companion');apple(f['out']/'._monitor_receipt.json');freeze(f)
 check(str(f['out']/'._monitor_receipt.json') in f['evidence']['all_output_sha256'],'registered_new_monitor_companion_frozen');m.check_frozen_worker_evidence(f['evidence']);cases.append(dict(name='valid_new_monitor_companion',pass_=True))
 f=base('valid_directory_companion');apple(f['out']/'._weighted_capture');freeze(f);m.check_frozen_worker_evidence(f['evidence']);check(str(f['out']/'._weighted_capture') in f['evidence']['all_output_sha256'],'registered_directory_companion_frozen');cases.append(dict(name='valid_directory_companion',pass_=True))
 reject('unknown_companion',lambda f:apple(f['out']/'._unregistered'))
 reject('fake_companion_bytes',lambda f:blob(f['out']/'._monitor_receipt.json',b'NOT_APPLEDOUBLE'))
 reject('companion_symlink',lambda f:(f['out']/'._monitor_receipt.json').symlink_to(f['monitor']))
 reject('unknown_ordinary_file',lambda f:blob(f['out']/'unknown'))
 reject('unknown_directory',lambda f:(f['out']/'unknown_dir').mkdir())
 reject('monitor_before_parse_changed',lambda f:f['monitor'].write_text(f['monitor'].read_text()+' '))
 def parse_swap(f):
  original=json.loads
  def loads(text,*args,**kw):
   r=original(text,*args,**kw)
   if isinstance(r,dict) and r.get('job_id')==f['jobid']:f['monitor'].write_text(f['monitor'].read_text()+' ')
   return r
  prior=m.json;m.json=types.SimpleNamespace(loads=loads)
  try:
   try:freeze(f);err=None
   except BaseException as e:err=repr(e)
  finally:m.json=prior
  check(err is not None,'monitor_after_parse_changed_reject');cases.append(dict(name='monitor_after_parse_changed',error=err,pass_=True))
 parse_swap(base('monitor_after_parse_changed'))
 f=base('companion_after_freeze');comp=apple(f['out']/'._monitor_receipt.json');freeze(f);comp.write_bytes(comp.read_bytes()[:-1]+b'X')
 try:m.check_frozen_worker_evidence(f['evidence']);err=None
 except BaseException as e:err=repr(e)
 check(err is not None,'frozen_companion_drift_reject');cases.append(dict(name='companion_after_freeze',error=err,pass_=True))
 f=base('arithmetic_external_route',True);apple(f['root']/'._independent_control_arithmetic.json');freeze(f);m.check_frozen_worker_evidence(f['evidence'])
 check(str(f['root']/'independent_control_arithmetic.json') in f['evidence']['workers'][f['jobid']]['sha256'],'outside_worker_report_bound_to_arithmetic_worker')
 check(str(f['root']/'._independent_control_arithmetic.json') in f['evidence']['all_output_sha256'],'outside_worker_report_companion_bound');cases.append(dict(name='arithmetic_external_route',pass_=True))
 f=base('arithmetic_deleted_marker_companion',True);apple(f['out']/'._arithmetic_pending.json');freeze(f);m.check_frozen_worker_evidence(f['evidence']);cases.append(dict(name='arithmetic_deleted_marker_companion_exact_registered_committed_marker',pass_=True))
 check(not (f['out']/'arithmetic_pending.json').exists(),'successful_private_pending_absent')
 f=base('epoch_budget');old=f['root'].parent/'native_canonical_calibration_v4_6';blob(old/'data',b'A'*31);blob(f['root']/'data',b'B'*17)
 value=c.canonical_epoch_bytes(f['root']);expected=sum(p.stat().st_size for d in [old,f['root']] for p in d.rglob('*') if p.is_file());check(value==expected,'all_epochs8GiB_meter_includes_old_failed_and_new');cases.append(dict(name='all_retained_epochs_budget',bytes=value,pass_=True))
 # Reproduce actual old failed transport case with read-only metadata/byte hashes.
 realplan=json.loads((P/'manifests/native_canonical_calibration_plan_v4_6.json').read_text());realjob=realplan['jobs'][0];realout=Path(realjob['output_dir']);monitor=realout/'monitor_receipt.json';realmon=json.loads(monitor.read_text())
 m.PLAN=P/'manifests/native_canonical_calibration_plan_v4_6.json';m.ADMISSION=P/'manifests/native_canonical_calibration_admission_v4_6.json'
 ev=dict(workers={},all_output_sha256={});m.freeze_worker_evidence(realplan,c.sha(m.ADMISSION),realjob['job_id'],realmon['command'],realout,ev,realjob,dict(path=str(monitor),sha256=c.sha(monitor)))
 check(str(realout/'._monitor_receipt.json') in ev['all_output_sha256'],'actual_failed_v6_monitor_companion_now_registered_and_frozen');m.check_frozen_worker_evidence(ev)
 cases.append(dict(name='actual_v6_transport_mismatch_readonly_metadata_hash_reproduction',pass_=True,no_new_fit=True,not_full16_completion=True))
 check(all(c.sha(p)==h for p,h in before.items()),'all_old_and_new_sources_unchanged_by_checks')
 receipt=dict(schema='author_canonical7_metadata_fixture_v1',status='AUTHOR_METADATA_SANITY_PASS_NOT_INDEPENDENT_ADMISSION',cases=cases,case_count=len(cases),checks=checks,check_count=len(checks),all_checks_pass=all(x['pass_'] for x in checks),source_sha256=before,no_workers_or_fits_or_mutex=True,no_original_input_or_reference_body_read=True,actual_v6_read_scope='Only already-authorized preserved output byte hashing and monitor/capture metadata; no numeric estimator outcomes inspected or recomputed.')
 c.write_intended_new(R/'canonical_calibration_transport_fixture_receipt_v4_7_2.json',receipt,c.intended_sha(receipt))
 print(json.dumps(dict(status=receipt['status'],cases=len(cases),checks=len(checks))))
if __name__=='__main__':main()
