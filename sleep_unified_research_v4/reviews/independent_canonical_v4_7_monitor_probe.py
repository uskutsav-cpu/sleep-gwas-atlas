#!/usr/bin/env python3
"""Tiny actual monitor-consumer contradiction probe; no worker, mutex or data."""
import importlib.util,json,sys,hashlib
from pathlib import Path
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_7 as common
spec=importlib.util.spec_from_file_location('_independent_canonical7_monitor',S/'42_prepare_native_canonical_calibration_v4_7.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
ROOT=P.parents[1]/'independent_canonical_v4_7_monitor_probe';OUT=R/'independent_canonical_v4_7_monitor_probe_receipt.json'
if ROOT.exists() or OUT.exists():raise RuntimeError('OWN_FRESH_NAMESPACE_REQUIRED')
ROOT.mkdir();cases=[]
for name in ['healthy','cleanup_errors','run_error','teardown_cleanup_error','first_teardown_exception']:
 d=ROOT/name;out=d/'native';out.mkdir(parents=True);(out/'weighted_capture').mkdir()
 m.PLAN=d/'plan.json';m.ADMISSION=d/'admission.json'
 common.write_intended_new(m.PLAN,{'own_fixture':name},common.intended_sha({'own_fixture':name}));common.write_intended_new(m.ADMISSION,{'own_fixture':name},common.intended_sha({'own_fixture':name}))
 job={'job_id':'OWN_METADATA','output_dir':str(out)};plan={'jobs':[job],'ssd_output_root':str(d)};command=['NO_PROCESS']
 cap=out/'native_capture_receipt.json';common.write_intended_new(cap,{'job':job,'plan_sha256':common.sha(m.PLAN)},common.intended_sha({'job':job,'plan_sha256':common.sha(m.PLAN)}))
 monitor=out/'monitor_receipt.json';record={'job_id':job['job_id'],'plan_sha256':common.sha(m.PLAN),'admission_sha256':common.sha(m.ADMISSION),'command':command,'returncode':0,'stop_reason':None,'process_group_teardown':{'remaining_group_members':[]},'final_guard_violations':[],'owned_group_disappearance_verified':True,'bulk_output_hashes_completed_before_final_guards':True,'output_hash_error':None,'deferred_termination_requests':[],'cleanup_errors':[],'run_error':None,'first_teardown_attempt':{'signal_zero_group_absent':True},'all_output_sha256':{str(cap):common.sha(cap)}}
 if name=='cleanup_errors':record['cleanup_errors']=['METADATA_PROOF_WRITE_FAILED']
 if name=='run_error':record['run_error']={'type':'RuntimeError','repr':'explicit failure'}
 if name=='teardown_cleanup_error':record['process_group_teardown']['cleanup_error']='explicit failure'
 if name=='first_teardown_exception':record['first_teardown_attempt']['exception']={'type':'RuntimeError','repr':'explicit failure'}
 expected=common.intended_sha(record);common.write_intended_new(monitor,record,expected)
 evidence={'workers':{},'all_output_sha256':{}};error=None
 try:m.freeze_worker_evidence(plan,common.sha(m.ADMISSION),job['job_id'],command,out,evidence,job,{'path':str(monitor),'sha256':expected});accepted=True
 except BaseException as e:accepted=False;error=repr(e)
 cases.append({'case':name,'accepted':accepted,'error':error,'monitor_path':str(monitor),'monitor_sha256':expected,'frozen_sha256':evidence['all_output_sha256'],'record':record})
receipt={'schema':'independent_canonical_v4_7_monitor_contradiction_probe','cases':cases,'status':'MEASURED_FOUR_CONTRADICTORY_FAILURE_ACCEPTANCES','actual_workers':0,'actual_fits':0,'actual_mutex':0,'production_body_reads':0,'source_sha256':{str(p):common.sha(p) for p in [S/'42_prepare_native_canonical_calibration_v4_7.py',S/'canonical_calibration_common_v4_7.py']},'qualification':'Deliberately contradictory monitor records; normal current producer sets stop_reason on these faults. No production failure inferred.'}
common.write_intended_new(OUT,receipt,common.intended_sha(receipt));print(json.dumps({'sha256':common.sha(OUT),'cases':[(r['case'],r['accepted']) for r in cases]}))
