#!/usr/bin/env python3
"""Independent frozen61 metadata and tiny controls; no input bodies or workers."""
import ast
from contextlib import contextmanager
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1]
S=P/'scripts'
SOURCE=S/'61_native_checkpoint_continuation.py'
PIN='53ba7ad9558a72dc1d1e9d9c7ccc357ed211a128bed6b2a898e406f5719b7a00'
PLAN=P/'manifests/native_extension_checkpoint_continuation_v4_2.json'
PLAN_PIN='6f21b9192780cd0490f722439c2ee7c7b62c37b7f0cfe96aa4515558f1483a75'
OUT=P/'statistical_validation/genomicsem_native_continuation_controls_receipt_v4_2.json'
sys.dont_write_bytecode=True
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_3 as common


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    code=SOURCE.read_bytes()
    if hashlib.sha256(code).hexdigest()!=PIN:
        raise RuntimeError('FROZEN61_CHANGED_BEFORE_REVIEW_FIXTURE')
    ast.parse(code)
    checks=[];observations=[]
    def check(value,label):
        if not value:raise RuntimeError('TINY_FIXTURE_FAILED: '+label)
        checks.append(label)
    check(common.sha(PLAN)==PLAN_PIN,'exact_frozen_continuation_plan_SHA')
    frozen=json.loads(PLAN.read_text());historical=json.loads(Path(frozen['historical_plan']).read_text())
    for p,h in frozen['dependencies_sha256'].items():check(common.sha(p)==h,'frozen_dependency_before '+p)
    check(frozen['original190_jobs']==historical['jobs'] and len(historical['jobs'])==190,'unchanged_original190_job_dictionaries')
    extension=[j for j in historical['jobs'] if j['stage']=='extension']
    check(frozen['completed_jobs']+frozen['remaining_jobs']==extension,'original112_extension_order_and_membership')
    check(len(frozen['completed_jobs'])==107 and sum(j['estimates'] for j in frozen['completed_jobs'] if j['kind']=='rg')==700,'exact107_completed_plus700_rg_checkpoint_scope')
    check([j['job_id'] for j in frozen['remaining_jobs']]==['extension_rg_snoring','extension_rg_sleep_apnea','extension_rg_sleep_efficiency','extension_rg_accel_sleep_duration','extension_rg_sleep_timing'],'exact5_original_remaining_batches')
    for p,h in frozen['completed_receipt_sha256'].items():check(common.sha(p)==h,'completed_receipt_metadata_before '+p)
    for r in frozen['interrupted_artifacts']:
        p=Path(r['original_path']);q=Path(r['preserved_path'])
        check(p.is_file() and not p.is_symlink() and p.stat().st_size==r['bytes'] and common.sha(p)==r['sha256'],'interrupted_native_artifact_hash_before '+p.name)
        check(p.parent==Path(historical['ssd_support_package'])/'native/extension_reproduction_v1' and q.parent==Path(frozen['namespace'])/'interrupted_snoring','exact_owned_source_destination_scope '+p.name)
    sidecars=sorted(Path(frozen['interrupted_artifacts'][0]['original_path']).parent.glob('._extension_rg_snoring*'))
    sidecar_inventory=[dict(path=str(p),bytes=p.stat().st_size,sha256=common.sha(p),regular_file=p.is_file(),symlink=p.is_symlink()) for p in sidecars]
    check(len(frozen['interrupted_artifacts'])==284,'284_biological_partial_artifact_map')
    check(len(sidecar_inventory)==284,'284_transport_sidecars_observed_separate_from_biological_map')
    output_verifier=load('_independent_real_baseline_output_guard',S/'sensitivity_executor_v4_3_1.py')
    with tempfile.TemporaryDirectory(prefix='genomicsem_continuation_draft_',dir=P/'statistical_validation') as d:
        tmp=Path(d);copyfile=tmp/'draft61.py';copyfile.write_bytes(code)
        m=load('_independent_unfrozen61_controls',copyfile)
        for label in ['healthy','admission_revoked_in_lock_wait','plan_changed_in_lock_wait','primary_monitor_write_failure','first_worker_failure_then_late_resource_failure','inconsistent_remaining_receipt','admission_revoked_after_primary_monitor','floor_crossed_after_primary_monitor']:
            base=tmp/label;base.mkdir();out=base/'continuation';hist=base/'historical.json';planfile=base/'plan.json';admfile=base/'admission.json';monitorfile=base/'monitor.json';reviewmd=base/'review.md';reviewjson=base/'review.json'
            reviewmd.write_text('CONTROL_ONLY\n');reviewjson.write_text('{}\n')
            support=base/'support';native=support/'native/extension_reproduction_v1';native.mkdir(parents=True)
            completed=[dict(job_id='CONTROL_H2_'+str(i),kind='h2',estimates=1,inputs=['CONTROL_INPUT'],options=[]) for i in range(100)]+[dict(job_id='CONTROL_RG_'+str(i),kind='rg',estimates=100,inputs=['CONTROL_INPUT','CONTROL_OTHER_INPUT'],options=[]) for i in range(7)]
            remaining=[dict(job_id='CONTROL_REMAINING_'+str(i),kind='rg',estimates=100,inputs=['CONTROL_INPUT','CONTROL_OTHER_INPUT'],options=[]) for i in range(5)]
            h=dict(ssd_support_package=str(support),inputs_verified=[dict(path=p,actual_sha256='CONTROL_SOURCE_HASH',match=True) for p in ['CONTROL_INPUT','CONTROL_OTHER_INPUT']])
            m.write_new(hist,h)
            expected_deps={'CONTROL_CODE':'CONTROL_CODE_HASH'}
            def write_job(j,corrupt=False):
                prefix=native/j['job_id'];full=Path(str(prefix)+'.full_precision.json');full.write_text('{"CONTROL_ONLY_NO_ESTIMATE":true}\n')
                command=['CONTROL_PY','-u',str(support/'scripts/native_ldsc_capture.py'),'--ldsc-dir',str(m.ROOT.parent/'ldsc-code'),'--'+j['kind'],','.join(j['inputs']),'--ref-ld-chr','CONTROL_REF/','--w-ld-chr','CONTROL_REF/','--print-delete-vals','--out',str(prefix)]+j['options']
                identities={p:'CONTROL_SOURCE_HASH' for p in j['inputs']}
                r=dict(job=j,command=command,returncode=-15 if corrupt else 0,scientific_cardinality_gate_pass=True,execution_identity_gate_pass=True,input_sha256=identities,input_sha256_after=identities,dependency_sha256_before=expected_deps,dependency_sha256_after=expected_deps,all_output_sha256={str(full):common.sha(full)},output_sha256=common.sha(full),fixture_only=True)
                m.write_new(native/(j['job_id']+'.execution_receipt.json'),r)
            for j in completed:write_job(j)
            plan=dict(historical_plan_sha256=common.sha(hist),completed_jobs=completed,remaining_jobs=remaining,prior_failed_monitor_sha256='CONTROL_ONLY_NO_HISTORICAL_FILE',guard={})
            plan['completed_receipt_sha256']={str(native/(j['job_id']+'.execution_receipt.json')):common.sha(native/(j['job_id']+'.execution_receipt.json')) for j in completed}
            m.write_new(planfile,plan);plan_sha=common.sha(planfile)
            adm=dict(execution_admitted=True,plan_sha256=plan_sha,executor_sha256=PIN,independent_review_sha256={str(reviewmd):common.sha(reviewmd),str(reviewjson):common.sha(reviewjson)})
            m.write_new(admfile,adm);adm_sha=common.sha(admfile)
            dispatched=[];preserved=[];late=[False]
            class FakeMonitor:
                @staticmethod
                def snapshot():return {'internal_free_bytes':0 if late[0] else 9<<30,'ssd_free_bytes':9<<30}
                @staticmethod
                def final_limits(state,*_):return 'INTERNAL_FULL_NATIVE_FLOOR_REACHED' if state['internal_free_bytes']<3<<30 else None
            class FakeProtected:
                TERMINATION_REQUEST=[]
                @staticmethod
                def stage_resource_gate(plan,started,where):
                    state,reason=protected.limits(plan,started,0)
                    if reason:raise RuntimeError('FINAL_RESOURCE_GATE_'+where+': '+reason)
                    return state
                @staticmethod
                def await_owned_cleanup(*_):return []
                baseline_output_binding_gate=staticmethod(output_verifier.baseline_output_binding_gate)
                @staticmethod
                def worker(command,prefix,receipt,*_):
                    dispatched.append(command)
                    if label=='first_worker_failure_then_late_resource_failure':raise RuntimeError('CONTROL_FIRST_WORKER_FAILURE')
                    for j in remaining:write_job(j,corrupt=label=='inconsistent_remaining_receipt' and j==remaining[0])
                    m.write_new(receipt,dict(status='WORKER_COMPLETE_VERIFIED',process_group_teardown={'remaining_group_members':[]},fixture_only=True))
            protected=FakeProtected()
            def mockload(path,_name):
                if Path(path).name=='30_prepare_and_run_ssd_native_campaign.py':return FakeMonitor()
                if Path(path).name=='sensitivity_executor_v4_3_1.py':return protected
                if Path(path).name=='58_run_finngen_feasibility_stage_v2.py':return SimpleNamespace(InheritedMutexSubprocess=lambda fd:SimpleNamespace(fd=fd))
                if Path(path).name=='04_native_reproduction_runner.py':return SimpleNamespace(PYTHON=Path('CONTROL_PY'),REF=Path('CONTROL_REF'))
                raise RuntimeError('UNEXPECTED_MOCK_IMPORT')
            @contextmanager
            def isolated_lock(*,before_release=None):
                fd=os.open(base/'isolated_lock',os.O_RDWR|os.O_CREAT,0o600);fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
                try:
                    if label=='admission_revoked_in_lock_wait':
                        changed=copy.deepcopy(adm);changed['execution_admitted']=False;admfile.write_text(json.dumps(changed))
                    if label=='plan_changed_in_lock_wait':planfile.write_text('{"CONTROL_CHANGED_PLAN":true}\n')
                    yield fd
                finally:
                    if before_release:before_release()
                    os.close(fd)
            def preserve(_plan):preserved.append(True);return {'CONTROL_ONLY_PRESERVATION':'CONTROL_SHA'}
            original_write=m.write_new
            def writer(path,value):
                if Path(path)==monitorfile and label=='primary_monitor_write_failure':raise OSError('CONTROL_PRIMARY_MONITOR_STORAGE_FAILURE')
                original_write(path,value)
                if Path(path)==monitorfile and label=='first_worker_failure_then_late_resource_failure':late[0]=True
                if Path(path)==monitorfile and label=='floor_crossed_after_primary_monitor':late[0]=True
                if Path(path)==monitorfile and label=='admission_revoked_after_primary_monitor':
                    changed=copy.deepcopy(adm);changed['execution_admitted']=False;admfile.write_text(json.dumps(changed))
            patches=[patch.object(m,'PLAN',planfile),patch.object(m,'HISTORICAL',hist),patch.object(m,'OUT',out),patch.object(m,'NEW_MONITOR',monitorfile),patch.object(m,'source_bindings',lambda *_:h),patch.object(m,'preserve_partial',preserve),patch.object(m,'load',mockload),patch.object(m,'write_new',writer),patch.object(common,'expected_historical_dependencies',lambda *_:expected_deps),patch.object(common,'exclusive_heavy_lock',isolated_lock),patch('subprocess.Popen',side_effect=RuntimeError('REAL_WORKERS_FORBIDDEN_IN_REVIEW'))]
            entered=[];error=None
            common.TERMINATION_REQUEST.clear()
            try:
                for p in patches:p.__enter__();entered.append(p)
                try:m.execute(plan_sha,admfile,adm_sha)
                except BaseException as e:error=type(e).__name__+': '+str(e)
            finally:
                for p in reversed(entered):p.__exit__(None,None,None)
                common.TERMINATION_REQUEST.clear()
            result=json.loads(monitorfile.read_text()) if monitorfile.exists() else None
            addendum=Path(str(monitorfile)+'.failure.json')
            observed=dict(case=label,real_workers=0,mocked_native_invoker_dispatches=len(dispatched),mock_preservation_calls=len(preserved),exception=error,monitor=result,postpersistence_addendum=json.loads(addendum.read_text()) if addendum.exists() else None)
            observations.append(observed)
            if label=='healthy':check(error is None and result['returncode']==0 and result['stop_reason'] is None,'healthy_mocked_scheduler_pass')
            elif label in ['admission_revoked_in_lock_wait','plan_changed_in_lock_wait']:
                check(len(dispatched)==len(preserved)==0 and error and 'PLAN_OR_ROOT_ADMISSION_CHANGED_BEFORE_ACTION' in result['stop_reason'],'corrected_'+label+'_blocks_before_preservation_and_dispatch')
            elif label=='primary_monitor_write_failure':check(error.startswith('SystemExit') and result is None and addendum.exists() and 'PRIMARY_MONITOR_SAVE_FAILURE' in observed['postpersistence_addendum']['stop_reason'],'corrected_primary_terminal_monitor_failure_preserves_fallback_receipt')
            elif label=='first_worker_failure_then_late_resource_failure':
                check('CONTROL_FIRST_WORKER_FAILURE' in result['stop_reason'] and addendum.exists() and observed['postpersistence_addendum']['stop_reason']==result['stop_reason'] and observed['postpersistence_addendum']['additional_failures'],'corrected_addendum_preserves_first_failure_and_records_late_failure')
            elif label=='inconsistent_remaining_receipt':check(error and 'CHECKPOINT_JOB_COMMAND_OR_SUCCESS_DIFFERS' in result['stop_reason'],'corrected_real_verify_completed_rejects_bad_new_receipt')
            else:check(error and result['returncode']==0 and result['stop_reason'] is None and observed['postpersistence_addendum']['stop_reason'].startswith('POST_PERSISTENCE_GATE:'),'postpersistence_failure_preserves_provisional_primary_plus_invalidating_addendum_'+label)
        # Genuine tiny same-volume file moves only; no actual native result moved.
        archive=tmp/'tiny_archive_success';archive.mkdir();source=tmp/'tiny_sources';source.mkdir()
        artifacts=[]
        for i in range(3):
            p=source/('CONTROL_ONLY_'+str(i));p.write_bytes(('CONTROL_ONLY_CONTENT_'+str(i)).encode())
            artifacts.append(dict(original_path=str(p),preserved_path=str(archive/'interrupted_snoring'/p.name),bytes=p.stat().st_size,sha256=common.sha(p)))
        with patch.object(m,'OUT',archive):proof=m.preserve_partial({'interrupted_artifacts':artifacts})
        check(all(not Path(r['original_path']).exists() and common.sha(r['preserved_path'])==r['sha256'] for r in artifacts),'tiny_same_volume_preservation_identity_no_deletion')
        events=[json.loads(x) for x in (archive/'interrupted_snoring_preservation.jsonl').read_text().splitlines()]
        check([x['event'] for x in events]==['PRESERVE_INTENT','PRESERVED_READBACK_VERIFIED']*3 and all(common.sha(p)==h for p,h in proof.items()),'tiny_fsynced_intent_commit_journal_and_receipt_hashes')
        # Local subprocess proxy forwards only inherited FD without launching.
        proxy=load('_independent_local58_proxy',S/'58_run_finngen_feasibility_stage_v2.py')
        captured=[]
        with patch.object(proxy.subprocess,'Popen',lambda *a,**k:captured.append(k) or 'CONTROL_NO_PROCESS'):
            check(proxy.InheritedMutexSubprocess(999).Popen(['CONTROL_ONLY'])=='CONTROL_NO_PROCESS' and captured[0]['pass_fds']==(999,),'local_proxy_passes_fd_no_global_subprocess_replacement')
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':peak*=1024
    check(peak<128*(1<<20),'stdlib_review_peak_RSS_below128MiB')
    for p,h in frozen['dependencies_sha256'].items():check(common.sha(p)==h,'frozen_dependency_after '+p)
    for p,h in frozen['completed_receipt_sha256'].items():check(common.sha(p)==h,'completed_receipt_metadata_after '+p)
    for r in frozen['interrupted_artifacts']:check(common.sha(r['original_path'])==r['sha256'],'interrupted_native_artifact_hash_after '+Path(r['original_path']).name)
    for r in sidecar_inventory:check(common.sha(r['path'])==r['sha256'],'transport_sidecar_hash_after '+Path(r['path']).name)
    check(common.sha(PLAN)==PLAN_PIN and common.sha(SOURCE)==PIN,'frozen_plan_and_code_unchanged_after_independent_controls')
    record=dict(status='FROZEN61_CORRECTED_TERMINAL_CONTROLS_PASS_TRANSPORT_SCOPE_REQUIRES_SEPARATE_PROOF',executor_sha256=PIN,plan_sha256=PLAN_PIN,checks=checks,check_count=len(checks),observations=observations,transport_sidecar_inventory=sidecar_inventory,transport_sidecars_in_frozen_plan=False,source_dependency_sha256=frozen['dependencies_sha256'],successful_receipt_metadata_sha256=frozen['completed_receipt_sha256'],interrupted_biological_artifact_map=frozen['interrupted_artifacts'],scope='Admission/scientific/native worker/real input gates mocked in scheduler witnesses; actual verify_completed arithmetic-free receipt guard with synthetic112receipts and real baseline output guard; tiny archive files and isolated real flock; no actual source or estimator validation',real_workers=0,raw_GWAS_body_reads=0,processed_GWAS_body_reads=0,estimator_calls=0,production_mutex_acquired=False,peak_RSS_bytes=peak,verifier_sha256=common.sha(__file__))
    with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(status=record['status'],executor_sha256=PIN,checks=len(checks),peak_RSS_bytes=peak,receipt_sha256=common.sha(OUT),real_workers=0),indent=2))


if __name__=='__main__':main()
