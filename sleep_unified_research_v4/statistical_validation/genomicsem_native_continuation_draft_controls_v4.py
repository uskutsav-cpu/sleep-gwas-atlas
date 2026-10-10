#!/usr/bin/env python3
"""Independent tiny controls of an unfrozen draft; no input bodies or workers."""
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
PIN='f7a4dd0116c439f17d5ed4fbda6217ac18b5cc4a816505819d182fe8500ae7c7'
OUT=P/'statistical_validation/genomicsem_native_continuation_draft_controls_receipt_v4.json'
sys.dont_write_bytecode=True
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_3 as common


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    code=SOURCE.read_bytes()
    if hashlib.sha256(code).hexdigest()!=PIN:
        raise RuntimeError('UNFROZEN_DRAFT_CHANGED_BEFORE_REVIEW_FIXTURE')
    ast.parse(code)
    checks=[];observations=[]
    def check(value,label):
        if not value:raise RuntimeError('TINY_FIXTURE_FAILED: '+label)
        checks.append(label)
    with tempfile.TemporaryDirectory(prefix='genomicsem_continuation_draft_',dir=P/'statistical_validation') as d:
        tmp=Path(d);copyfile=tmp/'draft61.py';copyfile.write_bytes(code)
        m=load('_independent_unfrozen61_controls',copyfile)
        for label in ['healthy','admission_revoked_in_lock_wait','plan_changed_in_lock_wait','primary_monitor_write_failure','first_worker_failure_then_late_resource_failure','inconsistent_remaining_receipt']:
            base=tmp/label;base.mkdir();out=base/'continuation';hist=base/'historical.json';planfile=base/'plan.json';admfile=base/'admission.json';monitorfile=base/'monitor.json';reviewmd=base/'review.md';reviewjson=base/'review.json'
            reviewmd.write_text('CONTROL_ONLY\n');reviewjson.write_text('{}\n')
            support=base/'support';native=support/'native/extension_reproduction_v1';native.mkdir(parents=True)
            completed=[dict(job_id='CONTROL_COMPLETED_'+str(i),kind='h2',estimates=1) for i in range(107)]
            remaining=[dict(job_id='CONTROL_REMAINING_'+str(i),kind='rg',estimates=100) for i in range(5)]
            m.write_new(hist,dict(ssd_support_package=str(support)))
            plan=dict(historical_plan_sha256=common.sha(hist),completed_jobs=completed,remaining_jobs=remaining,prior_failed_monitor_sha256='CONTROL_ONLY_NO_HISTORICAL_FILE',guard={})
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
                @staticmethod
                def worker(command,prefix,receipt,*_):
                    dispatched.append(command)
                    if label=='first_worker_failure_then_late_resource_failure':raise RuntimeError('CONTROL_FIRST_WORKER_FAILURE')
                    for j in completed+remaining:
                        r=dict(job=j,returncode=0,scientific_cardinality_gate_pass=True,execution_identity_gate_pass=True,fixture_only=True)
                        if label=='inconsistent_remaining_receipt' and j==remaining[0]:r['returncode']=-15
                        m.write_new(native/(j['job_id']+'.execution_receipt.json'),r)
                    m.write_new(receipt,dict(status='WORKER_COMPLETE_VERIFIED',process_group_teardown={'remaining_group_members':[]},fixture_only=True))
            protected=FakeProtected()
            def mockload(path,_name):
                if Path(path).name=='30_prepare_and_run_ssd_native_campaign.py':return FakeMonitor()
                if Path(path).name=='sensitivity_executor_v4_3_1.py':return protected
                if Path(path).name=='58_run_finngen_feasibility_stage_v2.py':return SimpleNamespace(InheritedMutexSubprocess=lambda fd:SimpleNamespace(fd=fd))
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
            patches=[patch.object(m,'PLAN',planfile),patch.object(m,'HISTORICAL',hist),patch.object(m,'OUT',out),patch.object(m,'NEW_MONITOR',monitorfile),patch.object(m,'verify_completed',lambda *_:{}),patch.object(m,'source_bindings',lambda *_:{}),patch.object(m,'preserve_partial',preserve),patch.object(m,'load',mockload),patch.object(m,'write_new',writer),patch.object(common,'exclusive_heavy_lock',isolated_lock),patch('subprocess.Popen',side_effect=RuntimeError('REAL_WORKERS_FORBIDDEN_IN_REVIEW'))]
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
                check(len(dispatched)==len(preserved)==1 and error and result['stop_reason'].endswith('TERMINAL_PLAN_OR_ADMISSION_CHANGED'),'witness_'+label+'_acts_before_terminal_rejection')
            elif label=='primary_monitor_write_failure':check(error.startswith('OSError') and result is None and not addendum.exists(),'witness_primary_terminal_monitor_failure_has_no_fallback_receipt')
            elif label=='first_worker_failure_then_late_resource_failure':
                check('CONTROL_FIRST_WORKER_FAILURE' in result['stop_reason'] and addendum.exists() and 'CONTROL_FIRST_WORKER_FAILURE' not in observed['postpersistence_addendum']['stop_reason'],'witness_addendum_stop_reason_replaces_first_monitor_failure')
            else:check(error is None and result['stop_reason'] is None,'witness_inconsistent_new_remaining_receipt_is_not_validated_at_terminal_gate')
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
    record=dict(status='UNFROZEN_DRAFT_REVIEW_CONTROL_FLOW_FINDINGS_NO_EXECUTION_ADMISSION',draft_sha256=PIN,current_source_sha256=common.sha(SOURCE),draft_may_have_been_corrected_after_snapshot=common.sha(SOURCE)!=PIN,checks=checks,observations=observations,scope='Entire admission/scientific/native worker/real input gates mocked only in scheduler witnesses; tiny archive files and isolated real flock; no actual source or estimator validation',real_workers=0,raw_GWAS_body_reads=0,processed_GWAS_body_reads=0,estimator_calls=0,production_mutex_acquired=False,peak_RSS_bytes=peak,verifier_sha256=common.sha(__file__))
    with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(status=record['status'],draft_sha256=PIN,checks=len(checks),peak_RSS_bytes=peak,receipt_sha256=common.sha(OUT),real_workers=0),indent=2))


if __name__=='__main__':main()
