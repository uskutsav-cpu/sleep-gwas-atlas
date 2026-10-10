#!/usr/bin/env python3
"""Pure software final-gate witness: no subprocess, GWAS, merge or fit."""
import importlib.util
import json
from pathlib import Path
import tempfile

PACKAGE=Path(__file__).resolve().parents[1]
SOURCE=PACKAGE/'scripts/sensitivity_executor_v4_3_1.py'
spec=importlib.util.spec_from_file_location('independent_sensitivity_worker_fixture',SOURCE)
executor=importlib.util.module_from_spec(spec);spec.loader.exec_module(executor)
original_sha=executor.sha
before=original_sha(SOURCE)
queried=[];resource_expired=[]

class Proc:
    pid=123456789
    returncode=0
    def poll(self):return 0

class Monitor:
    @staticmethod
    def terminate_owned(proc):return {'initial_group_members':[],'signals':[],'remaining_group_members':[]}

def limits(*args):
    queried.append(bool(resource_expired))
    return {'software_control':True,'elapsed_stage_seconds':7201 if resource_expired else 0},'STAGE_DEADLINE' if resource_expired else None

with tempfile.TemporaryDirectory(prefix='independent_sensitivity_hash_fault_',dir='/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/tmp') as folder:
    root=Path(folder);executor.SSD=root
    for child in ('logs_v4','receipts_v4','fits_v4'):(root/child).mkdir()
    plan=root/'mock.json';plan.write_text('{}\n');plan_sha=original_sha(plan)
    executor.subprocess.Popen=lambda *a,**k:Proc()
    executor.os.getpgid=lambda pid:pid
    executor.limits=limits
    def sha(path):
        if Path(path).name.endswith('.stdout.log'):resource_expired.append(True)
        return original_sha(path)
    executor.sha=sha
    receipt=root/'receipts_v4/mock.worker.json'
    try:
        executor.worker(['SOFTWARE_FIXTURE_ONLY'],root/'fits_v4/mock',receipt,{'guard':{'poll_seconds':2}},plan,plan_sha,0,Monitor,[True,[]])
    except RuntimeError as error:
        assert 'WORKER_STOPPED_PRESERVED' in str(error)
    else:
        raise AssertionError('POST_HASH_FAULT_ADMITTED_SUCCESS')
    result=json.loads(receipt.read_text())
    assert result['status']=='WORKER_FAILED_PRESERVED' and resource_expired and queried==[False,False,True]
    state,reason=limits()
    assert reason=='STAGE_DEADLINE'
assert original_sha(SOURCE)==before
proof={'status':'PASS_REPAIRED_POST_HASH_FINAL_RESOURCE_GATE','executor_sha256':before,
       'checker_sha256':original_sha(__file__),'worker_receipt_status':result['status'],
       'resource_queries_during_worker':len(queried)-1,'post_hash_violation_if_queried':reason,
       'real_subprocesses_or_estimators_or_audits':0,'GWAS_reads':0,'scope':'Software fixture changes the resource oracle during stdout hashing; no real resource exhaustion or source mutation.'}
out=PACKAGE/'reviews/independent_sensitivity_post_hash_repair_receipt_v4_3_1.json'
with out.open('x') as f:json.dump(proof,f,indent=2);f.write('\n')
print(json.dumps(proof))
