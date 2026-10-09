#!/usr/bin/env python3
"""Freeze dependency receipts and enforce a workload-specific resource plan."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent
SSD=Path('/Volumes/Extreme SSD')
OLD=SSD/'Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
PYTHON=OLD/'.ldsc-env/bin/python'
DEST=SSD/'sleep-unified-research-v1/recovery-2026-10-09/mvp_hm3_preprocessing_v3'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''): h.update(b)
    return h.hexdigest()

def save(path,obj):
    with path.open('x') as f: json.dump(obj,f,indent=2); f.write('\n')

def main():
    for name in ['logs','manifests','reviews','tables','tests']: (PACKAGE/name).mkdir(exist_ok=True)
    acquisition=ROOT/'sleep_unified_research_v1/logs/mvp_insomnia_acquisition_receipt_v1.json'
    source_receipt=json.loads(acquisition.read_text())
    if source_receipt['verification_status']!='PASS': raise RuntimeError('SOURCE_ACQUISITION_FAILED')
    worker=PACKAGE/'scripts/19_materialize_mvp_hm3.py'
    paths={'source':Path(source_receipt['path']),
           'hm3_reference':SSD/'sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz',
           'hm3_alleles':OLD/'work/track_b_completion/local_dependency_copies/ref/eur_w_ld_chr/w_hm3.snplist',
           'chain':OLD/'work/track_b_completion/local_dependency_copies/ref/hg38ToHg19.over.chain.gz',
           'historical_47':ROOT/'discovery_extension/scripts/47_stream_replication_sources.py',
           'historical_streaming_helper':ROOT/'discovery_extension/scripts/streaming_io.py',
           'chain_helper':ROOT/'scripts/liftover_chain.py','worker':worker,'runner':Path(__file__),
           'python_binary':PYTHON.resolve(),'protocol':PACKAGE/'FROZEN_MVP_PREPROCESSING_PROTOCOL_v4.md',
           'acquisition_receipt':acquisition,
           'original_schema_receipt':ROOT/'sleep_unified_research_v1/logs/mvp_insomnia_source_schema_v1.json',
           'definition_contract':ROOT/'sleep_unified_research_v2/source_definition_review/provenance_v2_definition_contract.json',
           'protocol_amendment':ROOT/'sleep_unified_research_v2/manifests/protocol_amendment_v3_receipt.json'}
    expected={'source':'46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835',
              'hm3_reference':'e6e4814d99a1eff91875014fe70c61f760182efdc1711ba6e155963cf5aa11f8',
              'hm3_alleles':'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed',
              'chain':'14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1'}
    before={k:sha(v) for k,v in paths.items()}
    if any(before[k]!=v for k,v in expected.items()): raise RuntimeError('PINNED_INPUT_HASH_GATE_FAILED')
    if paths['source'].stat().st_size!=575504178: raise RuntimeError('SOURCE_SIZE_GATE_FAILED')
    if DEST.exists(): raise RuntimeError('PRIOR_DESTINATION_PRESERVED')
    free=shutil.disk_usage(ROOT).free
    if free<512*1024**2 or shutil.disk_usage(SSD).free<5*1024**3: raise RuntimeError('LIGHTWEIGHT_RESOURCE_GATE_FAILED')
    DEST.mkdir(parents=True); (DEST/'tmp').mkdir()
    out=DEST/'mvp_GIA_insomnia.hm3.sumstats.gz'
    plan_path=PACKAGE/'manifests/mvp_preprocessing_execution_plan_v3.json'
    worker_receipt=PACKAGE/'logs/mvp_preprocessing_worker_receipt_v3.json'
    command=[str(PYTHON),str(worker),'--plan',str(plan_path)]
    plan={'frozen_utc':datetime.now(timezone.utc).isoformat(),'bound_paths':{k:str(v) for k,v in paths.items()},
          'inputs_sha256_before':before,'expected_registered_input_hashes':expected,
          'command':command,'python_version':subprocess.check_output([str(PYTHON),'--version'],text=True).strip(),
          'output_path':str(out),'worker_receipt_path':str(worker_receipt),
          'workers':1,'network_bytes':0,'internal_free_before':free,'internal_minimum_before':512*1024**2,
          'internal_emergency_stop_bytes':128*1024**2,'worker_observed_RSS_stop_bytes':768*1024**2,
          'maximum_seconds':3600,'combined_retained_output_limit_bytes':128*1024**2,
          'TMPDIR':str(DEST/'tmp'),'scope':'result-free preprocessing only',
          'full_native_LDSC_3GiB_internal_guard_unchanged':True,'pair_tests_admitted':0}
    save(plan_path,plan)
    frozen_plan_sha256=sha(plan_path)
    print(json.dumps({'status':'PLAN_FROZEN','path':str(plan_path),'sha256':sha(plan_path)}),flush=True)
    env=os.environ.copy(); env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(DEST/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    log_path=PACKAGE/'logs/mvp_preprocessing_native_v3.log'; began=time.time(); peak=0; stopped=None; escalated=False
    output_paths=[out,Path(str(out)+'.partial'),worker_receipt]
    def retained():
        n=0
        for p in output_paths:
            try:n+=p.stat().st_size
            except FileNotFoundError:pass
        return n
    with log_path.open('x') as log:
        proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            rss=subprocess.run(['ps','-p',str(proc.pid),'-o','rss='],capture_output=True,text=True).stdout.strip()
            if rss: peak=max(peak,int(rss)*1024)
            if peak>plan['worker_observed_RSS_stop_bytes']: stopped='WORKER_RSS_LIMIT'
            elif shutil.disk_usage(ROOT).free<plan['internal_emergency_stop_bytes']: stopped='INTERNAL_EMERGENCY_FLOOR'
            elif time.time()-began>plan['maximum_seconds']: stopped='TIME_LIMIT'
            elif retained()>plan['combined_retained_output_limit_bytes']: stopped='OUTPUT_SIZE_LIMIT'
            if stopped:
                proc.terminate()
                try:proc.wait(timeout=30)
                except subprocess.TimeoutExpired: escalated=True; proc.kill(); proc.wait()
                break
            time.sleep(2)
    after={k:sha(v) for k,v in paths.items()}
    if retained()>plan['combined_retained_output_limit_bytes']: stopped='OUTPUT_SIZE_LIMIT_AFTER_EXIT'
    elif shutil.disk_usage(ROOT).free<plan['internal_emergency_stop_bytes']: stopped='INTERNAL_EMERGENCY_FLOOR_AFTER_EXIT'
    elif time.time()-began>plan['maximum_seconds']: stopped='TIME_LIMIT_AFTER_EXIT'
    numerical=json.loads(worker_receipt.read_text()) if worker_receipt.exists() else None
    plan_unchanged=sha(plan_path)==frozen_plan_sha256
    success=proc.returncode==0 and not stopped and before==after and plan_unchanged and numerical is not None and numerical['status']=='PREPROCESSING_ONLY_PASS' and numerical['execution_plan_sha256']==frozen_plan_sha256
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'status':'RESULT_FREE_PREPROCESSING_PASS' if success else 'FAILED_PREPROCESSING_PRESERVED',
            'elapsed_seconds':time.time()-began,'exit_code':proc.returncode,'stop_reason':stopped,
            'termination_escalated':escalated,'worker_peak_observed_RSS_bytes':peak,
            'inputs_sha256_before':before,'inputs_sha256_after':after,'input_hashes_unchanged':before==after,
            'execution_plan_sha256':frozen_plan_sha256,'execution_plan_hash_unchanged':plan_unchanged,'worker_receipt_sha256':sha(worker_receipt) if numerical else None,
            'log_sha256':sha(log_path),'output_path':str(out),'output_sha256':sha(out) if out.exists() else None,
            'combined_output_bytes':retained(),'original_inputs_modified':False,'pair_tests_admitted':0,
            'complete_native_reproduction':False,'h2_or_rg_estimated':False}
    save(PACKAGE/'logs/mvp_preprocessing_execution_receipt_v3.json',result)
    print(json.dumps(result),flush=True)
    if not success: raise RuntimeError('FAILED_PREPROCESSING_PRESERVED')

if __name__=='__main__':main()
