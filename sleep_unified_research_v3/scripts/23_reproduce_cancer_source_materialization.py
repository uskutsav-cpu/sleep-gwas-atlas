#!/usr/bin/env python3
"""Replay unchanged source-specific native helpers into new SSD outputs."""
import argparse
import csv
from datetime import datetime,timezone
import gzip
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
PYTHON=SSD/'Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env/bin/python'
DEST=SSD/'sleep-unified-research-v1/recovery-2026-10-09/cancer_source_materialization_v3'
SCRIPTS={'breast_cancer':'15_materialize_bcac_breast.py','prostate_cancer':'20_materialize_practical_prostate.py','ovarian_cancer':'16_materialize_phelan_ovarian.py'}
BEGAN=None
OWNED_CHILD_STOPPED=None

def guard():
    if shutil.disk_usage(ROOT).free<128*1024**2:raise RuntimeError('INTERNAL_EMERGENCY_FLOOR')
    if shutil.disk_usage(DEST.parent).free<5*1024**3:raise RuntimeError('SSD_RUNTIME_RESERVE')
    if BEGAN is not None and time.monotonic()-BEGAN>3600:raise RuntimeError('ENTIRE_STAGE_TIME_LIMIT')

def sha(path):
    h=hashlib.sha256();checked=time.monotonic()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
            if time.monotonic()-checked>1:guard();checked=time.monotonic()
    return h.hexdigest()

def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')

def decoded(path):
    h=hashlib.sha256(); n=0;lines=0;checked=time.monotonic()
    with gzip.open(path,'rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b);n+=len(b);lines+=b.count(b'\n')
            if time.monotonic()-checked>1:
                guard()
                checked=time.monotonic()
    return {'decompressed_sha256':h.hexdigest(),'decompressed_bytes':n,'newline_count_including_header':lines,'full_gzip_crc_pass':True}

def main(trait):
    global BEGAN,OWNED_CHILD_STOPPED
    BEGAN=time.monotonic();guard()
    if shutil.disk_usage(ROOT).free<256*1024**2:raise RuntimeError('ROW_STREAMING_STAGE_LAUNCH_FLOOR')
    acquisition=PACKAGE/'logs'/(trait+'_acquisition_receipt_v3.json')
    receipt=json.loads(acquisition.read_text())
    if receipt['trait_id']!=trait or receipt['status']!='EXACT_HISTORICAL_SOURCE_REACQUIRED':raise RuntimeError('EXACT_ACQUISITION_GATE_FAILED')
    source=Path(receipt['path'])
    registry=ROOT/'config/public_gwas_sources.tsv'
    with registry.open(newline='') as f:registered=[r for r in csv.DictReader(f,delimiter='\t') if r['trait_ids']==trait]
    if len(registered)!=1:raise RuntimeError('UNIQUE_REGISTERED_SOURCE_GATE_FAILED')
    registered=registered[0]
    if receipt['expected_sha256']!=registered['archive_sha256'] or receipt['expected_bytes']!=int(registered['archive_bytes']) or receipt['url']!=registered['download_url']:
        raise RuntimeError('ACQUISITION_RECEIPT_REGISTRY_IDENTITY_GATE_FAILED')
    status_table=ROOT/'sleep_unified_research_v2/tables/CURRENT_CORE_SOURCE_STATUS.tsv'
    with status_table.open(newline='') as f:old=[r for r in csv.DictReader(f,delimiter='\t') if r['kind']=='core_raw' and r['trait_id']==trait]
    if len(old)!=1:raise RuntimeError('HISTORICAL_RAW_CARDINALITY_GATE_FAILED')
    old=old[0]; historical=Path(old['current_path']); script=ROOT/'scripts'/SCRIPTS[trait]
    blobs={}
    for p in [script,ROOT/'scripts/panel_guard.py',ROOT/'scripts/00_validate_panel.py']:
        rel=str(p.relative_to(ROOT)); current=subprocess.check_output(['git','rev-parse','HEAD:'+rel],cwd=ROOT,text=True).strip()
        original=subprocess.check_output(['git','rev-parse','659d01cf:'+rel],cwd=ROOT,text=True).strip()
        disk=subprocess.check_output(['git','hash-object',str(p)],cwd=ROOT,text=True).strip()
        if current!=original or disk!=original:raise RuntimeError('UNCHANGED_ORIGINAL_HELPER_GATE_FAILED '+rel)
        blobs[rel]=original
    paths={'original_container':source,'historical_raw':historical,'native_script':script,
           'panel_guard':ROOT/'scripts/panel_guard.py','panel_validator':ROOT/'scripts/00_validate_panel.py',
           'panel_manifest':ROOT/'config/analysis_panel.tsv','panel_lock':ROOT/'config/analysis_panel.lock.json',
           'native_python_binary':PYTHON.resolve(),'acquisition_receipt':acquisition,'source_status_table':status_table,'source_registry':registry,'runner':Path(__file__)}
    before={k:sha(v) for k,v in paths.items()}
    if before['original_container']!=receipt['expected_sha256'] or source.stat().st_size!=receipt['expected_bytes'] or before['historical_raw']!=old['expected_sha256']:
        raise RuntimeError('REGISTERED_ORIGINAL_INPUT_HASH_GATE_FAILED')
    outdir=DEST/trait
    if outdir.exists():raise RuntimeError('PRIOR_MATERIALIZATION_PRESERVED')
    if shutil.disk_usage(ROOT).free<256*1024**2 or shutil.disk_usage(DEST.parent).free<7*1024**3:
        raise RuntimeError('LIGHTWEIGHT_MATERIALIZATION_RESOURCE_GATE_FAILED')
    outdir.mkdir(parents=True);(outdir/'tmp').mkdir();output=outdir/(trait+'.txt.gz')
    command=[str(PYTHON),str(script),'--source',str(source),'--out',str(output)]
    if trait=='ovarian_cancer':command+=['--expected-sha256',receipt['expected_sha256']]
    plan={'frozen_utc':datetime.now(timezone.utc).isoformat(),'trait_id':trait,'command':command,
          'bound_paths':{k:str(v) for k,v in paths.items()},'inputs_sha256_before':before,
          'unchanged_original_main_helper_Git_blobs':blobs,'workers':1,'network_bytes':0,
          'internal_minimum_before_bytes':256*1024**2,'internal_free_before_bytes':shutil.disk_usage(ROOT).free,
          'internal_emergency_floor_bytes':128*1024**2,'SSD_minimum_before_bytes':7*1024**3,
          'SSD_runtime_floor_bytes':5*1024**3,'observed_worker_RSS_stop_bytes':128*1024**2,
          'output_final_and_partial_limit_bytes':2*1024**3,'maximum_entire_stage_seconds':3600,
          'resource_basis':'Unchanged stdlib helpers retain one source row and counters; no HM3 dictionary, pandas or estimator loaded',
          'source_id_from_frozen_registry':registered['source_id'],
          'comparison':'Exact historical decompressed bytes and line count; gzip bytes may differ from temporary filename/time metadata',
          'scope':'source-specific raw materialization only; no harmonization, munging, h2 or rg',
          'full_native_LDSC_3GiB_internal_guard_unchanged':True,'original_inputs_or_v1_v2_changed':False,
          'native_python_version':subprocess.check_output([str(PYTHON),'--version'],text=True).strip()}
    plan_path=PACKAGE/'manifests'/(trait+'_materialization_plan_v3.json');save(plan_path,plan);plan_sha=sha(plan_path)
    print(json.dumps({'status':'PLAN_FROZEN','trait':trait,'plan_sha256':plan_sha}),flush=True)
    env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(outdir/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    log_path=PACKAGE/'logs'/(trait+'_native_materialization_v3.log');began=time.monotonic();peak=0;stopped=None;escalated=False
    def retained():
        total=0
        for p in outdir.rglob('*'):
            try:
                if p.is_file():total+=p.stat().st_size
            except FileNotFoundError:pass
        return total
    proc=None
    try:
        with log_path.open('x') as log:
            proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
            while proc.poll() is None:
                guard()
                rss=subprocess.run(['ps','-p',str(proc.pid),'-o','rss='],capture_output=True,text=True).stdout.strip()
                if rss:peak=max(peak,int(rss)*1024)
                if peak>128*1024**2:stopped='WORKER_RSS_LIMIT'
                elif retained()>2*1024**3:stopped='OUTPUT_SIZE_LIMIT'
                if stopped:break
                time.sleep(2)
    finally:
        if proc is not None and proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=30)
            except subprocess.TimeoutExpired:escalated=True;proc.kill();proc.wait()
        OWNED_CHILD_STOPPED=proc is None or proc.poll() is not None
    if retained()>2*1024**3:stopped='OUTPUT_SIZE_LIMIT_AFTER_EXIT'
    elif shutil.disk_usage(ROOT).free<128*1024**2 or shutil.disk_usage(outdir).free<5*1024**3:stopped='RESOURCE_FLOOR_AFTER_EXIT'
    elif time.monotonic()-BEGAN>3600:stopped='ENTIRE_STAGE_TIME_LIMIT_AFTER_EXIT'
    old_content=new_content=None
    if proc.returncode==0 and not stopped and output.exists():old_content=decoded(historical);new_content=decoded(output)
    after={k:sha(v) for k,v in paths.items()}
    guard()
    success=(proc.returncode==0 and not stopped and before==after and sha(plan_path)==plan_sha and old_content is not None and old_content==new_content)
    output_sha=sha(output) if output.exists() else None
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'trait_id':trait,
            'status':'NATIVE_RAW_MATERIALIZATION_EXACT_CONTENT' if success else 'FAILED_OR_CONTENT_DIFFERENCE_PRESERVED',
            'exit_code':proc.returncode,'stop_reason':stopped,'worker_peak_observed_RSS_bytes':peak,
            'termination_escalated':escalated,'elapsed_seconds_including_comparison':time.monotonic()-began,
            'inputs_sha256_before':before,'inputs_sha256_after':after,'input_hashes_unchanged':before==after,
            'plan_sha256':plan_sha,'plan_hash_unchanged':sha(plan_path)==plan_sha,'log_sha256':sha(log_path),
            'historical_decompressed':old_content,'native_decompressed':new_content,
            'gzip_bytes_identical':output_sha==old['expected_sha256'],'output_path':str(output),'output_sha256':output_sha,
            'combined_new_output_bytes':retained(),'h2_or_rg_estimated':False,'full_QC_chain_reproduced':False,
            'original_inputs_modified':False}
    save(PACKAGE/'logs'/(trait+'_materialization_receipt_v3.json'),result)
    print(json.dumps({'trait':trait,'status':result['status'],'gzip_bytes_identical':result['gzip_bytes_identical']}),flush=True)
    if not success:raise RuntimeError('FAILED_OR_CONTENT_DIFFERENCE_PRESERVED')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--trait',required=True,choices=SCRIPTS);trait=parser.parse_args().trait
    try:main(trait)
    except BaseException as error:
        path=PACKAGE/'logs'/(trait+'_materialization_FAILURE_v3.json')
        if not path.exists():save(path,{'completed_utc':datetime.now(timezone.utc).isoformat(),'trait_id':trait,
            'status':'FAILED_NATIVE_STAGE_PRESERVED','error_type':type(error).__name__,'error':str(error),
            'owned_native_child_stopped':OWNED_CHILD_STOPPED,'elapsed_seconds':time.monotonic()-BEGAN if BEGAN is not None else None,
            'original_inputs_modified':False,'downstream_admission':False})
        raise
