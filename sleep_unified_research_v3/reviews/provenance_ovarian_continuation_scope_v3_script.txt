#!/usr/bin/env python3
"""Continue the unattempted ovarian object without retrying failed prostate."""
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import time

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent
helper_path=PACKAGE/'scripts/21_recover_remaining_core_containers.py'
approved_plan=json.loads((PACKAGE/'manifests/remaining_core_containers_resource_plan_v3.json').read_text())
if hashlib.sha256(helper_path.read_bytes()).hexdigest()!=approved_plan['script_sha256']:
    raise RuntimeError('APPROVED_HELPER_HASH_REQUIRED_BEFORE_IMPORT')
spec=importlib.util.spec_from_file_location('bounded_container_transfer',helper_path)
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
h.DEST=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/ovarian_only_continuation_v3')

def main():
    version=subprocess.check_output(['curl','--version'],text=True).splitlines()[0]
    if tuple(map(int,version.split()[1].split('.')))<(8,4,0):raise RuntimeError('CURL_STREAMING_CAP_REQUIRED')
    failure=PACKAGE/'logs/remaining_core_containers_FAILURE_v3.json'
    failed_source=PACKAGE/'logs/prostate_cancer_acquisition_receipt_v3.json'
    old_plan=PACKAGE/'manifests/remaining_core_containers_resource_plan_v3.json'
    if json.loads(failure.read_text())['error']!='FULL_HISTORICAL_SHA_FAILED_NO_SUBSTITUTION' or json.loads(failed_source.read_text())['status']!='SHA_MISMATCH_PRESERVED':
        raise RuntimeError('ORIGINAL_FAILURE_MUST_REMAIN_EXPLICIT')
    registry=ROOT/'config/public_gwas_sources.tsv'
    with registry.open(newline='') as f:sources=[r for r in csv.DictReader(f,delimiter='\t') if r['trait_ids']=='ovarian_cancer']
    if len(sources)!=1:raise RuntimeError('UNIQUE_FROZEN_OVARIAN_SOURCE_REQUIRED')
    s=sources[0];n=int(s['archive_bytes']);head_path=PACKAGE/'logs/ovarian_public_HEAD_v3.txt';identity=h.headers(head_path)
    original=json.loads(old_plan.read_text())
    if s not in original['sources'] or h.sha(registry)!=original['source_registry_sha256'] or h.sha(helper_path)!=original['script_sha256'] or h.sha(head_path)!=original['HEAD_sha256']['ovarian_cancer']:
        raise RuntimeError('ORIGINAL_FROZEN_INPUT_BINDING_FAILED')
    if identity.get('status_line','').split()[1:2]!=['200'] or identity.get('content-length')!=str(n) or not identity.get('etag'):raise RuntimeError('HEAD_IDENTITY_GATE_FAILED')
    original_dest=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/remaining_core_containers_v3/ovarian_cancer')
    if original_dest.exists() or h.DEST.exists():raise RuntimeError('PRIOR_OVARIAN_ATTEMPT_PRESERVED_REQUIRES_REVIEW')
    if shutil.disk_usage(ROOT).free<256*1024**2 or shutil.disk_usage(h.DEST.parent).free<2*n+5*1024**3:raise RuntimeError('CONTINUATION_RESOURCE_GATE_FAILED')
    h.DEST.mkdir();(h.DEST/s['trait_ids']/'chunks').mkdir(parents=True)
    plan={'frozen_utc':h.now(),'source':s,'network_bytes':n,'new_retention_bytes_full_and_chunks':2*n,
          'original_planned_total_network_bytes':original['network_bytes'],'unattempted_source_within_original_budget':True,
          'preserved_failed_prostate_bytes_used_for_science':False,'original_failure_sha256':h.sha(failure),
          'failed_source_receipt_sha256':h.sha(failed_source),'original_resource_plan_sha256':h.sha(old_plan),
          'helper_sha256':h.sha(helper_path),'wrapper_sha256':h.sha(__file__),'source_registry_sha256':h.sha(registry),
          'HEAD_sha256':h.sha(head_path),'curl_version':version,'workers':h.WORKERS,'chunk_bytes':h.CHUNK,'max_seconds_per_chunk':600,
          'maximum_total_seconds':5400,'deadline_clock':'monotonic','internal_launch_bytes':256*1024**2,
          'resource_basis':'Existing bounded curl transfers write bodies directly to SSD; observed curl RSS about8MiB each, small Python buffers',
          'internal_emergency_bytes':128*1024**2,'SSD_runtime_floor_with_in_flight_bytes':5*1024**3+h.WORKERS*h.CHUNK,
          'destination':str(h.DEST),'no_retry_or_source_substitution':True,'scientific_outcomes_accessed':False}
    plan_path=PACKAGE/'manifests/ovarian_only_continuation_plan_v3.json';h.save(plan_path,plan)
    print(json.dumps({'status':'OVARIAN_CONTINUATION_FROZEN','network_bytes':n,'plan_sha256':h.sha(plan_path)}),flush=True)
    h.BEGAN=time.monotonic();intervals=[(a,min(a+h.CHUNK,n)-1) for a in range(0,n,h.CHUNK)];results=[]
    with ThreadPoolExecutor(max_workers=h.WORKERS) as pool:
        pending=set();cursor=0
        try:
            while cursor<len(intervals) or pending:
                h.guard()
                while cursor<len(intervals) and len(pending)<h.WORKERS:
                    a,b=intervals[cursor];pending.add(pool.submit(h.fetch,s,identity,a,b));cursor+=1
                done,_=wait(pending,timeout=1,return_when=FIRST_COMPLETED)
                for future in done:
                    pending.remove(future);results.append(future.result());print(json.dumps({'chunks_complete':len(results),'chunks_total':len(intervals)}),flush=True)
        except BaseException:
            h.STOP.set()
            for future in pending:future.cancel()
            raise
    ordered=sorted(results,key=lambda r:r['start'])
    if [(r['start'],r['end']) for r in ordered]!=intervals or sum(r['actual_bytes'] for r in ordered)!=n:raise RuntimeError('RANGE_COVERAGE_FAILED')
    partial=h.DEST/s['trait_ids']/(s['archive_name']+'.assembling');final=h.DEST/s['trait_ids']/s['archive_name']
    import hashlib
    digest=hashlib.sha256()
    with partial.open('xb') as output:
        for r in ordered:
            h.guard()
            if h.sha(r['body_path'])!=r['body_sha256']:raise RuntimeError('CHUNK_HASH_CHANGED')
            with Path(r['body_path']).open('rb') as f:
                checked=time.monotonic()
                for b in iter(lambda:f.read(65536),b''):
                    output.write(b);digest.update(b)
                    if time.monotonic()-checked>1:h.guard();checked=time.monotonic()
    h.guard();passed=partial.stat().st_size==n and digest.hexdigest()==s['archive_sha256']
    if passed:partial.rename(final)
    r={'completed_utc':h.now(),'trait_id':s['trait_ids'],'source_id':s['source_id'],'url':s['download_url'],
       'expected_bytes':n,'actual_bytes':n,'expected_sha256':s['archive_sha256'],'actual_sha256':digest.hexdigest(),
       'path':str(final if passed else partial),'chunks':ordered,'resource_plan_sha256':h.sha(plan_path),
       'status':'EXACT_HISTORICAL_SOURCE_REACQUIRED' if passed else 'SHA_MISMATCH_PRESERVED','native_analysis_completed':False}
    h.save(PACKAGE/'logs/ovarian_cancer_acquisition_receipt_v3.json',r)
    print(json.dumps({k:v for k,v in r.items() if k!='chunks'}),flush=True)
    if not passed:raise RuntimeError('FULL_HISTORICAL_SHA_FAILED_NO_SUBSTITUTION')

if __name__=='__main__':
    try:main()
    except BaseException as error:
        h.STOP.set();p=PACKAGE/'logs/ovarian_continuation_FAILURE_v3.json'
        if not p.exists():h.save(p,{'completed_utc':h.now(),'status':'FAILED_CONTINUATION_PRESERVED','error':str(error),'downstream_admission':False,'original_inputs_modified':False})
        raise
