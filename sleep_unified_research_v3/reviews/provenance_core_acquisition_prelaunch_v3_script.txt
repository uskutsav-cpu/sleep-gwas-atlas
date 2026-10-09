#!/usr/bin/env python3
"""Recover three exact historical public containers with bounded byte ranges."""
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import threading
import time

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent
DEST=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/remaining_core_containers_v3')
CHUNK=64*1024**2
WORKERS=8
STOP=threading.Event()
BEGAN=None

def guard():
    if STOP.is_set():raise RuntimeError('TRANSFER_STOP_REQUESTED')
    if shutil.disk_usage(ROOT).free<128*1024**2:
        STOP.set();raise RuntimeError('INTERNAL_EMERGENCY_FLOOR')
    if shutil.disk_usage(DEST.parent).free<5*1024**3+WORKERS*CHUNK:
        STOP.set();raise RuntimeError('SSD_RESERVE_AND_IN_FLIGHT_FLOOR')
    if BEGAN is not None and time.monotonic()-BEGAN>5400:
        STOP.set();raise RuntimeError('TOTAL_TIME_LIMIT')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def now():return datetime.now(timezone.utc).isoformat()

def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')

def headers(path):
    text=Path(path).read_text(errors='replace').replace('\r\n','\n')
    block=next((b for b in reversed(text.strip().split('\n\n')) if b.startswith('HTTP/')),'')
    result={'status_line':block.splitlines()[0] if block else ''}
    for line in block.splitlines()[1:]:
        if ':' in line:
            k,v=line.split(':',1);result[k.lower()]=v.strip()
    return result

def fetch(s,identity,start,end):
    directory=DEST/s['trait_ids']/'chunks';stem='%012d-%012d'%(start,end)
    body=directory/(stem+'.bin');head=directory/(stem+'.headers');receipt=directory/(stem+'.json')
    if any(p.exists() for p in [body,head,receipt]):raise RuntimeError('PRIOR_CHUNK_PRESERVED')
    n=end-start+1; began=time.monotonic()
    conditional='If-Match: '+identity['etag'] if identity.get('etag') else 'If-Unmodified-Since: '+identity['last-modified']
    command=['curl','--silent','--show-error','--fail','--location','--max-redirs','3','--proto','=https','--proto-redir','=https',
             '--connect-timeout','20','--max-time','600','--max-filesize',str(n),
             '--range','%d-%d'%(start,end),'--header','Accept-Encoding: identity',
             '--header',conditional,'--dump-header',str(head),'--output',str(body),s['download_url']]
    guard()
    proc=subprocess.Popen(command,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
    stopped=False;escalated=False
    while proc.poll() is None:
        if STOP.wait(1):
            stopped=True;proc.terminate()
            try:proc.wait(timeout=30)
            except subprocess.TimeoutExpired:escalated=True;proc.kill();proc.wait()
            break
    _,stderr=proc.communicate()
    observed=headers(head) if head.exists() else {}
    length=body.stat().st_size if body.exists() else 0
    stable=observed.get('etag')==identity['etag'] if identity.get('etag') else observed.get('last-modified')==identity['last-modified']
    passed=(not stopped and not STOP.is_set() and proc.returncode==0 and length==n and observed.get('status_line','').split()[1:2]==['206']
            and observed.get('content-range')=='bytes %d-%d/%s'%(start,end,s['archive_bytes'])
            and observed.get('content-length')==str(n) and stable)
    r={'start':start,'end':end,'expected_bytes':n,'actual_bytes':length,'headers':observed,
       'curl_returncode':proc.returncode,'stderr':stderr,'elapsed_seconds':time.monotonic()-began,
       'shared_stop_observed':stopped,'owned_curl_termination_escalated':escalated,
       'body_path':str(body),'body_sha256':sha(body) if passed else None,
       'header_sha256':sha(head) if head.exists() else None,'status':'PASS_EXACT_RANGE' if passed else 'FAIL_PRESERVED'}
    save(receipt,r)
    if not passed:raise RuntimeError('RANGE_TRANSFER_FAILED '+str(receipt))
    return r

def main():
    global BEGAN
    registry=ROOT/'config/public_gwas_sources.tsv'
    with registry.open(newline='') as f:selected=[r for r in csv.DictReader(f,delimiter='\t') if r['trait_ids'] in {'breast_cancer','prostate_cancer','ovarian_cancer'}]
    if len(selected)!=3 or {s['trait_ids'] for s in selected}!={'breast_cancer','prostate_cancer','ovarian_cancer'}:raise RuntimeError('THREE_SOURCE_REGISTRY_GATE_FAILED')
    upstream={}; head_hashes={}
    for s in selected:
        short=s['trait_ids'].split('_')[0];path=PACKAGE/'logs'/(short+'_public_HEAD_v3.txt');h=headers(path)
        if h.get('status_line','').split()[1:2]!=['200'] or h.get('content-length')!=s['archive_bytes'] or h.get('accept-ranges')!='bytes':
            raise RuntimeError('HEAD_IDENTITY_OR_RANGE_GATE_FAILED '+s['trait_ids'])
        if not h.get('etag') and not h.get('last-modified'):raise RuntimeError('NO_CONDITIONAL_IDENTITY_GATE')
        upstream[s['trait_ids']]=h;head_hashes[s['trait_ids']]=sha(path)
    if DEST.exists():raise RuntimeError('PRIOR_DESTINATION_PRESERVED')
    total=sum(int(s['archive_bytes']) for s in selected)
    reserve=2*total+5*1024**3
    if shutil.disk_usage(DEST.parent).free<reserve or shutil.disk_usage(ROOT).free<512*1024**2:
        raise RuntimeError('RESOURCE_GATE_FAILED')
    version=subprocess.check_output(['curl','--version'],text=True).splitlines()[0]
    if tuple(map(int,version.split()[1].split('.')))<(8,4,0):raise RuntimeError('CURL_STREAMING_CAP_REQUIRED')
    DEST.mkdir(parents=True)
    plan={'frozen_utc':now(),'sources':selected,'upstream_headers':upstream,'HEAD_sha256':head_hashes,
          'source_registry_sha256':sha(registry),'script_sha256':sha(__file__),
          'network_bytes':total,'retention_bytes_full_and_chunks':2*total,'SSD_reserve_bytes':5*1024**3,
          'SSD_runtime_floor_with_in_flight_allowance_bytes':5*1024**3+WORKERS*CHUNK,'deadline_clock':'monotonic',
          'SSD_free_before':shutil.disk_usage(DEST).free,'workers':WORKERS,'chunk_bytes':CHUNK,
          'max_seconds_per_chunk':600,'maximum_all_sources_seconds':5400,'no_retries':True,
          'internal_source_bytes':0,'hash_buffer_bytes':65536,'internal_emergency_floor_bytes':128*1024**2,
          'curl_version':version,'TLS_verification_disabled':False,
          'breast_identity_limitation':'No ETag supplied; exact size/ranges plus Last-Modified conditional and final historical SHA are mandatory',
          'original_inputs_or_sealed_packages_modified':False,'scientific_outcomes_accessed':False}
    plan_path=PACKAGE/'manifests/remaining_core_containers_resource_plan_v3.json';save(plan_path,plan)
    print(json.dumps({'status':'RESOURCE_PLAN_FROZEN','network_bytes':total,'retention_bytes':2*total,'path':str(plan_path)}),flush=True)
    began=time.monotonic();BEGAN=began; acquired=[]
    for s in selected:
        n=int(s['archive_bytes']);(DEST/s['trait_ids']/'chunks').mkdir(parents=True)
        intervals=[(a,min(a+CHUNK,n)-1) for a in range(0,n,CHUNK)]; results=[]
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            pending=set();cursor=0
            try:
                while cursor<len(intervals) or pending:
                    guard()
                    while cursor<len(intervals) and len(pending)<WORKERS:
                        a,b=intervals[cursor];pending.add(pool.submit(fetch,s,upstream[s['trait_ids']],a,b));cursor+=1
                    done,_=wait(pending,timeout=1,return_when=FIRST_COMPLETED)
                    for future in done:
                        pending.remove(future);results.append(future.result())
                        print(json.dumps({'trait':s['trait_ids'],'chunks_complete':len(results),'chunks_total':len(intervals)}),flush=True)
            except BaseException:
                STOP.set()
                for future in pending:future.cancel()
                raise
        ordered=sorted(results,key=lambda r:r['start'])
        if [(r['start'],r['end']) for r in ordered]!=intervals or sum(r['actual_bytes'] for r in ordered)!=n:raise RuntimeError('RANGE_COVERAGE_FAILED')
        partial=DEST/s['trait_ids']/(s['archive_name']+'.assembling');final=DEST/s['trait_ids']/s['archive_name']
        digest=hashlib.sha256()
        with partial.open('xb') as output:
            for r in ordered:
                guard()
                if sha(r['body_path'])!=r['body_sha256']:raise RuntimeError('CHUNK_HASH_CHANGED')
                with Path(r['body_path']).open('rb') as f:
                    checked=time.monotonic()
                    for b in iter(lambda:f.read(65536),b''):
                        output.write(b);digest.update(b)
                        if time.monotonic()-checked>1:guard();checked=time.monotonic()
        guard()
        passed=partial.stat().st_size==n and digest.hexdigest()==s['archive_sha256']
        if passed:partial.rename(final)
        result={'completed_utc':now(),'trait_id':s['trait_ids'],'url':s['download_url'],
                'expected_bytes':n,'actual_bytes':partial.stat().st_size if not passed else final.stat().st_size,
                'expected_sha256':s['archive_sha256'],'actual_sha256':digest.hexdigest(),
                'path':str(final if passed else partial),'chunks':ordered,'resource_plan_sha256':sha(plan_path),
                'status':'EXACT_HISTORICAL_SOURCE_REACQUIRED' if passed else 'SHA_MISMATCH_PRESERVED',
                'native_analysis_completed':False}
        save(PACKAGE/'logs'/(s['trait_ids']+'_acquisition_receipt_v3.json'),result)
        if not passed:raise RuntimeError('FULL_HISTORICAL_SHA_FAILED_NO_SUBSTITUTION')
        acquired.append({k:v for k,v in result.items() if k!='chunks'})
        print(json.dumps({'trait':s['trait_ids'],'status':result['status'],'sha256':digest.hexdigest()}),flush=True)
    save(PACKAGE/'logs/remaining_core_containers_receipt_v3.json',{'completed_utc':now(),'elapsed_seconds':time.monotonic()-began,
         'sources':acquired,'all_exact_expected_hashes':True,'resource_plan_sha256':sha(plan_path),
         'network_bytes':total,'original_inputs_modified':False,'scientific_analyses_completed':False})

if __name__=='__main__':
    try:main()
    except BaseException as error:
        STOP.set()
        path=PACKAGE/'logs/remaining_core_containers_FAILURE_v3.json'
        if not path.exists():
            save(path,{'completed_utc':now(),'error_type':type(error).__name__,'error':str(error),
                       'elapsed_seconds':time.monotonic()-BEGAN if BEGAN is not None else None,
                       'status':'FAILED_TRANSFER_PRESERVED_NO_DOWNSTREAM_ADMISSION',
                       'active_task_owned_curl_stop_requested':True,'original_inputs_modified':False})
        raise
