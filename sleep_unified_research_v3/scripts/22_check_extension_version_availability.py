#!/usr/bin/env python3
"""HEAD-only availability check of all 100 immutable extension source versions."""
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def check(row):
    path=PACKAGE/'logs/extension_HEAD_v3'/(row['extension_trait_id']+'.txt')
    if path.exists():raise RuntimeError('PRIOR_HEAD_EVIDENCE_PRESERVED')
    began=time.time()
    p=subprocess.run(['curl','--head','--silent','--show-error','--fail','--location','--max-redirs','3',
                      '--connect-timeout','15','--max-time','45','--output',str(path),row['versioned_url']],capture_output=True,text=True)
    text=path.read_text(errors='replace').replace('\r\n','\n') if path.exists() else ''
    block=next((b for b in reversed(text.strip().split('\n\n')) if b.startswith('HTTP/')),'')
    h={}
    for line in block.splitlines()[1:]:
        if ':' in line:
            k,v=line.split(':',1);h[k.lower()]=v.strip()
    code=block.splitlines()[0].split()[1:2] if block else []
    passed=(p.returncode==0 and code==['200'] and h.get('content-length')==row['expected_size_bytes']
            and h.get('etag','').strip('"')==row['observed_etag'].strip('"')
            and h.get('x-amz-version-id')==row['s3_version_id'])
    return {'extension_trait_id':row['extension_trait_id'],'url':row['versioned_url'],
            'expected_bytes':int(row['expected_size_bytes']),'expected_checksum':row['expected_checksum'],
            'expected_version_id':row['s3_version_id'],'expected_ETag':row['observed_etag'],
            'observed_headers':h,'curl_returncode':p.returncode,'stderr':p.stderr,
            'elapsed_seconds':time.time()-began,'header_path':str(path.relative_to(PACKAGE)),
            'header_sha256':sha(path) if path.exists() else None,
            'status':'EXACT_VERSION_HEAD_MATCH_BODY_NOT_REVERIFIED' if passed else 'UNAVAILABLE_OR_CHANGED_NO_SUBSTITUTION'}

def main():
    snapshot=ROOT/'discovery_extension/provenance/panukbb/remote_object_snapshot.tsv'
    with snapshot.open(newline='') as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['object_role']=='phenotype_sumstats']
    if len(rows)!=100 or len({r['extension_trait_id'] for r in rows})!=100:raise RuntimeError('LOCKED_PANEL_CARDINALITY_FAILED')
    if any(not r['versioned_url'].startswith('https://pan-ukb-us-east-1.s3.amazonaws.com/') or 'versionId=' not in r['versioned_url'] for r in rows):raise RuntimeError('VERSION_PIN_REQUIRED')
    dest=PACKAGE/'logs/extension_HEAD_v3'
    if dest.exists():raise RuntimeError('PRIOR_HEAD_DIRECTORY_PRESERVED')
    dest.mkdir()
    began=time.time()
    with ThreadPoolExecutor(max_workers=8) as pool: results=list(pool.map(check,rows))
    report={'completed_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.time()-began,
            'snapshot_sha256':sha(snapshot),'code_sha256':sha(__file__),'workers':8,'maximum_seconds_per_HEAD':45,
            'summary_statistic_body_bytes_downloaded':0,'locked_source_total_bytes':sum(r['expected_bytes'] for r in results),
            'matched_HEAD_versions':sum(r['status']=='EXACT_VERSION_HEAD_MATCH_BODY_NOT_REVERIFIED' for r in results),
            'source_body_verification_or_native_reproduction_completed':False,'results':results}
    with (PACKAGE/'logs/extension_version_availability_v3.json').open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'}),flush=True)

if __name__=='__main__':main()
