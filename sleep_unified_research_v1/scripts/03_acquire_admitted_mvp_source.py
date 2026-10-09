#!/usr/bin/env python3
"""Acquire one metadata-admitted public sleep GWAS; never inspect rg outcomes."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/'reviews/provenance_v1_mvp_insomnia_source_contract.json'
DEST=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/new_validation_sources')

def main():
    c=json.loads(CONTRACT.read_text())
    protocol=ROOT/'FROZEN_NEW_ANALYSIS_PROTOCOL.md'
    if not protocol.is_file(): raise SystemExit('RESULT_FREE_PROTOCOL_REQUIRED')
    DEST.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(DEST).free < 5*1024**3: raise SystemExit('SSD_CAPACITY_GATE_FAILED')
    out=DEST/c['filename']; temp=out.with_name(out.name+'.downloading')
    if out.exists() or temp.exists(): raise SystemExit('PRESERVE_EXISTING_SOURCE; inspect receipt before any retry')
    h=hashlib.sha256(); md5=hashlib.md5(); n=0
    with urllib.request.urlopen(c['url'],timeout=120) as source, temp.open('xb') as destination:
        headers=dict(source.headers.items())
        while block:=source.read(4*1024*1024):
            destination.write(block); h.update(block); md5.update(block); n+=len(block)
    verified=n==c['bytes'] and md5.hexdigest()==c['md5']
    receipt={'completed_utc':datetime.now(timezone.utc).isoformat(),'url':c['url'],
        'path':str(out if verified else temp),'expected_bytes':c['bytes'],'actual_bytes':n,
        'expected_md5':c['md5'],'actual_md5':md5.hexdigest(),'actual_sha256':h.hexdigest(),
        'headers':headers,'verification_status':'PASS' if verified else 'FAIL',
        'protocol_sha256':hashlib.sha256(protocol.read_bytes()).hexdigest(),
        'source_contract_sha256':hashlib.sha256(CONTRACT.read_bytes()).hexdigest(),
        'genetic_correlation_outcomes_accessed':False,'original_input_files_modified':False}
    (ROOT/'logs/mvp_insomnia_acquisition_receipt_v1.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if not verified: raise SystemExit('SOURCE_HASH_FAILED_PARTIAL_PRESERVED')
    temp.rename(out)
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()
