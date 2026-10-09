#!/usr/bin/env python3
"""Seal small technical artifacts; GWAS bodies stay off Git."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import os
import re
import subprocess
import sys

PACKAGE=Path(__file__).resolve().parents[1];ROOT=PACKAGE.parent

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def main():
    manifest=PACKAGE/'hashes/ARTIFACT_SHA256_v3.tsv'
    if manifest.exists():raise RuntimeError('PRIOR_SEAL_PRESERVED')
    required=['FINAL_EVIDENCE_HANDOFF_V3.md','REVIEW_COMPLETION_V3.md','reviews/cancer_materialization_stream_review_v3.md',
              'reviews/provenance_cancer_source_final_forensics_v3.md','reviews/mvp_preprocessing_stream_review_v3.md',
              'manifests/figure_visual_QA_v3.json']
    if any(not (PACKAGE/p).is_file() for p in required):raise RuntimeError('REQUIRED_FINAL_EVIDENCE_MISSING')
    tests=[];env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    for name,n in [('sleep_unified_research_v1',10),('sleep_unified_research_v2',6),('sleep_unified_research_v3',8)]:
        command=[sys.executable,'-m','unittest','discover','-s',name+'/tests','-p','test_*.py']
        proc=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True)
        text=proc.stdout+proc.stderr
        path=PACKAGE/'logs'/('portable_'+name+'_tests_v3.log')
        with path.open('x') as f:f.write(text)
        passed=proc.returncode==0 and re.search(r'Ran '+str(n)+r' tests?\b',text) is not None
        tests.append({'package':name,'command':command,'returncode':proc.returncode,'expected_tests':n,'pass':passed,'log_sha256':sha(path)})
        if not passed:raise RuntimeError('PORTABLE_ARTIFACT_TEST_FAILED '+name)
    with (PACKAGE/'manifests/portable_test_receipt_v3.json').open('x') as f:
        json.dump({'completed_utc':datetime.now(timezone.utc).isoformat(),'tests':tests,'test_count':24,
                   'scope':'portable artifact invariants; not native scientific validation','sealer_sha256':sha(Path(__file__))},f,indent=2);f.write('\n')
    (PACKAGE/'hashes').mkdir(exist_ok=True)
    rows=[]
    for p in sorted(PACKAGE.rglob('*')):
        if p.is_symlink():raise RuntimeError('NO_SYMLINKED_SOURCE_ARTIFACTS')
        if not p.is_file():continue
        if p.suffix in {'.gz','.bgz','.bin','.zip','.pyc','.sqlite'}:raise RuntimeError('SOURCE_OR_RUNTIME_BODY_MUST_NOT_BE_SEALED '+str(p))
        if p.stat().st_size>5*1024**2:raise RuntimeError('UNREVIEWED_LARGE_ARTIFACT '+str(p))
        rows.append({'sha256':sha(p),'bytes':p.stat().st_size,'path':str(p.relative_to(PACKAGE))})
    with manifest.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sha256','bytes','path'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    print(json.dumps({'status':'SEALED_TECHNICAL_PACKAGE','files_excluding_manifest':len(rows),'bytes_excluding_manifest':sum(r['bytes'] for r in rows),
                      'manifest_sha256':sha(manifest),'portable_tests':24,'native_scientific_goal_complete':False}),flush=True)

if __name__=='__main__':main()
