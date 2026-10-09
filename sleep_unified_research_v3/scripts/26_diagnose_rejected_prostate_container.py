#!/usr/bin/env python3
"""Metadata-only row diagnostic of rejected bytes; no raw projection or result."""
from collections import Counter
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import shutil
import sys
import time

PACKAGE=Path(__file__).resolve().parents[1];ROOT=PACKAGE.parent
BEGAN=time.monotonic()

def guard():
    if shutil.disk_usage(ROOT).free<128*1024**2:raise RuntimeError('INTERNAL_EMERGENCY_FLOOR')
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>128*1024**2:raise RuntimeError('SELF_MAX_RSS_LIMIT')
    if time.monotonic()-BEGAN>1800:raise RuntimeError('TIME_LIMIT')

def sha(path):
    h=hashlib.sha256();checked=time.monotonic()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
            if time.monotonic()-checked>1:guard();checked=time.monotonic()
    return h.hexdigest()

def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')

def main():
    if shutil.disk_usage(ROOT).free<256*1024**2:raise RuntimeError('DIAGNOSTIC_LAUNCH_FLOOR')
    rpath=PACKAGE/'logs/prostate_cancer_acquisition_receipt_v3.json';r=json.loads(rpath.read_text())
    if r['status']!='SHA_MISMATCH_PRESERVED':raise RuntimeError('REJECTED_SOURCE_SCOPE_REQUIRED')
    source=Path(r['path']);historical=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/raw/prostate_cancer.txt.gz.provenance.json')
    native=ROOT/'scripts/20_materialize_practical_prostate.py';registry=ROOT/'config/public_gwas_sources.tsv'
    paths={'rejected_source':source,'receipt':rpath,'historical_provenance':historical,'unchanged_native_policy':native,'registry':registry,'code':Path(__file__)}
    before={k:sha(v) for k,v in paths.items()}
    if before['rejected_source']!=r['actual_sha256'] or source.stat().st_size!=r['expected_bytes']:raise RuntimeError('REJECTED_BYTES_CHANGED')
    plan={'frozen_utc':datetime.now(timezone.utc).isoformat(),'bound_paths':{k:str(v) for k,v in paths.items()},
          'inputs_sha256_before':before,'scope':'schema and source row-category counts only; no output projection, P/effect calculation, h2 or rg',
          'frozen_expected_source_gate_bypassed_for_science':False,'source_remains_rejected':True,
          'internal_launch_bytes':256*1024**2,'internal_emergency_bytes':128*1024**2,
          'self_peak_RSS_limit_bytes':128*1024**2,'maximum_seconds':1800,'network_bytes':0,
          'resource_basis':'one bounded source line at a time, counters and ten line-number examples'}
    plan_path=PACKAGE/'manifests/rejected_prostate_metadata_diagnostic_plan_v3.json';save(plan_path,plan)
    # Import constants only; never invoke native verify/materialize and never
    # override its historical SHA gate or frozen expected counts.
    sys.path.insert(0,str(ROOT/'scripts'))
    spec=importlib.util.spec_from_file_location('native_policy_constants',native);policy=importlib.util.module_from_spec(spec);spec.loader.exec_module(policy)
    c=Counter();examples=[];header_matches=False
    with source.open('rb') as f:
        header=f.readline().rstrip(b'\r\n');header_matches=header==policy.HEADER
        if not header_matches:raise RuntimeError('REJECTED_SOURCE_SCHEMA_CHANGED')
        for line_number,line in enumerate(f,start=2):
            if len(line)>128*1024:raise RuntimeError('BOUNDED_ROW_SIZE_EXCEEDED')
            c['source_rows']+=1;stripped=line.rstrip(b'\r\n');field_count=stripped.count(b'\t')+1
            if field_count!=policy.EXPECTED_FIELD_COUNT:
                c['malformed_rows']+=1
                if len(examples)<10:examples.append({'line_number':line_number,'field_count':field_count})
            else:
                fields=stripped.split(b'\t')
                if policy.RSID.fullmatch(fields[2]) is None:c['non_rsid_rows']+=1
                elif fields[5].upper() not in policy.VALID_ALLELES or fields[6].upper() not in policy.VALID_ALLELES:c['non_snp_rows']+=1
                elif fields[3] not in policy.AUTOSOMES:c['nonautosomal_rows']+=1
                else:c['retained_rows']+=1
            if c['source_rows']%100000==0:guard()
    old=json.loads(historical.read_text());keys=list(old['counts']);comparison={k:{'historical':old['counts'][k],'rejected_current':c[k],'difference':c[k]-old['counts'][k]} for k in keys}
    after={k:sha(v) for k,v in paths.items()};guard()
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'status':'REJECTED_SOURCE_METADATA_DIAGNOSTIC_ONLY',
            'header_matches_native_schema':header_matches,'counts':{k:c[k] for k in keys},'comparison_to_historical':comparison,
            'first_ten_malformed_line_numbers_and_field_counts':examples,'historical_first_malformed_examples':old['malformed_examples'],
            'exclusive_row_categories_reconcile':c['source_rows']==sum(c[k] for k in keys if k!='source_rows'),
            'inputs_sha256_before':before,'inputs_sha256_after':after,'input_hashes_unchanged':before==after,
            'execution_plan_sha256':sha(plan_path),'elapsed_seconds':time.monotonic()-BEGAN,
            'observed_self_max_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'source_still_rejected':True,'source_hash_discrepancy_resolved':False,'raw_projection_created':False,'h2_or_rg_estimated':False,
            'interpretation':'Counts can identify schema/category differences but cannot establish the cause or certify allele/effect data equality'}
    save(PACKAGE/'logs/rejected_prostate_metadata_diagnostic_v3.json',result)
    print(json.dumps({'status':result['status'],'comparison':comparison,'input_hashes_unchanged':before==after}),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException as error:
        p=PACKAGE/'logs/rejected_prostate_metadata_FAILURE_v3.json'
        if not p.exists():save(p,{'status':'FAILED_METADATA_DIAGNOSTIC_PRESERVED','error':str(error),'downstream_admission':False})
        raise
