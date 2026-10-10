#!/usr/bin/env python3
"""Full ordered original validation output/QC comparison, without estimators."""
import argparse
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
import math
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(v.is_symlink() for v in p.parents):
        raise RuntimeError('REGULAR_VALIDATION_REPLAY_EVIDENCE_REQUIRED: '+str(p))
    return p

def fixed_hashes(mapping):
    if not mapping:raise RuntimeError('NONEMPTY_FROZEN_REPLAY_MAP_REQUIRED')
    for path,digest in mapping.items():
        if sha(regular(path))!=digest:raise RuntimeError('FROZEN_VALIDATION_REPLAY_EVIDENCE_CHANGED: '+path)

def load_adapter(plan):
    spec=spec_from_file_location('_original_validation_local_adapter',plan['adapter_path'])
    m=module_from_spec(spec);spec.loader.exec_module(m);return m

def current_source(plan,member,adapter):
    adapter.completed_acquisition(plan)
    body=regular(member['body_path']);expected=member['original_source_identity']
    md5,digest=adapter.hashes(body)
    if body.stat().st_size!=expected['bytes'] or (md5,digest)!=(expected['md5'],expected['sha256']):
        raise RuntimeError('EXACT_FULL_VALIDATION_SOURCE_BODY_REQUIRED')
    return dict(path=str(body),bytes=body.stat().st_size,md5=md5,sha256=digest)

def qc_rows(path):
    with regular(path).open(newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        if r.fieldnames!=['metric','value','notes']:raise RuntimeError('ORIGINAL_VALIDATION_QC_SCHEMA_CHANGED')
        rows=list(r)
    if len({x['metric'] for x in rows})!=len(rows):raise RuntimeError('DUPLICATED_QC_METRIC')
    return rows

def literal_streams(candidate,original,expected_rows):
    hashes=[hashlib.sha256(),hashlib.sha256()];rows=0;seen=set();max_abs_z=0.
    with gzip.open(regular(candidate),'rb') as a,gzip.open(regular(original),'rb') as b:
        header_a,header_b=a.readline(),b.readline()
        if header_a!=b'SNP\tA1\tA2\tZ\tN\n' or header_a!=header_b:raise RuntimeError('ORIGINAL_ORDERED_HEADER_CHANGED')
        hashes[0].update(header_a);hashes[1].update(header_b)
        for left,right in zip_longest(a,b):
            if left is None or right is None or left!=right:raise RuntimeError('LITERAL_VALIDATION_STREAM_DIFFERS_AT_ROW_'+str(rows+1))
            for h,line in zip(hashes,[left,right]):h.update(line)
            fields=left.decode('utf-8').rstrip('\n').split('\t')
            if len(fields)!=5 or fields[0] in seen or fields[1] not in {'A','C','G','T'} or fields[2] not in {'A','C','G','T'} or fields[1]==fields[2] or set(fields[1:3]) in [{'A','T'},{'C','G'}]:
                raise RuntimeError('INVALID_OR_DUPLICATED_LITERAL_VALIDATION_ROW')
            seen.add(fields[0]);z,n=float(fields[3]),float(fields[4])
            if not math.isfinite(z) or not math.isfinite(n) or n<=0:raise RuntimeError('NONFINITE_VALIDATION_ROW')
            max_abs_z=max(max_abs_z,abs(z));rows+=1
    if rows!=expected_rows:raise RuntimeError('FULL_ORIGINAL_VALIDATION_ROW_COUNT_CHANGED')
    return dict(ordered_decompressed_SHA256=[h.hexdigest() for h in hashes],data_rows=rows,
        literal_header_SNP_allele_Z_N_order_equal=True,both_full_gzip_CRC_and_EOF_verified=True,
        finite_Z_N_rows=rows,literal_missing_Z_N_rows=0,max_abs_Z=max_abs_z)

def run(plan_path,expected_plan_sha,source_id,mode,out):
    if sha(regular(plan_path))!=expected_plan_sha:raise RuntimeError('FROZEN_VALIDATION_PLAN_CHANGED')
    plan=json.loads(Path(plan_path).read_text());fixed_hashes(plan['dependencies_sha256'])
    if plan['member_count']!=13 or plan['new_estimator_fits']!=0:raise RuntimeError('ORIGINAL_VALIDATION_PREPROCESSING_SCOPE_CHANGED')
    selected=[m for m in plan['members'] if m['source_id']==source_id]
    if len(selected)!=1:raise RuntimeError('ONE_REGISTERED_SOURCE_REQUIRED')
    member=selected[0];adapter=load_adapter(plan);body=current_source(plan,member,adapter)
    expected_out=member['source_gate_receipt' if mode=='source' else 'comparison_receipt']
    if str(out)!=expected_out:raise RuntimeError('EXACT_REGISTERED_VALIDATOR_OUTPUT_REQUIRED')
    frozen={str(plan_path):expected_plan_sha,member['original_munged']:member['original_munged_sha256'],
        member['original_qc']:member['original_qc_sha256'],member['original_receipt']:member['original_receipt_sha256']}
    if mode=='source':
        fixed_hashes(frozen)
        result=dict(status='EXACT_ORIGINAL_VALIDATION_SOURCE_GATE_PASS',source_id=source_id,
            plan_sha256=expected_plan_sha,source_body=body,acquisition_terminal_current=True,
            full_source_CRC_EOF_verified=False,raw_filter_replay_complete=False,new_estimator_fits=0)
    else:
        new_keys=['new_munged','new_qc','new_receipt','adapter_receipt','source_gate_receipt']
        frozen.update({member[k]:sha(regular(member[k])) for k in new_keys})
        fixed_hashes(frozen)
        source=json.loads(Path(member['source_gate_receipt']).read_text())
        if source['status']!='EXACT_ORIGINAL_VALIDATION_SOURCE_GATE_PASS' or source['plan_sha256']!=expected_plan_sha or source['source_body']!=body:
            raise RuntimeError('EXACT_PRIOR_CURRENT_SOURCE_GATE_REQUIRED')
        streams=literal_streams(member['new_munged'],member['original_munged'],member['expected_output_rows'])
        oldqc,newqc=qc_rows(member['original_qc']),qc_rows(member['new_qc'])
        old={r['metric']:r for r in oldqc};new={r['metric']:r for r in newqc}
        if list(old)!=list(new):raise RuntimeError('ORIGINAL_QC_METRIC_MEMBERSHIP_OR_ORDER_CHANGED')
        for metric,row in old.items():
            if metric=='munged_output_sha256':
                if row['value']!=member['original_munged_sha256'] or new[metric]['value']!=frozen[member['new_munged']] or row['notes']!=new[metric]['notes']:
                    raise RuntimeError('QC_COMPRESSED_OUTPUT_IDENTITY_CHANGED')
            elif row!=new[metric]:raise RuntimeError('FULL_ORIGINAL_VALIDATION_QC_DIFFERS: '+metric)
        if new['effective_sample_size']['value']!=member['effective_N_serialized_12g'] or int(new['output_rows']['value'])!=streams['data_rows']:
            raise RuntimeError('ORIGINAL_SAMPLE_N_OR_OUTPUT_COUNT_CHANGED')
        r=json.loads(Path(member['new_receipt']).read_text());sidecar=json.loads(Path(member['adapter_receipt']).read_text())
        identity=member['original_source_identity'];verification=r['source_verification']
        if r['pipeline_status']!='STREAM_HARMONIZE_DIRECT_MUNGE_PASS' or r['replication_source_id']!=source_id or r['source_url']!=member['frozen_source_metadata']['replication_source_url']:
            raise RuntimeError('EXACT_ORIGINAL_COLLECTOR_RESULT_IDENTITY_REQUIRED')
        if verification['verification_status']!='PASS' or verification['observed_sha256']!=identity['sha256'] or verification['observed_md5']!=identity['md5'] or verification['observed_size_bytes']!=identity['bytes'] or verification['source']!=member['body_path']:
            raise RuntimeError('EXACT_FULL_LOCAL_STREAM_VERIFICATION_REQUIRED')
        if r['munged_output_sha256']!=frozen[member['new_munged']] or r['harmonization_qc_sha256']!=frozen[member['new_qc']] or r['queue_sha256']!=plan['dependencies_sha256'][plan['queue']] or r['reference_sha256']!=plan['dependencies_sha256'][plan['reference']] or r['hm3_alleles_sha256']!=plan['dependencies_sha256'][plan['alleles']]:
            raise RuntimeError('EXACT_COLLECTOR_INPUT_OUTPUT_HASH_CLOSURE_REQUIRED')
        expected_code={path:plan['dependencies_sha256'][path] for path in [plan['original_collector'],plan['original_streaming_io']]}
        if r['pipeline_code_sha256']!=expected_code:raise RuntimeError('UNCHANGED_ORIGINAL_COLLECTOR_CODE_REQUIRED')
        expected_outputs={member[k]:frozen[member[k]] for k in ['new_munged','new_receipt','new_qc']}
        if sidecar['plan_sha256']!=expected_plan_sha or sidecar['source_id']!=source_id or sidecar['original_source_identity']!=identity or sidecar['output_sha256']!=expected_outputs or sidecar['full_local_source_retained'] is not True or sidecar['full_source_gzip_CRC_EOF_MD5_size_SHA256_verified'] is not True or sidecar['source_or_independent_replication_admitted'] is not False:
            raise RuntimeError('TRUTHFUL_LOCAL_REPLAY_SIDECAR_REQUIRED')
        result=dict(status='EXACT_ORIGINAL_VALIDATION_ORDERED_CONTENT_QC_REPLAY_PASS',source_id=source_id,
            plan_sha256=expected_plan_sha,source_body=body,source_gzip_CRC_EOF_verified=True,
            stream_comparison=streams,original_QC_metric_rows=len(oldqc),all_original_QC_metrics_and_N_strings_equal=True,
            compressed_byte_identity=frozen[member['new_munged']]==member['original_munged_sha256'],
            compressed_metadata_difference_cause_not_inferred=True,consumed_output_sha256=frozen,
            original_217_membership_and_41_estimates_unchanged=True,native_estimator_inputs_literal_equivalent=True,
            native190_already_executed_on_verified_original_processed_inputs=True,new_estimator_fits=0,
            independent_two_trait_replication_established=False,
            source_qualification='Historical rsid/allele projection and constant effective N preserved; INFO/per-variant N and shared-sleep/cohort/phenotype limitations are not cleared.')
    fixed_hashes(frozen);current_source(plan,member,adapter);fixed_hashes(plan['dependencies_sha256'])
    payload=json.dumps(result,indent=2,allow_nan=False)+'\n';intended=hashlib.sha256(payload.encode()).hexdigest()
    with Path(out).open('x') as f:f.write(payload)
    if sha(regular(out))!=intended:raise RuntimeError('INTENDED_VALIDATOR_RECEIPT_CHANGED')
    fixed_hashes(frozen);current_source(plan,member,adapter);fixed_hashes(plan['dependencies_sha256'])

def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha256',required=True);p.add_argument('--source-id',required=True);p.add_argument('--mode',choices=['source','compare'],required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.plan_sha256,a.source_id,a.mode,a.out)

if __name__=='__main__':main()
