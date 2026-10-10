#!/usr/bin/env python3
"""Full raw identity and literal recovered core-content/QC comparisons."""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
from extension_replay_common_v3 import sha,write_new,check_bindings,compare_munged


def compare_harmonized(new,old,expected_rows):
    digests=[hashlib.sha256(),hashlib.sha256()];counts=[0,0];unequal=0
    with gzip.open(new,'rb') as a,gzip.open(old,'rb') as b:
        heads=[a.readline(),b.readline()]
        if heads[0]!=heads[1] or heads[0].rstrip(b'\r\n')!=b'SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN':
            raise RuntimeError('HARMONIZED_CANONICAL_HEADER_DIFFERS')
        for h,v in zip(digests,heads):h.update(v)
        while True:
            lines=[a.readline(),b.readline()]
            if not any(lines):break
            for i,line in enumerate(lines):
                if line:
                    if line.rstrip(b'\r\n').count(b'\t')!=9:raise RuntimeError('HARMONIZED_ROW_WIDTH_DIFFERS')
                    counts[i]+=1;digests[i].update(line)
            unequal+=lines[0]!=lines[1]
    return {'status':'EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH' if unequal==0 and counts==[expected_rows,expected_rows] else 'HARMONIZED_CONTENT_MISMATCH_PRESERVED',
        'rows':counts,'expected_rows':expected_rows,'unequal_rows':unequal,
        'decompressed_sha256':[h.hexdigest() for h in digests],
        'compressed_sha256':[sha(new),sha(old)],'compressed_byte_identity_required':False,
        'literal_order_fields_missingness_and_numerical_serialization_compared':True,'both_gzip_CRC_and_EOF_verified':True}


def qc(path):
    metadata={};steps=[];in_steps=False
    with Path(path).open() as f:
        for line in f:
            fields=line.rstrip('\n').split('\t')
            if fields==['step','dropped','remaining']:in_steps=True;continue
            if len(fields)==3 and in_steps:steps.append({'reason':fields[0],'dropped':int(fields[1]),'remaining':int(fields[2])})
            elif len(fields)==2:metadata[fields[0]]=fields[1]
    return metadata,steps


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha',required=True);p.add_argument('--trait',required=True);p.add_argument('--mode',choices=['source','compare'],required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if sha(a.plan)!=a.plan_sha:raise RuntimeError('EXACT_CORE_REPLAY_PLAN_CHANGED')
    plan=json.loads(a.plan.read_text());check_bindings(plan)
    member=next(m for m in plan['members'] if m['trait_id']==a.trait);original=member['original_design']
    r={'status':'FAILED_PRESERVED_NO_AUTOMATIC_RETRY','trait':a.trait,'mode':a.mode,'plan_sha256':a.plan_sha,
       'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'estimator_calls':0}
    try:
        raw=Path(original['raw']['resolved_path']);expected=original['raw']
        if not raw.is_file() or raw.is_symlink() or raw.stat().st_size!=expected['bytes'] or sha(raw)!=expected['sealed_verified_sha256']:
            raise RuntimeError('ENTIRE_EXACT_CORE_RAW_IDENTITY_DIFFERS')
        r['raw_sha256']=expected['sealed_verified_sha256'];r['raw_bytes']=expected['bytes']
        if a.mode=='source':r['status']='EXACT_CORE_RAW_SOURCE_GATE_PASS'
        else:
            source=json.loads(Path(member['source_gate_receipt']).read_text())
            if source.get('status')!='EXACT_CORE_RAW_SOURCE_GATE_PASS' or source.get('trait')!=a.trait or source.get('plan_sha256')!=a.plan_sha or source.get('raw_sha256')!=r['raw_sha256']:
                raise RuntimeError('CORE_SOURCE_BEFORE_AFTER_PROOF_DIFFERS')
            prefilter=original['prefilter']
            if prefilter is not None:
                current=Path(member['harmonize_command'][member['harmonize_command'].index('--infile')+1])
                provenance=Path(member['harmonize_command'][member['harmonize_command'].index('--prefilter-provenance')+1])
                proof=json.loads(provenance.read_text())
                if current.stat().st_size!=prefilter['bytes'] or sha(current)!=prefilter['historical_expected_sha256']:
                    raise RuntimeError('EXACT_ORIGINAL_PREFILTER_OUTPUT_DIFFERS')
                if proof.get('input_sha256')!=r['raw_sha256'] or proof.get('input_bytes')!=r['raw_bytes'] or proof.get('output_sha256')!=sha(current) or proof.get('output_bytes')!=current.stat().st_size or proof.get('source_rows')!=original['raw_source_rows'] or proof.get('retained_rows')!=original['harmonizer_input_rows']:
                    raise RuntimeError('ORIGINAL_PREFILTER_COUNTS_OR_IDENTITIES_DIFFER')
                r['prefilter_proof_sha256']=sha(provenance)
            for key in ['harmonized','munged']:
                baseline=original[key]
                if sha(baseline['path'])!=baseline['sealed_sha256']:raise RuntimeError('RECOVERED_CORE_BASELINE_CHANGED')
            r['harmonized_comparison']=compare_harmonized(member['harmonized'],original['harmonized']['path'],original['harmonized_rows'])
            r['munged_comparison']=compare_munged(member['munged'],original['munged']['path'],original['historical_munge_log']['output_total_rows'])
            observed,steps=qc(member['harmonization_qc']);old=original['historical_QC']
            permitted_path_fields={'infile','variant_map','liftover_chain'}
            comparisons={k:observed.get(k)==v for k,v in old['metadata'].items() if k not in permitted_path_fields}
            comparisons['ordered_filter_steps']=steps==old['ordered_steps']
            r['QC_comparisons']=comparisons;r['execution_path_only_fields']=sorted(permitted_path_fields)
            if r['harmonized_comparison']['status']!='EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH' or r['munged_comparison']['status']!='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' or not all(comparisons.values()):
                raise RuntimeError('COMPLETE_CORE_SERIALIZED_CONTENT_OR_QC_DIFFERS')
            r['status']='EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS'
            r['historical_binary_runtime_identity_claimed']=False
            r['qualification']='Matching declared workflow versions; complete recovered-output content/QC agreement does not establish identical historical per-trait binaries or original compressed-byte identity.'
        check_bindings(plan)
        if sha(a.plan)!=a.plan_sha:raise RuntimeError('FINAL_CORE_REPLAY_PLAN_CHANGED')
    except BaseException as e:
        r['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY';r['error']=type(e).__name__+': '+str(e);raise
    finally:
        r['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();write_new(a.out,r)
    print(json.dumps({'trait':a.trait,'status':r['status'],'receipt':str(a.out),'sha256':sha(a.out)}),flush=True)


if __name__=='__main__':main()
