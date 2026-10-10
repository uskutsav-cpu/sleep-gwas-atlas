#!/usr/bin/env python3
"""Owned, bounded stdlib source gate or full processed-template comparison."""
import argparse
import datetime
import importlib.metadata
import json
from pathlib import Path
import sys

from extension_replay_common_v2 import sha,body_hashes,write_new,check_bindings,acquisition_receipt_gate,compare_munged,harmonized_receipt_gate,checkpoint_binding_gate


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha',required=True);p.add_argument('--admission',type=Path,required=True);p.add_argument('--admission-sha',required=True);p.add_argument('--trait',required=True);p.add_argument('--mode',choices=['source','compare'],required=True);p.add_argument('--checkpoint-binding',type=Path,required=True);p.add_argument('--checkpoint-binding-sha',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());member=next(x for x in plan['members'] if x['extension_trait_id']==a.trait)
    record=dict(schema='historical_extension_raw_pipeline_validation_v1',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                extension_trait_id=a.trait,mode=a.mode,plan_sha256=a.plan_sha,admission_sha256=a.admission_sha,status='FAILED_PRESERVED_NO_AUTOMATIC_RETRY')
    try:
        if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('PLAN_OR_ADMISSION_CHANGED')
        check_bindings(plan)
        admission=json.loads(a.admission.read_text())
        if not admission['execution_admitted'] or admission['plan_sha256']!=a.plan_sha:raise RuntimeError('ROOT_OPERATIONAL_ADMISSION_DIFFERS')
        if sha(a.checkpoint_binding)!=a.checkpoint_binding_sha:raise RuntimeError('DURABLE_SOURCE_CHECKPOINT_BINDING_CHANGED')
        binding=json.loads(a.checkpoint_binding.read_text())
        checkpoint_binding_gate(plan,member,binding,a.plan_sha,a.admission_sha)
        expected=binding['acquisition_receipt_sha256']
        record['checkpoint_binding_sha256']=a.checkpoint_binding_sha;record['acquisition_receipt_sha256']=expected
        environment={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy']}
        if not sys.version.startswith('3.9.23') or environment!={k:plan['environment'][k] for k in environment}:raise RuntimeError('PINNED_PYTHON_ENVIRONMENT_DIFFERS')
        record['environment_versions']=environment;record['python']=sys.version
        raw=body_hashes(member['raw'])
        if raw!={'sha256':member['raw_sha256'],'md5':member['raw_md5'],'bytes':member['raw_bytes']}:raise RuntimeError('ENTIRE_RAW_SHA_MD5_OR_BYTE_COUNT_DIFFERS')
        record['whole_raw_identity']=raw
        if a.mode=='source':
            record['gzip_CRC_and_EOF_scope']='Raw compressed hashes only here; unchanged harmonizer must subsequently decode to EOF and verify CRC.'
            record['status']='EXACT_RAW_SOURCE_GATE_PASS'
        else:
            before=json.loads(Path(member['source_gate_receipt']).read_text())
            if before['status']!='EXACT_RAW_SOURCE_GATE_PASS' or before['plan_sha256']!=a.plan_sha or before['whole_raw_identity']!=raw or before.get('checkpoint_binding_sha256')!=a.checkpoint_binding_sha or before.get('acquisition_receipt_sha256')!=expected or before.get('admission_sha256')!=a.admission_sha or before.get('extension_trait_id')!=a.trait:raise RuntimeError('RAW_BEFORE_AFTER_CHANGED_OR_UNCLEARED')
            if sha(member['archived_munged'])!=member['archived_munged_sha256']:raise RuntimeError('ARCHIVED_MUNGED_CHANGED')
            record['harmonization']=harmonized_receipt_gate(plan,member)
            comparison=compare_munged(member['munged'],member['archived_munged'],plan['template_rows'])
            record['processed_comparison']=comparison
            record['status']=comparison['status']
        check_bindings(plan)
        if sha(a.plan)!=a.plan_sha or sha(a.admission)!=a.admission_sha:raise RuntimeError('FINAL_PLAN_OR_ADMISSION_CHANGED')
        if sha(a.checkpoint_binding)!=a.checkpoint_binding_sha:raise RuntimeError('FINAL_SOURCE_CHECKPOINT_BINDING_CHANGED')
        checkpoint_binding_gate(plan,member,binding,a.plan_sha,a.admission_sha)
    except BaseException as e:
        record['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY';record['error']=type(e).__name__+': '+str(e)
        raise
    finally:
        record['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();write_new(a.out,record)
    if record['status'] not in ['EXACT_RAW_SOURCE_GATE_PASS','EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH']:raise SystemExit('REPLAY_COMPARISON_MISMATCH_PRESERVED')
    print(json.dumps(dict(status=record['status'],receipt=str(a.out),sha256=sha(a.out),estimator_calls=0)),flush=True)


if __name__=='__main__':main()
