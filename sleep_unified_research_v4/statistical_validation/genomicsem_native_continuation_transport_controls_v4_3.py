#!/usr/bin/env python3
"""Independent additive companion preservation controls; no original file moved."""
import ast
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
from unittest.mock import patch

P=Path(__file__).resolve().parents[1]
S=P/'scripts'
PLAN=P/'manifests/native_extension_checkpoint_continuation_v4_3.json'
PLAN_SHA='f087b74e1ebaf6b543513d69948d8dbb2232e5190360dd6fce3acc16e68ed5b3'
SCRIPT=S/'61_native_checkpoint_continuation_v3.py'
SCRIPT_SHA='b1309c0bd97736d7282c636956b859290dbe5430a5997137d8be0776e3160283'
OUT=P/'statistical_validation/genomicsem_native_continuation_transport_controls_receipt_v4_3.json'
PRIOR=P/'statistical_validation/genomicsem_native_continuation_controls_receipt_v4_2.json'
PRIOR_SHA='3fbadfcb0ba0ee3ace11ba1d3f09d746371024e1c826ee308c2486adb4b58e93'
sys.dont_write_bytecode=True
sys.path.insert(0,str(S))
import canonical_calibration_common_v4_3 as common


def load(path):
    spec=importlib.util.spec_from_file_location('_independent_transport_preservation',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    checks=[];observed=[]
    def check(value,label):
        if not value:raise RuntimeError('ADDITIVE_TRANSPORT_CONTROL_FAILED: '+label)
        checks.append(label)
    check(common.sha(PLAN)==PLAN_SHA and common.sha(SCRIPT)==SCRIPT_SHA,'exact_frozen_additive_v3_code_and_plan')
    check(common.sha(PRIOR)==PRIOR_SHA,'sealed1400_prior_terminal_and_source_controls_bound')
    prior=json.loads(PRIOR.read_text());v3=json.loads(PLAN.read_text());v2=json.loads((P/'manifests/native_extension_checkpoint_continuation_v4_2.json').read_text())
    for path,h in v3['dependencies_sha256'].items():check(common.sha(path)==h,'v3_frozen_dependency_before '+path)
    for key in ['original190_jobs','completed_jobs','remaining_jobs','completed_receipt_sha256','guard','historical_plan_sha256','prior_failed_monitor_sha256']:
        check(v3[key]==v2[key],'additive_transport_preserves_frozen_science_and_guards '+key)
    funcs=lambda path:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef)}
    a,b=funcs(S/'61_native_checkpoint_continuation.py'),funcs(SCRIPT)
    unchanged=['load','write_new','source_bindings','verify_completed','execute','main']
    for name in unchanged:check(a[name]==b[name],'AST_identical_prior_tested_function '+name)
    check({name for name in a if a[name]!=b[name]}=={'prepare','preserve_partial'},'only_preparation_and_preservation_functions_change')
    old_side={r['path']:r for r in prior['transport_sidecar_inventory']}
    check(len(v3['interrupted_artifacts'])==len(v3['interrupted_sidecars'])==len(old_side)==284,'exact284_biological_and284_transport_counts')
    check(len({r['original_path'] for r in v3['interrupted_artifacts']+v3['interrupted_sidecars']})==568,'unique568_original_paths')
    check(len({r['byte_backup_path'] for r in v3['interrupted_sidecars']})==284,'unique284_literal_backup_destinations')
    for i,(bio,side) in enumerate(zip(v3['interrupted_artifacts'],v3['interrupted_sidecars']),1):
        source=Path(bio['original_path']);target=Path(bio['preserved_path'])
        check(side['original_path']==str(source.with_name('._'+source.name)) and side['preserved_path']==str(target.with_name('._'+target.name)),'exact_paired_companion_mapping '+str(i))
        check((side['bytes'],side['sha256'])==(old_side[side['original_path']]['bytes'],old_side[side['original_path']]['sha256']) and side['biological_estimate'] is False,'original_transport_bytes_hash_and_nonbiological_scope '+str(i))
        check(Path(side['byte_backup_path'])==Path(v3['namespace'])/'interrupted_snoring_sidecar_original_bytes'/('sidecar_'+str(i).zfill(4)+'.bin'),'indexed_literal_backup_owned_namespace '+str(i))
    m=load(SCRIPT)
    with tempfile.TemporaryDirectory(prefix='genomicsem_transport_additive_',dir=P/'statistical_validation') as directory:
        tmp=Path(directory)
        for case in ['explicit_companion','implicit_companion','backup_write_failure','backup_readback_failure','transport_sidecar_hash_failure','ambiguous_companion','final_biological_readback_failure','final_backup_hash_failure']:
            base=tmp/case;base.mkdir();sources=base/'source';sources.mkdir();out=base/'archive';out.mkdir()
            biological=[];sidecars=[]
            for i in range(1,285):
                source=sources/('CONTROL_BIO_'+str(i));side=source.with_name('._'+source.name)
                source.write_bytes(('CONTROL_BIO_BYTES_'+str(i)).encode());side.write_bytes(('CONTROL_LITERAL_TRANSPORT_'+str(i)).encode())
                bio=dict(original_path=str(source),preserved_path=str(out/'interrupted_snoring'/source.name),bytes=source.stat().st_size,sha256=common.sha(source))
                companion=dict(original_path=str(side),preserved_path=str(out/'interrupted_snoring'/side.name),byte_backup_path=str(out/'interrupted_snoring_sidecar_original_bytes'/('sidecar_'+str(i).zfill(4)+'.bin')),bytes=side.stat().st_size,sha256=common.sha(side),biological_estimate=False)
                biological.append(bio);sidecars.append(companion)
            fixture=dict(interrupted_artifacts=biological,interrupted_sidecars=sidecars)
            original_rename=Path.rename;original_open=Path.open;original_sha=common.sha
            moves=[];first_move_backups=[];mutated=[False]
            def tiny_rename(source,target):
                target=Path(target)
                if source.name.startswith('CONTROL_BIO_'):
                    if not moves:first_move_backups.append(all(original_sha(r['byte_backup_path'])==r['sha256'] for r in sidecars))
                    moves.append(source.name)
                    result=original_rename(source,target)
                    ss=source.with_name('._'+source.name);st=target.with_name('._'+target.name)
                    if case in ['implicit_companion','transport_sidecar_hash_failure']:
                        original_rename(ss,st)
                    if case=='transport_sidecar_hash_failure' and len(moves)==1:st.write_bytes(b'CONTROL_CHANGED_TRANSPORT')
                    if case=='ambiguous_companion' and len(moves)==1:st.write_bytes(ss.read_bytes())
                    if len(moves)==284 and case=='final_biological_readback_failure':Path(biological[0]['preserved_path']).write_bytes(b'CONTROL_LATE_CORRUPTION')
                    if len(moves)==284 and case=='final_backup_hash_failure':Path(sidecars[0]['byte_backup_path']).write_bytes(b'CONTROL_LATE_BACKUP_CORRUPTION')
                    return result
                return original_rename(source,target)
            def tiny_open(path,mode='r',*args,**kwargs):
                if case=='backup_write_failure' and mode=='xb' and path.name=='sidecar_0142.bin':raise OSError('CONTROL_BACKUP_IO_FAILURE')
                return original_open(path,mode,*args,**kwargs)
            def tiny_sha(path):
                path=Path(path)
                if case=='backup_readback_failure' and path.name=='sidecar_0142.bin':return '0'*64
                return original_sha(path)
            error=None;proof=None
            with patch.object(m,'OUT',out),patch.object(Path,'rename',tiny_rename),patch.object(Path,'open',tiny_open),patch.object(common,'sha',tiny_sha),patch('subprocess.Popen',side_effect=RuntimeError('REAL_WORKERS_FORBIDDEN_IN_TRANSPORT_CONTROLS')):
                try:proof=m.preserve_partial(fixture)
                except BaseException as e:error=type(e).__name__+': '+str(e)
            events=[json.loads(line) for line in (out/'interrupted_snoring_preservation.jsonl').read_text().splitlines()]
            receipt=out/'interrupted_snoring_preservation_receipt.json'
            if case in ['explicit_companion','implicit_companion']:
                check(error is None and len(moves)==284 and first_move_backups==[True],'all284_literal_backups_verified_before_any_move_'+case)
                check(all(not Path(r['original_path']).exists() and original_sha(r['preserved_path'])==r['sha256'] for r in biological+sidecars),'final_all568_identity_and_source_absence_'+case)
                check(all(original_sha(r['byte_backup_path'])==r['sha256'] for r in sidecars),'final_all284_original_literal_backups_'+case)
                text=json.loads(receipt.read_text())
                check(text['sidecars_count_as_biological_estimates'] is False and len(text['artifacts'])==len(text['transport_sidecars'])==284,'explicit_nonbiological_sidecar_receipt_'+case)
                kinds=[r['event'] for r in events]
                check(kinds[:568]==['LITERAL_SIDECAR_BACKUP_INTENT','LITERAL_SIDECAR_BACKUP_READBACK_VERIFIED']*284 and kinds[568:]==['PAIRED_SAME_VOLUME_PRESERVATION_INTENT','PAIR_PRESERVED_READBACK_VERIFIED']*284,'fsynced_backup_and_paired_intent_commit_order_'+case)
                expected='EXPLICIT_COMPANION_RENAME' if case=='explicit_companion' else 'OBSERVED_IMPLICIT_COMPANION_RENAME'
                check({r['transport'] for r in events if r['event']=='PAIR_PRESERVED_READBACK_VERIFIED'}=={expected},'both_companion_transport_branches_'+case)
                check(all(original_sha(p)==h for p,h in proof.items()),'journal_and_exclusive_receipt_hash_bindings_'+case)
            elif case in ['backup_write_failure','backup_readback_failure']:
                check(error and not moves and not receipt.exists(),'backup_fault_stops_before_any_biological_move_'+case)
                check(all(original_sha(r['original_path'])==r['sha256'] for r in biological+sidecars),'all_original_fixture_bytes_retained_on_backup_fault_'+case)
            else:
                check(error and not receipt.exists() and first_move_backups==[True],'transport_or_final_fault_preserved_no_success_receipt_'+case)
                if case!='final_backup_hash_failure':check(all(original_sha(r['byte_backup_path'])==r['sha256'] for r in sidecars),'original284_literal_backups_retained_on_transport_fault_'+case)
            observed.append(dict(case=case,error=error,biological_moves=len(moves),all_backups_verified_before_first_move=first_move_backups,journal_event_count=len(events),success_receipt_exists=receipt.exists(),real_workers=0))
    for path,h in v3['dependencies_sha256'].items():check(common.sha(path)==h,'v3_frozen_dependency_after '+path)
    check(common.sha(PLAN)==PLAN_SHA and common.sha(SCRIPT)==SCRIPT_SHA and common.sha(PRIOR)==PRIOR_SHA,'frozen_versions_and_prior_controls_unchanged')
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':peak*=1024
    check(peak<128*(1<<20),'additive_stdlib_review_peak_RSS_below128MiB')
    record=dict(status='ADDITIVE_TRANSPORT_PRESERVATION_PREPARED_CODE_AND_TINY_CONTROLS_PASS',plan_sha256=PLAN_SHA,executor_sha256=SCRIPT_SHA,prior1400_control_sha256=PRIOR_SHA,prior_terminal_functions_reused_AST_identical=unchanged,checks=checks,check_count=len(checks),observations=observed,dependency_sha256=v3['dependencies_sha256'],sidecar_map_sha256=hashlib.sha256(json.dumps(v3['interrupted_sidecars'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),sidecar_count=284,biological_file_count=284,original_transport_backup_bytes=sum(r['bytes'] for r in v3['interrupted_sidecars']),frozen_success_receipt_map_unchanged=True,real_workers=0,original_failed_biological_or_sidecar_files_moved=0,GWAS_body_reads=0,estimator_calls=0,production_mutex_acquired=False,root_real_ExFAT_tiny_control_receipt_sha256=common.sha(P/'logs/native_sidecar_rename_control_v2.json'),scope='New sidecar/backup mappings and preservation function only;112/107 admission, cleanup and terminal controls reused from sealed1400check receipt by AST identity. Own tiny fixtures use explicit branch and mocked implicit branch on temporary files; root supplied genuine ExFAT control separately bound.',peak_RSS_bytes=peak,verifier_sha256=common.sha(__file__))
    with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(status=record['status'],checks=len(checks),prior_checks_reused=1400,peak_RSS_bytes=peak,receipt_sha256=common.sha(OUT),real_workers=0),indent=2))


if __name__=='__main__':main()
