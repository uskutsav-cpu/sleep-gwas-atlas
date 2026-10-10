#!/usr/bin/env python3
"""Independent metadata-only plan checks and tiny private controller/schema fixtures.

No actual source/reference body, scientific command, subprocess or estimator is used.
Controller tests replace all workers with explicitly non-scientific metadata stubs.
"""
import argparse
import collections
import copy
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
P = Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4')
PLAN = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/core_pipeline/original_small_input_replay_v1/core_original_small_replay_plan_v1.json')
PIN = '4b63e99da34f7ddcb05ab7109882c5caefea86d913c52ba417cf2e414e9d8b76'
FIX = P/'reviews/core_original_small_replay_controls_v1'
RECORDS = []
sys.path.insert(0, str(P/'scripts'))

def digest(q):return hashlib.sha256(Path(q).read_bytes()).hexdigest()
def new(q,v):
    with Path(q).open('x') as f:
        json.dump(v,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def module(q):
    s=importlib.util.spec_from_file_location('_private_small_review_'+str(len(RECORDS)),q)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def record(case,good,**detail):
    RECORDS.append(dict(case=case,control_pass=bool(good),**detail))
    assert good,RECORDS[-1]
def flag(cmd,key):return cmd[cmd.index(key)+1]

def metadata_checks():
    assert digest(PLAN)==PIN
    plan=json.loads(PLAN.read_text())
    design=json.loads((P/'reviews/independent_core_pipeline_replay_manifest_v4.json').read_text())
    rt=json.loads(Path(plan['qualified_runtime_receipt']).read_text())
    expected=[r for r in design['rows'] if r['prefilter'] is not None or r['trait_id']=='bmi']
    record('membership_original_order', [m['original_design'] for m in plan['members']]==expected
           and len(expected)==10 and len(design['rows'])==45,
           ordered_traits=[m['trait_id'] for m in plan['members']])
    expected_names=['ms','asthma','bmi','t2d','ldl','hdl','triglycerides','cad','telomere_length','melanoma']
    assert [m['trait_id'] for m in plan['members']]==expected_names
    rows=[]
    for m in plan['members']:
        d=m['original_design']; old=d['new_private_output_namespace']; newroot=str(PLAN.parent/m['trait_id'])
        transformed={}
        for label,cmd in d['command_templates'].items():
            transformed[label]=None if cmd is None else [rt['command'][0] if token=='{FROZEN_HARMONIZATION_PYTHON}'
                else newroot+token[len(old):] if token.startswith(old+'/') else token for token in cmd]
        record('command_exact_'+m['trait_id'], m['prefilter_command']==transformed['prefilter']
               and m['harmonize_command']==transformed['harmonize_original'] and m['munge_command']==transformed['munge_stock'])
        h=m['harmonize_command'];mu=m['munge_command']
        assert flag(h,'--trait')==m['trait_id'] and flag(h,'--source-build')==d['build']
        assert flag(h,'--config')==str(P.parent/'config/analysis_panel.tsv')
        assert flag(h,'--outdir')==str(Path(m['harmonized']).parent)
        assert flag(mu,'--sumstats')==m['harmonized'] and flag(mu,'--chunksize')=='500000' and flag(mu,'--out')==m['munged_prefix']
        assert d['raw']['historical_expected_sha256']==d['raw']['sealed_verified_sha256']
        qc=d['historical_QC'];rem=d['harmonizer_input_rows']
        for step in qc['ordered_steps']:
            rem-=step['dropped'];assert rem==step['remaining']
        assert rem==d['harmonized_rows']
        assert qc['metadata']['rows_in']==str(d['harmonizer_input_rows']) and qc['metadata']['rows_out']==str(rem)
        if d['prefilter'] is not None:
            assert d['prefilter']['sealed_sha256']==d['prefilter']['historical_expected_sha256']
            assert qc['metadata']['infile_sha256']==d['prefilter']['historical_expected_sha256']
            assert qc['metadata']['prefilter_source_sha256']==d['raw']['sealed_verified_sha256']
            assert qc['metadata']['prefilter_source_rows']==str(d['raw_source_rows'])
            assert qc['metadata']['prefilter_retained_rows']==str(d['harmonizer_input_rows'])
        rows.append(dict(trait=m['trait_id'],raw_bytes=d['raw']['bytes'],raw_source_rows=d['raw_source_rows'],
            harmonizer_input_rows=d['harmonizer_input_rows'],harmonized_rows=rem,
            munged_template_rows=d['historical_munge_log']['output_total_rows'],
            munged_nonmissing_rows=d['historical_munge_log']['output_nonmissing_rows'],
            prefilter_rerun=m['prefilter_command'] is not None,QC_step_count=len(qc['ordered_steps'])))
    record('scope_and_guard',sum(m['prefilter_command'] is not None for m in plan['members'])==6
           and plan['estimator_calls']==0 and plan['guard']['internal_floor_bytes']==3<<30
           and plan['guard']['SSD_floor_bytes']==5<<30 and plan['guard']['observed_aggregate_worker_RSS_limit_bytes']==2<<30
           and plan['guard']['new_output_limit_bytes']==16<<30 and plan['guard']['deadline_seconds']==96*3600
           and plan['guard']['per_worker_deadline_seconds']==2*3600 and plan['guard']['global_reservation_bytes']==300<<30)
    record('declared_runtime_qualification',rt['versions']['numpy']=='1.26.4' and rt['versions']['pandas']=='2.2.3'
           and rt['versions']['python'].startswith('3.11.11') and not rt['exact_historical_binary_recovery_claimed']
           and not rt['historical_per_trait_runtime_attestation_recovered'] and rt['regular_file_bytes']==533442552)
    ledger=json.loads(Path(plan['global_resource_ledger_path']).read_text())
    record('reservation_arithmetic',sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes']
           and ledger['ceiling_bytes']-ledger['reserved_total_bytes']==ledger['unallocated_margin_bytes']
           and ledger['reserved_total_bytes']<300<<30,
           reserved_total_bytes=ledger['reserved_total_bytes'],unallocated_margin_bytes=ledger['unallocated_margin_bytes'])
    checked={};excluded={}
    for path,wanted in plan['dependencies_sha256'].items():
        q=Path(path)
        if '/ref/' in path and not path.endswith('.provenance.json'):
            excluded[path]='reference body not opened';continue
        assert q.stat().st_size<=12*(1<<20),'metadata bound exceeded: '+path
        actual=digest(q);assert actual==wanted,path;checked[path]=actual
    record('metadata_dependencies',bool(checked),verified_count=len(checked),excluded_reference_count=len(excluded))
    for trait in ['ldl','hdl','triglycerides']:
        m=next(m for m in plan['members'] if m['trait_id']==trait);d=m['original_design']
        q=P.parent/'sleep_unified_research_v2/logs'/f'{trait}_native_prefilter_receipt_v2.json';pr=json.loads(q.read_text())
        prov=Path(flag(m['harmonize_command'],'--prefilter-provenance'));v=json.loads(prov.read_text())
        record('reused_lipid_proof_'+trait,pr['status']=='NATIVE_PREFILTER_EXACT_HISTORICAL_BYTES_AND_COUNTS'
            and pr['exit_code']==0 and pr['stop_reason'] is None and pr['input_hashes_unchanged']
            and pr['inputs_sha256_before']==pr['inputs_sha256_after']
            and pr['inputs_sha256_before']['raw_source']==d['raw']['sealed_verified_sha256']
            and pr['output_sha256']==d['prefilter']['historical_expected_sha256']==v['output_sha256']
            and pr['native_prefilter_provenance_sha256']==digest(prov)
            and v['source_rows']==d['raw_source_rows'] and v['retained_rows']==d['harmonizer_input_rows'])
    return dict(plan_sha256=PIN,members=rows,verified_metadata_code_sha256=checked,
                excluded_reference_bodies=excluded,original190_job_counts=dict(collections.Counter(j['stage'] for j in plan['original_190_jobs'])),
                actual_runtime_tree_rehashed_by_this_review=False)

def schema_controls():
    m=module(P/'scripts/72_validate_core_original_small_replay_v1.py')
    import extension_replay_common_v3 as common
    d=FIX/'schema';d.mkdir()
    def gz(name,payload,stamp=0):
        q=d/name;q.write_bytes(gzip.compress(payload,mtime=stamp));return q
    header=b'SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\n'
    lines=[b'fixture1\t1\t1\tA\tC\tNA\t0\t1\t1\t100\n',b'fixture2\t1\t2\tG\tA\t.2\t-0.0\t1\t1\t100\n']
    a=gz('harm_a.gz',header+b''.join(lines));b=gz('harm_b.gz',header+b''.join(lines),1)
    x=m.compare_harmonized(a,b,2)
    record('harm_exact_compressed_distinct',x['status']=='EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH'
           and x['compressed_sha256'][0]!=x['compressed_sha256'][1])
    for case,data,n in [('order',header+b''.join(reversed(lines)),2),('missing_literal',header+b''.join(lines).replace(b'NA',b'nan'),2),
                        ('numeric_literal',header+b''.join(lines).replace(b'-0.0',b'0.0'),2),('trailing_row',header+b''.join(lines)+lines[0],2),
                        ('expected_count',header+b''.join(lines),3)]:
        x=m.compare_harmonized(a,gz('harm_'+case+'.gz',data),n)
        record('harm_reject_'+case,x['status']=='HARMONIZED_CONTENT_MISMATCH_PRESERVED')
    bad=bytearray(b.read_bytes());bad[-8]^=1;(d/'harm_bad_crc.gz').write_bytes(bad)
    try:m.compare_harmonized(a,d/'harm_bad_crc.gz',2);caught=False
    except (gzip.BadGzipFile,EOFError):caught=True
    record('harm_CRC_failure',caught)
    head=b'SNP\tA1\tA2\tZ\tN\n';data=b'fixture1\tA\tC\t0.0\t100\nfixture2\tG\tA\tNA\tNA\n'
    a=gz('munge_a.gz',head+data);b=gz('munge_b.gz',head+data,1);x=common.compare_munged(a,b,2)
    record('munge_full_template_missing_count',x['status']=='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH'
           and x['counts_new']=={'finite_N_Z':1,'missing_or_nonfinite_N_Z':1})
    for case,payload,n in [('missing_literal',(head+data).replace(b'NA',b'nan'),2),('allele',(head+data).replace(b'fixture1\tA',b'fixture1\tT'),2),
                           ('numeric_literal',(head+data).replace(b'0.0',b'-0.0'),2),('expected_count',head+data,3)]:
        record('munge_reject_'+case,common.compare_munged(a,gz('munge_'+case+'.gz',payload),n)['status']=='DECOMPRESSED_TEMPLATE_MISMATCH_PRESERVED')
    bad=b.read_bytes()[:-3];(d/'munge_bad_eof.gz').write_bytes(bad)
    try:common.compare_munged(a,d/'munge_bad_eof.gz',2);caught=False
    except (gzip.BadGzipFile,EOFError):caught=True
    record('munge_EOF_failure',caught)
    q=d/'QC.txt';q.write_text('trait\tmetadata_fixture\n\nstep\tdropped\tremaining\nalpha\t2\t8\nbeta\t1\t7\n\nrows_out\t7\n')
    md,steps=m.qc(q)
    record('QC_order_parser',md=={'trait':'metadata_fixture','rows_out':'7'} and steps==[
        {'reason':'alpha','dropped':2,'remaining':8},{'reason':'beta','dropped':1,'remaining':7}])

def controller_control(case):
    d=FIX/('controller_'+case);d.mkdir();m=module(P/'scripts/71_run_core_original_small_replay_v1.py')
    m.SSD=d/'SSD';m.OUT=m.SSD/'core_pipeline/original_small_input_replay_v1';m.OUT.mkdir(parents=True)
    for sub in ['receipts_v4','logs_v4','tmp','cache']:(m.OUT/sub).mkdir()
    binding=d/'independent_review_metadata.json';new(binding,{'fixture_only':True})
    rtfile=d/'runtime_metadata.json';runtimefile=d/'runtime_file';runtimefile.write_text('metadata fixture runtime\n')
    new(rtfile,{'regular_files':{str(runtimefile):{'bytes':runtimefile.stat().st_size,'sha256':digest(runtimefile)}},'symlinks':{}})
    actual=json.loads(PLAN.read_text());plan={'guard':copy.deepcopy(actual['guard']),
        'qualified_runtime_receipt':str(rtfile),'qualified_runtime_receipt_sha256':digest(rtfile),
        'harmonization_python':'METADATA_STUB_NOT_AN_EXECUTABLE','dependencies_sha256':{str(binding):digest(binding)},
        'archived_input_sha256':{},'members':[]}
    plan['guard']['shared_heavy_worker_lock']=str(m.SSD/'native_heavy_worker.lock')
    for real in actual['members']:
        trait=real['trait_id'];td=m.OUT/trait;td.mkdir()
        for sub in ['harmonized','munged','prefilter','receipts']:(td/sub).mkdir()
        raw=td/'source_identity_metadata.json';new(raw,{'not_a_GWAS_body':True,'trait':trait})
        member={k:str(td/k) for k in ['harmonized','munged','source_gate_receipt','comparison_receipt','harmonization_qc']}
        member.update(trait_id=trait,munged_prefix=member['munged'],prefilter_command=None,
            harmonize_command=['METADATA_STUB','harmonize',member['harmonized']],
            munge_command=['METADATA_STUB','munge',member['munged']],
            original_design={'raw':{'resolved_path':str(raw),'sealed_verified_sha256':digest(raw)},'prefilter':None})
        plan['members'].append(member)
    planpath=d/'fixture_plan.json';new(planpath,plan);planhash=digest(planpath)
    admission=d/'fixture_admission.json';new(admission,{'execution_admitted':True,'plan_sha256':planhash,
        'executor_sha256':digest(P/'scripts/71_run_core_original_small_replay_v1.py'),
        'independent_review_sha256':{str(binding):digest(binding)},'metadata_fixture_only_no_actual_admission':True})
    a=SimpleNamespace(plan=planpath,plan_sha=planhash,admission=admission,admission_sha=digest(admission))
    class Stub:
        TERMINATION_REQUEST=[]
        def __init__(self):self.commands=[];self.cleanup_lock_verified=False;self.gates=[];self.postwrite=False
        def catchable_termination(self,signum,frame):self.TERMINATION_REQUEST.append(signum)
        def check_termination(self,where):
            if self.TERMINATION_REQUEST:raise RuntimeError('FIXTURE_DEFERRED_TERMINATION')
        def stage_resource_gate(self,pl,start,where):
            self.gates.append(where);self.check_termination(where)
            state,reason=self.limits(pl,start,0)
            if reason:raise RuntimeError(reason)
            if case=='posthash_resource' and where=='AFTER_CORE_OUTPUT_HASHES':raise RuntimeError('FIXTURE_POSTHASH_RESOURCE')
            if case=='postpersist_resource' and where=='POST_PERSISTENCE_CORE_TERMINAL':raise RuntimeError('FIXTURE_POSTPERSIST_RESOURCE')
            return state
        def baseline_gate(self,pl):return {'MOCK190_NOT_SCIENTIFIC_EVIDENCE':'metadata'}
        def safe_print(self,*args,**kwargs):pass
        def worker(self,command,prefix,receipt,pl,pp,ph,start,monitor,ownership):
            self.commands.append(command)
            if case in ['source_failure','owned_cleanup']:
                if case=='owned_cleanup':ownership[:]=[False,['METADATA_STUB_NO_PROCESS']]
                raise RuntimeError('FIXTURE_SOURCE_WORKER_FAILURE_NO_PROCESS')
            member=next(t for t in pl['members'] if str(prefix) in [t['source_gate_receipt'],t['comparison_receipt'],t['harmonized'],t['munged']])
            if str(prefix)==member['comparison_receipt']:
                new(prefix,{'status':'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS','fixture_only_not_scientific_evidence':True})
            else:new(prefix,{'fixture_only_not_scientific_evidence':True})
            if str(prefix)==member['harmonized']:new(member['harmonization_qc'],{'fixture_only_not_scientific_evidence':True})
            journal=Path(str(receipt)+'.failure_journal.jsonl');journal.write_text('{"fixture_only":true}\n')
            new(receipt,{'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':ph,'command':command,
                'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[]},
                'output_sha256':{str(prefix):digest(prefix),str(journal):digest(journal)}})
        def await_owned_cleanup(self,monitor,ownership,ph):
            fd=os.open(plan['guard']['shared_heavy_worker_lock'],os.O_RDWR)
            try:
                try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);self.cleanup_lock_verified=False
                except BlockingIOError:self.cleanup_lock_verified=True
            finally:os.close(fd)
            ownership[:]=[True,[]];return []
    stub=Stub();m.module=lambda name,path:stub if name=='_core_owned_worker' else SimpleNamespace()
    m.physical_mount=lambda:None
    class Disk:
        def disk_usage(self,path):return SimpleNamespace(free=1<<50)
    m.shutil=Disk();realwrite=m.write_new
    def fixture_write(path,value):
        if Path(path).name=='core_original_small_execution_receipt_v1.json':
            if case=='master_write_failure':Path(path).write_text('{"partial":');raise OSError('FIXTURE_MASTER_WRITE_FAILURE')
            realwrite(path,value)
            if case=='signal_during_master_write':stub.TERMINATION_REQUEST.append(signal.SIGTERM)
            if case=='output_mutation_during_master_write':Path(plan['members'][0]['harmonized']).write_text('MUTATED_METADATA_FIXTURE\n')
            return
        if case in ['master_write_failure','postpersist_resource','signal_during_master_write','output_mutation_during_master_write'] and str(path).endswith('.failure.json'):
            raise OSError('FIXTURE_SUPPLEMENTAL_WRITE_FAILURE')
        return realwrite(path,value)
    m.write_new=fixture_write
    # Module-local terminal class suppresses its supplemental writes for integration faults.
    if case in ['postpersist_resource','signal_during_master_write','output_mutation_during_master_write']:
        tm=module(P/'scripts/terminal_commit_common_v2.py');save=tm.save_new
        def termsave(path,value):
            if str(path).endswith('.failure.json'):raise OSError('FIXTURE_TERMINAL_SUPPLEMENTAL_WRITE_FAILURE')
            return save(path,value)
        tm.save_new=termsave;m.TerminalCommit=tm.TerminalCommit
    try:m.run(a);returned=True
    except SystemExit:returned=False
    import terminal_commit_common_v2 as consumer
    receipt=m.OUT/'core_original_small_execution_receipt_v1.json';pending=m.OUT/'core_pending_v1.json';seal=m.OUT/'core_terminal_seal_v1.json'
    bindings={'plan_sha256':planhash,'admission_sha256':a.admission_sha,'executor_sha256':digest(P/'scripts/71_run_core_original_small_replay_v1.py')}
    try:consumer.require_committed(pending,seal,bindings,{str(receipt):digest(receipt)});accepted=True
    except BaseException:accepted=False
    fd=os.open(plan['guard']['shared_heavy_worker_lock'],os.O_RDWR)
    try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);released=True
    except BlockingIOError:released=False
    finally:os.close(fd)
    expected=case=='success'
    record('controller_'+case,returned==expected and accepted==expected and released
        and (case!='owned_cleanup' or stub.cleanup_lock_verified),
        returned_success=returned,consumer_accepts=accepted,private_mutex_released=released,
        cleanup_private_mutex_held=stub.cleanup_lock_verified,metadata_stub_calls=len(stub.commands),
        actual_worker_calls=0,gate_sequence=stub.gates)

def main():
    FIX.mkdir(exist_ok=False);metadata=metadata_checks();schema_controls()
    for case in ['success','source_failure','owned_cleanup','posthash_resource','postpersist_resource',
                 'master_write_failure','signal_during_master_write','output_mutation_during_master_write']:
        controller_control(case)
    assert digest(PLAN)==PIN
    for path,h in metadata['verified_metadata_code_sha256'].items():assert digest(path)==h,path
    files={str(q):digest(q) for q in sorted(FIX.rglob('*')) if q.is_file() and not q.is_symlink()}
    out=P/'reviews/core_original_small_replay_controls_receipt_v1.json'
    new(out,{'schema':'independent_original_small_replay_metadata_controls_v1','plan_sha256':PIN,
        'control_script_sha256':digest(__file__),'control_count':len(RECORDS),'controls':RECORDS,
        'all_control_expectations_pass':all(q['control_pass'] for q in RECORDS),'metadata':metadata,
        'fixture_file_sha256':files,'GWAS_body_reads':0,'reference_body_reads':0,'real_data_decompressions':0,
        'tiny_synthetic_schema_gzip_fixtures_only':True,'actual_workers':0,'fits':0,
        'mock_controller_evidence_is_not_scientific_replay_evidence':True,'execution_admission_granted':False})
    print(json.dumps({'control_count':len(RECORDS),'receipt_sha256':digest(out),'actual_workers':0,'fits':0},indent=2))

if __name__=='__main__':main()
