#!/usr/bin/env python3
"""Bounded metadata and isolated software controls; never run replay workers."""
import argparse
import ast
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import resource
import signal
import sys
import tempfile
from types import SimpleNamespace

from extension_replay_common import sha,write_new,check_bindings,compare_munged,acquisition_receipt_gate

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v1'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan-sha',required=True);a=parser.parse_args()
    plan_path=OUT/'extension_pipeline_replay_plan_v1.json';plan=json.loads(plan_path.read_text());checks=[]
    def check(value,label):
        if not value:raise RuntimeError('SOFTWARE_CONTROL_FAILED: '+label)
        checks.append(label)
    check(sha(plan_path)==a.plan_sha,'exact_candidate_plan_SHA')
    check_bindings(plan,include_archived=True);check(True,'all_frozen_code_docs_receipt_reference_archived_processed_hashes')
    acq=json.loads(Path(plan['acquisition_plan']).read_text())
    check(len(plan['members'])==plan['member_count']==100,'exact100_panel_sources')
    check(plan['estimator_calls']==0 and plan['heavy_preprocessing_commands']==200 and plan['owned_stdlib_validation_commands']==200 and plan['total_owned_commands']==400,'exact400_owned_commands_zero_estimators')
    check(plan['reservation_arithmetic']['total_reserved_component_bytes']<300*(1<<30),'reservation_components_within_frozen300GiB')
    check(plan['guard']['internal_floor_bytes']==3*(1<<30) and plan['guard']['SSD_floor_bytes']==5*(1<<30) and plan['guard']['observed_aggregate_worker_RSS_limit_bytes']==2*(1<<30),'original3_5_2GiB_floors_RSS')
    check(plan['guard']['new_output_limit_bytes']==32*(1<<30) and plan['guard']['worker_count']==plan['guard']['BLAS_threads']==1,'32GiB_pipeline_cap_one_worker_BLAS1')
    check(plan['guard']['shared_heavy_worker_lock']==str(SSD/'native_heavy_worker.lock'),'exact_shared_exclusive_mutex_path')
    historical={m['extension_trait_id']:m for m in acq['members']}
    for m in plan['members']:
        t=m['extension_trait_id'];old=historical[t]
        check(m['acquisition_member']==old,'frozen_acquisition_member '+t)
        check((m['raw_sha256'],m['raw_md5'],m['raw_bytes'],m['raw'])==(old['expected_sha256'],old['expected_md5'],old['expected_bytes'],old['body_path']),'exact_original_source_hash_MD5_bytes_path '+t)
        expected=[plan['python'],'-u','-B',str(P.parent/'discovery_extension/scripts/10_harmonize_panukbb.py'),'--trait-id',t,'--source',m['raw'],'--panel',plan['panel'],'--reference',plan['reference'],'--out',m['harmonized'],'--qc-out',m['harmonization_qc'],'--receipt-out',m['harmonization_receipt']]
        check(m['harmonize_command']==expected,'exact_original_local_harmonizer_argv '+t)
        expected=[plan['python'],'-u','-B',str(Path(plan['ldsc_dir'])/'munge_sumstats.py'),'--sumstats',m['harmonized'],'--merge-alleles',plan['w_hm3'],'--snp','SNP','--a1','A1','--a2','A2','--frq','FRQ','--p','P','--N-col','N','--signed-sumstats','BETA,0','--chunksize','500000','--out',m['munged_prefix']]
        check(m['munge_command']==expected,'exact_original_munge_argv '+t)
        check(plan['archived_input_sha256'][m['archived_munged']]==m['archived_munged_sha256'],'exact_archival_processed_comparison_source '+t)
    names=['46_prepare_extension_pipeline_replay.py','47_run_extension_pipeline_replay.py','extension_replay_common.py','extension_replay_validate.py','extension_replay_fault_controls.py']
    for name in names:ast.parse((P/'scripts'/name).read_text(),filename=name);check(True,'Python_AST '+name)
    spec=importlib.util.spec_from_file_location('_isolated_pipeline_software_controls',P/'scripts/47_run_extension_pipeline_replay.py')
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    protected=runner.load_protected()
    saved=dict(out=runner.OUT,write=runner.write_new,argv=sys.argv,run=runner.run)
    try:
        with tempfile.TemporaryDirectory(prefix='extension_replay_controls_',dir=OUT/'tmp') as directory:
            tmp=Path(directory);runner.OUT=tmp
            # This is a labelled parser/software fixture, never scientific data.
            payload=b'SNP\tA1\tA2\tZ\tN\nCONTROL_ONLY_1\tA\tG\t1.234\t100.000\nCONTROL_ONLY_2\tC\tT\t\t\n'
            def gz(name,body,mtime):
                path=tmp/name
                with path.open('xb') as f:
                    with gzip.GzipFile(fileobj=f,filename=name,mode='wb',mtime=mtime) as z:z.write(body)
                return path
            archived=gz('archived_fixture.gz',payload,0);same=gz('new_fixture.gz',payload,1)
            r=compare_munged(same,archived,2)
            check(r['status']=='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' and r['compressed_sha256_new']!=r['compressed_sha256_archived'],'full_decompressed_identity_admitted_despite_gzip_serialization_difference')
            check(r['both_gzip_CRC_and_EOF_verified'] and r['counts_new']=={'finite_N_Z':1,'missing_or_nonfinite_N_Z':1},'entire_EOF_CRC_and_missing_N_Z_counts')
            for label,body in [
                ('SNP',payload.replace(b'CONTROL_ONLY_1',b'CONTROL_ONLY_CHANGED')),
                ('allele',payload.replace(b'\tA\tG\t',b'\tG\tA\t')),
                ('Z',payload.replace(b'1.234',b'1.235')),
                ('N',payload.replace(b'100.000',b'101.000')),
                ('literal_missingness',payload.replace(b'\t\t\n',b'\tNA\tNA\n')),
                ('order',payload.splitlines(keepends=True)[0]+b''.join(reversed(payload.splitlines(keepends=True)[1:]))),
                ('trailing_row',payload+b'CONTROL_ONLY_3\tA\tG\t2.000\t100.000\n'),
            ]:
                changed=gz(label+'.gz',body,0);r=compare_munged(changed,archived,2)
                check(r['status']=='DECOMPRESSED_TEMPLATE_MISMATCH_PRESERVED' and r['both_gzip_CRC_and_EOF_verified'],'reject_entire_stream_change_'+label)
            damaged=tmp/'CRC_damage.gz';bad=bytearray(same.read_bytes());bad[-8]^=1;damaged.write_bytes(bad)
            truncated=tmp/'truncated.gz';truncated.write_bytes(same.read_bytes()[:-4])
            for name,path in [('CRC',damaged),('truncated',truncated)]:
                try:compare_munged(path,archived,2)
                except (OSError,EOFError):check(True,'reject_'+name+'_without_false_EOF_proof')
                else:raise RuntimeError('corrupt gzip accepted')
            # Exact acquisition-gate metadata witness uses a 7-byte local fixture.
            # Its names/hashes are explicitly unrelated to GWAS or biological data.
            body=tmp/'source_fixture';body.write_bytes(b'control')
            import hashlib
            member=dict(index=2,extension_trait_id='SOFTWARE_ONLY',filename=body.name,body_path=str(body),expected_bytes=7,expected_md5=hashlib.md5(b'control').hexdigest(),expected_sha256=sha(body),expected_s3_version_id='CONTROL_VERSION',expected_etag='CONTROL_ETAG',url='https://invalid.example/SOFTWARE_ONLY')
            fixture_plan={'guard':plan['guard'],'acquisition_plan_sha256':'CONTROL_PLAN','acquisition_operational_identity':{'source_folder':str(tmp),'curl':'/usr/bin/curl','per_body_seconds_limit':7200,'resume_prefix_bytes':1}}
            receipt=tmp/'source_fixture_receipt.json';m={'acquisition_member':member,'acquisition_receipt':str(receipt)}
            command=['/usr/bin/curl','-q','--fail','--location','--max-redirs','5','--proto','=https','--proto-redir','=https','--tlsv1.2','--max-time','7200','--max-filesize','7','--dump-header',str(tmp/'logs/SOFTWARE_ONLY')+'.headers.txt','--output',str(body)+'.partial',member['url']]
            valid={'member':member,'plan_sha256':'CONTROL_PLAN','status':'EXACT_IMMUTABLE_SOURCE_ACQUIRED','returncode':0,'stop_reason':None,'teardown':{'teardown_verified':True,'remaining_group_members':[]},'post_cleanup_hash_resource_identity_gates_pass':True,'actual_size':7,'actual_md5':member['expected_md5'],'actual_sha256':member['expected_sha256'],'final_seal_md5':member['expected_md5'],'final_seal_sha256':member['expected_sha256'],'resume_offset':0,'command':command,'observed_headers':{'status':200,'x-amz-version-id':'CONTROL_VERSION','etag':'CONTROL_ETAG','content-length':'7'},'resource_before':{'internal_free_bytes':3*(1<<30),'ssd_free_bytes':5*(1<<30)},'resource_after':{'internal_free_bytes':3*(1<<30),'ssd_free_bytes':5*(1<<30)}}
            write_new(receipt,valid);acquisition_receipt_gate(fixture_plan,m,sha(receipt));check(True,'exact_successful_acquisition_metadata_fixture_admitted')
            for key,value in [('status','FAILED_PRESERVED_NO_AUTOMATIC_RETRY'),('actual_sha256','0'*64),('final_seal_sha256','0'*64),('plan_sha256','WRONG_PLAN'),('post_cleanup_hash_resource_identity_gates_pass',False)]:
                altered=copy.deepcopy(valid);altered[key]=value;path=tmp/(key+'.json');write_new(path,altered);bad=dict(m,acquisition_receipt=str(path))
                try:acquisition_receipt_gate(fixture_plan,bad,sha(path))
                except RuntimeError:check(True,'reject_acquisition_'+key)
                else:raise RuntimeError('uncleared acquisition accepted')
            root_fixture={'members':[{'acquisition_receipt':str(receipt)}]};admission={'execution_admitted':True,'plan_sha256':'P','executor_sha256':'E','acquisition_receipt_sha256':{str(receipt):sha(receipt)},'independent_review_sha256':{str(receipt):sha(receipt)}}
            runner.admission_gate(root_fixture,admission,'P','E');check(True,'explicit_root_admission_hash_bindings')
            for key,value in [('execution_admitted',False),('plan_sha256','WRONG'),('executor_sha256','WRONG'),('acquisition_receipt_sha256',{}),('independent_review_sha256',{})]:
                bad=copy.deepcopy(admission);bad[key]=value
                try:runner.admission_gate(root_fixture,bad,'P','E')
                except RuntimeError:check(True,'reject_root_admission_'+key)
                else:raise RuntimeError('invalid root admission accepted')
            invocations=[];runner.run=lambda _a:invocations.append(True);sys.argv=['software_driver']
            try:runner.main()
            except SystemExit as e:check('Prepared only' in str(e),'no_execute_flag_blocks_before_dispatch')
            check(not invocations,'no_real_worker_or_driver_dispatch_without_execute_flag')
            # Final persistence/late-signal/deadline controls use isolated files.
            protected.limits=lambda *_a:({'software_control':True},None)
            for label,write_failure,late_signal,deadline in [('healthy',False,False,False),('receipt_failure',True,False,False),('postseal_signal',False,True,False),('final_resource_failure',False,False,True)]:
                protected.TERMINATION_REQUEST.clear();target=tmp/(label+'.json');master={'status':'ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS'}
                protected.limits=lambda *_a:({'software_control':True},'STAGE_DEADLINE' if deadline else None)
                def controlled_write(path,data):
                    if Path(path)==target and write_failure:raise OSError('CONTROLLED_PRIMARY_STORAGE_FAILURE')
                    saved['write'](path,data)
                    if Path(path)==target and late_signal:protected.catchable_termination(signal.SIGTERM,None)
                runner.write_new=controlled_write
                runner.seal_master(protected,master,target,{},0,[True,[]])
                check((master['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS')==(label=='healthy'),'final_seal_status_'+label)
                if label!='healthy':check(Path(str(target)+'.failure.json').exists(),'failure_seal_preserved_'+label)
            check(not protected.TERMINATION_REQUEST,'no_pending_signal_after_control_reset')
    finally:
        runner.OUT=saved['out'];runner.write_new=saved['write'];sys.argv=saved['argv'];runner.run=saved['run'];protected.TERMINATION_REQUEST.clear()
    check_bindings(plan,include_archived=True);check(sha(plan_path)==a.plan_sha,'all_bindings_and_plan_unchanged_after_controls')
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    check(rss<128*(1<<20),'stdlib_control_RSS_under128MiB')
    record=dict(status='RESULT_FREE_BINDINGS_AND_SOFTWARE_FAULT_CONTROLS_PASS',plan_sha256=a.plan_sha,verifier_sha256=sha(__file__),checks=checks,check_count=len(checks),self_peak_RSS_bytes=rss,real_worker_subprocesses_launched=0,estimator_calls=0,real_harmonization_munging_or_merge_audits=0,raw_GWAS_body_reads=0,scope='Frozen metadata/hashes and isolated tiny parser/failure controls; not source replay or scientific validation.')
    out=OUT/'preparation_controls_v1.json';write_new(out,record)
    print(json.dumps(dict(status=record['status'],check_count=len(checks),sha256=sha(out),self_peak_RSS_bytes=rss,real_workers_launched=0),indent=2))


if __name__=='__main__':main()
