#!/usr/bin/env python3
"""Independent complete core10-v2 literal-content, QC and execution review.
Uses stdlib bounded streams, not the production comparison helper. No raw GWAS
body, reference body, native estimator, or runtime-tree rehash is performed.
"""
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import time

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
ROOT=SSD/'core_pipeline/original_small_input_replay_v2'
PLAN=ROOT/'core_original_small_replay_plan_v2.json'
MASTER=ROOT/'core_original_small_execution_receipt_v2.json'
TERMINAL=ROOT/'core_terminal_seal_v2.json'
REVIEW_PLAN=P/'reviews/independent_core_small_v2_whole_plan_v1.json'
OUT=P/'reviews/independent_core_small_v2_whole_receipt_v1.json'
T0=time.monotonic();B={};ROWS=[];WORKERS=[];OUTPUTS={};PREFILTERS=[]
GUARD=json.loads(REVIEW_PLAN.read_text())
RUNTIME_EXECUTABLES=set()
PROGRESS=P/'reviews/independent_core_small_v2_whole_progress_v1.jsonl'

def budget():
    assert time.monotonic()-T0<GUARD['deadline_seconds']
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<GUARD['maximum_observed_review_RSS_bytes']
    assert shutil.disk_usage('/System/Volumes/Data').free>=GUARD['internal_emergency_floor_bytes']
    assert shutil.disk_usage(SSD).free>=GUARD['SSD_emergency_floor_bytes']

def sha(path):
    h=hashlib.sha256();n=0
    with Path(path).open('rb') as f:
        for x in iter(lambda:f.read(65536),b''):
            h.update(x);n+=len(x)
            if n%(16*2**20)<65536:budget()
    return h.hexdigest()

def bind(path,expected=None):
    path=Path(path)
    assert path.is_file() and (not path.is_symlink() or str(path) in RUNTIME_EXECUTABLES),('regular file required',str(path))
    got=sha(path);assert expected is None or got==expected,(str(path),got,expected)
    B[str(path)]={'sha256':got,'bytes':path.stat().st_size}
    if path.is_symlink():B[str(path)].update(qualified_runtime_executable_symlink=True,literal_target=os.readlink(path),resolved_target=str(path.resolve()))
    return got

def read(path,expected=None):
    bind(path,expected);return json.loads(Path(path).read_text())

def parse_qc(path):
    lines=Path(path).read_text().splitlines();meta={};steps=[];in_steps=False
    for line in lines:
        if not line:continue
        fields=line.split('\t')
        if fields==['step','dropped','remaining']:in_steps=True;continue
        if len(fields)==3 and in_steps:
            steps.append({'reason':fields[0],'dropped':int(fields[1]),'remaining':int(fields[2])})
        elif len(fields)==2:
            in_steps=False;assert fields[0] not in meta,('duplicate QC key',path,fields[0]);meta[fields[0]]=fields[1]
        else:raise AssertionError(('QC width',path))
    remaining=int(meta['rows_in'])
    for step in steps:
        assert step['dropped']>=0;remaining-=step['dropped'];assert remaining==step['remaining']
    assert remaining==int(meta['rows_out'])
    return meta,steps

def meaningful_munge_log(path):
    # Capture every scientific/read/filter/count/summary line in original order;
    # only file path and timing lines are outside this exact equality check.
    text=Path(path).read_text();selected=[]
    for line in text.splitlines():
        if (line.startswith(('Removed ','Median value of ','Mean chi^2 = ','Lambda GC = ','Max chi^2 = ','WARNING:'))
            or re.fullmatch(r'Read \d+ SNPs (for allele merge\.|from --sumstats file\.)',line)
            or re.fullmatch(r'\d+ SNPs remain\.',line)
            or re.fullmatch(r'\d+ Genome-wide significant SNPs \(some may have been removed by filtering\)\.',line)):
            selected.append(line)
        elif line.startswith('Writing summary statistics for '):
            m=re.match(r'Writing summary statistics for (\d+) SNPs \((\d+) with nonmissing beta\) to ',line)
            assert m;selected.append('Writing summary statistics for '+m[1]+' SNPs ('+m[2]+' with nonmissing beta)')
    assert selected;return selected

def full_literal_pair(new,old,kind,expected):
    # Equal-sized decoded blocks establish entire literal ordered stream identity.
    # Parse only one equal stream, without storing source-row excerpts or sets.
    h=[hashlib.sha256(),hashlib.sha256()];sizes=[0,0];carry=b'';header=None;count=0
    finite=missing=0;nonpositive_N=0;empty_alleles=0
    with gzip.open(new,'rb') as a,gzip.open(old,'rb') as b:
        while True:
            blocks=[a.read(65536),b.read(65536)]
            if not any(blocks):break
            assert blocks[0]==blocks[1],('literal ordered stream mismatch',str(new),str(old),count)
            for i,x in enumerate(blocks):h[i].update(x);sizes[i]+=len(x)
            chunks=(carry+blocks[0]).split(b'\n');carry=chunks.pop();assert len(carry)<65536
            for line in chunks:
                fields=line.rstrip(b'\r').split(b'\t')
                if header is None:
                    header=fields
                    wanted=(b'SNP CHR BP A1 A2 FRQ BETA SE P N'.split() if kind=='harmonized' else b'SNP A1 A2 Z N'.split())
                    assert header==wanted,('literal canonical header',new)
                    continue
                assert len(fields)==len(header),('row width',new,count)
                count+=1
                if kind=='munged':
                    try:values=[float(fields[3]),float(fields[4])];ok=all(math.isfinite(v) for v in values)
                    except ValueError:ok=False
                    if ok:
                        finite+=1;nonpositive_N+=values[1]<=0;empty_alleles+=not fields[1] or not fields[2]
                    else:missing+=1
            budget()
    assert not carry,'missing final newline in decompressed stream'
    assert count==expected and h[0].hexdigest()==h[1].hexdigest()
    assert nonpositive_N==empty_alleles==0
    return dict(kind=kind,rows=count,header=[x.decode('ascii') for x in header],
                full_decompressed_bytes=sizes[0],full_decompressed_sha256=h[0].hexdigest(),
                ordered_literal_SNP_allele_N_Z_or_all_harmonized_fields_match=True,
                both_gzip_CRC_and_EOF_verified=True,finite_N_Z=finite if kind=='munged' else None,
                missing_or_nonfinite_N_Z=missing if kind=='munged' else None)

def expected_jobs(member,plan):
    trait=member['trait_id']
    def validator(mode,key):return [plan['harmonization_python'],'-B',str(P/'scripts/72_validate_core_original_small_replay_v2.py'),
        '--plan',str(PLAN),'--plan-sha',GUARD['production_plan_sha256'],'--trait',trait,'--mode',mode,'--out',member[key]]
    jobs=[('source',validator('source','source_gate_receipt'))]
    if member['prefilter_command'] is not None:jobs.append(('prefilter',member['prefilter_command']))
    jobs += [('harmonize',member['harmonize_command']),('munge',member['munge_command']),('compare',validator('compare','comparison_receipt'))]
    return jobs

def main():
    assert not OUT.exists() and not PROGRESS.exists();assert shutil.disk_usage('/System/Volumes/Data').free>=GUARD['internal_launch_floor_bytes']
    progress=PROGRESS.open('x')
    bind(REVIEW_PLAN)
    plan=read(PLAN,GUARD['production_plan_sha256']);master=read(MASTER,GUARD['master_sha256']);seal=read(TERMINAL,GUARD['terminal_sha256'])
    RUNTIME_EXECUTABLES.update([plan['python'],plan['harmonization_python']])
    assert master['status']=='ALL10_ORIGINAL_CORE_HARMONIZED_MUNGED_CONTENT_QC_REPLAY_PASS'
    assert master['estimator_calls']==0 and master['termination_requests']==[] and master['owned_cleanup_verified']
    assert len(master['completed_members'])==10 and len(master['worker_receipt_sha256'])==46
    executor=P/'scripts/71_run_core_original_small_replay_v2.py';bind(executor)
    admit_path=P/'manifests/core_original_small_replay_admission_v2.json';admit=read(admit_path,master['admission_sha256'])
    assert admit['execution_admitted'] and admit['plan_sha256']==sha(PLAN) and admit['executor_sha256']==sha(executor)
    for path,digest in admit['independent_review_sha256'].items():bind(path,digest)
    assert seal['status']=='REVIEWED_STAGE_TERMINAL_SEAL'
    assert seal['binding']==dict(plan_sha256=sha(PLAN),admission_sha256=sha(admit_path),executor_sha256=sha(executor))
    assert seal['result_receipt_sha256']=={str(MASTER):sha(MASTER)}
    assert seal['pending_path']==str(ROOT/'core_pending_v2.json') and seal['success_requires_absent_PENDING_and_all_failure_addenda']
    assert not (ROOT/'core_pending_v2.json').exists() and not (ROOT/'core_pending_v2.json').is_symlink()
    for p in [MASTER,TERMINAL]:assert not Path(str(p)+'.failure.json').exists() and not Path(str(p)+'.failure.json').is_symlink()
    oldfailure=read(plan['preserved_failed_v1_receipt'],plan['preserved_failed_v1_receipt_sha256'])
    assert oldfailure['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and oldfailure['owned_cleanup_verified']
    oldpending=Path(plan['preserved_failed_v1_receipt']).parent/'core_pending_v1.json'
    bind(oldpending,plan['dependencies_sha256'][str(oldpending)])
    for path,digest in master['historical190_gate_receipts'].items():bind(path,digest)
    assert len(plan['original_190_jobs'])==190
    runtime=read(plan['qualified_runtime_receipt'],plan['qualified_runtime_receipt_sha256'])
    assert master['runtime_before']==dict(runtime_receipt_sha256=sha(plan['qualified_runtime_receipt']),
        regular_files_verified=len(runtime['regular_files']),literal_symlinks_verified=len(runtime['symlinks']),historical_per_trait_binary_attestation=False)
    del runtime
    excluded=[]
    for path,digest in plan['dependencies_sha256'].items():
        f=Path(path)
        if f.name in ['hm3_grch37_variant_map.tsv.gz','hg38ToHg19.over.chain.gz','w_hm3.snplist']:
            excluded.append(dict(path=path,expected_sha256=digest));continue
        bind(f,digest)
    worker_paths=set();peaks=[];end_internal=[]
    for member,completed in zip(plan['members'],master['completed_members']):
        trait=member['trait_id'];d=member['original_design'];assert completed['trait']==trait
        for label,command in expected_jobs(member,plan):
            receipt=ROOT/'receipts_v4'/(trait+'__'+label+'.worker.json');worker_paths.add(str(receipt))
            r=read(receipt,master['worker_receipt_sha256'][str(receipt)])
            assert r['status']=='WORKER_COMPLETE_VERIFIED' and r['command']==command and r['plan_sha256']==sha(PLAN)
            assert r['returncode']==0 and r['stop_reason'] is None and r['plan_unchanged'] and r['owned_cleanup_verified']
            assert not r['metadata_errors'] and r['process_group_teardown']['remaining_group_members']==[]
            assert not Path(str(receipt)+'.failure.json').exists()
            final=r['resource_final'];g=plan['guard']
            assert final['internal_free_bytes']>=g['internal_floor_bytes'] and final['SSD_free_bytes']>=g['SSD_floor_bytes']
            assert final['observed_worker_RSS_bytes']<=g['observed_aggregate_worker_RSS_limit_bytes']
            assert final['core_namespace_bytes']<=g['new_output_limit_bytes'] and final['new_campaign_namespace_bytes']<=g['global_reservation_bytes']
            assert final['elapsed_stage_seconds']<=g['deadline_seconds'] and final['elapsed_worker_seconds']<=g['per_worker_deadline_seconds']
            peaks.append(r['peak_observed_owned_RSS_bytes']);end_internal.append(final['internal_free_bytes'])
            for path,digest in r['output_sha256'].items():
                if path in OUTPUTS:assert OUTPUTS[path]==digest
                OUTPUTS[path]=digest
            WORKERS.append(dict(trait=trait,label=label,receipt=str(receipt),receipt_sha256=sha(receipt),
                                exact_command_match=True,verified_empty_owned_group=True,peak_observed_RSS_bytes=peaks[-1]))
        source=read(member['source_gate_receipt'],completed['output_sha256'][member['source_gate_receipt']])
        assert source['status']=='EXACT_CORE_RAW_SOURCE_GATE_PASS' and source['trait']==trait and source['plan_sha256']==sha(PLAN)
        assert source['raw_sha256']==d['raw']['sealed_verified_sha256']==d['raw']['historical_expected_sha256'] and source['raw_bytes']==d['raw']['bytes']
        comparison=read(member['comparison_receipt'],completed['output_sha256'][member['comparison_receipt']])
        assert comparison['status']=='EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS' and comparison['plan_sha256']==sha(PLAN)
        assert comparison['raw_sha256']==source['raw_sha256'] and comparison['raw_bytes']==source['raw_bytes']
        pairs={}
        for kind in ['harmonized','munged']:
            new=Path(member[kind]);old=Path(d[kind]['path'])
            newsha=bind(new,completed['output_sha256'][str(new)]);oldsha=bind(old,d[kind]['sealed_sha256'])
            expected=d['harmonized_rows'] if kind=='harmonized' else d['historical_munge_log']['output_total_rows']
            pairs[kind]=full_literal_pair(new,old,kind,expected)
            pairs[kind].update(compressed_sha256_new=newsha,compressed_sha256_archived=oldsha,compressed_byte_identical=newsha==oldsha)
        assert pairs['munged']['finite_N_Z']==d['historical_munge_log']['output_nonmissing_rows']
        qc_old=Path(d['historical_QC']['path']);qc_new=Path(member['harmonization_qc'])
        bind(qc_old,d['historical_QC']['sha256']);bind(qc_new,completed['output_sha256'][str(qc_new)])
        om,oldsteps=parse_qc(qc_old);nm,newsteps=parse_qc(qc_new)
        assert om==d['historical_QC']['metadata'] and oldsteps==d['historical_QC']['ordered_steps'] and newsteps==oldsteps
        pathkeys={'infile','liftover_chain','variant_map'}
        assert all(nm[k]==v for k,v in om.items() if k not in pathkeys)
        assert nm['infile']==member['harmonize_command'][member['harmonize_command'].index('--infile')+1]
        assert int(nm['rows_out'])==pairs['harmonized']['rows'] and int(nm['rows_in'])==d['harmonizer_input_rows']
        pf={}
        if d['prefilter'] is not None:
            hcmd=member['harmonize_command'];current=Path(hcmd[hcmd.index('--infile')+1]);proof=Path(hcmd[hcmd.index('--prefilter-provenance')+1])
            pfr=read(proof,completed['consumed_prefilter_sha256'][str(proof)])
            bind(current,completed['consumed_prefilter_sha256'][str(current)])
            bind(d['prefilter']['path'],d['prefilter']['historical_expected_sha256'])
            assert sha(current)==d['prefilter']['historical_expected_sha256'] and pfr['output_sha256']==sha(current)
            assert pfr['source_rows']==d['raw_source_rows'] and pfr['retained_rows']==d['harmonizer_input_rows']
            assert pfr['input_sha256']==source['raw_sha256'] and pfr['input_bytes']==source['raw_bytes']
            assert pfr['allowlist_sha256']==nm['prefilter_allowlist_sha256'] and pfr['strategy']=='HAPMAP3_RSID_ALLOWLIST'
            pf=dict(trait=trait,replayed_now=member['prefilter_command'] is not None,
                    exact_historical_compressed_sha256=sha(current),source_rows=pfr['source_rows'],retained_rows=pfr['retained_rows'])
            PREFILTERS.append(pf)
        logold=Path(d['historical_munge_log']['path']);lognew=Path(member['munged_prefix']+'.log')
        bind(logold,d['historical_munge_log']['sha256']);bind(lognew)
        old_lines=meaningful_munge_log(logold);new_lines=meaningful_munge_log(lognew)
        assert new_lines==old_lines,('scientific munger log difference',trait)
        summary=dict(trait=trait,raw_bytes=source['raw_bytes'],raw_source_sha256_as_verified_by_executed_gates=source['raw_sha256'],
                     independent_raw_body_rehash=False,harmonizer_input_rows=d['harmonizer_input_rows'],
                     harmonized_rows=pairs['harmonized']['rows'],munged_template_rows=pairs['munged']['rows'],
                     munged_finite_N_Z=pairs['munged']['finite_N_Z'],munged_missing_N_Z=pairs['munged']['missing_or_nonfinite_N_Z'],
                     original_ordered_QC_steps=len(oldsteps),original_all_QC_steps_and_nonpath_metadata_match=True,
                     additional_QC_metadata_fields=sorted(set(nm)-set(om)),stock_munge_scientific_log_lines_match=len(old_lines),
                     compressed_harmonized_identical=pairs['harmonized']['compressed_byte_identical'],
                     compressed_munged_identical=pairs['munged']['compressed_byte_identical'],comparisons=pairs,prefilter=pf)
        ROWS.append(summary)
        progress.write(json.dumps(summary,sort_keys=True)+'\n');progress.flush();os.fsync(progress.fileno())
        print(json.dumps(dict(independent_complete_trait=trait,harmonized=summary['harmonized_rows'],munged_finite=summary['munged_finite_N_Z'])),flush=True)
    assert set(master['worker_receipt_sha256'])==worker_paths and len(worker_paths)==46
    progress.close()
    actual={str(p) for p in (ROOT/'receipts_v4').glob('*.worker.json') if not p.name.startswith('._')}
    assert actual==worker_paths
    for path,digest in OUTPUTS.items():bind(path,digest)
    for completed in master['completed_members']:
        for path,digest in completed['output_sha256'].items():bind(path,digest)
    for path,meta in B.items():
        assert sha(path)==meta['sha256'],('after-review file mutation',path)
        if meta.get('qualified_runtime_executable_symlink'):
            assert Path(path).is_symlink() and os.readlink(path)==meta['literal_target'] and str(Path(path).resolve())==meta['resolved_target']
    assert not (ROOT/'core_pending_v2.json').exists()
    budget()
    table=P/'reviews/independent_core_small_v2_whole_results_v1.tsv'
    fields=['trait','raw_bytes','raw_source_sha256_as_verified_by_executed_gates','independent_raw_body_rehash',
            'harmonizer_input_rows','harmonized_rows','munged_template_rows','munged_finite_N_Z','munged_missing_N_Z',
            'original_ordered_QC_steps','original_all_QC_steps_and_nonpath_metadata_match',
            'stock_munge_scientific_log_lines_match','compressed_harmonized_identical','compressed_munged_identical']
    with table.open('x') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(ROWS)
    result=dict(schema='independent_actual_core_small_v2_whole_content_QC_execution_v1',
                verdict='ALL10_ORDERED_HARMONIZED_MUNGED_CONTENT_QC_AND46_COMMANDS_INDEPENDENT_PASS',
                plan_sha256=sha(PLAN),master_sha256=sha(MASTER),terminal_sha256=sha(TERMINAL),
                rows=ROWS,workers=WORKERS,trait_count=10,command_count=46,prefilter_replayed=6,prefilter_reused=3,
                prefilt_records=PREFILTERS,all20_full_literal_decompressed_pairs_match=True,
                all20_pairs_both_gzip_CRC_and_EOF_verified=True,
                row_table_sha256=sha(table),code_metadata_and_derivative_before_after_sha256=B,
                aggregate_progress_sha256=sha(PROGRESS),
                reviewer_harness_corrections_preserved=['independent_core_small_v2_whole_harness_correction_v1.json','independent_core_small_v2_whole_harness_correction_v1_2.json'],
                unread_reference_dependencies=excluded,raw_GWAS_body_reads=0,fresh_raw_body_hashes_by_this_reviewer=False,
                actual_workers_or_estimators_launched_by_reviewer=0,production_stage_estimator_calls=0,
                original190_gate_metadata_count=len(master['historical190_gate_receipts']),
                complete_original190_numerical_reaudited=False,qualified_runtime_tree_rehashed_by_this_reviewer=False,
                exact_historical_per_trait_binaries_attested=False,compressed_byte_identity_claimed=False,
                old_failed_v1_preserved_with_PENDING=True,success_requires_exact_terminal_consumer=True,
                core_family_size=45,small_input_traits_replayed=10,other35_raw_chains_not_certified_by_this_review=True,
                fully_independent_two_trait_replication_claimed=False,scientific_rules_or_thresholds_modified=False,
                peak_observed_production_worker_RSS_bytes=max(peaks),minimum_production_worker_final_internal_free_bytes=min(end_internal),
                elapsed_review_seconds=time.monotonic()-T0,reviewer_max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                review_guard=GUARD,qualification='Full recovered-output ordered content/QC reproduction under qualified declared runtimes and pre-existing header-text compatibility. Raw byte identities use executed gates; no fresh raw/source-reference/runtime-tree body reread. Original scientific limitations remain.')
    with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(dict(verdict=result['verdict'],receipt_sha256=sha(OUT),elapsed=result['elapsed_review_seconds'],RSS=result['reviewer_max_RSS_bytes'])))

if __name__=='__main__':main()
