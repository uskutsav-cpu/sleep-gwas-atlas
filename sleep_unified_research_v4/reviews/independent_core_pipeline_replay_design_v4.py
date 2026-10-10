#!/usr/bin/env python3
"""Build proposed 45-trait replay provenance from code/metadata only."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
REVIEW=P/'reviews'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
NEW=SSD/'core_pipeline/independent_replay_proposal_v1'
ORIGINAL='659d01cf5b10a6debf1841404e7e01ef4750f755'
LDSC=ROOT.parent/'ldsc-code'
PYTHON=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env/bin/python')
META={}

def sha(path):
    path=Path(path);h=hashlib.sha256()
    assert path.stat().st_size<32<<20
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    META[str(path)]=dict(bytes=path.stat().st_size,sha256=h.hexdigest())
    return h.hexdigest()

def read_tsv(path):
    sha(path)
    with Path(path).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def record_path(row,kind='audit'):
    path=Path(row['path']);s=path.stat() if path.is_file() else None
    return dict(path=str(path),bytes=int(row['bytes'] or 0),sealed_sha256=row['actual_sha256'],
        historical_expected_sha256=row['expected_sha256'] or None,
        hash_evidence_status=row['status'],expected_hash_evidence=row['expected_hash_evidence'] or None,
        current_exists=bool(s),current_size_matches_record=bool(s and s.st_size==int(row['bytes'])),
        fresh_body_hash_verified_by_this_review=False)

def qc(path):
    sha(path);meta={};steps=[]
    for line in path.read_text().splitlines():
        v=line.split('\t')
        if len(v)==2:meta[v[0]]=v[1]
        elif len(v)==3 and v[1].isdigit() and v[2].isdigit():
            steps.append(dict(reason=v[0],dropped=int(v[1]),remaining=int(v[2])))
    ni,no=int(meta['rows_in']),int(meta['rows_out'])
    assert sum(x['dropped'] for x in steps)==ni-no
    return meta,steps

def main():
    started=time.monotonic()
    panel=read_tsv(ROOT/'config/analysis_panel.tsv')
    assert len(panel)==45 and len({r['trait_id'] for r in panel})==45
    inp=read_tsv(ROOT/'sleep_unified_research_v1/tables/native_input_hash_checks.tsv')
    read_tsv(ROOT/'sleep_unified_research_v1/tables/core_candidate_hash_audit.tsv')
    current=read_tsv(ROOT/'sleep_unified_research_v3/tables/CURRENT_CORE_SOURCE_STATUS_V3.tsv')
    read_tsv(ROOT/'sleep_unified_research_v2/tables/SOURCE_RECOVERY_UPDATES.tsv')
    pre={r['trait_id']:r for r in read_tsv(ROOT/'config/hm3_prefilter_plans.tsv')}
    maps={r['trait_id']:r for r in read_tsv(ROOT/'config/variant_mapping_plans.tsv')}
    lifts={r['trait_id']:r for r in read_tsv(ROOT/'config/liftover_plans.tsv')}
    schemas={r['trait_id']:r for r in read_tsv(ROOT/'config/gwas_schemas.tsv')}
    assert len(schemas)==45
    read_tsv(ROOT/'config/public_gwas_sources.tsv')
    by={(r['kind'],r['trait_id']):r for r in inp if r['trait_id']}
    source={(r['kind'],r['trait_id']):r for r in current}
    refs={r['kind']:r for r in inp if r['kind'] in ('variant_map','liftover_chain')}
    allow=next(r for r in inp if r['kind']=='official_ld_reference_member' and r['path'].endswith('w_hm3.snplist'))
    for r in list(refs.values())+[allow]:assert sha(r['path'])==r['actual_sha256']==r['expected_sha256']
    map_prov=Path(refs['variant_map']['path']+'.provenance.json');sha(map_prov)
    native_sha=sha(PYTHON)
    assert native_sha=='7dffb088cd3027e48f0127ced6f206e06111abed3b9457e580b6ebbf593c10ba'
    code={}
    for name in ['scripts/01_harmonize.py','scripts/21_prefilter_hm3.py','scripts/variant_map.py','scripts/liftover_chain.py',
                 'scripts/02_munge.sh','scripts/_common.sh','scripts/00_validate_panel.py']:
        h=sha(ROOT/name)
        original=subprocess.check_output(['git','show',ORIGINAL+':'+name],cwd=ROOT)
        blob=subprocess.check_output(['git','rev-parse',ORIGINAL+':'+name],cwd=ROOT,text=True).strip()
        code[name]=dict(path=str(ROOT/name),sha256=h,original_commit=ORIGINAL,original_blob=blob,
            original_sha256=hashlib.sha256(original).hexdigest(),exact_original_source_bytes=(ROOT/name).read_bytes()==original)
    for name in ['scripts/01_harmonize.py','scripts/21_prefilter_hm3.py','scripts/variant_map.py','scripts/liftover_chain.py']:
        assert code[name]['exact_original_source_bytes']
    assert sha(LDSC/'munge_sumstats.py')==hashlib.sha256(subprocess.check_output(['git','show','6c673952cee74bd5c57aef1555a03b1c015399a0:munge_sumstats.py'],cwd=LDSC)).hexdigest()
    for name in ['environment/workflow.yml','environment/tool_versions.tsv','sleep_unified_research_v1/sources/recovered/environment/tool_versions.tsv']:
        sha(ROOT/name)
    rows=[]
    for meta in panel:
        t=meta['trait_id'];raw=source['core_raw',t];oldraw=by['core_raw',t]
        assert raw['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE'
        assert raw['current_sha256']==raw['expected_sha256']==oldraw['expected_sha256']
        rp=Path(raw['current_path']);assert rp.is_file() and rp.stat().st_size==int(raw['current_bytes'])
        harm=record_path(by['core_harmonized',t]);munged=record_path(by['core_munged',t])
        assert harm['current_size_matches_record'] and munged['current_size_matches_record']
        qpath=Path(harm['path']).with_name(t+'.qc.txt');q,steps=qc(qpath)
        mlog=Path(munged['path']).with_name(t+'.log');sha(mlog);ml=mlog.read_text()
        read_match=re.search(r'Read ([0-9]+) SNPs from --sumstats file',ml)
        write_match=re.search(r'Writing summary statistics for ([0-9]+) SNPs \(([0-9]+) with nonmissing beta\)',ml)
        assert read_match and int(read_match.group(1))==int(q['rows_out'])
        assert write_match and int(write_match.group(1))==1217311
        pref=None;pref_input=str(rp);pref_prov=None;pref_command=None;pref_status='NOT_REQUIRED'
        if t in pre:
            pref=record_path(by['core_prefiltered',t]);assert pref['current_size_matches_record']
            assert pref['historical_expected_sha256']==pref['sealed_sha256']==q['infile_sha256']
            assert q['prefilter_source_sha256']==raw['current_sha256']
            pref_input=str(NEW/t/'prefilter'/f'{t}.hm3.tsv.gz')
            pref_prov=str(NEW/t/'prefilter'/f'{t}.hm3.provenance.json')
            pref_status='REPLAY_REQUIRED'
            pref_command=[str(PYTHON),'-B',str(ROOT/'scripts/21_prefilter_hm3.py'),'--input',str(rp),
                '--output',pref_input,'--provenance',pref_prov,'--snp-column',pre[t]['snp_column'],
                '--allowlist',allow['path'],'--expected-input-bytes',raw['current_bytes'],
                '--expected-input-sha256',raw['current_sha256']]
            if t in ['ldl','hdl','triglycerides']:
                receipt=P.parent/'sleep_unified_research_v2/logs'/f'{t}_native_prefilter_receipt_v2.json'
                sha(receipt);rr=json.loads(receipt.read_text())
                assert rr['status']=='NATIVE_PREFILTER_EXACT_HISTORICAL_BYTES_AND_COUNTS'
                assert rr['input_hashes_unchanged'] and rr['inputs_sha256_before']==rr['inputs_sha256_after']
                assert rr['output_sha256']==pref['sealed_sha256']
                assert rr['comparison']['source_rows']==[int(q['prefilter_source_rows'])]*2
                assert rr['comparison']['retained_rows']==[int(q['prefilter_retained_rows'])]*2
                assert rr['inputs_sha256_before']['native_script']==code['scripts/21_prefilter_hm3.py']['sha256']
                pref_input=rr['output_path'];pref_prov=rr['native_prefilter_provenance_path']
                prov=json.loads(Path(pref_prov).read_text());assert sha(pref_prov)==rr['native_prefilter_provenance_sha256']
                assert prov['output_sha256']==pref['sealed_sha256'] and prov['input_sha256']==raw['current_sha256']
                pref_status='REUSE_ALREADY_EXACT_NATIVE_REPLAY_NO_NEW_BODY_READ';pref_command=None
                pref.update(replay_receipt=str(receipt),replay_receipt_sha256=META[str(receipt)]['sha256'],
                    reuse_path=pref_input,reuse_provenance=pref_prov)
        else:assert q['infile_sha256']==raw['current_sha256']
        hcmd=['{FROZEN_HARMONIZATION_PYTHON}','-B',str(ROOT/'scripts/01_harmonize.py'),
            '--trait',t,'--config',str(ROOT/'config/analysis_panel.tsv'),'--infile',pref_input,
            '--outdir',str(NEW/t/'harmonized')]
        if meta['build'] in ('hg19','GRCh37','hg38','GRCh38'):hcmd+=['--source-build',meta['build']]
        if t in maps:
            v=maps[t];hcmd+=['--variant-map',refs['variant_map']['path'],'--variant-map-strategy',v['strategy'],
                '--expected-variant-map-bytes',v['map_bytes'],'--expected-variant-map-sha256',v['map_sha256']]
        if t in lifts:
            v=lifts[t];hcmd+=['--liftover-chain',refs['liftover_chain']['path'],
                '--expected-liftover-chain-bytes',v['chain_bytes'],'--expected-liftover-chain-sha256',v['chain_sha256']]
        if pref_prov:hcmd+=['--prefilter-provenance',pref_prov]
        ignore=any('FRQ column absent - source-level MAF QC must be documented' in x['reason'] for x in steps)
        mcmd=[str(PYTHON),'-B',str(LDSC/'munge_sumstats.py'),'--sumstats',str(NEW/t/'harmonized'/f'{t}.harmonized.tsv.gz'),
            '--merge-alleles',allow['path'],'--chunksize','500000']
        if ignore:mcmd+=['--ignore','FRQ']
        mcmd+=['--out',str(NEW/t/'munged'/t)]
        assert ('--ignore FRQ' in ml)==ignore
        nmode=next(x['reason'].split(': ',1)[1] for x in steps if x['reason'].startswith('sample-size mode:'))
        risk='LIKELY_EXCEEDS_2_GIB_FULL_DATAFRAME_REQUIRES_STREAMING_OR_JUSTIFIED_AMENDMENT' if int(q['rows_in'])>=5000000 else 'FULL_DATAFRAME_PEAK_UNMEASURED_REQUIRES_SUPERVISED_RESOURCE_PLAN'
        blockers=['EXECUTION_NOT_ADMITTED','HARMONIZATION_RUNTIME_BINDING_OR_EXPLICIT_VERSION_QUALIFICATION_REQUIRED',
                  'RESOURCE_SUPERVISOR_AND_IMMUTABLE_FAILURE_ORACLE_REQUIRED','BASELINE_DECOMPRESSED_CONTENT_HASHES_NOT_YET_ESTABLISHED_BY_THIS_REVIEW']
        if int(q['rows_in'])>=5000000:blockers.append('ORIGINAL_01_FULL_DATAFRAME_MEMORY_RISK')
        archive=source.get(('core_source_archive',t))
        rows.append(dict(trait_id=t,source_id=meta['source_id'],source_release=meta['dataset_version'],build=meta['build'],ancestry=meta['ancestry'],
            phenotype_type=meta['type'],phenotype_definition=meta['phenotype_definition'],source_schema=schemas[t],original_panel_metadata=meta,
            raw=dict(original_path=oldraw['path'],resolved_path=str(rp),bytes=int(raw['current_bytes']),
                historical_expected_sha256=raw['expected_sha256'],sealed_verified_sha256=raw['current_sha256'],proof=raw['proof'],
                status=raw['current_status'],fresh_body_hash_verified_by_this_review=False),source_archive=archive,
            prefilter=pref,prefilter_plan=pre.get(t),prefilter_replay_status=pref_status,
            raw_source_rows=int(q.get('prefilter_source_rows',q['rows_in'])),harmonizer_input_rows=int(q['rows_in']),
            harmonized_rows=int(q['rows_out']),harmonized=harm,munged=munged,
            historical_QC=dict(path=str(qpath),sha256=META[str(qpath)]['sha256'],metadata=q,ordered_steps=steps),
            historical_munge_log=dict(path=str(mlog),sha256=META[str(mlog)]['sha256'],output_total_rows=int(write_match.group(1)),
                output_nonmissing_rows=int(write_match.group(2)),warnings=[x for x in ml.splitlines() if 'WARNING' in x]),
            sample_size_mode=nmode,ignore_FRQ_only_with_original_absence_proof=ignore,
            mapping_plan=maps.get(t),liftover_plan=lifts.get(t),memory_risk=risk,method_blockers=blockers,
            new_private_output_namespace=str(NEW/t),new_harmonized_path=str(NEW/t/'harmonized'/f'{t}.harmonized.tsv.gz'),
            new_munged_path=str(NEW/t/'munged'/f'{t}.sumstats.gz'),
            command_templates=dict(prefilter=pref_command,harmonize_original=hcmd,munge_stock=mcmd),
            worker_environment=dict(TMPDIR=str(NEW/t/'tmp'),XDG_CACHE_HOME=str(NEW/t/'cache'),PYTHONDONTWRITEBYTECODE='1',
                OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1'),
            comparison_rule='EXACT_ORDERED_DECOMPRESSED_CONTENT_AND_QC_COUNTS_REQUIRED; COMPRESSED_SHA_REPORTED_SEPARATELY; NO_NEW_TOLERANCE_OR_QC_CHANGES'))
    assert len(rows)==45
    pending=[r['trait_id'] for r in rows if r['prefilter_replay_status']=='REPLAY_REQUIRED']
    assert set(pending)=={'ms','asthma','t2d','cad','telomere_length','melanoma'}
    before=dict(META)
    for path,r in before.items():assert sha(path)==r['sha256']
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<256<<20
    manifest=dict(schema='independent_proposed_core_pipeline_replay_design_v1',prepared_utc=datetime.now(timezone.utc).isoformat(),
        status='DESIGN_AND_METADATA_ONLY_EXECUTION_NOT_ADMITTED',original_main_commit=ORIGINAL,
        proposed_private_namespace=str(NEW),trait_count=45,raw_exact_sealed_identities=45,
        original_prefilters=9,already_native_reproduced_prefilters=['ldl','hdl','triglycerides'],remaining_prefilters=pending,
        source_archives=dict(total=len([r for r in current if r['kind']=='core_source_archive']),
            exact=len([r for r in current if r['kind']=='core_source_archive' and r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE']),
            unresolved=[r['trait_id'] for r in current if r['kind']=='core_source_archive' and r['current_status']!='EXACT_HISTORICAL_HASH_AVAILABLE']),
        code=code,ldsc_code_commit='6c673952cee74bd5c57aef1555a03b1c015399a0',native_python=dict(path=str(PYTHON),sha256=native_sha,
            declared_version='3.9.23',declared_numpy='1.21.5',declared_pandas='1.3.3',declared_scipy='1.7.3'),
        declared_original_workflow_runtime=dict(python='3.11.11',numpy='1.26.4',pandas='2.2.3',scipy='1.14.1',
            actual_per_trait_harmonization_runtime_not_proved_by_QC_ledgers=True,matching_binary_not_established_in_this_review=True),
        input_metadata_code_reference_hashes=before,all_metadata_hashes_unchanged=True,GWAS_body_reads=0,
        fresh_GWAS_hashes_or_decompression=0,workers_launched=0,rows=rows)
    target=REVIEW/'independent_core_pipeline_replay_manifest_v4.json'
    with target.open('x') as f:json.dump(manifest,f,indent=2);f.write('\n')
    fields=['trait_id','source_id','source_release','build','raw_path','raw_sha256','raw_bytes','prefilter_status','prefilter_baseline_sha256',
            'harmonized_baseline_sha256','munged_baseline_sha256','raw_source_rows','harmonizer_input_rows','harmonized_rows','munged_nonmissing_rows',
            'sample_size_mode','mapping','liftover','ignore_FRQ','memory_risk','new_output_namespace']
    with (REVIEW/'independent_core_pipeline_replay_manifest_v4.tsv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for r in rows:w.writerow(dict(trait_id=r['trait_id'],source_id=r['source_id'],source_release=r['source_release'],build=r['build'],
            raw_path=r['raw']['resolved_path'],raw_sha256=r['raw']['sealed_verified_sha256'],raw_bytes=r['raw']['bytes'],
            prefilter_status=r['prefilter_replay_status'],prefilter_baseline_sha256=r['prefilter']['sealed_sha256'] if r['prefilter'] else '',
            harmonized_baseline_sha256=r['harmonized']['sealed_sha256'],munged_baseline_sha256=r['munged']['sealed_sha256'],
            raw_source_rows=r['raw_source_rows'],harmonizer_input_rows=r['harmonizer_input_rows'],harmonized_rows=r['harmonized_rows'],
            munged_nonmissing_rows=r['historical_munge_log']['output_nonmissing_rows'],sample_size_mode=r['sample_size_mode'],
            mapping=r['mapping_plan']['strategy'] if r['mapping_plan'] else '',liftover=bool(r['liftover_plan']),
            ignore_FRQ=r['ignore_FRQ_only_with_original_absence_proof'],memory_risk=r['memory_risk'],new_output_namespace=r['new_private_output_namespace']))
    summary=dict(status=manifest['status'],completed_utc=manifest['prepared_utc'],checker_sha256=sha(__file__),
        manifest_sha256=sha(target),traits=45,raw_resolved=45,remaining_prefilters=pending,
        high_memory_risk_traits=[r['trait_id'] for r in rows if r['harmonizer_input_rows']>=5000000],
        historical_harmonized_or_munged_checkpoint_hashes_present=0,metadata_dependency_count=len(before),
        metadata_code_reference_bytes=sum(r['bytes'] for r in before.values()),before_after_metadata_hashes_pass=True,
        GWAS_body_reads=0,workers_launched=0,elapsed_seconds=time.monotonic()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with (REVIEW/'independent_core_pipeline_replay_design_receipt_v4.json').open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
    print(json.dumps(summary))

if __name__=='__main__':main()
