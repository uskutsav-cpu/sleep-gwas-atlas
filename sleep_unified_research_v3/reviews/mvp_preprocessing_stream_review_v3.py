#!/usr/bin/env python3
"""Bounded independent streaming audit; no worker imports or h2/rg estimates."""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import shutil
import sqlite3
import sys
import time

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
ROOT=PACKAGE.parent
AUDIT=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/mvp_independent_stream_review_v3b')
MAX_RSS=128*1024**2
MAX_FILES=96*1024**2
MAX_SECONDS=1800
INTERNAL_FLOOR=128*1024**2
BEGAN=time.monotonic()


def sha(path):
    before=path.stat()
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(65536),b''):h.update(block)
    after=path.stat()
    return {'path':str(path),'bytes':after.st_size,'sha256':h.hexdigest(),
            'size_mtime_stable':(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)}


def resources():
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak=peak if sys.platform=='darwin' else peak*1024
    size=sum(p.stat().st_size for p in AUDIT.iterdir() if p.is_file())
    if peak>MAX_RSS or size>MAX_FILES or time.monotonic()-BEGAN>MAX_SECONDS or shutil.disk_usage(ROOT).free<INTERNAL_FLOOR:
        raise RuntimeError('INDEPENDENT_AUDIT_RESOURCE_LIMIT')
    return peak,size


def main():
    if AUDIT.exists():raise RuntimeError('EXISTING_AUDIT_PRESERVED')
    if shutil.disk_usage(ROOT).free<256*1024**2 or shutil.disk_usage(AUDIT.parent).free<1024**3:
        raise RuntimeError('INDEPENDENT_AUDIT_LAUNCH_RESOURCE_GATE')
    AUDIT.mkdir()
    plan_path=PACKAGE/'manifests/mvp_preprocessing_execution_plan_v3.json'
    worker_path=PACKAGE/'logs/mvp_preprocessing_worker_receipt_v3.json'
    execution_path=PACKAGE/'logs/mvp_preprocessing_execution_receipt_v3.json'
    plan=json.loads(plan_path.read_text()); worker=json.loads(worker_path.read_text()); execution=json.loads(execution_path.read_text())
    paths={k:Path(v) for k,v in plan['bound_paths'].items()}
    output=Path(plan['output_path'])
    sampled_n=64
    policy={'started_utc':datetime.now(timezone.utc).isoformat(),'review_code_sha256':sha(Path(__file__))['sha256'],
            'maximum_observed_RSS_bytes':MAX_RSS,'maximum_retained_SSD_bytes':MAX_FILES,'maximum_seconds':MAX_SECONDS,
            'internal_launch_minimum_bytes':256*1024**2,'internal_emergency_floor_bytes':INTERNAL_FLOOR,
            'workers':1,'network_bytes':0,'hash_buffer_bytes':65536,'sqlite_cache_KiB':2048,
            'output_rows_all_streamed':True,'deterministic_output_sample_count':sampled_n,
            'sample_selection':'indices 1+floor((826027-1)*i/63), i=0..63; frozen from receipt before source rows accessed',
            'source_samples_and_SNP_level_audit_index_retained_only_on_task_SSD':True,
            'h2_or_rg_estimated':False,'pair_tests_admitted':0}
    with (AUDIT/'independent_audit_resource_plan.json').open('x') as f:json.dump(policy,f,indent=2);f.write('\n')
    hashes_before={k:sha(v) for k,v in paths.items()}
    receipt_before={'execution_plan':sha(plan_path),'worker_receipt':sha(worker_path),'execution_receipt':sha(execution_path),'output':sha(output),'native_log':sha(PACKAGE/'logs/mvp_preprocessing_native_v3.log')}
    checks=[]
    def check(name,actual,expected):
        checks.append({'name':name,'pass':actual==expected,'actual':actual,'expected':expected})
    check('worker_success',worker['status'],'PREPROCESSING_ONLY_PASS')
    check('execution_success',execution['status'],'RESULT_FREE_PREPROCESSING_PASS')
    check('plan_hash_worker',receipt_before['execution_plan']['sha256'],worker['execution_plan_sha256'])
    check('plan_hash_execution',receipt_before['execution_plan']['sha256'],execution['execution_plan_sha256'])
    check('worker_receipt_hash_execution',receipt_before['worker_receipt']['sha256'],execution['worker_receipt_sha256'])
    check('native_log_hash_execution',receipt_before['native_log']['sha256'],execution['log_sha256'])
    check('output_hash_worker',receipt_before['output']['sha256'],worker['output_sha256'])
    check('output_hash_execution',receipt_before['output']['sha256'],execution['output_sha256'])
    check('output_bytes_worker',receipt_before['output']['bytes'],worker['output_bytes'])
    check('all_actual_input_hashes_match_plan',{k:v['sha256'] for k,v in hashes_before.items()},plan['inputs_sha256_before'])
    check('all_actual_input_hashes_match_execution_before_after',execution['inputs_sha256_before']==execution['inputs_sha256_after']=={k:v['sha256'] for k,v in hashes_before.items()},True)
    check('recorded_worker_RSS_within_plan',execution['worker_peak_observed_RSS_bytes']<=plan['worker_observed_RSS_stop_bytes'],True)
    check('recorded_duration_within_plan',execution['elapsed_seconds']<=plan['maximum_seconds'],True)
    check('no_outcome_admission',execution['pair_tests_admitted']==0 and not execution['h2_or_rg_estimated'] and not execution['complete_native_reproduction'],True)
    resources()

    # All output rows live in a disk index, never a Python row dictionary/set.
    db_path=AUDIT/'output_identity_audit.sqlite'
    db=sqlite3.connect(db_path)
    db.execute('PRAGMA cache_size=-2048');db.execute('PRAGMA temp_store=FILE')
    db.execute('PRAGMA journal_mode=OFF');db.execute('PRAGMA synchronous=OFF')
    db.execute('CREATE TABLE snps(rs INTEGER PRIMARY KEY, a1 TEXT, a2 TEXT, ordinal INTEGER, allele_ok INTEGER DEFAULT 0, ref_ok INTEGER DEFAULT 0)')
    expected_rows=worker['counts']['output_rows']
    sample_positions={1+(expected_rows-1)*i//(sampled_n-1) for i in range(sampled_n)}
    samples={};output_counts=Counter();expected_n=format(float(Fraction(4*78566*329572,408138)),'.12g')
    with gzip.open(output,'rt',encoding='utf-8',newline='') as f:
        reader=csv.reader(f,delimiter='\t')
        check('output_header',next(reader),['SNP','A1','A2','Z','N'])
        for ordinal,row in enumerate(reader,1):
            output_counts['rows']+=1
            if len(row)!=5:raise RuntimeError('OUTPUT_COLUMN_COUNT')
            rs,a1,a2,z,n=row
            if not re.fullmatch('rs[1-9][0-9]*',rs):raise RuntimeError('OUTPUT_RSID_SYNTAX')
            if a1 not in 'ACGT' or a2 not in 'ACGT' or len(a1)!=1 or len(a2)!=1 or a1==a2 or {a1,a2} in [{'A','T'},{'C','G'}]:raise RuntimeError('OUTPUT_ALLELE_PAIR')
            if not math.isfinite(float(z)) or not math.isfinite(float(n)) or n!=expected_n:raise RuntimeError('OUTPUT_Z_N')
            if format(float(z),'.12g')!=z:raise RuntimeError('OUTPUT_Z_SERIALIZATION')
            try:db.execute('INSERT INTO snps(rs,a1,a2,ordinal) VALUES(?,?,?,?)',(int(rs[2:]),a1,a2,ordinal))
            except sqlite3.IntegrityError:raise RuntimeError('DUPLICATE_OUTPUT_RSID')
            if ordinal in sample_positions:samples[rs]={'output_ordinal':ordinal,'output':row}
            if ordinal%50000==0:db.commit();resources()
    db.commit()
    check('all_output_rows_count',output_counts['rows'],826027)
    check('output_count_matches_worker',output_counts['rows'],expected_rows)
    check('all_output_rows_unique',db.execute('SELECT COUNT(*) FROM snps').fetchone()[0],output_counts['rows'])
    check('deterministic_sample_count',len(samples),sampled_n)
    check('source_exclusion_reconciliation_from_receipt',worker['counts']['source_rows'],sum(worker['exclusive_exclusions'].values())+output_counts['rows'])
    sample_ids=set(samples)

    # Each output allele/coordinate is independently checked against pinned
    # references through the small disk index. Only the 64 sample reference
    # records are kept in a Python dictionary.
    reference_subset={};allele_counts=Counter()
    with paths['hm3_alleles'].open(newline='',encoding='utf-8') as f:
        for i,row in enumerate(csv.DictReader(f,delimiter='\t'),1):
            rs=row['SNP'];number=int(rs[2:])
            target=db.execute('SELECT a1,a2,allele_ok FROM snps WHERE rs=?',(number,)).fetchone()
            if target is not None:
                if target[2]:raise RuntimeError('DUPLICATE_OUTPUT_REFERENCE_ALLELE')
                if target[:2]!=(row['A1'].upper(),row['A2'].upper()):raise RuntimeError('OUTPUT_TARGET_ALLELE_MISMATCH')
                db.execute('UPDATE snps SET allele_ok=1 WHERE rs=?',(number,));allele_counts['matched']+=1
            if i%100000==0:db.commit();resources()
    db.commit()
    check('every_output_allele_matches_pinned_A1_A2',allele_counts['matched'],output_counts['rows'])
    reference_counts=Counter()
    with gzip.open(paths['hm3_reference'],'rt',encoding='utf-8',newline='') as f:
        for i,row in enumerate(csv.DictReader(f,delimiter='\t'),1):
            rs=row['SNP'];number=int(rs[2:])
            target=db.execute('SELECT ref_ok FROM snps WHERE rs=?',(number,)).fetchone()
            if target is not None:
                if target[0]:raise RuntimeError('DUPLICATE_OUTPUT_REFERENCE_COORDINATE')
                chromosome,position=int(row['CHR']),int(row['BP'])
                if not 1<=chromosome<=22 or position<1 or (chromosome==6 and 25000000<=position<=34000000):raise RuntimeError('OUTPUT_REFERENCE_GEOMETRY_MHC')
                db.execute('UPDATE snps SET ref_ok=1 WHERE rs=?',(number,));reference_counts['matched']+=1
                if rs in sample_ids:reference_subset[rs]=(chromosome,position)
            if i%100000==0:db.commit();resources()
    db.commit()
    check('every_output_has_pinned_GRCh37_reference_and_outside_MHC',reference_counts['matched'],output_counts['rows'])
    check('all_output_reference_flags_complete',db.execute('SELECT COUNT(*) FROM snps WHERE allele_ok!=1 OR ref_ok!=1').fetchone()[0],0)
    db.close()
    print(json.dumps({'stage':'OUTPUT_AND_REFERENCES_VERIFIED','rows':output_counts['rows'],'sample_count':len(samples)}),flush=True)

    # Stream all raw rows to CRC completion. Keep source excerpts only for the
    # prospectively selected 64 identities, with a strict candidate-count cap.
    source_counts=Counter();r2_min=math.inf;r2_max=-math.inf;source_candidates=defaultdict(list);candidate_count=0
    expected_header=json.loads(paths['original_schema_receipt'].read_text())['source_header']
    with gzip.open(paths['source'],'rt',encoding='utf-8',newline='') as f:
        reader=csv.reader(f,delimiter='\t');header=next(reader)
        check('source_header_exact',header,expected_header)
        fields={name:i for i,name in enumerate(header)}
        for ordinal,row in enumerate(reader,1):
            source_counts['rows']+=1
            if len(row)!=len(header):raise RuntimeError('SOURCE_ROW_COLUMN_COUNT')
            if (float(row[fields['n']]),float(row[fields['num_cases']]),float(row[fields['num_controls']]))!=(408138,78566,329572):raise RuntimeError('SOURCE_SAMPLE_COUNT')
            if row[fields['standard_error']].strip().lower() not in {'','na','#na','nan','null','none','.'}:raise RuntimeError('SOURCE_NONMISSING_SE')
            source_counts['literal_missing_SE']+=1
            r2=float(row[fields['r2']]);r2_min=min(r2_min,r2);r2_max=max(r2_max,r2);source_counts['R2_gt_1']+=r2>1
            raw_id=row[fields['rsid']]
            if raw_id in sample_ids:
                matches=[raw_id]
            elif raw_id.startswith('rs') and raw_id[2:].isdigit():
                matches=[]
            else:
                matches=[rs for rs in ('rs'+n for n in re.findall('rs([0-9]+)',raw_id)) if rs in sample_ids]
            for rs in matches:
                source_candidates[rs].append({'source_ordinal':ordinal,'source':dict(zip(header,row))})
                candidate_count+=1
            if candidate_count>512:raise RuntimeError('SOURCE_SAMPLE_CANDIDATE_CAP')
            if ordinal%1000000==0:
                resources();print(json.dumps({'stage':'SOURCE_STREAM','rows':ordinal}),flush=True)
    check('source_full_stream_rows',source_counts['rows'],19703815)
    check('source_literal_missing_SE_all_rows',source_counts['literal_missing_SE'],19703815)
    check('source_R2_range', [r2_min,r2_max],worker['R2_range'])
    check('source_R2_above_one_count',source_counts['R2_gt_1'],worker['whole_file_diagnostics']['R2_gt_1'])
    check('all_selected_ids_have_source_candidates',len(source_candidates),sampled_n)
    with (AUDIT/'source_candidate_checkpoint.json').open('x') as f:
        json.dump({'samples':samples,'source_candidates':source_candidates,'source_counts':dict(source_counts),'source_CRC_stream_completed':True},f)
        f.write('\n')

    # A separate chain parser stores only hits for the selected source points.
    # It imports neither the worker nor the implementation's liftover helper.
    point_by_chr=defaultdict(list)
    for rs,candidates in source_candidates.items():
        for c in candidates:
            chromosome=int(c['source']['chromosome']);position=int(c['source']['base_pair_location'])
            point_by_chr[chromosome].append((position,rs,c['source_ordinal']))
    for values in point_by_chr.values():values.sort()
    positions={k:[v[0] for v in values] for k,values in point_by_chr.items()}
    mappings=defaultdict(set);chain_counts=Counter();current=None;source_cursor=query_cursor=0
    with gzip.open(paths['chain'],'rt',encoding='ascii') as f:
        for line in f:
            parts=line.split()
            if not parts:current=None;continue
            if parts[0]=='chain':
                if len(parts)!=13 or parts[4]!='+' or parts[9] not in ['+','-']:raise RuntimeError('CHAIN_HEADER')
                source_chr=int(parts[2][3:]) if re.fullmatch('chr[0-9]+',parts[2]) else None
                query_chr=int(parts[7][3:]) if re.fullmatch('chr[0-9]+',parts[7]) else None
                current=(source_chr,query_chr,int(parts[8]),parts[9]);source_cursor=int(parts[5]);query_cursor=int(parts[10]);chain_counts['chains']+=1
                continue
            if current is None or len(parts) not in [1,3]:raise RuntimeError('CHAIN_BLOCK')
            size=int(parts[0]);dt=int(parts[1]) if len(parts)==3 else 0;dq=int(parts[2]) if len(parts)==3 else 0
            sc,qc,qsize,strand=current
            if sc in point_by_chr:
                low=bisect_left(positions[sc],source_cursor+1);high=bisect_right(positions[sc],source_cursor+size)
                for position,rs,ordinal in point_by_chr[sc][low:high]:
                    offset=position-source_cursor-1
                    bp=query_cursor+offset+1 if strand=='+' else qsize-query_cursor-offset
                    mappings[(rs,ordinal)].add((qc,bp,strand))
            source_cursor+=size+dt;query_cursor+=size+dq;chain_counts['blocks']+=1
    check('independent_chain_count',chain_counts['chains'],worker['chain_receipt']['chain_count'])
    check('independent_chain_block_count',chain_counts['blocks'],worker['chain_receipt']['block_count'])

    complement={'A':'T','T':'A','C':'G','G':'C'};sample_counts=Counter();sample_records=[]
    max_text_relative_difference=0.;max_text_last_digit_units=0.;max_binary_last_digit_units=0.
    for rs,sample in samples.items():
        eligible=[]
        for candidate in source_candidates[rs]:
            row=candidate['source'];hits=mappings[(rs,candidate['source_ordinal'])]
            if len(hits)!=1:continue
            chromosome,position,strand=next(iter(hits))
            if (chromosome,position)!=reference_subset[rs]:continue
            if chromosome==6 and 25000000<=position<=34000000:continue
            odds,low,high,p,af,r2=(float(row[k]) for k in ['odds_ratio','ci_lower','ci_upper','p_value','effect_allele_frequency','r2'])
            if not all(map(math.isfinite,[odds,low,high,p,af,r2])) or not(odds>0 and 0<low<high and 0<=p<=1 and 0<=af<=1 and min(af,1-af)>.01 and r2>.9):continue
            effect,other=row['effect_allele'].upper(),row['other_allele'].upper()
            if effect not in complement or other not in complement or effect==other or {effect,other} in [{'A','T'},{'C','G'}]:continue
            if strand=='-':effect,other=complement[effect],complement[other]
            a1,a2=sample['output'][1:3]
            pairs=[(effect,other),(complement[effect],complement[other])]
            sign=1 if (a1,a2) in pairs else -1 if (a2,a1) in pairs else None
            if sign is None:continue
            with localcontext() as ctx:
                ctx.prec=50
                log_or=Decimal(row['odds_ratio']).ln()
                uncertainty=(Decimal(row['ci_upper']).ln()-Decimal(row['ci_lower']).ln())/(Decimal(2)*Decimal('1.959963984540054'))
                exact_z=Decimal(sign)*log_or/uncertainty
                # Match the frozen input representation, independently taking
                # high-precision logs of the exact input binary doubles.
                binary_log_or=Decimal.from_float(odds).ln()
                binary_uncertainty=(Decimal.from_float(high).ln()-Decimal.from_float(low).ln())/Decimal.from_float(2*1.959963984540054)
                binary_z=Decimal(sign)*binary_log_or/binary_uncertainty
            actual=float(sample['output'][3]);estimate=float(binary_z);text_estimate=float(exact_z)
            quantum=10**(math.floor(math.log10(abs(estimate)))-11) if estimate else 1e-300
            within_precision=abs(actual-estimate)<=.5001*quantum+16*math.ulp(estimate)
            text_within_precision=abs(actual-text_estimate)<=.5001*quantum+16*math.ulp(text_estimate)
            eligible.append((candidate,sign,strand,exact_z,binary_z,within_precision,text_within_precision,quantum))
        if not eligible:raise RuntimeError('SAMPLE_NO_INDEPENDENT_ELIGIBLE_SOURCE')
        eligible.sort(key=lambda value:value[0]['source_ordinal'])
        chosen,sign,strand,z_decimal,binary_decimal,precision_ok,text_precision_ok,quantum=eligible[0]
        if not precision_ok:raise RuntimeError('SAMPLE_Z_OUTSIDE_DECLARED_PRECISION')
        if sample['output'][4]!=expected_n:raise RuntimeError('SAMPLE_NEFF')
        sample_counts['passed']+=1;sample_counts['orientation_sign_'+str(sign)]+=1;sample_counts['chain_'+strand]+=1
        sample_counts['eligible_duplicate_candidates']+=len(eligible)-1
        sample_counts['exact_text_decimal_formatted_Z_match']+=format(float(z_decimal),'.12g')==sample['output'][3]
        sample_counts['exact_binary_input_formatted_Z_match']+=format(float(binary_decimal),'.12g')==sample['output'][3]
        sample_counts['exact_text_decimal_within_12g_bound']+=text_precision_ok
        text_difference=abs(float(sample['output'][3])-float(z_decimal))
        binary_difference=abs(float(sample['output'][3])-float(binary_decimal))
        max_text_relative_difference=max(max_text_relative_difference,text_difference/max(abs(float(z_decimal)),1e-300))
        max_text_last_digit_units=max(max_text_last_digit_units,text_difference/quantum)
        max_binary_last_digit_units=max(max_binary_last_digit_units,binary_difference/quantum)
        sample_records.append({**sample,**chosen,'independent_GRCh37_reference':reference_subset[rs],
                               'independent_mapping':list(mappings[(rs,chosen['source_ordinal'])]),
                               'orientation_sign':sign,'independent_exact_text_decimal_Z':str(z_decimal),
                               'independent_exact_binary_input_decimal_Z':str(binary_decimal),
                               'matches_binary_input_12g_precision':precision_ok,'matches_exact_text_12g_precision':text_precision_ok,
                               'eligible_source_candidate_count':len(eligible)})
    check('deterministic_source_samples_pass_geometry_alleles_CI_Z_N',sample_counts['passed'],sampled_n)
    with (AUDIT/'deterministic_source_samples.jsonl').open('x') as f:
        for record in sample_records:f.write(json.dumps(record)+'\n')
    peak,size=resources()
    hashes_after={k:sha(v) for k,v in paths.items()}
    check('independent_bound_inputs_unchanged_during_audit',{k:v['sha256'] for k,v in hashes_before.items()},{k:v['sha256'] for k,v in hashes_after.items()})
    check('all_hashed_input_size_mtime_stable',all(v['size_mtime_stable'] for v in list(hashes_before.values())+list(hashes_after.values())),True)
    check('output_unchanged_after_audit',sha(output)['sha256'],receipt_before['output']['sha256'])
    audit_files={p.name:sha(p) for p in AUDIT.iterdir() if p.is_file()}
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'status':'INDEPENDENT_STREAMING_PREPROCESSING_PASS' if all(c['pass'] for c in checks) else 'INDEPENDENT_AUDIT_FAILURE',
            'review_code':sha(Path(__file__)),'input_hashes_before':hashes_before,'input_hashes_after':hashes_after,
            'receipt_and_output_hashes':receipt_before,'checks':checks,'check_count':len(checks),
            'all_output_rows_streamed_and_CRC_verified':True,'raw_source_full_stream_CRC_verified':True,
            'source_aggregate_counts':dict(source_counts),'output_aggregate_counts':dict(output_counts),
            'deterministic_sample_aggregate_counts':dict(sample_counts),'independent_chain_aggregate_counts':dict(chain_counts),
            'maximum_sample_Z_relative_difference_from_exact_source_text':max_text_relative_difference,
            'maximum_sample_Z_difference_from_exact_source_text_in_last_printed_digit_units':max_text_last_digit_units,
            'maximum_sample_Z_difference_from_exact_binary_inputs_in_last_printed_digit_units':max_binary_last_digit_units,
            'SSD_only_sample_and_audit_files':audit_files,'resource_plan':policy,'peak_observed_RSS_bytes':peak,
            'retained_SSD_bytes':size,'elapsed_seconds':time.monotonic()-BEGAN,'source_or_frozen_input_mutations':False,
            'h2_or_rg_estimated':False,'pair_tests_admitted':0,'full_native_reproduction':False,
            'limits':['All output structural/reference checks passed; source-to-output statistical and chain arithmetic independently checked for 64 predetermined rows using exact binary-float input values with independent 50-digit decimal logs, not all retained rows.',
                      'Printed 12g Z represents the frozen binary-float implementation. Exact original decimal source-text calculations can differ in final digits near OR=1; conditioning difference is quantified separately without calling it scientific uncertainty calibration.',
                      'CI-derived uncertainty is conditional on the frozen 95% log-scale Wald assumption and does not certify upstream effect semantics.',
                      'Source exclusions and every retained-row mapping/sign were not independently regenerated; no h2, rg, cohort overlap, 396/1200 completion or replication claim.',
                      'Hashes before/after and stat stability support unchanged inputs at observations, not a proof against transient external edits.']}
    with (HERE/'mvp_preprocessing_stream_review_v3.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'source_rows':source_counts['rows'],'output_rows':output_counts['rows'],'samples':dict(sample_counts),'peak_RSS':peak,'elapsed':result['elapsed_seconds']}),flush=True)
    if not all(c['pass'] for c in checks):raise SystemExit(1)


if __name__=='__main__':main()
