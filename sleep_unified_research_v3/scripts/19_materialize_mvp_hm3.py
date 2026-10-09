#!/usr/bin/env python3
"""Result-free, separately versioned MVP materialization and diagnostics."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
sys.path.insert(0, str(ROOT/'discovery_extension/scripts'))
spec = importlib.util.spec_from_file_location('historical47', ROOT/'discovery_extension/scripts/47_stream_replication_sources.py')
historical = importlib.util.module_from_spec(spec)
spec.loader.exec_module(historical)
sys.path.insert(0, str(ROOT/'scripts'))
from liftover_chain import load_chain

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536), b''): h.update(b)
    return h.hexdigest()

def finite(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError): return None

def bucket(value):
    return next((str(t) for t in [.001,.01,.05,.1,1,10] if value<=t), 'gt10')

def main(plan_path):
    plan=json.loads(Path(plan_path).read_text())
    paths={k:Path(v) for k,v in plan['bound_paths'].items()}
    hm3=historical.load_hm3(paths['hm3_reference'], paths['hm3_alleles'])
    chain, chain_receipt=load_chain(paths['chain'], plan['inputs_sha256_before']['chain'], 1246411)
    print('REFERENCE_READY nonambiguous_HM3='+str(len(hm3)), flush=True)
    counts=Counter(); excluded=Counter(); diagnostic=Counter(); p_buckets=Counter(); center_buckets=Counter()
    r2_range=[math.inf,-math.inf]; seen=set(); began=time.time()
    out=Path(plan['output_path']); tmp=Path(str(out)+'.partial')
    required={'chromosome','base_pair_location','effect_allele','other_allele','odds_ratio','standard_error','effect_allele_frequency','p_value','rsid','ci_upper','ci_lower','n','num_cases','num_controls','r2'}
    with gzip.open(paths['source'],'rt',newline='') as f, tmp.open('xb') as rawout:
        with gzip.GzipFile(filename='', mode='wb', fileobj=rawout, compresslevel=1, mtime=0) as compressed:
            import io
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as dest:
                writer=csv.writer(dest,delimiter='\t',lineterminator='\n'); writer.writerow(['SNP','A1','A2','Z','N'])
                reader=csv.DictReader(f,delimiter='\t'); header=reader.fieldnames
                expected_header=json.loads(paths['original_schema_receipt'].read_text())['source_header']
                if not required.issubset(header) or header!=expected_header: raise RuntimeError('SCHEMA_GATE_FAILED')
                neff=4/(1/78566+1/329572)
                for row in reader:
                    counts['source_rows']+=1
                    if counts['source_rows']%1000000==0:
                        print('ROWS %d OUTPUT %d SECONDS %.1f'%(counts['source_rows'],counts['output_rows'],time.time()-began),flush=True)
                    if (finite(row['n']),finite(row['num_cases']),finite(row['num_controls']))!=(408138,78566,329572):
                        raise RuntimeError('PER_VARIANT_SAMPLE_COUNTS_CHANGED')
                    se_text=row['standard_error'].strip()
                    supplied=finite(se_text)
                    if se_text.lower() in {'','na','#na','nan','null','none','.'}: diagnostic['literal_missing_supplied_SE']+=1
                    elif supplied is None: raise RuntimeError('UNEXPECTED_NONMISSING_INVALID_SUPPLIED_SE')
                    else: raise RuntimeError('UNEXPECTED_FINITE_SUPPLIED_SE')
                    odds,lo,hi=map(finite,(row['odds_ratio'],row['ci_lower'],row['ci_upper']))
                    good_ci=odds is not None and odds>0 and lo is not None and hi is not None and 0<lo<hi
                    if good_ci:
                        beta=math.log(odds); se=(math.log(hi)-math.log(lo))/(2*1.959963984540054)
                        good_ci=math.isfinite(se) and se>0
                        if good_ci:
                            diagnostic['OR_outside_CI']+=not(lo<=odds<=hi)
                            center_buckets[bucket(abs(beta-(math.log(lo)+math.log(hi))/2)/se)]+=1
                    r2=finite(row['r2'])
                    if r2 is not None:
                        r2_range[0]=min(r2_range[0],r2); r2_range[1]=max(r2_range[1],r2)
                        diagnostic['R2_gt_1']+=r2>1
                    ids=set(rs for rs in re.findall(r'rs[0-9]+',row['rsid']) if rs in hm3)
                    if not ids: excluded['not_nonambiguous_HM3']+=1; continue
                    if len(ids)!=1: excluded['multiple_distinct_HM3_rsIDs']+=1; continue
                    rs=next(iter(ids)); a1,a2,chr37,bp37=hm3[rs]
                    status,mapped=chain.map_point(row['chromosome'],row['base_pair_location'])
                    if status!='mapped': excluded['chain_'+status]+=1; continue
                    chr_target,bp_target,strand=mapped
                    if (chr_target,bp_target)!=(chr37,bp37): excluded['mapped_coordinate_HM3_mismatch']+=1; continue
                    if chr37==6 and 25000000<=bp37<=34000000: excluded['extended_MHC_GRCh37']+=1; continue
                    p=finite(row['p_value'])
                    if not good_ci or p is None or not 0<=p<=1 or r2 is None or r2<=.9:
                        excluded['invalid_OR_P_CI_or_R2_le_0_9']+=1; continue
                    af=finite(row['effect_allele_frequency'])
                    if af is None or not 0<=af<=1 or min(af,1-af)<=.01: excluded['invalid_AF_or_MAF_le_0_01']+=1; continue
                    ea,oa=row['effect_allele'].upper(),row['other_allele'].upper()
                    if strand=='-': ea,oa=ea.translate(historical.COMPLEMENT),oa.translate(historical.COMPLEMENT)
                    sign=historical.orientation(ea,oa,a1,a2)
                    if sign is None: excluded['allele_mismatch_or_ambiguous']+=1; continue
                    z=beta/se
                    if not math.isfinite(z): excluded['invalid_Z']+=1; continue
                    if rs in seen: excluded['duplicate_eligible_HM3_rsID']+=1; continue
                    p_derived=math.erfc(abs(z)/math.sqrt(2))
                    p_buckets[bucket(abs(p_derived-p)/max(p,1e-300))]+=1
                    counts['retained_chain_'+strand]+=1; counts['retained_orientation_sign_'+str(sign)]+=1
                    writer.writerow([rs,a1,a2,format(sign*z,'.12g'),format(neff,'.12g')])
                    seen.add(rs); counts['output_rows']+=1
    tmp.replace(out)
    if counts['source_rows']!=19703815: raise RuntimeError('SOURCE_ROW_CARDINALITY_CHANGED')
    reconciles=counts['source_rows']==sum(excluded.values())+counts['output_rows']
    success=counts['output_rows']>=700000 and reconciles
    result={'status':'PREPROCESSING_ONLY_PASS' if success else 'PREPROCESSING_FAILED_CARDINALITY',
            'counts':dict(counts),'exclusive_exclusions':dict(excluded),'source_row_reconciliation_pass':reconciles,
            'whole_file_diagnostics':dict(diagnostic),'whole_file_log_CI_centering_SE_unit_buckets':dict(center_buckets),
            'retained_P_CI_relative_difference_buckets':dict(p_buckets),'R2_range':r2_range,
            'reference_nonambiguous_HM3_rows':len(hm3),'chain_receipt':chain_receipt,'N_eff':neff,
            'header':header,'gzip_full_source_stream_crc_verified':True,'output_path':str(out),
            'output_sha256':sha(out),'output_bytes':out.stat().st_size,'worker_code_sha256':sha(__file__),
            'execution_plan_sha256':sha(plan_path),'elapsed_seconds':time.time()-began,
            'h2_or_correlation_estimated':False,'pair_tests_admitted':0,
            'limitations':['Historical-rule CI-derived Z assumes 95% log-scale Wald CI; supplied SE absent',
                          'Coordinate and allele concordance do not independently validate upstream effect semantics',
                          'h2, intercept, power and cohort-disjointness remain pending']}
    with Path(plan['worker_receipt_path']).open('x') as f: json.dump(result,f,indent=2); f.write('\n')
    print(json.dumps({'status':result['status'],'source_rows':counts['source_rows'],'output_rows':counts['output_rows']}),flush=True)
    if not success: raise RuntimeError('CARDINALITY_GATE_FAILED_OUTPUT_PRESERVED')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--plan',required=True); main(parser.parse_args().plan)
