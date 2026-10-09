#!/usr/bin/env python3
"""Constant-memory, whole-file source diagnostics; no pair rg or preprocessing."""
import csv
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import time

PACKAGE = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()

def finite(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError): return None

def main():
    receipt_path=PACKAGE/'logs/mvp_insomnia_acquisition_receipt_v1.json'
    acquisition=json.loads(receipt_path.read_text())
    source=Path(acquisition['path'])
    if acquisition['verification_status']!='PASS' or sha(source)!=acquisition['actual_sha256']:
        raise SystemExit('SOURCE_IDENTITY_GATE_FAILED')
    before=source.stat()
    counts=Counter(); ranges={}; discrepancies=Counter(); started=time.time()
    required={'chromosome','base_pair_location','effect_allele','other_allele','odds_ratio','standard_error','effect_allele_frequency','p_value','rsid','ci_upper','ci_lower','n','num_cases','num_controls','r2'}
    with gzip.open(source,'rt',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t'); header=reader.fieldnames
        if not required.issubset(header): raise SystemExit('SCHEMA_GATE_FAILED')
        for row in reader:
            counts['source_rows']+=1
            if counts['source_rows']%1000000==0:
                print(f"SOURCE_DIAGNOSTIC_ROWS {counts['source_rows']} elapsed={time.time()-started:.1f}",flush=True)
            v={k:finite(row[k]) for k in ['odds_ratio','standard_error','effect_allele_frequency','p_value','ci_upper','ci_lower','n','num_cases','num_controls','r2']}
            for k in ['n','num_cases','num_controls','r2']:
                x=v[k]
                if x is None: counts[k+'_missing']+=1
                else:
                    r=ranges.setdefault(k,[x,x]); r[0]=min(r[0],x); r[1]=max(r[1],x)
            if v['standard_error'] is None:counts['supplied_standard_error_missing']+=1
            else:counts['supplied_standard_error_present']+=1
            nc,nt,n=v['num_cases'],v['num_controls'],v['n']
            if nc and nt and nc>0 and nt>0:
                neff=4/(1/nc+1/nt); r=ranges.setdefault('N_eff',[neff,neff]);r[0]=min(r[0],neff);r[1]=max(r[1],neff)
                if nc!=78566 or nt!=329572:counts['variant_counts_differ_from_manifest']+=1
                if n!=nc+nt:counts['n_not_cases_plus_controls']+=1
            else:counts['invalid_variant_case_control_counts']+=1
            ea,oa=row['effect_allele'],row['other_allele']
            canonical=ea in 'ACGT' and oa in 'ACGT' and len(ea)==1 and len(oa)==1 and ea!=oa
            ambiguous=canonical and {ea,oa} in [{'A','T'},{'G','C'}]
            counts['noncanonical_alleles']+=not canonical
            counts['strand_ambiguous_alleles']+=ambiguous
            af,p,r2=v['effect_allele_frequency'],v['p_value'],v['r2']
            good_maf=af is not None and 0<=af<=1 and min(af,1-af)>0.01
            good_r2=r2 is not None and r2>0.9
            counts['historical_MAF_gt_0_01']+=good_maf
            counts['historical_R2_gt_0_9']+=good_r2
            odds,lo,hi=v['odds_ratio'],v['ci_lower'],v['ci_upper']
            good_ci=odds is not None and odds>0 and lo is not None and hi is not None and 0<lo<hi
            good_p=p is not None and 0<=p<=1
            counts['invalid_OR_CI']+=not good_ci
            counts['invalid_P']+=not good_p
            if not(good_ci and good_p):continue
            beta=math.log(odds); se=(math.log(hi)-math.log(lo))/(2*1.959963984540054)
            counts['CI_derived_logOR_SE_available']+=1
            if v['standard_error'] is not None:
                rel=abs(v['standard_error']-se)/se
                counts['supplied_SE_rel_diff_gt_0_01']+=rel>0.01
            source_p=math.erfc(abs(beta/se)/math.sqrt(2))
            # Continuous descriptive buckets, no post-hoc QC threshold or row deletion.
            rel=abs(source_p-p)/max(p,1e-300)
            bucket=next((str(t) for t in [0.001,0.01,0.05,0.1,1,10] if rel<=t),'gt10')
            discrepancies['P_relative_difference_le_'+bucket]+=1
            if good_maf and good_r2 and canonical and not ambiguous:
                counts['pre_HM3_allele_mapping_eligible_rows']+=1
    after=source.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
        raise SystemExit('SOURCE_CHANGED_DURING_DIAGNOSTIC')
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.time()-started,
            'source':str(source),'source_sha256':acquisition['actual_sha256'],
            'acquisition_receipt_sha256':sha(receipt_path),'diagnostic_code_sha256':sha(Path(__file__)),
            'protocol_sha256':sha(PACKAGE/'FROZEN_NEW_ANALYSIS_PROTOCOL.md'),
            'source_header':header,'counts':dict(counts),'ranges':ranges,'P_CI_descriptive_discrepancy_buckets':dict(discrepancies),
            'gzip_full_stream_crc_verified':True,'genetic_correlation_outcomes_accessed':False,
            'harmonization_or_h2_completed':False,'pair_tests_admitted_by_this_diagnostic':False,
            'limitations':['CI-derived SE follows historical MVP rule; missing supplied SE cannot be directly validated',
                          'P/CI buckets are descriptive and affected by upstream rounding; no new exclusion threshold',
                          'HM3 identity, signed allele orientation and GRCh38 mapping are not certified',
                          'Source h2/intercept, power and human phenotype/overlap review remain pending']}
    dest=PACKAGE/'logs/mvp_insomnia_source_schema_v1.json'
    if dest.exists():raise SystemExit('PRIOR_DIAGNOSTIC_PRESERVED; use a new version')
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
