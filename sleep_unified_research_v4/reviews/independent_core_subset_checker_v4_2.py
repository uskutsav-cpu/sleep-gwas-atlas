#!/usr/bin/env python3
"""Independent small-output checker; does not import the candidate collator/LDSC.

Only completed receipts/output logs and 200-element delete arrays are read.
Centered deletion variance avoids using the collator's pseudovalue variance code.
Normal tail bisection independently implements the liability scale conversion.
"""
import csv
import datetime
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import os

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SEEN={}
RTOL,ATOL,PATOL=1e-12,1e-15,1e-300
ARCHIVE=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def text(path):
    SEEN[str(path)]={'sha256':sha(path),'bytes':path.stat().st_size}
    return path.read_text()

def source_json(path):return json.loads(text(path))

def near(a,b,p=False):return abs(a-b)<=max(RTOL*max(abs(a),abs(b)),PATOL if p else ATOL)

def estimates(path):
    d=[float(v) for v in text(path).split()]
    assert len(d)==200 and all(math.isfinite(v) for v in d)
    return d

def centered_se(d):
    n=len(d);mean=math.fsum(d)/n
    # SE(delete jackknife)=(n-1) * sampleSD(delete) / sqrt(n).
    ss=math.fsum((v-mean)**2 for v in d)
    return (n-1)*math.sqrt(ss/(n*(n-1)))

def normal_threshold(k):
    assert 0<k<1
    lo,hi=-12.,12.
    for _ in range(90):
        mid=(lo+hi)/2
        tail=.5*math.erfc(mid/math.sqrt(2))
        if tail>k:lo=mid
        else:hi=mid
    return (lo+hi)/2

def scale_factor(sample,pop):
    if sample is None and pop is None:return 1.
    assert sample is not None and pop is not None
    a,b=float(sample),float(pop);assert 0<a<1 and 0<b<1
    threshold=normal_threshold(b)
    density=math.exp(-.5*threshold**2)/math.sqrt(2*math.pi)
    return (b*(1-b)/density)**2/(a*(1-a))

def compatible_token(value,token):
    d=Decimal(token);quantum=Decimal(10)**d.as_tuple().exponent
    return abs(Decimal(str(value))-d)<=quantum/2

def trait(path):return Path(path).name[:-len('.sumstats.gz')]

def correlation_tokens(path):
    lines=text(path).splitlines();at=lines.index('Summary of Genetic Correlation Results')
    header=lines[at+1].split();out={}
    for line in lines[at+2:]:
        if not line.strip():break
        parts=re.split(r'\.sumstats\.gz\s+',line,maxsplit=2)
        assert len(parts)==3
        a,b=parts[:2];v=parts[2].split();assert len(v)==len(header)-2
        out[(trait(a+'.sumstats.gz'),trait(b+'.sumstats.gz'))]=dict(zip(header[2:],v))
    # Extract scalar P values separately from the formatted summary table.
    chunks=re.split(r'Computing rg for phenotype \d+/\d+\n','\n'.join(lines[:at]))[1:]
    assert len(chunks)==len(out)
    for row,chunk in zip(out.values(),chunks):
        scalar=re.search(r'^P:\s*(\S+)',chunk,re.M);assert scalar
        row['scalar_P']=scalar[1]
        row['read_SNPs']=re.search(r'Read summary statistics for (\d+) SNPs\.',chunk)[1]
        row['merged_SNPs']=re.search(r'After merging with summary statistics, (\d+) SNPs remain\.',chunk)[1]
        row['valid_allele_SNPs']=re.search(r'(\d+) SNPs with valid alleles\.',chunk)[1]
    return out

def h2_tokens(path):
    out={}
    for line in text(path).splitlines():
        if line.startswith('Total ') and 'scale h2:' in line:
            m=re.match(r'Total (\w+) scale h2:\s*(\S+)\s+\(([^)]+)\)',line);assert m
            out.update(scale=m[1].lower(),h2=m[2],se=m[3])
        elif line.startswith('Intercept:'):
            m=re.match(r'Intercept:\s*(\S+)\s+\(([^)]+)\)',line);assert m
            out.update(intercept=m[1],intercept_se=m[2])
        elif line.startswith('Ratio:'):
            out['ratio_line']=line
        elif line.startswith('Ratio <'):
            out['ratio_line']=line
    assert all(k in out for k in ['scale','h2','se','intercept','intercept_se','ratio_line'])
    return out

def verify_job(job,prefix,plan):
    rp=Path(str(prefix)+'.execution_receipt.json');r=source_json(rp)
    assert r['job']==job and r['returncode']==0 and r['scientific_cardinality_gate_pass'] and r['execution_identity_gate_pass']
    assert r['input_sha256']==r['input_sha256_after']
    assert r['dependency_sha256_before']==r['dependency_sha256_after']
    assert list(r['input_sha256'])==job['inputs']
    assert r['all_output_sha256']
    for s,h in r['all_output_sha256'].items():
        q=Path(s);assert not q.name.startswith('._')
        SEEN[str(q)]={'sha256':sha(q),'bytes':q.stat().st_size}
        assert SEEN[str(q)]['sha256']==h
    capture=source_json(Path(str(prefix)+'.full_precision.json'))
    assert SEEN[str(prefix)+'.full_precision.json']['sha256']==r['output_sha256']
    assert len(capture['estimates'])==job['estimates']
    assert capture['arguments']['out']==str(prefix)
    assert capture['arguments']['n_blocks']==200
    # Verify provenance bindings by mapping the relocated wrapper to the plan entry.
    expected=dict(plan['dependencies_sha256'])
    old=str(ROOT/'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    moved=str(Path(plan['ssd_support_package'])/'scripts/native_ldsc_capture.py')
    expected[moved]=expected.pop(old)
    assert r['dependency_sha256_before']==expected
    if job['kind']=='rg':assert capture['arguments']['rg'].split(',')==job['inputs']
    else:
        assert capture['arguments']['h2']==job['inputs'][0]
        assert capture['estimates'][0]['input']==job['inputs'][0]
    return capture

def main():
    plan=source_json(P/'manifests/ssd_native_execution_plan_v4_3.json')
    native=Path(plan['ssd_support_package'])/'native/core_reproduction_v1'
    completed=[q for q in native.glob('*.execution_receipt.json') if not q.name.startswith('._')]
    sidecars=[q for q in native.glob('*.execution_receipt.json') if q.name.startswith('._')]
    frozen_path=ROOT/'sleep_unified_research_v1/tables/original_core_396.tsv'
    frozen=list(csv.DictReader(text(frozen_path).splitlines(),delimiter='\t'))
    oldidx={(r['sleep_trait'],r['disease_trait']):r for r in frozen}
    batch=next(j for j in plan['jobs'] if j['job_id']=='core_rg_insomnia')
    prefix=native/batch['job_id'];capture=verify_job(batch,prefix,plan)
    newtokens=correlation_tokens(Path(str(prefix)+'.log'))
    oldlogs={};pairs=[]
    mapping={'rg':('rg_ratio',None),'se':('rg_se',None),'z':('z',None),'p':('p',None),'h2_obs':('tot','hsq2'),'h2_obs_se':('tot_se','hsq2'),'h2_int':('intercept','hsq2'),'h2_int_se':('intercept_se','hsq2'),'gcov_int':('intercept','gencov'),'gcov_int_se':('intercept_se','gencov')}
    for e in capture['estimates']:
        key=(trait(e['p1']),trait(e['p2']));old=oldidx[key];oldpath=ARCHIVE/old['input_log']
        if str(oldpath) not in oldlogs:oldlogs[str(oldpath)]=correlation_tokens(oldpath)
        o,n=oldlogs[str(oldpath)][key],newtokens[key]
        assert set(o)==set(n)
        stem=str(prefix)+Path(e['p1']).name+'_'+Path(e['p2']).name
        deletions={c:estimates(Path(stem+'.'+c+'.delete')) for c in ['hsq1','hsq2','gencov']}
        checks={'printed_old_new_numeric_equal':all(float(o[k])==float(n[k]) for k in o)}
        for c,d in deletions.items():checks[c+'_se']=near(centered_se(d),e[c]['tot_se'])
        ratio=e['gencov']['tot']/math.sqrt(e['hsq1']['tot']*e['hsq2']['tot'])
        rd=[c/math.sqrt(a*b) for a,b,c in zip(deletions['hsq1'],deletions['hsq2'],deletions['gencov'])]
        se=centered_se(rd);bias=ratio+199*(ratio-math.fsum(rd)/200)
        z=ratio/se;p=math.erfc(abs(z)/math.sqrt(2))
        checks.update(rg=near(ratio,e['rg_ratio']),rg_jknife=near(bias,e['rg_jknife']),rg_se=near(se,e['rg_se']),z=near(z,e['z']),p=near(p,e['p'],True),positive_delete_denominators=all(a>0 and b>0 for a,b in zip(deletions['hsq1'],deletions['hsq2'])))
        for displayed,(field,component) in mapping.items():
            quantity=e[component][field] if component else e[field]
            checks['native_capture_to_display_'+displayed]=compatible_token(quantity,n[displayed])
            # Scalar P is used in frozen table; summary table P can be zero.
            token=o['scalar_P'] if displayed=='p' else o[displayed]
            checks['frozen_to_original_log_'+displayed]=compatible_token(float(token),old[displayed])
        checks['scalar_P_display']=compatible_token(e['p'],n['scalar_P'])
        pairs.append({'pair_id':'__'.join(key),'historical_tier':old['analysis_tier'],'old_log':str(oldpath),'checks':checks,'all_pass':all(checks.values()),'independent_estimates':{'rg':ratio,'bias_corrected_rg':bias,'se':se,'z':z,'p':p},'capture':{k:e[k] for k in ['rg_ratio','rg_jknife','rg_se','z','p']}})
    h2=[]
    for ident in ['insomnia','ms','melanoma','ldl','sleepdur','t2d','parkinson']:
        job=next(j for j in plan['jobs'] if j['job_id']=='core_h2_'+ident)
        prefix=native/job['job_id'];cap=verify_job(job,prefix,plan);e=cap['estimates'][0]
        old=h2_tokens(ARCHIVE/'results/logs'/('h2_'+ident+'.log'));new=h2_tokens(Path(str(prefix)+'.log'))
        factor=scale_factor(cap['arguments']['samp_prev'],cap['arguments']['pop_prev'])
        value,error=e['tot']*factor,e['tot_se']*factor
        checks={'printed_old_new_numeric_equal':all(old[k]==new[k] if k in ['scale','ratio_line'] else float(old[k])==float(new[k]) for k in old),
                'independent_observed_SE':near(centered_se(estimates(Path(str(prefix)+'.delete'))),e['tot_se']),
                'capture_to_reported_h2':compatible_token(value,new['h2']),'capture_to_reported_SE':compatible_token(error,new['se']),
                'capture_to_intercept':compatible_token(e['intercept'],new['intercept']),'capture_to_intercept_SE':compatible_token(e['intercept_se'],new['intercept_se']),
                'scale_matches_arguments':new['scale']==('observed' if cap['arguments']['samp_prev'] is None else 'liability'),
                'h2_z_scale_invariant':near(value/error,e['tot']/e['tot_se'])}
        h2.append({'trait':ident,'checks':checks,'all_pass':all(checks.values()),'captured_observed_h2':e['tot'],'captured_observed_SE':e['tot_se'],'independent_liability_factor':factor,'reported_h2':value,'reported_SE':error,'sample_prev_argument':cap['arguments']['samp_prev'],'population_prev_argument':cap['arguments']['pop_prev'],'reported_scale':new['scale'],'reported_log':new,'negative_ratio_special_case':e['ratio']<0,'liability_estimate_above_one':new['scale']=='liability' and value>1})
    assert len(pairs)==33
    after={s:sha(Path(s))==m['sha256'] for s,m in SEEN.items()}
    result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'First completed insomnia33 pair batch and seven selected standalone h2 jobs only; no whole core stage claim','receipt_snapshot_completed_biological_jobs':len(completed),'receipt_snapshot_AppleDouble_sidecars_excluded':len(sidecars),'monitor_terminal_receipt_present':(P/'logs/core_native_monitor_receipt_v4.json').exists(),'verified_rg_pairs':len(pairs),'verified_h2_jobs':len(h2),'all_subset_arithmetic_and_display_checks_pass':all(r['all_pass'] for r in pairs+h2),'all_small_consumed_artifact_hashes_unchanged':all(after.values()),'tolerances':{'rtol':RTOL,'atol_statistics':ATOL,'atol_p':PATOL},'independent_method':'Centered delete-value sample variance; direct ratio; erfc P; inverse-normal tail solved by90bisectionsteps; decimal token half-unit compatibility','pairs':pairs,'h2':h2,'inputs':SEEN,'checker_sha256':sha(Path(__file__)),'cross_fit_block_alignment_certified':False,'source_to_estimator_full_chain_certified':False,'new_biological_result':False,'candidate_collator_imported_or_executed':False,'initial_receipt_preserved':'reviews/independent_core_subset_receipt_v4.json','initial_oracle_qualification':'Only initial frozen-to-summary literal-equality assertions failed; frozen TSV is separately serialized. This version checks compatibility within the last written decimal unit, without asserting original unrounded identity.','GWAS_sources_read':False}
    dest=P/'reviews/independent_core_subset_receipt_v4_2.json'
    assert not dest.exists();dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pairs','h2','inputs']},indent=2))

if __name__=='__main__':main()
