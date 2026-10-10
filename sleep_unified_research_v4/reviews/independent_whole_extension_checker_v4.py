#!/usr/bin/env python3
"""Independent full extension verification; no candidate collator/LDSC import.

Reuses sealed reviewer centered-delete arithmetic and log token utilities.
Reads only metadata, code, reference hash dependencies, stock logs and small
captured/delete outputs. No raw or processed GWAS body reads, fits or estimators.
The pinned native interpreter is used solely for NumPy scalar formatting.
"""
import csv
import datetime
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
R = P/'reviews'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
RECOVERED = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs')
PLAN = P/'manifests/ssd_native_execution_plan_v4_3.json'
PLAN_SHA = 'f555dc441e8c93c528e153d2e8689edabb29ba17e4d8b0da68402ac4e4c98e88'
MONITOR = P/'logs/extension_native_monitor_receipt_v4_3.json'
MONITOR_SHA = 'a818244e5839714166df38b1e58b2520d5a8a629992fcc13783a78cec8b22ae2'
STEM = 'independent_whole_extension'


def module(path, name):
    spec = importlib.util.spec_from_file_location(name,path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


A = module(R/'independent_partial_extension_checkpoint_checker_v4_2.py','prior_partial_reviewer_math')
U = module(R/'independent_core_subset_checker_v4_3.py','prior_core_reviewer_tokens')
U.SEEN = A.SEEN
U.text, U.sha, U.source_json, U.estimates = A.read_text, A.sha, A.read_json, A.deletion


def register(path, expected=None):
    path = Path(path)
    assert not path.name.endswith(('.sumstats.gz','.bgz','.bgz.partial')), 'GWAS_BODY_READ_FORBIDDEN'
    return A.register(path,expected)


def tab(path):
    register(path)
    return list(csv.DictReader(Path(path).read_text().splitlines(),delimiter='\t'))


def verify_seal(path):
    register(path)
    for line in Path(path).read_text().splitlines():
        digest, name = line.split('  ',1)
        register(Path(path).parent/name,digest)


def h2_tokens(path):
    d = U.h2_tokens(path)
    body = A.read_text(path)
    for k,pat in [('lambda_gc',r'^Lambda GC:\s*(\S+)'),('mean_chisq',r'^Mean Chi\^2:\s*(\S+)'),
                  ('input_snps',r'^Read summary statistics for (\d+) SNPs\.'),
                  ('regression_snps',r'^After merging with regression SNP LD, (\d+) SNPs remain\.')]:
        d[k] = re.search(pat,body,re.M)[1]
    match = re.search(r'^Ratio:\s*([\d.eE+\-]+)\s+\(([^)]+)\)',body,re.M)
    if match:
        d.update(ratio=match[1],ratio_se=match[2])
    return d


def bh(ps):
    """Independent sorted rank products and suffix minimum, full family only."""
    assert len(ps) == 1200 and all(math.isfinite(p) and 0<=p<=1 for p in ps)
    order = sorted(range(len(ps)),key=lambda i:(ps[i],i))
    scaled = [min(1.0,ps[i]*len(ps)/(j+1)) for j,i in enumerate(order)]
    qs = [None]*len(ps)
    for j,i in enumerate(order):
        qs[i] = min(scaled[j:])
    return qs


def relative_error(a,b):
    return abs(a-b)/abs(b) if b else abs(a-b)


def write_table(path, rows):
    with Path(path).open('x',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return A.sha(path)


def main():
    started = time.monotonic()
    destinations = {k:R/(STEM+'_'+k+'_v4.tsv') for k in ('jobs','h2','rg','bh')}
    receipt_path = R/(STEM+'_receipt_v4.json')
    assert all(not p.exists() for p in [*destinations.values(),receipt_path])
    assert shutil.disk_usage('/System/Volumes/Data').free >= 256*1024**2
    register(Path(__file__))
    plan = A.read_json(PLAN)
    register(PLAN,PLAN_SHA)
    register(SSD/'manifests'/PLAN.name,PLAN_SHA)
    partial = A.read_json(R/'independent_partial_extension_checkpoint_receipt_v4_2.json')
    assert partial['checker_sha256'] == register(R/'independent_partial_extension_checkpoint_checker_v4_2.py')['sha256']
    verify_seal(R/'independent_partial_extension_checkpoint_v4_2.sha256')
    core = A.read_json(R/'independent_whole_core_receipt_v4_2.json')
    register(R/'independent_whole_core_receipt_v4_2.json','72732c58a0de3fd87995e748e111d3cf24cc4b157fdde443a59a19b252850fa0')
    register(R/'independent_core_subset_checker_v4_3.py',core['inputs'][str(R/'independent_core_subset_checker_v4_3.py')]['sha256'])
    original_path = ROOT/'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py'
    register(original_path,plan['original_runner_sha256'])
    original = module(original_path,'original_config_reader_only')
    original_plan_path = ROOT/'sleep_unified_research_v1/manifests/native_reproduction_jobs_v1.json'
    register(original_plan_path,plan['original_plan_sha256'])
    assert A.read_json(original_plan_path)['jobs'] == original.jobs() == plan['jobs'] and len(plan['jobs']) == 190
    for path,expected in {**plan['dependencies_sha256'],**plan['support_file_sha256']}.items():
        register(path,expected)
    A.REFERENCE_M = math.fsum(float(A.read_text(original.REF/(str(c)+'.l2.M_5_50')).strip()) for c in range(1,23))
    assert A.REFERENCE_M == 1173569
    support = Path(plan['ssd_support_package'])
    native = support/'native/extension_reproduction_v1'
    jobs = [j for j in plan['jobs'] if j['stage']=='extension']
    assert len(jobs)==112 and sum(j['kind']=='h2' for j in jobs)==100 and sum(j['estimates'] for j in jobs if j['kind']=='rg')==1200
    dependencies = dict(plan['dependencies_sha256'])
    old_wrapper = str(ROOT/'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    moved_wrapper = str(support/'scripts/native_ldsc_capture.py')
    assert old_wrapper in dependencies and moved_wrapper not in dependencies
    dependencies[moved_wrapper] = dependencies.pop(old_wrapper)
    admitted = {r['path']:r['actual_sha256'] for r in plan['inputs_verified']}
    assert len(admitted)==158 and all(r['match'] and r['expected_sha256']==r['actual_sha256'] for r in plan['inputs_verified'])
    ledger = {}
    for name,predicate in [('native_input_hash_checks.tsv',lambda r:r['kind']=='core_munged'),
                           ('archived_extension_recovery.tsv',lambda r:r['status']=='EXACT_RECEIPT_HASH_RECOVERED')]:
        ledger.update({r['path']:r['actual_sha256'] for r in tab(support/'tables'/name) if predicate(r)})
    assert all(ledger[path]==digest for path,digest in admitted.items())
    register(MONITOR,MONITOR_SHA)
    monitor = A.read_json(MONITOR)
    assert monitor['stage']=='extension' and monitor['plan_sha256']==PLAN_SHA
    assert monitor['returncode']==0 and monitor['stop_reason'] is None and monitor['owned_cleanup_verified']
    assert monitor['process_group_teardown']['remaining_group_members']==[]
    assert monitor['unchanged_successful_jobs_reused']==107 and monitor['remaining_frozen_jobs']==5
    assert monitor['old_failed_attempt_reclassified'] is False
    state = monitor['final_resource_snapshot']
    assert state['internal_free_bytes']>=plan['internal_floor_bytes'] and state['ssd_free_bytes']>=plan['ssd_floor_bytes']
    continuation_path = P/'manifests/native_extension_checkpoint_continuation_v4_3.json'
    register(continuation_path,monitor['continuation_plan_sha256'])
    continuation = A.read_json(continuation_path)
    assert continuation['historical_plan_sha256']==PLAN_SHA and continuation['original190_jobs']==plan['jobs']
    assert continuation['completed_jobs']+continuation['remaining_jobs']==jobs
    assert len(continuation['completed_jobs'])==107 and len(continuation['remaining_jobs'])==5
    for path,h in continuation['dependencies_sha256'].items():
        register(path,h)
    register(continuation['prior_failed_monitor'],monitor['prior_failed_monitor_sha256'])
    failed_monitor = A.read_json(continuation['prior_failed_monitor'])
    assert failed_monitor['returncode']==-15 and failed_monitor['stop_reason']=='INTERNAL_FULL_NATIVE_FLOOR_REACHED'
    assert failed_monitor['process_group_teardown']['remaining_group_members']==[]
    admission_path = P/'manifests/native_extension_checkpoint_continuation_admission_v4_3.json'
    register(admission_path,monitor['admission_sha256'])
    admission = A.read_json(admission_path)
    assert admission['execution_admitted'] is True and admission['plan_sha256']==monitor['continuation_plan_sha256']
    register(P/'scripts/61_native_checkpoint_continuation_v3.py',admission['executor_sha256'])
    for path,h in admission['independent_review_sha256'].items():
        register(path,h)
    register(monitor['worker_receipt'],monitor['worker_receipt_sha256'])
    worker = A.read_json(monitor['worker_receipt'])
    assert worker['status']=='WORKER_COMPLETE_VERIFIED' and worker['returncode']==0 and worker['stop_reason'] is None
    assert worker['plan_sha256']==monitor['continuation_plan_sha256'] and worker['owned_cleanup_verified']
    assert worker['process_group_teardown']['remaining_group_members']==[]
    for path,h in worker['output_sha256'].items():
        register(path,h)
    stdout_path = next(path for path in worker['output_sha256'] if path.endswith('.stdout.log'))
    stdout = A.read_text(stdout_path)
    skipped = re.findall(r'^CHECKPOINT_VERIFIED\s+\d+/112\s+(\S+)',stdout,re.M)
    started_jobs = re.findall(r'^NATIVE_START\s+\d+/112\s+(\S+)',stdout,re.M)
    completed_jobs = re.findall(r'^NATIVE_COMPLETE\s+(\S+)\s+elapsed_seconds=',stdout,re.M)
    assert skipped == [j['job_id'] for j in continuation['completed_jobs']]
    assert started_jobs == completed_jobs == [j['job_id'] for j in continuation['remaining_jobs']]
    preservation = None
    for path,h in monitor['interrupted_preservation_sha256'].items():
        register(path,h)
        if path.endswith('.json'):
            preservation = A.read_json(path)
    assert preservation['status']=='ALL284_BIOLOGICAL_AND284_SIDECAR_ARTIFACTS_PRESERVED_WITHOUT_DELETION'
    assert len(preservation['artifacts'])==len(preservation['transport_sidecars'])==284
    assert preservation['sidecars_count_as_biological_estimates'] is False
    for item in preservation['artifacts']:
        register(item['preserved_path'],item['sha256'])
        assert Path(item['preserved_path']).stat().st_size==item['bytes']
    # Sidecar contents are outside numerical evidence; the hash-bound specialist
    # transport proof is retained, and sidecars never count as jobs/estimates.
    receipts = [p for p in native.glob('*.execution_receipt.json') if not p.name.startswith('._')]
    sidecars = [p for p in native.glob('*.execution_receipt.json') if p.name.startswith('._')]
    assert len(receipts)==112 and set(map(str,receipts))==set(monitor['native_receipt_sha256'])
    assert set(receipts)=={Path(str(native/j['job_id'])+'.execution_receipt.json') for j in jobs}
    frozen_path = ROOT/'sleep_unified_research_v1/tables/original_extension_1200.tsv'
    frozen_h2_path = ROOT/'discovery_extension/results/ldsc/extension_trait_readiness.tsv'
    frozen = tab(frozen_path)
    old_idx = {(r['sleep_trait'],r['extension_trait_id']):r for r in frozen}
    old_h2 = {r['extension_trait_id']:r for r in tab(frozen_h2_path)}
    rg_table_path = P/'tables/extension_native_full_precision_rg_v4.tsv'
    h2_table_path = P/'tables/extension_native_full_precision_h2_v4.tsv'
    new_pairs = tab(rg_table_path)
    new_idx = {(r['sleep_trait'],r['outcome_trait']):r for r in new_pairs}
    new_h2_rows = tab(h2_table_path)
    h2_idx = {r['trait_id']:r for r in new_h2_rows}
    verification = tab(P/'tables/extension_native_job_verification_v4.tsv')
    verify_idx = {r['job_id']:r for r in verification}
    arithmetic = tab(P/'statistical_validation/extension_native_block_arithmetic_v4.tsv')
    ar_idx = {(r['kind'],r['identity']):r for r in arithmetic}
    assert len(frozen)==len(old_idx)==len(new_pairs)==len(new_idx)==1200 and set(old_idx)==set(new_idx)
    assert len(new_h2_rows)==len(h2_idx)==len(old_h2)==100 and set(old_h2)==set(h2_idx)
    assert len(verification)==len(verify_idx)==112 and len(arithmetic)==len(ar_idx)==1300
    job_rows,h2_rows,pair_rows,formats,values = [],[],[],[],[]
    logs = {}
    array_count = 0
    fields = {'rg':('rg_ratio',None),'se':('rg_se',None),'z':('z',None),'p':('p',None),
              'h2_obs':('tot','hsq2'),'h2_obs_se':('tot_se','hsq2'),'h2_int':('intercept','hsq2'),
              'h2_int_se':('intercept_se','hsq2'),'gcov_int':('intercept','gencov'),'gcov_int_se':('intercept_se','gencov')}
    for job in jobs:
        A.guard()
        prefix = native/job['job_id']
        rp = Path(str(prefix)+'.execution_receipt.json')
        register(rp,monitor['native_receipt_sha256'][str(rp)])
        if str(rp) in continuation['completed_receipt_sha256']:
            register(rp,continuation['completed_receipt_sha256'][str(rp)])
        r = A.read_json(rp)
        assert r['job']==job and r['returncode']==0 and r['scientific_cardinality_gate_pass'] is True
        assert r['execution_identity_gate_pass'] is True and r['original_outputs_modified'] is False
        assert r['command']==A.expected_command(job,prefix,original,support)
        assert r['input_sha256']==r['input_sha256_after']=={path:admitted[path] for path in job['inputs']}
        assert r['dependency_sha256_before']==r['dependency_sha256_after']==dependencies
        files = [p for p in prefix.parent.glob(prefix.name+'*') if p.is_file() and not p.name.startswith('._')]
        assert set(map(str,files))==set(r['all_output_sha256'])|{str(rp)}
        for path,h in r['all_output_sha256'].items():
            assert Path(path).parent==native and Path(path).name.startswith(prefix.name)
            register(path,h)
        cp = Path(str(prefix)+'.full_precision.json')
        register(cp,r['output_sha256'])
        cap = A.read_json(cp)
        args = cap['arguments']
        assert cap['ldsc_dir']==str(ROOT.parent/'ldsc-code')
        assert cap['python'].startswith('3.9.23') and cap['libraries']=={'numpy':'1.21.5','pandas':'1.3.3','scipy':'1.7.3'}
        assert cap['instrumentation']=='return-value capture; numerical estimator unmodified'
        assert args['out']==str(prefix) and args['n_blocks']==200 and args['print_delete_vals'] is True
        assert args['ref_ld_chr']==args['w_ld_chr']==str(original.REF)+'/'
        assert args['no_check_alleles'] is False and args['no_intercept'] is False
        assert args['samp_prev'] is None and args['pop_prev'] is None
        assert len(cap['estimates'])==job['estimates'] and original.finite_estimates(job['kind'],cap['estimates'],job['estimates'])
        arrays = [p for p in files if p.suffix in ('.delete','.part_delete')]
        count = 2 if job['kind']=='h2' else 3*job['estimates']
        assert len(arrays)==count
        array_count += count
        v = verify_idx[job['job_id']]
        assert v['kind']==job['kind'] and int(v['estimates'])==job['estimates'] and v['receipt_path']==str(rp)
        assert v['receipt_sha256']==A.SEEN[str(rp)]['sha256'] and v['full_precision_sha256']==r['output_sha256']
        assert all(v[k]=='True' for k in ['all_output_hashes_verified','input_identity_before_after','execution_identity_gate_pass'])
        assert float(v['elapsed_seconds'])==r['elapsed_seconds']
        job_rows.append({'job_id':job['job_id'],'kind':job['kind'],'estimates':job['estimates'],'receipt_sha256':A.SEEN[str(rp)]['sha256'],
                         'capture_sha256':r['output_sha256'],'command_input_dependency_output_bindings_pass':True,
                         'reused_success_preserved':str(rp) in continuation['completed_receipt_sha256'],'delete_arrays_200':count})
        if job['kind']=='h2':
            assert args['h2']==job['inputs'][0] and args['rg'] is None and cap['estimates'][0]['input']==job['inputs'][0]
            ident,e = U.trait(job['inputs'][0]),cap['estimates'][0]
            old,new = old_h2[ident],h2_idx[ident]
            old_path = RECOVERED/old['input_log']
            o,n = h2_tokens(old_path),h2_tokens(Path(str(prefix)+'.log'))
            d,part = A.deletion(str(prefix)+'.delete'),A.deletion(str(prefix)+'.part_delete')
            se = A.centered_se(d)
            factor = U.scale_factor(args['samp_prev'],args['pop_prev'])
            checks = {'printed_old_native_numeric_equal':all(o[k]==n[k] if k in ['scale','ratio_line'] else float(o[k])==float(n[k]) for k in o),
                      'observed_200_delete_se':A.near(se,e['tot_se']),
                      'partition_delete_times_M':all(A.near(a,b*A.REFERENCE_M) for a,b in zip(d,part)),
                      'observed_scale_no_prevalence_conversion':o['scale']==n['scale']==new['scale']=='observed' and factor==1.0,
                      'prevalence_arguments_table_empty':new['sample_prevalence_argument']==new['population_prevalence_argument']=='',
                      'source_identity_fields':new['original_log']==str(old_path) and new['original_log_sha256']==A.SEEN[str(old_path)]['sha256'] and new['native_full_precision_path']==str(cp) and new['native_full_precision_sha256']==r['output_sha256'] and new['original_table']==str(frozen_h2_path),
                      'historical_fields_preserved':all(new['original_'+k]==value for k,value in old.items()),
                      'reproduction_status':new['reproduction_status']=='PRINTED_PRECISION_CONCORDANT' and new['historical_full_precision_available']=='False'}
            expected = {'h2_observed_full_precision':e['tot'],'h2_observed_se_full_precision':e['tot_se'],
                        'reported_h2_full_precision':e['tot']*factor,'reported_h2_se_full_precision':e['tot_se']*factor,'liability_factor':factor,
                        'h2_z_full_precision':e['tot']/e['tot_se'],'intercept_full_precision':e['intercept'],'intercept_se_full_precision':e['intercept_se'],
                        'lambda_gc_full_precision':e['lambda_gc'],'mean_chisq_full_precision':e['mean_chisq'],'ratio_full_precision':e['ratio'],'ratio_se_full_precision':e['ratio_se']}
            for k,value in expected.items():
                checks['table_'+k] = new[k]=='' if value is None else float(new[k])==value
            checks['table_input_snps'] = int(new['input_snp_count'])==int(n['input_snps'])
            checks['table_regression_snps'] = int(new['regression_snp_count'])==int(n['regression_snps'])
            checks['table_h2_Z_gate'] = (new['h2_z_ge_4_full_precision']=='True')==(e['tot']/e['tot_se']>=4)
            checks['table_intercept_gate'] = (new['intercept_le_1_2_full_precision']=='True')==(e['intercept']<=1.2)
            checks['observed_physical_flag'] = new['physical_liability_point_in_0_1']=='NOT_APPLICABLE'
            ar = ar_idx[('h2',ident)]
            checks['collator_arithmetic_row'] = ar['status']=='INDEPENDENT_200_BLOCK_H2_ARITHMETIC_PASS' and A.near(se,float(ar['independent_block_se'])) and float(ar['captured_se'])==e['tot_se']
            fmt = {'h2':e['tot']*factor,'se':e['tot_se']*factor,'intercept':e['intercept'],'intercept_se':e['intercept_se'],'lambda_gc':e['lambda_gc'],'mean_chisq':e['mean_chisq']}
            if 'ratio' in n:
                fmt.update(ratio=e['ratio'],ratio_se=e['ratio_se'])
            else:
                checks['ratio_qualifier_consistent'] = ('NA' in n['ratio_line'] and e['mean_chisq']<=1) or ('Ratio <' in n['ratio_line'] and e['ratio']<0 and e['mean_chisq']>1)
            for field,value in fmt.items():
                values.append(value);formats.append({'identity':ident,'kind':'h2','field':field,'expected_token':n[field]})
            for field,token in [('h2','h2'),('h2_se','se'),('LDSC_intercept','intercept'),('LDSC_intercept_se','intercept_se'),('lambda_gc','lambda_gc'),('mean_chi2','mean_chisq')]:
                checks['frozen_4g_'+field] = format(float(old[field]),'.4g')==format(float(o[token]),'.4g')
            h2_rows.append({'trait_id':ident,'observed_h2':e['tot'],'independent_centered_delete_se':se,'captured_se':e['tot_se'],
                            'reported_scale':n['scale'],'liability_factor':factor,'checks_count':len(checks),'all_pass':all(checks.values()),
                            'failed_checks':','.join(k for k,v in checks.items() if not v)})
            continue
        assert args['rg'].split(',')==job['inputs'] and args['h2'] is None
        assert [(e['p1'],e['p2']) for e in cap['estimates']]==[(job['inputs'][0],p2) for p2 in job['inputs'][1:]]
        nt = U.correlation_tokens(Path(str(prefix)+'.log'))
        assert len(nt)==100
        for e in cap['estimates']:
            key = U.trait(e['p1']),U.trait(e['p2'])
            old,new = old_idx[key],new_idx[key]
            old_path = RECOVERED/old['input_log']
            if str(old_path) not in logs:
                logs[str(old_path)] = U.correlation_tokens(old_path)
            o,n = logs[str(old_path)][key],nt[key]
            stem = str(prefix)+Path(e['p1']).name+'_'+Path(e['p2']).name
            deletes = {c:A.deletion(stem+'.'+c+'.delete') for c in ('hsq1','hsq2','gencov')}
            checks = {'printed_old_native_numeric_equal':all(float(o[k])==float(n[k]) for k in o),
                      'positive_point_and_delete_denominators':e['hsq1']['tot']>0 and e['hsq2']['tot']>0 and all(a>0 and b>0 for a,b in zip(deletes['hsq1'],deletes['hsq2']))}
            assert checks['positive_point_and_delete_denominators']
            for c,d in deletes.items():
                checks[c+'_centered_delete_se'] = A.near(A.centered_se(d),e[c]['tot_se'])
            ratio = e['gencov']['tot']/math.sqrt(e['hsq1']['tot']*e['hsq2']['tot'])
            rd = [c/math.sqrt(a*b) for a,b,c in zip(deletes['hsq1'],deletes['hsq2'],deletes['gencov'])]
            se = A.centered_se(rd)
            bias = ratio+199*(ratio-math.fsum(rd)/200)
            z,p = ratio/se,math.erfc(abs(ratio/se)/math.sqrt(2))
            for k,value in {'rg_ratio':ratio,'rg_jknife':bias,'rg_se':se,'z':z,'p':p}.items():
                checks['independent_'+k] = A.near(value,e[k],k=='p')
            checks['P_from_captured_Z'] = A.near(math.erfc(abs(e['z'])/math.sqrt(2)),e['p'],True)
            for field,capfield in [('rg','rg_ratio'),('se','rg_se'),('z','z'),('p','p'),('rg_jknife_bias_corrected','rg_jknife')]:
                checks['table_'+field] = float(new[field])==e[capfield]
            checks['table_CI_lower'] = A.near(float(new['ci_lower_95']),e['rg_ratio']-1.959963984540054*e['rg_se'])
            checks['table_CI_upper'] = A.near(float(new['ci_upper_95']),e['rg_ratio']+1.959963984540054*e['rg_se'])
            for display,(field,component) in fields.items():
                value = e[component][field] if component else e[field]
                checks['exact_summary_4f_'+display] = format(value,'.4f')==n[display]
            frozen_map = {'rg':'rg','se':'se','z':'z','p':'scalar_P','extension_h2_observed':'h2_obs','extension_h2_observed_se':'h2_obs_se',
                          'extension_h2_intercept':'h2_int','extension_h2_intercept_se':'h2_int_se','cross_trait_LDSC_intercept':'gcov_int','cross_trait_LDSC_intercept_se':'gcov_int_se'}
            for field,token in frozen_map.items():
                checks['frozen_declared_4g_'+field] = format(float(old[field]),'.4g')==format(float(o[token]),'.4g')
            for c in ('hsq1','hsq2','gencov'):
                for k,value in e[c].items():
                    checks['table_'+c+'_'+k] = new[c+'_'+k]=='' if value is None else float(new[c+'_'+k])==value
            for dest,token in [('input_snp_count','read_SNPs'),('snp_overlap_after_merge','merged_SNPs'),('snp_overlap_valid_alleles','valid_allele_SNPs')]:
                checks['table_'+dest] = int(new[dest])==int(n[token])
            checks['historical_fields_preserved'] = all(new['original_'+k]==value for k,value in old.items())
            checks['source_identity_fields'] = new['original_log']==str(old_path) and new['original_log_sha256']==A.SEEN[str(old_path)]['sha256'] and new['native_full_precision_path']==str(cp) and new['native_full_precision_sha256']==r['output_sha256'] and new['original_table']==str(frozen_path)
            checks['status_flags'] = new['reproduction_status']=='PRINTED_PRECISION_CONCORDANT' and new['historical_full_precision_available']=='False' and new['captured_full_precision_to_native_summary_pass']=='True'
            checks['printed_field_counts'] = int(new['printed_fields_compared'])==int(new['printed_fields_concordant'])==14
            checks['pairwise_h2_Z_gate'] = (new['pairwise_h2_z_ge_4_both']=='True')==all(e[c]['tot']/e[c]['tot_se']>=4 for c in ('hsq1','hsq2'))
            checks['pairwise_intercept_gate'] = (new['pairwise_h2_intercept_le_1_2_both']=='True')==all(e[c]['intercept']<=1.2 for c in ('hsq1','hsq2'))
            for field,capfield in [('rg','rg_ratio'),('se','rg_se'),('z','z'),('p','p')]:
                checks['serialized_difference_'+field] = float(new[capfield+'_native_minus_original_serialized'])==e[capfield]-float(old[field])
            ar = ar_idx[('rg','__'.join(key))]
            checks['collator_arithmetic_row'] = ar['status']=='INDEPENDENT_200_BLOCK_RATIO_ARITHMETIC_PASS' and ar['capture_to_native_summary_pass']=='True' and ar['cross_estimator_boundaries_certified']=='False'
            for field in ('hsq1_se','hsq2_se','gencov_se','rg_ratio','rg_jknife','rg_se','z','p'):
                checks['collator_arithmetic_'+field] = ar[field]=='True'
            values.append(e['p']);formats.append({'identity':'__'.join(key),'kind':'rg','field':'scalar_P','expected_token':n['scalar_P']})
            pair_rows.append({'sleep_trait':key[0],'outcome_trait':key[1],'independent_ratio':ratio,'independent_bias_corrected_rg':bias,
                              'independent_centered_se':se,'captured_se':e['rg_se'],'independent_z':z,'independent_p':p,'captured_p':e['p'],
                              'p_relative_error':relative_error(p,e['p']),'checks_count':len(checks),'all_pass':all(checks.values()),
                              'failed_checks':','.join(k for k,v in checks.items() if not v)})
    assert len(job_rows)==112 and len(h2_rows)==100 and len(pair_rows)==1200 and array_count==3800
    formatter = 'import json,sys,numpy as np;np.set_printoptions(linewidth=1000,precision=4);x=json.load(sys.stdin);print(json.dumps(dict(numpy=np.__version__,values=[str(np.matrix(v)).replace("[", "").replace("]", "").strip() for v in x])))'
    env = {**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','TMPDIR':str(SSD/'tmp')}
    done = subprocess.run([str(original.PYTHON),'-B','-c',formatter],input=json.dumps(values),capture_output=True,text=True,check=True,timeout=30,env=env)
    rendered = json.loads(done.stdout)
    assert rendered['numpy']=='1.21.5' and len(rendered['values'])==len(formats)
    for row,token in zip(formats,rendered['values']):
        row.update(actual_token=token,exact_match=token==row['expected_token'])
    native_q = bh([r['captured_p'] for r in pair_rows])
    frozen_q = bh([float(old_idx[(r['sleep_trait'],r['outcome_trait'])]['p']) for r in pair_rows])
    bh_rows = []
    for pair,q,oq in zip(pair_rows,native_q,frozen_q):
        key = pair['sleep_trait'],pair['outcome_trait']
        old,new = old_idx[key],new_idx[key]
        old_recorded_q = float(old['extension_fdr'])
        native_pass,frozen_pass = q<.05,old_recorded_q<.05
        changed = native_pass!=frozen_pass
        checks = {'native_q_table':A.near(q,float(new['native_fdr_original_complete_family']),True),
                  'frozen_q_independent':A.near(oq,old_recorded_q,True),
                  'frozen_q_table':float(new['frozen_fdr_original_complete_family'])==old_recorded_q,
                  'native_status':(new['native_fdr_pass_0_05']=='True')==native_pass,
                  'frozen_status':(new['frozen_fdr_pass_0_05']=='True')==frozen_pass,
                  'boundary_change_status':(new['fdr_boundary_status_changed']=='True')==changed}
        bh_rows.append({'sleep_trait':key[0],'outcome_trait':key[1],'independent_native_q':q,'independent_frozen_q':oq,
                        'recorded_frozen_q':old_recorded_q,'native_positive':native_pass,'frozen_positive':frozen_pass,
                        'boundary_changed':changed,'all_pass':all(checks.values()),'failed_checks':','.join(k for k,v in checks.items() if not v)})
    comparison = A.read_json(P/'logs/extension_native_comparison_receipt_v4.json')
    register(P/'scripts/33_compare_native_campaign_v4_4.py',comparison['collator_sha256'])
    for path,info in comparison['sources_before'].items():
        register(path,info['sha256'])
        assert Path(path).stat().st_size==info['bytes']
    native_positive = sum(r['native_positive'] for r in bh_rows)
    frozen_positive = sum(r['frozen_positive'] for r in bh_rows)
    changes = sum(r['boundary_changed'] for r in bh_rows)
    independently_pass = all(r['all_pass'] for r in h2_rows+pair_rows+bh_rows) and all(r['exact_match'] for r in formats)
    assert comparison['job_count']==112 and comparison['h2_count']==100 and comparison['rg_count']==1200
    assert comparison['monitor_sha256']==MONITOR_SHA and comparison['plan_sha256']==PLAN_SHA
    assert comparison['native_fdr_positives']==native_positive and comparison['frozen_fdr_positives']==frozen_positive and comparison['fdr_boundary_status_changes']==changes
    assert comparison['arithmetic_tolerance']=={'relative':A.RTOL,'absolute_estimates':A.ATOL,'absolute_p':A.PATOL}
    assert not comparison['arithmetic_failures'] and not comparison['meaningful_printed_discrepancies']
    unchanged = all(A.sha(path)==info['sha256'] and Path(path).stat().st_size==info['bytes'] for path,info in A.SEEN.items())
    assert unchanged
    A.guard()
    tables = {key:{'path':str(destinations[key]),'sha256':write_table(destinations[key],rows)} for key,rows in [('jobs',job_rows),('h2',h2_rows),('rg',pair_rows),('bh',bh_rows)]}
    result = {'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'status':'PASS_WHOLE_EXTENSION_PROCESSED_INPUT_NATIVE_REPRODUCTION' if independently_pass and changes==0 else 'REVIEW_REQUIRED_EXTENSION_DISCREPANCY',
              'all_independent_checks_pass':independently_pass,'job_count':112,'standalone_h2_count':100,'rg_count':1200,'finite200_delete_array_count':array_count,
              'monitor_terminal_success_verified':True,'monitor_sha256':MONITOR_SHA,'plan_sha256':PLAN_SHA,
              '107_success_receipts_unchanged':True,'5_remaining_original_commands_completed':started_jobs,
              'continuation_stdout_skips':len(skipped),'continuation_stdout_starts':len(started_jobs),'continuation_stdout_completes':len(completed_jobs),
              'interrupted_snoring_biological_artifacts_preserved':284,'interrupted_snoring_accepted_partial_estimates':0,
              'excluded_native_receipt_AppleDouble_sidecars':len(sidecars),'biological_receipts':len(receipts),
              'historical_full_precision_available':False,'historical_log_numeric_concordance_h2':sum(r['all_pass'] for r in h2_rows),
              'historical_log_numeric_concordance_rg':sum(r['all_pass'] for r in pair_rows),
              'native_summary_exact_4f_checks':12000,'pinned_numpy_rg_scalar_P_checks':1200,
              'pinned_numpy_h2_scalar_checks':sum(r['kind']=='h2' for r in formats),'formatter_numpy_version':rendered['numpy'],
              'format_failures':[r for r in formats if not r['exact_match']],
              'numerical_or_table_failures':[r for r in h2_rows+pair_rows+bh_rows if not r['all_pass']],
              'native_complete1200_BH_positives':native_positive,'frozen_complete1200_BH_positives':frozen_positive,'BH_boundary_status_changes':changes,
              'highest_centered_reconstructed_P_relative_error':max(r['p_relative_error'] for r in pair_rows),
              'nearest_native_q_to_0_05':min(bh_rows,key=lambda r:abs(r['independent_native_q']-.05)),
              'all_h2_reported_scales_observed':all(r['reported_scale']=='observed' and r['liability_factor']==1 for r in h2_rows),
              'reference_M_5_50_total':A.REFERENCE_M,'tolerances':{'rtol':A.RTOL,'atol_statistics':A.ATOL,'atol_P':A.PATOL},
              'arithmetic_method':'Unchanged sealed reviewer centered 200-delete variance, direct covariance/sqrt(h2_1*h2_2), bias correction from deletion mean, erfc normal P; independent full-family rank/suffix-min BH.',
              'source_identity_scope':'158 admitted input identities reconciled to frozen ledgers and all112 per-job before/after hashes; no fresh raw or munged GWAS body hashes or decompression.',
              'limitations':['Printed original-log precision concordance does not establish unavailable archival unrounded identity.',
                             'Standalone extension h2 is observed-scale under original supplied N; effective N/ascertainment/absolute variance-fraction validity remain uncertified.',
                             'Output delete arrays verify within-fit arithmetic only; cross-fit SNP/mask/genomic alignment and corrected sampling covariance are not certified.',
                             'Raw-to-harmonized-to-munged source chain and scientific calibration remain separately assessed.',
                             'Core396, extension1200, validation41 and original217 families stay separate; validation completion and all190 admission are outside this review.',
                             'No new biological findings or execution admission arise from this reproducibility audit.'],
              'candidate_collator_imported_or_executed':False,'GWAS_bodies_read':False,'native_estimators_or_audits_launched':0,
              'formatter_subprocess_only':True,'formatter_code':formatter,
              'tables':tables,'consumed_artifacts':A.SEEN,'all_consumed_artifact_hashes_unchanged':unchanged,
              'resource_plan':{'hash_buffer_bytes':65536,'launch_internal_free_minimum':256*1024**2,'emergency_internal_floor':A.MIN_FREE,'maximum_RSS_bytes':A.MAX_RSS,'deadline_seconds':A.DEADLINE},
              'resource_usage':{'elapsed_seconds':time.monotonic()-started,'peak_self_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
              'checker_sha256':A.sha(Path(__file__))}
    with receipt_path.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({k:result[k] for k in ['status','all_independent_checks_pass','job_count','standalone_h2_count','rg_count','finite200_delete_array_count','native_complete1200_BH_positives','frozen_complete1200_BH_positives','BH_boundary_status_changes','highest_centered_reconstructed_P_relative_error','resource_usage']}))
    if not independently_pass or changes:
        raise SystemExit('INDEPENDENT_EXTENSION_REVIEW_REQUIRED; evidence preserved')


if __name__=='__main__':
    main()
