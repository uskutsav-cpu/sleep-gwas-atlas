#!/usr/bin/env python3
"""Independent full validation verification; no candidate collator/LDSC import.

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
MONITOR = P/'logs/validation_native_monitor_receipt_v4.json'
MONITOR_SHA = '751ae05e90048e638e43fa7c6f37f384faef88306c105f74ea618797e4c8e3d1'
STEM = 'independent_whole_validation'


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


def relative_error(a,b):
    return abs(a-b)/abs(b) if b else abs(a-b)


def write_table(path, rows):
    with Path(path).open('x',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return A.sha(path)


def exact_owned_outputs(job, prefix):
    # ARTHROSIS is a literal prefix of ARTHROSIS_KNEE: broad globs would
    # incorrectly attribute one valid job's outputs to the other job.
    paths = {Path(str(prefix)+suffix) for suffix in ('.log','.stdout.log','.full_precision.json','.execution_receipt.json')}
    if job['kind']=='h2':
        paths.update(Path(str(prefix)+suffix) for suffix in ('.delete','.part_delete'))
    else:
        for p2 in job['inputs'][1:]:
            stem = str(prefix)+Path(job['inputs'][0]).name+'_'+Path(p2).name
            paths.update(Path(stem+'.'+component+'.delete') for component in ('hsq1','hsq2','gencov'))
    assert all(path.is_file() for path in paths)
    return sorted(paths)


def main():
    started = time.monotonic()
    destinations = {k:R/(STEM+'_'+k+'_v4.tsv') for k in ('jobs','h2','rg','classification217','aggregate190')}
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
    native = support/'native/validation_reproduction_v1'
    jobs = [j for j in plan['jobs'] if j['stage']=='validation']
    assert len(jobs)==21 and sum(j['kind']=='h2' for j in jobs)==13 and sum(j['estimates'] for j in jobs if j['kind']=='rg')==41
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
    assert monitor['stage']=='validation' and monitor['plan_sha256']==PLAN_SHA
    assert monitor['returncode']==0 and monitor['stop_reason'] is None
    assert monitor['process_group_teardown']['remaining_group_members']==[]
    state = monitor['final_resource_snapshot']
    assert state['internal_free_bytes']>=plan['internal_floor_bytes'] and state['ssd_free_bytes']>=plan['ssd_floor_bytes']
    register(monitor['stdout_path'],monitor['stdout_sha256'])
    stdout = A.read_text(monitor['stdout_path'])
    started_jobs = re.findall(r'^NATIVE_START\s+\d+/21\s+(\S+)',stdout,re.M)
    completed_jobs = re.findall(r'^NATIVE_COMPLETE\s+(\S+)\s+elapsed_seconds=',stdout,re.M)
    assert started_jobs==completed_jobs==[j['job_id'] for j in jobs]
    ext_receipt_path = R/'independent_whole_extension_receipt_v4.json'
    register(ext_receipt_path,'9fe497d26046b9dc629cb18479886543f55f00912a493d44f1b5379662ab43aa')
    ext = A.read_json(ext_receipt_path)
    verify_seal(R/'independent_whole_extension_v4.sha256')
    assert core['all_independent_checks_pass'] and ext['all_independent_checks_pass']
    core_jobs = {j['job_id']:j for j in core['jobs']}
    ext_jobs = {j['job_id']:j for j in tab(R/'independent_whole_extension_jobs_v4.tsv')}
    aggregate = []
    for j in plan['jobs']:
        prefix = support/'native'/(j['stage']+'_reproduction_v1')/j['job_id']
        rp = Path(str(prefix)+'.execution_receipt.json')
        r = A.read_json(rp)
        proof = core_jobs.get(j['job_id']) if j['stage']=='core' else ext_jobs.get(j['job_id']) if j['stage']=='extension' else None
        if proof is not None:
            register(rp,proof['receipt_sha256'])
            register(str(prefix)+'.full_precision.json',proof['capture_sha256'])
        assert r['job']==j and r['returncode']==0 and r['scientific_cardinality_gate_pass'] is True and r['execution_identity_gate_pass'] is True
        assert r['original_outputs_modified'] is False and r['command']==A.expected_command(j,prefix,original,support)
        assert r['input_sha256']==r['input_sha256_after']=={path:admitted[path] for path in j['inputs']}
        assert r['dependency_sha256_before']==r['dependency_sha256_after']==dependencies
        files = exact_owned_outputs(j,prefix)
        assert set(map(str,files))==set(r['all_output_sha256'])|{str(rp)}
        for path,h in r['all_output_sha256'].items():
            register(path,h)
        cp = Path(str(prefix)+'.full_precision.json')
        register(cp,r['output_sha256'])
        cap = A.read_json(cp)
        assert len(cap['estimates'])==j['estimates'] and original.finite_estimates(j['kind'],cap['estimates'],j['estimates'])
        count = sum(p.suffix in ('.delete','.part_delete') for p in files)
        assert count==(2 if j['kind']=='h2' else 3*j['estimates'])
        aggregate.append({'job_id':j['job_id'],'stage':j['stage'],'kind':j['kind'],'estimates':j['estimates'],
                          'receipt_path':str(rp),'receipt_sha256':A.SEEN[str(rp)]['sha256'],'capture_sha256':r['output_sha256'],
                          'current_command_input_dependency_output_bindings_pass':True,'delete_arrays_200_bound':count,
                          'numerical_proof':'SEALED_WHOLE_CORE_CENTERED_ADJUDICATION' if j['stage']=='core' else 'SEALED_WHOLE_EXTENSION' if j['stage']=='extension' else 'CURRENT_WHOLE_VALIDATION'})
    assert len(aggregate)==190 and {stage:sum(r['stage']==stage for r in aggregate) for stage in ('core','extension','validation')}=={'core':57,'extension':112,'validation':21}
    actual = [p for p in (support/'native').rglob('*.execution_receipt.json') if not p.name.startswith('._')]
    assert len(actual)==190 and set(map(str,actual))=={r['receipt_path'] for r in aggregate}
    recorded = [x for x in monitor['native_receipt_paths'] if not Path(x).name.startswith('._')]
    monitor_sidecars = [x for x in monitor['native_receipt_paths'] if Path(x).name.startswith('._')]
    assert len(recorded)==190 and set(recorded)==set(map(str,actual))
    receipts = [p for p in native.glob('*.execution_receipt.json') if not p.name.startswith('._')]
    sidecars = [p for p in native.glob('*.execution_receipt.json') if p.name.startswith('._')]
    assert len(receipts)==21 and set(receipts)=={Path(str(native/j['job_id'])+'.execution_receipt.json') for j in jobs}
    frozen_path = ROOT/'discovery_extension/results/replication/replication_rg.tsv'
    frozen_h2_path = ROOT/'discovery_extension/results/replication/replication_source_h2.tsv'
    frozen = tab(frozen_path)
    old_idx = {(r['sleep_trait'],r['replication_source_id']):r for r in frozen}
    old_h2 = {r['replication_source_id']:r for r in tab(frozen_h2_path)}
    rg_table_path = P/'tables/validation_native_full_precision_rg_v4.tsv'
    h2_table_path = P/'tables/validation_native_full_precision_h2_v4.tsv'
    new_pairs = tab(rg_table_path)
    new_idx = {(r['sleep_trait'],r['outcome_trait']):r for r in new_pairs}
    new_h2_rows = tab(h2_table_path)
    h2_idx = {r['trait_id']:r for r in new_h2_rows}
    verification = tab(P/'tables/validation_native_job_verification_v4.tsv')
    verify_idx = {r['job_id']:r for r in verification}
    arithmetic = tab(P/'statistical_validation/validation_native_block_arithmetic_v4.tsv')
    ar_idx = {(r['kind'],r['identity']):r for r in arithmetic}
    assert len(frozen)==len(old_idx)==len(new_pairs)==len(new_idx)==41 and set(old_idx)==set(new_idx)
    assert len(new_h2_rows)==len(h2_idx)==len(old_h2)==13 and set(old_h2)==set(h2_idx)
    assert len(verification)==len(verify_idx)==21 and len(arithmetic)==len(ar_idx)==54
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
        r = A.read_json(rp)
        assert r['job']==job and r['returncode']==0 and r['scientific_cardinality_gate_pass'] is True
        assert r['execution_identity_gate_pass'] is True and r['original_outputs_modified'] is False
        assert r['command']==A.expected_command(job,prefix,original,support)
        assert r['input_sha256']==r['input_sha256_after']=={path:admitted[path] for path in job['inputs']}
        assert r['dependency_sha256_before']==r['dependency_sha256_after']==dependencies
        files = exact_owned_outputs(job,prefix)
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
                         'current_validation_receipt':True,'delete_arrays_200':count})
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
                checks['frozen_exact_printed_'+field] = float(old[field])==float(o[token])
            h2_rows.append({'trait_id':ident,'observed_h2':e['tot'],'independent_centered_delete_se':se,'captured_se':e['tot_se'],
                            'reported_scale':n['scale'],'liability_factor':factor,'checks_count':len(checks),'all_pass':all(checks.values()),
                            'failed_checks':','.join(k for k,v in checks.items() if not v)})
            continue
        assert args['rg'].split(',')==job['inputs'] and args['h2'] is None
        assert [(e['p1'],e['p2']) for e in cap['estimates']]==[(job['inputs'][0],p2) for p2 in job['inputs'][1:]]
        nt = U.correlation_tokens(Path(str(prefix)+'.log'))
        assert len(nt)==job['estimates']
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
            frozen_map = {'rg':'rg','se':'se','z':'z','p':'scalar_P','replication_h2_observed':'h2_obs','replication_h2_observed_se':'h2_obs_se',
                          'replication_h2_intercept':'h2_int','replication_h2_intercept_se':'h2_int_se','cross_trait_LDSC_intercept':'gcov_int','cross_trait_LDSC_intercept_se':'gcov_int_se'}
            for field,token in frozen_map.items():
                checks['frozen_exact_summary_'+field] = A.near(float(old[field]),math.erfc(abs(float(o['z']))/math.sqrt(2)),True) if field=='p' else float(old[field])==float(o[token])
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
    assert len(job_rows)==21 and len(h2_rows)==13 and len(pair_rows)==41 and array_count==149
    formatter = 'import json,sys,numpy as np;np.set_printoptions(linewidth=1000,precision=4);x=json.load(sys.stdin);print(json.dumps(dict(numpy=np.__version__,values=[str(np.matrix(v)).replace("[", "").replace("]", "").strip() for v in x])))'
    env = {**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','TMPDIR':str(SSD/'tmp')}
    done = subprocess.run([str(original.PYTHON),'-B','-c',formatter],input=json.dumps(values),capture_output=True,text=True,check=True,timeout=30,env=env)
    rendered = json.loads(done.stdout)
    assert rendered['numpy']=='1.21.5' and len(rendered['values'])==len(formats)
    for row,token in zip(formats,rendered['values']):
        row.update(actual_token=token,exact_match=token==row['expected_token'])
    alpha = .05/217
    manifest_path = ROOT/'discovery_extension/config/replication_manifest.tsv'
    manifest = tab(manifest_path)
    lock = A.read_json(ROOT/'discovery_extension/config/replication_manifest.lock.json')
    assert register(manifest_path)['sha256']==lock['manifest_sha256']
    assert lock['pair_count']==len(manifest)==217 and [r['pair_id'] for r in manifest]==lock['pair_ids_in_locked_order']
    assert format(alpha,'.12g')==format(lock['bonferroni_alpha'],'.12g')
    original217 = tab(ROOT/'sleep_unified_research_v1/REPLICATION_RESULTS.tsv')
    cohort = {r['pair_id']:r for r in tab(ROOT/'sleep_unified_research_v1/reviews/cohort_independence_v1.tsv')}
    original_results = {r['pair_id']:r for r in tab(ROOT/'discovery_extension/results/replication/replication_results.tsv')}
    register(ROOT/'discovery_extension/scripts/49_collate_replication_rg.py')
    register(ROOT/'discovery_extension/scripts/48_collate_replication_h2.py')
    register(ROOT/'discovery_extension/scripts/23_collate_replication.py')
    assert len(original217)==len(cohort)==len(original_results)==217
    assert [r['pair_id'] for r in original217]==lock['pair_ids_in_locked_order']
    original_idx = {r['pair_id']:r for r in original217}
    calculated = {(r['sleep_trait'],r['outcome_trait']):r for r in pair_rows}
    unavailable = set(lock['unavailable_pair_ids_in_locked_order'])
    testable = set(lock['testable_pair_ids_in_locked_order'])
    assert len(unavailable)==159 and len(testable)==58 and unavailable|testable==set(original_idx) and not unavailable&testable
    classification = []
    threshold_changes,direction_changes,h2_gate_changes = [],[],[]
    recomputed_historical_P_errors = []
    scalar_historical_P_differences = []
    for source_id,old in old_h2.items():
        native_h = h2_idx[source_id]
        oldpass = old['primary_status']=='PASS' and float(old['h2_z'])>=4 and float(old['LDSC_intercept'])<=1.2
        newpass = float(native_h['h2_z_full_precision'])>=4 and float(native_h['intercept_full_precision'])<=1.2
        if oldpass!=newpass:h2_gate_changes.append(source_id)
    for source in manifest:
        pair_id = source['pair_id']
        old,review,legacy = original_idx[pair_id],cohort[pair_id],original_results[pair_id]
        assert all(old[k]==v for k,v in legacy.items())
        assert old['current_class']==review['reviewed_class'] and old['current_QC_reason']==review['replication_outcome_h2_qc_reason']
        assert old['sleep_input_independent'].lower()==old['both_trait_independent_replication'].lower()=='false'
        assert old['discovery_rg']==source['discovery_rg'] and old['discovery_se']==source['discovery_se']
        present = pair_id not in unavailable and old_h2[source['replication_source_id']]['primary_status']=='PASS'
        nativep,nativerg,historicalp,direction,nativepositive,historicalpositive = None,None,None,None,None,None
        if pair_id in unavailable:
            nativeclass = 'NO_ELIGIBLE_EXTERNAL_SOURCE_NOT_A_NULL_RG_TEST'
            assert old['analysis_status']=='NO_INDEPENDENT_DATASET_PRE_RESULT_CLASSIFICATION' and old['replication_p']=='NA'
        else:
            source_id = source['replication_source_id']
            h = old_h2[source_id]
            native_h = h2_idx[source_id]
            native_h2pass = float(native_h['h2_z_full_precision'])>=4 and float(native_h['intercept_full_precision'])<=1.2
            historical_h2pass = h['primary_status']=='PASS' and float(h['h2_z'])>=4 and float(h['LDSC_intercept'])<=1.2
            assert present==historical_h2pass
            if not historical_h2pass:
                nativeclass = 'QC_INELIGIBLE_NOT_A_NULL_RG_TEST'
                assert old['analysis_status']=='REPLICATION_H2_FAILED_NO_PAIR_TEST' and old['replication_p']=='NA' and old['replication_h2_pass']=='False'
            else:
                key = source['sleep_trait'],source_id
                new,fit,frozenrow = new_idx[key],calculated[key],old_idx[key]
                assert new['original_pair_id']==frozenrow['pair_id']==pair_id
                nativep,nativerg = fit['captured_p'],float(new['rg'])
                historicalp = float(old['replication_p'])
                assert historicalp==float(frozenrow['p'])
                direction = float(source['discovery_rg'])*nativerg>0
                historical_direction = float(source['discovery_rg'])*float(old['replication_rg'])>0
                assert (old['direction_concordant']=='True')==historical_direction
                nativepositive,historicalpositive = nativep<alpha,historicalp<alpha
                assert nativepositive==(nativep<lock['bonferroni_alpha']) and historicalpositive==(historicalp<lock['bonferroni_alpha'])
                if nativepositive!=historicalpositive:threshold_changes.append(pair_id)
                if direction!=historical_direction:direction_changes.append(pair_id)
                nativeclass = 'QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION' if direction and nativepositive else 'DIRECTIONALLY_CONCORDANT_BELOW_FROZEN_THRESHOLD' if direction else 'DIRECTIONALLY_DISCORDANT_NO_VALIDATION'
                assert native_h2pass
                assert new['historical_classification']==old['current_class'] and new['current_interpretation']==old['current_interpretation']
                assert (new['native_p_pass_original_217_bonferroni']=='True')==nativepositive
                assert new['independent_two_trait_replication']=='False' and new['sampling_covariance_corrected']=='False'
                assert review['Dsleep_Vsleep_intersection']=='SAME_SLEEP_GWAS_REUSED_EXECUTED'
                validation_job = next(j for j in jobs if j['kind']=='rg' and U.trait(j['inputs'][0])==key[0])
                extension_job = next(j for j in plan['jobs'] if j['stage']=='extension' and j['kind']=='rg' and U.trait(j['inputs'][0])==key[0])
                assert validation_job['inputs'][0]==extension_job['inputs'][0]
                original_log = RECOVERED/frozenrow['input_log']
                token = logs[str(original_log)][key]
                hp = math.erfc(abs(float(token['z']))/math.sqrt(2))
                assert A.near(hp,historicalp,True)
                recomputed_historical_P_errors.append(relative_error(hp,historicalp))
                scalar_historical_P_differences.append(relative_error(float(token['scalar_P']),historicalp))
                dz = (float(old['replication_rg'])-float(old['discovery_rg']))/math.sqrt(float(old['discovery_se'])**2+float(old['replication_se'])**2)
                assert A.near(dz,float(old['heterogeneity_z'])) and A.near(dz*dz,float(old['heterogeneity_Q']))
                assert A.near(math.erfc(abs(dz)/math.sqrt(2)),float(old['heterogeneity_p']),True)
        assert nativeclass==old['current_class']
        classification.append({'pair_id':pair_id,'sleep_trait':source['sleep_trait'],'replication_source_id':source['replication_source_id'],
                               'native_estimate_present':present,'historical_class':old['current_class'],'independent_current_class':nativeclass,
                               'native_rg':nativerg,'native_p':nativep,'historical_p':historicalp,'exact_alpha_0_05_over_217':alpha,
                               'direction_concordant':direction,'native_bonferroni_positive':nativepositive,'historical_bonferroni_positive':historicalpositive,
                               'threshold_boundary_changed':nativepositive!=historicalpositive if present else None,
                               'independent_two_trait_replication':False,'all_classification_checks_pass':True})
    class_counts = {c:sum(r['independent_current_class']==c for r in classification) for c in set(r['independent_current_class'] for r in classification)}
    assert sum(r['native_estimate_present'] for r in classification)==41 and not h2_gate_changes and not direction_changes and not threshold_changes
    comparison = A.read_json(P/'logs/validation_native_comparison_receipt_v4.json')
    register(P/'scripts/33_compare_native_campaign_v4_4.py',comparison['collator_sha256'])
    for path,info in comparison['sources_before'].items():
        register(path,info['sha256'])
        assert Path(path).stat().st_size==info['bytes']
    independently_pass = all(r['all_pass'] for r in h2_rows+pair_rows) and all(r['exact_match'] for r in formats)
    assert comparison['job_count']==21 and comparison['h2_count']==13 and comparison['rg_count']==41
    assert comparison['monitor_sha256']==MONITOR_SHA and comparison['plan_sha256']==PLAN_SHA
    assert comparison['arithmetic_tolerance']=={'relative':A.RTOL,'absolute_estimates':A.ATOL,'absolute_p':A.PATOL}
    assert not comparison['arithmetic_failures'] and not comparison['meaningful_printed_discrepancies']
    unchanged = all(A.sha(path)==info['sha256'] and Path(path).stat().st_size==info['bytes'] for path,info in A.SEEN.items())
    assert unchanged
    A.guard()
    tables = {key:{'path':str(destinations[key]),'sha256':write_table(destinations[key],rows)} for key,rows in [('jobs',job_rows),('h2',h2_rows),('rg',pair_rows),('classification217',classification),('aggregate190',aggregate)]}
    result = {'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'status':'PASS_WHOLE_VALIDATION_AND_FULL190_PROCESSED_INPUT_REPRODUCTION' if independently_pass else 'REVIEW_REQUIRED_VALIDATION_DISCREPANCY',
              'all_independent_checks_pass':independently_pass,'validation_job_count':21,'standalone_h2_count':13,'rg_count':41,'finite200_validation_delete_array_count':array_count,
              'monitor_terminal_success_verified':True,'monitor_sha256':MONITOR_SHA,'plan_sha256':PLAN_SHA,
              'all190_current_job_receipt_command_input_dependency_output_bindings_verified':True,
              'completed_jobs_by_stage':{'core':57,'extension':112,'validation':21},'full190_delete_arrays_bound':sum(r['delete_arrays_200_bound'] for r in aggregate),
              'aggregate_capture_estimates':sum(r['estimates'] for r in aggregate),'monitor_biological_receipts':len(recorded),'monitor_AppleDouble_receipt_sidecars_excluded':len(monitor_sidecars),
              'sealed_core_numerical_receipt_sha256':A.SEEN[str(R/'independent_whole_core_receipt_v4_2.json')]['sha256'],
              'sealed_extension_numerical_receipt_sha256':A.SEEN[str(ext_receipt_path)]['sha256'],
              'current_validation_receipt_AppleDouble_sidecars_excluded':len(sidecars),
              'historical_full_precision_available':False,'native_summary_exact_4f_checks':410,'pinned_numpy_rg_scalar_P_checks':41,
              'pinned_numpy_h2_scalar_checks':sum(r['kind']=='h2' for r in formats),'formatter_numpy_version':rendered['numpy'],
              'format_failures':[r for r in formats if not r['exact_match']],
              'numerical_or_table_failures':[r for r in h2_rows+pair_rows if not r['all_pass']],
              'original217_classification_counts':class_counts,'original217_estimated_pairs':41,'independent_two_trait_replication_count':0,
              'exact_bonferroni_alpha':alpha,'historical_lock_serialized_alpha':lock['bonferroni_alpha'],
              'native_vs_historical_P_threshold_boundary_changes':threshold_changes,'native_vs_historical_direction_changes':direction_changes,'native_vs_historical_source_h2_gate_changes':h2_gate_changes,
              'historical_P_rule':'49_collate_replication_rg.py recomputes erfc(abs(four-decimal stock summary Z)/sqrt(2)); historical P is not direct scalar log P nor unavailable unrounded P.',
              'maximum_historical_erfc_printed_Z_relative_error':max(recomputed_historical_P_errors),
              'maximum_scalar_log_P_vs_historical_reconstructed_P_relative_difference':max(scalar_historical_P_differences),
              'highest_centered_reconstructed_P_relative_error':max(r['p_relative_error'] for r in pair_rows),
              'all_h2_reported_scales_observed':all(r['reported_scale']=='observed' and r['liability_factor']==1 for r in h2_rows),
              'reference_M_5_50_total':A.REFERENCE_M,'tolerances':{'rtol':A.RTOL,'atol_statistics':A.ATOL,'atol_P':A.PATOL},
              'arithmetic_method':'Unchanged sealed reviewer centered 200-delete variance, direct covariance/sqrt(h2_1*h2_2), bias correction, erfc normal P; exact original217 threshold/direction/QC classification, no BH on41.',
              'source_identity_scope':'158 admitted input identities and all190 before/after receipt bindings; no fresh raw/munged GWAS body reads or hashes.',
              'limitations':['Printed original-log agreement does not establish unavailable archival unrounded identity.',
                             '13 standalone validation h2 values are observed-scale under supplied N; absolute variance-fraction/ascertainment validity remains separate.',
                             '41 rg use the exact same discovery sleep inputs; qualified positives are external outcome-side validation, not independent two-trait replication.',
                             '17 QC-ineligible and159 unavailable rows are not null rg tests; original217 denominator retained.',
                             'Legacy heterogeneity arithmetic assumes zero covariance and remains an uncalibrated diagnostic; shared-covariance/genomic alignment/scientific calibration gates remain separate.',
                             'All190 operational/numerical reproduction does not establish the full raw-to-estimator chain, permissions, phenotype equivalence or new biological findings.'],
              'candidate_collator_imported_or_executed':False,'GWAS_bodies_read':False,'native_estimators_or_audits_launched':0,
              'formatter_subprocess_only':True,'formatter_code':formatter,'tables':tables,'consumed_artifacts':A.SEEN,
              'all_consumed_artifact_hashes_unchanged':unchanged,
              'resource_plan':{'hash_buffer_bytes':65536,'launch_internal_free_minimum':256*1024**2,'emergency_internal_floor':A.MIN_FREE,'maximum_RSS_bytes':A.MAX_RSS,'deadline_seconds':A.DEADLINE},
              'resource_usage':{'elapsed_seconds':time.monotonic()-started,'peak_self_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
              'checker_sha256':A.sha(Path(__file__))}
    with receipt_path.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({k:result[k] for k in ['status','all_independent_checks_pass','validation_job_count','standalone_h2_count','rg_count','finite200_validation_delete_array_count','original217_classification_counts','independent_two_trait_replication_count','highest_centered_reconstructed_P_relative_error','resource_usage']}))
    if not independently_pass:
        raise SystemExit('INDEPENDENT_VALIDATION_REVIEW_REQUIRED; evidence preserved')



if __name__=='__main__':
    main()
