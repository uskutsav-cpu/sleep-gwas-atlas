#!/usr/bin/env python3
"""Independent completed-core verification, reading small outputs only.

No candidate-collator or LDSC import. Own prior centered-delete utility reused.
Native NumPy is invoked only to format scalar numbers, never to estimate a fit.
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
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
SPEC = importlib.util.spec_from_file_location('independent_prior_utility', Path(__file__).with_name('independent_core_subset_checker_v4_3.py'))
U = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(U)
LIMIT = 128 * 1024**2
FLOOR = 128 * 1024**2
BEGAN = time.monotonic()


def guard():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        peak *= 1024
    assert peak < LIMIT, ('REVIEWER_RSS_CAP', peak)
    assert shutil.disk_usage(ROOT).free >= FLOOR, 'REVIEWER_INTERNAL_EMERGENCY_FLOOR'
    assert time.monotonic() - BEGAN < 600, 'REVIEWER_TIME_CAP'
    return peak


def tab(path):
    return list(csv.DictReader(U.text(path).splitlines(), delimiter='\t'))


def extra_h2_tokens(path):
    d = U.h2_tokens(path)
    body = U.text(path)
    for k, pat in [('lambda_gc', r'^Lambda GC:\s*(\S+)'), ('mean_chisq', r'^Mean Chi\^2:\s*(\S+)'),
                   ('input_snps', r'^Read summary statistics for (\d+) SNPs\.'),
                   ('regression_snps', r'^After merging with regression SNP LD, (\d+) SNPs remain\.')]:
        d[k] = re.search(pat, body, re.M)[1]
    m = re.search(r'^Ratio:\s*([\d.eE+\-]+)\s+\(([^)]+)\)', body, re.M)
    if m:
        d.update(ratio=m[1], ratio_se=m[2])
    return d


def independent_bh(ps):
    # Sorted rank products, followed by an explicit suffix minimum.
    assert all(math.isfinite(p) and 0 <= p <= 1 for p in ps)
    ordered = sorted(enumerate(ps), key=lambda item: (item[1], item[0]))
    scaled = [min(1., p * len(ps) / (i+1)) for i, (_, p) in enumerate(ordered)]
    q = [None] * len(ps)
    for i, (identity, _) in enumerate(ordered):
        q[identity] = min(scaled[i:])
    return q


def relerr(a, b):
    return abs(a-b)/abs(b) if b else abs(a-b)


def main():
    assert shutil.disk_usage(ROOT).free >= 256 * 1024**2, 'REVIEWER_LAUNCH_FLOOR'
    U.text(Path(__file__))
    U.text(Path(__file__).with_name('independent_core_subset_checker_v4_3.py'))
    U.text(P/'reviews/independent_whole_core_plan_v4.json')
    U.text(ROOT/'scripts/05_collate.py')
    U.text(ROOT/'scripts/102_analyze_phase1_atlas.py')
    U.text(P/'scripts/33_compare_native_campaign_original_pseudovariance_v1.py')
    plan = U.source_json(P/'manifests/ssd_native_execution_plan_v4_3.json')
    monitor = U.source_json(P/'logs/core_native_monitor_receipt_v4.json')
    original_comparison = U.source_json(P/'logs/core_native_comparison_receipt_v4.json')
    assert monitor['returncode'] == 0 and monitor['stop_reason'] is None
    assert monitor['process_group_teardown']['remaining_group_members'] == []
    assert monitor['plan_sha256'] == U.sha(P/'manifests/ssd_native_execution_plan_v4_3.json')
    assert original_comparison['monitor_sha256'] == U.sha(P/'logs/core_native_monitor_receipt_v4.json')
    assert len(original_comparison['arithmetic_failures']) == 1
    assert original_comparison['arithmetic_failures'][0]['identity'] == 'snoring__bmi'
    native = Path(plan['ssd_support_package'])/'native/core_reproduction_v1'
    jobs = [j for j in plan['jobs'] if j['stage'] == 'core']
    receipts = [p for p in native.glob('*.execution_receipt.json') if not p.name.startswith('._')]
    sidecars = [p for p in native.glob('*.execution_receipt.json') if p.name.startswith('._')]
    assert len(jobs) == len(receipts) == 57
    assert set(receipts) == {Path(str(native/j['job_id'])+'.execution_receipt.json') for j in jobs}
    monitor_biological = [Path(s) for s in monitor['native_receipt_paths'] if not Path(s).name.startswith('._')]
    assert set(monitor_biological) == set(receipts) and len(monitor_biological) == 57
    input_idx = {r['path']:r for r in plan['inputs_verified']}
    for s, expected in {**plan['dependencies_sha256'], **plan['support_file_sha256']}.items():
        q = Path(s)
        U.SEEN[s] = {'sha256': U.sha(q), 'bytes': q.stat().st_size}
        assert U.SEEN[s]['sha256'] == expected
    U.SEEN[monitor['stdout_path']] = {'sha256':U.sha(Path(monitor['stdout_path'])), 'bytes':Path(monitor['stdout_path']).stat().st_size}
    assert U.SEEN[monitor['stdout_path']]['sha256'] == monitor['stdout_sha256']
    frozen = tab(ROOT/'sleep_unified_research_v1/tables/original_core_396.tsv')
    oldidx = {(r['sleep_trait'],r['disease_trait']):r for r in frozen}
    native_pairs = tab(P/'tables/core_native_full_precision_rg_v4.tsv')
    newidx = {(r['sleep_trait'],r['outcome_trait']):r for r in native_pairs}
    native_h2 = tab(P/'tables/core_native_full_precision_h2_v4.tsv')
    hidx = {r['trait_id']:r for r in native_h2}
    tab(P/'statistical_validation/core_native_block_arithmetic_v4.tsv')
    tab(P/'tables/core_native_job_verification_v4.tsv')
    assert len(oldidx) == len(newidx) == 396 and len(hidx) == 45
    assert set(oldidx) == set(newidx)
    oldlogs, pair_rows, h2_rows, format_values, format_checks, job_records = {}, [], [], [], [], []
    fields = {'rg':('rg_ratio',None), 'se':('rg_se',None), 'z':('z',None), 'p':('p',None),
              'h2_obs':('tot','hsq2'), 'h2_obs_se':('tot_se','hsq2'), 'h2_int':('intercept','hsq2'),
              'h2_int_se':('intercept_se','hsq2'), 'gcov_int':('intercept','gencov'), 'gcov_int_se':('intercept_se','gencov')}
    adjudication = None
    for job in jobs:
        guard()
        prefix = native/job['job_id']
        capture = U.verify_job(job,prefix,plan)
        r = U.source_json(Path(str(prefix)+'.execution_receipt.json'))
        assert all(r['input_sha256'][s] == input_idx[s]['expected_sha256'] == input_idx[s]['actual_sha256'] for s in job['inputs'])
        job_records.append({'job_id':job['job_id'], 'kind':job['kind'], 'estimates':job['estimates'], 'all_output_hashes_verified':True,
                            'receipt_sha256':U.sha(Path(str(prefix)+'.execution_receipt.json')),
                            'capture_sha256':U.sha(Path(str(prefix)+'.full_precision.json')),
                            'input_hash_binding_matches_frozen_plan':True})
        if job['kind'] == 'h2':
            ident = U.trait(job['inputs'][0]); e = capture['estimates'][0]
            old = extra_h2_tokens(U.ARCHIVE/'results/logs'/('h2_'+ident+'.log'))
            n = extra_h2_tokens(Path(str(prefix)+'.log'))
            factor = U.scale_factor(capture['arguments']['samp_prev'],capture['arguments']['pop_prev'])
            se = U.centered_se(U.estimates(Path(str(prefix)+'.delete')))
            value, error = e['tot']*factor, e['tot_se']*factor
            checks = {'old_native_printed_numeric_equal':all(old[k] == n[k] if k in ['scale','ratio_line'] else float(old[k]) == float(n[k]) for k in old),
                      'observed_200_delete_se':U.near(se,e['tot_se']),
                      'reported_scale_matches_prev_arguments':n['scale'] == ('observed' if capture['arguments']['samp_prev'] is None else 'liability'),
                      'h2_z_scale_invariant':U.near(value/error,e['tot']/e['tot_se'])}
            for k,v in {'h2_observed_full_precision':e['tot'], 'h2_observed_se_full_precision':e['tot_se'],
                        'reported_h2_full_precision':value, 'reported_h2_se_full_precision':error, 'liability_factor':factor,
                        'h2_z_full_precision':e['tot']/e['tot_se'], 'intercept_full_precision':e['intercept'],
                        'intercept_se_full_precision':e['intercept_se'], 'lambda_gc_full_precision':e['lambda_gc'],
                        'mean_chisq_full_precision':e['mean_chisq'], 'ratio_full_precision':e['ratio'], 'ratio_se_full_precision':e['ratio_se']}.items():
                checks['native_table_'+k] = U.near(v,float(hidx[ident][k]))
            fmt = {'h2':value, 'se':error, 'intercept':e['intercept'], 'intercept_se':e['intercept_se'], 'lambda_gc':e['lambda_gc'], 'mean_chisq':e['mean_chisq']}
            if 'ratio' in n:
                fmt.update(ratio=e['ratio'],ratio_se=e['ratio_se'])
            else:
                checks['ratio_qualifier_consistent'] = ('NA' in n['ratio_line'] and e['mean_chisq'] <= 1) or ('Ratio <' in n['ratio_line'] and e['ratio'] < 0 and e['mean_chisq'] > 1)
            for k,v in fmt.items():
                format_values.append(v); format_checks.append({'identity':ident,'kind':'h2','field':k,'expected_token':n[k]})
            h2_rows.append({'trait_id':ident, 'checks':checks, 'all_pass':all(checks.values()), 'independent_centered_delete_se':se,
                            'captured_observed_h2':e['tot'], 'captured_observed_se':e['tot_se'], 'independent_liability_factor':factor,
                            'reported_h2':value, 'reported_se':error, 'reported_scale':n['scale'], 'sample_prev_argument':capture['arguments']['samp_prev'],
                            'population_prev_argument':capture['arguments']['pop_prev'], 'liability_point_above_one':n['scale'] == 'liability' and value > 1,
                            'observed_fit_variance_fraction_scientifically_certified':False})
            continue
        newtokens = U.correlation_tokens(Path(str(prefix)+'.log'))
        assert len(newtokens) == 33
        for e in capture['estimates']:
            key = U.trait(e['p1']),U.trait(e['p2']); old,new = oldidx[key],newidx[key]
            oldpath = U.ARCHIVE/old['input_log']
            if str(oldpath) not in oldlogs:
                oldlogs[str(oldpath)] = U.correlation_tokens(oldpath)
            o,n = oldlogs[str(oldpath)][key],newtokens[key]
            stem = str(prefix)+Path(e['p1']).name+'_'+Path(e['p2']).name
            dels = {c:U.estimates(Path(stem+'.'+c+'.delete')) for c in ('hsq1','hsq2','gencov')}
            checks = {'printed_old_native_numeric_equal':all(float(o[k]) == float(n[k]) for k in o),
                      'positive_200_delete_denominators':all(a > 0 and b > 0 for a,b in zip(dels['hsq1'],dels['hsq2']))}
            assert checks['positive_200_delete_denominators']
            for c,d in dels.items():
                checks[c+'_centered_delete_se'] = U.near(U.centered_se(d),e[c]['tot_se'])
            ratio = e['gencov']['tot']/math.sqrt(e['hsq1']['tot']*e['hsq2']['tot'])
            ds = [c/math.sqrt(a*b) for a,b,c in zip(dels['hsq1'],dels['hsq2'],dels['gencov'])]
            se = U.centered_se(ds)
            bias = ratio + 199*(ratio-math.fsum(ds)/200)
            z = ratio/se; p = math.erfc(abs(z)/math.sqrt(2))
            checks.update(rg=U.near(ratio,e['rg_ratio']),rg_jknife=U.near(bias,e['rg_jknife']),se=U.near(se,e['rg_se']),
                          z=U.near(z,e['z']),p=U.near(p,e['p'],True),p_from_captured_z=U.near(math.erfc(abs(e['z'])/math.sqrt(2)),e['p'],True))
            for k in ('rg','se','z','p'):
                value = e['rg_ratio'] if k == 'rg' else e['rg_se'] if k == 'se' else e[k]
                checks['native_table_'+k] = float(new[k]) == value
            checks['native_table_bias_corrected_rg'] = float(new['rg_jknife_bias_corrected']) == e['rg_jknife']
            checks['native_table_ci_lower'] = U.near(float(new['ci_lower_95']),ratio-1.959963984540054*e['rg_se'])
            checks['native_table_ci_upper'] = U.near(float(new['ci_upper_95']),ratio+1.959963984540054*e['rg_se'])
            for displayed,(field,component) in fields.items():
                v = e[component][field] if component else e[field]
                checks['exact_summary_4f_'+displayed] = format(v,'.4f') == n[displayed]
                token = o['scalar_P'] if displayed == 'p' else o[displayed]
                checks['frozen_log_declared_4g_'+displayed] = format(float(token),'.4g') == format(float(old[displayed]),'.4g')
            for component in ('hsq1','hsq2','gencov'):
                for k,v in e[component].items():
                    checks['native_table_'+component+'_'+k] = float(new[component+'_'+k]) == v
            for dest,token in [('input_snp_count','read_SNPs'),('snp_overlap_after_merge','merged_SNPs'),('snp_overlap_valid_alleles','valid_allele_SNPs')]:
                checks['native_table_'+dest] = int(new[dest]) == int(n[token])
            format_values.append(e['p']);format_checks.append({'identity':'__'.join(key),'kind':'rg','field':'scalar_P','expected_token':n['scalar_P']})
            if key == ('snoring','bmi'):
                # Diagnostic reconstruction of the first collator's exact variance route.
                pseudo = [200*ratio-199*x for x in ds]
                pseudo_se = math.sqrt(statistics.variance(pseudo)/200)
                pseudo_p = math.erfc(abs(ratio/pseudo_se)/math.sqrt(2))
                native_z_p = math.erfc(abs(e['z'])/math.sqrt(2))
                adjudication = {'pair_id':'snoring__bmi','native_se':e['rg_se'],'centered_se':se,'pseudo_se':pseudo_se,
                                'native_z':e['z'],'native_p':e['p'],'centered_reconstructed_p':p,'pseudo_reconstructed_p':pseudo_p,
                                'erfc_captured_z':native_z_p,'centered_p_relative_error':relerr(p,e['p']),
                                'pseudo_p_relative_error':relerr(pseudo_p,e['p']),'captured_z_erfc_relative_error':relerr(native_z_p,e['p']),
                                'pseudo_p_fixed_tolerance_pass':U.near(pseudo_p,e['p'],True),'centered_p_fixed_tolerance_pass':U.near(p,e['p'],True),
                                'captured_z_erfc_fixed_tolerance_pass':U.near(native_z_p,e['p'],True),'tolerance_relaxed':False,
                                'initial_receipt_preserved':True,'interpretation':'Equivalent jackknife variance routes differ by binary64 roundoff, amplified by the extreme normal tail. The initial pseudovariance failure is retained. Centered reconstruction and erfc(captured Z) pass the original fixed tolerance.'}
            pair_rows.append({'pair_id':'__'.join(key),'historical_tier':old['analysis_tier'],'all_pass':all(checks.values()),'checks':checks,
                              'independent_centered_se':se,'independent_ratio':ratio,'independent_bias_corrected_rg':bias,'independent_z':z,'independent_p':p,
                              'p_relative_error':relerr(p,e['p']),'native_p':e['p'],'original_log':str(oldpath)})
    assert len(pair_rows) == 396 and len(h2_rows) == 45
    # Pinned native NumPy scalar rendering, no LDSC import/data read/fit.
    python = next(k for k in plan['dependencies_sha256'] if k.endswith('/.ldsc-env/bin/python'))
    formatter = 'import json,sys,numpy as np;np.set_printoptions(linewidth=1000,precision=4);x=json.load(sys.stdin);print(json.dumps(dict(numpy=np.__version__,values=[str(np.matrix(v)).replace("[", "").replace("]", "").strip() for v in x])))'
    env = {**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','TMPDIR':str(Path(plan['ssd_support_package']).parents[1]/'tmp')}
    done = subprocess.run([python,'-B','-c',formatter],input=json.dumps(format_values),capture_output=True,text=True,check=True,timeout=120,env=env)
    rendered = json.loads(done.stdout)
    assert rendered['numpy'] == '1.21.5' and len(rendered['values']) == len(format_checks)
    for row,token in zip(format_checks,rendered['values']):
        row.update(actual_token=token,exact_match=token == row['expected_token'])
    qs = independent_bh([r['native_p'] for r in pair_rows])
    bh_rows = []
    for pair,q in zip(pair_rows,qs):
        a,b = pair['pair_id'].split('__'); key = a,b
        n,o = newidx[key],oldidx[key]
        bh_rows.append({'pair_id':pair['pair_id'],'tier':pair['historical_tier'],'q':q,'table_q':float(n['native_fdr_original_complete_family']),
                        'q_table_matches':U.near(q,float(n['native_fdr_original_complete_family'])),'native_positive':q < .05,
                        'frozen_positive':float(o['fdr']) < .05,'boundary_changed':(q < .05) != (float(o['fdr']) < .05)})
    peak = guard()
    hashes_unchanged = all(U.sha(Path(s)) == m['sha256'] for s,m in U.SEEN.items())
    all_pass = all(r['all_pass'] for r in pair_rows+h2_rows) and all(r['exact_match'] for r in format_checks) and all(r['q_table_matches'] and not r['boundary_changed'] for r in bh_rows) and hashes_unchanged
    result = {'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Entire completed core processed-input native campaign:57jobs,45h2,396rg; no raw GWAS rerun, extension or validation completion claim.',
              'job_count':len(jobs),'h2_count':len(h2_rows),'rg_count':len(pair_rows),'AppleDouble_receipt_sidecars_excluded':len(sidecars),
              'monitor_terminal_success_verified':True,'original_comparison_failure_count_preserved':len(original_comparison['arithmetic_failures']),
              'all_independent_checks_pass':all_pass,'fixed_tolerances':{'rtol':U.RTOL,'atol_statistics':U.ATOL,'atol_p':U.PATOL},
              'independent_method':'Centered delete variance; direct covariance/h2 ratio; erfc P; bias correction from deletion mean; inverse-normal tail bisection; suffix-min BH; pinned NumPy formatting only.',
              'exact_summary_4f_count':3960,'exact_numpy_rg_scalar_p_count':396,'exact_numpy_h2_scalar_count':sum(r['kind']=='h2' for r in format_checks),
              'native_family_bh_positives':sum(r['native_positive'] for r in bh_rows),'frozen_family_bh_positives':sum(r['frozen_positive'] for r in bh_rows),
              'native_primary_bh_positives':sum(r['native_positive'] and r['tier']=='PRIMARY_PHASE1' for r in bh_rows),
              'native_sensitivity_bh_positives':sum(r['native_positive'] and r['tier']=='QC_FAILED_SENSITIVITY' for r in bh_rows),
              'bh_boundary_changes':sum(r['boundary_changed'] for r in bh_rows),'highest_centered_reconstructed_p_relative_error':max(r['p_relative_error'] for r in pair_rows),
              'reported_liability_points_above_one':[r['trait_id'] for r in h2_rows if r['liability_point_above_one']],
              'snoring_bmi_adjudication':adjudication,'pairs':pair_rows,'h2':h2_rows,'format_checks':format_checks,'bh':bh_rows,'jobs':job_records,
              'small_output_dependency_support_hashes_before_and_after_verified':hashes_unchanged,'inputs':U.SEEN,
              'GWAS_source_bodies_read':False,'current_input_hashes_independently_rehashed':False,'input_hash_scope':'All57before/after execution receipt bindings independently match frozen plan; large processed inputs not rehashed again in this review.',
              'candidate_collator_imported_or_executed':False,'historical_full_precision_available':False,'frozen_serialization_precision':'Original collator %.4g, subsequent atlas CSV %.15g may preserve float conversion tails; never infer original unrounded identity.',
              'cross_fit_block_alignment_certified':False,'source_to_estimator_complete_raw_chain_certified':False,'sampling_covariance_corrected':False,
              'absolute_h2_scientific_validity_certified':False,'new_biological_finding':False,'resource_limits':{'reviewer_rss_bytes':LIMIT,'launch_internal_bytes':256*1024**2,'emergency_internal_bytes':FLOOR,'wall_seconds':600,'hash_buffer_bytes':65536},
              'resource_usage':{'peak_observed_reviewer_rss_bytes':peak,'elapsed_seconds':time.monotonic()-BEGAN},'formatter_numpy_version':rendered['numpy'],'formatter_code':formatter}
    dest = P/'reviews/independent_whole_core_receipt_v4.json'
    assert not dest.exists();dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pairs','h2','format_checks','bh','jobs','inputs','formatter_code']},indent=2))
    if not all_pass:
        raise SystemExit('INDEPENDENT_CORE_REVIEW_FAILED; receipt preserved')


if __name__ == '__main__':
    main()
