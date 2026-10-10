#!/usr/bin/env python3
"""Collate completed native stages; retain archival precision and scientific QC.

No LDSC imports or execution. Standard-library block arithmetic is a numerical
check of each estimator separately, never evidence of cross-fit block alignment.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import native_stage_completion_v4_3 as stage_completion

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
ARCHIVE = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')
RECOVERED = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs')
INPUTS = {}
RTOL, ATOL, P_ATOL = 1e-12, 1e-15, 1e-300


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def record(path):
    path = Path(path)
    observed = {'sha256': sha(path), 'bytes': path.stat().st_size}
    prior = INPUTS.get(str(path))
    if prior is not None and observed != prior:
        raise RuntimeError('CONSUMED_SOURCE_CHANGED_SINCE_FIRST_OBSERVATION: ' + str(path))
    if prior is None:
        INPUTS[str(path)] = observed
    return path


def load(path):
    return json.loads(record(path).read_text())


def rows(path):
    with record(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def name(path):
    return Path(path).name.removesuffix('.sumstats.gz')


def output(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise RuntimeError('VERSIONED_OUTPUT_EXISTS ' + str(path))
    if path.suffix == '.json':
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    else:
        assert data
        fields = list(dict.fromkeys(k for r in data for k in r))
        with path.open('x', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)


def summary(path):
    text = record(path).read_text()
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line == 'Summary of Genetic Correlation Results')
    header = lines[start + 1].split()
    data = []
    for i in range(start + 2, len(lines)):
        if not lines[i].strip():
            break
        m = re.match(r'^(.*?)\.sumstats\.gz\s+(.*?)\.sumstats\.gz\s+(.+)$', lines[i])
        assert m, (path, i)
        values = [m[1] + '.sumstats.gz', m[2] + '.sumstats.gz', *m[3].split()]
        assert len(values) == len(header)
        data.append(dict(zip(header, values), log_line=i + 1))
    chunks = re.split(r'^Computing rg for phenotype \d+/\d+\s*$', '\n'.join(lines[:start]), flags=re.M)[1:]
    assert len(chunks) == len(data), path
    for row, chunk in zip(data, chunks):
        p = re.search(r'^P:\s*([-\d.eE+]+)', chunk, re.M)
        assert p, path
        row['scalar_p'] = p[1]
        for field, pattern in (
            ('input_snp_count', r'^Read summary statistics for (\d+) SNPs\.'),
            ('snp_overlap_after_merge', r'^After merging with summary statistics, (\d+) SNPs remain\.'),
            ('snp_overlap_valid_alleles', r'^(\d+) SNPs with valid alleles\.'),
        ):
            m = re.search(pattern, chunk, re.M)
            assert m, (path, field)
            row[field] = int(m[1])
    return data


def h2log(path):
    text = record(path).read_text()
    hm = re.search(r'^Total (Observed|Liability) scale h2:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)', text, re.M)
    im = re.search(r'^Intercept:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)', text, re.M)
    assert hm and im, path
    out = {'scale': hm[1].lower(), 'h2': float(hm[2]), 'h2_se': float(hm[3]),
           'intercept': float(im[1]), 'intercept_se': float(im[2])}
    for key, pattern in (
        ('lambda_gc', r'^Lambda GC:\s*([-\d.eE+]+)'),
        ('mean_chisq', r'^Mean Chi\^2:\s*([-\d.eE+]+)'),
        ('input_snp_count', r'^Read summary statistics for (\d+) SNPs\.'),
        ('regression_snp_count', r'^After merging with regression SNP LD, (\d+) SNPs remain\.'),
    ):
        m = re.search(pattern, text, re.M)
        assert m, (path, key)
        out[key] = float(m[1])
    m = re.search(r'^Ratio:\s*([-\d.eE+]+) \(([-\d.eE+]+)\)', text, re.M)
    if m:
        out['ratio'], out['ratio_se'] = map(float, m.groups())
    elif 'Ratio: NA (mean chi^2 < 1)' in text:
        out['ratio'], out['ratio_se'] = 'NA_MEAN_CHISQ_LT_1', 'NOT_ESTIMATED'
    else:
        assert 'Ratio < 0 (usually indicates GC correction).' in text
        out['ratio'], out['ratio_se'] = 'LT_ZERO', 'NOT_REPORTED'
    return out


def bh(ps):
    assert all(math.isfinite(p) and 0 <= p <= 1 for p in ps)
    order = sorted(range(len(ps)), key=lambda i: (ps[i], i))
    out, running = [None] * len(ps), 1.0
    for rank in range(len(ps), 0, -1):
        i = order[rank - 1]
        running = min(running, len(ps) * ps[i] / rank)
        out[i] = running
    return out


def close(a, b, field):
    return math.isclose(a, b, rel_tol=RTOL, abs_tol=P_ATOL if field == 'p' else ATOL)


def delete(path):
    vals = [float(line) for line in record(path).read_text().splitlines()]
    assert len(vals) == 200 and all(math.isfinite(v) for v in vals), path
    return vals


def block_se(total, deletes):
    n = len(deletes)
    # Algebraically identical to Var(pseudovalues)/B, avoiding subtraction
    # of the shared large B*point term from every pseudovalue.
    center = math.fsum(deletes) / n
    return math.sqrt((n - 1) / n * math.fsum((x - center) ** 2 for x in deletes))


def check_rg(prefix, e):
    stem = str(prefix) + Path(e['p1']).name + '_' + Path(e['p2']).name
    d = {key: delete(stem + '.' + key + '.delete') for key in ('hsq1', 'hsq2', 'gencov')}
    checks = {}
    for key in d:
        v = block_se(e[key]['tot'], d[key])
        checks[key + '_se'] = close(v, e[key]['tot_se'], key)
    ratio = e['gencov']['tot'] / math.sqrt(e['hsq1']['tot'] * e['hsq2']['tot'])
    checks['rg_ratio'] = close(ratio, e['rg_ratio'], 'rg_ratio')
    if all(a > 0 and b > 0 for a, b in zip(d['hsq1'], d['hsq2'])):
        ratios = [c / math.sqrt(a * b) for a, b, c in zip(d['hsq1'], d['hsq2'], d['gencov'])]
        pseudo = [200 * ratio - 199 * x for x in ratios]
        jack, se = math.fsum(pseudo) / 200, block_se(ratio, ratios)
        checks['rg_jknife'] = close(jack, e['rg_jknife'], 'rg_jknife')
        checks['rg_se'] = close(se, e['rg_se'], 'rg_se')
        z, p = ratio / se, math.erfc(abs(ratio / se) / math.sqrt(2))
        checks['z'], checks['p'] = close(z, e['z'], 'z'), close(p, e['p'], 'p')
        status = 'INDEPENDENT_200_BLOCK_RATIO_ARITHMETIC_PASS' if all(checks.values()) else 'ARITHMETIC_DISCREPANCY'
    else:
        status = 'NONPOSITIVE_DELETE_DENOMINATOR_REQUIRES_REVIEW'
    return {'checks': checks, 'status': status, 'cross_estimator_boundaries_certified': False}


def archived_log(stage, kind, ident):
    if stage == 'core':
        return ARCHIVE / 'results/logs' / (kind + '_' + ident + '.log')
    if stage == 'extension':
        return RECOVERED / 'discovery_extension/logs' / kind / (kind + '_' + ident + '.log')
    return RECOVERED / 'discovery_extension/logs/replication' / kind / (kind + ('_replication_' if kind == 'rg' else '_') + ident + '.log')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['core', 'extension', 'validation'], required=True)
    args = parser.parse_args()
    stage = args.stage
    plan = load(PLAN)
    record(Path(stage_completion.__file__))
    monitor = stage_completion.stage_monitor(stage, record=record)
    assert monitor['returncode'] == 0 and monitor['stop_reason'] is None
    assert monitor['plan_sha256'] == sha(PLAN)
    assert monitor['process_group_teardown']['remaining_group_members'] == []
    jobs = [j for j in plan['jobs'] if j['stage'] == stage]
    native = Path(plan['ssd_support_package']) / 'native' / (stage + '_reproduction_v1')
    receipt_paths = [p for p in native.glob('*.execution_receipt.json') if not p.name.startswith('._')]
    assert len(receipt_paths) == len(jobs)
    if stage == 'core':
        frozen_path = ROOT / 'sleep_unified_research_v1/tables/original_core_396.tsv'
        frozen = rows(frozen_path)
        index = {(r['sleep_trait'], r['disease_trait']): r for r in frozen}
        frozen_h2_path = ROOT / 'sleep_unified_research_v1/sources/recovered/results/tables/h2_summary.tsv'
        hindex = {r['trait']: r for r in rows(frozen_h2_path)}
    elif stage == 'extension':
        frozen_path = ROOT / 'sleep_unified_research_v1/tables/original_extension_1200.tsv'
        frozen = rows(frozen_path)
        index = {(r['sleep_trait'], r['extension_trait_id']): r for r in frozen}
        frozen_h2_path = ROOT / 'discovery_extension/results/ldsc/extension_trait_readiness.tsv'
        hindex = {r['extension_trait_id']: r for r in rows(frozen_h2_path)}
    else:
        frozen_path = ROOT / 'discovery_extension/results/replication/replication_rg.tsv'
        frozen = rows(frozen_path)
        index = {(r['sleep_trait'], r['replication_source_id']): r for r in frozen}
        frozen_h2_path = ROOT / 'discovery_extension/results/replication/replication_source_h2.tsv'
        hindex = {r['replication_source_id']: r for r in rows(frozen_h2_path)}
    assert len(index) == len(frozen)
    pairs, heritabilities, arithmetic, receipts, differences = [], [], [], [], []
    for job in jobs:
        prefix = native / job['job_id']
        rp = Path(str(prefix) + '.execution_receipt.json')
        receipt = load(rp)
        assert receipt['job'] == job and receipt['returncode'] == 0
        assert receipt['scientific_cardinality_gate_pass'] and receipt['execution_identity_gate_pass']
        assert receipt['input_sha256'] == receipt['input_sha256_after']
        assert receipt['dependency_sha256_before'] == receipt['dependency_sha256_after']
        assert receipt['all_output_sha256']
        for path, expected in receipt['all_output_sha256'].items():
            assert sha(record(path)) == expected, path
        full_path = Path(str(prefix) + '.full_precision.json')
        capture = load(full_path)
        assert sha(full_path) == receipt['output_sha256']
        assert len(capture['estimates']) == job['estimates']
        receipts.append({'job_id': job['job_id'], 'kind': job['kind'], 'estimates': job['estimates'],
                         'receipt_path': str(rp), 'receipt_sha256': sha(rp), 'full_precision_sha256': sha(full_path),
                         'all_output_hashes_verified': True, 'input_identity_before_after': True,
                         'execution_identity_gate_pass': True, 'elapsed_seconds': receipt['elapsed_seconds']})
        if job['kind'] == 'h2':
            ident = name(job['inputs'][0])
            old = hindex[ident]
            oldpath = archived_log(stage, 'h2', ident)
            oldlog, newlog = h2log(oldpath), h2log(Path(str(prefix) + '.log'))
            concordant = oldlog == newlog
            e = capture['estimates'][0]
            P, K = capture['arguments']['samp_prev'], capture['arguments']['pop_prev']
            factor = 1.0
            if P is not None and K is not None:
                P, K = float(P), float(K)
                threshold = -statistics.NormalDist().inv_cdf(K)
                density = math.exp(-threshold * threshold / 2) / math.sqrt(2 * math.pi)
                factor = K * K * (1 - K) ** 2 / (P * (1 - P) * density * density)
            row = {'trait_id': ident, 'stage': stage, 'h2_observed_full_precision': e['tot'],
                   'h2_observed_se_full_precision': e['tot_se'], 'reported_h2_full_precision': e['tot'] * factor,
                   'reported_h2_se_full_precision': e['tot_se'] * factor, 'scale': newlog['scale'],
                   'observed_fit_interpretation': 'Fitted LDSC total under supplied N conventions; population/sample variance-fraction interpretation not independently certified.',
                   'sample_prevalence_argument': P, 'population_prevalence_argument': K, 'liability_factor': factor,
                   'h2_z_full_precision': e['tot'] / e['tot_se'], 'intercept_full_precision': e['intercept'],
                   'intercept_se_full_precision': e['intercept_se'], 'lambda_gc_full_precision': e['lambda_gc'],
                   'mean_chisq_full_precision': e['mean_chisq'], 'ratio_full_precision': e['ratio'],
                   'ratio_se_full_precision': e['ratio_se'], 'input_snp_count': int(newlog['input_snp_count']),
                   'regression_snp_count': int(newlog['regression_snp_count']),
                   'h2_z_ge_4_full_precision': e['tot'] / e['tot_se'] >= 4,
                   'intercept_le_1_2_full_precision': e['intercept'] <= 1.2,
                   'physical_liability_point_in_0_1': 0 <= e['tot'] * factor <= 1 if newlog['scale'] == 'liability' else 'NOT_APPLICABLE',
                   'reproduction_status': 'PRINTED_PRECISION_CONCORDANT' if concordant else 'NUMERICALLY_DIFFERENT',
                   'original_log': str(oldpath), 'original_log_sha256': sha(oldpath),
                   'native_full_precision_path': str(full_path), 'native_full_precision_sha256': sha(full_path),
                   'original_table': str(frozen_h2_path), 'historical_full_precision_available': False}
            row.update({'original_' + k: v for k, v in old.items()})
            heritabilities.append(row)
            se = block_se(e['tot'], delete(str(prefix) + '.delete'))
            ar = {'identity': ident, 'kind': 'h2', 'independent_block_se': se, 'captured_se': e['tot_se'],
                  'status': 'INDEPENDENT_200_BLOCK_H2_ARITHMETIC_PASS' if close(se, e['tot_se'], 'se') else 'ARITHMETIC_DISCREPANCY'}
            arithmetic.append(ar)
            if not concordant:
                differences.append({'kind': 'h2', 'identity': ident, 'original_printed': oldlog, 'native_printed': newlog})
            continue
        sleep = name(job['inputs'][0])
        # Core archival primary and sensitivity rows have separate stock logs.
        # Respect each frozen row's recorded input_log instead of assuming one batch.
        historical = {}
        for key, old in index.items():
            if key[0] != sleep:
                continue
            oldpath = (ARCHIVE if stage == 'core' else RECOVERED) / old['input_log']
            if oldpath not in historical:
                historical[oldpath] = {(name(r['p1']), name(r['p2'])): r for r in summary(oldpath)}
        newlog = {(name(r['p1']), name(r['p2'])): r for r in summary(Path(str(prefix) + '.log'))}
        assert len(newlog) == job['estimates']
        for e in capture['estimates']:
            key = (name(e['p1']), name(e['p2']))
            old = index[key]
            oldpath = (ARCHIVE if stage == 'core' else RECOVERED) / old['input_log']
            a, b = historical[oldpath][key], newlog[key]
            fields = [k for k in a if k not in {'p1', 'p2', 'log_line'}]
            comparison = {k: float(a[k]) == float(b[k]) for k in fields}
            concordant = all(comparison.values())
            # Pinned ldsc.py uses pandas display.float_format='{:.4f}'.
            # Bind JSON capture to the native summary, including its separate
            # fixed-width P representation (which can legitimately print zero).
            captured_summary = {'rg': e['rg_ratio'], 'se': e['rg_se'], 'z': e['z'], 'p': e['p'],
                                'h2_obs': e['hsq2']['tot'], 'h2_obs_se': e['hsq2']['tot_se'],
                                'h2_int': e['hsq2']['intercept'], 'h2_int_se': e['hsq2']['intercept_se'],
                                'gcov_int': e['gencov']['intercept'], 'gcov_int_se': e['gencov']['intercept_se']}
            display_checks = {field: float(f'{value:.4f}') == float(b[field])
                              for field, value in captured_summary.items()}
            row = {'stage': stage, 'sleep_trait': key[0], 'outcome_trait': key[1], 'rg': e['rg_ratio'],
                   'rg_jknife_bias_corrected': e['rg_jknife'], 'se': e['rg_se'], 'z': e['z'], 'p': e['p'],
                   'ci_lower_95': e['rg_ratio'] - 1.959963984540054 * e['rg_se'],
                   'ci_upper_95': e['rg_ratio'] + 1.959963984540054 * e['rg_se'],
                   'reproduction_status': 'PRINTED_PRECISION_CONCORDANT' if concordant else 'NUMERICALLY_DIFFERENT',
                   'printed_fields_compared': len(fields), 'printed_fields_concordant': sum(comparison.values()),
                   'captured_full_precision_to_native_summary_pass': all(display_checks.values()),
                   'input_snp_count': b['input_snp_count'], 'snp_overlap_after_merge': b['snp_overlap_after_merge'],
                   'snp_overlap_valid_alleles': b['snp_overlap_valid_alleles'],
                   'original_log': str(oldpath), 'original_log_sha256': sha(oldpath), 'original_log_line': a['log_line'],
                   'native_full_precision_path': str(full_path), 'native_full_precision_sha256': sha(full_path),
                   'original_table': str(frozen_path), 'historical_full_precision_available': False,
                   'pairwise_h2_z_ge_4_both': all(e[k]['tot'] / e[k]['tot_se'] >= 4 for k in ('hsq1', 'hsq2')),
                   'pairwise_h2_intercept_le_1_2_both': all(e[k]['intercept'] <= 1.2 for k in ('hsq1', 'hsq2'))}
            for component in ('hsq1', 'hsq2', 'gencov'):
                row.update({component + '_' + k: v for k, v in e[component].items()})
            row.update({'original_' + k: v for k, v in old.items()})
            for field, quantity in (('rg', 'rg_ratio'), ('se', 'rg_se'), ('z', 'z'), ('p', 'p')):
                row[quantity + '_native_minus_original_serialized'] = e[quantity] - float(old[field])
            pairs.append(row)
            ar = check_rg(prefix, e)
            if not all(display_checks.values()):
                ar['status'] = 'CAPTURE_TO_NATIVE_DISPLAY_DISCREPANCY'
            arithmetic.append({'identity': '__'.join(key), 'kind': 'rg', 'status': ar['status'],
                               **ar['checks'], 'capture_to_native_summary_pass': all(display_checks.values()),
                               'cross_estimator_boundaries_certified': False})
            if not concordant:
                differences.append({'kind': 'rg', 'identity': key, 'nonconcordant_fields': [k for k in comparison if not comparison[k]],
                                    'original_printed': a, 'native_printed': b})
    expected = {'core': (45, 396, 57), 'extension': (100, 1200, 112), 'validation': (13, 41, 21)}[stage]
    assert (len(heritabilities), len(pairs), len(receipts)) == expected
    assert len({(r['sleep_trait'], r['outcome_trait']) for r in pairs}) == len(pairs)
    if stage in ('core', 'extension'):
        for row, q in zip(pairs, bh([r['p'] for r in pairs])):
            row['native_fdr_original_complete_family'] = q
            row['native_fdr_pass_0_05'] = q < .05
            field = 'fdr' if stage == 'core' else 'extension_fdr'
            row['frozen_fdr_original_complete_family'] = float(row['original_' + field])
            row['frozen_fdr_pass_0_05'] = float(row['original_' + field]) < .05
            row['fdr_boundary_status_changed'] = row['native_fdr_pass_0_05'] != row['frozen_fdr_pass_0_05']
    else:
        classified = rows(ROOT / 'sleep_unified_research_v1/REPLICATION_RESULTS.tsv')
        by_pair = {r['pair_id']: r for r in classified}
        for row in pairs:
            original = by_pair[row['original_pair_id']]
            row['historical_classification'] = original['current_class']
            row['native_p_pass_original_217_bonferroni'] = row['p'] < .05 / 217
            row['independent_two_trait_replication'] = False
            row['current_interpretation'] = original['current_interpretation']
            row['sampling_covariance_corrected'] = False
    arithmetic_failures = [r for r in arithmetic if r['status'] not in (
        'INDEPENDENT_200_BLOCK_RATIO_ARITHMETIC_PASS', 'INDEPENDENT_200_BLOCK_H2_ARITHMETIC_PASS')]
    result = {'recorded_utc': datetime.now(timezone.utc).isoformat(), 'stage': stage,
              'scope': 'New full-precision processed-input native historical reproduction; raw-to-estimator replay remains separately assessed.',
              'monitor_sha256': sha(stage_completion.stage_monitor_path(stage)),
              'plan_sha256': sha(PLAN), 'job_count': len(jobs), 'h2_count': len(heritabilities), 'rg_count': len(pairs),
              'collator_sha256': sha(Path(__file__)),
              'block_variance_arithmetic': 'Centered deletion variance (B-1)/B*sum((delete-mean)^2), equivalent to pseudovalue sample variance/B; initial core pseudo-arithmetic precision flag preserved separately.',
              'reproduction_status_counts_h2': dict(Counter(r['reproduction_status'] for r in heritabilities)),
              'reproduction_status_counts_rg': dict(Counter(r['reproduction_status'] for r in pairs)),
              'arithmetic_status_counts': dict(Counter(r['status'] for r in arithmetic)),
              'arithmetic_tolerance': {'relative': RTOL, 'absolute_estimates': ATOL, 'absolute_p': P_ATOL},
              'printed_comparison': 'Exact numeric equality of the literal stock-log representations; archival full precision unavailable. Frozen serialized values preserved separately.',
              'liability_conversion': 'Raw observed fit retained; pinned h2_obs_to_liab formula using exact logged P,K arguments, implemented independently with NormalDist.',
              'absolute_h2_scope': 'Reported-scale arithmetic reproduces supplied conventions; effective-N/sample-fraction consistency, ascertainment and population prevalence remain separately reviewed.',
              'historical_tables_modified': False, 'cross_estimator_boundaries_certified': False,
              'meaningful_printed_discrepancies': differences, 'arithmetic_failures': arithmetic_failures,
              'native_fdr_positives': sum(r.get('native_fdr_pass_0_05', False) for r in pairs),
              'frozen_fdr_positives': sum(r.get('frozen_fdr_pass_0_05', False) for r in pairs),
              'fdr_boundary_status_changes': sum(r.get('fdr_boundary_status_changed', False) for r in pairs),
              'sources_before': INPUTS}
    stage_completion.stage_monitor(stage, record=record)
    assert all(sha(path) == meta['sha256'] for path, meta in INPUTS.items()), 'CONSUMED_SOURCE_CHANGED'
    result['all_consumed_hashes_unchanged_after'] = True
    paths = [PACKAGE / 'tables' / (stage + '_native_full_precision_rg_v4.tsv'),
             PACKAGE / 'tables' / (stage + '_native_full_precision_h2_v4.tsv'),
             PACKAGE / 'statistical_validation' / (stage + '_native_block_arithmetic_v4.tsv'),
             PACKAGE / 'tables' / (stage + '_native_job_verification_v4.tsv'),
             PACKAGE / 'logs' / (stage + '_native_comparison_receipt_v4.json')]
    assert not any(p.exists() for p in paths), 'VERSIONED_STAGE_OUTPUT_EXISTS'
    for p, data in zip(paths, [pairs, heritabilities, arithmetic, receipts, result]):
        output(p, data)
    print(json.dumps({k: v for k, v in result.items() if k not in {'sources_before', 'meaningful_printed_discrepancies', 'arithmetic_failures'}}, indent=2))
    if differences or arithmetic_failures or result['fdr_boundary_status_changes']:
        raise SystemExit('NUMERICAL_REVIEW_REQUIRED; all results and discrepancies preserved')


if __name__ == '__main__':
    main()
