#!/usr/bin/env python3
"""Read-only frozen native checkpoint audit; no GWAS bodies or estimators.

Independently uses 64-KiB hashes and centered 200-delete variance.  Reads only
code, metadata, stock logs, captured estimates and small output delete arrays.
Never imports the candidate collator or LDSC; original runner import calls only
its config-reading jobs() function.  Writes new immutable review artifacts.
"""
import csv
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import time

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
REVIEW = P / 'reviews'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = P / 'manifests/ssd_native_execution_plan_v4_3.json'
EXPECTED_PLAN = 'f555dc441e8c93c528e153d2e8689edabb29ba17e4d8b0da68402ac4e4c98e88'
STEM = 'independent_partial_extension_checkpoint'
SEEN = {}
RTOL, ATOL, PATOL = 1e-12, 1e-15, 1e-300
START = time.monotonic()
MIN_FREE, MAX_RSS, DEADLINE = 128 * 1024**2, 128 * 1024**2, 600


def guard():
    assert shutil.disk_usage('/System/Volumes/Data').free >= MIN_FREE, 'AUDIT_INTERNAL_EMERGENCY_FLOOR'
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss <= MAX_RSS, 'AUDIT_RSS_LIMIT'
    assert time.monotonic() - START <= DEADLINE, 'AUDIT_DEADLINE'


def sha(path):
    guard()
    path = Path(path)
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def register(path, expected=None):
    path = Path(path)
    assert not path.name.startswith('._'), 'APPLEDOUBLE_IS_NOT_A_BIOLOGICAL_ARTIFACT'
    actual = sha(path)
    if expected is not None:
        assert actual == expected, 'HASH_MISMATCH ' + str(path)
    record = {'sha256': actual, 'bytes': path.stat().st_size}
    if str(path) in SEEN:
        assert SEEN[str(path)] == record, 'CONSUMED_ARTIFACT_CHANGED_DURING_AUDIT'
    SEEN[str(path)] = record
    return record


def read_text(path):
    register(path)
    return Path(path).read_text()


def read_json(path):
    return json.loads(read_text(path))


def near(a, b, p=False):
    return abs(a - b) <= max(RTOL * max(abs(a), abs(b)), PATOL if p else ATOL)


def centered_se(d):
    n = len(d)
    mean = math.fsum(d) / n
    ss = math.fsum((value - mean)**2 for value in d)
    return (n - 1) * math.sqrt(ss / (n * (n - 1)))


def deletion(path):
    values = [float(x) for x in read_text(path).split()]
    assert len(values) == 200 and all(math.isfinite(x) for x in values), 'INVALID_200_DELETE_ARRAY ' + str(path)
    return values


def summary_tokens(path):
    lines = read_text(path).splitlines()
    at = lines.index('Summary of Genetic Correlation Results')
    header = lines[at + 1].split()
    rows = []
    for line in lines[at + 2:]:
        if not line.strip():
            break
        parts = re.split(r'\.sumstats\.gz\s+', line, maxsplit=2)
        assert len(parts) == 3
        values = parts[2].split()
        assert len(values) == len(header) - 2
        rows.append((parts[0].strip() + '.sumstats.gz', parts[1].strip() + '.sumstats.gz', dict(zip(header[2:], values))))
    return rows


def expected_command(job, prefix, original, support):
    return [str(original.PYTHON), '-u', str(support / 'scripts/native_ldsc_capture.py'),
            '--ldsc-dir', str(ROOT.parent / 'ldsc-code'), '--' + job['kind'], ','.join(job['inputs']),
            '--ref-ld-chr', str(original.REF) + '/', '--w-ld-chr', str(original.REF) + '/',
            '--print-delete-vals', '--out', str(prefix)] + job['options']


def fit_checks(job, prefix, cap):
    rows = []
    if job['kind'] == 'h2':
        d = deletion(str(prefix) + '.delete')
        part = deletion(str(prefix) + '.part_delete')
        e = cap['estimates'][0]
        checks = {'standalone_observed_SE': near(centered_se(d), e['tot_se']),
                  'one_annotation_partition_matches_total': all(near(a, b) for a, b in zip(d, part))}
        if job['stage'] == 'extension':
            checks['observed_scale_no_prevalence_arguments'] = cap['arguments']['samp_prev'] is None and cap['arguments']['pop_prev'] is None
        rows.append({'fit_id': job['job_id'], 'checks': checks, 'all_pass': all(checks.values())})
        return rows, 2
    displayed = summary_tokens(str(prefix) + '.log')
    assert len(displayed) == job['estimates']
    fields = {'rg': ('rg_ratio', None), 'se': ('rg_se', None), 'z': ('z', None), 'p': ('p', None),
              'h2_obs': ('tot', 'hsq2'), 'h2_obs_se': ('tot_se', 'hsq2'),
              'h2_int': ('intercept', 'hsq2'), 'h2_int_se': ('intercept_se', 'hsq2'),
              'gcov_int': ('intercept', 'gencov'), 'gcov_int_se': ('intercept_se', 'gencov')}
    for e, (p1, p2, tokens) in zip(cap['estimates'], displayed):
        assert (p1, p2) == (e['p1'], e['p2'])
        stem = str(prefix) + Path(p1).name + '_' + Path(p2).name
        d = {name: deletion(stem + '.' + name + '.delete') for name in ['hsq1', 'hsq2', 'gencov']}
        checks = {name + '_total_SE': near(centered_se(values), e[name]['tot_se']) for name, values in d.items()}
        valid = e['hsq1']['tot'] > 0 and e['hsq2']['tot'] > 0 and all(a > 0 and b > 0 for a, b in zip(d['hsq1'], d['hsq2']))
        checks['positive_h2_point_and_delete_denominators'] = valid
        computed = {}
        if valid:
            rg = e['gencov']['tot'] / math.sqrt(e['hsq1']['tot'] * e['hsq2']['tot'])
            rd = [c / math.sqrt(a * b) for a, b, c in zip(d['hsq1'], d['hsq2'], d['gencov'])]
            se = centered_se(rd)
            bias = rg + 199 * (rg - math.fsum(rd) / 200)
            z = rg / se
            p = math.erfc(abs(z) / math.sqrt(2))
            computed = {'rg_ratio': rg, 'rg_jknife': bias, 'rg_se': se, 'z': z, 'p': p}
            checks.update({field: near(value, e[field], field == 'p') for field, value in computed.items()})
            checks['P_from_captured_Z'] = near(math.erfc(abs(e['z']) / math.sqrt(2)), e['p'], True)
        for display, (field, component) in fields.items():
            value = e[component][field] if component else e[field]
            checks['summary_exact_4f_' + display] = format(value, '.4f') == tokens[display]
        fit = {'fit_id': Path(p1).name[:-12] + '__' + Path(p2).name[:-12], 'checks': checks, 'all_pass': all(checks.values())}
        if not fit['all_pass']:
            fit['independent_arithmetic'] = computed
            fit['captured_arithmetic'] = {field: e[field] for field in computed}
        rows.append(fit)
    return rows, 3 * len(rows)


def main():
    outputs = [REVIEW / (STEM + '_receipt_v4.json'), REVIEW / (STEM + '_v4.tsv')]
    assert all(not p.exists() for p in outputs), 'IMMUTABLE_REVIEW_TARGET_EXISTS'
    assert shutil.disk_usage('/System/Volumes/Data').free >= 256 * 1024**2
    plan = read_json(PLAN)
    assert SEEN[str(PLAN)]['sha256'] == EXPECTED_PLAN
    original_path = ROOT / 'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py'
    register(original_path, plan['original_runner_sha256'])
    register(P / 'scripts/30_prepare_and_run_ssd_native_campaign.py', plan['launcher_sha256'])
    register(P / 'scripts/31_invoke_original_native_runner.py', plan['invoker_sha256'])
    original_plan = read_json(ROOT / 'sleep_unified_research_v1/manifests/native_reproduction_jobs_v1.json')
    assert SEEN[str(ROOT / 'sleep_unified_research_v1/manifests/native_reproduction_jobs_v1.json')]['sha256'] == plan['original_plan_sha256']
    assert original_plan['jobs'] == plan['jobs'] and len(plan['jobs']) == 190
    spec = importlib.util.spec_from_file_location('audit_original_config_reader', original_path)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    assert original.jobs() == plan['jobs']
    for config in [ROOT / 'config/analysis_panel.tsv', ROOT / 'discovery_extension/config/candidate_traits.tsv',
                   ROOT / 'discovery_extension/config/replication_manifest.tsv', ROOT / 'discovery_extension/results/replication/replication_rg_jobs.tsv']:
        register(config)
    for path, expected in plan['dependencies_sha256'].items():
        register(path, expected)
    for path, expected in plan['support_file_sha256'].items():
        register(path, expected)
    support = Path(plan['ssd_support_package'])
    expected_deps = dict(plan['dependencies_sha256'])
    old_wrapper = str(ROOT / 'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    new_wrapper = str(support / 'scripts/native_ldsc_capture.py')
    assert old_wrapper in expected_deps and new_wrapper not in expected_deps
    expected_deps[new_wrapper] = expected_deps.pop(old_wrapper)
    admitted = {r['path']: r['actual_sha256'] for r in plan['inputs_verified']}
    assert len(admitted) == 158 and all(r['match'] is True and r['expected_sha256'] == r['actual_sha256'] for r in plan['inputs_verified'])
    ledger = {}
    for name, filter_function in [('native_input_hash_checks.tsv', lambda r: r['kind'] == 'core_munged'),
                                   ('archived_extension_recovery.tsv', lambda r: r['status'] == 'EXACT_RECEIPT_HASH_RECOVERED')]:
        rows = csv.DictReader(read_text(support / 'tables' / name).splitlines(), delimiter='\t')
        ledger.update({r['path']: r['actual_sha256'] for r in rows if filter_function(r)})
    assert all(ledger[path] == digest for path, digest in admitted.items())
    monitor_path = P / 'logs/extension_native_monitor_receipt_v4.json'
    monitor = read_json(monitor_path)
    register(SSD / 'logs' / monitor_path.name, SEEN[str(monitor_path)]['sha256'])
    assert monitor['stage'] == 'extension' and monitor['plan_sha256'] == EXPECTED_PLAN
    assert monitor['returncode'] == -15 and monitor['stop_reason'] == 'INTERNAL_FULL_NATIVE_FLOOR_REACHED'
    assert monitor['process_group_teardown']['remaining_group_members'] == []
    assert monitor['final_resource_snapshot']['internal_free_bytes'] < plan['internal_floor_bytes'] == 3 * 1024**3
    register(monitor['stdout_path'], monitor['stdout_sha256'])
    runner_stdout = read_text(monitor['stdout_path'])
    core_proof = read_json(REVIEW / 'independent_whole_core_receipt_v4_2.json')
    assert SEEN[str(REVIEW / 'independent_whole_core_receipt_v4_2.json')]['sha256'] == '72732c58a0de3fd87995e748e111d3cf24cc4b157fdde443a59a19b252850fa0'
    # Core arithmetic is reused from this independent complete seal; current
    # core receipt/capture/output hashes are nevertheless checked again below.
    table, fits, seen_receipts, array_count = [], [], [], 0
    for index, job in enumerate(plan['jobs'], 1):
        guard()
        prefix = support / 'native' / (job['stage'] + '_reproduction_v1') / job['job_id']
        receipt_path = Path(str(prefix) + '.execution_receipt.json')
        files = sorted(p for p in prefix.parent.glob(prefix.name + '*') if p.is_file() and not p.name.startswith('._'))
        row = {'frozen_order': index, 'job_id': job['job_id'], 'stage': job['stage'], 'kind': job['kind'],
               'expected_estimates': job['estimates'], 'captured_estimates': 0,
               'status': 'UNSTARTED', 'receipt_path': '', 'receipt_sha256': '',
               'capture_sha256': '', 'output_files': len(files), 'delete_arrays_200': 0,
               'arithmetic_checks_pass': '', 'original_command_verified': False,
               'input_hash_metadata_verified': False, 'current_output_hashes_verified': False}
        if receipt_path.exists():
            r = read_json(receipt_path)
            assert r['job'] == job and r['returncode'] == 0 and r['scientific_cardinality_gate_pass'] is True
            assert r['execution_identity_gate_pass'] is True and r['original_outputs_modified'] is False
            assert r['command'] == expected_command(job, prefix, original, support)
            expected_inputs = {path: admitted[path] for path in job['inputs']}
            assert r['input_sha256'] == r['input_sha256_after'] == expected_inputs
            assert r['dependency_sha256_before'] == r['dependency_sha256_after'] == expected_deps
            assert set(map(str, files)) == set(r['all_output_sha256']) | {str(receipt_path)}
            for path, digest in r['all_output_sha256'].items():
                assert Path(path).parent == prefix.parent and Path(path).name.startswith(prefix.name)
                register(path, digest)
            capture_path = Path(str(prefix) + '.full_precision.json')
            register(capture_path, r['output_sha256'])
            cap = read_json(capture_path)
            args = cap['arguments']
            assert cap['ldsc_dir'] == str(ROOT.parent / 'ldsc-code')
            assert cap['python'].startswith('3.9.23') and cap['libraries'] == {'numpy': '1.21.5', 'pandas': '1.3.3', 'scipy': '1.7.3'}
            assert cap['instrumentation'] == 'return-value capture; numerical estimator unmodified'
            assert args['out'] == str(prefix) and args['n_blocks'] == 200 and args['print_delete_vals'] is True
            assert args['ref_ld_chr'] == args['w_ld_chr'] == str(original.REF) + '/'
            assert args['no_check_alleles'] is False and args['no_intercept'] is False
            assert len(cap['estimates']) == job['estimates'] and original.finite_estimates(job['kind'], cap['estimates'], job['estimates'])
            if job['kind'] == 'h2':
                assert args['h2'] == job['inputs'][0] and cap['estimates'][0]['input'] == job['inputs'][0]
                assert args['rg'] is None
            else:
                assert args['rg'].split(',') == job['inputs'] and args['h2'] is None
                assert [(e['p1'], e['p2']) for e in cap['estimates']] == [(job['inputs'][0], p2) for p2 in job['inputs'][1:]]
            native_arrays = sorted(p for p in files if p.suffix == '.delete' or p.suffix == '.part_delete')
            if job['stage'] == 'extension':
                fitted, count = fit_checks(job, prefix, cap)
                fits.extend(fitted)
                assert len(native_arrays) == count
                checks_pass = all(fit['all_pass'] for fit in fitted)
            else:
                for path in native_arrays:
                    deletion(path)
                count = len(native_arrays)
                assert count == (2 if job['kind'] == 'h2' else 3 * job['estimates'])
                checks_pass = 'PRIOR_INDEPENDENT_WHOLE_CORE_SEAL'
            array_count += count
            seen_receipts.append(str(receipt_path))
            row.update(status='COMPLETED_SEALED', captured_estimates=len(cap['estimates']),
                       receipt_path=str(receipt_path), receipt_sha256=SEEN[str(receipt_path)]['sha256'],
                       capture_sha256=r['output_sha256'], delete_arrays_200=count,
                       arithmetic_checks_pass=checks_pass, original_command_verified=True,
                       input_hash_metadata_verified=True, current_output_hashes_verified=True)
        elif files:
            assert job['job_id'] == 'extension_rg_snoring'
            row.update(status='STOPPED_PARTIAL_UNSEALED', receipt_path=str(receipt_path))
        table.append(row)
    counts = {stage: sum(row['stage'] == stage and row['status'] == 'COMPLETED_SEALED' for row in table) for stage in ['core', 'extension', 'validation']}
    assert counts == {'core': 57, 'extension': 107, 'validation': 0}
    recorded_receipts = [path for path in monitor['native_receipt_paths'] if not Path(path).name.startswith('._')]
    recorded_sidecars = [path for path in monitor['native_receipt_paths'] if Path(path).name.startswith('._')]
    assert set(recorded_receipts) == set(seen_receipts) and len(recorded_sidecars) == 164
    native = support / 'native'
    actual_receipts = sorted(str(path) for path in native.rglob('*.execution_receipt.json') if not path.name.startswith('._'))
    assert actual_receipts == sorted(seen_receipts)
    snoring = next(job for job in plan['jobs'] if job['job_id'] == 'extension_rg_snoring')
    prefix = native / 'extension_reproduction_v1/extension_rg_snoring'
    partial = sorted(p for p in prefix.parent.glob(prefix.name + '*') if p.is_file() and not p.name.startswith('._'))
    partial_files = []
    for path in partial:
        info = register(path)
        item = {'path': str(path), **info, 'native_delete_cardinality': None}
        if path.suffix == '.delete':
            item['native_delete_cardinality'] = len(deletion(path))
        partial_files.append(item)
    assert len(partial_files) == 284 and sum(p['native_delete_cardinality'] == 200 for p in partial_files) == 282
    expected_partial = {str(prefix) + '.log', str(prefix) + '.stdout.log'}
    expected_partial.update(str(prefix) + Path(snoring['inputs'][0]).name + '_' + Path(p2).name + '.' + component + '.delete'
                            for p2 in snoring['inputs'][1:95] for component in ['hsq1', 'hsq2', 'gencov'])
    assert {p['path'] for p in partial_files} == expected_partial
    partial_stdout = read_text(str(prefix) + '.stdout.log')
    assert 'Summary of Genetic Correlation Results' not in partial_stdout
    assert 'Computing rg for phenotype 96/101' in partial_stdout
    assert not Path(str(prefix) + '.full_precision.json').exists()
    assert not Path(str(prefix) + '.execution_receipt.json').exists()
    starts = re.findall(r'^NATIVE_START\s+\d+/112\s+(\S+)$', runner_stdout, re.M)
    completes = re.findall(r'^NATIVE_COMPLETE\s+(\S+)\s+elapsed_seconds=', runner_stdout, re.M)
    assert len(starts) == 108 and len(completes) == 107 and starts[-1] == 'extension_rg_snoring'
    assert set(completes) == {row['job_id'] for row in table if row['stage'] == 'extension' and row['status'] == 'COMPLETED_SEALED'}
    after = {path: sha(path) == info['sha256'] for path, info in SEEN.items()}
    assert all(after.values())
    guard()
    incomplete = [job for job, row in zip(plan['jobs'], table) if row['status'] != 'COMPLETED_SEALED']
    failures = [fit for fit in fits if not fit['all_pass']]
    result = {
        'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'status': 'PARTIAL_CHECKPOINT_PROVENANCE_VERIFIED' if not failures else 'PARTIAL_CHECKPOINT_PROVENANCE_VERIFIED_NUMERICAL_REVIEW_REQUIRED',
        'scope': 'Terminal original extension checkpoint only; completed extension within-fit arithmetic; current core outputs rehashed against previous complete independent core proof.',
        'frozen_jobs': 190, 'completed_biological_jobs': 164, 'completed_jobs_by_stage': counts,
        'completed_extension_h2': 100, 'completed_extension_rg_batches': 7, 'completed_extension_rg_pairs': 700,
        'completed_extension_delete_arrays_200': 2300, 'all_completed_delete_arrays_200': array_count,
        'partial_snoring_artifacts': partial_files, 'partial_snoring_complete_delete_triplets': 94,
        'partial_snoring_accepted_estimates': 0,
        'completed_extension_batch_names': [row['job_id'] for row in table if row['stage'] == 'extension' and row['kind'] == 'rg' and row['status'] == 'COMPLETED_SEALED'],
        'incomplete_original_jobs': incomplete,
        'numerical_fits_checked': len(fits), 'all_completed_extension_within_fit_arithmetic_and_exact_summary_display_pass': not failures,
        'numerical_failures': failures,
        'tolerances': {'rtol': RTOL, 'atol_statistics': ATOL, 'atol_P': PATOL},
        'arithmetic_method': 'Centered 200-delete variance; full point covariance / sqrt(h2_1*h2_2); bias corrected ratio; erfc two-sided normal P. No collator import. No tolerance relaxation.',
        'summary_display_method': 'Exact native .4f summary tokens for 10 fields per completed extension rg row. Scalar P display and historical extension comparison are outside this checkpoint audit.',
        'source_hash_verification_scope': '158 current-plan input identities and completed per-job before/after input SHA receipts reconciled to admitted ledgers; no fresh GWAS body read or hash in this review.',
        'core_arithmetic_proof_sha256': SEEN[str(REVIEW / 'independent_whole_core_receipt_v4_2.json')]['sha256'],
        'monitor': {key: monitor[key] for key in ['stage', 'started_utc', 'completed_utc', 'returncode', 'stop_reason', 'process_group_teardown', 'final_resource_snapshot', 'peak_observed_owned_rss_bytes', 'elapsed_seconds', 'stdout_path', 'stdout_sha256']},
        'monitor_biological_receipt_paths': len(recorded_receipts), 'monitor_AppleDouble_sidecars_excluded': len(recorded_sidecars),
        'original_runner_stdout_starts': len(starts), 'original_runner_stdout_completions': len(completes),
        'partial_snoring_stock_log_bytes': Path(str(prefix) + '.log').stat().st_size,
        'partial_snoring_stdout_last_attempted_phenotype_number': 96,
        'full_1200_extension_family_complete': False, 'all_190_prerequisites_complete': False,
        'validation_41_pair_family_complete': False, 'cross_fit_genomic_alignment_certified': False,
        'raw_to_harmonized_to_munged_chain_certified': False, 'new_biological_claim_admitted': False,
        'execution_authorized_by_this_audit': False, 'GWAS_bodies_read': False, 'native_estimators_or_audits_launched': 0,
        'continuation_requirements': [
            'Use a distinct additive operational plan/controller/monitor/stdout version bound to the failed original monitor, this audit and an immutable snoring interruption archive receipt.',
            'The unchanged original04/31 may skip107 extension successes only after exact original command, input, dependency, full capture and complete output hash checks; preserve every success byte and receipt.',
            'Freeze all284 snoring biological artifact source/destination paths, SHA and byte identities before same-volume owned rename; journal each action, retain full readback and rollback proof, never delete or overwrite. Sidecars are transport metadata, not estimates.',
            'Replay the entire snoring100-pair batch plus four unstarted original extension batches in the frozen order/commands/options;94 partial delete triplets are insufficient for completed batch admission.',
            'Retain unchanged3GiB internal/5GiB SSD floors, one estimator/one BLAS thread, resource caps, exclusive shared heavy lock, deferred catchable signals and finally-safe owned-group cleanup/quarantine.',
            'Seal continuation separately and bind both old failed monitor and new successful monitor; never overwrite or relabel the failed stage receipt.',
            'Validation21 commands remain unstarted and require separate original-command execution and receipts before all190 prerequisites can clear.',
            'Do not calculate replacement family BH from700 partial pairs or certify1200/41/190 completion from partial receipts.'
        ],
        'consumed_artifacts': SEEN, 'all_consumed_artifact_hashes_unchanged': True,
        'resource_plan': {'hash_buffer_bytes': 65536, 'launch_internal_free_minimum': 256 * 1024**2,
                          'emergency_internal_floor': MIN_FREE, 'maximum_RSS_bytes': MAX_RSS,
                          'deadline_seconds': DEADLINE, 'observed_peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          'elapsed_seconds': time.monotonic() - START},
        'checker_sha256': sha(Path(__file__))
    }
    assert len(fits) == 800 and len(incomplete) == 26 and array_count == 3578
    with outputs[1].open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(table)
    result['table_sha256'] = sha(outputs[1])
    with outputs[0].open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['partial_snoring_artifacts', 'consumed_artifacts', 'incomplete_original_jobs', 'continuation_requirements', 'monitor', 'numerical_failures']}, indent=2))
    print(json.dumps({'numerical_failures': failures}, indent=2))


if __name__ == '__main__':
    main()
