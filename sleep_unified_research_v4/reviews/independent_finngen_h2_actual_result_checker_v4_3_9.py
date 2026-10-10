"""Current verified report consumption and independent native-delete arithmetic."""
import csv
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import resource
import sys
import time

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
sys.path.insert(0, str(P / 'scripts'))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
load = lambda p: json.loads(Path(p).read_text())
started = time.monotonic()
plan_path = P / 'manifests/finngen_observed_h2_plan_v4_3_9.json'
plan_sha = '8167144ca070f30738d8671ce951d71078ae001926ba6c5abd89d3aaaa2a54a1'
report_path = P / 'statistical_validation/finngen_observed_h2_actual_diagnostic_v4_3_9.json'
assert sha(plan_path) == plan_sha
assert sha(report_path) == '6b5a10a4401b52075d38574a20e9ba6919d7b89ddaa0335f789c44a0840604a5'
consumer_path = P / 'scripts/60_verify_finngen_h2_diagnostic_v5.py'
assert sha(consumer_path) == 'e9e85f9a79d6cc07da4da9b41b3a0d6bd9b95a969b4165d31a0217e0950719e2'
spec = importlib.util.spec_from_file_location('_actual_h2_verified_report_reader', consumer_path)
consumer = importlib.util.module_from_spec(spec); spec.loader.exec_module(consumer)
consumption = consumer.require_verified_report(plan_path, plan_sha, report_path)
# Read-only consumer performs its mandatory current compressed-derivative hashes;
# it does not reopen raw source/reference bodies, decompress, fit or acquire locks.
plan, report = load(plan_path), load(report_path)
job = plan['jobs'][0]
capture = load(job['result_receipt'])
stage, worker = load(plan['stage_receipt']), load(job['worker_receipt'])
assert len(plan['jobs']) == len(capture['estimates']) == len(capture['final_intersections']) == 1
assert capture['plan_sha256'] == worker['plan_sha256'] == stage['plan_sha256'] == plan_sha
assert capture['job_id'] == job['job_id']
assert worker['command'] == [plan_sha if v == '{PLAN_SHA256}' else v for v in job['command_template']]
assert worker['metadata_errors'] == [] and worker['plan_unchanged'] is True
assert worker['status'] == 'WORKER_COMPLETE_VERIFIED' and worker['returncode'] == 0
assert worker['stop_reason'] is None and worker['owned_cleanup_verified'] is True
assert worker['process_group_teardown']['remaining_group_members'] == []
assert not worker['process_group_teardown'].get('cleanup_error')
x, intersection = capture['estimates'][0], capture['final_intersections'][0]
assert capture['input_sha256_before'] == capture['input_sha256_after'] == plan['input_sha256']
assert x['n_blocks'] == intersection['n_blocks'] == 200 and x['n_annot'] == 1
assert intersection['final_ordered_SNP_count'] == 1156359
assert intersection['final_Z_filter_removed'] == 0 and intersection['effective_two_step'] == 30
assert intersection['two_step_hsq_mask']['count'] == 1156315
assert intersection['two_step_hsq_mask']['total'] - intersection['two_step_hsq_mask']['count'] == x['twostep_filtered'] == 44
assert capture['warnings_and_errors'] == report['warning_lines'] == []
def scalar_vector(values):
    out = [v[0] if isinstance(v, list) and len(v) == 1 else v for v in values]
    assert len(out) == 200 and all(type(v) in [int, float] and math.isfinite(v) for v in out)
    return out
def decimal_variance(values):
    with localcontext() as context:
        context.prec = 70
        data = [Decimal.from_float(v) for v in values]
        mean = sum(data) / Decimal(len(data))
        var = Decimal(len(data) - 1) / Decimal(len(data)) * sum((v - mean) ** 2 for v in data)
        return float(var), float(var.sqrt())
vectors = {key:scalar_vector(x[key]) for key in ['tot_delete_values', 'part_delete_values', 'intercept_delete_values']}
checks = {}
for key, field in [('tot_delete_values', 'tot_se'), ('intercept_delete_values', 'intercept_se')]:
    variance, se = decimal_variance(vectors[key])
    assert math.isclose(se, x[field], rel_tol=1e-12, abs_tol=1e-15)
    checks[field] = dict(native=x[field], independent_decimal_SE=se, independent_decimal_variance=variance,
                        difference=se-x[field], fixed_tolerance_pass=True)
assert math.isclose(checks['tot_se']['independent_decimal_variance'], x['tot_cov'], rel_tol=1e-12, abs_tol=1e-15)
stock_checks = {}
for key, suffix in [('tot_delete_values', '.delete'), ('part_delete_values', '.part_delete')]:
    path = Path(job['out_prefix'] + suffix)
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    assert len(lines) == 200 and all(len(v) == 1 for v in lines)
    stock = [float(v[0]) for v in lines]
    assert all(math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15) for a, b in zip(stock, vectors[key]))
    stock_checks[key] = dict(rows=200, maximum_absolute_difference=max(abs(a-b) for a,b in zip(stock,vectors[key])),
                             stock_sha256=sha(path), fixed_tolerance_pass=True)
M, N = x['M'][0][0], float(format(plan['assumed_effective_N'], '.12g'))
assert math.isclose(x['tot'], x['coef'][0] * M, rel_tol=1e-12, abs_tol=1e-15)
assert math.isclose(x['tot_se'], x['coef_se'][0] * M, rel_tol=1e-12, abs_tol=1e-15)
assert math.isclose(x['coef'][0], x['jknife']['est'][0][0] / N, rel_tol=1e-12, abs_tol=1e-15)
assert all(math.isclose(a, b*M, rel_tol=1e-12, abs_tol=1e-15) for a,b in zip(vectors['tot_delete_values'],vectors['part_delete_values']))
assert all(math.isclose(v, j[1], rel_tol=1e-12, abs_tol=1e-15) for v,j in zip(vectors['intercept_delete_values'],x['jknife']['delete_values']))
d = report['diagnostic']
z = x['tot'] / x['tot_se']
ci = [x['tot'] - 1.959963984540054*x['tot_se'], x['tot'] + 1.959963984540054*x['tot_se']]
ratio, ratio_se = (x['intercept']-1)/(x['mean_chisq']-1), x['intercept_se']/(x['mean_chisq']-1)
for actual, expected in [(d['h2'],x['tot']), (d['h2_SE'],x['tot_se']), (d['h2_Z'],z),
    (d['h2_CI_lower'],ci[0]), (d['h2_CI_upper'],ci[1]), (d['intercept'],x['intercept']),
    (d['intercept_SE'],x['intercept_se']), (d['mean_chi_square'],x['mean_chisq']),
    (d['lambda_GC'],x['lambda_gc']), (x['ratio'],ratio), (x['ratio_se'],ratio_se)]:
    assert math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-15)
assert d['necessary_h2_Z_ge4'] is (z >= 4) and d['necessary_intercept_le1_2'] is (x['intercept'] <= 1.2)
assert d['source_INFO_gate_verified'] is d['per_variant_N_verified'] is d['scientific_source_admitted'] is False
assert d['fully_independent_two_trait_replication'] is False and d['new_rg_estimates'] == 0
assert report['common_cross_trait_boundaries_established'] is report['original_families_changed'] is False
table_path = report_path.with_suffix('.tsv')
table = list(csv.DictReader(table_path.open(), delimiter='\t'))
assert len(table) == 1 and set(table[0]) == set(d)
assert all(table[0][k] == str(v) for k, v in d.items())
log = Path(job['out_prefix'] + '.log').read_text()
assert 'After merging with regression SNP LD, 1156359 SNPs remain.' in log
assert 'Using two-step estimator with cutoff at 30.' in log
assert not any('WARNING' in line or 'ERROR' in line for line in log.splitlines())
receipt = dict(schema='independent_actual_finngen_h2_result_arithmetic_v4_3_9', status='PASS',
    plan_sha256=plan_sha, report_sha256=sha(report_path), consumption=consumption,
    full_precision_capture_sha256=sha(job['result_receipt']), native_stage_sha256=sha(plan['stage_receipt']),
    independent_decimal_delete_checks=checks, independent_stock_delete_checks=stock_checks,
    captured_M=M, serialized_assumed_N=N, point_coefficient_M_identity_pass=True,
    total_partition_intercept_delete_transform_identity_pass=True, scalar_table_log_concordance_pass=True,
    diagnostic=d, captured_ratio=ratio, captured_ratio_SE=ratio_se,
    final_intersection=intersection, warning_lines=[], stock_delete_export_scope='Total and partition available; intercept verified from captured native array only.',
    numerical_tolerances=dict(relative=1e-12,absolute=1e-15), actual_fits_run_by_reviewer=0,
    raw_reference_body_reads=0, derivative_decompression_or_old_suites=0,
    mandatory_consumer_current_compressed_derivative_hashes='Executed unchanged read-only consumer, including its two prior preprocessing gates.',
    elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
out = R / 'independent_finngen_h2_actual_result_arithmetic_receipt_v4_3_9.json'
with out.open('x') as f: json.dump(receipt,f,indent=2); f.write('\n')
print(json.dumps(dict(status=receipt['status'], h2=x['tot'], SE=x['tot_se'], Z=z, CI=ci,
    intercept=x['intercept'], intercept_SE=x['intercept_se'], ratio=ratio, stock_delete_rows=400,
    captured_intercept_delete_rows=200, max_RSS_bytes=receipt['max_RSS_bytes'],elapsed_seconds=receipt['elapsed_seconds'])))
