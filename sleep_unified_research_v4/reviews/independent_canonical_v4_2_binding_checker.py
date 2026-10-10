#!/usr/bin/env python3
"""Read-only independent metadata/hash snapshot, no executor imports or fits."""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import shutil
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
PLAN = PACKAGE / 'manifests/native_canonical_calibration_plan_v4_2.json'
EXPECTED = '5954e930b567b25b009769a1010009d494385714a6e7a2f47b7dd67c434d3917'
OUTPUT = PACKAGE / 'reviews/independent_canonical_v4_2_binding_receipt.json'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while True:
            data = handle.read(65536)
            if not data: break
            h.update(data)
    return h.hexdigest()

def load(path):
    return json.loads(Path(path).read_text())

def check(condition, label):
    if not condition: raise AssertionError(label)

start = time.monotonic()
check(shutil.disk_usage('/System/Volumes/Data').free >= 256 * 1024**2, 'review internal floor')
check(digest(PLAN) == EXPECTED, 'plan identity')
plan = load(PLAN)
old = load(PACKAGE / 'manifests/native_canonical_calibration_plan_v4.json')
check(len(plan['jobs']) == len(old['jobs']) == 16, '16 controls')
check(plan['allowed_control_ids'] == old['allowed_control_ids'], 'controls')
for before, after in zip(old['jobs'], plan['jobs']):
    for key in ('job_id', 'control_id', 'method', 'partition', 'duplicate', 'inputs', 'input_sha256'):
        check(before[key] == after[key], 'control changed: ' + key)
    left, right = before['ldsc_args'][:], after['ldsc_args'][:]
    left[left.index('--out')+1] = '<versioned-output>'
    right[right.index('--out')+1] = '<versioned-output>'
    check(left == right, 'native CLI changed')
for key in ('arithmetic', 'reference_sha256', 'coordinate_map_sha256', 'interval_tsv_sha256', 'ordered_all_reference_records_sha256', 'ordered_construction_records_sha256', 'inputs_current_sha256'):
    check(plan[key] == old[key], 'scientific policy/identity changed: ' + key)
check(plan['output_limit_bytes'] == 8*1024**3 and plan['worker_count'] == plan['blas_threads'] == 1, 'resource policy')
admission = load(PACKAGE / 'manifests/native_canonical_calibration_admission_v4_2.json')
check(admission['execution_admitted'] is False and admission['plan_sha256'] == EXPECTED, 'admission state')
identities = {}
for name, expected in {**plan['dependencies_sha256'], **plan['inputs_current_sha256']}.items():
    actual = digest(name)
    check(actual == expected, 'bound dependency/input changed: ' + name)
    identities[name] = {'sha256':actual, 'bytes':Path(name).stat().st_size}
    if name.endswith('.py'): ast.parse(Path(name).read_text(), filename=name)
fixtures = []
for stem in ('canonical_native_calibration_fixture_receipt_v4_3', 'canonical_operational_correction_fixture_receipt_v4_2', 'canonical_cleanup_and_final_guard_fixture_receipt_v4_2'):
    path = PACKAGE / 'statistical_validation' / (stem + '.json')
    record = load(path)
    fixture_keys = [k for k in record if k.endswith('fixture_sha256')]
    for key in fixture_keys:
        matching = [name for name in plan['dependencies_sha256'] if name.endswith('.py') and plan['dependencies_sha256'][name] == record[key]]
        check(bool(matching), 'fixture source not bound: ' + stem)
    fixtures.append({'path':str(path), 'sha256':digest(path), 'status':record.get('status'), 'recorded_utc':record.get('recorded_utc')})
historical_path = Path(plan['historical_plan_path'])
historical = load(historical_path)
check(digest(historical_path) == plan['historical_plan_sha256'], 'historical plan')
expected_deps = dict(historical['dependencies_sha256'])
original = str(ROOT / 'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
relocated = str(Path(historical['ssd_support_package']) / 'scripts/native_ldsc_capture.py')
value = expected_deps.pop(original)
check(digest(original) == digest(relocated) == historical['support_file_sha256'][relocated] == value, 'sole relocation identity')
check(relocated not in expected_deps, 'duplicate relocated dependency')
expected_deps[relocated] = value
historical_inputs = {r['path']:r['actual_sha256'] for r in historical['inputs_verified'] if r['match']}
counts = {stage:{'expected':0,'present_verified':0,'missing':[]} for stage in ('core','extension','validation')}
receipts, output_count, output_bytes = {}, 0, 0
for job in historical['jobs']:
    stage = job['stage']
    counts[stage]['expected'] += 1
    directory = Path(historical['ssd_support_package']) / 'native' / (stage+'_reproduction_v1')
    path = directory / (job['job_id'] + '.execution_receipt.json')
    if not path.exists():
        counts[stage]['missing'].append(job['job_id'])
        continue
    first = digest(path)
    record = load(path)
    check(record['job'] == job and record['returncode'] == 0 and record['scientific_cardinality_gate_pass'] is True and record['execution_identity_gate_pass'] is True, 'historical job success')
    check(record['dependency_sha256_before'] == record['dependency_sha256_after'] == expected_deps, 'historical exact relocated dependency map')
    expected_inputs = {name:historical_inputs[name] for name in job['inputs']}
    check(record['input_sha256'] == record['input_sha256_after'] == expected_inputs, 'historical input receipt maps')
    outputs = record['all_output_sha256']
    full = str(directory / (job['job_id']+'.full_precision.json'))
    check(outputs.get(full) == record['output_sha256'], 'full precision required')
    for name, expected in outputs.items():
        p = Path(name)
        check(p.parent.resolve() == directory.resolve() and p.name.startswith(job['job_id']), 'output scope')
        check(digest(p) == expected, 'output hash')
        output_count += 1
        output_bytes += p.stat().st_size
    check(digest(path) == first, 'receipt changed during inspection')
    receipts[str(path)] = first
    counts[stage]['present_verified'] += 1
stage_receipts = {}
for stage in counts:
    p = PACKAGE / 'logs' / (stage+'_native_monitor_receipt_v4.json')
    if not p.exists():
        stage_receipts[stage] = {'present':False}
        continue
    r = load(p)
    check(r['stage'] == stage and r['plan_sha256'] == plan['historical_plan_sha256'] and r['returncode'] == 0 and r['stop_reason'] is None and r['process_group_teardown']['remaining_group_members'] == [], 'stage monitor completion binding')
    stage_receipts[stage] = {'present':True,'sha256':digest(p)}
for name, item in identities.items():
    check(digest(name) == item['sha256'], 'bound file changed during review')
check(digest(PLAN) == EXPECTED, 'plan changed')
result = {'status':'PASS_BOUND_IDENTITIES_BUT_EXECUTION_BLOCKED_OPERATIONAL_FINDINGS',
    'completed_utc':datetime.now(timezone.utc).isoformat(), 'plan_sha256':EXPECTED,
    'checker_sha256':digest(__file__), 'dependencies_count':len(plan['dependencies_sha256']),
    'input_compressed_identity_count':len(plan['inputs_current_sha256']), 'before_after_bound_hashes_pass':True,
    'bound_identities':identities, 'fixtures':fixtures, 'exact_control_matrix_and_scientific_identity_preserved':True,
    'historical_counts_snapshot':counts, 'historical_stage_receipts_snapshot':stage_receipts,
    'present_historical_receipt_sha256':receipts, 'historical_output_hash_count':output_count,
    'historical_output_hash_bytes':output_bytes, 'all190_completed_claimed':False,
    'admission_execution_admitted':False, 'native_fits_launched_by_reviewer':0,
    'GWAS_outcome_columns_read_by_reviewer':False, 'raw_or_giant_sources_rehashed_by_reviewer':False,
    'remaining_operational_blockers':['Catchable supervisor termination is not deferred/converted to protected failure.', 'Cleanup-lock retry diagnostics can throw and escape retention.', 'Final resource/deadline gates precede all-output hashing.'],
    'cross_arm_numerical_rule_status':'Accepted as a frozen-tolerance implementation comparison; unequal arrays are not asserted bitwise identical or a proof that roundoff is their sole cause.',
    'scientific_covariance_or_realistic_LD_calibration_admitted':False,
    'elapsed_seconds':time.monotonic()-start, 'peak_rss_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
with OUTPUT.open('x') as handle:
    json.dump(result,handle,indent=2);handle.write('\n')
print(json.dumps({k:result[k] for k in ('status','dependencies_count','historical_counts_snapshot','historical_output_hash_count','historical_output_hash_bytes','elapsed_seconds')}))
