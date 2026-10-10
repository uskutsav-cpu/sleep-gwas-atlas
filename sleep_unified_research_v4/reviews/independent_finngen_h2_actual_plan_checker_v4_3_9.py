"""Actual one-job/producer metadata check; no body rereads or science worker."""
import hashlib
import importlib.util
import json
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
load = lambda p: json.loads(Path(p).read_text())
def regular(p):
    p = Path(p)
    assert p.is_file() and not p.is_symlink() and not any(q.is_symlink() for q in p.parents)
    return p
plan_path = P / 'manifests/finngen_observed_h2_plan_v4_3_9.json'
expected = '8167144ca070f30738d8671ce951d71078ae001926ba6c5abd89d3aaaa2a54a1'
assert sha(regular(plan_path)) == expected
plan = load(plan_path)
candidate = load(R / 'independent_finngen_h2_candidate_review_v4_3_9.json')
assert sha(R / 'independent_finngen_h2_candidate_review_seal_v4_3_9.json') == '401c5e664561ca93e30138e4c41d691b56837b9903a583df722cd67de8b3e555'
assert candidate['conditional_candidate_pass'] and len(candidate['controls']) == 22
assert all(c['passed'] for c in candidate['controls'])
for path, digest in candidate['fixed_draft_producer_consumer_sha256'].items():
    assert sha(regular(path)) == digest
assert plan['stage'] == 'observed_h2' and plan['scope'] == 'FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY'
assert plan['new_rg_commands'] == 0 and plan['scientific_source_admitted'] is False
assert plan['independent_replication_established'] is False and len(plan['jobs']) == 1
assert plan['all_190_original_commands_and_numerical_reviews_required'] is True
job = plan['jobs'][0]
assert job['kind'] == 'h2' and job['estimates'] == 1 and job['inputs'] == [plan['derivative']]
assert job['out_prefix'] == job['output_prefix']
assert job['result_receipt'] == job['output_prefix'] + '.full_precision.json'
base = load(plan['baseline_gate_plan'])
assert sha(plan['baseline_gate_plan']) == plan['baseline_gate_plan_sha256'] == '05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
assert plan['environment'] == base['environment'] and plan['python'] == base['python']
assert plan['ldsc_dir'] == base['ldsc_dir']
args = ['--h2', plan['derivative'], '--ref-ld-chr', base['reference_prefix'], '--w-ld-chr',
        base['reference_prefix'], '--n-blocks', '200', '--print-delete-vals', '--out', job['output_prefix']]
assert job['ldsc_args'] == args
assert job['command_template'] == [plan['python'], '-u', str(P / 'scripts/38_sensitivity_ldsc_capture.py'),
    '--plan', str(plan_path), '--plan-sha', '{PLAN_SHA256}', '--job-id', job['job_id']] + args
assert plan['guard'] == dict(worker_count=1, BLAS_threads=1, internal_floor_bytes=3 << 30,
    SSD_floor_bytes=5 << 30, observed_aggregate_worker_RSS_limit_bytes=2 << 30,
    new_output_limit_bytes=2 << 30, deadline_seconds=7200, poll_seconds=2)
assert plan['global_reservation_bytes'] == 300 << 30
assert plan['assumed_effective_N'] == 4 / (1 / 51643 + 1 / 446273)
assert len(plan['required_independent_review_paths']) == 7
assert plan['required_independent_review_paths'][-3:] == [str(R / n) for n in [
    'independent_finngen_h2_prelaunch_v4_3_9.md', 'independent_finngen_h2_prelaunch_v4_3_9.json',
    'independent_finngen_h2_prelaunch_seal_v4_3_9.json']]
proof = plan['prior_preprocessing_terminal']
assert proof['plan_sha256'] == '78f0faea56024be8eca65421bf655c0f7c224afe6f029ed3ca0cea22491781b8'
for path, digest in proof['metadata_sha256'].items():
    assert sha(regular(path)) == digest
prior = load(proof['plan_path'])
stage = load(proof['stage_receipt'])
assert stage['plan_sha256'] == proof['plan_sha256']
assert stage['status'] == 'QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED'
assert stage['new_rg_commands'] == 0 and stage['owned_cleanup_verified'] is True
assert stage['biological_or_replication_claim_admitted'] is False
assert stage['result_sha256'] == proof['result_sha256']
audit = load(R / 'independent_finngen_actual_derivative_receipt_v4_3_8.json')
assert audit['status'] == 'PASS' and audit['derivative'] == plan['derivative']
assert plan['input_sha256'] == {plan['derivative']: audit['derivative_sha256']}
assert proof['result_sha256'][plan['derivative']] == audit['derivative_sha256']
for path, digest in proof['result_sha256'].items():
    if path != plan['derivative']: assert sha(regular(path)) == digest
worker_path = prior['jobs'][0]['worker_receipt']
worker = load(worker_path)
cmd = prior['jobs'][0]['command_template']
fd_index = cmd.index('--inherited-heavy-lock-fd') + 1
fd = worker['command'][fd_index]
assert fd.isdecimal()
assert worker['command'] == [proof['plan_sha256'] if x == '{PLAN_SHA256}' else fd if x == '{HEAVY_LOCK_FD}' else x for x in cmd]
assert worker['status'] == 'WORKER_COMPLETE_VERIFIED' and worker['plan_sha256'] == proof['plan_sha256']
assert worker['returncode'] == 0 and worker['stop_reason'] is None and worker['metadata_errors'] == []
assert worker['plan_unchanged'] is True and worker['owned_cleanup_verified'] is True
assert worker['process_group_teardown']['remaining_group_members'] == []
assert not worker['process_group_teardown'].get('cleanup_error')
assert worker['output_sha256'] == {k:v for k,v in proof['result_sha256'].items() if k != worker_path}
for path in [worker_path, proof['stage_receipt'], prior['terminal_pending']]:
    marker = Path(str(path) + '.failure.json')
    assert not marker.exists() and not marker.is_symlink()
spec = importlib.util.spec_from_file_location('_read_only_terminal2_actual', P / 'scripts/terminal_commit_common_v2.py')
terminal = importlib.util.module_from_spec(spec); spec.loader.exec_module(terminal)
binding = dict(plan_sha256=proof['plan_sha256'], admission_sha256=stage['admission_sha256'],
               executor_sha256=proof['executor_sha256'])
terminal.require_committed(prior['terminal_pending'], prior['terminal_seal'], binding,
                          {proof['stage_receipt']:proof['stage_receipt_sha256']})
assert sha(prior['terminal_seal']) == proof['terminal_seal_sha256']
critical = ['38_sensitivity_ldsc_capture.py', 'sensitivity_capture_common.py', 'sensitivity_executor_v4_4.py',
    '30_prepare_and_run_ssd_native_campaign.py', 'canonical_calibration_common_v4_5.py',
    'canonical_calibration_common_v4_3.py', 'native_stage_completion_v4_3.py', 'extension_replay_common_v4.py',
    'terminal_commit_common_v2.py', '58_run_finngen_feasibility_stage_v3_9.py',
    '59_prepare_finngen_feasibility_stage_v3_9.py', '60_verify_finngen_h2_diagnostic_v5.py']
for name in critical:
    path = P / 'scripts' / name
    assert sha(regular(path)) == plan['dependencies_sha256'][str(path)]
assert all(plan['dependencies_sha256'][p] == d for p, d in base['dependencies_sha256'].items())
assert len(plan['dependencies_sha256']) == 2353
assert all(Path(p).is_absolute() for p in plan['dependencies_sha256'])
counts = [sum('/' + key + '/' in k for k in stage['historical_baseline_receipt_sha256'] if k.endswith('.execution_receipt.json'))
          for key in ['core_reproduction_v1', 'extension_reproduction_v1', 'validation_reproduction_v1']]
assert counts == [57, 112, 21]
assert sha(plan_path) == expected
print(json.dumps(dict(status='ACTUAL_ONE_H2_PLAN_AND_CURRENT_PRODUCER_PASS', plan_sha256=expected,
    original_baseline_member_counts=counts, current_prior_terminal_verified=True,
    current_prior_result_artifacts=5, inherited_candidate_helper_controls=22,
    critical_executable_dependencies_checked=len(critical), derivative_reopened=False,
    raw_reference_runtime_bodies_or_censuses_or_workers_or_locks=0)))
