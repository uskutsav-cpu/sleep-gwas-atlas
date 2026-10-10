#!/usr/bin/env python3
"""Independent metadata/hash preflight: no GWAS input reads or native imports."""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import time

PACKAGE = Path(__file__).resolve().parents[1]
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
PLAN = SSD / 'sensitivity_operational_plan_v4_3_1.json'
EXPECTED = '05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
CONTROL = SSD / 'proofs/sensitivity_correction_controls_v4_3_1.json'
EXPECTED_CONTROL = 'ea86dcfd8f9e6dedb6368d6ac198288cadbba6e0e0e6f140fc6d51e210c33be8'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while True:
            data = f.read(65536)
            if not data: break
            h.update(data)
    return h.hexdigest()

def load(path): return json.loads(Path(path).read_text())
def check(ok, label):
    if not ok: raise AssertionError(label)

start = time.monotonic()
check(digest(PLAN) == EXPECTED and digest(CONTROL) == EXPECTED_CONTROL, 'frozen plan/control')
plan = load(PLAN)
old = load(plan['superseded_operational_plan'])
member = load(plan['scientific_member_manifest'])
check(digest(plan['scientific_member_manifest']) == plan['scientific_member_manifest_sha256'], 'frozen scientific matrix')
check(len(plan['jobs']) == len(member['jobs']) == 26 and sum(j['estimates'] for j in plan['jobs']) == 62, '26/62')
check(len(plan['audit_jobs']) == 52 and sum(j['expected_identities'] for j in plan['audit_jobs']) == 124, '52/124')
for kind in ('jobs', 'audit_jobs'):
    check(len(plan[kind]) == len(old[kind]), 'matrix lengths')
    for new, prior in zip(plan[kind], old[kind]):
        for key, value in prior.items():
            if key not in ('out_prefix', 'ldsc_args'): check(new[key] == value, 'matrix member changed: ' + key)
        argv = new['ldsc_args'][:]
        argv[argv.index('--out')+1] = prior['out_prefix']
        check(argv == prior['ldsc_args'], 'native options changed')
for new, frozen in zip(plan['jobs'], member['jobs']):
    for key in ('job_id', 'kind', 'inputs', 'estimates'): check(new[key] == frozen[key], 'scientific member')
    args = new['ldsc_args']
    expected_args = ['--'+new['kind'], ','.join(new['inputs']), '--ref-ld-chr',plan['reference_prefix'], '--w-ld-chr',plan['reference_prefix'], '--n-blocks','200','--print-delete-vals', '--out',new['out_prefix']] + frozen['options']
    check(args == expected_args, 'exact scientific options')
for key in ('input_sha256','guard','baseline_dependency_sha256','original_190_jobs','native_command_count','new_fitted_estimates','stock_merge_only_command_count','stock_merge_only_intersection_identities'):
    check(plan[key] == old[key], 'frozen policy changed: '+key)
check(not (set(plan['dependencies_sha256']) & set(plan['input_sha256'])), 'GWAS dependency overlap')
bound = {}
for name, expected in plan['dependencies_sha256'].items():
    actual = digest(name)
    check(actual == expected, 'bound identity changed')
    bound[name] = {'sha256':actual,'bytes':Path(name).stat().st_size}
    if name.endswith('.py'): ast.parse(Path(name).read_text(), filename=name)
proof = load(plan['derivative_proof'])
check(digest(plan['derivative_proof']) == plan['derivative_proof_sha256'], 'derivative proof identity')
check(proof['manifest_sha256'] == plan['scientific_member_manifest_sha256'] and proof['status'] == 'RESULT_FREE_ONLY_FINITE_N_CHANGED' and proof['estimator_calls'] == proof['scientific_filter_changes'] == 0, 'derivative proof scope')
check(len(proof['results']) == 2, 'two derivatives')
for row in proof['results']:
    check(row['source_sha256_before'] == row['source_sha256_after'] == plan['baseline_input_sha256'][row['source']], 'original derivative source identity')
    check(row['derivative_sha256'] == plan['input_sha256'][row['derivative']], 'derivative bound identity')
    check(row['counts']['template_rows'] == row['full_stream_verified_rows'] == 1217311, 'full template rows')
    check(row['counts']['finite_N_changed_rows'] + row['counts']['nonfinite_N_preserved_rows'] == 1217311, 'N reconciliation')
    check(row['original_and_derivative_invariant_match'] and row['literal_missingness_preserved'] and row['gzip']['independent_second_compression_matches'], 'derivative proof flags')
    expected_n = 376169 if Path(row['source']).name == 'ms.sumstats.gz' else 290130
    check(row['new_finite_N'] == expected_n, 'N only convention')
normal = dict(plan['baseline_dependency_sha256'])
relocation = plan['baseline_dependency_relocation']
normal[relocation['relocated_path']] = normal.pop(relocation['original_path'])
check(normal == plan['baseline_relocated_dependency_sha256'] and len(normal) == 56 and relocation['changed_path_keys'] == 1, 'sole wrapper relocation')
binding = plan['core_precision_adjudication']
for name, expected in binding['evidence_sha256'].items(): check(digest(name) == expected, 'precision evidence hash')
comparison = load(binding['comparison'])
check(comparison['arithmetic_failures'] == [binding['exact_preserved_failure']], 'only exact snoring BMI failure')
check(binding['exact_preserved_failure']['identity'] == 'snoring__bmi' and binding['exact_preserved_failure']['p'] is False, 'exact exception identity')
adjudication = load(binding['adjudication'])
check(adjudication['fixed_tolerance_unchanged'] and adjudication['initial_collator_failure_retained'] and not adjudication['native_meaningful_numerical_discrepancy_identified'], 'precision qualification')
control = load(CONTROL)
verifier = PACKAGE/'scripts/verify_sensitivity_corrections_v4_3_1.py'
check(control['plan_sha256'] == EXPECTED and control['verifier_sha256'] == digest(verifier), 'control binding')
check(control['status'] == 'CORRECTION_BINDINGS_AND_FAILURE_CONTROLS_PASS' and control['check_count'] == len(control['checks']) == 1082, 'control counts')
check(control['real_worker_subprocesses_launched'] == control['estimator_calls'] == control['stock_merge_audits_launched'] == 0, 'control scope')
counts = {s: {'expected':0,'present':0,'terminal_monitor_present':(PACKAGE/'logs'/(s+'_native_monitor_receipt_v4.json')).exists()} for s in ('core','extension','validation')}
for job in plan['original_190_jobs']:
    counts[job['stage']]['expected'] += 1
    path = Path(plan['baseline_support_package'])/'native'/(job['stage']+'_reproduction_v1')/(job['job_id']+'.execution_receipt.json')
    counts[job['stage']]['present'] += int(path.exists())
for name, r in bound.items(): check(digest(name) == r['sha256'], 'changed during review')
check(digest(PLAN) == EXPECTED and digest(CONTROL) == EXPECTED_CONTROL, 'plan/control changed during review')
receipt = {'status':'PASS_OPERATIONAL_CORRECTION_PENDING_ALL190_BASELINE_PREREQUISITES',
    'completed_utc':datetime.now(timezone.utc).isoformat(), 'plan_sha256':EXPECTED, 'control_sha256':EXPECTED_CONTROL,
    'checker_sha256':digest(__file__), 'before_after_dependency_hashes_pass':True,'dependencies_count':len(bound),'dependencies':bound,
    'scientific_matrix':{'native_commands':26,'new_estimates':62,'audit_commands':52,'audit_identity_instances':124,'sole_option_change_from_old_plan':'output namespace'},
    'input_identities_checked_as_metadata_only':True,'GWAS_input_body_reads':0,'derivative_proof_independently_stream_reproduced':False,
    'derivative_proof_hash_and_code_review_only':True,'derivative_metadata':proof['results'],
    'exact_snoring_BMI_precision_exception_evidence_pass':True,'control_check_count':1082,'controls_rerun_by_reviewer':False,
    'actual_baseline_presence_snapshot':counts,'all190_prerequisite_admitted':False,
    'remaining_operational_blockers':[],
    'native_workers_or_audits_launched_by_reviewer':0,'scientific_validity_or_cross_fit_alignment_claimed':False,
    'elapsed_seconds':time.monotonic()-start,'peak_RSS_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
out = PACKAGE/'reviews/independent_sensitivity_v4_3_1_binding_receipt.json'
with out.open('x') as f: json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({k:receipt[k] for k in ('status','dependencies_count','actual_baseline_presence_snapshot','elapsed_seconds')}))
