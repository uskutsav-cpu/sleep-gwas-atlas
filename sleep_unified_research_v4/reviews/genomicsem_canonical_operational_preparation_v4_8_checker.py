"""Author preparation audit only; no worker, numerical import, or admission."""
import ast
import difflib
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
S = P / 'scripts'
NAMES = ['42_prepare_native_canonical_calibration', 'canonical_calibration_common',
         'canonical_calibration_capture', 'canonical_calibration_arithmetic']
EXPECTED_PLAN = 'c9ad3842526a61ab24ff8fd39222c408fb611ba2bbf26e7ef0959f1e86d8de51'
EXPECTED_FAILURE = 'c80990a119349b0deda5b4d1c6ee15ceb1d16f3c79d0f9f457894b8998674215'


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def write_new(path, value):
    with path.open('x') as f:
        f.write(value)


def main():
    checks = []
    def require(label, condition):
        checks.append({'label': label, 'pass': bool(condition)})
        if not condition:
            raise AssertionError(label)

    plan_path = P / 'manifests/native_canonical_calibration_plan_v4_8.json'
    require('frozen_prepared_plan_sha', sha(plan_path) == EXPECTED_PLAN)
    plan = json.loads(plan_path.read_text())
    old_path = Path(plan['preserved_v4_7_plan'])
    require('old_plan_sha', sha(old_path) == plan['preserved_v4_7_plan_sha256'])
    old = json.loads(old_path.read_text())
    failure_path = Path(plan['preserved_failed_attempt_receipt'])
    require('failed_preserved_inventory_sha', sha(failure_path) == EXPECTED_FAILURE)
    failure = json.loads(failure_path.read_text())
    require('failure_is_preserved_not_credited', failure['status'] == 'FAILED_PRESERVED'
            and failure['completed_native_controls_credited'] == 0
            and failure['method_calibration_pass'] is False)
    require('exact_preserved_inventory_count', failure['file_count'] == 92
            and len(failure['all_preserved_file_identity']) == 92)
    total = 0
    for path, identity in failure['all_preserved_file_identity'].items():
        p = Path(path)
        require('preserved_regular:' + path, p.is_file() and not p.is_symlink()
                and not any(q.is_symlink() for q in p.parents))
        require('preserved_size_sha:' + path,
                p.stat().st_size == identity['bytes'] and sha(p) == identity['sha256'])
        total += identity['bytes']
    require('exact_preserved_bytes', total == failure['preserved_bytes'] == 127954369)
    require('new_plan_preservation_map_matches',
            plan['preserved_v4_6_file_sha256'] == {
                k: v['sha256'] for k, v in failure['all_preserved_file_identity'].items()})
    require('prior21_review_map_matches',
            plan['preserved_v4_6_review_artifact_sha256'] == failure['old_review_artifact_sha256']
            and len(failure['old_review_artifact_sha256']) == 21)
    for path, digest in failure['old_review_artifact_sha256'].items():
        require('prior21_review_sha:' + path, sha(path) == digest)
    for label, path_key, hash_key, map_key, count in [
            ('prior_v7_complete_rejection', 'preserved_v4_7_complete_rejection_seal', 'preserved_v4_7_complete_rejection_seal_sha256', 'preserved_v4_7_rejection_artifact_sha256', 25),
            ('prior_v7_author_preparation', 'preserved_v4_7_author_seal', 'preserved_v4_7_author_seal_sha256', 'preserved_v4_7_author_artifact_sha256', 21)]:
        seal_path = Path(plan[path_key])
        require(label+'_seal_sha', sha(seal_path) == plan[hash_key])
        seal = json.loads(seal_path.read_text())
        require(label+'_exact_artifact_map', seal['file_sha256'] == plan[map_key] and len(seal['file_sha256']) == count)
        for path, digest in seal['file_sha256'].items():
            require(label+'_artifact_sha:' + path, sha(path) == digest and plan['dependencies_sha256'].get(path) == digest)
        require(label+'_seal_in_live_dependencies', plan['dependencies_sha256'].get(str(seal_path)) == plan[hash_key])
    require('prior23_rejection_seal_retained',
            str(R/'independent_canonical_executor_prelaunch_review_v4_7_seal.json') in plan['preserved_v4_7_rejection_artifact_sha256'])
    inherited = json.loads(Path(plan['preparation_deferred_identity_addendum']).read_text())
    require('exact_sealed_deferred_addendum_sha', sha(Path(plan['preparation_deferred_identity_addendum'])) == plan['preparation_deferred_identity_addendum_sha256'])
    require('exact50_inherited_not_freshly_hashed', len(plan['preparation_deferred_dependency_sha256']) == plan['preparation_deferred_dependency_count'] == 50
            and plan['preparation_deferred_dependency_sha256'] == inherited['remaining50_source_reference_runtime_dependency_sha256'])
    for path, digest in plan['preparation_deferred_dependency_sha256'].items():
        require('deferred_identity_retained_without_body_read:' + path, plan['dependencies_sha256'].get(path) == old['dependencies_sha256'].get(path) == digest)
    require('future_admission_and_190_gates_unchanged_AST',
            ast.dump(next(n for n in ast.parse((S/'canonical_calibration_common_v4_7.py').read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'admit'), include_attributes=False)
            == ast.dump(next(n for n in ast.parse((S/'canonical_calibration_common_v4_8.py').read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'admit'), include_attributes=False))
    for map_path, map_field in [(R/'genomicsem_canonical_operational_preparation_v4_7_fixture_map.json', None),
                                 (R/'independent_canonical_v4_7_preserved_fixture_map.json', 'regular_file_sha256')]:
        old_fixture_map = json.loads(map_path.read_text())
        old_fixture_map = old_fixture_map[map_field] if map_field else {k: v['sha256'] for k, v in old_fixture_map.items()}
        for path, digest in old_fixture_map.items():
            require('preserved_v7_tiny_witness_file_sha:' + path, sha(path) == digest)
    for path, digest in old['dependencies_sha256'].items():
        require('original_dependency_identity_retained:' + path,
                plan['dependencies_sha256'].get(path) == digest)
    require('exact16_same_jobs_only_output_relocated', len(old['jobs']) == len(plan['jobs']) == 16
            and json.loads(json.dumps(old['jobs']).replace(old['ssd_output_root'], plan['ssd_output_root'])) == plan['jobs'])
    for key in ['scope', 'allowed_control_ids', 'historical_plan_path', 'historical_plan_sha256',
                'block_receipt_path', 'block_receipt_sha256', 'build', 'interval_tsv_sha256',
                'ordered_all_reference_records_sha256', 'ordered_construction_records_sha256',
                'coordinate_map_path', 'coordinate_map_sha256', 'partition_amendment',
                'static_construction_exclusions', 'reference_sha256', 'inputs_current_sha256',
                'python', 'ldsc_dir', 'environment', 'worker_count', 'blas_threads',
                'shared_lock_path', 'shared_lock_policy', 'internal_floor_bytes', 'ssd_floor_bytes',
                'worker_rss_limit_bytes', 'output_limit_bytes', 'per_worker_seconds_limit',
                'runtime_poll_seconds', 'resource_storage_assumption', 'operational_resource_amendment',
                'arithmetic', 'cross_arm_identity_rule', 'allow_41_covariance_outcomes',
                'allow_calibrated_biological_p_values', 'realistic_LD_sampling_calibration_pass',
                'remaining_method_gate', 'global_resource_ledger_path', 'global_resource_ledger_sha256',
                'global_reservation_bytes']:
        require('unchanged_scientific_or_resource_field:' + key, plan[key] == old[key])
    require('unchanged_explicit_caps', plan['output_limit_bytes'] == 8*(1<<30)
            and plan['all_retained_canonical_epochs_output_cap_bytes'] == 8*(1<<30)
            and plan['global_reservation_bytes'] == 300*(1<<30)
            and plan['internal_floor_bytes'] == 3*(1<<30)
            and plan['ssd_floor_bytes'] == 5*(1<<30)
            and plan['worker_rss_limit_bytes'] == 2*(1<<30)
            and plan['per_worker_seconds_limit'] == 7200)
    require('no_execution_admission_or_namespace', plan['preparation_only'] is True
            and plan['execution_admitted'] is False and plan['native_jobs_launched'] == 0
            and not (P/'manifests/native_canonical_calibration_admission_v4_8.json').exists()
            and not Path(plan['ssd_output_root']).exists())
    require('old_pending_remains', Path(failure['pending_path']).is_file())
    fixture_path = R/'canonical_calibration_monitor_consumer_fixture_receipt_v4_8_2.json'
    fixture = json.loads(fixture_path.read_text())
    require('author_fixture16cases63assertions', fixture['all_checks_pass']
            and fixture['case_count'] == 16 and fixture['check_count'] == 63
            and fixture['workers_launched'] == 0 and fixture['mutex_acquisitions'] == 0
            and fixture['GWAS_reference_runtime_body_reads'] == 0)
    sources = {}
    diff = []
    for name in NAMES:
        a = S/(name+'_v4_7.py'); b = S/(name+'_v4_8.py')
        for f in [a, b]:
            sources[str(f)] = sha(f)
            require('author_fixture_exact_source:' + f.name,
                    fixture['source_sha256'][str(f)] == sources[str(f)])
            compile(f.read_text(), str(f), 'exec')
        require('new_source_bound_in_plan:' + b.name,
                plan['dependencies_sha256'][str(b)] == sources[str(b)])
        diff.extend(difflib.unified_diff(a.read_text().splitlines(True), b.read_text().splitlines(True),
                                      fromfile=str(a), tofile=str(b)))
    def func(path, name):
        return next(n for n in ast.parse(path.read_text()).body
                    if isinstance(n, ast.FunctionDef) and n.name == name)
    require('arithmetic_validate_exact_AST',
            ast.dump(func(S/'canonical_calibration_arithmetic_v4_7.py', 'validate'), include_attributes=False)
            == ast.dump(func(S/'canonical_calibration_arithmetic_v4_8.py', 'validate'), include_attributes=False))
    require('capture_exact_import_only',
            (S/'canonical_calibration_capture_v4_7.py').read_text().replace(
                'canonical_calibration_common_v4_7', 'canonical_calibration_common_v4_8')
            == (S/'canonical_calibration_capture_v4_8.py').read_text())
    diff_path = R/'genomicsem_canonical_operational_preparation_v4_8.diff'
    write_new(diff_path, ''.join(diff))
    fixture_map = {str(f): {'bytes': f.stat().st_size, 'sha256': sha(f)}
                   for f in sorted((R/'canonical_calibration_monitor_consumer_fixtures_v4_8_2').rglob('*'))
                   if f.is_file() and not f.is_symlink()}
    fixture_map_path = R/'genomicsem_canonical_operational_preparation_v4_8_fixture_map.json'
    write_new(fixture_map_path, json.dumps(fixture_map, indent=2, allow_nan=False)+'\n')
    receipt = {'schema': 'author_canonical_operational_preparation_audit_v4_8',
               'completed_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'PREPARED_FOR_SEPARATE_INDEPENDENT_REVIEW_NOT_ADMITTED',
               'plan_path': str(plan_path), 'plan_sha256': sha(plan_path),
               'failed_preserved_inventory_path': str(failure_path), 'failed_preserved_inventory_sha256': sha(failure_path),
               'source_sha256': sources, 'checks': checks, 'check_count': len(checks),
               'all_checks_pass': all(c['pass'] for c in checks),
               'fixture_receipt_path': str(fixture_path), 'fixture_receipt_sha256': sha(fixture_path),
               'fixture_map_path': str(fixture_map_path), 'fixture_map_sha256': sha(fixture_map_path),
               'diff_path': str(diff_path), 'diff_sha256': sha(diff_path),
               'workers_launched': 0, 'heavy_lock_acquired': False,
               'new_outcome_analysis': False, 'independent_review_or_root_admission': False,
               'preserved92_bytes_rehashed_only': total,
               'deferred50_dependencies_not_freshly_hashed': plan['preparation_deferred_dependency_sha256'],
               'preparation_dependency_qualification': plan['preparation_dependency_qualification'],
               'preserved_v7_complete_rejection_seal_sha256': plan['preserved_v4_7_complete_rejection_seal_sha256'],
               'preserved_v7_author_preparation_seal_sha256': plan['preserved_v4_7_author_seal_sha256'],
               'actual_producer_stop_guarantees_unchanged': True,
               'ordinary_production_false_success_claim': False,
               'author_harness_failure_preserved': 'Initial tiny final-consumer healthy fixture omitted required fake plan.python; retained initial checker/log/fixtures; corrected fixture _v4_8_2 adds only that fake key and new private output paths. Candidate sources unchanged.',
               'no_original_GWAS_or_reference_body_read_in_this_checker': True,
               'lock_probe_qualification': 'Prior read-only NB probe BUSY_NO_LOCK_ACQUIRED; parent reports other independently admitted corelarge35 stage. Own failed group97558 absence was established in frozen inventory; no new process or lock probe here.',
               'remaining_gates': ['Separate independent exact v4_8 source/plan review',
                   'Root admission and fresh unchanged guards/mutex availability before any execution',
                   'Actual16 control workers plus17th arithmetic success and Terminal2 completion',
                   'Independent native-control numerical verification',
                   'Realistic LD-aware sampling calibration before any biological41 covariance release']}
    out = R/'genomicsem_canonical_operational_preparation_v4_8.json'
    write_new(out, json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': receipt['status'], 'check_count': len(checks), 'all_checks_pass': True,
                      'receipt': str(out), 'receipt_sha256': sha(out), 'diff_sha256': sha(diff_path)}))


if __name__ == '__main__':
    main()
