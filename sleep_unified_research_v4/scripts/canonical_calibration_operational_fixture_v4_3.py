#!/usr/bin/env python3
"""No-fit fixtures for corrected history, teardown, lock, and method guards."""
import ast
import copy
import fcntl
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace

import canonical_calibration_common_v4_3 as common


def load_executor():
    path = common.PACKAGE / 'scripts/42_prepare_native_canonical_calibration_v4_3.py'
    spec = importlib.util.spec_from_file_location('canonical_mock_operational_executor', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def method_fixture():
    """Execute exactly two pinned functions with every I/O/fit operation mocked."""
    source_path = common.ROOT.parent / 'ldsc-code/ldscore/sumstats.py'
    native = json.loads((common.PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json').read_text())
    assert common.sha(source_path) == native['dependencies_sha256'][str(source_path)]
    source = ast.parse(source_path.read_text())
    functions = [node for node in source.body if isinstance(node, ast.FunctionDef) and node.name in ('estimate_rg', '_split_or_none')]
    observations = []
    def pair(rows, args, log, M, refs, weights, index):
        observations.append({'two_step': args.two_step, 'intercept_h2': args.intercept_h2, 'intercept_gencov': args.intercept_gencov})
        return 'MOCK_RETURN_NO_FIT'
    scope = {'copy': copy, '_parse_rg': lambda paths: (['fixture_A', 'fixture_B'], ['A', 'B']), '_check_arg_len': lambda *args: None,
             '_read_ld_sumstats': lambda *args, **kw: (SimpleNamespace(shape=(1, 1)), 'W', ['L2'], 'MOCK_ROWS', None),
             '_read_other_sumstats': lambda *args: 'MOCK_FINAL_ROWS', '_rg': pair,
             '_print_gencor': lambda *args: None, '_get_rg_table': lambda *args: 'MOCK_TABLE_NO_STATISTICS'}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source_path), 'exec'), scope)
    log = SimpleNamespace(log=lambda message: None)
    for twostep in (None, 30):
        args = SimpleNamespace(rg='fixture_A,fixture_B', intercept_h2=None, intercept_gencov=None, samp_prev=None, pop_prev=None,
                               no_intercept=False, out='fixture', two_step=twostep, print_cov=False, print_delete_vals=False)
        result = scope['estimate_rg'](args, log)
        assert result == ['MOCK_RETURN_NO_FIT']
        assert observations[-1] == {'two_step': twostep, 'intercept_h2': [None, None], 'intercept_gencov': [None, None]}
    return {'exact_pinned_function_source_sha256': common.sha(source_path), 'observations': observations,
            'one_step_default_condition_false_after_list_normalization': True, 'free_intercepts_preserved': True, 'estimator_calls': 0}


def history_fixture(root):
    original_package = common.PACKAGE
    package, support = root / 'history_package', root / 'history_support'
    (package / 'logs').mkdir(parents=True)
    (support / 'scripts').mkdir(parents=True)
    original = common.ROOT / 'sleep_unified_research_v1/scripts/native_ldsc_capture.py'
    capture = support / 'scripts/native_ldsc_capture.py'
    capture.write_bytes(original.read_bytes())
    dependency_value = common.sha(original)
    deps = {str(original): dependency_value, 'other_exact_dependency': 'b'*64}
    relocated = {str(capture): dependency_value, 'other_exact_dependency': 'b'*64}
    jobs = []
    for i in range(190):
        stage = ('core', 'extension', 'validation')[i % 3]
        jobs.append({'job_id': 'fixture_' + str(i), 'stage': stage, 'inputs': ['fixture_input']})
    historical = {'jobs': jobs, 'ssd_support_package': str(support), 'dependencies_sha256': deps,
                  'support_file_sha256': {str(capture): dependency_value},
                  'inputs_verified': [{'path': 'fixture_input', 'actual_sha256': 'c'*64, 'match': True}]}
    native_path = root / 'mock_history_plan.json'
    common.write_new(native_path, historical)
    plan = {'historical_plan_path': str(native_path), 'historical_plan_sha256': common.sha(native_path)}
    first_output = None
    for job in jobs:
        directory = support / 'native' / (job['stage'] + '_reproduction_v1')
        directory.mkdir(parents=True, exist_ok=True)
        output = directory / (job['job_id'] + '.full_precision.json')
        common.write_new(output, {'mock': True, 'no_estimates': True})
        record = {'job': job, 'returncode': 0, 'scientific_cardinality_gate_pass': True, 'execution_identity_gate_pass': True,
                  'dependency_sha256_before': relocated, 'dependency_sha256_after': relocated,
                  'input_sha256': {'fixture_input': 'c'*64}, 'input_sha256_after': {'fixture_input': 'c'*64},
                  'all_output_sha256': {str(output): common.sha(output)}, 'output_sha256': common.sha(output)}
        common.write_new(directory / (job['job_id'] + '.execution_receipt.json'), record)
        if first_output is None: first_output = output
    for stage in ('core', 'extension', 'validation'):
        common.write_new(package / 'logs' / (stage + '_native_monitor_receipt_v4.json'),
                         {'stage': stage, 'returncode': 0, 'stop_reason': None, 'plan_sha256': plan['historical_plan_sha256'],
                          'process_group_teardown': {'remaining_group_members': []}})
    common.PACKAGE = package
    try:
        receipts = common.historical_completion(plan)
        assert len(receipts) == 193
        content = first_output.read_bytes()
        first_output.write_text('fixture_tamper')
        try: common.historical_completion(plan)
        except RuntimeError as error: assert str(error).startswith('BOUND_IDENTITY_CHANGED')
        else: raise AssertionError('Tampered historical output admitted')
        first_output.write_bytes(content)
        monitor = package / 'logs/core_native_monitor_receipt_v4.json'
        record = json.loads(monitor.read_text())
        record['plan_sha256'] = 'd'*64
        monitor.write_text(json.dumps(record))
        try: common.historical_completion(plan)
        except RuntimeError as error: assert str(error) == 'HISTORICAL_STAGE_MONITOR_PLAN_BINDING_MISMATCH: core'
        else: raise AssertionError('Wrong historical monitor plan admitted')
    finally:
        common.PACKAGE = original_package
    real = json.loads((common.PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json').read_text())
    mapped = common.expected_historical_dependencies(real)
    current = Path(real['ssd_support_package']) / 'native/core_reproduction_v1/core_h2_insomnia.execution_receipt.json'
    observed = json.loads(current.read_text())
    assert mapped == observed['dependency_sha256_before'] == observed['dependency_sha256_after']
    return {'synthetic_historical_receipts': 193, 'exact_relocation_pass': True, 'tampered_output_rejected': True,
            'wrong_monitor_plan_rejected': True, 'existing_core_dependency_map_exactly_matches': True,
            'real_historical_campaign_complete_claimed': False}


def monitor_fixtures(root):
    executor = load_executor()
    executor.PLAN = root / 'mock_monitor_plan.json'
    executor.ADMISSION = root / 'mock_monitor_admission.json'
    executor.PLAN.write_text('{}')
    outcomes = []
    for case in ('success', 'helper_exception', 'residual_children', 'post_teardown_floor', 'post_teardown_output', 'post_teardown_admission'):
        output_root = root / case
        output_root.mkdir()
        out = output_root / 'worker'
        out.mkdir()
        executor.ADMISSION.write_text('{}')
        admission_sha = common.sha(executor.ADMISSION)
        plan = {'ssd_output_root': str(output_root), 'internal_floor_bytes': 100, 'ssd_floor_bytes': 100,
                'worker_rss_limit_bytes': 1000, 'output_limit_bytes': 1000000, 'per_worker_seconds_limit': 10, 'runtime_poll_seconds': 2}
        calls = {'teardown': 0, 'snapshot': 0, 'size': 0}
        class FakeProcess:
            pid, returncode = 888888, 0
            def poll(self): return self.returncode
            def wait(self, timeout=None): return self.returncode
        proc = FakeProcess()
        def teardown(process):
            calls['teardown'] += 1
            if case == 'helper_exception' and calls['teardown'] == 1:
                raise RuntimeError('FIXTURE_HELPER_TEARDOWN_FAILURE')
            if case == 'post_teardown_admission': executor.ADMISSION.write_text('{"changed":true}')
            return {'initial_group_members': [{'pid': 999999, 'state': 'S'}] if case == 'residual_children' else [], 'signals': [], 'remaining_group_members': []}
        helper = SimpleNamespace(terminate_owned=teardown, owned_rss=lambda pid: (0, []))
        def snapshot(path):
            calls['snapshot'] += 1
            return {'internal_free_bytes': 50 if case == 'post_teardown_floor' and calls['snapshot'] >= 3 else 10000,
                    'ssd_free_bytes': 10000}
        def size(path):
            calls['size'] += 1
            return 1000001 if case == 'post_teardown_output' and calls['size'] >= 2 else 0
        executor.monitor_helpers = lambda: helper
        executor.snapshot, executor.output_bytes = snapshot, size
        executor.owned_group_exists = lambda pid: False
        executor.subprocess = SimpleNamespace(Popen=lambda *args, **kwargs: proc, STDOUT=-2)
        ownership = []
        failed = False
        try:
            executor.monitored_worker(plan, admission_sha, case, ['MOCK_COMMAND_NO_PROCESS'], out, 0, ownership)
        except RuntimeError:
            failed = True
        receipt = json.loads((out / 'monitor_receipt.json').read_text())
        assert ownership == [] and receipt['owned_group_disappearance_verified']
        assert failed == (case != 'success'), (case, receipt)
        if case == 'helper_exception':
            assert receipt['stop_reason'] == 'OWNED_GROUP_TEARDOWN_EXCEPTION'
            assert receipt['first_teardown_attempt']['exception']['repr'] == "RuntimeError('FIXTURE_HELPER_TEARDOWN_FAILURE')"
            assert (out / 'monitor_failure_receipt.json').exists() and (out / 'owned_group_cleanup_recovery.json').exists()
        if case == 'residual_children': assert receipt['stop_reason'] == 'RESIDUAL_OWNED_GROUP_AFTER_PARENT_EXIT'
        if case == 'post_teardown_floor': assert receipt['final_guard_violations'] == ['INTERNAL_RESOURCE_FLOOR']
        if case == 'post_teardown_output': assert receipt['final_guard_violations'] == ['OUTPUT_LIMIT']
        if case == 'post_teardown_admission': assert receipt['final_guard_violations'] == ['ADMISSION_CHANGED_DURING_RUN']
        outcomes.append({'case': case, 'correct_success_or_failure': True, 'monitor_receipt_written': True, 'actual_processes_launched': 0})
    return outcomes


def lock_fixture(root):
    original_lock, original_time = common.SHARED_LOCK, common.time
    common.SHARED_LOCK = root / 'mock_shared.lock'
    common.time = SimpleNamespace(sleep=lambda seconds: None)
    attempts = []
    def cleanup():
        with common.SHARED_LOCK.open('a+') as competing:
            try: fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError: pass
            else: raise AssertionError('Shared lock released before cleanup')
        attempts.append('lock_still_held')
        if len(attempts) == 1: raise RuntimeError('FIXTURE_CLEANUP_RETRY')
    try:
        with common.exclusive_heavy_lock(before_release=cleanup): pass
        assert len(attempts) == 2
        with common.SHARED_LOCK.open('a+') as competing:
            fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    finally:
        common.SHARED_LOCK, common.time = original_lock, original_time
    return {'retained_through_cleanup_exception': True, 'released_after_cleanup_confirmation': True, 'mock_cleanup_attempts': 2}


def main():
    with tempfile.TemporaryDirectory(prefix='canonical_operational_mock_') as temporary:
        root = Path(temporary)
        result = {'method': method_fixture(), 'history': history_fixture(root), 'monitor_cases': monitor_fixtures(root), 'shared_lock': lock_fixture(root)}
    result.update({'schema': 'canonical_operational_correction_fixture_v1', 'recorded_utc': common.utc(), 'fixture_sha256': common.sha(__file__),
                   'status': 'RESULT_FREE_OPERATIONAL_MOCKS_PASS_NOT_NATIVE_CALIBRATION', 'native_fits_launched': 0, 'GWAS_rows_read': False})
    common.write_new(common.PACKAGE / 'statistical_validation/canonical_operational_correction_fixture_receipt_v4_3.json', result)
    print(json.dumps(result))


if __name__ == '__main__': main()
