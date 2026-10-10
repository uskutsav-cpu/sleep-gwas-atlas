#!/usr/bin/env python3
"""Extra no-process fixtures for final guards and unresolved-group lock holding."""
import fcntl
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace

import canonical_calibration_common_v4_2 as common


def main():
    source = common.PACKAGE / 'scripts/42_prepare_native_canonical_calibration_v4_2.py'
    spec = importlib.util.spec_from_file_location('canonical_cleanup_mock', source)
    executor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(executor)
    with tempfile.TemporaryDirectory(prefix='canonical_cleanup_mock_') as temporary:
        root = Path(temporary)
        executor.PLAN, executor.ADMISSION = root / 'plan.json', root / 'admission.json'
        executor.PLAN.write_text('{}')
        executor.ADMISSION.write_text('{}')
        plan_hash, admission_hash = common.sha(executor.PLAN), common.sha(executor.ADMISSION)
        plan = {'internal_floor_bytes': 100, 'ssd_floor_bytes': 100, 'output_limit_bytes': 1000, 'per_worker_seconds_limit': 10}
        healthy = {'internal_free_bytes': 1000, 'ssd_free_bytes': 1000}
        assert executor.final_violations(plan, healthy, 0, 11, admission_hash, plan_hash) == ['WORKER_DEADLINE']
        assert executor.final_violations(plan, {'internal_free_bytes': 1000, 'ssd_free_bytes': 50}, 0, 0, admission_hash, plan_hash) == ['SSD_RESOURCE_FLOOR']
        executor.PLAN.write_text('{"changed":true}')
        assert executor.final_violations(plan, healthy, 0, 0, admission_hash, plan_hash) == ['PLAN_CHANGED_DURING_RUN']
        calls = {'helper': 0, 'kills': 0, 'lock_held_while_unresolved': 0}
        alive = {'group': True}
        class FakeProcess:
            pid = 888889
            def wait(self, timeout=None): return 0
        def helper_teardown(proc):
            calls['helper'] += 1
            raise RuntimeError('MOCK_UNRESOLVED_TEARDOWN')
        def kill(pid, sig):
            assert pid == 888889
            calls['kills'] += 1
            if calls['kills'] == 2: alive['group'] = False
        original_lock = common.SHARED_LOCK
        common.SHARED_LOCK = root / 'temporary_shared.lock'
        def pause(seconds):
            assert alive['group']
            with common.SHARED_LOCK.open('a+') as competing:
                try: fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError: calls['lock_held_while_unresolved'] += 1
                else: raise AssertionError('Lock released while mock group alive')
        executor.monitor_helpers = lambda: SimpleNamespace(terminate_owned=helper_teardown)
        executor.owned_group_exists = lambda pid: alive['group']
        executor.os = SimpleNamespace(killpg=kill)
        executor.time = SimpleNamespace(sleep=pause)
        output = root / 'mock_worker'
        output.mkdir()
        entry = {'proc': FakeProcess(), 'output_dir': output, 'cleanup_errors': []}
        ownership = [entry]
        try:
            with common.exclusive_heavy_lock(before_release=lambda: executor.await_owned_cleanup(ownership, None) if ownership else None):
                executor.await_owned_cleanup(ownership, None)
            assert ownership == [] and not alive['group'] and calls == {'helper': 2, 'kills': 2, 'lock_held_while_unresolved': 1}
            with common.SHARED_LOCK.open('a+') as competing:
                fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            proof = json.loads((output / 'owned_group_cleanup_recovery.json').read_text())
            assert proof['group_disappearance_verified_before_shared_lock_release'] and proof['cleanup_errors']
        finally:
            common.SHARED_LOCK = original_lock
    record = {'schema': 'canonical_cleanup_and_final_guard_fixture_v1', 'recorded_utc': common.utc(), 'fixture_sha256': common.sha(__file__),
              'deadline_failure_pass': True, 'SSD_floor_failure_pass': True, 'plan_change_failure_pass': True,
              'unresolved_group_recovery_pass': True, 'mock_counts': calls,
              'lock_retained_until_mock_group_disappearance': True, 'actual_processes_launched': 0, 'native_fits_launched': 0,
              'status': 'OPERATIONAL_MOCK_PASS_ONLY_NO_NATIVE_CALIBRATION'}
    common.write_new(common.PACKAGE / 'statistical_validation/canonical_cleanup_and_final_guard_fixture_receipt_v4_2.json', record)
    print(json.dumps(record))


if __name__ == '__main__': main()
