#!/usr/bin/env python3
"""Catchable-signal/diagnostic/post-hash fault fixtures; no estimator or worker."""
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import signal
import tempfile
from types import SimpleNamespace

import canonical_calibration_common_v4_3 as common


def executor_module():
    path = common.PACKAGE / 'scripts/42_prepare_native_canonical_calibration_v4_3.py'
    spec = importlib.util.spec_from_file_location('canonical_signal_mock', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FaultStream:
    def __init__(self, error): self.error = error
    def write(self, value): raise self.error
    def flush(self): raise self.error


def diagnostic_lock_cases(root):
    old_lock, old_time, old_stderr = common.SHARED_LOCK, common.time, common.sys.stderr
    result = []
    try:
        common.time = SimpleNamespace(sleep=lambda seconds: None)
        for number, error in enumerate((BrokenPipeError('MOCK_BROKEN_PIPE'), KeyboardInterrupt('MOCK_DIAGNOSTIC_INTERRUPT'))):
            common.SHARED_LOCK = root / ('fault_lock_' + str(number))
            common.sys.stderr = FaultStream(error)
            attempts = []
            def cleanup():
                with common.SHARED_LOCK.open('a+') as competing:
                    try: fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    except BlockingIOError: pass
                    else: raise AssertionError('Diagnostic failure released lock')
                attempts.append(True)
                if len(attempts) == 1: raise RuntimeError('MOCK_UNVERIFIED_CLEANUP')
            with common.exclusive_heavy_lock(before_release=cleanup): pass
            assert len(attempts) == 2
            # Formatting failures belong inside the same safe boundary.
            common.safe_diagnostic(lambda: (_ for _ in ()).throw(error))
            result.append({'diagnostic_fault': type(error).__name__, 'lock_retained_until_second_cleanup': True, 'formatting_fault_contained': True})
    finally:
        common.SHARED_LOCK, common.time, common.sys.stderr = old_lock, old_time, old_stderr
    return result


def real_deferred_signals(root):
    old_lock = common.SHARED_LOCK
    common.SHARED_LOCK = root / 'real_signal_fixture.lock'
    common.TERMINATION_REQUEST.clear()
    numbers = (signal.SIGTERM, signal.SIGHUP, signal.SIGINT)
    previous = {number: signal.getsignal(number) for number in numbers}
    try:
        with common.deferred_termination_signals():
            with common.exclusive_heavy_lock():
                for number in numbers:
                    os.kill(os.getpid(), number)
                assert common.TERMINATION_REQUEST == list(numbers)
        assert {number: signal.getsignal(number) for number in numbers} == previous
        return {'real_signals_to_fixture_process': list(numbers), 'handler_requests_deferred_without_raising': True, 'original_handlers_restored': True}
    finally:
        common.SHARED_LOCK = old_lock
        common.TERMINATION_REQUEST.clear()


def guardian_diagnostic_case(root):
    executor = executor_module()
    old_lock, old_stderr = common.SHARED_LOCK, common.sys.stderr
    common.SHARED_LOCK = root / 'guardian_diagnostic.lock'
    common.sys.stderr = FaultStream(KeyboardInterrupt('MOCK_GUARDIAN_DIAGNOSTIC'))
    alive = {'group': True}
    count = {'kills': 0, 'lock_observations': 0}
    class Process:
        pid = 888890
        def wait(self, timeout=None): return 0
    def helper(proc): raise RuntimeError('MOCK_TEARDOWN_HELPER')
    def kill(pid, sig):
        count['kills'] += 1
        if count['kills'] == 2: alive['group'] = False
    def pause(seconds):
        with common.SHARED_LOCK.open('a+') as competing:
            try: fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError: count['lock_observations'] += 1
            else: raise AssertionError('Guardian diagnostic released lock')
    executor.monitor_helpers = lambda: SimpleNamespace(terminate_owned=helper)
    executor.os, executor.time = SimpleNamespace(killpg=kill), SimpleNamespace(sleep=pause)
    executor.owned_group_exists = lambda pid: alive['group']
    output = root / 'guardian_worker'
    output.mkdir()
    ownership = [{'proc': Process(), 'output_dir': output, 'cleanup_errors': []}]
    try:
        with common.exclusive_heavy_lock(before_release=lambda: executor.await_owned_cleanup(ownership, None) if ownership else None):
            executor.await_owned_cleanup(ownership, None)
        assert not ownership and not alive['group'] and count == {'kills': 2, 'lock_observations': 1}
        return {'guardian_diagnostic_interrupt_contained': True, 'lock_retained_while_mock_group_alive': True, 'counts': count}
    finally:
        common.SHARED_LOCK, common.sys.stderr = old_lock, old_stderr


def registration_and_hash_cases(root):
    results = []
    for case in ('signal_at_registration', 'floor_after_hash', 'deadline_after_hash', 'admission_after_hash', 'output_after_hash'):
        executor = executor_module()
        folder = root / case
        folder.mkdir()
        out = folder / 'worker'
        out.mkdir()
        executor.PLAN, executor.ADMISSION = folder / 'plan.json', folder / 'admission.json'
        executor.PLAN.write_text('{}')
        executor.ADMISSION.write_text('{}')
        admission_sha = common.sha(executor.ADMISSION)
        state = {'hash_started': False, 'clock': 0.}
        class Process:
            pid, returncode = 888891, 0
            def poll(self): return 0
            def wait(self, timeout=None): return 0
        def create(*args, **kwargs):
            if case == 'signal_at_registration': common.catchable_termination(signal.SIGTERM, None)
            return Process()
        helper = SimpleNamespace(terminate_owned=lambda proc: {'initial_group_members': [], 'signals': [], 'remaining_group_members': []},
                                 owned_rss=lambda pid: (0, []))
        def hashing(path):
            path = Path(path)
            if path.parent == out:
                state['hash_started'] = True
                if case == 'deadline_after_hash': state['clock'] = 11.
                if case == 'admission_after_hash': executor.ADMISSION.write_text('{"changed":true}')
            return common.sha(path)
        executor.sha = hashing
        executor.time = SimpleNamespace(monotonic=lambda: state['clock'], sleep=lambda seconds: None)
        executor.snapshot = lambda path: {'internal_free_bytes': 1000, 'ssd_free_bytes': 50 if case == 'floor_after_hash' and state['hash_started'] else 1000}
        executor.output_bytes = lambda path: 1001 if case == 'output_after_hash' and state['hash_started'] else 0
        executor.monitor_helpers = lambda: helper
        executor.owned_group_exists = lambda pid: False
        executor.subprocess = SimpleNamespace(Popen=create, STDOUT=-2)
        plan = {'ssd_output_root': str(folder), 'internal_floor_bytes': 100, 'ssd_floor_bytes': 100, 'worker_rss_limit_bytes': 1000,
                'output_limit_bytes': 1000, 'per_worker_seconds_limit': 10, 'runtime_poll_seconds': 2}
        common.TERMINATION_REQUEST.clear()
        ownership = []
        try:
            executor.monitored_worker(plan, admission_sha, case, ['MOCK_NO_PROCESS'], out, 0, ownership)
        except RuntimeError: pass
        else: raise AssertionError('Fault unexpectedly admitted success: ' + case)
        receipt = json.loads((out / 'monitor_receipt.json').read_text())
        assert not ownership and receipt['bulk_output_hashes_completed_before_final_guards']
        if case == 'signal_at_registration':
            assert receipt['deferred_termination_requests'] == [signal.SIGTERM]
            assert receipt['stop_reason'].startswith('DEFERRED_TERMINATION_SIGNAL')
            assert receipt['owned_group_disappearance_verified']
        else:
            expected = {'floor_after_hash': 'SSD_RESOURCE_FLOOR', 'deadline_after_hash': 'WORKER_DEADLINE',
                        'admission_after_hash': 'ADMISSION_CHANGED_DURING_RUN', 'output_after_hash': 'OUTPUT_LIMIT'}[case]
            assert receipt['stop_reason'] == expected and expected in receipt['final_guard_violations']
        results.append({'case': case, 'failure_receipt_written': True, 'failure_gate_pass': True, 'actual_processes_launched': 0})
        common.TERMINATION_REQUEST.clear()
    return results


def main():
    with tempfile.TemporaryDirectory(prefix='canonical_signal_mock_') as temporary:
        root = Path(temporary)
        result = {'diagnostic_lock_cases': diagnostic_lock_cases(root), 'signals': real_deferred_signals(root),
                  'guardian': guardian_diagnostic_case(root), 'registration_and_hash_cases': registration_and_hash_cases(root)}
    result.update({'schema': 'canonical_signal_diagnostic_post_hash_fixture_v1', 'recorded_utc': common.utc(), 'fixture_sha256': common.sha(__file__),
                   'status': 'SIGNAL_DIAGNOSTIC_AND_POST_HASH_FAULT_CONTROLS_PASS_NO_NATIVE_CALIBRATION', 'native_fits_launched': 0,
                   'estimator_modules_imported': False, 'GWAS_rows_read': False})
    common.write_new(common.PACKAGE / 'statistical_validation/canonical_signal_diagnostic_post_hash_fixture_receipt_v4_3.json', result)
    print(json.dumps(result))


if __name__ == '__main__': main()
