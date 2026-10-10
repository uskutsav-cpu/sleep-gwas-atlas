#!/usr/bin/env python3
"""Real mutex contention and owned-group failure injection; no network/estimators."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types

PACKAGE = Path(__file__).resolve().parents[1]
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
source = PACKAGE / 'scripts/52_acquire_extension_raw_sources_v4.py'
spec = importlib.util.spec_from_file_location('acquisition_fault_control', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
checks = []
before_hash = module.hashes(source)[1]
plan_hash = module.hashes(module.PLAN)[1]

with tempfile.TemporaryDirectory(prefix='acquisition_fault_controls_', dir=SSD/'tmp') as temporary:
    folder = Path(temporary)
    lock = folder / 'exclusive.lock'
    child_code = "import fcntl,sys; f=open(sys.argv[1],'a+b'); fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)"
    with module.family_lock(lock):
        result = subprocess.run([sys.executable, '-B', '-c', child_code, str(lock)], capture_output=True, text=True)
        assert result.returncode != 0 and 'BlockingIOError' in result.stderr
        checks.append({'control': 'real_second_process_mutex_contention', 'pass': True})
    result = subprocess.run([sys.executable, '-B', '-c', child_code, str(lock)], capture_output=True, text=True)
    assert result.returncode == 0
    checks.append({'control': 'mutex_release_after_owned_scope', 'pass': True})
    proc = subprocess.Popen([sys.executable, '-B', '-c', 'import time;time.sleep(60)'], start_new_session=True)
    assert os.getpgid(proc.pid) == proc.pid
    def fail_cleanup(proc):
        raise RuntimeError('INJECTED_SHARED_TEARDOWN_FAILURE')
    monitor = types.SimpleNamespace(terminate_owned=fail_cleanup)
    cleanup = module.cleanup_proof(monitor, proc)
    assert cleanup['teardown_verified'] and 'INJECTED_SHARED_TEARDOWN_FAILURE' in cleanup['cleanup_error']
    assert proc.returncode == -9 and cleanup['remaining_group_members'] == []
    checks.append({'control': 'cleanup_exception_independent_owned_pgid_kill_and_reap', 'pass': True, 'proof': cleanup})
    immutable = folder / 'immutable.json'
    module.write_new(immutable, {'original': True})
    try:
        module.write_new(immutable, {'overwrite': True})
        raise AssertionError('IMMUTABLE_WRITE_ACCEPTED')
    except FileExistsError:
        assert json.loads(immutable.read_text()) == {'original': True}
    checks.append({'control': 'exclusive_receipt_write_refuses_overwrite', 'pass': True})
    # Invoke the real source-finally path with a failed worker and failed shared
    # cleanup helper. There is no HTTP request: /usr/bin/false is the worker.
    module.FOLDER = folder / 'attempt'
    module.PACKAGE = folder / 'package'
    module.CURL = Path('/usr/bin/false')
    module.assert_bindings = lambda *args: None
    for child in ('raw', 'logs', 'receipts'):
        (module.FOLDER / child).mkdir(parents=True)
    state = {'internal_free_bytes': 4*1024**3, 'ssd_free_bytes': 10*1024**3}
    monitor = types.SimpleNamespace(snapshot=lambda: dict(state), terminate_owned=fail_cleanup,
                                    owned_rss=lambda pid:(0,[pid]), final_limits=lambda *args:None,
                                    group_members=lambda pid:[])
    plan = {'explicit_resume': None, 'internal_floor_bytes': 3*1024**3, 'ssd_floor_bytes': 5*1024**3,
            'per_body_seconds_limit': 10, 'family_seconds_limit': 20, 'runtime_poll_seconds': .01}
    member = {'index': 2, 'extension_trait_id': 'fault_fixture', 'body_path': str(module.FOLDER/'raw/body.bin'),
              'expected_bytes': 10, 'expected_md5': '0'*32, 'expected_sha256': '0'*64,
              'url': 'https://invalid.example/never-requested'}
    import time
    try:
        module.acquire(plan, member, monitor, '0'*64, time.monotonic())
        raise AssertionError('FAILED_WORKER_WAS_ADMITTED')
    except RuntimeError as exc:
        assert 'SOURCE_STOP_PRESERVED' in str(exc)
    source_receipt = module.FOLDER / 'receipts/fault_fixture.json'
    receipt = json.loads(source_receipt.read_text())
    assert receipt['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
    assert 'PROCESS_TEARDOWN_FAILURE' in receipt['stop_reason']
    assert receipt['teardown']['teardown_verified']
    assert not Path(member['body_path']).exists()
    copied = module.PACKAGE / 'source_provenance/extension_raw_acquisition_v4_4/fault_fixture.json'
    assert source_receipt.read_bytes() == copied.read_bytes()
    checks.append({'control': 'worker_failure_plus_cleanup_exception_always_preserves_failure_receipts', 'pass': True,
                   'worker_returncode': receipt['returncode'], 'receipt_sha256': module.hashes(source_receipt)[1]})

    import signal
    old_handlers = {number:signal.getsignal(number) for number in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT)}
    try:
        for number in old_handlers:
            signal.signal(number,module.catchable_termination)
            os.kill(os.getpid(),number)
        assert module.TERMINATION_REQUEST == list(old_handlers)
        checks.append({'control':'real_SIGTERM_SIGHUP_SIGINT_handlers_defer_without_raising','pass':True})
    finally:
        for number,handler in old_handlers.items():signal.signal(number,handler)
        module.TERMINATION_REQUEST.clear()
    original_cleanup = module.cleanup_proof
    original_write = module.write_new
    original_retain = module.retain_lock_until_gone
    quarantine_calls=[]
    def forced_unverified(*args):
        return {'teardown_verified':False,'cleanup_error':'INJECTED_UNVERIFIED_TEARDOWN'}
    def persistence_failure(path,value):
        if Path(path).name=='fault_fixture.json':raise OSError('INJECTED_RECEIPT_PERSISTENCE_FAILURE')
        return original_write(path,value)
    def checked_retain(proc,receipt):
        competitor=subprocess.run([sys.executable,'-B','-c',child_code,str(lock)],capture_output=True,text=True)
        assert competitor.returncode!=0 and 'BlockingIOError' in competitor.stderr
        quarantine_calls.append(proc.pid)
        return original_retain(proc,receipt)
    module.cleanup_proof=forced_unverified
    module.write_new=persistence_failure
    module.retain_lock_until_gone=checked_retain
    module.FOLDER=folder/'persistence_failure_attempt'
    for child in ('raw','logs','receipts'):(module.FOLDER/child).mkdir(parents=True)
    member={**member,'body_path':str(module.FOLDER/'raw/body.bin')}
    with module.family_lock(lock):
        try:
            module.acquire(plan,member,monitor,'0'*64,time.monotonic())
            raise AssertionError('PERSISTENCE_FAILURE_WAS_ACCEPTED')
        except OSError as exc:
            assert 'INJECTED_RECEIPT_PERSISTENCE_FAILURE' in str(exc)
    assert quarantine_calls
    checks.append({'control':'receipt_write_failure_cannot_bypass_owned_quarantine_and_flock','pass':True,
                   'quarantine_calls':len(quarantine_calls)})
    module.write_new=original_write
    quarantine_calls.clear()
    module.FOLDER=folder/'status_failure_attempt'
    for child in ('raw','logs','receipts'):(module.FOLDER/child).mkdir(parents=True)
    member={**member,'body_path':str(module.FOLDER/'raw/body.bin')}
    def status_failure(*args):raise OSError('INJECTED_PRE_PERSISTENCE_STATUS_FAILURE')
    module.safe_state=status_failure
    with module.family_lock(lock):
        try:
            module.acquire(plan,member,monitor,'0'*64,time.monotonic())
            raise AssertionError('STATUS_FAILURE_WAS_ACCEPTED')
        except OSError as exc:
            assert 'INJECTED_PRE_PERSISTENCE_STATUS_FAILURE' in str(exc)
    assert quarantine_calls
    checks.append({'control':'early_cleanup_status_failure_still_enforces_outer_owned_quarantine','pass':True,
                   'quarantine_calls':len(quarantine_calls)})

assert module.hashes(source)[1] == before_hash
assert module.hashes(PACKAGE/'manifests/extension_raw_acquisition_plan_v4_4.json')[1] == plan_hash
output = PACKAGE / 'logs/extension_acquisition_fault_controls_v4_4.json'
module.write_new(output, {'completed_utc': datetime.now(timezone.utc).isoformat(), 'executor_sha256': before_hash,
                          'plan_sha256': plan_hash, 'checks': checks, 'all_pass': all(c['pass'] for c in checks),
                          'network_requests': 0, 'native_estimators_launched': 0,
                          'test_script_sha256': module.hashes(Path(__file__))[1]})
print(json.dumps({'all_pass': True, 'controls': len(checks), 'receipt': str(output)}))
