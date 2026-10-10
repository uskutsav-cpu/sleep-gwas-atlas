#!/usr/bin/env python3
"""Independent helper-only repaired diagnostic/signal control, real flock."""
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import tempfile
from types import SimpleNamespace

PACKAGE=Path(__file__).resolve().parents[1]
SOURCE=PACKAGE/'scripts/canonical_calibration_common_v4_3.py'
digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
before=digest(SOURCE)
spec=importlib.util.spec_from_file_location('independent_canonical_common_repair',SOURCE)
common=importlib.util.module_from_spec(spec);spec.loader.exec_module(common)
results=[]
with tempfile.TemporaryDirectory(prefix='independent_canonical_repair_',dir='/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/tmp') as folder:
    for index,fault in enumerate((BrokenPipeError('INDEPENDENT_DIAGNOSTIC_FAULT'),KeyboardInterrupt('INDEPENDENT_DIAGNOSTIC_FAULT'))):
        common.SHARED_LOCK=Path(folder)/('lock'+str(index))
        common.time=SimpleNamespace(sleep=lambda seconds:None)
        def diagnostic(*args,_error=fault,**kwargs):raise _error
        common.print=diagnostic
        attempts=[]
        def cleanup():
            with common.SHARED_LOCK.open('a+') as competing:
                try:fcntl.flock(competing.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:pass
                else:raise AssertionError('LOCK_RELEASED_WITH_UNVERIFIED_CLEANUP')
            attempts.append(True)
            if len(attempts)==1:raise RuntimeError('INDEPENDENT_UNVERIFIED_FIRST_CLEANUP')
        with common.exclusive_heavy_lock(before_release=cleanup):pass
        assert len(attempts)==2
        with common.SHARED_LOCK.open('a+') as competing:fcntl.flock(competing.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        results.append({'diagnostic_fault':type(fault).__name__,'cleanup_attempts':2,'lock_held_until_second_cleanup':True,'released_after_certified_callback':True})
    numbers=(signal.SIGTERM,signal.SIGHUP,signal.SIGINT)
    prior={number:signal.getsignal(number) for number in numbers}
    with common.deferred_termination_signals():
        for number in numbers:os.kill(os.getpid(),number)
        assert common.TERMINATION_REQUEST==list(numbers)
        try:common.assert_no_termination()
        except RuntimeError as e:assert 'DEFERRED_TERMINATION_SIGNAL' in str(e)
        else:raise AssertionError('PENDING_SIGNAL_NOT_ACTED_ON')
    assert {number:signal.getsignal(number) for number in numbers}==prior
assert digest(SOURCE)==before
record={'status':'PASS_REPAIRED_DIAGNOSTIC_LOCK_AND_DEFERRED_SIGNAL_HELPERS_ONLY',
        'common_helper_sha256':before,'checker_sha256':digest(__file__),'diagnostic_controls':results,
        'real_self_signals':list(numbers),'signals_deferred_then_protected_assert_failed':True,'signal_handlers_restored':True,
        'native_or_other_workers_launched':0,'estimator_imports':0,'GWAS_reads':0,
        'limit':'Helpers and real temporary flock only; no live native worker or empirical calibration tested.'}
out=PACKAGE/'reviews/independent_canonical_repaired_guardian_receipt_v4_3.json'
with out.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps(record))
