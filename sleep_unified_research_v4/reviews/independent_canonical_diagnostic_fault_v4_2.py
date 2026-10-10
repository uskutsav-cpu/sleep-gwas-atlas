#!/usr/bin/env python3
"""Bound helper-only fault demonstration: real flock, no workers/estimators."""
import fcntl
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile

PACKAGE = Path(__file__).resolve().parents[1]
SOURCE = PACKAGE / 'scripts/canonical_calibration_common_v4_2.py'
before = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
spec = importlib.util.spec_from_file_location('reviewed_common_only', SOURCE)
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
attempts = []

def cleanup():
    attempts.append('unverified_cleanup')
    raise RuntimeError('REVIEW_INJECTED_UNVERIFIED_CLEANUP')

def broken_diagnostic(*args, **kwargs):
    raise BrokenPipeError('REVIEW_INJECTED_STDERR_FAILURE')

with tempfile.TemporaryDirectory(prefix='independent_canonical_diagnostic_', dir='/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/tmp') as folder:
    common.SHARED_LOCK = Path(folder) / 'owned.lock'
    common.print = broken_diagnostic
    escaped = False
    try:
        with common.exclusive_heavy_lock(before_release=cleanup):
            pass
    except BrokenPipeError:
        escaped = True
    assert escaped and len(attempts) == 1
    with common.SHARED_LOCK.open('a+') as competitor:
        fcntl.flock(competitor.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        reacquired_after_unverified_cleanup = True
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == before
receipt = {'status':'CONFIRMED_SUPERVISOR_LOCK_RETENTION_FAILURE_ON_DIAGNOSTIC_ERROR',
           'helper_sha256':before, 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'cleanup_remained_unverified':True, 'cleanup_attempts':len(attempts),
           'diagnostic_exception_escaped':escaped,
           'separate_file_description_reacquired_after_unverified_cleanup':reacquired_after_unverified_cleanup,
           'owned_workers_launched':0, 'native_estimators_launched':0,
           'limit':'No worker descriptor exists in this helper-only fixture. A real inherited worker description may keep the mutex held; supervisor cleanup/retry nonetheless escaped.'}
out = PACKAGE / 'reviews/independent_canonical_diagnostic_fault_receipt_v4_2.json'
with out.open('x') as handle:
    json.dump(receipt, handle, indent=2); handle.write('\n')
print(json.dumps(receipt))
