"""Small identity/serialization helpers; no numerical library or estimator import."""
import hashlib
import json
import fcntl
import os
import sys
import time
import signal
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
PLAN = PACKAGE / 'manifests/native_canonical_calibration_plan_v4_3.json'
ADMISSION = PACKAGE / 'manifests/native_canonical_calibration_admission_v4_3.json'
SHARED_LOCK = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/native_heavy_worker.lock')
TERMINATION_REQUEST = []


def catchable_termination(signum, frame):
    # Do not raise between process creation and ownership registration, or
    # inside cleanup/persistence. The protected monitor interprets requests.
    if len(TERMINATION_REQUEST) < 64:
        TERMINATION_REQUEST.append(signum)


@contextmanager
def deferred_termination_signals():
    numbers = (signal.SIGTERM, signal.SIGHUP, signal.SIGINT)
    previous = {number: signal.getsignal(number) for number in numbers}
    try:
        for number in numbers:
            signal.signal(number, catchable_termination)
        yield
    finally:
        for number, handler in previous.items():
            signal.signal(number, handler)


def assert_no_termination():
    if TERMINATION_REQUEST:
        raise RuntimeError('DEFERRED_TERMINATION_SIGNAL: ' + str(TERMINATION_REQUEST))


def safe_diagnostic(value_factory):
    # Formatting and emission BOTH belong inside the non-throwing boundary.
    try:
        value = value_factory() if callable(value_factory) else value_factory
        print(value, file=sys.stderr, flush=True)
    except BaseException:
        pass


@contextmanager
def exclusive_heavy_lock(*, before_release=None):
    """Shared by all native sensitivity/calibration/raw-pipeline executors."""
    with SHARED_LOCK.open('a+') as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('SHARED_HEAVY_WORKER_LOCK_BUSY_NO_WORKER_LAUNCHED') from error
        try:
            yield handle.fileno()
        finally:
            # A teardown exception may never release the mutex around a live
            # group. Retry the cleanup callback while retaining the open lock.
            # This also catches interruption of cleanup itself. SIGKILL of the
            # supervisor cannot be handled by Python and needs host review.
            if before_release is not None:
                while True:
                    try:
                        before_release()
                        break
                    except BaseException as error:
                        safe_diagnostic(lambda: 'SHARED_HEAVY_WORKER_LOCK_RETAINED: ' + repr(error))
                        try:
                            time.sleep(2)
                        except BaseException:
                            pass
            # Closing after verified cleanup releases our description. Avoid an
            # explicit LOCK_UN that would also unlock inherited descriptions.


def validate_inherited_lock(fd):
    stat, current = os.fstat(fd), SHARED_LOCK.stat()
    if (stat.st_dev, stat.st_ino) != (current.st_dev, current.st_ino):
        raise RuntimeError('INHERITED_HEAVY_LOCK_IDENTITY_MISMATCH')
    # Re-acquiring flock on the inherited open file description preserves the
    # parent's ownership; an unrelated simultaneous worker cannot acquire it.
    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(record, f, indent=2, allow_nan=False)
        f.write('\n')


def check_hashes(expected):
    actual = {p: sha(p) for p in expected}
    if actual != expected:
        raise RuntimeError('BOUND_IDENTITY_CHANGED: ' + ', '.join(p for p in expected if actual[p] != expected[p]))
    return actual


def clean(value, np):
    if isinstance(value, dict):
        return {str(k): clean(v, np) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v, np) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist(), np)
    if isinstance(value, np.generic):
        return clean(value.item(), np)
    if isinstance(value, float) and not np.isfinite(value):
        return {'nonfinite': str(value)}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError('UNSUPPORTED_CAPTURE_TYPE: ' + str(type(value)))


def expected_historical_dependencies(historical):
    """One exact, source-bound PACKAGE relocation; no hash/basename matching."""
    expected = dict(historical['dependencies_sha256'])
    original = str(ROOT / 'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    relocated = str(Path(historical['ssd_support_package']) / 'scripts/native_ldsc_capture.py')
    if original not in expected or relocated in expected:
        raise RuntimeError('HISTORICAL_CAPTURE_RELOCATION_SCHEMA_CHANGED')
    value = expected.pop(original)
    if historical['support_file_sha256'].get(relocated) != value or sha(original) != value or sha(relocated) != value:
        raise RuntimeError('HISTORICAL_CAPTURE_RELOCATION_IDENTITY_CHANGED')
    expected[relocated] = value
    return expected


def historical_completion(plan):
    """Require all 190 original jobs and three completed, reaped stage monitors."""
    historical_path = Path(plan['historical_plan_path'])
    if sha(historical_path) != plan['historical_plan_sha256']:
        raise RuntimeError('HISTORICAL_PLAN_IDENTITY_CHANGED')
    historical = json.loads(historical_path.read_text())
    if len(historical['jobs']) != 190:
        raise RuntimeError('HISTORICAL_CARDINALITY_CHANGED')
    expected_dependencies = expected_historical_dependencies(historical)
    receipts = {}
    identities = {r['path']: r['actual_sha256'] for r in historical['inputs_verified'] if r['match']}
    for job in historical['jobs']:
        p = Path(historical['ssd_support_package']) / 'native' / (job['stage'] + '_reproduction_v1') / (job['job_id'] + '.execution_receipt.json')
        r = json.loads(p.read_text())
        if r.get('job') != job or r.get('returncode') != 0 or r.get('scientific_cardinality_gate_pass') is not True or r.get('execution_identity_gate_pass') is not True:
            raise RuntimeError('HISTORICAL_JOB_NOT_COMPLETED: ' + job['job_id'])
        if r.get('dependency_sha256_before') != expected_dependencies or r.get('dependency_sha256_after') != expected_dependencies:
            raise RuntimeError('HISTORICAL_DEPENDENCY_RECEIPT_MISMATCH: ' + job['job_id'])
        expected_inputs = {p: identities[p] for p in job['inputs']}
        if r.get('input_sha256') != expected_inputs or r.get('input_sha256_after') != expected_inputs:
            raise RuntimeError('HISTORICAL_INPUT_RECEIPT_MISMATCH: ' + job['job_id'])
        outputs = r.get('all_output_sha256')
        output_dir = Path(historical['ssd_support_package']) / 'native' / (job['stage'] + '_reproduction_v1')
        full_precision = str(output_dir / (job['job_id'] + '.full_precision.json'))
        if not isinstance(outputs, dict) or not outputs or outputs.get(full_precision) != r.get('output_sha256'):
            raise RuntimeError('HISTORICAL_FULL_PRECISION_OUTPUT_BINDING_MISSING: ' + job['job_id'])
        for output in outputs:
            op = Path(output)
            if op.parent.resolve() != output_dir.resolve() or not op.name.startswith(job['job_id']):
                raise RuntimeError('HISTORICAL_OUTPUT_SCOPE_MISMATCH: ' + job['job_id'])
        check_hashes(outputs)
        receipts[str(p)] = sha(p)
    for stage in ('core', 'extension', 'validation'):
        p = PACKAGE / 'logs' / (stage + '_native_monitor_receipt_v4.json')
        r = json.loads(p.read_text())
        if r.get('returncode') != 0 or r.get('stop_reason') is not None or r.get('process_group_teardown', {}).get('remaining_group_members') != []:
            raise RuntimeError('HISTORICAL_STAGE_NOT_COMPLETED_AND_REAPED: ' + stage)
        if r.get('stage') != stage or r.get('plan_sha256') != plan['historical_plan_sha256']:
            raise RuntimeError('HISTORICAL_STAGE_MONITOR_PLAN_BINDING_MISMATCH: ' + stage)
        receipts[str(p)] = sha(p)
    return receipts


def admit(plan_path=PLAN, admission_path=ADMISSION):
    """Separate review must explicitly bind exact prepared code/plan and history."""
    plan_path, admission_path = Path(plan_path), Path(admission_path)
    plan = json.loads(plan_path.read_text())
    a = json.loads(admission_path.read_text())
    if a.get('execution_admitted') is not True or a.get('scope') != 'METHOD_CALIBRATION_ONLY':
        raise RuntimeError('NATIVE_CALIBRATION_EXECUTION_NOT_ADMITTED')
    if a.get('allow_41_covariance_outcomes') is not False or a.get('allow_calibrated_biological_p_values') is not False:
        raise RuntimeError('BIOLOGICAL_ANALYSIS_SCOPE_NOT_ADMITTED')
    if plan.get('shared_lock_path') != str(SHARED_LOCK):
        raise RuntimeError('SHARED_HEAVY_WORKER_LOCK_BINDING_REQUIRED')
    if a.get('plan_sha256') != sha(plan_path) or a.get('independent_binding_review_pass') is not True or a.get('resource_plan_review_pass') is not True:
        raise RuntimeError('INDEPENDENT_REVIEW_BINDING_REQUIRED')
    review_path = Path(a['independent_review_path'])
    if sha(review_path) != a.get('independent_review_sha256'):
        raise RuntimeError('INDEPENDENT_REVIEW_CHANGED')
    if a.get('historical_campaign_receipt_sha256') != historical_completion(plan):
        raise RuntimeError('HISTORICAL_COMPLETION_NOT_BOUND_BY_REVIEW')
    check_hashes(plan['dependencies_sha256'])
    return plan, a
