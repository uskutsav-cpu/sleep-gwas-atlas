"""Explicit stage receipt selection after the preserved extension floor stop.

No estimator import, worker launch or source-body read. The caller records each
consumed metadata identity, and still verifies every scientific job separately.
"""
import hashlib
import json
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
HISTORICAL = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
HISTORICAL_SHA = 'f555dc441e8c93c528e153d2e8689edabb29ba17e4d8b0da68402ac4e4c98e88'
CONTINUATION = PACKAGE / 'manifests/native_extension_checkpoint_continuation_v4_3.json'
CONTINUATION_SHA = 'f087b74e1ebaf6b543513d69948d8dbb2232e5190360dd6fce3acc16e68ed5b3'
ADMISSION = PACKAGE / 'manifests/native_extension_checkpoint_continuation_admission_v4_3.json'
ADMISSION_SHA = '0d50eb33161c12bb7e6801bb7310c9a7e26fa75d78c842df6001e10bf3d11135'
FAILED = PACKAGE / 'logs/extension_native_monitor_receipt_v4.json'
FAILED_SHA = '0e05b99fc6913168214e57abf7c4ba5ee2aa9b53e3af360b2114a1f4b3c3c264'
CONTROLLER_SHA = 'b1309c0bd97736d7282c636956b859290dbe5430a5997137d8be0776e3160283'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def stage_monitor_path(stage):
    if stage not in ('core', 'extension', 'validation'):
        raise ValueError('EXACT_HISTORICAL_STAGE_REQUIRED')
    suffix = '_v4_3.json' if stage == 'extension' else '_v4.json'
    return PACKAGE / 'logs' / (stage + '_native_monitor_receipt' + suffix)


def stage_monitor(stage, record=None):
    consumed = {}

    def read(path, expected=None):
        path = Path(path)
        if path.is_symlink() or not path.is_file():
            raise RuntimeError('STAGE_EVIDENCE_NOT_REGULAR_FILE: ' + str(path))
        actual = sha(path)
        if expected is not None and actual != expected:
            raise RuntimeError('STAGE_EVIDENCE_IDENTITY_CHANGED: ' + str(path))
        consumed[str(path)] = actual
        if record is not None:
            record(path)
        return json.loads(path.read_text()) if path.suffix == '.json' else None

    historical = read(HISTORICAL, HISTORICAL_SHA)
    path = stage_monitor_path(stage)
    if Path(str(path) + '.failure.json').exists():
        raise RuntimeError('TERMINAL_FAILURE_INVALIDATES_PROVISIONAL_STAGE_MONITOR')
    monitor = read(path)
    if monitor.get('stage') != stage or monitor.get('plan_sha256') != HISTORICAL_SHA:
        raise RuntimeError('EXACT_STAGE_AND_HISTORICAL_PLAN_REQUIRED')
    if monitor.get('returncode') != 0 or monitor.get('stop_reason') is not None:
        raise RuntimeError('STAGE_NOT_COMPLETE')
    if monitor.get('process_group_teardown', {}).get('remaining_group_members') != []:
        raise RuntimeError('STAGE_OWNED_GROUP_DISAPPEARANCE_UNPROVEN')
    if stage == 'extension':
        plan = read(CONTINUATION, CONTINUATION_SHA)
        admission = read(ADMISSION, ADMISSION_SHA)
        failed = read(FAILED, FAILED_SHA)
        if failed.get('stop_reason') != 'INTERNAL_FULL_NATIVE_FLOOR_REACHED' or failed.get('returncode') != -15:
            raise RuntimeError('PRESERVED_PREDECESSOR_FAILURE_IDENTITY_REQUIRED')
        if monitor.get('continuation_plan_sha256') != CONTINUATION_SHA or monitor.get('admission_sha256') != ADMISSION_SHA:
            raise RuntimeError('EXACT_CONTINUATION_ADMISSION_REQUIRED')
        if monitor.get('prior_failed_monitor_sha256') != FAILED_SHA or monitor.get('old_failed_attempt_reclassified') is not False:
            raise RuntimeError('PREDECESSOR_FAILURE_MUST_REMAIN_FAILED')
        if (monitor.get('unchanged_successful_jobs_reused'), monitor.get('remaining_frozen_jobs')) != (107, 5):
            raise RuntimeError('EXACT107_PLUS5_CHECKPOINT_CONTINUATION_REQUIRED')
        if monitor.get('owned_cleanup_verified') is not True or monitor.get('termination_requests') != [] or monitor.get('additional_failures'):
            raise RuntimeError('CONTINUATION_CLEANUP_OR_TERMINATION_FAILURE')
        if admission.get('execution_admitted') is not True or admission.get('plan_sha256') != CONTINUATION_SHA or admission.get('executor_sha256') != CONTROLLER_SHA:
            raise RuntimeError('ROOT_ADMISSION_SCHEMA_OR_CONTROLLER_CHANGED')
        read_controller = PACKAGE / 'scripts/61_native_checkpoint_continuation_v3.py'
        if sha(read_controller) != CONTROLLER_SHA:
            raise RuntimeError('ADMITTED_CONTINUATION_CONTROLLER_CHANGED')
        consumed[str(read_controller)] = CONTROLLER_SHA
        if record is not None:
            record(read_controller)
        for review, expected in admission['independent_review_sha256'].items():
            read(review, expected)
        for dependency, expected in plan['dependencies_sha256'].items():
            actual = sha(dependency)
            if actual != expected:
                raise RuntimeError('CONTINUATION_FROZEN_DEPENDENCY_CHANGED: ' + dependency)
            consumed[dependency] = actual
            if record is not None:
                record(dependency)
        preservation = monitor.get('interrupted_preservation_sha256')
        if not isinstance(preservation, dict) or len(preservation) != 2:
            raise RuntimeError('PAIRED_INTERRUPTED_PRESERVATION_PROOF_REQUIRED')
        preservation_receipts = []
        for proof, expected in preservation.items():
            if sha(proof) != expected:
                raise RuntimeError('INTERRUPTED_PRESERVATION_PROOF_CHANGED')
            consumed[proof] = expected
            if record is not None:
                record(proof)
            if proof.endswith('.json'):
                preservation_receipts.append(json.loads(Path(proof).read_text()))
        if len(preservation_receipts) != 1:
            raise RuntimeError('EXACT_PRESERVATION_RECEIPT_REQUIRED')
        pr = preservation_receipts[0]
        if pr.get('status') != 'ALL284_BIOLOGICAL_AND284_SIDECAR_ARTIFACTS_PRESERVED_WITHOUT_DELETION' or pr.get('artifacts') != plan['interrupted_artifacts'] or pr.get('transport_sidecars') != plan['interrupted_sidecars'] or pr.get('sidecars_count_as_biological_estimates') is not False:
            raise RuntimeError('EXACT568_PAIRED_PRESERVATION_RECEIPT_REQUIRED')
        worker = read(monitor['worker_receipt'], monitor['worker_receipt_sha256'])
        if worker.get('status') != 'WORKER_COMPLETE_VERIFIED' or worker.get('process_group_teardown', {}).get('remaining_group_members') != []:
            raise RuntimeError('ORIGINAL_RUNNER_CONTINUATION_WORKER_NOT_COMPLETE')
        jobs = [j for j in historical['jobs'] if j['stage'] == 'extension']
        expected_paths = {str(Path(historical['ssd_support_package']) / 'native/extension_reproduction_v1' / (j['job_id'] + '.execution_receipt.json')) for j in jobs}
        captured = monitor.get('native_receipt_sha256')
        if len(expected_paths) != 112 or not isinstance(captured, dict) or set(captured) != expected_paths:
            raise RuntimeError('EXACT112_EXTENSION_RECEIPT_BINDINGS_REQUIRED')
        for receipt, expected in captured.items():
            read(receipt, expected)
        state = monitor.get('final_resource_snapshot', {})
        if state.get('internal_free_bytes', -1) < 3 * 1024**3 or state.get('ssd_free_bytes', -1) < 5 * 1024**3:
            raise RuntimeError('CONTINUATION_FINAL_FROZEN_RESOURCE_FLOORS_FAILED')
    if any(sha(p) != expected for p, expected in consumed.items()):
        raise RuntimeError('STAGE_EVIDENCE_CHANGED_DURING_GATE')
    if Path(str(path) + '.failure.json').exists():
        raise RuntimeError('LATE_TERMINAL_FAILURE_INVALIDATES_STAGE')
    return monitor
