#!/usr/bin/env python3
"""Independent metadata-only controls for the frozen downstream selector.

No campaign monitor is consumed and neither collator main nor any estimator is
invoked. All stage evidence below is synthetic and lives in temporary storage.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import tempfile
from datetime import datetime, timezone

P = Path(__file__).resolve().parents[1]
S = P / 'scripts'
SELECTOR = S / 'native_stage_completion_v4_3.py'
COLLATOR = S / '33_compare_native_campaign_v4_3.py'
PINNED = {
    str(SELECTOR): 'ed3b4d7af808d004447c955bf1ebafe6f70951ce9a178cab000ba865e13d487c',
    str(COLLATOR): 'abf7799709309f86af415badb244f810145de511d91f99f231bdce0a0356c47f',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, sort_keys=True) + '\n')
    return sha(path)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fixture(mod, root):
    mod.PACKAGE = root
    root.joinpath('scripts').mkdir(parents=True)
    controller = root / 'scripts/61_native_checkpoint_continuation_v3.py'
    controller.write_text('# Synthetic controller identity; never executed.\n')
    mod.CONTROLLER_SHA = sha(controller)
    jobs = [{'stage': 'extension', 'job_id': 'synthetic_%03d' % i} for i in range(112)]
    historical = root / 'historical.json'
    mod.HISTORICAL = historical
    mod.HISTORICAL_SHA = dump(historical, {'jobs': jobs, 'ssd_support_package': str(root / 'support')})
    failed = root / 'failed.json'
    mod.FAILED = failed
    mod.FAILED_SHA = dump(failed, {'stop_reason': 'INTERNAL_FULL_NATIVE_FLOOR_REACHED', 'returncode': -15})
    review = root / 'independent_review.md'
    review.write_text('Synthetic review only.\n')
    dep = root / 'dependency.py'
    dep.write_text('# Synthetic dependency.\n')
    artifacts = [{'synthetic_biological_identity': i} for i in range(284)]
    sidecars = [{'synthetic_transport_identity': i} for i in range(284)]
    plan = root / 'continuation.json'
    mod.CONTINUATION = plan
    mod.CONTINUATION_SHA = dump(plan, {'dependencies_sha256': {str(dep): sha(dep)},
                                      'interrupted_artifacts': artifacts,
                                      'interrupted_sidecars': sidecars})
    admission = root / 'admission.json'
    mod.ADMISSION = admission
    mod.ADMISSION_SHA = dump(admission, {'execution_admitted': True,
                                        'plan_sha256': mod.CONTINUATION_SHA,
                                        'executor_sha256': mod.CONTROLLER_SHA,
                                        'independent_review_sha256': {str(review): sha(review)}})
    journal = root / 'preservation.jsonl'
    journal.write_text('{"event":"SYNTHETIC_ONLY"}\n')
    preservation = root / 'preservation.json'
    pr = {'status': 'ALL284_BIOLOGICAL_AND284_SIDECAR_ARTIFACTS_PRESERVED_WITHOUT_DELETION',
          'artifacts': artifacts, 'transport_sidecars': sidecars,
          'sidecars_count_as_biological_estimates': False, 'journal_sha256': sha(journal)}
    dump(preservation, pr)
    worker = root / 'worker.json'
    wr = {'status': 'WORKER_COMPLETE_VERIFIED', 'process_group_teardown': {'remaining_group_members': []}}
    dump(worker, wr)
    receipts = {}
    for j in jobs:
        r = root / 'support/native/extension_reproduction_v1' / (j['job_id'] + '.execution_receipt.json')
        receipts[str(r)] = dump(r, {'synthetic': True, 'job_id': j['job_id'], 'audit_note': 'first'})
    mon = {'stage': 'extension', 'plan_sha256': mod.HISTORICAL_SHA, 'returncode': 0,
           'stop_reason': None, 'process_group_teardown': {'remaining_group_members': []},
           'continuation_plan_sha256': mod.CONTINUATION_SHA, 'admission_sha256': mod.ADMISSION_SHA,
           'prior_failed_monitor_sha256': mod.FAILED_SHA, 'old_failed_attempt_reclassified': False,
           'unchanged_successful_jobs_reused': 107, 'remaining_frozen_jobs': 5,
           'owned_cleanup_verified': True, 'termination_requests': [],
           'interrupted_preservation_sha256': {str(preservation): sha(preservation), str(journal): sha(journal)},
           'worker_receipt': str(worker), 'worker_receipt_sha256': sha(worker),
           'native_receipt_sha256': receipts,
           'final_resource_snapshot': {'internal_free_bytes': 3 * 1024**3, 'ssd_free_bytes': 5 * 1024**3},
           'audit_note': 'first'}
    mon_path = mod.stage_monitor_path('extension')
    dump(mon_path, mon)
    return dict(mon=mon, mon_path=mon_path, pr=pr, preservation=preservation,
                journal=journal, wr=wr, worker=worker, dep=dep, review=review,
                historical=historical, admission=admission, failed=failed)


def run():
    for path, expected in PINNED.items():
        if sha(path) != expected:
            raise RuntimeError('REVIEWED_SOURCE_CHANGED: ' + path)
    sys.path.insert(0, str(S))
    mod = module(SELECTOR, '_stage_selector_review_v4_3')
    col = module(COLLATOR, '_stage_collator_review_v4_3')
    checks = []
    with tempfile.TemporaryDirectory(prefix='native-stage-review-') as td:
        root = Path(td)

        def case(name, mutate=None, expected=None, callback=None):
            f = fixture(mod, root / ('case_%02d' % len(checks)))
            if mutate:
                mutate(f)
                dump(f['mon_path'], f['mon'])
            try:
                mod.stage_monitor('extension', record=(lambda path: callback(f, path)) if callback else None)
            except Exception as exc:
                if expected is None or expected not in str(exc):
                    raise
                checks.append({'case': name, 'pass': True, 'rejected': str(exc).split(':')[0]})
            else:
                if expected is not None:
                    raise AssertionError('EXPECTED_REJECTION_MISSING: ' + name)
                checks.append({'case': name, 'pass': True, 'accepted_synthetic_metadata': True})
            return f

        case('complete_synthetic_extension_metadata')
        changes = [
            ('wrong_stage', lambda f: f['mon'].update(stage='validation'), 'EXACT_STAGE_AND_HISTORICAL_PLAN_REQUIRED'),
            ('wrong_historical_plan', lambda f: f['mon'].update(plan_sha256='bad'), 'EXACT_STAGE_AND_HISTORICAL_PLAN_REQUIRED'),
            ('nonzero_returncode', lambda f: f['mon'].update(returncode=-15), 'STAGE_NOT_COMPLETE'),
            ('recorded_stop', lambda f: f['mon'].update(stop_reason='SYNTHETIC_STOP'), 'STAGE_NOT_COMPLETE'),
            ('remaining_group_member', lambda f: f['mon']['process_group_teardown'].update(remaining_group_members=[123]), 'STAGE_OWNED_GROUP_DISAPPEARANCE_UNPROVEN'),
            ('wrong_continuation_plan', lambda f: f['mon'].update(continuation_plan_sha256='bad'), 'EXACT_CONTINUATION_ADMISSION_REQUIRED'),
            ('wrong_admission', lambda f: f['mon'].update(admission_sha256='bad'), 'EXACT_CONTINUATION_ADMISSION_REQUIRED'),
            ('predecessor_reclassified', lambda f: f['mon'].update(old_failed_attempt_reclassified=True), 'PREDECESSOR_FAILURE_MUST_REMAIN_FAILED'),
            ('wrong_predecessor_identity', lambda f: f['mon'].update(prior_failed_monitor_sha256='bad'), 'PREDECESSOR_FAILURE_MUST_REMAIN_FAILED'),
            ('wrong_reuse_count', lambda f: f['mon'].update(unchanged_successful_jobs_reused=106), 'EXACT107_PLUS5_CHECKPOINT_CONTINUATION_REQUIRED'),
            ('wrong_remaining_count', lambda f: f['mon'].update(remaining_frozen_jobs=6), 'EXACT107_PLUS5_CHECKPOINT_CONTINUATION_REQUIRED'),
            ('cleanup_unverified', lambda f: f['mon'].update(owned_cleanup_verified=False), 'CONTINUATION_CLEANUP_OR_TERMINATION_FAILURE'),
            ('termination_request', lambda f: f['mon'].update(termination_requests=[15]), 'CONTINUATION_CLEANUP_OR_TERMINATION_FAILURE'),
            ('additional_failure', lambda f: f['mon'].update(additional_failures=['synthetic']), 'CONTINUATION_CLEANUP_OR_TERMINATION_FAILURE'),
            ('missing_preservation_proof', lambda f: f['mon']['interrupted_preservation_sha256'].pop(str(f['journal'])), 'PAIRED_INTERRUPTED_PRESERVATION_PROOF_REQUIRED'),
            ('wrong_preservation_identity', lambda f: f['mon']['interrupted_preservation_sha256'].update({str(f['preservation']): 'bad'}), 'INTERRUPTED_PRESERVATION_PROOF_CHANGED'),
            ('missing_native_receipt', lambda f: f['mon']['native_receipt_sha256'].pop(next(iter(f['mon']['native_receipt_sha256']))), 'EXACT112_EXTENSION_RECEIPT_BINDINGS_REQUIRED'),
            ('wrong_native_receipt_identity', lambda f: f['mon']['native_receipt_sha256'].update({next(iter(f['mon']['native_receipt_sha256'])): 'bad'}), 'STAGE_EVIDENCE_IDENTITY_CHANGED'),
            ('internal_floor_crossed', lambda f: f['mon']['final_resource_snapshot'].update(internal_free_bytes=3 * 1024**3 - 1), 'CONTINUATION_FINAL_FROZEN_RESOURCE_FLOORS_FAILED'),
            ('ssd_floor_crossed', lambda f: f['mon']['final_resource_snapshot'].update(ssd_free_bytes=5 * 1024**3 - 1), 'CONTINUATION_FINAL_FROZEN_RESOURCE_FLOORS_FAILED'),
            ('changed_dependency', lambda f: f['dep'].write_text('# changed\n'), 'CONTINUATION_FROZEN_DEPENDENCY_CHANGED'),
            ('changed_independent_review', lambda f: f['review'].write_text('changed\n'), 'STAGE_EVIDENCE_IDENTITY_CHANGED'),
            ('changed_failed_predecessor', lambda f: f['failed'].write_text('{}\n'), 'STAGE_EVIDENCE_IDENTITY_CHANGED'),
            ('initial_failure_addendum', lambda f: Path(str(f['mon_path']) + '.failure.json').write_text('{}\n'), 'TERMINAL_FAILURE_INVALIDATES_PROVISIONAL_STAGE_MONITOR'),
        ]
        for name, mutate, expected in changes:
            case(name, mutate, expected)
        def preservation_change(f, **changes):
            f['pr'].update(**changes)
            f['mon']['interrupted_preservation_sha256'][str(f['preservation'])] = dump(f['preservation'], f['pr'])
        case('wrong_preservation_status', lambda f: preservation_change(f, status='INCOMPLETE'), 'EXACT568_PAIRED_PRESERVATION_RECEIPT_REQUIRED')
        case('wrong_biological_map', lambda f: preservation_change(f, artifacts=[]), 'EXACT568_PAIRED_PRESERVATION_RECEIPT_REQUIRED')
        case('wrong_sidecar_map', lambda f: preservation_change(f, transport_sidecars=[]), 'EXACT568_PAIRED_PRESERVATION_RECEIPT_REQUIRED')
        case('sidecar_counted_as_biology', lambda f: preservation_change(f, sidecars_count_as_biological_estimates=True), 'EXACT568_PAIRED_PRESERVATION_RECEIPT_REQUIRED')
        def worker_change(f, **changes):
            f['wr'].update(**changes)
            f['mon']['worker_receipt_sha256'] = dump(f['worker'], f['wr'])
        case('failed_worker', lambda f: worker_change(f, status='WORKER_FAILED_PRESERVED'), 'ORIGINAL_RUNNER_CONTINUATION_WORKER_NOT_COMPLETE')
        case('worker_residual_member', lambda f: worker_change(f, process_group_teardown={'remaining_group_members': [123]}), 'ORIGINAL_RUNNER_CONTINUATION_WORKER_NOT_COMPLETE')
        def late_failure(f, path):
            if str(path) == list(f['mon']['native_receipt_sha256'])[-1]:
                Path(str(f['mon_path']) + '.failure.json').write_text('{}\n')
        case('late_failure_addendum', expected='LATE_TERMINAL_FAILURE_INVALIDATES_STAGE', callback=late_failure)
        def late_drift(f, path):
            if str(path) == list(f['mon']['native_receipt_sha256'])[-1]:
                f['review'].write_text('changed after review read\n')
        case('evidence_changed_within_one_gate', expected='STAGE_EVIDENCE_CHANGED_DURING_GATE', callback=late_drift)
        f = fixture(mod, root / 'original_stage_paths')
        for stage in ('core', 'validation'):
            p = mod.stage_monitor_path(stage)
            assert p.name == stage + '_native_monitor_receipt_v4.json'
            dump(p, {'stage': stage, 'plan_sha256': mod.HISTORICAL_SHA, 'returncode': 0,
                     'stop_reason': None, 'process_group_teardown': {'remaining_group_members': []}})
            mod.stage_monitor(stage)
            checks.append({'case': 'original_' + stage + '_monitor_path', 'pass': True})
        try:
            mod.stage_monitor_path('extension_v4_3')
        except ValueError:
            checks.append({'case': 'invalid_stage', 'pass': True})
        else:
            raise AssertionError('INVALID_STAGE_NOT_REJECTED')
        witnesses = []
        for target in ('selected_monitor', 'native_job_receipt'):
            f = fixture(mod, root / ('rebinding_' + target))
            col.INPUTS = {}
            mod.stage_monitor('extension', record=col.record)
            before = {p: dict(v) for p, v in col.INPUTS.items()}
            if target == 'native_job_receipt':
                victim = Path(next(iter(f['mon']['native_receipt_sha256'])))
                r = json.loads(victim.read_text())
                r['audit_note'] = 'changed between gate calls'
                f['mon']['native_receipt_sha256'][str(victim)] = dump(victim, r)
            else:
                victim = f['mon_path']
                f['mon']['audit_note'] = 'changed between gate calls'
            dump(f['mon_path'], f['mon'])
            mod.stage_monitor('extension', record=col.record)
            original_final_assertion_passes = all(col.sha(p) == m['sha256'] for p, m in col.INPUTS.items())
            assert before[str(victim)]['sha256'] != col.INPUTS[str(victim)]['sha256']
            assert original_final_assertion_passes
            witnesses.append({'case': target + '_first_observation_overwritten',
                              'witness_confirmed': True,
                              'initial_identity_changed': True,
                              'first_and_second_selector_gates_passed': True,
                              'original_final_unchanged_assertion_passed': True,
                              'real_campaign_evidence_consumed': False})
    assert all(sha(p) == s for p, s in PINNED.items())
    return {'recorded_utc': datetime.now(timezone.utc).isoformat(),
            'status': 'SYNTHETIC_SELECTOR_CONTROLS_PASS_WITH_CONFIRMED_COLLATOR_REBINDING_BLOCKER',
            'reviewed_sha256': PINNED, 'helper_sha256': sha(__file__),
            'passing_control_count': len(checks), 'controls': checks,
            'confirmed_blocker_witnesses': witnesses,
            'peak_review_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'real_workers_launched': 0, 'collator_main_called': False,
            'real_stage_monitor_consumed': False, 'GWAS_body_or_genotype_reads': 0,
            'scientific_or_covariance_outcomes_computed': False,
            'actual_extension_completion_claim': False,
            'source_identities_unchanged_after': True}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError('DISTINCT_REVIEW_OUTPUT_REQUIRED')
    result = run()
    with args.output.open('x') as out:
        out.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'passing_control_count', 'peak_review_rss_bytes', 'real_workers_launched')}))
