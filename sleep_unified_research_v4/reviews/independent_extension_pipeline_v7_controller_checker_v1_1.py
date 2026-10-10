"""Whole actual controller, fabricated metadata workers/proof/resource/lock IO.

No subprocess/estimator/curl, real mutex, source/reference/runtime body reader,
admission or preparer is invoked. Six final-oracle fixtures deliberately replace
one stored worker/hash at the mocked family boundary to isolate final oracle
coverage after healthy immediate consumption; this is not a real producer claim.
"""
import copy
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
from types import SimpleNamespace

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
S = P / 'scripts'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
ACTUAL_PLAN = SSD / 'extension_pipeline_replay_v7/extension_pipeline_replay_plan_v7.json'
OUT = R / 'independent_extension_pipeline_v7_control_fixtures_v1_1'
sys.path.insert(0, str(S))
import terminal_commit_common_v2 as terminal_helper

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def mutate(record, kind):
    if kind == 'metadata_errors': record['metadata_errors'] = ['PRIVATE_MOCK_ERROR']
    elif kind == 'plan_unchanged': record['plan_unchanged'] = False
    elif kind == 'cleanup_error': record['process_group_teardown']['cleanup_error'] = 'PRIVATE_MOCK_ERROR'
    elif kind == 'wrong_command': record['command'] += ['UNFROZEN_ARGUMENT']
    elif kind == 'returncode': record['returncode'] = 17
    elif kind == 'stop_reason': record['stop_reason'] = 'PRIVATE_MOCK_STOP'
    elif kind == 'returncode_bool': record['returncode'] = False
    elif kind == 'empty_outputs': record['output_sha256'] = {}
    else: raise AssertionError(kind)

def run():
    started = time.monotonic()
    assert sha(ACTUAL_PLAN) == 'ad0c08bc0b0a62c8a1bc070221689f2733b0eefe8f25de8cc1bd8689a81b1ba5'
    actual = json.loads(ACTUAL_PLAN.read_text())
    actual_out = ACTUAL_PLAN.parent
    candidate = S / '47_run_extension_pipeline_replay_v7.py'
    candidate_sha = sha(candidate)
    assert candidate_sha == '6abd57c9d743ea1c33b7ccbed0f88c51a7af0972db24655284f3e279fa483f00'
    OUT.mkdir(exist_ok=False)
    six = ['metadata_errors', 'plan_unchanged', 'cleanup_error', 'wrong_command', 'returncode', 'stop_reason']
    labels = ['healthy400'] + ['first_' + x for x in six] + ['final_' + x for x in six] + [
        'receipt_after_parse', 'receipt_during_output_consumption', 'failure_addendum', 'receipt_samebytes_symlink',
        'output_samebytes_symlink', 'empty_outputs', 'returncode_bool', 'master_postwrite_mutation',
        'final_missing_worker_key', 'old_output_after_master_write', 'postpersist_internal_floor',
        'postpersist_deferred_signal', 'failed_supplemental_writes', 'owned_cleanup_before_release']
    results = []
    for label in labels:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 128 << 20 or time.monotonic() - started > 180:
            raise RuntimeError('INDEPENDENT_METADATA_REVIEW_RESOURCE_BOUND')
        folder = OUT / label
        campaign = folder / 'campaign'
        namespace = campaign / 'extension_pipeline_replay_v7'
        for child in ['harmonized', 'munged', 'qc', 'receipts', 'logs_v4', 'receipts_v4', 'tmp', 'cache']:
            (namespace / child).mkdir(parents=True, exist_ok=True)
        runner = load('_independent_extension7_' + label, candidate)
        runner.SSD, runner.OUT = campaign, namespace
        # Compact invented operational plan only; separate audit binds full real plan.
        top_keys=['guard','members','python','acquisition_execution_binding','acquisition_plan_sha256','checkpoint_receipt_binding_policy']
        plan = copy.deepcopy({k:actual[k] for k in top_keys})
        member_keys=['extension_trait_id','index','acquisition_receipt','harmonized','harmonization_qc','harmonization_receipt','munged','munged_prefix','source_gate_receipt','comparison_receipt','harmonize_command','munge_command','archived_munged_sha256']
        plan['members']=[{k:m[k] for k in member_keys} for m in plan['members']]
        binding=plan['acquisition_execution_binding']
        plan['acquisition_execution_binding']={k:binding[k] for k in ['root_admission','executor_sha256','sha256']}
        plan['acquisition_execution_binding']['sha256']={binding['root_admission']:binding['sha256'][binding['root_admission']]}
        plan['guard']['shared_heavy_worker_lock'] = str(campaign / 'native_heavy_worker.lock')
        for member in plan['members']:
            for key, value in list(member.items()):
                if isinstance(value, str) and value.startswith(str(actual_out) + '/'):
                    member[key] = value.replace(str(actual_out) + '/', str(namespace) + '/', 1)
                elif key in ['harmonize_command', 'munge_command']:
                    member[key] = [x.replace(str(actual_out) + '/', str(namespace) + '/', 1) if x.startswith(str(actual_out) + '/') else x for x in value]
        plan_path, admission_path = folder / 'plan.json', folder / 'admission.json'
        save(plan_path, plan)
        ph = sha(plan_path)
        review = folder / 'review.txt'
        review.write_text('PRIVATE MOCK ROOT REVIEW IDENTITY\n')
        acq = plan['acquisition_execution_binding']
        admission = dict(execution_admitted=True, plan_sha256=ph, executor_sha256=candidate_sha,
            acquisition_receipt_sha256={}, checkpoint_receipt_binding_policy=plan['checkpoint_receipt_binding_policy'],
            acquisition_execution_plan_sha256=plan['acquisition_plan_sha256'], acquisition_executor_sha256=acq['executor_sha256'],
            acquisition_root_admission_sha256=acq['sha256'][acq['root_admission']], independent_review_sha256={str(review): sha(review)})
        save(admission_path, admission)
        args = SimpleNamespace(plan=plan_path, plan_sha=ph, admission=admission_path, admission_sha=sha(admission_path))
        flags, counters = {}, dict(worker=0, acquire=0, close=0, cleanup=0, current_fd=None, baseline=0, checkpoint=0)
        target_worker = namespace / 'receipts_v4' / (plan['members'][0]['extension_trait_id'] + '__source.worker.json')
        target_output = Path(plan['members'][0]['source_gate_receipt'])
        lookup = {str(namespace / 'receipts_v4' / (m['extension_trait_id'] + '__' + kind + '.worker.json')): (m, kind)
            for m in plan['members'] for kind in ['source', 'harmonize', 'munge', 'compare']}
        baseline = {'PRIVATE_MOCK_COMPLETE190': 'a' * 64}
        family = {'PRIVATE_MOCK_COMPLETE100': True}
        protected = SimpleNamespace(TERMINATION_REQUEST=[])
        protected.catchable_termination = lambda n, frame: protected.TERMINATION_REQUEST.append(n)
        protected.check_termination = lambda where: (_ for _ in ()).throw(RuntimeError('PRIVATE_DEFERRED_TERMINATION_' + where)) if protected.TERMINATION_REQUEST else None
        def baseline_gate(_plan):
            counters['baseline'] += 1
            return baseline
        protected.baseline_gate = baseline_gate
        def stage_resource_gate(_plan, start, where):
            protected.check_termination(where)
            state, reason = protected.limits(_plan, start, 0)
            if reason: raise RuntimeError(reason)
            return state
        protected.stage_resource_gate = stage_resource_gate
        protected.safe_print = lambda *a, **k: None
        def await_cleanup(monitor, ownership, plan_sha):
            assert counters['current_fd'] is not None
            counters['cleanup'] += 1
            ownership[:] = [True, []]
            return []
        protected.await_owned_cleanup = await_cleanup
        def worker(command, prefix, receipt, _plan, path, plan_sha, family_start, monitor, ownership):
            assert counters['current_fd'] is not None
            counters['worker'] += 1
            member, kind = lookup[str(receipt)]
            keys = {'source': ['source_gate_receipt'], 'harmonize': ['harmonized', 'harmonization_qc', 'harmonization_receipt'],
                'munge': ['munged'], 'compare': ['comparison_receipt']}[kind]
            for key in keys:
                out = Path(member[key])
                if key == 'comparison_receipt': save(out, dict(status='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH', private_fixture_only=True))
                else: out.write_text('PRIVATE METADATA MARKER; NO GWAS CONTENT\n')
            journal = Path(str(receipt) + '.failure_journal.jsonl')
            journal.write_text('{"private_mock_metadata_worker":true}\n')
            outputs = {member[key]: sha(member[key]) for key in keys}
            outputs[str(journal)] = sha(journal)
            record = dict(status='WORKER_COMPLETE_VERIFIED', plan_sha256=ph, command=list(command), returncode=0,
                stop_reason=None, metadata_errors=[], plan_unchanged=True, owned_cleanup_verified=True,
                process_group_teardown=dict(remaining_group_members=[], cleanup_error=None, verified=True), output_sha256=outputs)
            if Path(receipt) == target_worker:
                if label.startswith('first_'): mutate(record, label[len('first_'):])
                elif label in ['empty_outputs', 'returncode_bool']: mutate(record, label)
            save(receipt, record)
            if Path(receipt) == target_worker:
                if label == 'failure_addendum': save(Path(str(receipt) + '.failure.json'), {'private_failure': True})
                if label in ['receipt_samebytes_symlink', 'output_samebytes_symlink']:
                    target = Path(receipt) if label == 'receipt_samebytes_symlink' else target_output
                    copy_path = folder / 'same_bytes_copy.txt'
                    copy_path.write_bytes(target.read_bytes())
                    target.unlink()
                    target.symlink_to(copy_path)
                if label == 'owned_cleanup_before_release':
                    ownership[:] = [False, ['PRIVATE_FAKE_OWNED_MEMBER']]
                    raise RuntimeError('PRIVATE_WORKER_FAILURE_WITH_MOCK_OWNERSHIP')
        protected.worker = worker
        runner.load_protected = lambda: protected
        runner.physical_mount = lambda: {'PRIVATE_MOCK_PHYSICAL_SSD': True}
        runner.check_bindings = lambda *a, **k: None
        runner.acquisition_operational_gate = lambda *a, **k: None
        runner.checkpoint_binding_gate = lambda *a, **k: None
        runner.acquisition_family_gate = lambda *a, **k: family
        def wait_checkpoint(_plan, member, _a, _admission, _protected, start):
            assert counters['current_fd'] is None
            counters['checkpoint'] += 1
            return dict(plan_sha256=ph, pipeline_root_admission_sha256=args.admission_sha,
                extension_trait_id=member['extension_trait_id'], acquisition_receipt=member['acquisition_receipt'],
                acquisition_receipt_sha256='b' * 64, acquisition_execution_binding_sha256=_plan['acquisition_execution_binding']['sha256'])
        runner.wait_checkpoint = wait_checkpoint
        def wait_family(_plan, _protected, start, receipts):
            assert counters['current_fd'] is None and len(receipts) == 100
            master = inspect.currentframe().f_back.f_locals['master']
            if label.startswith('final_') and label != 'final_missing_worker_key':
                record = json.loads(target_worker.read_text())
                mutate(record, label[len('final_'):])
                save(target_worker, record)
                master['worker_receipt_sha256'][str(target_worker)] = sha(target_worker)
                flags['final_oracle_fixture_injection'] = True
            if label == 'final_missing_worker_key': master['worker_receipt_sha256'].pop(str(target_worker))
            return family
        runner.wait_family_commit = wait_family
        def acquire(*a):
            assert counters['current_fd'] is None
            counters['acquire'] += 1
            counters['current_fd'] = 917
            return 917
        runner.acquire_source_mutex = acquire
        def close(fd):
            assert fd == counters['current_fd'] == 917
            counters['close'] += 1
            counters['current_fd'] = None
        runner.os = SimpleNamespace(getpid=os.getpid, close=close)
        runner.signal = SimpleNamespace(SIGINT=signal.SIGINT, SIGTERM=signal.SIGTERM, SIGHUP=signal.SIGHUP, signal=lambda *a: None)
        runner.shutil = SimpleNamespace(disk_usage=lambda path: SimpleNamespace(free=2 << 30 if flags.get('low_internal') and str(path) == '/System/Volumes/Data' else 900 << 30))
        runner.namespace_bytes = lambda path: 1 << 20
        runner.time = SimpleNamespace(monotonic=time.monotonic, sleep=lambda seconds: (_ for _ in ()).throw(AssertionError('NO_WAIT_IN_FIXTURE')))
        original_writer, original_helper_save = runner.write_new, terminal_helper.save_new
        def writer(path, value):
            if label == 'failed_supplemental_writes' and str(path).endswith('.failure.json'): raise OSError('PRIVATE_SECONDARY_WRITE_FAILURE')
            original_writer(path, value)
            if Path(path) == namespace / 'pipeline_execution_receipt_v7.json':
                if label == 'master_postwrite_mutation': Path(path).write_text(Path(path).read_text() + ' ')
                if label == 'old_output_after_master_write': target_output.write_text('PRIVATE POSTPERSIST DRIFT\n')
                if label in ['postpersist_internal_floor', 'failed_supplemental_writes']: flags['low_internal'] = True
                if label == 'postpersist_deferred_signal': protected.TERMINATION_REQUEST.append(signal.SIGTERM)
        runner.write_new = writer
        def helper_save(path, value):
            if label == 'failed_supplemental_writes' and str(path).endswith('.failure.json'): raise OSError('PRIVATE_TERMINAL_SECONDARY_WRITE_FAILURE')
            original_helper_save(path, value)
        terminal_helper.save_new = helper_save
        def parsed(payload):
            value = json.loads(payload)
            if isinstance(value, dict) and value.get('status') == 'WORKER_COMPLETE_VERIFIED' and not flags.get('parsed_target'):
                flags['parsed_target'] = True
                if label == 'receipt_after_parse': target_worker.write_text(target_worker.read_text() + ' ')
                if label == 'receipt_during_output_consumption': flags['consume_mutation_ready'] = True
            return value
        runner.json = SimpleNamespace(loads=parsed, dumps=json.dumps)
        def current_sha(path):
            value = sha(path)
            if Path(path) == target_output and flags.get('consume_mutation_ready'):
                target_worker.write_text(target_worker.read_text() + ' ')
                flags['consume_mutation_ready'] = False
            return value
        runner.sha = current_sha
        accepted, error, consumed = False, None, False
        try:
            try: runner.run(args); accepted = True
            except BaseException as exc: error = type(exc).__name__ + ': ' + str(exc)
            master_path = namespace / 'pipeline_execution_receipt_v7.json'
            if master_path.exists():
                try:
                    terminal_helper.require_committed(namespace / 'pipeline_pending_v7.json', namespace / 'pipeline_terminal_seal_v7.json',
                        {'plan_sha256': ph, 'admission_sha256': args.admission_sha, 'executor_sha256': candidate_sha}, {str(master_path): sha(master_path)})
                    consumed = True
                except BaseException: pass
            record = json.loads(master_path.read_text()) if master_path.exists() else {}
            expected = label == 'healthy400'
            passed = accepted == expected and consumed == expected and counters['current_fd'] is None and counters['acquire'] == counters['close']
            if not expected: passed = passed and (namespace / 'pipeline_pending_v7.json').exists()
            if label == 'healthy400': passed = passed and counters['worker'] == 400 and len(record['completed_members']) == 100 and len(record['worker_receipt_sha256']) == 400
            if label == 'owned_cleanup_before_release': passed = passed and counters['cleanup'] == 1
            results.append(dict(label=label, passed=passed, execute_success=accepted, Terminal2_consumer_success=consumed,
                expected_success=expected, error=error, preserved_primary_status=record.get('status'),
                pending_retained=(namespace / 'pipeline_pending_v7.json').exists(), metadata_worker_calls=counters['worker'],
                completed_members=len(record.get('completed_members', [])), worker_map_size=len(record.get('worker_receipt_sha256', {})),
                mock_mutex_acquire_close_counts=[counters['acquire'], counters['close']], mock_cleanup_calls=counters['cleanup'],
                final_oracle_isolation_state_injection=flags.get('final_oracle_fixture_injection', False), actual_worker_lock_calls=0))
        finally: terminal_helper.save_new = original_helper_save
        print(json.dumps({'control': label, 'passed': results[-1]['passed'], 'metadata_workers': counters['worker']}), flush=True)
    receipt = dict(schema='independent_extension7_whole_controller_metadata_controls_v1_1', controls=results,
        all_expected_results=all(row['passed'] for row in results), controls_count=len(results), reviewed_controller_sha256=candidate_sha,
        controller_unchanged=sha(candidate) == candidate_sha, actual_plan_sha256=sha(ACTUAL_PLAN),
        fixture_regular_sha256={str(path): sha(path) for path in sorted(OUT.rglob('*')) if path.is_file() and not path.is_symlink()},
        fixture_symlinks={str(path): os.readlink(path) for path in OUT.rglob('*') if path.is_symlink()},
        actual_worker_estimator_curl_mutex_body_runtime_read_calls=0, old402_1729_or190_reruns=0,
        qualifications='Science/source/acquisition/full190/physical resource probes and lock/process operations are mocked; current exact controller, strict receipt oracle, cardinality and Terminal2 persistence behavior execute on new private metadata fixtures. Final six use explicit private state injection to isolate final oracle, not a genuine producer failure claim.',
        elapsed_seconds=time.monotonic() - started, max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    receipt['fixture_bytes'] = sum(Path(path).stat().st_size for path in receipt['fixture_regular_sha256'])
    assert receipt['fixture_bytes'] < 64 << 20
    target = R / 'independent_extension_pipeline_v7_controller_receipt_v1_1.json'
    with target.open('x') as stream: json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({'all_pass': receipt['all_expected_results'], 'controls': len(results), 'peak_RSS': receipt['max_RSS_bytes'], 'fixture_bytes': receipt['fixture_bytes']}))

if __name__ == '__main__':
    run()
