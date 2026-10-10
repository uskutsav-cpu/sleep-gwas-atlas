#!/usr/bin/env python3
"""Narrow v3 frozen-plan and synthetic terminal-integration review.

Imports define helpers only. Every worker, source/body operation and baseline in
the scheduling witness is labelled synthetic; real subprocess creation is banned.
"""
import ast
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

P = Path(__file__).resolve().parents[1]
S = P / 'scripts'
V = P / 'statistical_validation'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = SSD / 'extension_pipeline_replay_v3/extension_pipeline_replay_plan_v3.json'
PRIOR = SSD / 'extension_pipeline_replay_v2/extension_pipeline_replay_plan_v2.json'
RAW = P / 'manifests/extension_raw_acquisition_plan_v4_6.json'
PINNED = {
    str(PLAN): '565b34e997b2cab841f3a101cb4f8649c4c4660e56f5b1421ce4f68f29afe4ac',
    str(PRIOR): '77f008725b8d87ebe5e08c85f8667c8b55f02a7065ff3a6e46a7c45d3fd38cb3',
    str(RAW): '393af7ffc3231ceddf829ae669558001a7411bb4a7917b950a713c41546237fe',
    str(S / '47_run_extension_pipeline_replay_v3.py'): 'd207e1f703a032408c08e5adcf16bace35fab1e691f377156499bb760dda3c90',
    str(S / 'extension_replay_common_v3.py'): 'ca5eaf792eb4a06ef04d323249472452774eada0e54109b645d1b95df6d5bee4',
    str(S / 'extension_replay_validate_v3.py'): 'acfdaf61e7597aebf10aecb5cdfe6d0f2c0fd02b409ca989306536c60d6828e5',
    str(S / '46_prepare_extension_pipeline_replay_v3.py'): 'b804f4fbc31385a13d9a8bf051fa52d3d39035dc799b761225ed72b256e549f1',
    str(S / 'terminal_commit_common_v2.py'): '9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd',
    str(S / 'sensitivity_executor_v4_4.py'): '7ac802596f29adfacb017f406676e6a30ee00bf890904976afda7f88a3ffa9c5',
    str(P / 'reviews/genomicsem_extension_checkpoint_review_seal_v2.json'): 'f7bf85183737644bc49b3760fffd65d73f63dd365325c7985afdd09f8239bc9d',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def functions(path):
    return {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef)}


def run():
    checks = []
    def check(value, label):
        if not value:
            raise AssertionError(label)
        checks.append({'case': label, 'pass': True})
    for p, h in PINNED.items():
        check(sha(p) == h, 'pinned_identity ' + p)
    old_seal = json.loads((P / 'reviews/genomicsem_extension_checkpoint_review_seal_v2.json').read_text())
    for p, meta in old_seal['artifacts_sha256'].items():
        digest = meta['sha256'] if isinstance(meta, dict) else meta
        check(sha(p) == digest, 'preserved_v2_artifact ' + p)
    old, plan, raw = [json.loads(p.read_text()) for p in [PRIOR, PLAN, RAW]]
    check(len(old['members']) == len(plan['members']) == len(raw['members']) == 100, 'exact100_three_plan_cardinality')
    for a, b, r in zip(old['members'], plan['members'], raw['members']):
        expected = copy.deepcopy(a)
        for k, v in list(expected.items()):
            if isinstance(v, str) and v.startswith(str(PRIOR.parent) + '/'):
                expected[k] = v.replace(str(PRIOR.parent) + '/', str(PLAN.parent) + '/', 1)
            elif k in ['harmonize_command', 'munge_command']:
                expected[k] = [x.replace(str(PRIOR.parent) + '/', str(PLAN.parent) + '/', 1) if x.startswith(str(PRIOR.parent) + '/') else r['body_path'] if x == a['raw'] else x for x in v]
        expected['raw'], expected['acquisition_member'] = r['body_path'], r
        reuse = raw['reused_checkpoints'].get(str(r['index']))
        expected['acquisition_receipt'] = reuse['receipt_path'] if reuse else str(Path(raw['pending_path']).parent / 'receipts' / (r['extension_trait_id'] + '.json'))
        expected['acquisition_origin_plan_sha256'] = reuse['original_acquisition_plan_sha256'] if reuse else PINNED[str(RAW)]
        check(expected == b, 'exact_path_origin_only_member_transformation ' + b['extension_trait_id'])
        check({k: v for k, v in a['acquisition_member'].items() if k != 'body_path'} == {k: v for k, v in r.items() if k != 'body_path'}, 'original_source_science_identity ' + b['extension_trait_id'])
    scientific_keys = ['reference', 'reference_sha256', 'w_hm3', 'template_rows', 'environment', 'archived_input_sha256', 'original_190_jobs', 'baseline_execution_plan_sha256', 'baseline_input_sha256', 'baseline_dependency_sha256', 'baseline_relocated_dependency_sha256']
    for k in scientific_keys:
        check(old[k] == plan[k], 'unchanged_scientific_plan_field ' + k)
    check(all(plan['dependencies_sha256'].get(k) == v for k, v in old['dependencies_sha256'].items()), 'all373_prior_dependencies_retained')
    check(all(plan['dependencies_sha256'].get(k) == v for k, v in raw['bound_sources'].items()), 'all137_v6_bound_sources_retained')
    check(sum(m['acquisition_origin_plan_sha256'] == old['acquisition_plan_sha256'] for m in plan['members']) == 2, 'two_preserved_original_receipt_origins')
    check(sum(m['acquisition_origin_plan_sha256'] == PINNED[str(RAW)] for m in plan['members']) == 98, '98_new_v6_receipt_origins')
    expected_guard = copy.deepcopy(old['guard']); expected_guard['family_terminal_assembly_seconds'] = 600
    check(plan['guard'] == expected_guard, 'only_guard_addition600_second_family_commit_wait')
    for p, h in PINNED.items():
        if '/scripts/' in p:
            check(plan['dependencies_sha256'].get(p) == h, 'prepared_plan_exact_code_binding ' + p)
    old_f, new_f = functions(S / 'extension_replay_common_v2.py'), functions(S / 'extension_replay_common_v3.py')
    for name in old_f:
        if name not in ('acquisition_receipt_gate', 'checkpoint_binding_gate'):
            check(old_f[name] == new_f[name], 'unchanged_common_function ' + name)
    check(functions(S / 'extension_replay_validate_v2.py') == functions(S / 'extension_replay_validate_v3.py'), 'validator_functions_unchanged_only_common_import_changes')
    prot_old, prot_new = functions(S / 'sensitivity_executor_v4_3_1.py'), functions(S / 'sensitivity_executor_v4_4.py')
    for name in ('worker', 'await_owned_cleanup', 'baseline_monitor_binding_gate', 'baseline_output_binding_gate', 'admit_exact_core_adjudication'):
        check(prot_old[name] == prot_new[name], 'protected_function_unchanged ' + name)
    sys.path.insert(0, str(S))
    common = module(S / 'extension_replay_common_v3.py', '_common_narrow_v3')
    runner = module(S / '47_run_extension_pipeline_replay_v3.py', '_runner_narrow_v3')
    protected = module(S / 'sensitivity_executor_v4_4.py', '_protected_narrow_v3')
    terminal_mod = module(S / 'terminal_commit_common_v2.py', '_terminal_narrow_v3')
    observations = []
    with tempfile.TemporaryDirectory(prefix='extension-v3-review-') as td:
        tmp = Path(td)
        # Integration only: generic terminal algorithm receives real isolated
        # markers, a tiny receipt and labelled gate callbacks; no full re-audit.
        for label in ['healthy', 'post_persistence_resource_failure', 'late_signal', 'receipt_failure']:
            folder = tmp / label; folder.mkdir()
            target = folder / 'receipt.json'; pending = folder / 'pending.json'; seal = folder / 'seal.json'
            terminal = terminal_mod.TerminalCommit(pending, seal, {'fixture_only': True})
            protected.TERMINATION_REQUEST.clear(); events = []; late = [False]
            master = {'status': 'ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS'}
            def limits(*_):
                events.append('resource_query_late=' + str(late[0]))
                return {'fixture_only': True}, ('SSD_SPACE_GUARD' if late[0] and label == 'post_persistence_resource_failure' else None)
            def writer(path, obj):
                if Path(path) == target and label == 'receipt_failure':
                    raise OSError('CONTROLLED_FINAL_RECEIPT_IO_FAILURE')
                common.write_new(path, obj)
                if Path(path) == target:
                    events.append('persisted'); late[0] = True
                    if label == 'late_signal':
                        protected.catchable_termination(signal.SIGTERM, None)
            with patch.object(protected, 'limits', limits), patch.object(runner, 'write_new', writer):
                runner.seal_master(protected, master, target, {}, 0, [True, []], terminal, lambda: events.append('identity_query'))
            if label == 'healthy':
                check(master['status'] == 'ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS' and not pending.exists() and seal.is_file(), 'healthy_actual_terminal_integration_commits')
                terminal_mod.require_committed(pending, seal, {'fixture_only': True}, {str(target): sha(target)})
                check(events.count('resource_query_late=True') == 2, 'healthy_two_post_persistence_resource_queries')
            else:
                check(master['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and pending.exists() and Path(str(target) + '.failure.json').exists(), 'pending_and_failure_addendum_retained ' + label)
            if label == 'post_persistence_resource_failure':
                check('resource_query_late=True' in events, 'v2_post_persistence_false_PASS_witness_now_rejected')
                observations.append({'witness': 'v2_post_persistence_resource_boundary', 'v3_status': master['status'], 'events': events, 'durable_pending_retained': pending.exists()})
        protected.TERMINATION_REQUEST.clear()
        # Original400-command scheduler witness, preserving its exact timing:
        # family fails inside the final synthetic compare, after final discovery.
        schedule = tmp / 'scheduler'; schedule.mkdir()
        for n in ['receipts', 'receipts_v4', 'harmonized', 'munged', 'qc', 'logs_v4', 'tmp', 'cache']:
            (schedule / n).mkdir()
        members = []
        for i in range(1, 101):
            t = 'CONTROL_ONLY_' + str(i)
            members.append(dict(index=i, extension_trait_id=t, acquisition_receipt=str(schedule / (t + '.acquired.json')), harmonized=str(schedule / 'harmonized' / (t + '.txt.gz')), harmonization_qc=str(schedule / 'qc' / (t + '.tsv')), harmonization_receipt=str(schedule / 'receipts' / (t + '.harm.json')), munged_prefix=str(schedule / 'munged' / t), munged=str(schedule / 'munged' / (t + '.sumstats.gz')), source_gate_receipt=str(schedule / 'receipts' / (t + '.source.json')), comparison_receipt=str(schedule / 'receipts' / (t + '.compare.json')), archived_munged_sha256='FIXTURE_ONLY_NO_ARCHIVE', harmonize_command=['FIXTURE_ONLY_HARMONIZER_NEVER_EXECUTED'], munge_command=['FIXTURE_ONLY_MUNGE_NEVER_EXECUTED']))
        family = schedule / 'acquisition_primary.json'
        fp = dict(guard=copy.deepcopy(plan['guard']), members=members, python='FIXTURE_ONLY_NO_RUNTIME', acquisition_family_terminal_contract={'primary_receipt_paths': [str(family), str(schedule / 'acquisition_second_copy.json')], 'pending_path': str(schedule / 'acquisition_pending.json'), 'terminal_seal_path': str(schedule / 'acquisition_seal.json')})
        fixture_plan, fixture_admission = schedule / 'plan.json', schedule / 'admission.json'
        common.write_new(fixture_plan, fp); common.write_new(fixture_admission, {'fixture_only': True})
        args = SimpleNamespace(plan=fixture_plan, plan_sha=sha(fixture_plan), admission=fixture_admission, admission_sha=sha(fixture_admission))
        calls = []
        class FakeProtected:
            TERMINATION_REQUEST = []
            catchable_termination = staticmethod(lambda *_: None)
            check_termination = staticmethod(lambda *_: None)
            safe_print = staticmethod(lambda *_a, **_k: None)
            baseline_gate = staticmethod(lambda *_: {'fixture_only': True})
            stage_resource_gate = staticmethod(lambda *_: {'fixture_only': True})
            await_owned_cleanup = staticmethod(lambda *_: (_ for _ in ()).throw(RuntimeError('NO_OWNED_FIXTURE_PROCESSES_EXPECTED')))
            @staticmethod
            def worker(command, prefix, *_):
                calls.append(command)
                if '--mode' in command:
                    mode = command[command.index('--mode') + 1]
                    common.write_new(prefix, {'status': 'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' if mode == 'compare' else 'EXACT_RAW_SOURCE_GATE_PASS', 'fixture_only': True})
                    if len(calls) == 400:
                        common.write_new(family, {'status': 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT', 'completed_source_count': 100, 'stop_reason': 'FIXTURE_FINAL_FAMILY_GATE_FAILURE'})
                elif command == ['FIXTURE_ONLY_HARMONIZER_NEVER_EXECUTED']:
                    prefix.write_bytes(b'FIXTURE_ONLY')
                elif command == ['FIXTURE_ONLY_MUNGE_NEVER_EXECUTED']:
                    Path(str(prefix) + '.sumstats.gz').write_bytes(b'FIXTURE_ONLY')
                else:
                    raise RuntimeError('UNRECOGNIZED_SYNTHETIC_WORKER')
        def fake_wait(_p, m, *_):
            return dict(plan_sha256=args.plan_sha, pipeline_root_admission_sha256=args.admission_sha, extension_trait_id=m['extension_trait_id'], acquisition_receipt=m['acquisition_receipt'], acquisition_receipt_sha256=hashlib.sha256(m['extension_trait_id'].encode()).hexdigest(), fixture_only=True)
        patches = [patch.object(runner, 'OUT', schedule), patch.object(runner, 'load_protected', lambda: FakeProtected()), patch.object(runner, 'admission_gate', lambda *_: None), patch.object(runner, 'check_bindings', lambda *_a, **_k: None), patch.object(runner, 'physical_mount', lambda: {'fixture_only': True}), patch.object(runner, 'wait_checkpoint', fake_wait), patch.object(runner, 'checkpoint_binding_gate', lambda *_: None), patch.object(runner, 'acquire_source_mutex', lambda *_: os.open(schedule / 'isolated_descriptor', os.O_RDWR | os.O_CREAT, 0o600)), patch('subprocess.Popen', side_effect=RuntimeError('FORBIDDEN_REAL_WORKER_LAUNCH'))]
        entered = []
        try:
            for p in patches:
                p.__enter__(); entered.append(p)
            try:
                runner.run(args)
            except SystemExit as e:
                check('PIPELINE_STOP_PRESERVED_REQUIRES_REVIEW' in str(e), 'scheduler_reports_preserved_failure')
        finally:
            for p in reversed(entered):
                p.__exit__(None, None, None)
        result = json.loads((schedule / 'pipeline_execution_receipt_v3.json').read_text())
        check(len(calls) == 400 and result['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and 'SOURCE_ACQUISITION_FAMILY_FAILED_PRESERVED' in result['error'], 'v2_late_failed_family_false_PASS_witness_now_rejected')
        check((schedule / 'pipeline_pending_v3.json').exists() and not (schedule / 'pipeline_terminal_seal_v3.json').exists(), 'failed_family_cannot_commit_pipeline_terminal')
        observations.append({'witness': 'v2_late_acquisition_family_failure', 'fixture_dispatches': len(calls), 'real_worker_dispatches': 0, 'v3_status': result['status'], 'error': result['error'], 'durable_pending_retained': True})
    for p, h in PINNED.items():
        check(sha(p) == h, 'pinned_identity_after ' + p)
    return {'recorded_utc': datetime.now(timezone.utc).isoformat(), 'status': 'NARROW_V3_EXPLICIT_PLAN_SOURCE_AND_TERMINAL_INTEGRATION_PASS', 'prepared_plan': str(PLAN), 'prepared_plan_sha256': PINNED[str(PLAN)], 'pinned_sha256_before_and_after': PINNED, 'helper_sha256': sha(__file__), 'checks': checks, 'check_count': len(checks), 'v2_false_PASS_witnesses_now_rejected': observations, 'prior1397_controls_reused_by_unchanged_science_functions_and_seal': True, 'source_origins_preserved': 2, 'source_origins_v6': 98, 'interface_qualification': 'Runner default plan filename is stale v2; root must supply explicit --plan with exact v3 path/SHA. Wrong omitted default fails before workers.', 'execution_admitted': False, 'real_worker_subprocesses_launched': 0, 'estimator_calls': 0, 'raw_or_processed_GWAS_body_reads': 0, 'genotype_decode_calls': 0, 'production_mutex_acquired': False, 'actual190_completion_claim': False, 'generic_terminal_full_controls_repeated': False, 'self_peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


if __name__ == '__main__':
    out = V / 'genomicsem_extension_checkpoint_review_controls_v3.json'
    if out.exists():
        raise RuntimeError('DISTINCT_REVIEW_OUTPUT_REQUIRED')
    result = run()
    with out.open('x') as f:
        f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'check_count', 'real_worker_subprocesses_launched', 'self_peak_RSS_bytes']}))
