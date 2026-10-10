#!/usr/bin/env python3
"""Metadata-only v5 binding review and three bounded finalization witnesses.

Does not import or run any estimator, curl, body decoder or real acquire().
Only execute() family finalization is exercised with 100 synthetic checkpoints,
mocked transfer/identity helpers, real tiny immutable JSON writes and own flock.
"""
import ast
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import shutil
import tempfile
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
REVIEW = P / 'reviews'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = P / 'manifests/extension_raw_acquisition_plan_v4_5.json'
EXPECTED_PLAN = '48e3a26a518e382f093638ebc1de91faa4af519615723e4d512fd623c4794d13'
EXECUTOR = P / 'scripts/52_acquire_extension_raw_sources_v5.py'
OUT = REVIEW / 'independent_acquisition_v5_binding_receipt_v4.json'
SEEN = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def register(path, expected=None):
    path = Path(path)
    assert not path.name.endswith(('.bgz', '.bgz.partial', '.gz', '.gz.partial')), 'NO_BODY_READ_ALLOWED'
    info = {'sha256': sha(path), 'bytes': path.stat().st_size}
    if expected is not None:
        assert info['sha256'] == expected
    if str(path) in SEEN:
        assert SEEN[str(path)] == info
    SEEN[str(path)] = info
    return info


def jread(path):
    register(path)
    return json.loads(Path(path).read_text())


def module_from(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def finalization_witness(plan, scenario):
    """Use exact frozen execute() but make any real source/network work impossible."""
    m = module_from(EXECUTOR, 'independent_v5_' + scenario)
    writes = []
    with tempfile.TemporaryDirectory(prefix='acquisition_v5_metadata_', dir=SSD / 'tmp') as tmp:
        temp = Path(tmp)
        m.PACKAGE = temp / 'repo'
        m.FOLDER = temp / 'ssd'
        m.PLAN = temp / 'plan.json'
        m.ADMISSION = temp / 'admission.json'
        local_plan = json.loads(json.dumps(plan))
        local_plan['explicit_resume']['original_partial_path'] = str(temp / 'invented-prefix')
        # An invented three-byte fixture, never a study source or source excerpt.
        (temp / 'invented-prefix').write_bytes(b'abc')
        local_plan['explicit_resume']['prefix_bytes'] = 3
        (m.PLAN).write_text(json.dumps(local_plan))
        admission = {'execution_admitted': True, 'plan_sha256': EXPECTED_PLAN,
                     'executor_sha256': 'e' * 64,
                     'independent_review_artifact_sha256': {'synthetic.md': 'a' * 64, 'synthetic.json': 'b' * 64}}
        m.ADMISSION.write_text(json.dumps(admission))
        old_family_lock, old_write = m.family_lock, m.write_new
        @contextmanager
        def own_lock(ignored):
            with old_family_lock(temp / 'fixture.lock') as handle:
                yield handle
        m.family_lock = own_lock
        calls = {'snapshot': 0, 'transfer': 0, 'reuse': 0}
        def snapshot():
            calls['snapshot'] += 1
            # Family snapshots: initial; final before publication; safe_state
            # during publication; post-publication final guard.
            fail = scenario.startswith('post_guard_') and calls['snapshot'] == 4
            return {'internal_free_bytes': 0 if fail else 8 * 1024**3,
                    'ssd_free_bytes': 500 * 1024**3}
        monitor = SimpleNamespace(INTERNAL_FLOOR=plan['internal_floor_bytes'], SSD_FLOOR=plan['ssd_floor_bytes'],
                                  RSS_LIMIT=plan['maximum_owned_transfer_rss_bytes'], OUTPUT_LIMIT=plan['monitor_output_limit_bytes'], snapshot=snapshot)
        m.load_monitor = lambda: monitor
        m.physical_mount = lambda: None
        m.assert_bindings = lambda *args: None
        m.global_namespace_gate = lambda *args: 1
        m.subprocess = SimpleNamespace(run=lambda *args, **kwargs: SimpleNamespace(stdout=plan['curl_version']))
        def fake_hash(path):
            path = Path(path)
            if path == m.PLAN:
                return '0' * 32, EXPECTED_PLAN
            if path == EXECUTOR:
                return '0' * 32, 'e' * 64
            if path == temp / 'invented-prefix':
                return '0' * 32, local_plan['explicit_resume']['prefix_sha256']
            return '0' * 32, 'f' * 64
        m.hashes = fake_hash
        def checkpoint(member, key):
            calls[key] += 1
            return {'path': str(temp / ('synthetic-source-' + str(member['index']) + '.json')), 'sha256': 'f' * 64}
        m.reused_checkpoint_gate = lambda plan, member, monitor: checkpoint(member, 'reuse')
        m.acquire = lambda plan, member, *args: checkpoint(member, 'transfer')
        m.safe_print = lambda *args: None
        def write(path, value):
            path = Path(path)
            writes.append(str(path.relative_to(temp)))
            if scenario == 'second_success_copy_write_failure' and path == m.FOLDER / 'extension_raw_acquisition_family_receipt_v4_5.json':
                raise OSError('INDEPENDENT_SECOND_COPY_WRITE_FAILURE')
            if scenario == 'post_guard_first_invalidator_write_failure' and path == m.PACKAGE / 'logs/extension_raw_acquisition_family_receipt_v4_5.json.failure.json':
                raise OSError('INDEPENDENT_INVALIDATOR_WRITE_FAILURE')
            return old_write(path, value)
        m.write_new = write
        exception = None
        try:
            m.execute(EXPECTED_PLAN)
        except BaseException as exc:
            exception = type(exc).__name__ + ': ' + str(exc)
        paths = {'primary_repo': m.PACKAGE / 'logs/extension_raw_acquisition_family_receipt_v4_5.json',
                 'primary_ssd': m.FOLDER / 'extension_raw_acquisition_family_receipt_v4_5.json',
                 'invalidator_repo': m.PACKAGE / 'logs/extension_raw_acquisition_family_receipt_v4_5.json.failure.json',
                 'invalidator_ssd': m.FOLDER / 'extension_raw_acquisition_family_receipt_v4_5.json.failure.json'}
        observed = {key: {'exists': path.exists(), 'status': json.loads(path.read_text())['status'] if path.exists() else None}
                    for key, path in paths.items()}
        assert exception is not None and calls['transfer'] == 98
        unsafe = observed['primary_repo']['status'] == 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED' and not observed['invalidator_repo']['exists'] and not observed['invalidator_ssd']['exists']
        if scenario == 'post_guard_failure_addenda_pass':
            assert not unsafe and observed['invalidator_repo']['status'] == observed['invalidator_ssd']['status'] == 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
        else:
            assert unsafe
        return {'scenario': scenario, 'calls': calls, 'observed_artifacts': observed, 'exception': exception,
                'write_order': writes, 'provisional_success_without_failure_oracle_observed': unsafe,
                'real_curl_or_native_workers': 0, 'real_GWAS_body_reads': 0, 'invented_prefix_bytes': 3}


def main():
    start = time.monotonic()
    assert not OUT.exists()
    assert shutil.disk_usage('/System/Volumes/Data').free >= 256 * 1024**2
    plan = jread(PLAN)
    assert SEEN[str(PLAN)]['sha256'] == EXPECTED_PLAN
    register(SSD / 'manifests' / PLAN.name, EXPECTED_PLAN)
    for path, expected in plan['bound_sources'].items():
        register(path, expected)
    checkpoint = jread(REVIEW / 'independent_acquisition_checkpoint_receipt_v4_4.json')
    old = jread(P / 'manifests/extension_raw_acquisition_plan_v4_4.json')
    assert len(plan['members']) == len(checkpoint['proposed_mixed_origin_members']) == 100
    for original, proposed, member in zip(old['members'], checkpoint['proposed_mixed_origin_members'], plan['members']):
        assert all(member[key] == value for key, value in original.items() if key != 'body_path')
        assert member['body_path'] == proposed['body_path']
    assert set(plan['reused_checkpoints']) == {'1', '2'} and plan['network_sources_to_transfer'] == 98
    for key, fixed in plan['reused_checkpoints'].items():
        expected = checkpoint['source_receipt_proofs'][int(key) - 1]
        assert fixed == {'receipt_path': expected['receipt_path'], 'receipt_sha256': expected['receipt_sha256'],
                         'original_acquisition_plan_sha256': old['version'] and checkpoint['plan_sha256']}
    assert plan['explicit_resume']['source_index'] == 3
    resume = checkpoint['source3_resume']
    for key in ['source_trait', 'original_partial_path', 'prefix_bytes', 'prefix_sha256']:
        assert plan['explicit_resume'][key] == resume[key]
    assert plan['additional_preserved_prefix_bytes'] == checkpoint['budget']['all_preserved_extra_prefix_bytes']
    assert plan['new_network_bytes'] == checkpoint['budget']['remaining_network_bytes_after_two_reused_bodies_and_source3_resume']
    assert plan['SSD_reservation_bytes'] == 300 * 1024**3
    assert plan['transfer_worker_count'] == 1 and plan['automatic_retry'] is False
    assert not (P / 'manifests/extension_raw_acquisition_admission_v4_5.json').exists()
    assert not (SSD / 'extension_raw_replay_v5/family_ownership.json').exists()
    old_ast = ast.parse((P / 'scripts/52_acquire_extension_raw_sources_v4.py').read_text())
    new_ast = ast.parse(EXECUTOR.read_text())
    names = ['safe_print', 'catchable_termination', 'write_new', 'family_lock', 'cleanup_proof', 'retain_lock_until_gone', 'header_fields', 'safe_state']
    ast_equal = {}
    for name in names:
        a = next(n for n in old_ast.body if isinstance(n, ast.FunctionDef) and n.name == name)
        b = next(n for n in new_ast.body if isinstance(n, ast.FunctionDef) and n.name == name)
        ast_equal[name] = ast.dump(a, include_attributes=False) == ast.dump(b, include_attributes=False)
    assert all(ast_equal.values())
    controls = jread(P / 'logs/extension_acquisition_fault_controls_v4_4.json')
    register(P / 'scripts/53_verify_acquisition_fault_controls_v4.py', controls['test_script_sha256'])
    assert controls['executor_sha256'] == sha(P / 'scripts/52_acquire_extension_raw_sources_v4.py')
    assert controls['all_pass'] is True and len(controls['checks']) == 8 and all(c['pass'] for c in controls['checks'])
    # Actual existing gate: only receipt/header content and source-body stat are
    # accessed. It never hashes/decompresses a GWAS body. The v5 wrapper's body
    # hashing is not invoked; its checkpoint body proof is already sealed.
    gate_path = P / 'scripts/extension_replay_common_v2.py'
    gate = module_from(gate_path, 'independent_exact_common_v2')
    gate_plan = jread(plan['reused_checkpoint_gate_plan'])
    reused = []
    for member in plan['members'][:2]:
        old_member = next(m for m in gate_plan['members'] if m['index'] == member['index'])
        r = gate.acquisition_receipt_gate(gate_plan, old_member, plan['reused_checkpoints'][str(member['index'])]['receipt_sha256'])
        assert r['member'] == member and r['plan_sha256'] == checkpoint['plan_sha256']
        reused.append(member['extension_trait_id'])
    witnesses = [finalization_witness(plan, scenario) for scenario in ['second_success_copy_write_failure', 'post_guard_first_invalidator_write_failure', 'post_guard_failure_addenda_pass']]
    assert all(sha(path) == info['sha256'] for path, info in SEEN.items())
    result = {'status': 'REJECT_V5_FAMILY_FINALIZATION_PERSISTENCE_FAILURE', 'completed_utc': datetime.now(timezone.utc).isoformat(),
              'plan_sha256': EXPECTED_PLAN, 'executor_sha256': sha(EXECUTOR),
              'ordered100_original_and_checkpoint_map_match': True, 'exact2_reused_gate_metadata_pass': reused,
              'new98_target_paths_match': True, 'source3_exact_resume_identity_match': True,
              'global_budget_and_old_prefix_totals_match': True, 'old_cleanup_signal_helpers_AST_identical': ast_equal,
              'prior8_control_receipt_sha256': SEEN[str(P / 'logs/extension_acquisition_fault_controls_v4_4.json')]['sha256'],
              'control_scope': 'Prior eight real helper/signal controls verified by source hashes and exact helper AST identity; not rerun. Three new pure family metadata witnesses, no transfer or source decompression.',
              'family_finalization_witnesses': witnesses,
              'material_findings': [{'id': 'F1', 'priority': 1, 'lines': [503, 516],
                                     'issue': 'Success persistence precedes fallible post-persistence invalidation. Second success-copy write failure escapes before invalidation; invalidator-write failure also escapes with one/two ALL100 success files and no failure oracle.',
                                     'required_correction': 'Preserve prepared v5. Use distinct version with an authoritative durable terminal commit/receipt protocol that cannot admit provisional successes when either persistence or invalidation fails. Cover both replica writes and all invalidator writes under fail-closed terminal ownership. Consumers must require complete terminal proof and veto invalidation; add isolated fixtures for both demonstrated failures.'}],
              'remaining_design_limits': ['No extra GWAS body reads: actual first2 body hashes and source3 partial proof are reused from sealed two-pass checkpoint audit.',
                                         'Global SSD namespace cap is measured before/after source execution and during transfer; final family caps are rechecked after published metadata. Future scratch/metadata and first2 mixed-origin proofs remain explicit.',
                                         'No source/gzip/scientific or completed100 claim follows from preparation or synthetic checkpoints.'],
              'metadata_bindings': SEEN, 'before_after_metadata_hashes_pass': True,
              'root_admission_absent': True, 'execution_clearance': False,
              'network_requests': 0, 'real_transfer_or_native_workers': 0, 'GWAS_body_reads': 0, 'GWAS_decompressions': 0,
              'checker_sha256': sha(Path(__file__)), 'elapsed_seconds': time.monotonic() - start,
              'peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    with OUT.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({key: result[key] for key in ['status', 'ordered100_original_and_checkpoint_map_match', 'exact2_reused_gate_metadata_pass', 'old_cleanup_signal_helpers_AST_identical', 'family_finalization_witnesses', 'elapsed_seconds', 'peak_RSS_bytes']}, indent=2))


if __name__ == '__main__':
    main()
