#!/usr/bin/env python3
"""Independent metadata-only v6 binding audit and tiny lifecycle controls.

Never starts curl, native workers, estimators or a GWAS body reader. Exact
execute()/finalize_family() run with 100 invented receipt entries, mocked
acquisition and identity prerequisites, real tiny JSON/fsync operations and
an independently owned temporary flock. Exact acquire() is exercised only
up to rejection of its post-prefix launch fence, using invented b'abc'.
"""
import ast
from contextlib import contextmanager, redirect_stdout
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
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
PLAN = P / 'manifests/extension_raw_acquisition_plan_v4_6.json'
EXPECTED_PLAN = '393af7ffc3231ceddf829ae669558001a7411bb4a7917b950a713c41546237fe'
EXECUTOR = P / 'scripts/52_acquire_extension_raw_sources_v6.py'
EXPECTED_EXECUTOR = 'd408b245c362c274e95f8c33727ba016367ed14b1e0cd0a77adfd29ae577bf8d'
OUT = REVIEW / 'independent_acquisition_v6_binding_receipt_v4.json'
SEEN = {}


def digest(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            md5.update(block)
            sha.update(block)
    return md5.hexdigest(), sha.hexdigest()


def register(path, expected=None):
    path = Path(path)
    assert not path.name.endswith(('.bgz', '.bgz.partial', '.gz', '.gz.partial')), 'NO_GWAS_BODY_READ_ALLOWED'
    info = {'sha256': digest(path)[1], 'bytes': path.stat().st_size}
    if expected is not None:
        assert info['sha256'] == expected, str(path)
    if str(path) in SEEN:
        assert SEEN[str(path)] == info
    SEEN[str(path)] = info
    return info


def jread(path):
    register(path)
    return json.loads(Path(path).read_text())


def module_from(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def consumer_oracle(primaries, terminal, pending, expected_plan):
    """Separate strict commit oracle, deliberately not the executor's status flag."""
    if pending.exists() or any(Path(str(p)+'.failure.json').exists() for p in [*primaries, terminal]):
        return False
    if not all(p.is_file() for p in [*primaries, terminal]):
        return False
    a, b, s = [json.loads(p.read_text()) for p in [*primaries, terminal]]
    actual = {str(p): digest(p)[1] for p in primaries}
    return (a == b and a['status'] == 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED'
            and a['plan_sha256'] == expected_plan
            and a['completed_source_count'] == len(a['source_receipts']) == 100
            and a['reused_exact_source_count'] == 2 and a['new_exact_source_count'] == 98
            and s['status'] == 'ALL100_FAMILY_TERMINAL_SEAL'
            and s['plan_sha256'] == expected_plan and s['executor_sha256'] == 'e'*64
            and s['primary_receipt_sha256'] == actual
            and len(set(actual.values())) == 1 and s['source_receipts'] == a['source_receipts']
            and s['pending_path'] == str(pending)
            and s['full_gzip_EOF_or_pipeline_replay_certified'] is False)


def family_control(plan, scenario):
    m = module_from(EXECUTOR, 'independent_acquisition_v6_' + scenario)
    writes = []
    with tempfile.TemporaryDirectory(prefix='acquisition_v6_metadata_', dir=SSD/'tmp') as tmp:
        temp = Path(tmp)
        m.PACKAGE, m.FOLDER = temp/'repo', temp/'ssd'
        m.PLAN, m.ADMISSION = temp/'plan.json', temp/'admission.json'
        local = json.loads(json.dumps(plan))
        prefix = temp/'invented-prefix'
        prefix.write_bytes(b'abc')
        local['explicit_resume'].update(original_partial_path=str(prefix), prefix_bytes=3)
        pending, terminal = m.FOLDER/'family_pending_v4_6.json', m.FOLDER/'family_terminal_seal_v4_6.json'
        local.update(pending_path=str(pending), terminal_seal_path=str(terminal))
        m.PLAN.write_text(json.dumps(local))
        m.ADMISSION.write_text(json.dumps({'execution_admitted': True, 'plan_sha256': EXPECTED_PLAN,
                                          'executor_sha256': 'e'*64,
                                          'independent_review_artifact_sha256': {'synthetic.md': 'a'*64, 'synthetic.json': 'b'*64}}))
        original_lock, original_write, original_fsync = m.family_lock, m.write_new, m.fsync_directory
        @contextmanager
        def own_lock(ignored):
            with original_lock(temp/'fixture.lock') as handle:
                yield handle
        m.family_lock = own_lock
        calls = {'snapshot': 0, 'synthetic_acquire': 0, 'synthetic_reuse': 0, 'fsync_directory': 0}
        phase = {'terminal_persisted': False, 'clock': 0.0}
        def snapshot():
            calls['snapshot'] += 1
            bad = scenario == 'post_terminal_floor_and_all_invalidators_fail' and phase['terminal_persisted']
            return {'internal_free_bytes': 0 if bad else 8*1024**3, 'ssd_free_bytes': 500*1024**3}
        monitor = SimpleNamespace(INTERNAL_FLOOR=plan['internal_floor_bytes'], SSD_FLOOR=plan['ssd_floor_bytes'],
                                  RSS_LIMIT=plan['maximum_owned_transfer_rss_bytes'], OUTPUT_LIMIT=plan['monitor_output_limit_bytes'], snapshot=snapshot)
        m.load_monitor = lambda: monitor
        # These prerequisites are explicitly mocked in lifecycle controls.
        # Their actual frozen files are independently checked in main().
        m.physical_mount = lambda: None
        m.assert_bindings = lambda *args: None
        m.global_namespace_gate = lambda *args: 1
        m.subprocess = SimpleNamespace(run=lambda *args, **kwargs: SimpleNamespace(stdout=plan['curl_version']))
        m.time = SimpleNamespace(monotonic=lambda: phase['clock'])
        def hashes(path):
            path = Path(path)
            if path == m.PLAN:
                return '0'*32, EXPECTED_PLAN
            if path == EXECUTOR:
                return '0'*32, 'e'*64
            if path == prefix:
                return '0'*32, local['explicit_resume']['prefix_sha256']
            if path.name.startswith('synthetic-source-'):
                return '0'*32, 'f'*64
            assert path.is_relative_to(temp), 'NO_REAL_BODY_READ_IN_FIXTURE'
            return digest(path)
        m.hashes = hashes
        def checkpoint(member, label):
            calls[label] += 1
            return {'path': str(temp/('synthetic-source-'+str(member['index'])+'.json')), 'sha256': 'f'*64}
        m.reused_checkpoint_gate = lambda plan, member, monitor: checkpoint(member, 'synthetic_reuse')
        m.acquire = lambda plan, member, *args: checkpoint(member, 'synthetic_acquire')
        m.safe_print = lambda *args: None
        primaries = [m.PACKAGE/'logs/extension_raw_acquisition_family_receipt_v4_6.json', m.FOLDER/'extension_raw_acquisition_family_receipt_v4_6.json']
        def write(path, value):
            path = Path(path)
            writes.append(str(path.relative_to(temp)))
            if scenario in ('second_copy_and_all_invalidators_fail', 'post_terminal_floor_and_all_invalidators_fail', 'terminal_write_and_all_invalidators_fail') and path.name.endswith('.failure.json'):
                raise OSError('INDEPENDENT_ALL_INVALIDATOR_WRITES_FAIL')
            if scenario == 'second_copy_and_all_invalidators_fail' and path == primaries[1]:
                raise OSError('INDEPENDENT_SECOND_PRIMARY_COPY_WRITE_FAIL')
            if scenario == 'terminal_write_and_all_invalidators_fail' and path == terminal:
                raise OSError('INDEPENDENT_TERMINAL_SEAL_WRITE_FAIL')
            original_write(path, value)
            if path == terminal:
                phase['terminal_persisted'] = True
                if scenario == 'deferred_signal_after_terminal':
                    m.catchable_termination(15, None)
                if scenario == 'primary_tamper_after_terminal':
                    with primaries[0].open('a') as f:
                        f.write('\n')
                if scenario == 'deadline_after_terminal':
                    phase['clock'] = plan['family_seconds_limit'] + 1.0
        m.write_new = write
        def fsync_dir(path):
            calls['fsync_directory'] += 1
            if scenario == 'pending_directory_fsync_failure' and calls['fsync_directory'] == 1:
                raise OSError('INDEPENDENT_PENDING_DIRECTORY_FSYNC_FAIL')
            return original_fsync(path)
        m.fsync_directory = fsync_dir
        if scenario == 'pending_unlink_failure':
            class UnlinkFaultPath(type(Path())):
                def unlink(self, *args, **kwargs):
                    if str(self) == str(pending):
                        raise OSError('INDEPENDENT_PENDING_UNLINK_FAIL')
                    return super().unlink(*args, **kwargs)
            m.Path = UnlinkFaultPath
        exception = None
        try:
            with redirect_stdout(io.StringIO()):
                m.execute(EXPECTED_PLAN)
        except BaseException as exc:
            exception = type(exc).__name__ + ': ' + str(exc)
        expected_success = scenario in ('healthy_commit', 'duplicate_after_healthy_commit')
        admitted = consumer_oracle(primaries, terminal, pending, EXPECTED_PLAN)
        assert admitted == expected_success
        assert (exception is None) == expected_success
        if scenario == 'pending_directory_fsync_failure':
            assert calls['synthetic_acquire'] == 0 and not (m.FOLDER/'family_ownership.json').exists()
        else:
            assert calls['synthetic_acquire'] == 98 and calls['synthetic_reuse'] == 4
        assert pending.exists() != expected_success
        duplicate = None
        if scenario == 'duplicate_after_healthy_commit':
            before = {str(p): digest(p)[1] for p in [*primaries, terminal]}
            transfer_before = calls['synthetic_acquire']
            try:
                m.execute(EXPECTED_PLAN)
            except BaseException as exc:
                duplicate = type(exc).__name__ + ': ' + str(exc)
            assert duplicate and 'PRIOR_ATTEMPT_REQUIRES_NEW_VERSION_AND_EXPLICIT_AUDIT' in duplicate
            assert calls['synthetic_acquire'] == transfer_before
            assert {str(p): digest(p)[1] for p in [*primaries, terminal]} == before
        paths = {'primary_repo': primaries[0], 'primary_ssd': primaries[1], 'terminal': terminal, 'pending': pending}
        observed = {k: {'exists': p.exists(), 'status': json.loads(p.read_text())['status'] if p.exists() else None,
                        'failure_marker': Path(str(p)+'.failure.json').exists()} for k,p in paths.items()}
        return {'scenario': scenario, 'pass': True, 'strict_consumer_admitted': admitted,
                'exception': exception, 'duplicate_exception': duplicate, 'calls': calls,
                'observed_artifacts': observed, 'write_order': writes,
                'real_network_or_native_workers': 0, 'real_GWAS_body_reads': 0, 'invented_prefix_bytes': 3}


def prefix_fence_control(plan, scenario):
    m = module_from(EXECUTOR, 'independent_v6_prefix_' + scenario)
    with tempfile.TemporaryDirectory(prefix='acquisition_v6_prefix_', dir=SSD/'tmp') as tmp:
        temp = Path(tmp)
        m.PACKAGE, m.FOLDER = temp/'repo', temp/'ssd'
        for child in ('raw','logs','receipts'):
            (m.FOLDER/child).mkdir(parents=True)
        prefix = temp/'invented-prefix'
        prefix.write_bytes(b'abc')
        local = json.loads(json.dumps(plan))
        local['explicit_resume'].update(original_partial_path=str(prefix), prefix_bytes=3, prefix_sha256=digest(prefix)[1])
        member = dict(plan['members'][2], body_path=str(m.FOLDER/'raw'/'invented-output'))
        calls = {'snapshot': 0, 'binding': 0, 'Popen': 0}
        def snapshot():
            calls['snapshot'] += 1
            bad = scenario == 'post_copy_floor_failure' and calls['snapshot'] == 2
            return {'internal_free_bytes': 0 if bad else 8*1024**3, 'ssd_free_bytes': 500*1024**3}
        monitor = SimpleNamespace(snapshot=snapshot)
        def binding(*args):
            calls['binding'] += 1
            if scenario == 'post_copy_binding_failure' and calls['binding'] == 2:
                raise AssertionError('INDEPENDENT_POST_COPY_BINDING_CHANGED')
        def impossible_popen(*args, **kwargs):
            calls['Popen'] += 1
            raise AssertionError('REAL_WORKER_FORBIDDEN')
        def hashes(path):
            path = Path(path)
            assert path.is_relative_to(temp), 'NO_REAL_BODY_READ_IN_PREFIX_CONTROL'
            return digest(path)
        m.assert_bindings = binding
        m.physical_mount = lambda: None
        m.global_namespace_gate = lambda *args: 1
        m.subprocess = SimpleNamespace(Popen=impossible_popen, STDOUT=-2)
        m.hashes = hashes
        exception = None
        try:
            with redirect_stdout(io.StringIO()):
                m.acquire(local, member, monitor, EXPECTED_PLAN, time.monotonic())
        except BaseException as exc:
            exception = type(exc).__name__ + ': ' + str(exc)
        partial = Path(member['body_path']+'.partial')
        receipt = json.loads((m.FOLDER/'receipts'/(member['extension_trait_id']+'.json')).read_text())
        assert exception and calls['Popen'] == 0
        assert partial.read_bytes() == b'abc'
        assert receipt['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        assert receipt['teardown'] == {'worker_launched': False, 'remaining_group_members': [], 'teardown_verified': True}
        return {'scenario': scenario, 'pass': True, 'calls': calls, 'exception': exception,
                'source_receipt_status': receipt['status'], 'retained_invented_prefix_bytes': partial.stat().st_size,
                'real_network_or_native_workers': 0, 'real_GWAS_body_reads': 0}


def main():
    start = time.monotonic()
    assert not OUT.exists()
    assert shutil.disk_usage('/System/Volumes/Data').free >= 256*1024**2
    register(EXECUTOR, EXPECTED_EXECUTOR)
    plan = jread(PLAN)
    assert SEEN[str(PLAN)]['sha256'] == EXPECTED_PLAN
    register(SSD/'manifests'/PLAN.name, EXPECTED_PLAN)
    for path, expected in plan['bound_sources'].items():
        register(path, expected)
    checkpoint = jread(REVIEW/'independent_acquisition_checkpoint_receipt_v4_4.json')
    old = jread(P/'manifests/extension_raw_acquisition_plan_v4_4.json')
    v5 = jread(P/'manifests/extension_raw_acquisition_plan_v4_5.json')
    rejected = jread(REVIEW/'independent_acquisition_v5_binding_receipt_v4.json')
    assert rejected['status'] == 'REJECT_V5_FAMILY_FINALIZATION_PERSISTENCE_FAILURE'
    assert len(plan['members']) == len(checkpoint['proposed_mixed_origin_members']) == 100
    assert [m['index'] for m in plan['members']] == list(range(1,101))
    assert len({m['extension_trait_id'] for m in plan['members']}) == len({m['filename'] for m in plan['members']}) == 100
    for original, member in zip(old['members'], plan['members']):
        assert all(member[k] == v for k,v in original.items() if k != 'body_path')
        expected_body = original['body_path'] if member['index'] <= 2 else str(SSD/'extension_raw_replay_v6/raw'/member['filename'])
        assert member['body_path'] == expected_body
    assert plan['reused_checkpoints'] == v5['reused_checkpoints'] and set(plan['reused_checkpoints']) == {'1','2'}
    assert plan['network_sources_to_transfer'] == 98 and plan['explicit_resume'] == v5['explicit_resume']
    resume = checkpoint['source3_resume']
    for key in ['source_trait','original_partial_path','prefix_bytes','prefix_sha256']:
        assert plan['explicit_resume'][key] == resume[key]
    assert plan['additional_preserved_prefix_bytes'] == checkpoint['budget']['all_preserved_extra_prefix_bytes'] == 2405433344
    assert plan['new_network_bytes'] == checkpoint['budget']['remaining_network_bytes_after_two_reused_bodies_and_source3_resume'] == 221795203885
    assert sum(m['expected_bytes'] for m in plan['members']) == plan['compressed_network_bytes'] == 227610388647
    for key in ['SSD_reservation_bytes','internal_floor_bytes','ssd_floor_bytes','maximum_owned_transfer_rss_bytes','per_body_seconds_limit','family_seconds_limit','transfer_worker_count','automatic_retry','runtime_poll_seconds','curl_version','monitor_path','curl_path','monitor_output_limit_bytes']:
        assert plan[key] == v5[key]
    assert plan['SSD_reservation_bytes'] == 300*1024**3 and plan['transfer_worker_count'] == 1 and plan['automatic_retry'] is False
    ledger = jread(plan['global_resource_ledger_path'])
    assert SEEN[plan['global_resource_ledger_path']]['sha256'] == plan['global_resource_ledger_sha256']
    assert sum(ledger['component_bytes'].values()) == ledger['reserved_total_bytes'] == 318518274622
    assert ledger['ceiling_bytes'] == plan['SSD_reservation_bytes']
    assert ledger['unallocated_margin_bytes'] == ledger['ceiling_bytes'] - ledger['reserved_total_bytes'] == 3604272578
    assert ledger['component_bytes']['raw_bytes'] == plan['compressed_network_bytes']
    assert ledger['component_bytes']['preserved_prefix_bytes'] == plan['additional_preserved_prefix_bytes']
    assert ledger['includes_AppleDouble_transport_sidecars_in_actual_meter'] is True
    assert not (P/'manifests/extension_raw_acquisition_admission_v4_6.json').exists()
    assert not (SSD/'extension_raw_replay_v6/family_ownership.json').exists()
    assert not Path(plan['pending_path']).exists() and not Path(plan['terminal_seal_path']).exists()
    durability = jread(P/'logs/raw_family_directory_durability_control_v1.json')
    assert durability['directory_fsync'] == 'PASS'
    old_ast, new_ast = [ast.parse(p.read_text()) for p in [P/'scripts/52_acquire_extension_raw_sources_v4.py', EXECUTOR]]
    names = ['safe_print','catchable_termination','write_new','family_lock','cleanup_proof','retain_lock_until_gone','header_fields','safe_state']
    ast_equal = {}
    for name in names:
        a = next(n for n in old_ast.body if isinstance(n,ast.FunctionDef) and n.name == name)
        b = next(n for n in new_ast.body if isinstance(n,ast.FunctionDef) and n.name == name)
        ast_equal[name] = ast.dump(a,include_attributes=False) == ast.dump(b,include_attributes=False)
    assert all(ast_equal.values())
    controls = jread(P/'logs/extension_acquisition_fault_controls_v4_4.json')
    register(P/'scripts/53_verify_acquisition_fault_controls_v4.py', controls['test_script_sha256'])
    assert controls['executor_sha256'] == digest(P/'scripts/52_acquire_extension_raw_sources_v4.py')[1]
    assert controls['all_pass'] is True and len(controls['checks']) == 8 and all(c['pass'] for c in controls['checks'])
    gate = module_from(P/'scripts/extension_replay_common_v2.py', 'independent_v6_unchanged_gate')
    gate_plan = jread(plan['reused_checkpoint_gate_plan'])
    reused = []
    for member in plan['members'][:2]:
        original = next(m for m in gate_plan['members'] if m['index'] == member['index'])
        r = gate.acquisition_receipt_gate(gate_plan,original,plan['reused_checkpoints'][str(member['index'])]['receipt_sha256'])
        assert r['member'] == member and r['plan_sha256'] == checkpoint['plan_sha256']
        reused.append(member['extension_trait_id'])
    scenarios = ['second_copy_and_all_invalidators_fail','post_terminal_floor_and_all_invalidators_fail',
                 'terminal_write_and_all_invalidators_fail','pending_unlink_failure','pending_directory_fsync_failure',
                 'deferred_signal_after_terminal','primary_tamper_after_terminal','deadline_after_terminal',
                 'healthy_commit','duplicate_after_healthy_commit']
    family = [family_control(plan,s) for s in scenarios]
    prefix = [prefix_fence_control(plan,s) for s in ['post_copy_binding_failure','post_copy_floor_failure']]
    assert all(digest(path)[1] == info['sha256'] and Path(path).stat().st_size == info['bytes'] for path,info in SEEN.items())
    result = {'status': 'PASS_V6_BINDING_AND_BOUNDED_LIFECYCLE_PREFLIGHT',
              'completed_utc': datetime.now(timezone.utc).isoformat(), 'plan_sha256': EXPECTED_PLAN,
              'executor_sha256': EXPECTED_EXECUTOR, 'bound_dependency_count': len(plan['bound_sources']),
              'original100_identities_unchanged': True, 'reused2_unchanged_metadata_gate_pass': reused,
              'new98_private_v6_paths_match': True, 'source3_exact_resume_identity_match': True,
              'old_checkpoint_body_hashes_reused_without_new_body_reads': True,
              'global_ledger': {'sha256': plan['global_resource_ledger_sha256'], 'reserved_total_bytes': ledger['reserved_total_bytes'],
                                'ceiling_bytes': ledger['ceiling_bytes'], 'unallocated_margin_bytes': ledger['unallocated_margin_bytes'],
                                'both_preserved_prefix_bytes': plan['additional_preserved_prefix_bytes']},
              'old_cleanup_signal_helpers_AST_identical': ast_equal,
              'prior8_real_helper_control_receipt_sha256': SEEN[str(P/'logs/extension_acquisition_fault_controls_v4_4.json')]['sha256'],
              'durability_control_sha256': SEEN[str(P/'logs/raw_family_directory_durability_control_v1.json')]['sha256'],
              'family_lifecycle_controls': family, 'post_prefix_launch_fence_controls': prefix,
              'control_count': len(family)+len(prefix), 'all_controls_pass': True, 'material_blockers': [],
              'required_consumer_contract': ['Require exact ALL100 terminal seal and both primary copies with their bound current hashes.',
                                             'Require 100 ordered original source receipts, exact plan/executor identities, reused2/new98 mapping.',
                                             'Reject any PENDING or primary/terminal failure marker; never admit an ALL100 primary alone.'],
              'qualifications': ['Prepared executor/plan review only. Root admission absent; this review launches no transfer.',
                                 'No additional GWAS body reads, gzip/EOF validation, source pipeline replay or native/statistical inference.',
                                 'Ten new real tiny-file lifecycle controls mock acquisition/identity prerequisites. Two new prefix controls copy only invented three-byte content; no Popen occurs.',
                                 'Eight prior real helper/signal controls apply by pinned receipt hashes and exact helper AST identity; they were not rerun.',
                                 'Runtime consumers must implement the complete terminal contract before accepting future receipts; future consumer code is not certified here.',
                                 'Final marker unlink is last commit operation. A crash restoring the unfsynced unlink conservatively rejects completion; no power-loss simulation was performed.',
                                 '300 GiB actual namespace ceiling and 5 GiB physical SSD free floor are distinct guards; the reservation ledger includes both preserved prefixes and planned namespaces.'],
              'resource': {'elapsed_seconds': time.monotonic()-start, 'peak_self_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                           'internal_free_bytes_at_finish': shutil.disk_usage('/System/Volumes/Data').free,
                           'SSD_free_bytes_at_finish': shutil.disk_usage(SSD).free,
                           'GWAS_body_reads': 0, 'real_curl_or_native_worker_count': 0},
              'metadata_before_after': SEEN}
    with OUT.open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'status': result['status'], 'controls': result['control_count'], 'bound_dependencies': len(plan['bound_sources']),
                      'metadata_files': len(SEEN), 'elapsed_seconds': result['resource']['elapsed_seconds'],
                      'peak_self_RSS_bytes': result['resource']['peak_self_RSS_bytes']}))


if __name__ == '__main__':
    main()
