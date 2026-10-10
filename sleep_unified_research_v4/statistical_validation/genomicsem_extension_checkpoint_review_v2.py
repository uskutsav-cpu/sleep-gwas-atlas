#!/usr/bin/env python3
"""Independent metadata and tiny-fixture review; never read GWAS inputs or launch workers."""
import ast
import copy
import csv
import datetime
import fcntl
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import signal
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
SCRIPTS = P / 'scripts'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT = SSD / 'extension_pipeline_replay_v2'
PLAN = OUT / 'extension_pipeline_replay_plan_v2.json'
PLAN_SHA = '77f008725b8d87ebe5e08c85f8667c8b55f02a7065ff3a6e46a7c45d3fd38cb3'
CONTROL_SHA = 'b666e7626eb26f8b338340bdc4948bc286e85c37db3d68831d8a3a3ce2c5d2e0'
RESULT = P / 'statistical_validation/genomicsem_extension_checkpoint_review_controls_v2.json'
sys.dont_write_bytecode = True
sys.path.insert(0, str(SCRIPTS))
import extension_replay_common_v2 as common


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    item = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(item)
    return item


def main():
    checks = []
    observations = []
    def check(value, label):
        if not value:
            raise RuntimeError('INDEPENDENT_FIXTURE_CHECK_FAILED: ' + label)
        checks.append(label)
    check(common.sha(PLAN) == PLAN_SHA, 'exact_reviewed_v2_plan')
    check(common.sha(OUT / 'preparation_controls_v2.json') == CONTROL_SHA, 'exact_preparer553_control_receipt')
    plan = json.loads(PLAN.read_text())
    forbidden = set(plan['archived_input_sha256']) | set(plan['baseline_input_sha256']) | {m['raw'] for m in plan['members']}
    admitted_metadata = {p: h for p, h in plan['dependencies_sha256'].items() if p not in forbidden}
    for path, digest in admitted_metadata.items():
        check(common.sha(path) == digest, 'non_GWAS_dependency_before ' + path)
    check(len(plan['members']) == len({m['extension_trait_id'] for m in plan['members']}) == len({m['raw'] for m in plan['members']}) == 100, '100_unique_members_and_raw_paths')
    acq = json.loads(Path(plan['acquisition_plan']).read_text())
    historical = {m['extension_trait_id']: m for m in acq['members']}
    with Path(plan['panel']).open() as f:
        panel = list(csv.DictReader(f, delimiter='\t'))
    check([m['extension_trait_id'] for m in plan['members']] == [r['extension_trait_id'] for r in panel] == [m['extension_trait_id'] for m in acq['members']], 'original_panel_acquisition_replay_order_identical')
    check([m['index'] for m in plan['members']] == list(range(1, 101)), 'original_indices1_through100_no_skips')
    panel_by_id = {r['extension_trait_id']: r for r in panel}
    for m in plan['members']:
        t = m['extension_trait_id']; a = historical[t]; row = panel_by_id[t]
        check(m['acquisition_member'] == a, 'exact_acquisition_member ' + t)
        check((m['raw'], m['raw_sha256'], m['raw_md5'], m['raw_bytes']) == (a['body_path'], a['expected_sha256'], a['expected_md5'], a['expected_bytes']), 'exact_source_tuple ' + t)
        old = json.loads(Path(a['historical_streaming_receipt']).read_text()); v = old['source_verification']
        check(old['pipeline_status'] == 'STREAM_HARMONIZE_MUNGE_PASS' and old['extension_trait_id'] == t and old['source_url'] == v['source'] == a['url'] and v['verification_status'] == 'PASS', 'exact_original_streaming_identity ' + t)
        check((v['observed_sha256'], v['observed_md5'], v['observed_size_bytes'], v['http_version_id'], v['http_etag'].strip('"')) == (a['expected_sha256'], a['expected_md5'], a['expected_bytes'], a['expected_s3_version_id'], a['expected_etag']), 'original_SHA_MD5_bytes_version_ETag ' + t)
        check(old['reference_sha256'] == plan['reference_sha256'] and old['munged_output_sha256'] == m['archived_munged_sha256'] == plan['archived_input_sha256'][m['archived_munged']], 'original_reference_and_processed_receipt_metadata ' + t)
        expected = [plan['python'], '-u', '-B', str(ROOT/'discovery_extension/scripts/10_harmonize_panukbb.py'), '--trait-id', t, '--source', m['raw'], '--panel', plan['panel'], '--reference', plan['reference'], '--out', m['harmonized'], '--qc-out', m['harmonization_qc'], '--receipt-out', m['harmonization_receipt']]
        check(m['harmonize_command'] == expected, 'unchanged_local_harmonization_argv ' + t)
        expected = [plan['python'], '-u', '-B', str(Path(plan['ldsc_dir'])/'munge_sumstats.py'), '--sumstats', m['harmonized'], '--merge-alleles', plan['w_hm3'], '--snp', 'SNP', '--a1', 'A1', '--a2', 'A2', '--frq', 'FRQ', '--p', 'P', '--N-col', 'N', '--signed-sumstats', 'BETA,0', '--chunksize', '500000', '--out', m['munged_prefix']]
        check(m['munge_command'] == expected, 'unchanged_stock_munge_argv ' + t)
        sample_n = 4.0 / (1.0 / int(row['cases']) + 1.0 / int(row['controls'])) if row['binary_or_continuous'] == 'binary' else float(row['sample_size'])
        check(math.isfinite(sample_n) and sample_n > 0 and row['ancestry'] == a['ancestry'] == 'EUR' and row['build'] == a['build'] == 'GRCh37', 'metadata_N_ancestry_build_eligible ' + t)
    original = json.loads(Path(plan['preserved_v1_plan']).read_text())
    old_namespace = str(Path(plan['preserved_v1_plan']).parent)
    expected_members = copy.deepcopy(original['members'])
    for m in expected_members:
        for k, v in list(m.items()):
            if isinstance(v, str) and v.startswith(old_namespace + '/'):
                m[k] = v.replace(old_namespace+'/', str(OUT)+'/', 1)
            elif k in ['harmonize_command', 'munge_command']:
                m[k] = [x.replace(old_namespace+'/', str(OUT)+'/', 1) for x in v]
    check(plan['members'] == expected_members, 'v1_v2_member_change_only_output_namespace')
    for k in ['reference', 'reference_sha256', 'w_hm3', 'template_rows', 'environment', 'archived_input_sha256', 'original_190_jobs', 'baseline_execution_plan_sha256', 'baseline_input_sha256', 'reservation_arithmetic']:
        check(plan[k] == original[k], 'v1_v2_frozen_science_and_guards ' + k)
    check(len(plan['original_190_jobs']) == 190, 'original190_prerequisite_is_complete_declared_job_family')
    check((plan['heavy_preprocessing_commands'], plan['owned_stdlib_validation_commands'], plan['total_owned_commands'], plan['estimator_calls']) == (200, 200, 400, 0), '400_owned_commands_zero_native_estimators')
    for n in ['46_prepare_extension_pipeline_replay_v2.py', '47_run_extension_pipeline_replay_v2.py', 'extension_replay_common_v2.py', 'extension_replay_validate_v2.py', 'extension_replay_fault_controls_v2.py']:
        ast.parse((SCRIPTS/n).read_text(), filename=n)
        check(True, 'independent_AST ' + n)
    runner = module('_independent_extension_checkpoint_runner_v2', SCRIPTS/'47_run_extension_pipeline_replay_v2.py')
    protected = runner.load_protected()
    streaming = module('_independent_original_streaming_reader', ROOT/'discovery_extension/scripts/streaming_io.py')
    with tempfile.TemporaryDirectory(prefix='genomicsem_extension_review_', dir=P/'statistical_validation') as d:
        tmp = Path(d)
        def gz(name, payload, stamp=0):
            p = tmp/name
            with p.open('xb') as f:
                with gzip.GzipFile(fileobj=f, filename=name, mode='wb', mtime=stamp) as z:
                    z.write(payload)
            return p
        payload = b'SNP\tA1\tA2\tZ\tN\nCONTROL_A\tA\tG\t1.000\t100.000\nCONTROL_B\tC\tT\t\t\n'
        archive = gz('archive_fixture.gz', payload); same = gz('new_fixture.gz', payload, 1)
        r = common.compare_munged(same, archive, 2)
        check(r['status'] == 'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' and r['compressed_sha256_new'] != r['compressed_sha256_archived'], 'decompressed_identity_not_gzip_container_identity')
        check(r['counts_new'] == {'finite_N_Z': 1, 'missing_or_nonfinite_N_Z': 1} and r['both_gzip_CRC_and_EOF_verified'], 'finite_missingness_counts_and_EOF_proof')
        changes = {'SNP': payload.replace(b'CONTROL_A', b'CONTROL_X'), 'A1': payload.replace(b'\tA\tG\t', b'\tG\tG\t'), 'A2': payload.replace(b'\tA\tG\t', b'\tA\tC\t'), 'Z': payload.replace(b'1.000', b'2.000'), 'N': payload.replace(b'100.000', b'101.000'), 'numeric_alias': payload.replace(b'1.000', b'1.0'), 'missingness_alias': payload.replace(b'\t\t\n', b'\tNA\tNA\n'), 'order': b''.join([payload.splitlines(keepends=True)[0], *reversed(payload.splitlines(keepends=True)[1:])]), 'extra_row': payload+b'CONTROL_C\tA\tG\t0.000\t100.000\n', 'missing_row': b''.join(payload.splitlines(keepends=True)[:2])}
        for label, changed in changes.items():
            item = gz(label+'.gz', changed)
            rr = common.compare_munged(item, archive, 2)
            check(rr['status'] == 'DECOMPRESSED_TEMPLATE_MISMATCH_PRESERVED' and rr['both_gzip_CRC_and_EOF_verified'], 'entire_template_rejects_'+label)
        changed = gz('early_mismatch_bad_tail.gz', changes['SNP'])
        damaged = bytearray(changed.read_bytes()); damaged[-8] ^= 1; changed.write_bytes(damaged)
        truncated = tmp/'truncated_fixture.gz'; truncated.write_bytes(same.read_bytes()[:-4])
        for label, p in [('CRC_after_early_mismatch', changed), ('truncated_EOF', truncated)]:
            try:
                common.compare_munged(p, archive, 2)
            except (OSError, EOFError):
                check(True, 'reject_'+label+'_no_false_EOF_pass')
            else:
                raise RuntimeError('GZIP_DAMAGE_WAS_ADMITTED')
        missing = gz('all_missing_fixture.gz', b'SNP\tA1\tA2\tZ\tN\nCONTROL_A\t\t\t\t\n')
        rr = common.compare_munged(missing, missing, 1)
        check(rr['status'] == 'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' and rr['counts_new']['finite_N_Z'] == 0, 'fixture_documents_template_identity_does_not_assert_finite_fit_universe')
        with streaming.open_verified_gzip_text(local_path=same, expected_md5='md5:'+hashlib.md5(same.read_bytes()).hexdigest(), expected_size_bytes=same.stat().st_size) as (f, proof):
            check(f.read().encode() == payload, 'original_streaming_reader_entire_tiny_fixture')
        check(proof['verification_status'] == 'PASS' and proof['observed_sha256'] == common.sha(same), 'original_reader_full_compressed_hash_after_EOF')
        for label, kw in [('MD5', dict(expected_md5='md5:'+'0'*32, expected_size_bytes=same.stat().st_size)), ('size', dict(expected_md5='md5:'+hashlib.md5(same.read_bytes()).hexdigest(), expected_size_bytes=same.stat().st_size+1))]:
            try:
                with streaming.open_verified_gzip_text(local_path=same, **kw) as (f, _):
                    f.read()
            except ValueError:
                check(True, 'original_reader_rejects_'+label)
            else:
                raise RuntimeError('ORIGINAL_READER_IDENTITY_FAILURE_ADMITTED')
        for label in ['healthy', 'late_signal', 'receipt_failure', 'post_persistence_resource_failure']:
            protected.TERMINATION_REQUEST.clear(); events=[]; late=[False]
            target = tmp/(label+'.json'); master={'status':'ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS'}
            def limits(*_):
                events.append('resource_query_late='+str(late[0]))
                return {'fixture':True}, ('SSD_SPACE_GUARD' if late[0] and label=='post_persistence_resource_failure' else None)
            protected.limits = limits
            def writer(path, value):
                if Path(path)==target and label=='receipt_failure':
                    raise OSError('CONTROLLED_FINAL_RECEIPT_IO_FAILURE')
                common.write_new(path, value);events.append('persisted')
                if Path(path)==target:
                    late[0]=True
                    if label=='late_signal':
                        protected.catchable_termination(signal.SIGTERM, None)
            with patch.object(runner, 'write_new', writer):
                runner.seal_master(protected, master, target, {}, 0, [True, []])
            if label in ['late_signal','receipt_failure']:
                check(master['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and Path(str(target)+'.failure.json').exists(), 'final_seal_fails_and_preserves_'+label)
            else:
                check(master['status'] == 'ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS', 'final_seal_observed_status_'+label)
            if label=='post_persistence_resource_failure':
                check(events == ['resource_query_late=False','persisted'], 'witness_no_resource_query_after_terminal_receipt_persistence')
                observations.append({'finding':'terminal_receipt_persistence_resource_boundary', 'events':events, 'observed_status':master['status'], 'injected_condition':'SSD floor crossed inside final receipt writer; no real floor changed'})
        protected.TERMINATION_REQUEST.clear()
        # Call the original scheduling function only with all source/worker,
        # acquisition, filesystem-resource and baseline gates replaced by labelled
        # tiny-fixture mocks. subprocess creation is forbidden explicitly.
        # This checks terminal control flow, not science or actual admission.
        schedule = tmp/'scheduler_fixture';schedule.mkdir()
        for name in ['receipts','receipts_v4','harmonized','munged','qc','logs_v4','tmp','cache']:
            (schedule/name).mkdir()
        members=[]
        for i in range(1,101):
            t='CONTROL_ONLY_'+str(i)
            members.append(dict(index=i, extension_trait_id=t, acquisition_receipt=str(schedule/(t+'.acquired.json')), harmonized=str(schedule/'harmonized'/(t+'.txt.gz')), harmonization_qc=str(schedule/'qc'/(t+'.tsv')), harmonization_receipt=str(schedule/'receipts'/(t+'.harm.json')), munged_prefix=str(schedule/'munged'/t), munged=str(schedule/'munged'/(t+'.sumstats.gz')), source_gate_receipt=str(schedule/'receipts'/(t+'.source.json')), comparison_receipt=str(schedule/'receipts'/(t+'.compare.json')), archived_munged_sha256='FIXTURE_ONLY_NO_ARCHIVE', harmonize_command=['FIXTURE_ONLY_HARMONIZER_NEVER_EXECUTED'], munge_command=['FIXTURE_ONLY_MUNGE_NEVER_EXECUTED']))
        family=schedule/'acquisition_family_fixture.json'
        fp=dict(guard=copy.deepcopy(plan['guard']), members=members, python='FIXTURE_ONLY_NO_NATIVE_RUNTIME', acquisition_family_receipt=str(family))
        fixture_plan=schedule/'plan.json';common.write_new(fixture_plan, fp)
        fixture_admission=schedule/'admission.json';common.write_new(fixture_admission, {'fixture_only':True})
        args=SimpleNamespace(plan=fixture_plan, plan_sha=common.sha(fixture_plan), admission=fixture_admission, admission_sha=common.sha(fixture_admission))
        calls=[]
        class FakeProtected:
            TERMINATION_REQUEST=[]
            catchable_termination=staticmethod(lambda *_:None)
            check_termination=staticmethod(lambda *_:None)
            safe_print=staticmethod(lambda *_a, **_k:None)
            baseline_gate=staticmethod(lambda *_: {'fixture_only':True})
            stage_resource_gate=staticmethod(lambda *_: {'fixture_only':True})
            await_owned_cleanup=staticmethod(lambda *_: (_ for _ in ()).throw(RuntimeError('NO_OWNED_FIXTURE_PROCESSES_EXPECTED')))
            @staticmethod
            def worker(command, prefix, *_args):
                calls.append(command)
                if '--mode' in command:
                    mode=command[command.index('--mode')+1]
                    common.write_new(prefix, {'status':'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' if mode=='compare' else 'EXACT_RAW_SOURCE_GATE_PASS', 'fixture_only':True})
                    if len(calls)==400:
                        common.write_new(family, {'status':'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT', 'completed_source_count':100, 'stop_reason':'FIXTURE_FINAL_FAMILY_GATE_FAILURE'})
                elif command==['FIXTURE_ONLY_HARMONIZER_NEVER_EXECUTED']:
                    prefix.write_bytes(b'FIXTURE_ONLY')
                elif command==['FIXTURE_ONLY_MUNGE_NEVER_EXECUTED']:
                    Path(str(prefix)+'.sumstats.gz').write_bytes(b'FIXTURE_ONLY')
                else:
                    raise RuntimeError('UNRECOGNIZED_FAKE_WORKER')
        def fake_wait(_plan, m, _a, _admission, _protected, _started):
            return dict(plan_sha256=args.plan_sha, pipeline_root_admission_sha256=args.admission_sha, extension_trait_id=m['extension_trait_id'], acquisition_receipt=m['acquisition_receipt'], acquisition_receipt_sha256=hashlib.sha256(m['extension_trait_id'].encode()).hexdigest(), fixture_only=True)
        patches=[patch.object(runner,'OUT',schedule),patch.object(runner,'load_protected',lambda:FakeProtected()),patch.object(runner,'admission_gate',lambda *_:None),patch.object(runner,'check_bindings',lambda *_a,**_k:None),patch.object(runner,'physical_mount',lambda:{'fixture_only':True}),patch.object(runner,'wait_checkpoint',fake_wait),patch.object(runner,'checkpoint_binding_gate',lambda *_:None),patch.object(runner,'acquire_source_mutex',lambda *_:os.open(schedule/'mock_lock_descriptor',os.O_RDWR|os.O_CREAT,0o600)),patch('subprocess.Popen',side_effect=RuntimeError('FORBIDDEN_REAL_WORKER_LAUNCH'))]
        entered=[]
        try:
            for p in patches:p.__enter__();entered.append(p)
            runner.run(args)
        finally:
            for p in reversed(entered):p.__exit__(None,None,None)
        terminal=json.loads((schedule/'pipeline_execution_receipt_v2.json').read_text())
        check(len(calls)==400 and terminal['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS' and json.loads(family.read_text())['status']=='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT', 'witness_late_failed_acquisition_family_not_rechecked_in_final_scheduler_path')
        observations.append({'finding':'late_acquisition_family_failure', 'fixture_dispatches':len(calls), 'real_worker_dispatches':0, 'source_family_status':json.loads(family.read_text())['status'], 'pipeline_status':terminal['status'], 'source_family_failure_published':'inside final mock compare after final source discovery', 'scope':'all scientific/admission gates deliberately mocked; exact original scheduler terminal control flow'})
        # Real isolated flock, no shared production mutex, no process launch.
        fd=os.open(tmp/'isolated_flock_fixture',os.O_RDWR|os.O_CREAT,0o600)
        fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ownership=[False,[SimpleNamespace(pid=777777, returncode=0)]]; attempts=[]
        class FakeMonitor:
            @staticmethod
            def terminate_owned(_proc):
                other=os.open(tmp/'isolated_flock_fixture',os.O_RDWR)
                try:
                    try:fcntl.flock(other,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    except BlockingIOError:attempts.append('lock_retained')
                    else:raise RuntimeError('ISOLATED_LOCK_RELEASED_TOO_EARLY')
                finally:os.close(other)
                if len(attempts)==1:raise RuntimeError('CONTROLLED_FIRST_TEARDOWN_QUERY_FAILURE')
                return dict(initial_group_members=[],remaining_group_members=[],signals=[])
        with patch.object(protected,'safe_print',side_effect=BrokenPipeError('CONTROLLED_DIAGNOSTIC_FAILURE')),patch.object(protected.time,'sleep',lambda _:None),patch.object(protected,'SSD',tmp),patch.object(protected,'save',lambda *_:None):
            rr=protected.await_owned_cleanup(FakeMonitor(), ownership, PLAN_SHA)
        check(ownership==[True,[]] and attempts==['lock_retained','lock_retained'] and rr==[], 'shared_cleanup_helper_retains_isolated_real_flock_across_teardown_and_diagnostic_failure')
        os.close(fd)
    for path,digest in admitted_metadata.items():
        check(common.sha(path)==digest, 'non_GWAS_dependency_after '+path)
    check(common.sha(PLAN)==PLAN_SHA and common.sha(OUT/'preparation_controls_v2.json')==CONTROL_SHA, 'prepared_plan_and_preparer_controls_unchanged')
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':peak*=1024
    check(peak<128*(1<<20), 'review_self_peak_RSS_below128MiB')
    record=dict(status='INDEPENDENT_METADATA_AND_TINY_FIXTURES_COMPLETE_WITH_TERMINAL_FINDINGS', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), plan_sha256=PLAN_SHA, preparer_control_sha256=CONTROL_SHA, verifier_sha256=common.sha(__file__), non_GWAS_dependency_sha256=admitted_metadata, dependency_hashes_checked_before_and_after=len(admitted_metadata), archived_GWAS_dependencies_deliberately_unread={p:h for p,h in plan['dependencies_sha256'].items() if p in forbidden}, raw_GWAS_body_reads=0, processed_GWAS_body_reads=0, genotype_decode_calls=0, real_worker_subprocesses_launched=0, estimator_calls=0, no_production_mutex_acquired=True, check_count=len(checks), checks=checks, observed_terminal_findings=observations, self_peak_RSS_bytes=peak, scope='Source/code/receipt metadata and isolated literal parser/control-flow/real-flock fixtures only; original helper imports, no numerical engine import, no real preprocessing, no source acquisition, no historical190 execution completion claim.')
    common.write_new(RESULT,record)
    print(json.dumps(dict(status=record['status'],check_count=len(checks),metadata_dependencies=len(admitted_metadata),deliberately_unread_GWAS_dependencies=len(record['archived_GWAS_dependencies_deliberately_unread']),real_workers=0,peak_RSS_bytes=peak,receipt_sha256=common.sha(RESULT)),indent=2))


if __name__=='__main__':main()
