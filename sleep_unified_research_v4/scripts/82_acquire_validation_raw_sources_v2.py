#!/usr/bin/env python3
"""Versioned, exclusively owned replay of the historical 100 source bodies."""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import plistlib
import re
import shutil
import signal
import subprocess
import sys
import time
from terminal_commit_common_v2 import TerminalCommit


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = PACKAGE / 'manifests/validation_raw_acquisition_plan_v4_2.json'
ADMISSION = PACKAGE / 'manifests/validation_raw_acquisition_admission_v4_2.json'
FOLDER = SSD / 'validation_raw_replay_v2'
MONITOR = PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py'
CURL = Path('/usr/bin/curl')

TERMINATION_REQUEST = []

def utc():
    return datetime.now(timezone.utc).isoformat()

def safe_print(value):
    try:
        print(json.dumps(value), flush=True)
    except BaseException:
        pass

def catchable_termination(signum, frame):
    # Defer catchable signals through registration, cleanup and quarantine.
    # Only the protected worker/seal path turns the request into a failure.
    TERMINATION_REQUEST.append(signum)

def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with Path(path).open('rb') as handle:
        for data in iter(lambda: handle.read(4 * 1024**2), b''):
            md5.update(data)
            sha.update(data)
    return md5.hexdigest(), sha.hexdigest()

def write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())

@contextmanager
def family_lock(path):
    # The persistent file is never removed; locks belong to the live descriptor.
    with Path(path).open('a+b') as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield handle
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

def load_monitor():
    spec = importlib.util.spec_from_file_location('raw_replay_monitor_v4', MONITOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def physical_mount():
    info = plistlib.loads(subprocess.run(
        ['diskutil', 'info', '-plist', '/Volumes/Extreme SSD'],
        capture_output=True, check=True).stdout)
    assert os.path.ismount('/Volumes/Extreme SSD')
    assert info['VolumeUUID'].upper() == '77FD98CC-B09E-3DAA-ACC0-82FF2B776B28'
    assert not info.get('ReadOnlyVolume', False)

def global_namespace_gate(plan):
    total=sum(p.stat().st_size for p in SSD.rglob('*') if p.is_file() and not p.is_symlink())
    assert total<=plan['SSD_reservation_bytes'], '300GIB_GLOBAL_NEW_CAMPAIGN_NAMESPACE_LIMIT'
    return total

def assert_bindings(plan, expected_plan_sha):
    assert hashes(PLAN)[1] == expected_plan_sha, 'FROZEN_PLAN_CHANGED'
    for path, expected in plan['bound_sources'].items():
        assert hashes(path)[1] == expected, 'BOUND_DEPENDENCY_CHANGED: ' + path
    if 'reviewed_admission_sha256' in plan:
        assert hashes(ADMISSION)[1] == plan['reviewed_admission_sha256'], 'ADMISSION_CHANGED'
        for path, expected in plan['reviewed_audit_bindings'].items():
            assert hashes(path)[1] == expected, 'INDEPENDENT_REVIEW_CHANGED'

def header_fields(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith('HTTP/'):
            result = {'status': int(line.split()[1])}
        elif ':' in line:
            key, value = line.split(':', 1)
            result[key.strip().lower()] = value.strip().strip('"')
    return result

def safe_state(monitor):
    try:
        return monitor.snapshot()
    except BaseException as exc:
        return {'snapshot_error': type(exc).__name__ + ': ' + str(exc)}

def cleanup_proof(monitor, proc):
    if proc is None:
        return {'worker_launched': False, 'remaining_group_members': [], 'teardown_verified': True}
    try:
        proof = monitor.terminate_owned(proc)
        assert not monitor.group_members(proc.pid)
        return {**proof, 'teardown_verified': True}
    except BaseException as exc:
        proof = {'cleanup_error': type(exc).__name__ + ': ' + str(exc), 'teardown_verified': False}
        # A failing shared helper cannot suppress the failure receipt. Use the
        # exact owned PGID, independently check disappearance and keep ownership
        # if even this fallback cannot certify it.
        signals = []
        try:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
                signals.append('SIGKILL')
            except ProcessLookupError:
                pass
            proc.wait(timeout=10)
            deadline = time.monotonic() + 10
            while time.monotonic() < deadline:
                try:
                    os.killpg(proc.pid, 0)
                except ProcessLookupError:
                    proof.update({'teardown_verified': True, 'remaining_group_members': [],
                                  'independent_fallback_signals': signals})
                    break
                time.sleep(.25)
        except BaseException as fallback_error:
            proof['fallback_error'] = type(fallback_error).__name__ + ': ' + str(fallback_error)
        if not proof['teardown_verified']:
            proof['group_disappearance_unverified'] = True
        return proof

def retain_lock_until_gone(proc, receipt_path):
    # This exceptional quarantine remains inside the family flock. A survivor
    # cannot outlive ownership merely because a cleanup helper raised.
    attempts = 0
    while True:
        attempts += 1
        gone = False
        try:
            proc.poll()
            os.killpg(proc.pid, 0)
        except ProcessLookupError:
            gone = True
        except BaseException as exc:
            if attempts % 30 == 1:
                safe_print({'OWNED_GROUP_QUARANTINED_LOCK_HELD': proc.pid, 'error': str(exc)})
        if gone:
            try:
                proc.wait(timeout=10)
            except BaseException:
                gone = False
            if gone:
                try:
                    write_new(Path(str(receipt_path) + '.quarantine_resolution.json'),
                              {'completed_utc': utc(), 'owned_pgid': proc.pid, 'group_disappearance_verified': True,
                               'attempts': attempts, 'primary_failure_receipt_preserved': True})
                except BaseException as exc:
                    safe_print({'GROUP_GONE_RESOLUTION_RECEIPT_ERROR': str(exc), 'owned_pgid': proc.pid})
                return
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            continue
        except BaseException:
            pass
        if attempts % 30 == 1:
            safe_print({'OWNED_GROUP_QUARANTINED_LOCK_HELD': proc.pid})
        # Interrupts do not release the lock while an owned group survives.
        try:
            time.sleep(2)
        except BaseException:
            pass

def acquire(plan, member, monitor, expected_plan_sha, family_start):
    trait = member['source_id']
    final = Path(member['body_path'])
    partial = Path(str(final) + '.partial')
    receipt_path = FOLDER / 'receipts' / (trait + '.json')
    header = FOLDER / 'logs' / (trait + '.headers.txt')
    stdout = FOLDER / 'logs' / (trait + '.curl.log')
    assert all(not p.exists() and not p.is_symlink() for p in (final, partial, receipt_path, header, stdout))
    assert all(not parent.is_symlink() for p in (final, partial, receipt_path, header, stdout) for parent in p.parents)
    before = monitor.snapshot()
    global_namespace_gate(plan)
    assert before['internal_free_bytes'] >= plan['internal_floor_bytes']
    assert before['ssd_free_bytes'] >= plan['ssd_floor_bytes']
    assert_bindings(plan, expected_plan_sha)
    resume = None
    offset = resume['prefix_bytes'] if resume else 0
    command = [str(CURL), '-q', '--fail', '--location', '--max-redirs', '5', '--proto', '=https',
               '--proto-redir', '=https', '--tlsv1.2', '--max-time', str(plan['per_body_seconds_limit']),
               '--max-filesize', str(member['expected_bytes']), '--dump-header', str(header),
               '--output', str(partial)]
    if offset:
        command += ['--continue-at', str(offset)]
    command.append(member['url'])
    receipt = {'member': member, 'started_utc': utc(), 'command': command, 'plan_sha256': expected_plan_sha,
               'resource_before': before, 'resume_offset': offset, 'status': 'RUNNING',
               'original_source_sha256_required': member['expected_sha256']}
    proc = None
    reason = None
    samples = []
    peak = 0
    start = time.monotonic()
    print(json.dumps({'source_start': member['index'], 'total': 13, 'trait': trait, 'resume_offset': offset}), flush=True)
    try:
        if resume:
            original = Path(resume['original_partial_path'])
            assert original.stat().st_size == offset and hashes(original)[1] == resume['prefix_sha256']
            with original.open('rb') as source, partial.open('xb') as target:
                shutil.copyfileobj(source, target, 4 * 1024**2)
                target.flush()
                os.fsync(target.fileno())
            assert partial.stat().st_size == offset and hashes(partial)[1] == resume['prefix_sha256']
        assert_bindings(plan,expected_plan_sha)
        physical_mount()
        global_namespace_gate(plan)
        launch_state=monitor.snapshot()
        assert launch_state['internal_free_bytes']>=plan['internal_floor_bytes'] and launch_state['ssd_free_bytes']>=plan['ssd_floor_bytes']
        assert time.monotonic()-start<=plan['per_body_seconds_limit'] and time.monotonic()-family_start<=plan['family_seconds_limit']
        assert not TERMINATION_REQUEST, 'DEFERRED_TERMINATION_BEFORE_WORKER'
        with stdout.open('x') as log:
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                    env={**os.environ, 'TMPDIR': str(SSD / 'tmp')})
            assert os.getpgid(proc.pid) == proc.pid
            while proc.poll() is None:
                state = monitor.snapshot()
                global_namespace_gate(plan)
                rss, pids = monitor.owned_rss(proc.pid)
                peak = max(peak, rss)
                size = partial.stat().st_size if partial.exists() else 0
                reason = monitor.final_limits(state, rss, size, time.monotonic() - start)
                if size > member['expected_bytes']:
                    reason = 'SOURCE_EXCEEDS_PINNED_SIZE'
                if time.monotonic() - start > plan['per_body_seconds_limit']:
                    reason = 'PER_BODY_DEADLINE'
                if time.monotonic() - family_start > plan['family_seconds_limit']:
                    reason = 'FAMILY_DEADLINE'
                samples.append({**state, 'source_bytes': size, 'owned_rss_bytes': rss, 'owned_pids': pids})
                if TERMINATION_REQUEST:
                    reason = 'DEFERRED_TERMINATION_SIGNAL: ' + str(TERMINATION_REQUEST)
                if reason:
                    raise RuntimeError(reason)
                if len(samples) % 30 == 0:
                    assert hashes(PLAN)[1] == expected_plan_sha
                    assert hashes(Path(__file__))[1] == plan['bound_sources'][str(Path(__file__))]
                    assert hashes(MONITOR)[1] == plan['bound_sources'][str(MONITOR)]
                    print(json.dumps({'source_index': member['index'], 'downloaded_bytes': size,
                                      'elapsed_seconds': round(time.monotonic()-start, 1)}), flush=True)
                time.sleep(plan['runtime_poll_seconds'])
        assert not TERMINATION_REQUEST, 'DEFERRED_TERMINATION_AFTER_WORKER'
        assert proc.returncode == 0, 'CURL_NONZERO_EXIT: ' + str(proc.returncode)
        assert not monitor.group_members(proc.pid), 'UNEXPECTED_TRANSFER_DESCENDANTS'
        observed = header_fields(header)
        assert observed['status'] == (206 if offset else 200)
        assert_headers(member, observed)
        assert int(observed['content-length']) == member['expected_bytes'] - offset
        if offset:
            assert observed['content-range'] == f"bytes {offset}-{member['expected_bytes']-1}/{member['expected_bytes']}"
        assert partial.stat().st_size == member['expected_bytes']
        md5, sha256 = hashes(partial)
        assert (md5, sha256) == (member['expected_md5'], member['expected_sha256'])
        assert_bindings(plan, expected_plan_sha)
        physical_mount()
        assert monitor.final_limits(monitor.snapshot(), peak, partial.stat().st_size, time.monotonic()-start) is None
        assert time.monotonic()-start <= plan['per_body_seconds_limit']
        assert time.monotonic()-family_start <= plan['family_seconds_limit']
        receipt.update({'status': 'EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED', 'actual_size': partial.stat().st_size,
                        'actual_md5': md5, 'actual_sha256': sha256, 'observed_headers': observed,
                        'gzip_full_stream_verified': False, 'source_body_committed': False})
    except BaseException as exc:
        reason = reason or type(exc).__name__ + ': ' + str(exc)
        receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
    finally:
        teardown = {'teardown_verified': False, 'cleanup_not_yet_certified': True}
        try:
            teardown = cleanup_proof(monitor, proc)
            if not teardown['teardown_verified'] or teardown.get('cleanup_error'):
                reason = 'PROCESS_TEARDOWN_FAILURE: ' + teardown.get('cleanup_error', '')
                receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
            # Cleanup failures cannot bypass the immutable failure receipt.
            receipt.update({'completed_utc': utc(), 'elapsed_seconds': time.monotonic()-start,
                            'returncode': proc.returncode if proc else None, 'stop_reason': reason,
                            'teardown': teardown, 'peak_observed_owned_rss_bytes': peak,
                            'resource_after': safe_state(monitor), 'samples': samples})
            try:
                for label, path in [('headers', header), ('retained_partial', partial)]:
                    try:
                        md5_now, sha_now = hashes(path) if path.exists() else (None, None)
                        receipt[label + '_sha256'] = sha_now
                        receipt[label + '_bytes'] = path.stat().st_size if path.exists() else None
                        if label == 'retained_partial' and receipt['status'] == 'EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED':
                            assert (md5_now, sha_now, path.stat().st_size) == (
                                member['expected_md5'], member['expected_sha256'], member['expected_bytes'])
                            receipt['final_seal_md5'] = md5_now
                            receipt['final_seal_sha256'] = sha_now
                    except BaseException as exc:
                        receipt[label + '_inspection_error'] = str(exc)
                        receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
                        receipt['stop_reason'] = 'FINAL_BODY_INSPECTION_FAILURE: ' + str(exc)
                if receipt['status'] == 'EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED':
                    try:
                        if resume:
                            original = Path(resume['original_partial_path'])
                            assert original.stat().st_size == offset and hashes(original)[1] == resume['prefix_sha256']
                        assert_bindings(plan, expected_plan_sha)
                        physical_mount()
                        final_state = monitor.snapshot()
                        global_namespace_gate(plan)
                        final_reason = monitor.final_limits(final_state, peak, partial.stat().st_size, time.monotonic()-start)
                        assert final_reason is None, final_reason
                        assert time.monotonic()-start <= plan['per_body_seconds_limit'], 'FINAL_PER_BODY_DEADLINE'
                        assert time.monotonic()-family_start <= plan['family_seconds_limit'], 'FINAL_FAMILY_DEADLINE'
                        assert not TERMINATION_REQUEST, 'DEFERRED_TERMINATION_BEFORE_SEAL'
                        receipt['resource_after'] = final_state
                        receipt['post_cleanup_hash_resource_identity_gates_pass'] = True
                        assert not final.exists() and not final.is_symlink()
                        partial.rename(final)
                    except BaseException as exc:
                        receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
                        receipt['stop_reason'] = 'FINAL_SEAL_FAILURE: ' + str(exc)
                write_new(receipt_path, receipt)
                write_new(PACKAGE / 'source_provenance/validation_raw_acquisition_v4_2' / receipt_path.name, receipt)
            finally:
                # Persistence, stdout and catchable termination failures cannot
                # bypass quarantine while group disappearance is unverified.
                if not teardown['teardown_verified'] and proc is not None:
                    retain_lock_until_gone(proc, receipt_path)
        finally:
            if not teardown['teardown_verified'] and proc is not None:
                retain_lock_until_gone(proc, receipt_path)
    if receipt['status'] != 'EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED':
        raise RuntimeError('SOURCE_STOP_PRESERVED: ' + str(receipt.get('stop_reason')))
    print(json.dumps({'source_complete': member['index'], 'sha256': receipt['actual_sha256'],
                      'elapsed_seconds': receipt['elapsed_seconds']}), flush=True)
    return {'path': str(receipt_path), 'sha256': hashes(receipt_path)[1]}


def assert_headers(member, observed):
    if observed['status'] != 200 or int(observed['content-length']) != member['expected_bytes']:
        raise RuntimeError('EXACT_FULL_BODY_HTTP_STATUS_AND_LENGTH_REQUIRED')
    for key, expected in member['expected_transport_headers'].items():
        if observed.get(key) != expected:
            raise RuntimeError('FROZEN_TRANSPORT_HEADER_CHANGED: ' + key)


def prepare():
    design_path = PACKAGE/'statistical_validation/validation_raw_replay_design_v1.json'
    if hashes(design_path)[1] != 'f5291abd8989e6444f9aabe9ce2b0dd15258412acaa850dac1c59e72d86c7403':
        raise RuntimeError('EXACT_ORIGINAL13_DESIGN_CHANGED')
    design = json.loads(design_path.read_text())
    seal_path = PACKAGE/'statistical_validation/validation_raw_replay_design_seal_v1.json'
    seal = json.loads(seal_path.read_text())
    bound = dict(design['metadata_dependencies_sha256'])
    bound.update(seal['file_sha256'])
    files = [design_path,seal_path,Path(__file__),MONITOR,CURL,
             PACKAGE/'scripts/52_acquire_extension_raw_sources_v8.py',
             PACKAGE/'scripts/81_acquire_validation_raw_sources_v1.py',
             PACKAGE/'manifests/validation_raw_acquisition_plan_v4_1.json',
             PACKAGE/'scripts/terminal_commit_common_v2.py',
             PACKAGE/'manifests/global_SSD_resource_reservation_v4_4.json']
    for path in files:bound[str(path)] = hashes(path)[1]
    if bound[str(files[-1])] != '04e2ea61d9a06281468674c1a8326367ffb7e507857814365e63108ee38ebbc7':
        raise RuntimeError('CONSISTENT_300GIB_LEDGER_REQUIRED')
    prior_seal_path=PACKAGE/'reviews/validation13_acquisition_independent_review_seal_v1.json'
    if hashes(prior_seal_path)[1]!='d40c0524acfe924c3ae86cbdfec1fe8d7df30afee3c5aa32cb377e6361223de9':raise RuntimeError('PRESERVED_V1_REJECTION_CHANGED')
    bound[str(prior_seal_path)]=hashes(prior_seal_path)[1]
    prior_seal=json.loads(prior_seal_path.read_text())
    for path,artifact in prior_seal.get('file_sha256',prior_seal.get('artifacts',{})).items():
        digest=artifact['sha256'] if isinstance(artifact,dict) else artifact
        if hashes(path)[1]!=digest:raise RuntimeError('PRESERVED_V1_REVIEW_ARTIFACT_CHANGED')
        bound[path]=digest
    ledger = json.loads(files[-1].read_text())
    if ledger['component_bytes']['proposed_validation13_raw_bytes'] != 10058648185 or ledger['reserved_total_bytes'] >= ledger['ceiling_bytes']:
        raise RuntimeError('EXACT_ORIGINAL13_RESERVATION_REQUIRED')
    members = []
    for m in design['members']:
        historical = m['historical_source_body']
        probe = m['HEAD_probe']
        fields = probe['headers']
        if probe['status'] != 200 or int(fields['content-length']) != historical['observed_size_bytes']:
            raise RuntimeError('FROZEN_PUBLIC_SOURCE_LENGTH_CHANGED')
        transport = {'etag':fields['etag'].strip('"'),'last-modified':fields['last-modified']}
        if m['source_id'].startswith('finngen_'):
            transport.update({key:fields[key] for key in ['x-goog-generation','x-goog-hash','x-goog-stored-content-length']})
            if '?generation='+fields['x-goog-generation'] not in historical['source'] or probe['identity_checks'] != {'etag':True,'generation':True,'length':True}:
                raise RuntimeError('EXACT_GCS_GENERATION_REQUIRED')
        elif m['current_official_yaml_MD5_if_MVP'] != historical['observed_md5']:
            raise RuntimeError('MVP_CURRENT_OFFICIAL_AND_ORIGINAL_BODY_MD5_DIFFER')
        filename = Path(m['proposed_body_path']).name
        if Path(m['proposed_body_path']) != SSD/'validation_raw_replay_v1/raw'/filename:
            raise RuntimeError('EXACT_FROZEN_PRIVATE_SOURCE_ROUTE_REQUIRED')
        members.append(dict(index=m['index'],source_id=m['source_id'],filename=filename,
                            body_path=str(FOLDER/'raw'/filename),url=historical['source'],
                            expected_bytes=historical['observed_size_bytes'],expected_md5=historical['observed_md5'],
                            expected_sha256=historical['observed_sha256'],expected_transport_headers=transport,
                            original_http_etag=historical['http_etag'],
                            transport_qualification='Pinned GCS generation and original full SHA/MD5' if m['source_id'].startswith('finngen_') else 'Current EBI HEAD ETag differs from historical transport metadata; original complete SHA/MD5/size and current official YAML MD5 must match. No changed body is accepted.'))
    if len(members) != 13 or sum(m['expected_bytes'] for m in members) != 10058648185:
        raise RuntimeError('EXACT_ORIGINAL13_CARDINALITY_AND_BYTES_REQUIRED')
    for path,digest in bound.items():
        if Path(path).is_symlink() or not Path(path).is_file() or hashes(path)[1] != digest:
            raise RuntimeError('FROZEN_DESIGN_DEPENDENCY_CHANGED: '+path)
    physical_mount()
    state = load_monitor().snapshot()
    if state['internal_free_bytes'] < 3<<30 or state['ssd_free_bytes'] < 300<<30:
        raise RuntimeError('UNCHANGED_STORAGE_RESERVATION_GUARD_FAILS')
    for folder in ['raw','logs','receipts']:(FOLDER/folder).mkdir(parents=True,exist_ok=True)
    plan = dict(schema='frozen_original13_validation_source_acquisition_v2',version=2,prepared_utc=utc(),
                member_count=13,members=members,compressed_network_bytes=10058648185,
                private_namespace=str(FOLDER),bound_sources=bound,executor_path=str(Path(__file__)),
                internal_floor_bytes=3<<30,ssd_floor_bytes=5<<30,SSD_reservation_bytes=300<<30,
                maximum_owned_transfer_rss_bytes=2<<30,monitor_output_limit_bytes=4<<30,
                per_body_seconds_limit=7200,family_seconds_limit=96*3600,runtime_poll_seconds=2,
                exclusive_family_lock_path=str(SSD/'extension_raw_acquisition_family.lock'),
                pending_path=str(FOLDER/'family_pending_v4_2.json'),terminal_seal_path=str(FOLDER/'family_terminal_seal_v4_2.json'),
                terminal_protocol='REVIEWED_TERMINAL2_WITH_TWO_EXACT_PRIMARY_COPIES;FULL_CURRENT13_BODY_HASHES_AND_ALL_RECEIPT_COPIES_CHECKED_BEFORE_AND_AFTER_PERSISTENCE',
                curl_version=subprocess.run([str(CURL),'-q','--version'],capture_output=True,text=True,check=True).stdout,
                global_resource_ledger_path=str(files[-1]),global_resource_ledger_sha256=bound[str(files[-1])],
                resource_at_plan=state,source_bodies_read_in_preparation=0,workers_launched=0,
                automatic_retry=False,source_substitution=False,scientific_membership_changes=False,
                independent_two_trait_replication_claimed=False,shared_single_transfer_family=True,
                human_or_restricted_access_requested=False)
    write_new(PLAN,plan)
    write_new(SSD/'manifests'/PLAN.name,plan)
    safe_print(dict(plan_sha256=hashes(PLAN)[1],sources=13,source_body_bytes=10058648185,workers_launched=0))


def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(parent.is_symlink() for parent in p.parents):
        raise RuntimeError('EXACT_REGULAR_ACQUISITION_EVIDENCE_REQUIRED: '+str(p))
    return p


def execute(expected_plan_sha):
    if not re.fullmatch(r'[0-9a-f]{64}',expected_plan_sha) or hashes(regular(PLAN))[1] != expected_plan_sha:
        raise RuntimeError('EXACT_FROZEN_ORIGINAL13_PLAN_REQUIRED')
    plan=json.loads(PLAN.read_text());admission=json.loads(regular(ADMISSION).read_text())
    if admission.get('execution_admitted') is not True or admission.get('plan_sha256')!=expected_plan_sha or admission.get('executor_sha256')!=hashes(Path(__file__))[1]:
        raise RuntimeError('EXACT_ROOT_ORIGINAL13_ADMISSION_REQUIRED')
    reviews=admission['independent_review_artifact_sha256']
    if not reviews or not any(p.endswith('.md') for p in reviews) or not any(p.endswith('.json') for p in reviews):
        raise RuntimeError('INDEPENDENT_OPERATIONAL_REVIEW_REQUIRED')
    plan['reviewed_admission_sha256']=hashes(ADMISSION)[1];plan['reviewed_audit_bindings']=reviews
    if plan['member_count']!=13 or len(plan['members'])!=13 or len({m['source_id'] for m in plan['members']})!=13 or len({m['filename'] for m in plan['members']})!=13:
        raise RuntimeError('EXACT_ORIGINAL13_MEMBERSHIP_REQUIRED')
    if [m['index'] for m in plan['members']] != list(range(1,14)) or sum(m['expected_bytes'] for m in plan['members'])!=10058648185 or plan['compressed_network_bytes']!=10058648185:
        raise RuntimeError('EXACT_ORIGINAL13_ORDER_AND_BYTES_REQUIRED')
    assert_bindings(plan,expected_plan_sha)
    monitor=load_monitor()
    if (monitor.INTERNAL_FLOOR,monitor.SSD_FLOOR,monitor.RSS_LIMIT,monitor.OUTPUT_LIMIT)!=(plan['internal_floor_bytes'],plan['ssd_floor_bytes'],plan['maximum_owned_transfer_rss_bytes'],plan['monitor_output_limit_bytes']):
        raise RuntimeError('UNCHANGED_TRANSFER_MONITOR_LIMITS_REQUIRED')
    if subprocess.run([str(CURL),'-q','--version'],capture_output=True,text=True,check=True).stdout!=plan['curl_version']:
        raise RuntimeError('FROZEN_CURL_RUNTIME_CHANGED')
    physical_mount()
    with family_lock(plan['exclusive_family_lock_path']):
        if (FOLDER/'family_ownership.json').exists() or (FOLDER/'family_ownership.json').is_symlink():
            raise RuntimeError('PRIOR_ACQUISITION_ATTEMPT_PRESERVED_NO_RETRY')
        start=time.monotonic();master=dict(status='RUNNING',started_utc=utc(),plan_sha256=expected_plan_sha,source_receipts=[],owner_pid=os.getpid())
        binding=dict(plan_sha256=expected_plan_sha,admission_sha256=hashes(ADMISSION)[1],executor_sha256=hashes(Path(__file__))[1])
        fixed_primary_hashes={}
        terminal=TerminalCommit(plan['pending_path'],plan['terminal_seal_path'],binding)
        write_new(FOLDER/'family_ownership.json',dict(owner_pid=os.getpid(),started_utc=utc(),plan_sha256=expected_plan_sha,immutable_attempt=2))
        def resources():
            state=monitor.snapshot();global_namespace_gate(plan)
            if state['internal_free_bytes']<plan['internal_floor_bytes'] or state['ssd_free_bytes']<plan['ssd_floor_bytes'] or time.monotonic()-start>plan['family_seconds_limit']:
                raise RuntimeError('FINAL_SOURCE_FAMILY_RESOURCE_OR_DEADLINE_GUARD')
            return state
        def termination():
            if TERMINATION_REQUEST:raise RuntimeError('DEFERRED_TERMINATION_SIGNAL')
        def identity():
            assert_bindings(plan,expected_plan_sha);physical_mount()
            regular(plan['pending_path'])
            if Path(plan['terminal_seal_path']).exists():regular(plan['terminal_seal_path'])
            for path,digest in fixed_primary_hashes.items():
                if hashes(regular(path))[1]!=digest:raise RuntimeError('FROZEN_INTENDED_FAMILY_PRIMARY_CHANGED')
            if len(master['source_receipts'])!=13:raise RuntimeError('ALL13_EXACT_BODY_PROOFS_REQUIRED')
            for member,item in zip(plan['members'],master['source_receipts']):
                path=regular(item['path']);mirror=regular(PACKAGE/'source_provenance/validation_raw_acquisition_v4_2'/path.name)
                if hashes(path)[1]!=item['sha256'] or hashes(mirror)[1]!=item['sha256']:raise RuntimeError('EXACT_BOTH_SOURCE_RECEIPT_COPIES_REQUIRED')
                receipt=json.loads(path.read_text());body=regular(member['body_path'])
                if receipt['member']!=member or receipt['status']!='EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED' or receipt['plan_sha256']!=expected_plan_sha or receipt['returncode']!=0 or receipt['stop_reason'] is not None or not receipt['teardown']['teardown_verified'] or receipt['teardown'].get('remaining_group_members') or receipt['teardown'].get('cleanup_error'):
                    raise RuntimeError('EXACT_COMPLETED_REAPED_SOURCE_PRODUCER_REQUIRED')
                if (body.stat().st_size,*hashes(body))!=(member['expected_bytes'],member['expected_md5'],member['expected_sha256']):
                    raise RuntimeError('CURRENT_FULL_ORIGINAL_SOURCE_IDENTITY_CHANGED')
                assert_headers(member,receipt['observed_headers'])
                for value in [path,mirror]:
                    if Path(str(value)+'.failure.json').exists() or Path(str(value)+'.failure.json').is_symlink():raise RuntimeError('SOURCE_FAILURE_VETOES_FAMILY_SUCCESS')
            termination();resources()
        try:
            resources();termination()
            for member in plan['members']:
                master['source_receipts'].append(acquire(plan,member,monitor,expected_plan_sha,start))
            identity()
            master.update(status='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED',full_current_body_SHA_MD5_size_verified=True,
                          scientific_replication_admitted=False,gzip_full_EOF_or_pipeline_replay_certified=False)
        except BaseException as error:
            master.update(status='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT',stop_reason=type(error).__name__+': '+str(error))
        master.update(completed_utc=utc(),elapsed_seconds=time.monotonic()-start,completed_source_count=len(master['source_receipts']),resource_after=safe_state(monitor))
        primaries=[PACKAGE/'logs/validation_raw_acquisition_family_receipt_v4_2.json',FOLDER/'validation_raw_acquisition_family_receipt_v4_2.json']
        try:
            intended_primary_hash=hashlib.sha256((json.dumps(master,indent=2)+'\n').encode()).hexdigest()
            fixed_primary_hashes={str(path):intended_primary_hash for path in primaries}
            for path in primaries:write_new(path,master)
            primary_hashes={str(path):hashes(regular(path))[1] for path in primaries}
            if primary_hashes!=fixed_primary_hashes:raise RuntimeError('TWO_PRIMARY_FAMILY_RECEIPTS_DIFFER_FROM_FROZEN_INTENT')
            if master['status']=='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED':
                if not terminal.commit(primary_hashes,identity_gate=identity,resource_gate=resources,termination_gate=termination):
                    raise RuntimeError('SOURCE_FAMILY_TERMINAL2_COMMIT_FAILED')
        except BaseException as error:
            for path in primaries:
                try:write_new(Path(str(path)+'.failure.json'),dict(status='FAMILY_FINALIZATION_FAILED_PRESERVED',error=repr(error)))
                except BaseException:pass
            raise
        if master['status']!='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED':
            raise SystemExit('ACQUISITION_STOP_PRESERVED_REQUIRES_REVIEW')


if __name__=='__main__':
    for name in ['SIGINT','SIGTERM','SIGHUP']:signal.signal(getattr(signal,name),catchable_termination)
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--execute',action='store_true')
    p.add_argument('--expected-plan-sha256');args=p.parse_args()
    if args.prepare:prepare()
    else:
        if not args.expected_plan_sha256:p.error('Exact independently reviewed plan SHA required')
        execute(args.expected_plan_sha256)
