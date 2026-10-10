#!/usr/bin/env python3
"""Acquire the exact 13 historical validation bodies with frozen terminal intent."""
import argparse
import copy
import errno
import stat
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
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
PLAN = PACKAGE / 'manifests/validation_raw_acquisition_plan_v4_4.json'
ADMISSION = PACKAGE / 'manifests/validation_raw_acquisition_admission_v4_4.json'
FOLDER = SSD / 'validation_raw_replay_v4'
ORIGINAL_PLAN = PACKAGE / 'manifests/validation_raw_acquisition_plan_v4_3.json'
ORIGINAL_PLAN_SHA = '8e93d548163bda284b6292b9780962060bbc3ff4afde44e5a916c80f3b69b708'
STOP = PACKAGE / 'logs/validation_raw13_v3_actual_stop_root_receipt_v1.json'
STOP_SHA = 'e1142ad40f71b73c0af1e0cd2526902628b614d74df04561b2ed46a1261fd606'
LEDGER = PACKAGE / 'manifests/global_SSD_resource_reservation_v4_6.json'
PREFIX_BYTES = 57671680
PREFIX_SHA = '3537a1b14c7e3560c43d38898dfadd6ba280738498ef5de4c21195c962a1a415'
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
    """Non-atomic live regular-byte meter; only disappearing entries may vanish.

    Each entry has one no-follow stat. Directory descriptors prevent following
    replacement symlinks. Root absence, permissions, EIO and all other errors
    fail the stage. Bound scientific/evidence files are separately fail-closed.
    """
    total = 0
    def disappeared(path):
        safe_print({'GLOBAL_METER_CONCURRENT_ENOENT_REMOVAL': str(path)})
    def walk(fd, path):
        nonlocal total
        try:
            with os.scandir(fd) as entries:
                for entry in entries:
                    child = path / entry.name
                    try:
                        info = entry.stat(follow_symlinks=False)
                    except OSError as exc:
                        if exc.errno != errno.ENOENT:
                            raise
                        disappeared(child)
                        continue
                    if stat.S_ISREG(info.st_mode):
                        total += info.st_size
                        if total > plan['SSD_reservation_bytes']:
                            raise RuntimeError('300GIB_GLOBAL_NEW_CAMPAIGN_NAMESPACE_LIMIT')
                    elif stat.S_ISDIR(info.st_mode):
                        try:
                            child_fd = os.open(entry.name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
                        except OSError as exc:
                            if exc.errno != errno.ENOENT:
                                raise
                            disappeared(child)
                            continue
                        try:
                            walk(child_fd, child)
                        finally:
                            os.close(child_fd)
        except OSError as exc:
            if exc.errno != errno.ENOENT:
                raise
            disappeared(path)
    root_fd = os.open(SSD, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        walk(root_fd, SSD)
    finally:
        os.close(root_fd)
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

def transfer_command(plan, member, header, partial, offset=0):
    command = [str(CURL), '-q', '--fail', '--location', '--max-redirs', '5', '--proto', '=https',
               '--proto-redir', '=https', '--tlsv1.2', '--max-time', str(plan['per_body_seconds_limit']),
               '--max-filesize', str(member['expected_bytes']), '--dump-header', str(header),
               '--output', str(partial)]
    if offset:
        command += ['--continue-at', str(offset), '--header', 'If-Match: ' + member['original_http_etag']]
    command.append(member['url'])
    return command


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
    resume = plan['resume_source12'] if member['index'] == 12 else None
    offset = resume['prefix_bytes'] if resume else 0
    command = transfer_command(plan,member,header,partial,offset)
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
            original = regular(resume['original_partial_path'])
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
        assert_headers(member, observed, offset)
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
                            original = regular(resume['original_partial_path'])
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
                write_new(PACKAGE / 'source_provenance/validation_raw_acquisition_v4_4' / receipt_path.name, receipt)
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


def assert_headers(member, observed, offset=0):
    if observed['status'] != (206 if offset else 200) or int(observed['content-length']) != member['expected_bytes'] - offset:
        raise RuntimeError('EXACT_BODY_HTTP_STATUS_AND_LENGTH_REQUIRED')
    if offset and observed.get('content-range') != f"bytes {offset}-{member['expected_bytes']-1}/{member['expected_bytes']}":
        raise RuntimeError('EXACT_PINNED_RANGE_REQUIRED')
    for key, expected in member['expected_transport_headers'].items():
        if observed.get(key) != expected:
            raise RuntimeError('FROZEN_TRANSPORT_HEADER_CHANGED: ' + key)


def source_receipt_gate(plan, member, item, current_plan_sha, verify_body=True):
    origin = plan['receipt_origins'][member['source_id']]
    path, mirror = regular(item['path']), regular(origin['mirror_path'])
    if str(path) != origin['primary_path']:
        raise RuntimeError('EXACT_MIXED_ORIGIN_SOURCE_RECEIPT_REQUIRED')
    if origin['reused_v3'] and item['sha256'] != origin['sha256']:
        raise RuntimeError('ORIGINAL_V3_RECEIPT_IDENTITY_CHANGED')
    if hashes(path)[1] != item['sha256'] or hashes(mirror)[1] != item['sha256']:
        raise RuntimeError('EXACT_BOTH_SOURCE_RECEIPT_COPIES_REQUIRED')
    receipt = json.loads(path.read_text())
    expected_origin = ORIGINAL_PLAN_SHA if origin['reused_v3'] else current_plan_sha
    offset = 0 if origin['reused_v3'] or member['index'] != 12 else plan['resume_source12']['prefix_bytes']
    if (receipt['member'] != member or receipt['status'] != 'EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED'
            or receipt['plan_sha256'] != expected_origin or receipt['returncode'] != 0
            or receipt['stop_reason'] is not None or receipt['teardown']['teardown_verified'] is not True
            or receipt['teardown'].get('remaining_group_members') != [] or receipt['teardown'].get('cleanup_error')
            or receipt['resume_offset'] != offset or receipt.get('post_cleanup_hash_resource_identity_gates_pass') is not True):
        raise RuntimeError('EXACT_COMPLETED_REAPED_SOURCE_PRODUCER_REQUIRED')
    if (receipt['actual_size'], receipt['actual_md5'], receipt['actual_sha256'],
            receipt['final_seal_md5'], receipt['final_seal_sha256']) != (
            member['expected_bytes'], member['expected_md5'], member['expected_sha256'],
            member['expected_md5'], member['expected_sha256']):
        raise RuntimeError('EXACT_ORIGINAL_RECEIPT_BODY_IDENTITY_REQUIRED')
    if receipt['command'] != transfer_command(plan,member,origin['header_path'],member['body_path']+'.partial',offset):
        raise RuntimeError('EXACT_ORIGIN_TRANSFER_COMMAND_REQUIRED')
    header = regular(origin['header_path'])
    if hashes(header)[1] != receipt['headers_sha256'] or header_fields(header) != receipt['observed_headers']:
        raise RuntimeError('EXACT_ORIGIN_HEADER_EVIDENCE_REQUIRED')
    assert_headers(member, receipt['observed_headers'], offset)
    for value in [path, mirror]:
        failure = Path(str(value) + '.failure.json')
        if failure.exists() or failure.is_symlink():
            raise RuntimeError('SOURCE_FAILURE_VETOES_FAMILY_SUCCESS')
    if hashes(path)[1] != item['sha256'] or hashes(mirror)[1] != item['sha256']:
        raise RuntimeError('SOURCE_RECEIPT_CHANGED_DURING_CONSUMPTION')
    body = regular(member['body_path'])
    if body.stat().st_size != member['expected_bytes']:
        raise RuntimeError('CURRENT_ORIGINAL_SOURCE_SIZE_CHANGED')
    if verify_body and hashes(body) != (member['expected_md5'], member['expected_sha256']):
        raise RuntimeError('CURRENT_FULL_ORIGINAL_SOURCE_IDENTITY_CHANGED')
    return receipt


def prepare(diagnosis_seal, diagnosis_seal_sha, diagnosis_root, diagnosis_root_sha):
    """Freeze only this continuation metadata; no source-body reads/SSD writes."""
    if hashes(regular(ORIGINAL_PLAN))[1] != ORIGINAL_PLAN_SHA or hashes(regular(STOP))[1] != STOP_SHA:
        raise RuntimeError('EXACT_STOPPED_V3_CHECKPOINT_REQUIRED')
    old = json.loads(ORIGINAL_PLAN.read_text())
    stop = json.loads(STOP.read_text())
    if (stop['source_bodies_completed'], stop['failed_source_index'], stop['failed_partial_bytes'],
            stop['failed_partial_sha256'], stop['whole_family_terminal_success'], stop['pending_retained'],
            stop['source13_started'], stop['owned_teardown_verified']) != (11,12,PREFIX_BYTES,PREFIX_SHA,False,True,False,True):
        raise RuntimeError('EXACT_PRESERVED11_PREFIX12_CHECKPOINT_REQUIRED')
    for path, digest in [(diagnosis_seal,diagnosis_seal_sha),(diagnosis_root,diagnosis_root_sha)]:
        if not re.fullmatch('[0-9a-f]{64}', digest or '') or hashes(regular(path))[1] != digest:
            raise RuntimeError('FORMAL_STAT_DIAGNOSIS_AND_ROOT_CONSUMPTION_REQUIRED')
    ledger = json.loads(regular(LEDGER).read_text())
    if ledger['ceiling_bytes'] != 300<<30 or ledger['component_bytes']['proposed_validation13_raw_bytes'] != 10058648185 or ledger['preserved_validation12_ENOENT_partial_bytes'] != PREFIX_BYTES or ledger['reserved_total_bytes'] >= ledger['ceiling_bytes']:
        raise RuntimeError('EXACT_ADDITIVE_PREFIX_RESERVATION_REQUIRED')
    original_family = json.loads(regular(SSD/'validation_raw_replay_v3/validation_raw_acquisition_family_receipt_v4_3.json').read_text())
    if original_family['status'] != 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT' or original_family['completed_source_count'] != 11 or len(original_family['source_receipts']) != 11:
        raise RuntimeError('ORIGINAL_FAILED_FAMILY_NOT_COMPLETE_REQUIRED')
    original_start = datetime.fromisoformat(original_family['started_utc'])
    original_deadline = original_start + timedelta(seconds=old['family_seconds_limit'])
    plan = copy.deepcopy(old)
    for member in plan['members'][11:]:
        member['body_path'] = str(FOLDER/'raw'/member['filename'])
    bound = dict(old['bound_sources'])
    # Bind failed-family metadata and source receipts, never the source bodies.
    for path,digest in stop['file_sha256'].items():
        if path.endswith('.partial'):
            continue
        if hashes(regular(path))[1] != digest:
            raise RuntimeError('PRESERVED_STOP_METADATA_CHANGED')
        bound[path] = digest
    files = [ORIGINAL_PLAN,STOP,Path(__file__),LEDGER,Path(diagnosis_seal),Path(diagnosis_root),
             PACKAGE/'reviews/validation13_acquisition_independent_review_seal_v3.json']
    for path in files:
        bound[str(path)] = hashes(regular(path))[1]
    origins = {}
    for member in plan['members']:
        sid = member['source_id']; reused = member['index'] <= 11
        folder = SSD/'validation_raw_replay_v3' if reused else FOLDER
        mirror = PACKAGE/('source_provenance/validation_raw_acquisition_v4_3' if reused else 'source_provenance/validation_raw_acquisition_v4_4')/(sid+'.json')
        primary = folder/'receipts'/(sid+'.json')
        origins[sid] = dict(reused_v3=reused,primary_path=str(primary),mirror_path=str(mirror),header_path=str(folder/'logs'/(sid+'.headers.txt')))
        if reused:
            origins[sid]['sha256'] = stop['file_sha256'][str(primary)]
    prefix = old['members'][11]['body_path']+'.partial'
    plan.update(schema='frozen_original13_validation_source_acquisition_v4',version=4,prepared_utc=utc(),
        members=plan['members'],private_namespace=str(FOLDER),executor_path=str(Path(__file__)),bound_sources=bound,
        compressed_network_bytes=sum(m['expected_bytes'] for m in plan['members'][11:])-PREFIX_BYTES,
        original_full_body_bytes=10058648185,reused_source_count=11,new_transfer_count=2,receipt_origins=origins,
        resume_source12=dict(original_partial_path=prefix,prefix_bytes=PREFIX_BYTES,prefix_sha256=PREFIX_SHA,
            original_failed_receipt_path=str(SSD/'validation_raw_replay_v3/receipts/finngen_r13_K11_CHOLELITH.json')),
        pending_path=str(FOLDER/'family_pending_v4_4.json'),terminal_seal_path=str(FOLDER/'family_terminal_seal_v4_4.json'),
        global_resource_ledger_path=str(LEDGER),global_resource_ledger_sha256=bound[str(LEDGER)],
        original_v3_plan_path=str(ORIGINAL_PLAN),original_v3_plan_sha256=ORIGINAL_PLAN_SHA,
        preserved_v3_stop_root_receipt_path=str(STOP),preserved_v3_stop_root_receipt_sha256=STOP_SHA,
        formal_diagnosis_seal_path=str(diagnosis_seal),formal_diagnosis_root_path=str(diagnosis_root),
        source_bodies_read_in_preparation=0,workers_launched=0,
        preparation_qualification='CURRENT_METADATA_ONLY;NO_FRESH11_BODY_HASHES_RUNTIME_CENSUS_MOUNT_OR_LOCK;ROOT_ADMISSION_REQUIRED',
        resource_meter_qualification='NON_ATOMIC_NOFOLLOW_SINGLE_STAT_CENSUS;ONLY_CONCURRENT_ENOENT_ENTRY_REMOVAL_SKIPPED_WITH_DIAGNOSTIC;BOUND_EVIDENCE_REMAINS_FAIL_CLOSED',
        old_failed_family_never_relabelled=True,automatic_retry=False,
        original_family_started_utc=original_start.isoformat(),original_family_deadline_utc=original_deadline.isoformat())
    for member in plan['members'][:11]:
        origin=origins[member['source_id']]
        source_receipt_gate(plan,member,dict(path=origin['primary_path'],sha256=origin['sha256']),None,verify_body=False)
    write_new(PLAN,plan)
    safe_print(dict(plan_sha256=hashes(PLAN)[1],reused_sources=11,new_transfers=2,source12_prefix_bytes=PREFIX_BYTES,workers_launched=0))



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
    if [m['index'] for m in plan['members']] != list(range(1,14)) or sum(m['expected_bytes'] for m in plan['members'])!=10058648185 or plan['original_full_body_bytes']!=10058648185 or plan['compressed_network_bytes']!=1557744023 or plan['reused_source_count']!=11 or plan['new_transfer_count']!=2:
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
        for folder in ['raw','logs','receipts']:(FOLDER/folder).mkdir(parents=True,exist_ok=True)
        continuation_start=time.monotonic()
        remaining=(datetime.fromisoformat(plan['original_family_deadline_utc'])-datetime.now(timezone.utc)).total_seconds()
        if remaining <= 0:raise RuntimeError('ORIGINAL96H_FAMILY_DEADLINE_EXPIRED')
        start=continuation_start-(plan['family_seconds_limit']-remaining)
        master=dict(status='RUNNING',started_utc=utc(),plan_sha256=expected_plan_sha,source_receipts=[],owner_pid=os.getpid(),
                    original_family_started_utc=plan['original_family_started_utc'],original_family_deadline_utc=plan['original_family_deadline_utc'])
        binding=dict(plan_sha256=expected_plan_sha,admission_sha256=hashes(ADMISSION)[1],executor_sha256=hashes(Path(__file__))[1])
        fixed_primary_hashes={}
        fixed_terminal_intent_sha=None
        terminal=TerminalCommit(plan['pending_path'],plan['terminal_seal_path'],binding)
        write_new(FOLDER/'family_ownership.json',dict(owner_pid=os.getpid(),started_utc=utc(),plan_sha256=expected_plan_sha,immutable_attempt=4))
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
            if Path(plan['terminal_seal_path']).exists() or Path(plan['terminal_seal_path']).is_symlink():
                if fixed_terminal_intent_sha is None or hashes(regular(plan['terminal_seal_path']))[1]!=fixed_terminal_intent_sha:
                    raise RuntimeError('FROZEN_PRIVATE_TERMINAL_SEAL_INTENT_CHANGED')
            for path,digest in fixed_primary_hashes.items():
                if hashes(regular(path))[1]!=digest:raise RuntimeError('FROZEN_INTENDED_FAMILY_PRIMARY_CHANGED')
            if len(master['source_receipts'])!=13:raise RuntimeError('ALL13_EXACT_BODY_PROOFS_REQUIRED')
            for member,item in zip(plan['members'],master['source_receipts']):
                source_receipt_gate(plan,member,item,expected_plan_sha)
            termination();resources()
        try:
            resources();termination()
            for member in plan['members']:
                if member['index'] <= 11:
                    origin=plan['receipt_origins'][member['source_id']]
                    item=dict(path=origin['primary_path'],sha256=origin['sha256'])
                    source_receipt_gate(plan,member,item,expected_plan_sha,verify_body=False)
                else:
                    item=acquire(plan,member,monitor,expected_plan_sha,start)
                master['source_receipts'].append(item)
            identity()
            master.update(status='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED',full_current_body_SHA_MD5_size_verified=True,
                          scientific_replication_admitted=False,gzip_full_EOF_or_pipeline_replay_certified=False,
                          reused_v3_source_count=11,newly_acquired_source_count=2,old_failed_family_never_relabelled=True,
                          receipt_origins=plan['receipt_origins'])
        except BaseException as error:
            master.update(status='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT',stop_reason=type(error).__name__+': '+str(error))
        master.update(completed_utc=utc(),elapsed_seconds=time.monotonic()-continuation_start,elapsed_original_family_seconds=time.monotonic()-start,completed_source_count=len(master['source_receipts']),resource_after=safe_state(monitor))
        primaries=[PACKAGE/'logs/validation_raw_acquisition_family_receipt_v4_4.json',FOLDER/'validation_raw_acquisition_family_receipt_v4_4.json']
        try:
            intended_primary_hash=hashlib.sha256((json.dumps(master,indent=2)+'\n').encode()).hexdigest()
            fixed_primary_hashes={str(path):intended_primary_hash for path in primaries}
            for path in primaries:write_new(path,master)
            primary_hashes={str(path):hashes(regular(path))[1] for path in primaries}
            if primary_hashes!=fixed_primary_hashes:raise RuntimeError('TWO_PRIMARY_FAMILY_RECEIPTS_DIFFER_FROM_FROZEN_INTENT')
            if master['status']=='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED':
                intended_seal=dict(status='REVIEWED_STAGE_TERMINAL_SEAL',binding=binding,
                    result_receipt_sha256=fixed_primary_hashes,pending_path=str(plan['pending_path']),
                    pending_sha256=terminal.pending_sha,success_requires_absent_PENDING_and_all_failure_addenda=True)
                fixed_terminal_intent_sha=hashlib.sha256((json.dumps(intended_seal,indent=2,allow_nan=False)+'\n').encode()).hexdigest()
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
    p.add_argument('--expected-plan-sha256')
    p.add_argument('--diagnosis-seal');p.add_argument('--diagnosis-seal-sha256')
    p.add_argument('--diagnosis-root');p.add_argument('--diagnosis-root-sha256');args=p.parse_args()
    if args.prepare:prepare(args.diagnosis_seal,args.diagnosis_seal_sha256,args.diagnosis_root,args.diagnosis_root_sha256)
    else:
        if not args.expected_plan_sha256:p.error('Exact independently reviewed plan SHA required')
        execute(args.expected_plan_sha256)
