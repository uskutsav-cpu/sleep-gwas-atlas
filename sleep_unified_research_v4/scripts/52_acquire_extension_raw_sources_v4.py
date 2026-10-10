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

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = PACKAGE / 'manifests/extension_raw_acquisition_plan_v4_4.json'
ADMISSION = PACKAGE / 'manifests/extension_raw_acquisition_admission_v4_4.json'
FOLDER = SSD / 'extension_raw_replay_v4'
MONITOR = PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py'
CURL = Path('/usr/bin/curl')
OLD_PLAN = PACKAGE / 'manifests/extension_raw_acquisition_plan_v4.json'
OLD_PLAN_SHA = 'd8d294a6ff8d7b4f613db7e8b87853df975a32df972fed1999789b39a7f10788'


def utc():
    return datetime.now(timezone.utc).isoformat()


def safe_print(value):
    try:
        print(json.dumps(value), flush=True)
    except BaseException:
        pass


TERMINATION_REQUEST = []


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


def prepare():
    assert hashes(OLD_PLAN)[1] == OLD_PLAN_SHA
    old = json.loads(OLD_PLAN.read_text())
    old_family = PACKAGE / 'logs/extension_raw_acquisition_family_receipt_v4.json'
    failed = json.loads(old_family.read_text())
    assert failed['status'] == 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
    assert failed['completed_source_count'] == 0
    first_id = old['members'][0]['extension_trait_id']
    old_receipt = SSD / 'extension_raw_replay/receipts' / (first_id + '.json')
    failure = json.loads(old_receipt.read_text())
    assert failure['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
    assert failure['teardown']['remaining_group_members'] == []
    original_partial = Path(old['members'][0]['body_path'] + '.partial')
    prefix_size = original_partial.stat().st_size
    assert 0 < prefix_size < old['members'][0]['expected_bytes']
    prefix_sha = hashes(original_partial)[1]
    members = []
    bound = dict(old['bound_sources'])
    for p in (OLD_PLAN, old_family, old_receipt, MONITOR, Path(__file__), CURL,
              PACKAGE / 'scripts/45_acquire_extension_raw_sources.py',
              PACKAGE / 'scripts/48_acquire_extension_raw_sources_v2.py',
              PACKAGE / 'scripts/50_acquire_extension_raw_sources_v3.py',
              PACKAGE / 'manifests/extension_raw_acquisition_plan_v4_3.json',
              PACKAGE / 'manifests/extension_raw_acquisition_plan_v4_2.json',
              PACKAGE / 'logs/extension_raw_acquisition_explicit_interruption_v4.json',
              PACKAGE / 'FROZEN_EXTENSION_RAW_REPLAY_OPERATIONAL_AMENDMENT_v4.md'):
        bound[str(p)] = hashes(p)[1]
    for row in old['members']:
        original = ROOT / 'discovery_extension/provenance/streaming_receipts' / (row['extension_trait_id'] + '.json')
        record = json.loads(original.read_text())
        source = record['source_verification']
        assert source['verification_status'] == 'PASS'
        assert source['source'] == row['url'] == record['source_url']
        assert source['expected_md5'] == source['observed_md5'] == row['expected_md5']
        assert source['expected_size_bytes'] == source['observed_size_bytes'] == row['expected_bytes']
        assert source['http_version_id'] == row['expected_s3_version_id']
        assert source['http_etag'].strip('"') == row['expected_etag']
        expected_sha = source['observed_sha256']
        assert re.fullmatch('[0-9a-f]{64}', expected_sha)
        bound[str(original)] = hashes(original)[1]
        members.append({**row, 'body_path': str(FOLDER / 'raw' / row['filename']),
                        'expected_sha256': expected_sha,
                        'historical_streaming_receipt': str(original)})
    curl_version = subprocess.run([str(CURL), '-q', '--version'], capture_output=True, text=True, check=True).stdout
    assert curl_version.startswith('curl 8.7.1 ')
    monitor = load_monitor()
    assert (monitor.INTERNAL_FLOOR, monitor.SSD_FLOOR, monitor.RSS_LIMIT, monitor.OUTPUT_LIMIT) == (3*1024**3, 5*1024**3, 2*1024**3, 4*1024**3)
    physical_mount()
    state = monitor.snapshot()
    assert state['internal_free_bytes'] >= old['internal_floor_bytes']
    assert state['ssd_free_bytes'] >= old['SSD_reservation_bytes']
    plan = {k:v for k,v in old.items() if k not in ('members', 'bound_sources', 'prepared_utc', 'resource_at_plan')}
    plan.update({'prepared_utc': utc(), 'members': members, 'bound_sources': bound,
                 'curl_path': str(CURL), 'curl_version': curl_version,
                 'executor_path': str(Path(__file__)), 'monitor_path': str(MONITOR),
                 'monitor_output_limit_bytes': monitor.OUTPUT_LIMIT,
                 'runtime_poll_seconds': 2, 'resource_at_plan': state,
                 'exclusive_family_lock_path': str(SSD / 'extension_raw_acquisition_family.lock'),
                 'explicit_resume': {'source_trait': first_id, 'original_partial_path': str(original_partial),
                                     'prefix_bytes': prefix_size, 'prefix_sha256': prefix_sha,
                                     'prior_failed_receipt': str(old_receipt),
                                     'method': 'COPY_PRESERVED_PREFIX_THEN_PINNED_HTTP206_RANGE_RESUME_FULL_BODY_MD5_AND_ORIGINAL_SHA256'},
                 'version': 4, 'prior_attempt_preserved': True,
                 'additional_preserved_prefix_bytes': prefix_size,
                 'execution_admitted_only_after_independent_operational_review': True})
    write_new(PLAN, plan)
    write_new(SSD / 'manifests' / PLAN.name, plan)
    print(json.dumps({'plan_path': str(PLAN), 'plan_sha256': hashes(PLAN)[1], 'prefix_bytes': prefix_size,
                      'members': len(members), 'all100_original_sha256_bound': True}), flush=True)


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
    trait = member['extension_trait_id']
    final = Path(member['body_path'])
    partial = Path(str(final) + '.partial')
    receipt_path = FOLDER / 'receipts' / (trait + '.json')
    header = FOLDER / 'logs' / (trait + '.headers.txt')
    stdout = FOLDER / 'logs' / (trait + '.curl.log')
    assert all(not p.exists() for p in (final, partial, receipt_path, header, stdout))
    before = monitor.snapshot()
    assert before['internal_free_bytes'] >= plan['internal_floor_bytes']
    assert before['ssd_free_bytes'] >= plan['ssd_floor_bytes']
    assert_bindings(plan, expected_plan_sha)
    resume = plan['explicit_resume'] if member['index'] == 1 else None
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
    print(json.dumps({'source_start': member['index'], 'total': 100, 'trait': trait, 'resume_offset': offset}), flush=True)
    try:
        if resume:
            original = Path(resume['original_partial_path'])
            assert original.stat().st_size == offset and hashes(original)[1] == resume['prefix_sha256']
            with original.open('rb') as source, partial.open('xb') as target:
                shutil.copyfileobj(source, target, 4 * 1024**2)
                target.flush()
                os.fsync(target.fileno())
            assert partial.stat().st_size == offset and hashes(partial)[1] == resume['prefix_sha256']
        assert not TERMINATION_REQUEST, 'DEFERRED_TERMINATION_BEFORE_WORKER'
        with stdout.open('x') as log:
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                    env={**os.environ, 'TMPDIR': str(SSD / 'tmp')})
            assert os.getpgid(proc.pid) == proc.pid
            while proc.poll() is None:
                state = monitor.snapshot()
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
        assert observed['x-amz-version-id'] == member['expected_s3_version_id']
        assert observed['etag'] == member['expected_etag']
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
        receipt.update({'status': 'EXACT_IMMUTABLE_SOURCE_ACQUIRED', 'actual_size': partial.stat().st_size,
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
                        if label == 'retained_partial' and receipt['status'] == 'EXACT_IMMUTABLE_SOURCE_ACQUIRED':
                            assert (md5_now, sha_now, path.stat().st_size) == (
                                member['expected_md5'], member['expected_sha256'], member['expected_bytes'])
                            receipt['final_seal_md5'] = md5_now
                            receipt['final_seal_sha256'] = sha_now
                    except BaseException as exc:
                        receipt[label + '_inspection_error'] = str(exc)
                        receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
                        receipt['stop_reason'] = 'FINAL_BODY_INSPECTION_FAILURE: ' + str(exc)
                if receipt['status'] == 'EXACT_IMMUTABLE_SOURCE_ACQUIRED':
                    try:
                        if resume:
                            original = Path(resume['original_partial_path'])
                            assert original.stat().st_size == offset and hashes(original)[1] == resume['prefix_sha256']
                        assert_bindings(plan, expected_plan_sha)
                        physical_mount()
                        final_state = monitor.snapshot()
                        final_reason = monitor.final_limits(final_state, peak, partial.stat().st_size, time.monotonic()-start)
                        assert final_reason is None, final_reason
                        assert time.monotonic()-start <= plan['per_body_seconds_limit'], 'FINAL_PER_BODY_DEADLINE'
                        assert time.monotonic()-family_start <= plan['family_seconds_limit'], 'FINAL_FAMILY_DEADLINE'
                        assert not TERMINATION_REQUEST, 'DEFERRED_TERMINATION_BEFORE_SEAL'
                        receipt['resource_after'] = final_state
                        receipt['post_cleanup_hash_resource_identity_gates_pass'] = True
                        assert not final.exists()
                        partial.rename(final)
                    except BaseException as exc:
                        receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
                        receipt['stop_reason'] = 'FINAL_SEAL_FAILURE: ' + str(exc)
                write_new(receipt_path, receipt)
                write_new(PACKAGE / 'source_provenance/extension_raw_acquisition_v4_4' / receipt_path.name, receipt)
            finally:
                # Persistence, stdout and catchable termination failures cannot
                # bypass quarantine while group disappearance is unverified.
                if not teardown['teardown_verified'] and proc is not None:
                    retain_lock_until_gone(proc, receipt_path)
        finally:
            if not teardown['teardown_verified'] and proc is not None:
                retain_lock_until_gone(proc, receipt_path)
    if receipt['status'] != 'EXACT_IMMUTABLE_SOURCE_ACQUIRED':
        raise RuntimeError('SOURCE_STOP_PRESERVED: ' + str(receipt.get('stop_reason')))
    print(json.dumps({'source_complete': member['index'], 'sha256': receipt['actual_sha256'],
                      'elapsed_seconds': receipt['elapsed_seconds']}), flush=True)
    return {'path': str(receipt_path), 'sha256': hashes(receipt_path)[1]}


def execute(expected_plan_sha):
    assert not TERMINATION_REQUEST, 'DEFERRED_TERMINATION_BEFORE_FAMILY'
    assert re.fullmatch('[0-9a-f]{64}', expected_plan_sha)
    assert hashes(PLAN)[1] == expected_plan_sha
    plan = json.loads(PLAN.read_text())
    admission = json.loads(ADMISSION.read_text())
    assert admission['execution_admitted'] is True
    assert admission['plan_sha256'] == expected_plan_sha
    assert admission['executor_sha256'] == hashes(Path(__file__))[1]
    review_bindings = admission['independent_review_artifact_sha256']
    assert isinstance(review_bindings, dict) and len(review_bindings) >= 2, 'NONEMPTY_INDEPENDENT_REVIEW_BINDINGS_REQUIRED'
    assert any(path.endswith('.json') for path in review_bindings)
    assert any(path.endswith('.md') for path in review_bindings)
    assert all(re.fullmatch('[0-9a-f]{64}', value) for value in review_bindings.values())
    plan['reviewed_admission_sha256'] = hashes(ADMISSION)[1]
    plan['reviewed_audit_bindings'] = admission['independent_review_artifact_sha256']
    assert len(plan['members']) == len({m['extension_trait_id'] for m in plan['members']}) == len({m['filename'] for m in plan['members']}) == 100
    assert sum(m['expected_bytes'] for m in plan['members']) == plan['compressed_network_bytes']
    assert_bindings(plan, expected_plan_sha)
    monitor = load_monitor()
    assert (monitor.INTERNAL_FLOOR, monitor.SSD_FLOOR, monitor.RSS_LIMIT, monitor.OUTPUT_LIMIT) == (
        plan['internal_floor_bytes'], plan['ssd_floor_bytes'], plan['maximum_owned_transfer_rss_bytes'], plan['monitor_output_limit_bytes'])
    assert subprocess.run([str(CURL), '-q', '--version'], capture_output=True, text=True, check=True).stdout == plan['curl_version']
    physical_mount()
    with family_lock(plan['exclusive_family_lock_path']):
        # A duplicate invocation rejects before any source, receipt or stdout write.
        for child in ('raw', 'logs', 'receipts'):
            (FOLDER / child).mkdir(parents=True, exist_ok=True)
        assert not (FOLDER / 'family_ownership.json').exists(), 'PRIOR_ATTEMPT_REQUIRES_NEW_VERSION_AND_EXPLICIT_AUDIT'
        state = monitor.snapshot()
        assert state['internal_free_bytes'] >= plan['internal_floor_bytes']
        assert state['ssd_free_bytes'] >= plan['SSD_reservation_bytes']
        write_new(FOLDER / 'family_ownership.json', {'owner_pid': os.getpid(), 'started_utc': utc(),
                                                  'plan_sha256': expected_plan_sha, 'immutable_attempt': 4})
        master = {'started_utc': utc(), 'plan_sha256': expected_plan_sha, 'resource_start': state,
                  'source_receipts': [], 'status': 'RUNNING', 'owner_pid': os.getpid()}
        start = time.monotonic()
        try:
            for member in plan['members']:
                assert time.monotonic()-start <= plan['family_seconds_limit']
                master['source_receipts'].append(acquire(plan, member, monitor, expected_plan_sha, start))
            assert_bindings(plan, expected_plan_sha)
            resume = plan['explicit_resume']
            original = Path(resume['original_partial_path'])
            assert original.stat().st_size == resume['prefix_bytes'] and hashes(original)[1] == resume['prefix_sha256']
            physical_mount()
            final_state = monitor.snapshot()
            assert final_state['internal_free_bytes'] >= plan['internal_floor_bytes']
            assert final_state['ssd_free_bytes'] >= plan['ssd_floor_bytes']
            assert time.monotonic()-start <= plan['family_seconds_limit']
            for source_receipt in master['source_receipts']:
                assert hashes(source_receipt['path'])[1] == source_receipt['sha256']
            master['final_resource_and_identity_gates_pass'] = True
            master['status'] = 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED'
        except BaseException as exc:
            master['status'] = 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
            master['stop_reason'] = type(exc).__name__ + ': ' + str(exc)
        finally:
            master.update({'completed_utc': utc(), 'elapsed_seconds': time.monotonic()-start,
                           'resource_after': safe_state(monitor), 'completed_source_count': len(master['source_receipts'])})
            write_new(PACKAGE / 'logs/extension_raw_acquisition_family_receipt_v4_4.json', master)
            write_new(FOLDER / 'extension_raw_acquisition_family_receipt_v4_4.json', master)
        if master['status'] != 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED':
            raise SystemExit('ACQUISITION_STOP_PRESERVED_REQUIRES_REVIEW')


if __name__ == '__main__':
    for name in ('SIGTERM', 'SIGHUP', 'SIGINT'):
        signal.signal(getattr(signal, name), catchable_termination)
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--execute', action='store_true')
    parser.add_argument('--expected-plan-sha256')
    args = parser.parse_args()
    if args.prepare:
        prepare()
    else:
        assert args.expected_plan_sha256, 'EXPLICIT_REVIEWED_PLAN_SHA_REQUIRED'
        execute(args.expected_plan_sha256)
