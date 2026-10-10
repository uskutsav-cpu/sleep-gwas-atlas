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
PLAN = PACKAGE / 'manifests/extension_raw_acquisition_plan_v4_7.json'
ADMISSION = PACKAGE / 'manifests/extension_raw_acquisition_admission_v4_7.json'
FOLDER = SSD / 'extension_raw_replay_v7'
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
    """Freeze one manual continuation after the documented source4 TCP failure."""
    old_path=PACKAGE/'manifests/extension_raw_acquisition_plan_v4_6.json'
    assert hashes(old_path)[1]=='393af7ffc3231ceddf829ae669558001a7411bb4a7917b950a713c41546237fe'
    old=json.loads(old_path.read_text())
    failed_path=PACKAGE/'logs/extension_raw_acquisition_family_receipt_v4_6.json'
    failed=json.loads(failed_path.read_text())
    assert failed['status']=='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT' and failed['completed_source_count']==3
    assert failed['stop_reason']=='RuntimeError: SOURCE_STOP_PRESERVED: AssertionError: CURL_NONZERO_EXIT: 56'
    fourth=old['members'][3]
    failed_source=SSD/'extension_raw_replay_v6/receipts'/(fourth['extension_trait_id']+'.json')
    r=json.loads(failed_source.read_text())
    assert r['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and r['returncode']==56
    assert r['stop_reason']=='AssertionError: CURL_NONZERO_EXIT: 56'
    assert r['teardown']['teardown_verified'] and r['teardown']['remaining_group_members']==[]
    prefix=Path(fourth['body_path']+'.partial')
    assert not prefix.is_symlink() and prefix.stat().st_size==r['retained_partial_bytes']==1126899820
    assert 0<prefix.stat().st_size<fourth['expected_bytes'] and hashes(prefix)[1]==r['retained_partial_sha256']
    resource_ledger=PACKAGE/'manifests/global_SSD_resource_reservation_v4_3.json'
    assert hashes(resource_ledger)[1]=='6bd5c7cec9d9880ae630e5bfd431b548bd4cc2f3837d3ce12935fda326a9c590'
    ledger=json.loads(resource_ledger.read_text())
    assert ledger['reserved_total_bytes']==319645174442<ledger['ceiling_bytes']==old['SSD_reservation_bytes']
    gate_plan=SSD/'extension_pipeline_replay_v3/extension_pipeline_replay_plan_v3.json'
    assert hashes(gate_plan)[1]=='565b34e997b2cab841f3a101cb4f8649c4c4660e56f5b1421ce4f68f29afe4ac'
    review=[PACKAGE/'reviews'/n for n in ['independent_acquisition_v6_preflight_v4.md','independent_acquisition_v6_binding_receipt_v4.json','independent_acquisition_v6_metadata_controls_v4.py','independent_acquisition_v6_preflight_v4.sha256']]
    bound=dict(old['bound_sources'])
    files=[old_path,failed_path,failed_source,Path(old['pending_path']),resource_ledger,gate_plan,PACKAGE/'scripts/extension_replay_common_v3.py',PACKAGE/'scripts/52_acquire_extension_raw_sources_v6.py',Path(__file__),*review]
    for name in ['genomicsem_extension_checkpoint_review_v3.md','genomicsem_extension_checkpoint_review_v3.json']:files.append(PACKAGE/'reviews'/name)
    for f in files:bound[str(f)]=hashes(f)[1]
    reused={};members=[];successful={item['path']:item['sha256'] for item in failed['source_receipts']}
    for m in old['members']:
        if m['index']<=3:
            origin=old['reused_checkpoints'].get(str(m['index']))
            path=Path(origin['receipt_path']) if origin else SSD/'extension_raw_replay_v6/receipts'/(m['extension_trait_id']+'.json')
            receipt=json.loads(path.read_text());assert receipt['status']=='EXACT_IMMUTABLE_SOURCE_ACQUIRED'
            assert successful[str(path)]==hashes(path)[1]
            reused[str(m['index'])]=dict(receipt_path=str(path),receipt_sha256=hashes(path)[1],original_acquisition_plan_sha256=receipt['plan_sha256'])
            bound[str(path)]=hashes(path)[1];members.append(m)
        else:members.append({**m,'body_path':str(FOLDER/'raw'/m['filename'])})
    plan=dict(old)
    plan.update(prepared_utc=utc(),version=7,members=members,bound_sources=bound,executor_path=str(Path(__file__)),reused_checkpoints=reused,reused_checkpoint_gate_plan=str(gate_plan),
        explicit_resume=dict(source_index=4,source_trait=fourth['extension_trait_id'],original_partial_path=str(prefix),prefix_bytes=prefix.stat().st_size,prefix_sha256=r['retained_partial_sha256'],prior_failed_receipt=str(failed_source),method='COPY_PRESERVED_SOURCE4_PREFIX_THEN_PINNED_HTTP206_FULL_ORIGINAL_SHA_MD5'),
        additional_preserved_prefix_bytes=old['additional_preserved_prefix_bytes']+prefix.stat().st_size,
        prior_attempt_preserved=True,network_sources_to_transfer=97,
        new_network_bytes=sum(m['expected_bytes'] for m in members[3:])-prefix.stat().st_size,
        resource_at_plan=load_monitor().snapshot(),global_resource_ledger_path=str(resource_ledger),global_resource_ledger_sha256=hashes(resource_ledger)[1],
        pending_path=str(FOLDER/'family_pending_v4_7.json'),terminal_seal_path=str(FOLDER/'family_terminal_seal_v4_7.json'),
        prior_v6_plan_sha256=hashes(old_path)[1],prior_v6_failed_family_sha256=hashes(failed_path)[1],
        root_manual_continuation_reason='Actual curl56 receive timeout; identical pinned version and prefix retained; no automatic retry or source substitution')
    physical_mount()
    assert plan['resource_at_plan']['internal_free_bytes']>=plan['internal_floor_bytes']
    assert plan['resource_at_plan']['ssd_free_bytes']>=plan['SSD_reservation_bytes']
    for m in members[:3]:reused_checkpoint_gate(plan,m,load_monitor())
    write_new(PLAN,plan);write_new(SSD/'manifests'/PLAN.name,plan)
    safe_print(dict(plan_sha256=hashes(PLAN)[1],original_sources=100,reused_sources=3,network_sources=97,source4_resume_prefix_bytes=prefix.stat().st_size,workers_launched=0))


def global_namespace_gate(plan):
    total=sum(p.stat().st_size for p in SSD.rglob('*') if p.is_file() and not p.is_symlink())
    assert total<=plan['SSD_reservation_bytes'], '300GIB_GLOBAL_NEW_CAMPAIGN_NAMESPACE_LIMIT'
    return total


def reused_checkpoint_gate(plan,member,monitor):
    """Reuse only the two exact receipts; no receipt or source bytes are edited."""
    fixed=plan['reused_checkpoints'][str(member['index'])]
    assert hashes(fixed['receipt_path'])[1]==fixed['receipt_sha256']
    gate_path=Path(plan['reused_checkpoint_gate_plan'])
    assert hashes(gate_path)[1]==plan['bound_sources'][str(gate_path)]
    spec=importlib.util.spec_from_file_location('_unchanged_v4_checkpoint_gate',PACKAGE/'scripts/extension_replay_common_v3.py')
    gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
    old_plan=json.loads(gate_path.read_text())
    old_member=next(m for m in old_plan['members'] if m['index']==member['index'])
    r=gate.acquisition_receipt_gate(old_plan,old_member,fixed['receipt_sha256'])
    assert r['member']==member and r['plan_sha256']==fixed['original_acquisition_plan_sha256']
    body=Path(member['body_path'])
    assert (body.stat().st_size,*hashes(body))==(member['expected_bytes'],member['expected_md5'],member['expected_sha256'])
    state=monitor.snapshot()
    global_namespace_gate(plan)
    assert state['internal_free_bytes']>=plan['internal_floor_bytes'] and state['ssd_free_bytes']>=plan['ssd_floor_bytes']
    return dict(path=fixed['receipt_path'],sha256=fixed['receipt_sha256'])


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
    global_namespace_gate(plan)
    assert before['internal_free_bytes'] >= plan['internal_floor_bytes']
    assert before['ssd_free_bytes'] >= plan['ssd_floor_bytes']
    assert_bindings(plan, expected_plan_sha)
    resume = plan['explicit_resume'] if member['index'] == plan['explicit_resume']['source_index'] else None
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
                        global_namespace_gate(plan)
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
                write_new(PACKAGE / 'source_provenance/extension_raw_acquisition_v4_7' / receipt_path.name, receipt)
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


def fsync_directory(path):
    fd=os.open(path,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)


def finalize_family(plan,expected_plan_sha,master,monitor,start,pending,pending_sha):
    """All success copies are provisional while the durable PENDING exists.

    The final unlink is the commit boundary. No fallible work/check/diagnostic
    follows it. A crash that reintroduces the marker causes conservative rejection.
    Failure receipts are supplemental: their persistence is never the veto oracle.
    """
    primaries=[PACKAGE/'logs/extension_raw_acquisition_family_receipt_v4_7.json',FOLDER/'extension_raw_acquisition_family_receipt_v4_7.json']
    terminal=Path(plan['terminal_seal_path'])
    committed=False
    try:
        master.update(completed_utc=utc(),elapsed_seconds=time.monotonic()-start,resource_after=safe_state(monitor),completed_source_count=len(master['source_receipts']),terminal_protocol=plan['terminal_protocol'],pending_path=str(pending),terminal_seal_path=str(terminal))
        assert pending.is_file() and hashes(pending)[1]==pending_sha, 'DURABLE_PENDING_CHANGED_OR_MISSING'
        for path in primaries:write_new(path,master);fsync_directory(path.parent)
        if master['status']!='ALL100_EXACT_SOURCE_BODIES_ACQUIRED':return False
        assert len(master['source_receipts'])==100
        assert master['reused_exact_source_count']==len(plan['reused_checkpoints']) and master['new_exact_source_count']==100-len(plan['reused_checkpoints'])
        primary_sha={str(path):hashes(path)[1] for path in primaries}
        assert len(set(primary_sha.values()))==1, 'PRIMARY_RECEIPT_COPIES_DIFFER'
        for item in master['source_receipts']:assert hashes(item['path'])[1]==item['sha256']
        seal=dict(status='ALL100_FAMILY_TERMINAL_SEAL',plan_sha256=expected_plan_sha,executor_sha256=hashes(Path(__file__))[1],primary_receipt_sha256=primary_sha,source_receipts=master['source_receipts'],pending_path=str(pending),pending_sha256=pending_sha,created_utc=utc(),full_gzip_EOF_or_pipeline_replay_certified=False)
        write_new(terminal,seal);fsync_directory(terminal.parent)
        terminal_sha=hashes(terminal)[1]
        assert_bindings(plan,expected_plan_sha);physical_mount()
        global_namespace_gate(plan)
        for item in master['source_receipts']:assert hashes(item['path'])[1]==item['sha256']
        for path,h in primary_sha.items():assert hashes(path)[1]==h
        assert hashes(terminal)[1]==terminal_sha and hashes(pending)[1]==pending_sha
        for path in [*primaries,terminal]:assert not Path(str(path)+'.failure.json').exists(), 'TERMINAL_FAILURE_ALREADY_EXISTS'
        state=monitor.snapshot()
        assert state['internal_free_bytes']>=plan['internal_floor_bytes'] and state['ssd_free_bytes']>=plan['ssd_floor_bytes']
        assert time.monotonic()-start<=plan['family_seconds_limit'] and not TERMINATION_REQUEST
        # Commit only after post-persistence identity/resource/deadline checks.
        # Do not fsync/query/print after this boundary: a crash may restore the
        # PENDING marker and reject an otherwise complete transfer conservatively.
        pending.unlink()
        committed=True
    except BaseException as exc:
        master['status']='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
        master.setdefault('stop_reason','FAMILY_FINALIZATION_FAILURE: '+type(exc).__name__+': '+str(exc))
        master['finalization_error']=type(exc).__name__+': '+str(exc)
        for path in [*primaries,terminal]:
            try:write_new(Path(str(path)+'.failure.json'),master)
            except BaseException as secondary:safe_print(dict(event='SUPPLEMENTAL_FAILURE_RECEIPT_SAVE_ERROR',path=str(path),error=repr(secondary),pending_must_veto_all_provisional_success=True))
    return committed


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
        assert_bindings(plan,expected_plan_sha);physical_mount();global_namespace_gate(plan)
        pending=Path(plan['pending_path'])
        terminal=Path(plan['terminal_seal_path'])
        assert not pending.exists() and not terminal.exists(), 'PRIOR_TERMINAL_ATTEMPT_PRESERVED'
        write_new(pending,{'status':'PENDING_NOT_ADMISSIBLE','plan_sha256':expected_plan_sha,'owner_pid':os.getpid(),'started_utc':utc()})
        fsync_directory(pending.parent)
        pending_sha=hashes(pending)[1]
        write_new(FOLDER / 'family_ownership.json', {'owner_pid': os.getpid(), 'started_utc': utc(),
                                                  'plan_sha256': expected_plan_sha, 'immutable_attempt': 6})
        master = {'started_utc': utc(), 'plan_sha256': expected_plan_sha, 'resource_start': state,
                  'source_receipts': [], 'status': 'RUNNING', 'owner_pid': os.getpid()}
        start = time.monotonic()
        try:
            for member in plan['members']:
                assert time.monotonic()-start <= plan['family_seconds_limit']
                if str(member['index']) in plan['reused_checkpoints']:
                    master['source_receipts'].append(reused_checkpoint_gate(plan,member,monitor))
                    safe_print(dict(source_reused=member['index'],original_exact_receipt_preserved=True))
                else:master['source_receipts'].append(acquire(plan, member, monitor, expected_plan_sha, start))
            assert_bindings(plan, expected_plan_sha)
            resume = plan['explicit_resume']
            original = Path(resume['original_partial_path'])
            assert original.stat().st_size == resume['prefix_bytes'] and hashes(original)[1] == resume['prefix_sha256']
            physical_mount()
            final_state = monitor.snapshot()
            master['final_global_campaign_namespace_bytes']=global_namespace_gate(plan)
            assert final_state['internal_free_bytes'] >= plan['internal_floor_bytes']
            assert final_state['ssd_free_bytes'] >= plan['ssd_floor_bytes']
            assert time.monotonic()-start <= plan['family_seconds_limit']
            for source_receipt in master['source_receipts']:
                assert hashes(source_receipt['path'])[1] == source_receipt['sha256']
            assert len(master['source_receipts'])==100
            for m in plan['members'][:len(plan['reused_checkpoints'])]:reused_checkpoint_gate(plan,m,monitor)
            master['reused_exact_source_count']=len(plan['reused_checkpoints'])
            master['new_exact_source_count']=100-len(plan['reused_checkpoints'])
            master['final_resource_and_identity_gates_pass'] = True
            master['status'] = 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED'
        except BaseException as exc:
            master['status'] = 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
            master['stop_reason'] = type(exc).__name__ + ': ' + str(exc)
        finally:
            committed=finalize_family(plan,expected_plan_sha,master,monitor,start,pending,pending_sha)
        if not committed or master['status'] != 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED':
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
