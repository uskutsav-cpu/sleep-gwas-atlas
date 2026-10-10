#!/usr/bin/env python3
"""Sequential exact100 pinned-body acquisition with immutable checkpoints."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import plistlib
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN = PACKAGE / 'manifests/extension_raw_acquisition_plan_v4.json'


def utc():
    return datetime.now(timezone.utc).isoformat()


def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(b); sha.update(b)
    return md5.hexdigest(), sha.hexdigest()


def headers(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith('HTTP/'):
            result = {'status': int(line.split()[1])}
        elif ':' in line:
            k, v = line.split(':', 1); result[k.strip().lower()] = v.strip().strip('"')
    return result


def main():
    plan = json.loads(PLAN.read_text())
    assert len(plan['members']) == len({r['filename'] for r in plan['members']}) == 100
    assert sum(r['expected_bytes'] for r in plan['members']) == plan['compressed_network_bytes']
    info = plistlib.loads(subprocess.run(['diskutil', 'info', '-plist', '/Volumes/Extreme SSD'], capture_output=True, check=True).stdout)
    assert info['VolumeUUID'].upper() == '77FD98CC-B09E-3DAA-ACC0-82FF2B776B28'
    assert os.path.ismount('/Volumes/Extreme SSD') and not info.get('ReadOnlyVolume', False)
    spec = importlib.util.spec_from_file_location('monitor', PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py')
    monitor = importlib.util.module_from_spec(spec); spec.loader.exec_module(monitor)
    assert monitor.snapshot()['internal_free_bytes'] >= plan['internal_floor_bytes']
    assert monitor.snapshot()['ssd_free_bytes'] >= plan['SSD_reservation_bytes']
    folder = SSD / 'extension_raw_replay'
    for child in ('raw', 'logs', 'receipts'):
        (folder / child).mkdir(parents=True, exist_ok=True)
    copied = PACKAGE / 'source_provenance/extension_raw_acquisition_v4'
    copied.mkdir(exist_ok=True)
    master_path = PACKAGE / 'logs/extension_raw_acquisition_family_receipt_v4.json'
    assert not master_path.exists(), 'PRIOR_FAMILY_TERMINAL_REQUIRES_EXPLICIT_AUDIT'
    master = {'started_utc': utc(), 'plan_sha256': hashes(PLAN)[1], 'executor_sha256': hashes(Path(__file__))[1],
              'resource_start': monitor.snapshot(), 'source_receipts': [], 'status': 'RUNNING'}
    started_family = time.monotonic()
    try:
        for member in plan['members']:
            assert all(hashes(p)[1] == h for p, h in plan['bound_sources'].items())
            assert time.monotonic() - started_family <= plan['family_seconds_limit']
            rp = folder / 'receipts' / (member['extension_trait_id'] + '.json')
            final = Path(member['body_path']); partial = Path(str(final) + '.partial')
            if rp.exists():
                previous = json.loads(rp.read_text())
                assert previous['status'] == 'EXACT_IMMUTABLE_SOURCE_ACQUIRED', 'FAILED_RECEIPT_REQUIRES_EXPLICIT_AUDIT'
                assert previous['member'] == member and final.stat().st_size == member['expected_bytes']
                assert hashes(final) == (member['expected_md5'], previous['actual_sha256'])
                master['source_receipts'].append({'path': str(rp), 'sha256': hashes(rp)[1]})
                print(json.dumps({'checkpoint_verified': member['index'], 'trait': member['extension_trait_id']}), flush=True)
                continue
            assert not final.exists() and not partial.exists(), 'UNSEALED_SOURCE_PRESERVED'
            guard = monitor.snapshot()
            assert guard['internal_free_bytes'] >= plan['internal_floor_bytes'] and guard['ssd_free_bytes'] >= plan['ssd_floor_bytes']
            header = folder / 'logs' / (member['extension_trait_id'] + '.headers.txt')
            stdout = folder / 'logs' / (member['extension_trait_id'] + '.curl.log')
            command = ['curl', '--fail', '--location', '--proto', '=https', '--tlsv1.2', '--max-time',
                       str(plan['per_body_seconds_limit']), '--dump-header', str(header), '--output', str(partial), member['url']]
            receipt = {'member': member, 'started_utc': utc(), 'command': command, 'plan_sha256': hashes(PLAN)[1],
                       'executor_sha256': master['executor_sha256'], 'resource_before': guard}
            proc = None; reason = None; samples = []; peak = 0; teardown = None; started = time.monotonic()
            print(json.dumps({'source_start': member['index'], 'total': 100, 'trait': member['extension_trait_id'], 'expected_bytes': member['expected_bytes']}), flush=True)
            try:
                with stdout.open('x') as log:
                    proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                            env={**os.environ, 'TMPDIR': str(SSD / 'tmp')})
                    assert os.getpgid(proc.pid) == proc.pid
                    while proc.poll() is None:
                        state = monitor.snapshot(); rss, pids = monitor.owned_rss(proc.pid); peak = max(peak, rss)
                        size = partial.stat().st_size if partial.exists() else 0
                        reason = monitor.final_limits(state, rss, size, time.monotonic() - started)
                        if size > member['expected_bytes']:
                            reason = 'SOURCE_EXCEEDS_PINNED_SIZE'
                        if time.monotonic() - started > plan['per_body_seconds_limit']:
                            reason = 'PER_BODY_DEADLINE'
                        if time.monotonic() - started_family > plan['family_seconds_limit']:
                            reason = 'FAMILY_DEADLINE'
                        samples.append({**state, 'source_bytes': size, 'owned_rss_bytes': rss, 'owned_pids': pids})
                        if reason:
                            teardown = monitor.terminate_owned(proc); raise RuntimeError(reason)
                        if len(samples) % 30 == 0:
                            print(json.dumps({'source_index': member['index'], 'downloaded_bytes': size, 'elapsed_seconds': round(time.monotonic() - started, 1)}), flush=True)
                        time.sleep(2)
                teardown = monitor.terminate_owned(proc)
                assert proc.returncode == 0 and not teardown['initial_group_members']
                observed = headers(header)
                assert observed['status'] == 200 and observed['x-amz-version-id'] == member['expected_s3_version_id']
                assert observed['etag'] == member['expected_etag'] and int(observed['content-length']) == member['expected_bytes']
                assert partial.stat().st_size == member['expected_bytes']
                md5, sha256 = hashes(partial)
                assert md5 == member['expected_md5']
                assert monitor.final_limits(monitor.snapshot(), peak, partial.stat().st_size, time.monotonic() - started) is None
                assert all(hashes(p)[1] == h for p, h in plan['bound_sources'].items())
                partial.rename(final)
                receipt.update({'status': 'EXACT_IMMUTABLE_SOURCE_ACQUIRED', 'actual_size': final.stat().st_size,
                                'actual_md5': md5, 'actual_sha256': sha256, 'observed_headers': observed,
                                'gzip_full_stream_verified': False, 'source_body_committed': False})
            except BaseException as exc:
                reason = reason or type(exc).__name__ + ': ' + str(exc)
                receipt['status'] = 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
                raise
            finally:
                if proc is not None and (proc.poll() is None or monitor.group_members(proc.pid)):
                    teardown = monitor.terminate_owned(proc)
                receipt.update({'completed_utc': utc(), 'elapsed_seconds': time.monotonic() - started,
                                'returncode': proc.returncode if proc else None, 'stop_reason': reason,
                                'teardown': teardown, 'peak_observed_owned_rss_bytes': peak,
                                'resource_after': monitor.snapshot(), 'samples': samples,
                                'headers_sha256': hashes(header)[1] if header.exists() else None})
                rp.write_text(json.dumps(receipt, indent=2) + '\n')
                (copied / rp.name).write_bytes(rp.read_bytes())
            master['source_receipts'].append({'path': str(rp), 'sha256': hashes(rp)[1]})
            print(json.dumps({'source_complete': member['index'], 'sha256': receipt['actual_sha256'], 'elapsed_seconds': receipt['elapsed_seconds']}), flush=True)
        master['status'] = 'ALL100_EXACT_SOURCE_BODIES_ACQUIRED'
    except BaseException as exc:
        master['status'] = 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
        master['stop_reason'] = type(exc).__name__ + ': ' + str(exc)
        raise
    finally:
        master.update({'completed_utc': utc(), 'elapsed_seconds': time.monotonic() - started_family,
                       'resource_after': monitor.snapshot(), 'completed_source_count': len(master['source_receipts'])})
        master_path.write_text(json.dumps(master, indent=2) + '\n')
        (folder / master_path.name).write_bytes(master_path.read_bytes())


if __name__ == '__main__':
    main()
