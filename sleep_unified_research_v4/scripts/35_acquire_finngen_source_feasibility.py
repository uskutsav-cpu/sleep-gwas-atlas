#!/usr/bin/env python3
"""Acquire one bounded, generation-pinned public source; no scientific fitting."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import plistlib
import shutil
import subprocess
import time
import importlib.util

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
URL = 'https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_INSOMNIA.gz?generation=1777989563097164'
SIZE = 809346932
MD5 = 'f16c21acbf8b9ebfeb0f26a19f0ecfe3'
GENERATION = '1777989563097164'
VOLUME = '77FD98CC-B09E-3DAA-ACC0-82FF2B776B28'


def utc():
    return datetime.now(timezone.utc).isoformat()


def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(b); sha.update(b)
    return md5.hexdigest(), sha.hexdigest()


def main():
    spec = importlib.util.spec_from_file_location('native_monitor', PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py')
    monitor = importlib.util.module_from_spec(spec); spec.loader.exec_module(monitor)
    info = plistlib.loads(subprocess.run(['diskutil', 'info', '-plist', '/Volumes/Extreme SSD'], capture_output=True, check=True).stdout)
    assert info['VolumeUUID'].upper() == VOLUME and os.path.ismount('/Volumes/Extreme SSD')
    assert not info.get('ReadOnlyVolume', False)
    folder = SSD / 'new_source_feasibility/finngen_R13_F5_INSOMNIA'
    folder.mkdir(parents=True, exist_ok=False)
    body = folder / 'finngen_R13_F5_INSOMNIA.gz.partial'
    final = folder / 'finngen_R13_F5_INSOMNIA.gz'
    header = folder / 'download_headers.txt'
    stdout = folder / 'curl_stdout.log'
    receipt_path = PACKAGE / 'logs/finngen_R13_insomnia_acquisition_receipt_v4.json'
    assert not receipt_path.exists()
    sources = [PACKAGE / 'FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md',
               PACKAGE / 'tables/sleep_clinician_finngen_MVP_admissibility_v4.tsv',
               PACKAGE / 'source_provenance/sleep_clinician_finngen_R13_insomnia_HEAD_v4.json', Path(__file__)]
    receipt = {'started_utc': utc(), 'url': URL, 'expected_size': SIZE, 'expected_md5': MD5,
               'expected_generation': GENERATION, 'volume_uuid': VOLUME, 'maximum_new_source_bytes': 2 * 1024**3,
               'scope': 'PUBLIC_SOURCE_ACQUISITION_ONLY_NO_ESTIMATION_NO_PAIR_INFERENCE',
               'bound_sources': {str(p): hashes(p)[1] for p in sources}, 'resource_before': monitor.snapshot()}
    command = ['curl', '--fail', '--location', '--proto', '=https', '--tlsv1.2', '--max-time', '3600',
               '--dump-header', str(header), '--output', str(body), URL]
    receipt['command'] = command
    proc = None; samples = []; reason = None; teardown = None; peak = 0
    started = time.monotonic()
    try:
        assert monitor.final_limits(monitor.snapshot(), 0, 0, 0) is None
        with stdout.open('x') as log:
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                    env={**os.environ, 'TMPDIR': str(SSD / 'tmp')})
            assert os.getpgid(proc.pid) == proc.pid
            while proc.poll() is None:
                state = monitor.snapshot(); rss, pids = monitor.owned_rss(proc.pid); peak = max(peak, rss)
                size = body.stat().st_size if body.exists() else 0
                samples.append({**state, 'owned_rss_bytes': rss, 'owned_pids': pids, 'source_bytes': size})
                reason = monitor.final_limits(state, rss, size, time.monotonic() - started)
                if size > SIZE:
                    reason = 'SOURCE_EXCEEDS_FROZEN_BYTE_COUNT'
                if time.monotonic() - started > 3600:
                    reason = 'SOURCE_DOWNLOAD_DEADLINE'
                if reason:
                    teardown = monitor.terminate_owned(proc); raise RuntimeError(reason)
                if len(samples) % 30 == 0:
                    print(json.dumps({'downloaded_bytes': size, 'expected_bytes': SIZE, 'elapsed_seconds': round(time.monotonic() - started, 1)}), flush=True)
                time.sleep(2)
        teardown = monitor.terminate_owned(proc)
        assert proc.returncode == 0 and not teardown['initial_group_members']
        assert body.stat().st_size == SIZE
        text = header.read_text()
        assert 'x-goog-generation: ' + GENERATION in text.lower()
        md5, sha256 = hashes(body)
        assert md5 == MD5
        assert monitor.final_limits(monitor.snapshot(), peak, body.stat().st_size, time.monotonic() - started) is None
        assert all(hashes(p)[1] == h for p, h in receipt['bound_sources'].items())
        body.rename(final)
        receipt.update({'status': 'EXACT_GENERATION_SIZE_MD5_ACQUIRED', 'source_path': str(final),
                        'actual_size': final.stat().st_size, 'actual_md5': md5, 'actual_sha256': sha256,
                        'gzip_full_stream_verified': False, 'scientific_source_admitted': False})
    except BaseException as exc:
        reason = reason or type(exc).__name__ + ': ' + str(exc)
        receipt['status'] = 'FAILED_PRESERVED_NO_RETRY'
        raise
    finally:
        if proc is not None and (proc.poll() is None or monitor.group_members(proc.pid)):
            teardown = monitor.terminate_owned(proc)
        receipt.update({'completed_utc': utc(), 'elapsed_seconds': time.monotonic() - started,
                        'returncode': proc.returncode if proc else None, 'stop_reason': reason,
                        'teardown': teardown, 'peak_observed_owned_rss_bytes': peak,
                        'resource_after': monitor.snapshot(), 'samples': samples,
                        'headers_sha256': hashes(header)[1] if header.exists() else None})
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        (folder / receipt_path.name).write_bytes(receipt_path.read_bytes())
    print(json.dumps({k: receipt[k] for k in ('status', 'actual_size', 'actual_md5', 'actual_sha256', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
