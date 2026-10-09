#!/usr/bin/env python3
"""Three predeclared capped byte probes; rejected source never admitted."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

OUT = Path(__file__).resolve().parent
PACKAGE = OUT.parent
ROOT = OUT.parents[1]
SSD = Path('/Volumes/Extreme SSD')
DEST = SSD / 'sleep-unified-research-v1/recovery-2026-10-09/provenance_review_prostate_1MiB_probes_v3'
BUFFER = 65536
PROBE = 1024**2
BEGAN = time.monotonic()
PEAK = 0


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(BUFFER), b''):
            h.update(block)
    return h.hexdigest()


def headers(path):
    responses = [b for b in re.split(r'\r?\n\r?\n', path.read_text().strip()) if b.startswith('HTTP/')]
    if not responses:
        return None, {}
    lines = responses[-1].splitlines()
    fields = {}
    for line in lines[1:]:
        if ':' in line:
            k, v = line.split(':', 1)
            fields[k.lower()] = v.strip()
    return int(lines[0].split()[1]), fields


def guard():
    global PEAK
    if shutil.disk_usage(ROOT).free < 128 * 1024**2 or shutil.disk_usage(SSD).free < 5 * 1024**3:
        raise RuntimeError('PROBE_RESOURCE_FLOOR')
    if time.monotonic() - BEGAN > 300:
        raise RuntimeError('PROBE_TOTAL_TIME_LIMIT')
    sample = subprocess.check_output(['ps', '-p', str(os.getpid()), '-o', 'rss='], text=True).strip()
    if not sample:
        raise RuntimeError('PROBE_SELF_RSS_UNOBSERVABLE')
    PEAK = max(PEAK, int(sample) * 1024)
    if PEAK > 128 * 1024**2:
        raise RuntimeError('PROBE_SELF_RSS_LIMIT')


def main():
    guard()
    if shutil.disk_usage(ROOT).free < 256 * 1024**2 or DEST.exists():
        raise RuntimeError('PROBE_LAUNCH_OR_PRIOR_DESTINATION_GATE')
    forensic_path = OUT / 'provenance_core_body_forensics_v3.json'
    forensic = json.loads(forensic_path.read_text())
    rejected = next(r for r in forensic['sources'] if r['trait_id'] == 'prostate_cancer')
    if rejected['classification'] != 'SHA_MISMATCH_PRESERVED' or not rejected['body_equals_chunk_concatenation']:
        raise RuntimeError('INDEPENDENT_REJECTED_BODY_REQUIRED')
    receipt_path = PACKAGE / 'logs/prostate_cancer_acquisition_receipt_v3.json'
    receipt = json.loads(receipt_path.read_text())
    if sha(receipt_path) != rejected['receipt_sha256']:
        raise RuntimeError('REJECTED_RECEIPT_CHANGED')
    body = Path(rejected['body_path'])
    body_before = body.stat()
    n = rejected['bytes']
    head_path = PACKAGE / 'logs/prostate_public_HEAD_v3.txt'
    status, frozen = headers(head_path)
    etag = frozen.get('etag')
    if status != 200 or not etag or etag.startswith('W/') or sha(head_path) != rejected['HEAD_sha256']:
        raise RuntimeError('STRONG_FROZEN_ETAG_REQUIRED')
    version = subprocess.check_output(['curl', '--version'], text=True).splitlines()[0]
    if tuple(map(int, version.split()[1].split('.'))) < (8, 4, 0):
        raise RuntimeError('CURL_STREAMING_CAP_REQUIRED')
    intervals = [('early', 0), ('middle', n // 2), ('late', n - PROBE)]
    plan = {'frozen_utc': now(), 'script_sha256': sha(__file__), 'forensic_result_sha256': sha(forensic_path),
            'receipt_sha256': sha(receipt_path), 'HEAD_sha256': sha(head_path),
            'public_url': receipt['url'], 'conditional_ETag': etag,
            'rejected_full_body_sha256': rejected['independent_body_sha256'],
            'offset_selection': 'fixed zero, integer half-size, and final 1MiB; declared before any probe response',
            'probes': [{'label': label, 'start': start, 'end': start + PROBE - 1} for label, start in intervals],
            'maximum_response_body_network_bytes': 3 * PROBE,
            'maximum_retained_probe_body_bytes': 3 * PROBE,
            'workers': 1, 'no_retry': True, 'curl_version': version,
            'max_seconds_per_request': 60, 'max_total_seconds': 300,
            'hash_and_comparison_buffer_bytes': BUFFER,
            'internal_launch_bytes': 256 * 1024**2, 'internal_emergency_bytes': 128 * 1024**2,
            'sampled_reviewer_RSS_stop_bytes': 128 * 1024**2,
            'SSD_runtime_floor_bytes': 5 * 1024**3, 'destination': str(DEST),
            'original_files_modified': False, 'association_values_or_h2_rg_parsed': False,
            'scientific_source_admission': False}
    plan_path = OUT / 'provenance_prostate_public_probes_plan_v3.json'
    save(plan_path, plan)
    DEST.mkdir()
    results = []
    for label, start in intervals:
        guard()
        end = start + PROBE - 1
        target, header = DEST / (label + '.bin'), DEST / (label + '.headers')
        command = ['curl', '--silent', '--show-error', '--fail', '--location', '--max-redirs', '3',
                   '--proto', '=https', '--proto-redir', '=https', '--connect-timeout', '20',
                   '--max-time', '60', '--max-filesize', str(PROBE), '--range', '%d-%d' % (start, end),
                   '--header', 'Accept-Encoding: identity', '--header', 'If-Match: ' + etag,
                   '--dump-header', str(header), '--output', str(target), receipt['url']]
        proc = None
        try:
            proc = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            while proc.poll() is None:
                guard()
                time.sleep(1)
        finally:
            if proc is not None and proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
        _, stderr = proc.communicate()
        actual = target.stat().st_size if target.exists() else 0
        http, observed = headers(header) if header.exists() else (None, {})
        valid = (proc.returncode == 0 and actual == PROBE and http == 206
                 and observed.get('content-length') == str(PROBE)
                 and observed.get('content-range') == 'bytes %d-%d/%d' % (start, end, n)
                 and observed.get('etag') == etag)
        reference_digest = hashlib.sha256()
        agrees = None
        if valid:
            agrees = True
            with body.open('rb') as reference, target.open('rb') as sample:
                reference.seek(start)
                for _ in range(PROBE // BUFFER):
                    a, b = reference.read(BUFFER), sample.read(BUFFER)
                    reference_digest.update(a)
                    agrees = agrees and a == b
        result = {'label': label, 'start': start, 'end': end, 'expected_probe_bytes': PROBE,
                  'actual_probe_bytes': actual, 'curl_returncode': proc.returncode, 'stderr': stderr,
                  'HTTP_status': http, 'response_headers': observed, 'exact_range_checks_pass': valid,
                  'probe_sha256': sha(target) if target.exists() else None,
                  'rejected_body_slice_sha256': reference_digest.hexdigest() if valid else None,
                  'probe_bytes_equal_rejected_body_slice': agrees,
                  'header_sha256': sha(header) if header.exists() else None,
                  'scientific_source_admission': False}
        save(DEST / (label + '.json'), result)
        results.append(result)
        print(json.dumps({'probe': label, 'exact_range_pass': valid, 'same_bytes': agrees}), flush=True)
    body_after = body.stat()
    unchanged = ((body_before.st_size, body_before.st_mtime_ns) == (body_after.st_size, body_after.st_mtime_ns)
                 and sha(receipt_path) == rejected['receipt_sha256'])
    guard()
    result = {'completed_utc': now(), 'plan_sha256': sha(plan_path), 'results': results,
              'all_three_exact_ranges_match_rejected_body': all(r['exact_range_checks_pass'] and r['probe_bytes_equal_rejected_body_slice'] for r in results),
              'rejected_body_stat_and_receipt_unchanged': unchanged,
              'response_body_bytes_downloaded': sum(r['actual_probe_bytes'] for r in results),
              'observed_reviewer_peak_RSS_bytes': PEAK, 'elapsed_seconds': time.monotonic() - BEGAN,
              'scientific_source_admission': False, 'original_files_modified': False,
              'interpretation_limit': 'Three sampled current public byte ranges cannot establish full remote-object identity or prove why the historical digest differs. Same HEAD/ETag is not body identity proof.'}
    save(OUT / 'provenance_prostate_public_probes_v3.json', result)
    print(json.dumps({'status': 'PUBLIC_PROBES_COMPLETE', 'all_three_match': result['all_three_exact_ranges_match_rejected_body']}), flush=True)


if __name__ == '__main__':
    try:
        main()
    except BaseException as error:
        path = OUT / 'provenance_prostate_public_probes_FAILURE_v3.json'
        if not path.exists():
            save(path, {'completed_utc': now(), 'error': str(error), 'original_files_modified': False,
                        'scientific_source_admission': False})
        raise
