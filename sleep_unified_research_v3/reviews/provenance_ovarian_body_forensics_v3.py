#!/usr/bin/env python3
"""Independent byte-only source review; no projection or scientific analysis."""
import csv
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
ROOT = OUT.parents[1]
PACKAGE = OUT.parent
SSD = Path('/Volumes/Extreme SSD')
SOURCE_BASE = SSD / 'sleep-unified-research-v1/recovery-2026-10-09/ovarian_only_continuation_v3'
PREFIX = 'provenance_ovarian_body_forensics_v3'
BEGAN = time.monotonic()
PEAK = 0
BUFFER = 65536


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def guard():
    global PEAK
    if shutil.disk_usage(ROOT).free < 128 * 1024**2:
        raise RuntimeError('INTERNAL_EMERGENCY_FLOOR')
    if shutil.disk_usage(SSD).free < 5 * 1024**3:
        raise RuntimeError('SSD_RUNTIME_FLOOR')
    if time.monotonic() - BEGAN > 3600:
        raise RuntimeError('REVIEW_TIME_LIMIT')
    sample = subprocess.run(['ps', '-p', str(os.getpid()), '-o', 'rss='],
                            capture_output=True, text=True, check=True).stdout.strip()
    if not sample:
        raise RuntimeError('SELF_RSS_UNOBSERVABLE')
    PEAK = max(PEAK, int(sample) * 1024)
    if PEAK > 128 * 1024**2:
        raise RuntimeError('SELF_RSS_LIMIT')
    if sum(p.stat().st_size for p in OUT.glob(PREFIX + '*') if p.is_file()) > 1024**2:
        raise RuntimeError('TECHNICAL_REPORT_RETENTION_LIMIT')


def digest(path, concatenated=None):
    h = hashlib.sha256()
    checked = time.monotonic()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(BUFFER), b''):
            h.update(block)
            if concatenated is not None:
                concatenated.update(block)
            if time.monotonic() - checked >= 1:
                guard()
                checked = time.monotonic()
    return h.hexdigest()


def parse_headers(path):
    blocks = re.split(r'\r?\n\r?\n', Path(path).read_text(errors='strict').strip())
    responses = [block for block in blocks if block.startswith('HTTP/')]
    if not responses:
        raise RuntimeError('NO_HTTP_HEADER_BLOCK')
    lines = responses[-1].splitlines()
    fields = {}
    for line in lines[1:]:
        if ':' in line:
            key, value = line.split(':', 1)
            fields[key.strip().lower()] = value.strip()
    return int(lines[0].split()[1]), fields


def main():
    guard()
    if shutil.disk_usage(ROOT).free < 256 * 1024**2:
        raise RuntimeError('REVIEW_LAUNCH_FLOOR')
    registry = ROOT / 'config/public_gwas_sources.tsv'
    if digest(registry) != '4b194d25641c2308a84428a7d127969e2eeec0b1884626ab0b8f5a49e1b7cc93':
        raise RuntimeError('FROZEN_REGISTRY_CHANGED')
    with registry.open(newline='') as handle:
        sources = {r['trait_ids']: r for r in csv.DictReader(handle, delimiter='\t')
                   if r['trait_ids'] in {'ovarian_cancer'}}
    if set(sources) != {'ovarian_cancer'}:
        raise RuntimeError('SOURCE_SCOPE_CARDINALITY')
    plan = {'frozen_utc': now(), 'script_sha256': digest(__file__),
            'scope': 'byte-only independent ovarian body and saved chunks review; no projection or association values',
            'network_bytes': 0, 'workers': 1, 'hash_buffer_bytes': BUFFER,
            'maximum_large_bytes_read': 2 * sum(int(r['archive_bytes']) for r in sources.values()),
            'internal_launch_bytes': 256 * 1024**2, 'internal_emergency_bytes': 128 * 1024**2,
            'sampled_reviewer_RSS_stop_bytes': 128 * 1024**2, 'SSD_runtime_floor_bytes': 5 * 1024**3,
            'maximum_technical_report_bytes': 1024**2, 'deadline_seconds': 3600,
            'deadline_clock': 'monotonic', 'source_files_modified': False,
            'h2_rg_or_association_values_parsed': False,
            'internal_free_before_bytes': shutil.disk_usage(ROOT).free,
            'SSD_free_before_bytes': shutil.disk_usage(SSD).free}
    plan_path = OUT / (PREFIX + '_resource_plan.json')
    save(plan_path, plan)
    print(json.dumps({'status': 'BYTE_REVIEW_PLAN_FROZEN', 'plan_sha256': digest(plan_path)}), flush=True)
    summary, rows = [], []
    for trait in ['ovarian_cancer']:
        s = sources[trait]
        n = int(s['archive_bytes'])
        receipt_path = PACKAGE / 'logs' / (trait + '_acquisition_receipt_v3.json')
        receipt_hash = digest(receipt_path)
        receipt = json.loads(receipt_path.read_text())
        if (receipt['trait_id'] != trait or receipt['expected_sha256'] != s['archive_sha256']
                or receipt['expected_bytes'] != n or receipt['url'] != s['download_url']):
            raise RuntimeError('RECEIPT_REGISTRY_IDENTITY_MISMATCH')
        body = Path(receipt['path'])
        expected_name = s['archive_name'] + ('.assembling' if trait == 'prostate_cancer' else '')
        if body.parent != SOURCE_BASE / trait or body.name != expected_name:
            raise RuntimeError('BODY_PATH_SCOPE_MISMATCH')
        body_before = body.stat()
        if body_before.st_size != n:
            raise RuntimeError('BODY_EXACT_SIZE_FAILED')
        full_sha = digest(body)
        ordered = sorted(receipt['chunks'], key=lambda r: r['start'])
        concatenated = hashlib.sha256()
        cursor = 0
        head_path = PACKAGE / 'logs' / (trait.split('_')[0] + '_public_HEAD_v3.txt')
        head_status, identity = parse_headers(head_path)
        if head_status != 200 or identity.get('content-length') != str(n):
            raise RuntimeError('FROZEN_HEAD_SIZE_MISMATCH')
        for i, chunk in enumerate(ordered):
            start, end = chunk['start'], chunk['end']
            if start != cursor or not start <= end < n or end - start + 1 > 64 * 1024**2:
                raise RuntimeError('CHUNK_COVERAGE_INVALID')
            length = end - start + 1
            stem = '%012d-%012d' % (start, end)
            chunk_body = SOURCE_BASE / trait / 'chunks' / (stem + '.bin')
            chunk_headers = chunk_body.with_suffix('.headers')
            chunk_json = chunk_body.with_suffix('.json')
            if Path(chunk['body_path']) != chunk_body or json.loads(chunk_json.read_text()) != chunk:
                raise RuntimeError('CHUNK_RECEIPT_BINDING_INVALID')
            status, observed = parse_headers(chunk_headers)
            validator = 'etag' if identity.get('etag') else 'last-modified'
            actual_chunk_sha = digest(chunk_body, concatenated)
            header_sha = digest(chunk_headers)
            passed = (chunk['status'] == 'PASS_EXACT_RANGE' and chunk['curl_returncode'] == 0
                      and chunk_body.stat().st_size == length == chunk['actual_bytes']
                      and chunk['expected_bytes'] == length and status == 206
                      and observed.get('content-range') == 'bytes %d-%d/%d' % (start, end, n)
                      and observed.get('content-length') == str(length)
                      and observed.get(validator) == identity.get(validator)
                      and actual_chunk_sha == chunk['body_sha256']
                      and header_sha == chunk['header_sha256'])
            if not passed:
                raise RuntimeError('INDEPENDENT_CHUNK_BINDING_FAILED')
            rows.append({'trait_id': trait, 'start': start, 'end': end, 'bytes': length,
                         'chunk_sha256': actual_chunk_sha, 'header_sha256': header_sha,
                         'chunk_receipt_sha256': digest(chunk_json), 'HTTP_status': status,
                         'identity_validator': validator, 'all_independent_checks_pass': True})
            cursor = end + 1
            if (i + 1) % 8 == 0:
                print(json.dumps({'trait': trait, 'independent_chunks_verified': i + 1}), flush=True)
        if cursor != n or concatenated.hexdigest() != full_sha:
            raise RuntimeError('ASSEMBLED_BODY_NOT_EXACT_CHUNK_CONCATENATION')
        body_after = body.stat()
        if ((body_before.st_size, body_before.st_mtime_ns) != (body_after.st_size, body_after.st_mtime_ns)
                or digest(receipt_path) != receipt_hash):
            raise RuntimeError('REVIEWED_INPUT_CHANGED')
        matches = full_sha == s['archive_sha256']
        expected_state = 'EXACT_HISTORICAL_SOURCE_REACQUIRED' if matches else 'SHA_MISMATCH_PRESERVED'
        if receipt['actual_sha256'] != full_sha or receipt['status'] != expected_state:
            raise RuntimeError('INDEPENDENT_BODY_RECEIPT_DISAGREEMENT')
        summary.append({'trait_id': trait, 'source_id': s['source_id'], 'body_path': str(body),
                        'bytes': n, 'frozen_expected_sha256': s['archive_sha256'],
                        'independent_body_sha256': full_sha,
                        'independent_concatenated_chunk_sha256': concatenated.hexdigest(),
                        'chunk_count': len(ordered), 'all_chunk_ranges_and_bindings_pass': True,
                        'body_equals_chunk_concatenation': True, 'matches_frozen_historical_SHA': matches,
                        'classification': expected_state, 'receipt_sha256': receipt_hash,
                        'HEAD_sha256': digest(head_path), 'independent_source_admission': matches})
        print(json.dumps({'trait': trait, 'classification': expected_state, 'SHA': full_sha}), flush=True)
    historical_support = None
    guard()
    result = {'completed_utc': now(), 'scope': plan['scope'], 'plan_sha256': digest(plan_path),
              'sources': summary, 'historical_raw_review_performed': False,
              'observed_reviewer_peak_RSS_bytes': PEAK, 'elapsed_seconds': time.monotonic() - BEGAN,
              'hash_buffer_bytes': BUFFER, 'network_bytes': 0, 'source_files_modified': False,
              'association_values_or_h2_rg_parsed': False,
              'source_identity_interpretation': 'Complete body SHA, chunk concatenation SHA, exact ranges and frozen historical SHA must agree.'}
    save(OUT / (PREFIX + '.json'), result)
    with (OUT / (PREFIX + '_chunks.tsv')).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({'status': 'INDEPENDENT_BYTE_REVIEW_COMPLETE', 'report': PREFIX + '.json',
                      'peak_RSS_bytes': PEAK, 'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    try:
        main()
    except BaseException as error:
        path = OUT / (PREFIX + '_FAILURE.json')
        if not path.exists():
            save(path, {'completed_utc': now(), 'error': str(error),
                        'source_files_modified': False, 'source_admission': False})
        raise
