#!/usr/bin/env python3
"""Pinned byte-range transfer after a preserved slow sequential attempt.

No phenotype analyses are performed. Chunk files and failure evidence remain
in a new task-owned SSD directory; admission requires original full-file SHA.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
DEST = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/exact_lipid_sources_ranged_v2')
TRAITS = {'hdl', 'ldl', 'triglycerides'}
CHUNK = 64*1024**2
WORKERS = 8

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024**2), b''):
            h.update(b)
    return h.hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')

def headers(path):
    text = path.read_text(errors='replace').replace('\r\n', '\n')
    block = next((b for b in reversed(text.strip().split('\n\n')) if b.startswith('HTTP/')), '')
    out = {'status_line': block.splitlines()[0] if block else ''}
    for line in block.splitlines()[1:]:
        if ':' in line:
            k, v = line.split(':', 1)
            out[k.lower()] = v.strip()
    return out

def fetch_chunk(source, etag, start, end):
    directory = DEST/source['trait_ids']/'chunks'
    stem = f'{start:012d}-{end:012d}'
    body = directory/(stem+'.bin')
    header = directory/(stem+'.headers')
    receipt = directory/(stem+'.json')
    if any(p.exists() for p in [body, header, receipt]):
        raise RuntimeError('PRIOR_CHUNK_PRESERVED_NO_OVERWRITE')
    expected = end-start+1
    began = time.time()
    process = subprocess.run(['curl', '--silent', '--show-error', '--fail', '--location',
                              '--connect-timeout', '20', '--max-time', '600',
                              '--max-filesize', str(expected), '--range', f'{start}-{end}',
                              '--header', 'Accept-Encoding: identity', '--header', 'If-Match: '+etag,
                              '--dump-header', str(header), '--output', str(body), source['download_url']],
                             text=True, capture_output=True)
    n = body.stat().st_size if body.exists() else 0
    h = headers(header) if header.exists() else {}
    passed = (process.returncode == 0 and n == expected and
              h.get('status_line', '').split()[1:2] == ['206'] and
              h.get('content-range') == f"bytes {start}-{end}/{source['archive_bytes']}" and
              h.get('content-length') == str(expected) and h.get('etag') == etag)
    result = {'start': start, 'end': end, 'expected_bytes': expected, 'actual_bytes': n,
              'curl_returncode': process.returncode, 'stderr': process.stderr, 'headers': h,
              'body_path': str(body), 'body_sha256': sha(body) if passed else None,
              'header_sha256': sha(header) if header.exists() else None,
              'elapsed_seconds': time.time()-began, 'status': 'PASS_EXACT_RANGE' if passed else 'FAIL_PRESERVED'}
    save(receipt, result)
    if not passed:
        raise RuntimeError('RANGE_TRANSFER_OR_IDENTITY_FAILED: '+str(receipt))
    return result

def main():
    registry = ROOT/'config/public_gwas_sources.tsv'
    with registry.open(newline='') as f:
        selected = [r for r in csv.DictReader(f, delimiter='\t') if r['trait_ids'] in TRAITS]
    if len(selected) != 3 or {s['trait_ids'] for s in selected} != TRAITS:
        raise SystemExit('EXACT_UNIQUE_THREE_SOURCE_GATE_FAILED')
    head_path = PACKAGE/'logs/exact_lipid_upstream_head_v2.json'
    upstream = json.loads(head_path.read_text())
    identity = lambda s: (s['trait_ids'], s['download_url'], int(s['archive_bytes']), s['archive_sha256'])
    actual = [(s['trait_id'], s['url'], s['expected_bytes'], s['expected_sha256']) for s in upstream['sources']]
    if len(actual) != 3 or set(actual) != {identity(s) for s in selected}:
        raise SystemExit('HEAD_RECEIPT_IDENTITY_GATE_FAILED')
    for r in upstream['sources']:
        p = PACKAGE/'logs'/(r['trait_id']+'_upstream_head_v2.txt')
        if r['returncode'] or not r['size_matches'] or sha(p) != r['header_sha256'] or r['headers'].get('accept-ranges') != 'bytes':
            raise SystemExit('HEAD_HEADER_HASH_OR_RANGE_GATE_FAILED')
    search_path = PACKAGE/'logs/missing_core_ssd_filename_search_v2.json'
    search = json.loads(search_path.read_text())
    if any(e['path'] != '/Volumes/Extreme SSD/.TemporaryItems' or 'Operation not permitted' not in e['error'] for e in search['errors']):
        raise SystemExit('UNREVIEWED_SEARCH_COVERAGE_ERROR')
    if any(Path(r['path']).name in {s['archive_name'] for s in selected} | {s['raw_files'] for s in selected} for r in search['candidates']):
        raise SystemExit('LOCAL_LIPID_CANDIDATE_REQUIRES_REVIEW')
    probe_header = PACKAGE/'logs/lipid_range_probe_headers_v2.txt'
    probe_body = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/range_probe_v2/hdl_first_1MiB.bin')
    probe = headers(probe_header)
    hdl = next(r for r in upstream['sources'] if r['trait_id'] == 'hdl')
    if (probe.get('status_line', '').split()[1:2] != ['206'] or probe.get('content-range') != 'bytes 0-1048575/2257306867' or
            probe.get('etag') != hdl['headers']['etag'] or probe_body.stat().st_size != 1048576):
        raise SystemExit('RANGE_PROBE_GATE_FAILED')
    sequential_path = PACKAGE/'logs/hdl_source_acquisition_receipt_v2.json'
    sequential = json.loads(sequential_path.read_text())
    if sequential['status'] != 'FAILED_PARTIAL_OR_MISMATCH_PRESERVED' or sequential['actual_sha256'] is not None:
        raise SystemExit('PREVIOUS_TRANSFER_STATE_REQUIRES_REVIEW')
    previous_partial = Path(sequential['path'])
    if previous_partial.stat().st_size != sequential['actual_bytes']:
        raise SystemExit('PREVIOUS_PARTIAL_STATE_CHANGED')
    curl_version = subprocess.check_output(['curl', '--version'], text=True).splitlines()[0]
    if tuple(map(int, curl_version.split()[1].split('.'))) < (8,4,0):
        raise SystemExit('CURL_STREAMING_SIZE_CAP_SUPPORT_REQUIRED')
    if DEST.exists():
        raise SystemExit('PRIOR_RANGED_DESTINATION_PRESERVED_NO_OVERWRITE')
    DEST.mkdir(parents=True)
    total = sum(int(s['archive_bytes']) for s in selected)
    footprint = 2*total + previous_partial.stat().st_size + probe_body.stat().st_size
    if shutil.disk_usage(DEST).free < footprint + 5*1024**3:
        raise SystemExit('SSD_RESOURCE_GATE_FAILED')
    plan = {'frozen_utc': now(), 'expected_new_network_bytes': total,
            'previous_transfer_and_probe_bytes': previous_partial.stat().st_size+probe_body.stat().st_size,
            'retention_footprint_bytes_including_chunks_and_previous_partial': footprint,
            'workers': WORKERS, 'chunk_bytes': CHUNK, 'maximum_seconds_per_chunk': 600,
            'source_order': [s['trait_ids'] for s in selected], 'sources': selected,
            'curl_version': curl_version, 'TLS_verification_disabled': False,
            'destination': str(DEST), 'SSD_free_bytes': shutil.disk_usage(DEST).free,
            'internal_source_bytes': 0, 'hash_and_assembly_buffer_bytes': 4*1024**2,
            'script_sha256': sha(__file__), 'source_registry_sha256': sha(registry),
            'HEAD_receipt_sha256': sha(head_path), 'filename_search_receipt_sha256': sha(search_path),
            'accepted_search_coverage_gaps': search['errors'], 'probe_header_sha256': sha(probe_header),
            'probe_body_sha256': sha(probe_body), 'sequential_failure_receipt_sha256': sha(sequential_path),
            'original_SSD_and_sealed_v1_modified': False, 'native_analyses_or_pair_outcomes_accessed': False}
    plan_path = PACKAGE/'manifests/exact_lipid_ranged_resource_plan_v2.json'
    save(plan_path, plan)
    print(json.dumps({'plan': str(plan_path), 'network_bytes': total, 'workers': WORKERS, 'retention_bytes': footprint}), flush=True)
    acquired = []
    began = time.time()
    for source in selected:
        directory = DEST/source['trait_ids']
        (directory/'chunks').mkdir(parents=True)
        n = int(source['archive_bytes'])
        etag = next(r['headers']['etag'] for r in upstream['sources'] if r['trait_id'] == source['trait_ids'])
        intervals = [(a, min(a+CHUNK, n)-1) for a in range(0, n, CHUNK)]
        results = []
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futures = [pool.submit(fetch_chunk, source, etag, a, b) for a,b in intervals]
            for future in as_completed(futures):
                results.append(future.result())
                print(json.dumps({'trait': source['trait_ids'], 'chunks_complete': len(results), 'chunks_total': len(intervals)}), flush=True)
        ordered = sorted(results, key=lambda r: r['start'])
        if [(r['start'], r['end']) for r in ordered] != intervals or sum(r['actual_bytes'] for r in ordered) != n:
            raise SystemExit('CHUNK_COMPLETENESS_OR_ORDER_GATE_FAILED')
        partial = directory/(source['archive_name']+'.assembling')
        final = directory/source['archive_name']
        if partial.exists() or final.exists():
            raise SystemExit('PRIOR_ASSEMBLY_PRESERVED')
        digest = hashlib.sha256()
        with partial.open('xb') as output:
            for chunk in ordered:
                path = Path(chunk['body_path'])
                if sha(path) != chunk['body_sha256']:
                    raise SystemExit('CHUNK_CHANGED_AFTER_RECEIPT')
                with path.open('rb') as f:
                    for block in iter(lambda: f.read(4*1024**2), b''):
                        output.write(block)
                        digest.update(block)
        passed = partial.stat().st_size == n and digest.hexdigest() == source['archive_sha256']
        if passed:
            partial.rename(final)
        result = {'completed_utc': now(), 'trait_id': source['trait_ids'], 'source_id': source['source_id'],
                  'url': source['download_url'], 'expected_bytes': n, 'actual_bytes': n,
                  'expected_sha256': source['archive_sha256'], 'actual_sha256': digest.hexdigest(),
                  'path': str(final if passed else partial), 'chunks': ordered,
                  'resource_plan_sha256': sha(plan_path), 'status': 'EXACT_HISTORICAL_SOURCE_REACQUIRED' if passed else 'FULL_SHA_MISMATCH_PRESERVED',
                  'native_preprocessing_or_estimation_completed': False}
        save(PACKAGE/'logs'/(source['trait_ids']+'_ranged_acquisition_receipt_v2.json'), result)
        print(json.dumps({'trait': source['trait_ids'], 'status': result['status'], 'sha256': result['actual_sha256']}), flush=True)
        if not passed:
            raise SystemExit('FULL_ORIGINAL_SHA_GATE_FAILED_NO_DOWNSTREAM_USE')
        acquired.append({k:v for k,v in result.items() if k != 'chunks'})
    save(PACKAGE/'logs/exact_lipid_ranged_acquisition_receipt_v2.json',
         {'completed_utc': now(), 'elapsed_seconds': time.time()-began, 'sources': acquired,
          'all_exact_expected_hashes': True, 'total_source_bytes': total, 'resource_plan_sha256': sha(plan_path),
          'original_SSD_and_v1_outputs_modified': False, 'native_reproduction_completed': False})

if __name__ == '__main__':
    main()
