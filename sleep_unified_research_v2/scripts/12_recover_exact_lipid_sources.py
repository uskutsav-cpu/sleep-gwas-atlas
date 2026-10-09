#!/usr/bin/env python3
"""Versioned, bounded recovery of exact historical public lipid source bytes.

Modes run sequentially: metadata search; upstream HEAD; acquisition.
Original SSD files, sealed v1 outputs, and statistical families are read-only.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
SSD = Path('/Volumes/Extreme SSD')
DEST = SSD / 'sleep-unified-research-v1/recovery-2026-10-09/exact_lipid_sources_v2'
TRAITS = {'hdl', 'ldl', 'triglycerides'}

def now():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()

def read(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')

def sources():
    return [r for r in read(ROOT / 'config/public_gwas_sources.tsv') if r['trait_ids'] in TRAITS]

def search_ssd():
    receipt = PACKAGE / 'logs/missing_core_ssd_filename_search_v2.json'
    if receipt.exists():
        raise SystemExit('PRIOR_SEARCH_PRESERVED')
    plan = json.loads((PACKAGE / 'manifests/missing_core_archive_search_plan_v2.json').read_text())
    names = plan['target_names']
    excluded = {'.git', 'node_modules', '.Trashes', '.Spotlight-V100', '.fseventsd',
                '__pycache__', 'site-packages', '.Trash'}
    errors = []
    matches = []
    directories = files = 0
    started = time.time()
    def on_error(exc):
        errors.append({'path': exc.filename, 'error': str(exc)})
    for base, dirs, filenames in os.walk(SSD, followlinks=False, onerror=on_error):
        directories += 1
        dirs[:] = [d for d in dirs if d not in excluded and not d.startswith('._')]
        for name in filenames:
            files += 1
            if name.startswith('._') or name not in names:
                continue
            path = Path(base) / name
            st = path.lstat()
            matches.append({'path': str(path), 'bytes': st.st_size,
                            'expected_sha256': names[name], 'is_symlink': path.is_symlink(),
                            'status': 'FILENAME_CANDIDATE_ONLY_NOT_HASH_VERIFIED'})
        if directories % 20000 == 0:
            print(json.dumps({'directories': directories, 'files': files, 'matches': len(matches)}), flush=True)
    result = {'completed_utc': now(), 'elapsed_seconds': time.time() - started,
              'root': str(SSD), 'directories_scanned': directories, 'files_scanned': files,
              'excluded_directory_names': sorted(excluded), 'follow_symlinks': False,
              'target_names': names, 'candidates': matches, 'errors': errors,
              'archive_members_not_searched_by_this_mode': True,
              'no_content_or_original_files_modified': True, 'script_sha256': sha(__file__)}
    save(receipt, result)
    print(json.dumps({'receipt': str(receipt), 'candidates': matches, 'errors': len(errors)}), flush=True)

def headers(path):
    text = path.read_text(errors='replace')
    blocks = text.replace('\r\n', '\n').strip().split('\n\n')
    block = next((x for x in reversed(blocks) if x.startswith('HTTP/')), '')
    out = {}
    for line in block.splitlines()[1:]:
        if ':' in line:
            key, value = line.split(':', 1)
            out[key.lower()] = value.strip()
    return out

def head():
    receipt = PACKAGE / 'logs/exact_lipid_upstream_head_v2.json'
    if receipt.exists():
        raise SystemExit('PRIOR_HEAD_PRESERVED')
    rows = []
    for source in sources():
        path = PACKAGE / 'logs' / (source['trait_ids'] + '_upstream_head_v2.txt')
        if path.exists():
            raise SystemExit('PRIOR_HEADER_PRESERVED')
        process = subprocess.run(['curl', '--silent', '--show-error', '--fail', '--location',
                                  '--head', '--connect-timeout', '20', '--max-time', '60',
                                  '--header', 'Accept-Encoding: identity',
                                  '--output', str(path), source['download_url']],
                                 text=True, capture_output=True)
        h = headers(path) if path.exists() else {}
        row = {'trait_id': source['trait_ids'], 'url': source['download_url'],
               'expected_bytes': int(source['archive_bytes']), 'expected_sha256': source['archive_sha256'],
               'returncode': process.returncode, 'stderr': process.stderr,
               'headers': h, 'header_sha256': sha(path) if path.exists() else None,
               'size_matches': h.get('content-length') == source['archive_bytes'],
               'source_release': source['source_id']}
        rows.append(row)
        print(json.dumps({'trait': source['trait_ids'], 'returncode': process.returncode,
                          'size_matches': row['size_matches']}), flush=True)
    save(receipt, {'completed_utc': now(), 'sources': rows, 'TLS_verification_disabled': False,
                   'large_response_body_downloaded': False, 'script_sha256': sha(__file__)})

def download():
    selected = sources()
    if len(selected) != 3:
        raise SystemExit('EXACT_THREE_SOURCE_GATE_FAILED')
    search = json.loads((PACKAGE / 'logs/missing_core_ssd_filename_search_v2.json').read_text())
    # Preserve the failed-directory receipt. The only accepted gap is the
    # OS-controlled temporary namespace, never a research source directory.
    accepted_search_errors = [e for e in search['errors']
                              if e['path'] == str(SSD / '.TemporaryItems')
                              and 'Operation not permitted' in e['error']]
    if len(accepted_search_errors) != len(search['errors']) or any(Path(r['path']).name in {s['archive_name'] for s in selected}
                               or Path(r['path']).name in {s['raw_files'] for s in selected}
                               for r in search['candidates']):
        raise SystemExit('LOCAL_CANDIDATES_OR_SEARCH_ERRORS_REQUIRE_REVIEW_BEFORE_DOWNLOAD')
    upstream = json.loads((PACKAGE / 'logs/exact_lipid_upstream_head_v2.json').read_text())
    expected_identity = {(s['trait_ids'], s['download_url'], int(s['archive_bytes']), s['archive_sha256']) for s in selected}
    actual_identity = {(r['trait_id'], r['url'], r['expected_bytes'], r['expected_sha256']) for r in upstream['sources']}
    if len(upstream['sources']) != 3 or actual_identity != expected_identity or not all(
            r['returncode'] == 0 and r['size_matches'] and
            r['headers'].get('content-length') == str(r['expected_bytes'])
            for r in upstream['sources']):
        raise SystemExit('UPSTREAM_IDENTITY_GATE_FAILED')
    curl_version = subprocess.check_output(['curl', '--version'], text=True).splitlines()[0]
    version_tuple = tuple(int(v) for v in curl_version.split()[1].split('.'))
    if version_tuple < (8, 4, 0):
        raise SystemExit('CURL_STREAMING_MAX_FILESIZE_SUPPORT_REQUIRED')
    DEST.mkdir(parents=True, exist_ok=True)
    total = sum(int(r['archive_bytes']) for r in selected)
    free = shutil.disk_usage(DEST).free
    if free < total + 5 * 1024**3:
        raise SystemExit('SSD_RESOURCE_GATE_FAILED')
    resource_plan = {'frozen_utc': now(), 'network_bytes_expected': total,
                     'maximum_retained_source_bytes': total, 'SSD_free_bytes_before': free,
                     'SSD_reserve_bytes': 5 * 1024**3, 'destination': str(DEST),
                     'workers': 1, 'streaming_hash_buffer_bytes': 4 * 1024**2,
                     'maximum_seconds_per_download': 1200, 'body_retry_count': 0,
                     'curl_version': curl_version, 'response_size_cap_per_source': 'exact expected source byte count (--max-filesize)',
                     'accepted_search_coverage_gaps': accepted_search_errors,
                     'source_absence_scope': 'readable SSD namespaces and named archive; OS temporary namespace inaccessible',
                     'internal_disk_source_bytes': 0, 'source_registry_sha256': sha(ROOT/'config/public_gwas_sources.tsv'),
                     'HEAD_receipt_sha256': sha(PACKAGE/'logs/exact_lipid_upstream_head_v2.json'),
                     'filename_search_receipt_sha256': sha(PACKAGE/'logs/missing_core_ssd_filename_search_v2.json'),
                     'script_sha256': sha(__file__), 'sources': selected,
                     'original_SSD_and_sealed_v1_modified': False,
                     'no_preprocessing_estimation_or_new_pair_outcomes': True}
    save(PACKAGE / 'manifests/exact_lipid_download_resource_plan_v2.json', resource_plan)
    print(json.dumps({'resource_plan_frozen': True, 'expected_network_bytes': total, 'SSD_free_bytes': free}), flush=True)
    rows = []
    started = time.time()
    for source in selected:
        trait = source['trait_ids']
        final = DEST / source['archive_name']
        partial = DEST / (source['archive_name'] + '.acquiring')
        if final.exists() or partial.exists():
            raise SystemExit('PRIOR_SOURCE_OR_PARTIAL_PRESERVED_NO_OVERWRITE')
        header_path = PACKAGE / 'logs' / (trait + '_download_headers_v2.txt')
        step_started = time.time()
        process = subprocess.run(['curl', '--silent', '--show-error', '--fail', '--location',
                                  '--connect-timeout', '20', '--max-time', '1200',
                                  '--max-filesize', source['archive_bytes'],
                                  '--header', 'Accept-Encoding: identity',
                                  '--dump-header', str(header_path), '--output', str(partial),
                                  source['download_url']], text=True, capture_output=True)
        n = partial.stat().st_size if partial.exists() else 0
        h = sha(partial) if process.returncode == 0 and n == int(source['archive_bytes']) else None
        passed = process.returncode == 0 and n == int(source['archive_bytes']) and h == source['archive_sha256']
        row = {'trait_id': trait, 'source_id': source['source_id'], 'url': source['download_url'],
               'expected_bytes': int(source['archive_bytes']), 'actual_bytes': n,
               'expected_sha256': source['archive_sha256'], 'actual_sha256': h,
               'curl_returncode': process.returncode, 'stderr': process.stderr,
               'elapsed_seconds': time.time() - step_started, 'headers': headers(header_path),
               'header_sha256': sha(header_path), 'status': 'EXACT_HISTORICAL_SOURCE_REACQUIRED' if passed else 'FAILED_PARTIAL_OR_MISMATCH_PRESERVED',
               'path': str(final if passed else partial)}
        if passed:
            partial.rename(final)
        save(PACKAGE / 'logs' / (trait + '_source_acquisition_receipt_v2.json'), row)
        rows.append(row)
        print(json.dumps(row), flush=True)
        if not passed:
            raise SystemExit('FAILED_SOURCE_RECEIPT_WRITTEN; NO_DOWNSTREAM_USE')
    save(PACKAGE / 'logs/exact_lipid_source_acquisition_receipt_v2.json',
         {'completed_utc': now(), 'elapsed_seconds': time.time() - started,
          'sources': rows, 'total_retained_bytes': total, 'all_exact_expected_hashes': True,
          'original_SSD_and_v1_outputs_modified': False, 'native_reproduction_completed': False})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['search-ssd', 'head', 'download'])
    mode = parser.parse_args().mode
    {'search-ssd': search_ssd, 'head': head, 'download': download}[mode]()
