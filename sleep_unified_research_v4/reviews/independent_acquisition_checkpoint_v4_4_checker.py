#!/usr/bin/env python3
"""Independent stopped-source checkpoint audit; no network or decompression.

Body access is allowlisted to two completed v4 bodies and source3's v4 partial.
Everything else is code, provenance metadata, headers, logs, or stat inventory.
No candidate executor import; only bounded stdlib 64-KiB hash streams.
"""
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import time
from urllib.parse import parse_qs, urlparse

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
REVIEW = P / 'reviews'
PLAN = P / 'manifests/extension_raw_acquisition_plan_v4_4.json'
PLAN_SHA = 'c318a18bd4912dd15ba3ffcdbd8fe9bf5c97aa535591834a2c065238bd23bc9c'
OUT = REVIEW / 'independent_acquisition_checkpoint_receipt_v4_4.json'
PROPOSED = SSD / 'extension_raw_replay_v5'
METADATA = {}
BODY_ALLOWLIST = set()
BODY = {}
START = time.monotonic()
READ_BYTES = 0


def state():
    return {'recorded_utc': datetime.now(timezone.utc).isoformat(),
            'internal_free_bytes': shutil.disk_usage('/System/Volumes/Data').free,
            'ssd_free_bytes': shutil.disk_usage(SSD).free,
            'peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def guard():
    s = state()
    assert s['internal_free_bytes'] >= 3 * 1024**3, 'REVIEW_INTERNAL_3GiB_FLOOR'
    assert s['ssd_free_bytes'] >= 5 * 1024**3, 'REVIEW_SSD_5GiB_FLOOR'
    assert s['peak_RSS_bytes'] <= 128 * 1024**2, 'REVIEW_RSS_CEILING'
    assert time.monotonic() - START < 1800, 'REVIEW_DEADLINE'
    return s


def file_identity(path):
    s = Path(path).stat()
    return {'bytes': s.st_size, 'device': s.st_dev, 'inode': s.st_ino,
            'mtime_ns': s.st_mtime_ns, 'ctime_ns': s.st_ctime_ns}


def stream_hash(path, body=False):
    global READ_BYTES
    path = Path(path)
    if body:
        assert str(path) in BODY_ALLOWLIST, 'BODY_OUTSIDE_EXPLICIT_AUTHORIZATION'
    else:
        assert not path.name.endswith(('.bgz', '.bgz.partial', '.gz', '.gz.partial')), 'BODY_IN_METADATA_HASH_PATH'
    before = file_identity(path)
    md5, sha = hashlib.md5(), hashlib.sha256()
    total = 0
    with path.open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            md5.update(block)
            sha.update(block)
            total += len(block)
            READ_BYTES += len(block)
            if total % (16 * 1024**2) == 0:
                guard()
            if body and total % (1024**3) == 0:
                print(json.dumps({'independent_body_hash': path.name, 'bytes_read': total}), flush=True)
    assert total == before['bytes'] and file_identity(path) == before, 'FILE_CHANGED_DURING_HASH'
    return {'md5': md5.hexdigest(), 'sha256': sha.hexdigest(), **before}


def register(path, expected=None):
    path = Path(path)
    info = stream_hash(path)
    assert not path.name.startswith('._')
    if expected is not None:
        assert info['sha256'] == expected, 'METADATA_HASH_MISMATCH ' + str(path)
    if str(path) in METADATA:
        assert METADATA[str(path)] == info
    METADATA[str(path)] = info
    return info


def read_json(path):
    register(path)
    return json.loads(Path(path).read_text())


def read_tsv(path):
    register(path)
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def headers(path):
    # Parse final HTTP response block independently of executor helper.
    raw = Path(path).read_text()
    blocks = [block for block in re.split(r'\n\s*\n', raw) if block.startswith('HTTP/')]
    assert blocks
    lines = blocks[-1].splitlines()
    result = {'status': int(lines[0].split()[1])}
    for line in lines[1:]:
        if ':' in line:
            key, value = line.split(':', 1)
            key = key.lower().strip()
            assert key not in result, 'DUPLICATE_FINAL_HTTP_HEADER'
            result[key] = value.strip().strip('"')
    return result


def group_snapshot(pgids, owner):
    result = subprocess.run(['ps', '-axo', 'pid=,pgid=,stat='], capture_output=True, text=True, check=True)
    rows = [line.split() for line in result.stdout.splitlines()]
    matches = [{'pid': int(a), 'pgid': int(b), 'state': c} for a, b, c in rows
               if int(b) in pgids or int(a) == owner]
    assert matches == [], 'OLD_ACQUISITION_OWNER_OR_GROUP_STILL_PRESENT'
    return {'checked_utc': datetime.now(timezone.utc).isoformat(), 'owned_pgids': sorted(pgids),
            'family_owner_pid': owner, 'matching_live_processes': matches}


def directory_inventory(folder):
    return [{'path': str(path), 'bytes': path.stat().st_size,
             'AppleDouble_transport_metadata': path.name.startswith('._')}
            for path in sorted(folder.rglob('*')) if path.is_file()]


def main():
    assert not OUT.exists()
    start_state = guard()
    plan = read_json(PLAN)
    assert METADATA[str(PLAN)]['sha256'] == PLAN_SHA
    assert plan['SSD_reservation_bytes'] == 300 * 1024**3
    assert plan['internal_floor_bytes'] == 3 * 1024**3 and plan['ssd_floor_bytes'] == 5 * 1024**3
    assert plan['transfer_worker_count'] == 1 and plan['automatic_retry'] is False
    for path, expected in plan['bound_sources'].items():
        register(path, expected)
    register(SSD / 'manifests' / PLAN.name, PLAN_SHA)
    admission = read_json(P / 'manifests/extension_raw_acquisition_admission_v4_4.json')
    assert admission['execution_admitted'] is True and admission['plan_sha256'] == PLAN_SHA
    assert admission['executor_sha256'] == plan['bound_sources'][plan['executor_path']]
    assert len(admission['independent_review_artifact_sha256']) >= 2
    for path, expected in admission['independent_review_artifact_sha256'].items():
        register(path, expected)
    family_path = P / 'logs/extension_raw_acquisition_family_receipt_v4_4.json'
    family = read_json(family_path)
    folder = Path(plan['members'][0]['body_path']).parent.parent
    register(folder / family_path.name, METADATA[str(family_path)]['sha256'])
    assert family['plan_sha256'] == PLAN_SHA
    assert family['status'] == 'FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT' and family['completed_source_count'] == 2
    assert len(family['source_receipts']) == 2 and 'INTERNAL_FULL_NATIVE_FLOOR_REACHED' in family['stop_reason']
    assert family['resource_after']['internal_free_bytes'] < plan['internal_floor_bytes']
    ownership = read_json(folder / 'family_ownership.json')
    assert ownership['plan_sha256'] == PLAN_SHA and ownership['immutable_attempt'] == 4 and ownership['owner_pid'] == family['owner_pid']
    old_plan = read_json(P / 'manifests/extension_raw_acquisition_plan_v4.json')
    panel = read_tsv(ROOT / 'discovery_extension/config/candidate_traits.tsv')
    remote = read_tsv(ROOT / 'discovery_extension/provenance/panukbb/remote_object_snapshot.tsv')
    assert len(plan['members']) == len(old_plan['members']) == len(panel) == 100
    assert len({m['extension_trait_id'] for m in plan['members']}) == len({m['filename'] for m in plan['members']}) == 100
    remote_rows = [row for row in remote if row['object_role'] == 'phenotype_sumstats']
    assert len(remote_rows) == len({row['extension_trait_id'] for row in remote_rows}) == 100
    remote_index = {row['extension_trait_id']: row for row in remote_rows}
    mixed, historical = [], {}
    for original, member, row in zip(old_plan['members'], plan['members'], panel):
        assert all(member[key] == value for key, value in original.items() if key != 'body_path')
        assert member['extension_trait_id'] == row['extension_trait_id']
        assert member['filename'] == row['source_filename'] and member['expected_bytes'] == int(row['source_file_size_bytes'])
        assert member['expected_md5'] == row['checksum'].split(':')[1]
        assert member['build'] == row['build'] and member['ancestry'] == row['ancestry']
        assert member['phenotype'] == row['phenotype_name'] and member['source_release'] == row['study_accession']
        snap = remote_index[member['extension_trait_id']]
        assert member['url'] == snap['versioned_url'] and member['expected_s3_version_id'] == snap['s3_version_id']
        assert member['expected_etag'] == snap['observed_etag'] and member['expected_bytes'] == int(snap['expected_size_bytes']) == int(snap['observed_content_length'])
        assert member['expected_md5'] == snap['expected_checksum'].split(':')[1]
        parsed = urlparse(member['url'])
        assert parsed.scheme == 'https' and parse_qs(parsed.query) == {'versionId': [member['expected_s3_version_id']]}
        h = read_json(member['historical_streaming_receipt'])
        source = h['source_verification']
        assert source['verification_status'] == 'PASS'
        assert source['source'] == h['source_url'] == member['url']
        assert source['expected_md5'] == source['observed_md5'] == member['expected_md5']
        assert source['expected_size_bytes'] == source['observed_size_bytes'] == member['expected_bytes']
        assert source['observed_sha256'] == member['expected_sha256']
        assert source['http_version_id'] == member['expected_s3_version_id'] and source['http_etag'].strip('"') == member['expected_etag']
        historical[member['extension_trait_id']] = source
        reuse = member['index'] <= 2
        body_path = Path(member['body_path']) if reuse else PROPOSED / 'raw' / member['filename']
        receipt_path = (folder if reuse else PROPOSED) / 'receipts' / (member['extension_trait_id'] + '.json')
        if not reuse:
            assert not body_path.exists() and not Path(str(body_path) + '.partial').exists() and not receipt_path.exists()
        mixed.append({**member, 'legacy_v4_body_path': member['body_path'], 'body_path': str(body_path),
                      'receipt_path': str(receipt_path), 'origin_attempt': 'v4_4_REUSE_UNCHANGED' if reuse else 'v5_NEW_EXPLICIT_ATTEMPT',
                      'body_acquired_now': reuse, 'reused_receipt_sha256': family['source_receipts'][member['index'] - 1]['sha256'] if reuse else None,
                      'continuation_method': 'UNCHANGED_EXACT_BODY_AND_RECEIPT' if reuse else ('COPY_PRESERVED_SOURCE3_PREFIX_THEN_EXACT_HTTP206_RANGE' if member['index'] == 3 else 'PINNED_FULL_HTTP200_BODY')})
    assert sum(m['expected_bytes'] for m in plan['members']) == plan['compressed_network_bytes']
    receipts, pgids = [], set()
    for member in plan['members'][:3]:
        trait = member['extension_trait_id']
        rp = folder / 'receipts' / (trait + '.json')
        r = read_json(rp)
        register(P / 'source_provenance/extension_raw_acquisition_v4_4' / rp.name, METADATA[str(rp)]['sha256'])
        assert r['member'] == member and r['plan_sha256'] == PLAN_SHA and r['original_source_sha256_required'] == member['expected_sha256']
        assert r['teardown']['remaining_group_members'] == [] and r['teardown']['teardown_verified'] is True
        worker_ids = {pid for sample in r['samples'] for pid in sample['owned_pids']}
        assert len(worker_ids) == 1
        pgids.update(worker_ids)
        header = folder / 'logs' / (trait + '.headers.txt')
        curl_stdout = folder / 'logs' / (trait + '.curl.log')
        register(header, r['headers_sha256'])
        register(curl_stdout)
        assert header.stat().st_size == r['headers_bytes']
        parsed_header = headers(header)
        offset = plan['explicit_resume']['prefix_bytes'] if member['index'] == 1 else 0
        expected_command = ['/usr/bin/curl', '-q', '--fail', '--location', '--max-redirs', '5', '--proto', '=https',
                            '--proto-redir', '=https', '--tlsv1.2', '--max-time', '7200', '--max-filesize', str(member['expected_bytes']),
                            '--dump-header', str(header), '--output', member['body_path'] + '.partial']
        if offset:
            expected_command += ['--continue-at', str(offset)]
        expected_command += [member['url']]
        assert r['command'] == expected_command and r['resume_offset'] == offset
        assert parsed_header['status'] == (206 if offset else 200)
        assert parsed_header['x-amz-version-id'] == member['expected_s3_version_id'] and parsed_header['etag'] == member['expected_etag']
        assert int(parsed_header['content-length']) == member['expected_bytes'] - offset
        assert parsed_header['last-modified'] == historical[trait]['http_last_modified']
        if offset:
            assert parsed_header['content-range'] == 'bytes %d-%d/%d' % (offset, member['expected_bytes'] - 1, member['expected_bytes'])
        if member['index'] <= 2:
            assert family['source_receipts'][member['index'] - 1] == {'path': str(rp), 'sha256': METADATA[str(rp)]['sha256']}
            assert r['status'] == 'EXACT_IMMUTABLE_SOURCE_ACQUIRED' and r['returncode'] == 0 and r['stop_reason'] is None
            assert r['post_cleanup_hash_resource_identity_gates_pass'] is True and r['observed_headers'] == parsed_header
            assert r['actual_size'] == r['retained_partial_bytes'] == member['expected_bytes']
            assert r['actual_md5'] == r['final_seal_md5'] == member['expected_md5']
            assert r['actual_sha256'] == r['final_seal_sha256'] == r['retained_partial_sha256'] == member['expected_sha256']
            assert r['gzip_full_stream_verified'] is False and r['source_body_committed'] is False
            body_path = Path(member['body_path'])
            assert body_path.is_file() and not Path(str(body_path) + '.partial').exists()
        else:
            assert r['status'] == 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and r['returncode'] == -15
            assert r['stop_reason'] == 'INTERNAL_FULL_NATIVE_FLOOR_REACHED'
            assert r['resource_after']['internal_free_bytes'] < plan['internal_floor_bytes']
            body_path = Path(member['body_path'] + '.partial')
            assert not Path(member['body_path']).exists() and body_path.stat().st_size == r['retained_partial_bytes'] == 1276116992
            assert 'actual_sha256' not in r and 'actual_size' not in r
        BODY_ALLOWLIST.add(str(body_path))
        receipts.append({'index': member['index'], 'trait': trait, 'receipt_path': str(rp),
                         'receipt_sha256': METADATA[str(rp)]['sha256'], 'status': r['status'],
                         'headers_path': str(header), 'headers_sha256': r['headers_sha256'],
                         'body_path': str(body_path), 'expected_current_sha256': r['retained_partial_sha256'],
                         'expected_current_bytes': r['retained_partial_bytes'], 'teardown': r['teardown'],
                         'returncode': r['returncode'], 'stop_reason': r['stop_reason'],
                         'parsed_final_HTTP_headers': parsed_header})
    group_before = group_snapshot(pgids, family['owner_pid'])
    with Path(plan['exclusive_family_lock_path']).open('rb') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            for pass_number in [1, 2]:
                for row in receipts:
                    info = stream_hash(row['body_path'], body=True)
                    assert info['sha256'] == row['expected_current_sha256'] and info['bytes'] == row['expected_current_bytes']
                    if row['index'] <= 2:
                        assert info['md5'] == plan['members'][row['index'] - 1]['expected_md5']
                    if pass_number == 1:
                        BODY[row['body_path']] = info
                    else:
                        assert BODY[row['body_path']] == info, 'AUTHORIZED_BODY_CHANGED_BETWEEN_HASH_PASSES'
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    old_prefix = Path(plan['explicit_resume']['original_partial_path'])
    old_prefix_before = file_identity(old_prefix)
    assert old_prefix_before['bytes'] == plan['explicit_resume']['prefix_bytes'] == 1129316352
    inventories = {str(folder): directory_inventory(folder) for folder in sorted(SSD.glob('extension_raw_replay*')) if folder.is_dir()}
    assert sum(not Path(row['path']).name.startswith('._') for row in inventories[str(folder)] if Path(row['path']).parent.name == 'raw') == 3
    # Budget includes first2 reused bodies exactly once, the original source1
    # failed prefix and the v4 source3 failed prefix retained after a v5 copy.
    preserved_prefix_bytes = old_prefix_before['bytes'] + receipts[2]['expected_current_bytes']
    existing_metadata_bytes = sum(row['bytes'] for entries in inventories.values() for row in entries
                                  if row['path'] not in BODY_ALLOWLIST and row['path'] != str(old_prefix))
    final_body_storage_bytes = plan['compressed_network_bytes'] + preserved_prefix_bytes
    minimum_storage_with_current_metadata = final_body_storage_bytes + existing_metadata_bytes
    assert minimum_storage_with_current_metadata + plan['ssd_floor_bytes'] < plan['SSD_reservation_bytes']
    body_stat_after = {path: file_identity(path) == {key: info[key] for key in ['bytes', 'device', 'inode', 'mtime_ns', 'ctime_ns']} for path, info in BODY.items()}
    metadata_after = {path: stream_hash(path) == info for path, info in METADATA.items()}
    assert all(body_stat_after.values()) and all(metadata_after.values()) and file_identity(old_prefix) == old_prefix_before
    group_after = group_snapshot(pgids, family['owner_pid'])
    final_state = guard()
    result = {
        'status': 'PASS_TWO_EXACT_BODIES_AND_SOURCE3_PRESERVED_PARTIAL',
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'review_scope': 'Stopped acquisition checkpoint and proposed mixed-source continuation only; no execution admission.',
        'plan_sha256': PLAN_SHA, 'family_receipt_sha256': METADATA[str(family_path)]['sha256'],
        'family_failed_status_preserved': True, 'body_successes': 2, 'failed_partial_sources': 1, 'unstarted_sources': 97,
        'source_receipt_proofs': receipts, 'authorized_body_hashes': BODY, 'body_hash_passes': 2,
        'authorized_body_bytes_per_pass': sum(info['bytes'] for info in BODY.values()),
        'before_after_body_hashes_and_stats_match': True, 'before_after_metadata_hashes_match': True,
        'current_source3_partial_full_source_identity_admitted': False,
        'source3_resume': {'source_trait': plan['members'][2]['extension_trait_id'], 'original_partial_path': receipts[2]['body_path'],
                           'prefix_bytes': receipts[2]['expected_current_bytes'], 'prefix_sha256': receipts[2]['expected_current_sha256'],
                           'prefix_md5_diagnostic_only': BODY[receipts[2]['body_path']]['md5'],
                           'new_partial_path': mixed[2]['body_path'] + '.partial', 'new_final_path': mixed[2]['body_path'],
                           'remaining_network_bytes': plan['members'][2]['expected_bytes'] - receipts[2]['expected_current_bytes'],
                           'required_HTTP_status': 206,
                           'required_Content_Range': 'bytes %d-%d/%d' % (receipts[2]['expected_current_bytes'], plan['members'][2]['expected_bytes'] - 1, plan['members'][2]['expected_bytes']),
                           'required_versionId': plan['members'][2]['expected_s3_version_id'],
                           'required_ETag': plan['members'][2]['expected_etag'], 'required_whole_MD5': plan['members'][2]['expected_md5'],
                           'required_whole_SHA256': plan['members'][2]['expected_sha256']},
        'original_source1_prefix': {'path': str(old_prefix), 'bytes': old_prefix_before['bytes'],
                                    'frozen_SHA256': plan['explicit_resume']['prefix_sha256'], 'body_rehashed_here': False,
                                    'scope': 'Current stat identity and unchanged frozen plan/failure receipt; excluded from explicit body-read authorization.',
                                    'stat_before_after_match': True},
        'proposed_mixed_origin_members': mixed,
        'mixed_origin_cardinality': {'ordered_members': 100, 'distinct_traits': 100, 'distinct_filenames': 100,
                                    'reused_v4_successes': 2, 'new_v5_sources': 98},
        'budget': {'global_reservation_bytes': plan['SSD_reservation_bytes'],
                   'complete100_body_bytes_with_first2_counted_once': plan['compressed_network_bytes'],
                   'source1_original_preserved_prefix_bytes': old_prefix_before['bytes'],
                   'source3_v4_preserved_prefix_bytes': receipts[2]['expected_current_bytes'],
                   'all_preserved_extra_prefix_bytes': preserved_prefix_bytes,
                   'final_body_and_preserved_prefix_bytes': final_body_storage_bytes,
                   'existing_acquisition_metadata_and_transport_bytes': existing_metadata_bytes,
                   'minimum_final_storage_including_existing_metadata': minimum_storage_with_current_metadata,
                   'remaining_within300GiB_before_SSD_floor_and_new_metadata': plan['SSD_reservation_bytes'] - minimum_storage_with_current_metadata,
                   'remaining_network_bytes_after_two_reused_bodies_and_source3_resume': sum(m['expected_bytes'] for m in mixed[2:]) - receipts[2]['expected_current_bytes'],
                   'limit_scope': 'Global acquisition namespaces, every old/new final/partial body, transport metadata, logs/receipts and working copies count. Future scratch/new metadata require frozen allocation plus runtime total-byte gates; this is a minimum projected footprint, not a complete runtime reservation proof.'},
        'current_acquisition_namespace_stat_inventory': inventories,
        'old_owned_group_proofs': {'before': group_before, 'after': group_after, 'exclusive_family_flock_acquired_and_released_read_only': True},
        'metadata_bindings': METADATA,
        'source_identity_status': 'All100 unchanged original panel/snapshot/source URL/version/ETag/bytes/MD5/SHA provenance reconciled; actual complete bytes independently verified only for first2.',
        'continuation_requirements': [
            'Freeze a newv5 plan/code/monitor/admission/ownership namespace, binding this review, all old plans/receipts/failures/ownership, current successful bodies/receipts and retained prefixes; do not mutate v4.',
            'Use exact mixed100 identity and origin paths: reuse first2 v4 body/receipt/header identities unchanged, acquire remaining98 only at newv5 paths; no success rewriting or re-download.',
            'Copy source3 current prefix into an exclusive newv5 partial, prove original and copy byte/SHA identity, fsync, preserve original; then request exact pinned HTTPS version206 suffix and Content-Range/Content-Length/ETag/versionId.',
            'Final source3 admission requires full original bytecount+MD5+SHA; prefix MD5 is diagnostic and cannot substitute. No gzip/EOF/scientific chain claim follows from hash equality.',
            'All subsequent97 sources require exact pinned200 headers/bytes/originalMD5/SHA and immutable success/failure receipts, one worker, no automatic retries or source substitutions.',
            'Enforce global300GiB accounting including both extra retained prefixes and all old/new namespaces/copies/metadata; retain unchanged3GiB internal and5GiB SSD floors, original RSS/time/max-filesize limits and terminal/posthash gates.',
            'Hold the exclusive family flock through immutable receipts and authoritative owned-group teardown/quarantine with deferred catchable signals; fresh root admission must bind nonempty exact independent MD/JSON reviews.',
            'Reconcile final100 membership/order/identity and both old reused proofs plus98 new proofs. Preserve failedv4 family receipt and write distinctv5 all100 receipt; no prior failure relabel or partial-family promotion.'
        ],
        'network_requests': 0, 'body_writes': 0, 'GWAS_decompressions': 0, 'native_workers_or_estimators_launched': 0,
        'scientific_reproduction_or_results_admitted': False, 'execution_authorized': False,
        'resource_plan': {'hash_buffer_bytes': 65536, 'internal_floor_bytes': 3 * 1024**3, 'ssd_floor_bytes': 5 * 1024**3,
                          'maximum_RSS_bytes': 128 * 1024**2, 'deadline_seconds': 1800, 'body_workers': 0,
                          'initial': start_state, 'final': final_state, 'elapsed_seconds': time.monotonic() - START,
                          'read_bytes_including_metadata_and_two_body_passes': READ_BYTES},
        'checker_sha256': stream_hash(Path(__file__))['sha256']
    }
    with OUT.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({key: result[key] for key in ['status', 'body_successes', 'failed_partial_sources', 'unstarted_sources', 'authorized_body_bytes_per_pass', 'mixed_origin_cardinality', 'budget', 'resource_plan']}, indent=2))


if __name__ == '__main__':
    main()
