#!/usr/bin/env python3
"""Independent small-buffer metadata audit; no candidate import or network IO."""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
PLAN = ROOT / 'manifests/extension_raw_acquisition_plan_v4_4.json'
EXPECT_PLAN = 'c318a18bd4912dd15ba3ffcdbd8fe9bf5c97aa535591834a2c065238bd23bc9c'
OUT = ROOT / 'reviews/independent_acquisition_v4_binding_receipt.json'

def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while True:
            chunk = handle.read(65536)
            if not chunk:
                break
            result.update(chunk)
    return result.hexdigest()

def read(path):
    return json.loads(Path(path).read_text())

def check(value, message):
    if not value:
        raise AssertionError(message)

start = time.monotonic()
check(digest(PLAN) == EXPECT_PLAN, 'plan identity')
plan = read(PLAN)
old_path = ROOT / 'manifests/extension_raw_acquisition_plan_v4.json'
old = read(old_path)
check(len(plan['members']) == len(old['members']) == 100, '100 members')
check(len({r['extension_trait_id'] for r in plan['members']}) == 100, 'unique trait')
check(len({r['filename'] for r in plan['members']}) == 100, 'unique filename')
identities = []
for old_row, new_row in zip(old['members'], plan['members']):
    for key, value in old_row.items():
        if key != 'body_path':
            check(new_row[key] == value, 'historical member changed: ' + key)
    receipt_path = Path(new_row['historical_streaming_receipt'])
    expected_path = REPO / 'discovery_extension/provenance/streaming_receipts' / (new_row['extension_trait_id'] + '.json')
    check(receipt_path == expected_path, 'original receipt path')
    original = read(receipt_path)
    source = original['source_verification']
    check(source['verification_status'] == 'PASS', 'original receipt status')
    check(source['source'] == original['source_url'] == new_row['url'], 'url')
    check(source['observed_sha256'] == new_row['expected_sha256'], 'original whole-body sha')
    check(re.fullmatch('[0-9a-f]{64}', new_row['expected_sha256']) is not None, 'sha grammar')
    check(source['expected_md5'] == source['observed_md5'] == new_row['expected_md5'], 'md5')
    check(source['expected_size_bytes'] == source['observed_size_bytes'] == new_row['expected_bytes'], 'size')
    check(source['http_version_id'] == new_row['expected_s3_version_id'], 'S3 version')
    check(source['http_etag'].strip('"') == new_row['expected_etag'], 'etag')
    check(receipt_path.as_posix() in plan['bound_sources'], 'original receipt not bound')
    check(Path(new_row['body_path']).parent.name == 'raw' and Path(new_row['body_path']).parents[1].name == 'extension_raw_replay_v4', 'new attempt namespace')
    identities.append({'index': new_row['index'], 'trait': new_row['extension_trait_id'], 'original_sha256': new_row['expected_sha256']})
for key, value in old.items():
    if key not in ('members', 'bound_sources', 'prepared_utc', 'resource_at_plan'):
        check(plan[key] == value, 'original policy changed: ' + key)
check(sum(row['expected_bytes'] for row in plan['members']) == plan['compressed_network_bytes'], 'sum bytes')
check(max(row['expected_bytes'] for row in plan['members']) == plan['largest_body_bytes'], 'largest bytes')
check(plan['transfer_worker_count'] == 1 and plan['automatic_retry'] is False, 'one worker/no retry')
check(plan['SSD_reservation_bytes'] > plan['compressed_network_bytes'] + plan['additional_preserved_prefix_bytes'] + plan['ssd_floor_bytes'], 'reservation budget')
bound = {}
for name, expected in plan['bound_sources'].items():
    observed = digest(name)
    check(observed == expected, 'bound hash mismatch: ' + name)
    bound[name] = {'sha256': observed, 'bytes': Path(name).stat().st_size}
executor = Path(plan['executor_path'])
control_source = ROOT / 'scripts/53_verify_acquisition_fault_controls_v4.py'
control_path = ROOT / 'logs/extension_acquisition_fault_controls_v4_4.json'
controls = read(control_path)
check(controls['executor_sha256'] == digest(executor), 'control executor identity')
check(controls['plan_sha256'] == EXPECT_PLAN, 'control plan identity')
check(controls['test_script_sha256'] == digest(control_source), 'control code identity')
check(controls['all_pass'] is True and len(controls['checks']) == 8 and all(c['pass'] is True for c in controls['checks']), 'eight controls')
check(controls['network_requests'] == controls['native_estimators_launched'] == 0, 'control scope')
for path in (executor, control_source):
    ast.parse(path.read_text(), filename=str(path))
ssd_plan = Path(plan['exclusive_family_lock_path']).parent / 'manifests' / PLAN.name
check(digest(ssd_plan) == EXPECT_PLAN, 'SSD plan copy')
partial = Path(plan['explicit_resume']['original_partial_path'])
check(partial.stat().st_size == plan['explicit_resume']['prefix_bytes'], 'preserved prefix stat')
check(not (ROOT / 'manifests/extension_raw_acquisition_admission_v4_4.json').exists(), 'unexpected execution admission')
for name, record in bound.items():
    check(digest(name) == record['sha256'], 'bound file changed during audit: ' + name)
check(digest(PLAN) == EXPECT_PLAN, 'plan changed during audit')
receipt = {
    'status': 'PASS_BOUND_METADATA_AND_CONTROL_IDENTITY_NO_NETWORK_NO_NATIVE',
    'completed_utc': datetime.now(timezone.utc).isoformat(),
    'independent_checker_sha256': digest(__file__),
    'plan_sha256': EXPECT_PLAN,
    'executor_sha256': digest(executor),
    'control_receipt_sha256': digest(control_path),
    'control_script_sha256': digest(control_source),
    'bound_dependencies_count': len(bound),
    'bound_dependencies_bytes_per_pass': sum(r['bytes'] for r in bound.values()),
    'before_after_dependency_hashes_pass': True,
    'bound_dependencies': bound,
    'ordered_original_100_identity_pass': True,
    'original_members': identities,
    'compressed_network_bytes': plan['compressed_network_bytes'],
    'preserved_prefix_bytes': partial.stat().st_size,
    'prefix_body_rehashed_by_this_review': False,
    'prior_prefix_hash_scope': 'Frozen plan and bound failure proof; existing 1.129GB prefix stat only in this review.',
    'observed_controls': controls['checks'],
    'control_coverage_limit': 'Read and hash verify previously run controls; did not rerun. Unverified-quarantine branch controls use a completed false worker, while a separate actual live child exercises SIGKILL/reap fallback.',
    'admission_present_at_review': False,
    'source_bodies_acquired_by_this_review': 0,
    'estimators_launched_by_this_review': 0,
    'gzip_and_scientific_reproduction_admitted': False,
    'elapsed_seconds': time.monotonic() - start,
    'peak_rss_platform_units': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
}
with OUT.open('x') as handle:
    json.dump(receipt, handle, indent=2)
    handle.write('\n')
print(json.dumps({k: receipt[k] for k in ('status', 'plan_sha256', 'bound_dependencies_count', 'bound_dependencies_bytes_per_pass', 'elapsed_seconds')}))
