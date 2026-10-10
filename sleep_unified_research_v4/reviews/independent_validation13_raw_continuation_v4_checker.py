"""Narrow independent metadata/AST/arithmetic review; never open GWAS bodies."""
import ast
import hashlib
import json
from datetime import datetime, timedelta
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
load = lambda p: json.loads(Path(p).read_text())
code = P / 'scripts/108_acquire_validation_raw_sources_v4.py'
plan_path = P / 'manifests/validation_raw_acquisition_plan_v4_4.json'
ledger_path = P / 'manifests/global_SSD_resource_reservation_v4_6.json'
assert sha(code) == '6fec56a75ad71d841379436a461ac63d0d8cbcbe0ac62cf4aa960ead6a783d88'
assert sha(plan_path) == '3be0806d08c650bb733fcfb6d96f3f695aad2667f03d10d89de6a80ee8b715ee'
assert sha(ledger_path) == 'd5d52a07477bbd7cac4dd6f0bc6d1bf7f6da2d965564cb360e8b660e299d2c87'
plan = load(plan_path)
old = load(P / 'manifests/validation_raw_acquisition_plan_v4_3.json')
ledger = load(ledger_path)
prior_ledger = load(P / 'manifests/global_SSD_resource_reservation_v4_5.json')
controls = load(R / 'validation13_raw_continuation_author_controls_v4.json')
author_seal_path = R / 'validation13_raw_continuation_author_preparation_seal_v4.json'
assert sha(author_seal_path) == 'c9e8236d259ecdeb61c39aaa6cb497c295425523a50fbe1fd04f37c55c6ea464'
author_seal = load(author_seal_path)
assert len(author_seal['file_sha256']) == 10
for path, digest in author_seal['file_sha256'].items():
    assert sha(path) == digest  # Ten code/metadata files only; no nested seal traversal.
assert plan['member_count'] == len(plan['members']) == len(plan['receipt_origins']) == 13
assert [m['index'] for m in plan['members']] == list(range(1, 14))
assert len({m['source_id'] for m in plan['members']}) == 13
for before, after in zip(old['members'], plan['members']):
    expected = dict(before)
    if after['index'] >= 12:
        expected['body_path'] = str(Path(plan['private_namespace']) / 'raw' / before['filename'])
    assert expected == after
    origin = plan['receipt_origins'][after['source_id']]
    assert origin['reused_v3'] is (after['index'] <= 11)
    epoch = 'v3' if origin['reused_v3'] else 'v4'
    assert str(Path(origin['primary_path']).parent.parent).endswith('validation_raw_replay_' + epoch)
    assert Path(origin['header_path']).name == after['source_id'] + '.headers.txt'
    if origin['reused_v3']:
        assert plan['bound_sources'][origin['primary_path']] == origin['sha256']
        assert plan['bound_sources'][origin['mirror_path']] == origin['sha256']
assert len(plan['bound_sources']) == 120 and all(Path(p).is_absolute() for p in plan['bound_sources'])
assert all(plan['bound_sources'][p] == d for p, d in old['bound_sources'].items())
keys = ['member_count', 'internal_floor_bytes', 'ssd_floor_bytes', 'maximum_owned_transfer_rss_bytes',
        'monitor_output_limit_bytes', 'SSD_reservation_bytes', 'per_body_seconds_limit',
        'family_seconds_limit', 'runtime_poll_seconds', 'exclusive_family_lock_path', 'curl_version']
assert all(plan[k] == old[k] for k in keys)
assert (plan['internal_floor_bytes'], plan['ssd_floor_bytes'], plan['maximum_owned_transfer_rss_bytes'],
        plan['per_body_seconds_limit'], plan['family_seconds_limit'], plan['SSD_reservation_bytes']) == (
        3 << 30, 5 << 30, 2 << 30, 7200, 345600, 300 << 30)
prefix = plan['resume_source12']
assert prefix['prefix_bytes'] == 57671680
assert prefix['prefix_sha256'] == '3537a1b14c7e3560c43d38898dfadd6ba280738498ef5de4c21195c962a1a415'
assert prefix['original_partial_path'] == old['members'][11]['body_path'] + '.partial'
assert sum(m['expected_bytes'] for m in plan['members']) == plan['original_full_body_bytes'] == 10058648185
assert plan['members'][11]['expected_bytes'] - prefix['prefix_bytes'] == 752143553
assert plan['members'][12]['expected_bytes'] == 805600470
assert sum(m['expected_bytes'] for m in plan['members'][11:]) - prefix['prefix_bytes'] == plan['compressed_network_bytes'] == 1557744023
assert datetime.fromisoformat(plan['original_family_deadline_utc']) == (
    datetime.fromisoformat(plan['original_family_started_utc']) + timedelta(seconds=345600))
assert plan['original_family_deadline_utc'] == '2026-10-14T02:58:07.094803+00:00'
# If continuation monotonic time is M and remaining original budget is R,
# start=M-(L-R); the unchanged elapsed>L guard expires exactly R seconds later.
for remaining in [0.25, 3600, 345600]:
    M, L = 1000000, 345600
    start = M - (L - remaining)
    assert (M + remaining) - start == L
expected_components = dict(prior_ledger['component_bytes'])
expected_components['preserved_prefix_bytes'] += prefix['prefix_bytes']
expected_components['validation_raw13_v4_continuation_metadata_cap_bytes'] = 8 << 20
assert ledger['component_bytes'] == expected_components
assert sum(expected_components.values()) == ledger['reserved_total_bytes'] == 319465399037
assert ledger['ceiling_bytes'] == 300 << 30
assert ledger['ceiling_bytes'] - ledger['reserved_total_bytes'] == ledger['unallocated_margin_bytes'] == 2657148163
assert ledger['reserved_total_bytes'] - prior_ledger['reserved_total_bytes'] == 66060288
def functions(path):
    return {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(Path(path).read_text()).body
            if isinstance(n, ast.FunctionDef)}
current = functions(code)
tested = functions(R / 'validation13_raw_continuation_tested_source_v4.py')
assert current.keys() == tested.keys()
assert all(current[n] == tested[n] for n in current if n != 'prepare')
donor = functions(P / 'scripts/90_acquire_validation_raw_sources_v3.py')
assert len(controls['unchanged_helper_AST']) == 14
assert all(current[n] == donor[n] for n in controls['unchanged_helper_AST'])
assert controls['code_sha256'] == sha(R / 'validation13_raw_continuation_tested_source_v4.py')
assert controls['control_count'] == len(controls['controls']) == 18 and controls['assertion_count'] == 83
assert all(c['status'] in ['PASS', 'EXPECTED_REJECTION'] for c in controls['controls'])
print(json.dumps(dict(status='QUALIFIED_PRELAUNCH_DELTA_PASS', plan_sha256=sha(plan_path),
    executor_sha256=sha(code), ledger_sha256=sha(ledger_path), original_members_unchanged=13,
    bound_absolute_metadata_paths=120, old_dependency_hash_values_preserved=True,
    inherited_author_controls=18, inherited_author_assertions=83,
    focused_controls_reexecuted=0, production_body_reads=0, real_workers_network_locks_or_census=0,
    network_bytes=1557744023, reservation_bytes=319465399037, margin_bytes=2657148163,
    original_deadline_utc=plan['original_family_deadline_utc'], root_execution_admission_granted=False)))
