"""Consume the genuine13-source v4 producer, preserving11 v3 origins and206."""
import importlib.util
import json
from pathlib import Path

from extension_replay_common_v4 import sha
from terminal_commit_common_v2 import require_committed


def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(v.is_symlink() for v in p.parents):
        raise RuntimeError('REGULAR_ACQUISITION_PROFILE_EVIDENCE_REQUIRED: '+str(p))
    return p


def fixed_hashes(mapping):
    if not mapping:raise RuntimeError('NONEMPTY_ACQUISITION_PROFILE_MAP_REQUIRED')
    for path,digest in mapping.items():
        if sha(regular(path))!=digest:raise RuntimeError('ACQUISITION_PROFILE_CHANGED: '+path)


def completed_acquisition(plan):
    contract=plan['acquisition_terminal'];fixed_hashes(contract['metadata_sha256'])
    path=regular(contract['producer_plan_path']);raw=json.loads(path.read_text())
    if sha(path)!=contract['binding']['plan_sha256'] or raw['receipt_origins']!=contract['receipt_origins']:
        raise RuntimeError('EXACT_CURRENT_V4_ORIGIN_PLAN_REQUIRED')
    if (raw['reused_source_count'],raw['new_transfer_count'],len(raw['members']),raw['original_full_body_bytes'])!=(11,2,13,10058648185):
        raise RuntimeError('EXACT11_REUSED2_NEW_ORIGINAL13_FAMILY_REQUIRED')
    if [m['index'] for m in raw['members']]!=list(range(1,14)) or raw['resume_source12']['prefix_bytes']!=57671680:
        raise RuntimeError('EXACT_SOURCE_ORDER_AND_SOURCE12_RESUME_REQUIRED')
    pairs=contract['source_receipt_pairs']
    if len(pairs)!=13 or [x['member'] for x in pairs]!=raw['members'] or len({x['member']['source_id'] for x in pairs})!=13:
        raise RuntimeError('EXACT13_ORDERED_DISTINCT_PRODUCER_MEMBERS_REQUIRED')
    paths=[x[k] for x in pairs for k in ['primary_path','mirror_path']]
    if len(paths)!=26 or len(set(paths))!=26 or paths!=contract['individual_receipt_paths']:
        raise RuntimeError('EXACT26_DISTINCT_FROZEN_RECEIPT_COPIES_REQUIRED')
    executor=regular(contract['producer_executor_path'])
    if str(executor)!=raw['executor_path'] or sha(executor)!=contract['binding']['executor_sha256']:
        raise RuntimeError('EXACT108_PRODUCER_IMPLEMENTATION_REQUIRED')
    spec=importlib.util.spec_from_file_location('_reviewed_validation108_profile',executor)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    expected_sources=[]
    for pair in pairs:
        origin=raw['receipt_origins'][pair['member']['source_id']]
        if (pair['primary_path'],pair['mirror_path'])!=(origin['primary_path'],origin['mirror_path']):
            raise RuntimeError('EXACT_ORIGIN_RECEIPT_ROUTES_REQUIRED')
        for value in [pair['primary_path'],pair['mirror_path']]:
            if contract['metadata_sha256'].get(value)!=pair['sha256'] or not Path(value).is_absolute():
                raise RuntimeError('EXACT_FROZEN_DOUBLE_RECEIPT_BINDING_REQUIRED')
        item=dict(path=pair['primary_path'],sha256=pair['sha256'])
        # Existing108 checks origin-specific plan, command, status, finalSHA/MD5,
        # source12 exact206 range, headers and cleanup. No body hashing here;
        # unchanged96/97 current_source and local EOF paths perform that work.
        producer.source_receipt_gate(raw,pair['member'],item,contract['binding']['plan_sha256'],verify_body=False)
        expected_sources.append(item)
    primary=contract['primary_receipt_sha256']
    if len(primary)!=2 or len(set(primary.values()))!=1:
        raise RuntimeError('TWO_BYTE_IDENTICAL_GENUINE_FULL_FAMILY_COPIES_REQUIRED')
    for path,digest in primary.items():
        if contract['metadata_sha256'].get(path)!=digest:raise RuntimeError('EXACT_FULL_FAMILY_PRIMARY_MAP_REQUIRED')
        family=json.loads(regular(path).read_text())
        if (family.get('status')!='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED'
            or family.get('completed_source_count')!=13 or family.get('plan_sha256')!=contract['binding']['plan_sha256']
            or family.get('source_receipts')!=expected_sources or family.get('receipt_origins')!=raw['receipt_origins']
            or family.get('full_current_body_SHA_MD5_size_verified') is not True
            or family.get('reused_v3_source_count')!=11 or family.get('newly_acquired_source_count')!=2
            or family.get('old_failed_family_never_relabelled') is not True):
            raise RuntimeError('EXACT_GENUINE_TERMINAL13_SOURCE_FAMILY_REQUIRED')
    require_committed(contract['pending_path'],contract['seal_path'],contract['binding'],primary)
    fixed_hashes(contract['metadata_sha256'])
    return dict(status='GENUINE13_TERMINAL_MIXED_ORIGIN_PRODUCER_METADATA_PASS',sources=13,
                inherited_v3_sources=11,new_v4_sources=2,source12_resume_offset=57671680,
                fresh_source_body_hashes_performed=False,gzip_EOF_or_scientific_replay_claimed=False)
