#!/usr/bin/env python3
"""Bind the exact locked100 immutable public bodies and bounded storage budget."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    panel_path = ROOT / 'discovery_extension/config/candidate_traits.tsv'
    snapshot_path = ROOT / 'discovery_extension/provenance/panukbb/remote_object_snapshot.tsv'
    panel, snapshot = rows(panel_path), rows(snapshot_path)
    remote = {r['extension_trait_id']: r for r in snapshot if r['object_role'] == 'phenotype_sumstats'}
    assert len(panel) == len(remote) == 100
    out = PACKAGE / 'manifests/extension_raw_acquisition_plan_v4.json'
    assert not out.exists()
    members = []
    for index, r in enumerate(panel, 1):
        q = remote[r['extension_trait_id']]
        assert int(r['source_file_size_bytes']) == int(q['expected_size_bytes'])
        assert r['checksum'] == q['expected_checksum'] and r['checksum'].startswith('md5:')
        assert r['source_filename'] == q['filename']
        assert 'versionId=' in q['versioned_url'] and q['s3_version_id'] in q['versioned_url']
        assert q['head_verification_status'].startswith('HEAD_SIZE_VERSION')
        members.append({'index': index, 'extension_trait_id': r['extension_trait_id'],
                        'url': q['versioned_url'], 'filename': q['filename'],
                        'expected_bytes': int(q['expected_size_bytes']), 'expected_md5': r['checksum'][4:],
                        'expected_s3_version_id': q['s3_version_id'], 'expected_etag': q['observed_etag'],
                        'build': r['build'], 'ancestry': r['ancestry'], 'phenotype': r['phenotype_name'],
                        'source_release': r['study_accession'], 'body_path': str(SSD / 'extension_raw_replay/raw' / q['filename'])})
    expected = sum(r['expected_bytes'] for r in members)
    assert expected == 227610388647
    reserve = 300 * 1024**3
    internal, free = shutil.disk_usage('/System/Volumes/Data').free, shutil.disk_usage(SSD).free
    assert internal >= 3 * 1024**3 and free >= reserve
    source_paths = [panel_path, snapshot_path, PACKAGE / 'FROZEN_EXTENSION_RAW_REPLAY_RESOURCE_PLAN_v1.md', Path(__file__)]
    plan = {'prepared_utc': datetime.now(timezone.utc).isoformat(), 'members': members,
            'compressed_network_bytes': expected, 'SSD_reservation_bytes': reserve, 'largest_body_bytes': max(r['expected_bytes'] for r in members),
            'internal_floor_bytes': 3 * 1024**3, 'ssd_floor_bytes': 5 * 1024**3,
            'maximum_owned_transfer_rss_bytes': 2 * 1024**3, 'per_body_seconds_limit': 7200,
            'family_seconds_limit': 96 * 3600, 'transfer_worker_count': 1, 'automatic_retry': False,
            'retained_raw_bodies': True, 'protected_GWAS_bodies_committed': False,
            'bound_sources': {str(p): sha(p) for p in source_paths},
            'resource_at_plan': {'internal_free_bytes': internal, 'ssd_free_bytes': free},
            'scientific_membership_or_threshold_changes': False, 'scope': 'EXACT_HISTORICAL_RAW_PROVENANCE_REPLAY_ONLY'}
    out.write_text(json.dumps(plan, indent=2) + '\n')
    (SSD / 'manifests' / out.name).write_bytes(out.read_bytes())
    print(json.dumps({'members': len(members), 'network_GiB': expected / 1024**3, 'SSD_budget_GiB': 300, 'plan_sha256': sha(out)}))


if __name__ == '__main__':
    main()
