#!/usr/bin/env python3
"""Restore exact DIRECT_TSV source bytes from verified historical gzip copies."""
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
DEST = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/direct_tsv_source_recovery_v2')
TRAITS = {'colorectal_cancer', 'lung_cancer'}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024**2), b''):
            h.update(block)
    return h.hexdigest()

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def save(path, data):
    with path.open('x') as f:
        json.dump(data, f, indent=2)
        f.write('\n')

def main():
    sources = [r for r in read(ROOT/'config/public_gwas_sources.tsv') if r['trait_ids'] in TRAITS]
    checks = read(ROOT/'sleep_unified_research_v1/tables/native_input_hash_checks.tsv')
    if len(sources) != 2 or {r['trait_ids'] for r in sources} != TRAITS or any(r['archive_member'] != 'DIRECT_TSV' for r in sources):
        raise SystemExit('EXACT_DIRECT_TSV_FAMILY_GATE_FAILED')
    if DEST.exists():
        raise SystemExit('PRIOR_RECOVERY_DESTINATION_PRESERVED')
    DEST.mkdir(parents=True)
    total = sum(int(r['archive_bytes']) for r in sources)
    if shutil.disk_usage(DEST).free < total + 5*1024**3:
        raise SystemExit('SSD_RESOURCE_GATE_FAILED')
    inputs = []
    for source in sources:
        matches = [r for r in checks if r['kind'] == 'core_raw' and r['trait_id'] == source['trait_ids']]
        if len(matches) != 1 or matches[0]['status'] != 'MATCH_EXPECTED_SHA256':
            raise SystemExit('HISTORICAL_RAW_IDENTITY_GATE_FAILED')
        inputs.append(dict(source=source, raw_check=matches[0]))
    plan = {'frozen_utc': datetime.now(timezone.utc).isoformat(), 'inputs': inputs,
            'maximum_restored_source_bytes': total, 'network_bytes': 0,
            'stream_buffer_bytes': 4*1024**2, 'workers': 1,
            'destination': str(DEST), 'source_registry_sha256': sha(ROOT/'config/public_gwas_sources.tsv'),
            'original_hash_ledger_sha256': sha(ROOT/'sleep_unified_research_v1/tables/native_input_hash_checks.tsv'),
            'script_sha256': sha(__file__), 'original_SSD_and_v1_modified': False}
    plan_path = PACKAGE/'manifests/direct_tsv_source_recovery_plan_v2.json'
    save(plan_path, plan)
    began = time.time()
    results = []
    for item in inputs:
        source = item['source']
        check = item['raw_check']
        raw = Path(check['path'])
        before = raw.stat()
        if before.st_size != int(check['bytes']) or sha(raw) != check['expected_sha256']:
            raise SystemExit('CURRENT_RAW_HASH_GATE_FAILED')
        expected_n = int(source['archive_bytes'])
        partial = DEST/(source['archive_name']+'.restoring')
        final = DEST/source['archive_name']
        n = 0
        digest = hashlib.sha256()
        with gzip.open(raw, 'rb') as f, partial.open('xb') as out:
            for block in iter(lambda: f.read(4*1024**2), b''):
                n += len(block)
                if n > expected_n:
                    raise SystemExit('RESTORATION_SIZE_CAP_EXCEEDED_PARTIAL_PRESERVED')
                digest.update(block)
                out.write(block)
        after = raw.stat()
        stable = (before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino)
        passed = stable and n == expected_n and digest.hexdigest() == source['archive_sha256']
        if passed:
            partial.rename(final)
        row = {'trait_id': source['trait_ids'], 'source_id': source['source_id'], 'raw_path': str(raw),
               'raw_sha256_before': check['expected_sha256'], 'raw_stable_during_restore': stable,
               'expected_source_bytes': expected_n, 'actual_source_bytes': n,
               'expected_source_sha256': source['archive_sha256'], 'actual_source_sha256': digest.hexdigest(),
               'path': str(final if passed else partial),
               'status': 'EXACT_SOURCE_BYTES_RESTORED_FROM_HISTORICAL_GZIP' if passed else 'MISMATCH_PRESERVED_NOT_ADMITTED',
               'recovery_distinction': 'new decompressed copy; original missing source path remains absent'}
        results.append(row)
        save(PACKAGE/'logs'/(source['trait_ids']+'_direct_source_restoration_v2.json'), row)
        print(json.dumps(row), flush=True)
        if not passed:
            raise SystemExit('EXACT_SOURCE_RESTORATION_FAILED')
    save(PACKAGE/'logs/direct_tsv_source_recovery_receipt_v2.json',
         {'completed_utc': datetime.now(timezone.utc).isoformat(), 'elapsed_seconds': time.time()-began,
          'results': results, 'resource_plan_sha256': sha(plan_path),
          'original_SSD_and_v1_outputs_modified': False, 'native_harmonization_reproduced': False})

if __name__ == '__main__':
    main()
