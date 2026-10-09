#!/usr/bin/env python3
"""Copy and verify a concrete relocation candidate; never remove original files."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SOURCE = Path('/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03/data/raw.local-preserved/.archives')
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
TARGET = SSD / 'storage_relocation_candidate' / 'archives'

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def guard():
    if shutil.disk_usage('/System/Volumes/Data').free < 128 * 1024**2:
        raise RuntimeError('INTERNAL_COPY_EMERGENCY_FLOOR')
    if shutil.disk_usage(SSD).free < 5 * 1024**3:
        raise RuntimeError('SSD_COPY_RESERVE')

def main():
    start = time.monotonic()
    registry = ROOT / 'config/public_gwas_sources.tsv'
    with registry.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    expected = {}
    for row in rows:
        if (SOURCE / row['archive_name']).is_file():
            expected[row['archive_name']] = row['archive_sha256']
    companion = 'neale_ukbb_round2_variants.tsv.bgz'
    companion_sha = 'e7f035e9264536b416e7b7a514dd863bd508dc380f40e48f12ffbc7c06891ea5'
    assert companion_sha in registry.read_text(), 'COMPANION_NOT_BOUND_TO_REGISTRY'
    expected[companion] = companion_sha
    actual_names = {p.name for p in SOURCE.iterdir() if p.is_file()}
    assert set(expected) == actual_names and len(expected) == 7
    total = sum((SOURCE / name).stat().st_size for name in expected)
    assert total == 4355618590
    if shutil.disk_usage(SSD).free < total + 5 * 1024**3:
        raise RuntimeError('COPY_CAPACITY_GATE')
    plan = {
        'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'source_directory': str(SOURCE), 'target_directory': str(TARGET),
        'files': expected, 'total_bytes': total,
        'registry_sha256': sha(registry), 'script_sha256': sha(Path(__file__)),
        'worker_count': 1, 'maximum_buffer_bytes': 4 * 1024**2,
        'network_bytes': 0, 'ssd_reserve_bytes': 5 * 1024**3,
        'internal_emergency_floor_bytes': 128 * 1024**2,
        'original_changes_authorized_by_this_script': False,
        'relocation_requires_explicit_permission': 'Preserved files live in a checkout used by another research project; no deletion or symlink replacement here.',
        'frozen_full_native_internal_floor_bytes': 3 * 1024**3,
        'full_native_guard_relaxed': False,
    }
    plan_path = PACKAGE / 'manifests/internal_archive_relocation_candidate_plan_v4.json'
    if plan_path.exists():
        raise RuntimeError('PLAN_ALREADY_EXISTS_NO_SILENT_RETRY')
    plan_path.write_text(json.dumps(plan, indent=2) + '\n')
    TARGET.mkdir(parents=True, exist_ok=False)
    proof = []
    for name in sorted(expected):
        guard()
        if time.monotonic() - start > 1800:
            raise RuntimeError('COPY_DEADLINE')
        source, target = SOURCE / name, TARGET / name
        stat_before = source.stat()
        h = hashlib.sha256()
        with source.open('rb') as inp, target.open('xb') as out:
            for b in iter(lambda: inp.read(4 * 1024**2), b''):
                guard()
                out.write(b)
                h.update(b)
            out.flush()
            os.fsync(out.fileno())
        shutil.copystat(source, target)
        copied = h.hexdigest()
        target_sha = sha(target)
        source_after_sha = sha(source)
        stat_after = source.stat()
        assert copied == target_sha == source_after_sha == expected[name], name
        assert (stat_before.st_size, stat_before.st_mtime_ns) == (stat_after.st_size, stat_after.st_mtime_ns)
        proof.append({'name': name, 'source_path': str(source), 'target_path': str(target),
                      'bytes': stat_before.st_size, 'expected_sha256': expected[name],
                      'copy_stream_sha256': copied, 'target_sha256': target_sha,
                      'source_after_sha256': source_after_sha, 'all_match': True})
        print(json.dumps({'verified': name, 'bytes': stat_before.st_size}), flush=True)
    receipt = {'completed_utc': datetime.now(timezone.utc).isoformat(),
               'elapsed_seconds': time.monotonic() - start, 'plan_sha256': sha(plan_path),
               'files': proof, 'files_verified': len(proof), 'copied_bytes': total,
               'original_files_removed': 0, 'original_directory_replaced': False,
               'internal_bytes_freed': 0, 'native_jobs_launched': 0,
               'proposed_atomic_action': 'Only after explicit authorization, relocate the verified seven-file directory to the SSD and leave a symlink at its original path; preserve all Git worktrees and original bytes.',
               'internal_free_bytes': shutil.disk_usage('/System/Volumes/Data').free}
    path = PACKAGE / 'logs/internal_archive_relocation_candidate_receipt_v4.json'
    path.write_text(json.dumps(receipt, indent=2) + '\n')
    (SSD / 'manifests/internal_archive_relocation_candidate_receipt_v4.json').write_bytes(path.read_bytes())
    print(json.dumps({'files_verified': len(proof), 'copied_bytes': total, 'original_changes': False}), flush=True)

if __name__ == '__main__':
    main()
