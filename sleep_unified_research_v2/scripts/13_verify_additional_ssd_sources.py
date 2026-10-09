#!/usr/bin/env python3
"""Hash only exact missing-source filename candidates without changing them."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

PACKAGE = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4*1024**2), b''):
            h.update(b)
    return h.hexdigest()

def main():
    receipt = PACKAGE/'logs/additional_ssd_sources_verification_v2.json'
    table = PACKAGE/'tables/additional_ssd_source_recovery_v2.tsv'
    if receipt.exists() or table.exists():
        raise SystemExit('PRIOR_HASH_AUDIT_PRESERVED')
    inventory = PACKAGE/'logs/missing_core_ssd_filename_search_v2.json'
    plan = PACKAGE/'manifests/missing_core_archive_search_plan_v2.json'
    search = json.loads(inventory.read_text())
    archive = json.loads(plan.read_text())
    started = time.time()
    rows = []
    for candidate in search['candidates']:
        path = Path(candidate['path'])
        before = path.lstat()
        expected_size = archive['expected_sizes'][candidate['expected_sha256']]
        row = dict(candidate, expected_bytes=expected_size, actual_sha256='', stable_during_hash=False)
        if path.is_symlink() or before.st_size != expected_size:
            row['status'] = 'SYMLINK_OR_DIFFERENT_SIZE_NOT_ADMITTED'
        else:
            row['actual_sha256'] = sha(path)
            after = path.lstat()
            row['stable_during_hash'] = (before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino)
            row['status'] = 'EXACT_EXPECTED_SHA256_READ_ONLY_SOURCE_RECOVERED' if row['stable_during_hash'] and row['actual_sha256'] == row['expected_sha256'] else 'UNSTABLE_OR_MISMATCH_NOT_ADMITTED'
        rows.append(row)
        print(json.dumps({'name': path.name, 'status': row['status']}), flush=True)
    with table.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    hashes = {r['expected_sha256'] for r in rows if r['status'].startswith('EXACT_EXPECTED')}
    result = {'completed_utc': datetime.now(timezone.utc).isoformat(), 'elapsed_seconds': time.time()-started,
              'unique_missing_source_hashes_recovered': len(hashes), 'candidate_paths_verified': len(rows),
              'ledger_entries_satisfied': sum(len(archive['targets_by_sha256'][h]) for h in hashes),
              'rows': rows, 'search_receipt_sha256': sha(inventory), 'search_plan_sha256': sha(plan),
              'table_sha256': sha(table), 'script_sha256': sha(Path(__file__)),
              'read_buffer_bytes': 4*1024**2, 'other_projects_and_original_inputs_modified': False,
              'source_bytes_copied_or_committed': False, 'native_preprocessing_or_estimation_completed': False}
    with receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({'receipt': str(receipt), 'unique_hashes': len(hashes), 'ledger_entries_satisfied': result['ledger_entries_satisfied']}), flush=True)

if __name__ == '__main__':
    main()
