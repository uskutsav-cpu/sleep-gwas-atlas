"""Build a new canonical, zero-preserving AVI bank; preserve all older results."""
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import json
import sqlite3
import sys
import threading
import time
from pathlib import Path
import numpy as np
from audit_permitted_AVI_alignment import ROOT, BASE, DATA, OUT, sha

sys.path.append(str(ROOT / 'work/stage5_env/lib/python3.13/site-packages'))
import pysam

URL = 'https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'


def main():
    start = time.monotonic()
    initial_path = OUT / 'permitted_AVI_alignment_receipt.json'
    initial = json.loads(initial_path.read_text())
    identity_path = OUT / 'permitted_AVI_exact_bank_identity_audit.csv.gz'
    assert sha(identity_path) == initial['identity_audit_table_sha256']
    resolution_path = OUT / 'canonical_reference_role_resolution_receipt.json'
    resolution = json.loads(resolution_path.read_text())
    assert resolution['reference_resolution_complete_without_network_errors']
    assert resolution['source_alignment_audit_sha256'] == sha(initial_path)
    role_path = OUT / 'canonical_reference_role_resolution.csv'
    assert sha(role_path) == resolution['resolution_table_sha256']
    with role_path.open() as f:
        roles = {int(r['bank_id']): r for r in csv.DictReader(f)}
    assert len(roles) == 601 and all(r['full_construct_reference_verified'] == 'True' for r in roles.values())
    with gzip.open(identity_path, 'rt') as f:
        identities = list(csv.DictReader(f))
    assert len(identities) == initial['used_physical_pairs_audited']
    wanted = {r['canonical_genomic_variant'] for r in roles.values()}
    index = BASE / 'sleep_frailty_ml_v7/sources/AVI_SNVS.tsv.gz.tbi'
    assert hashlib.md5(index.read_bytes()).hexdigest() == '911ee04b6e8fba1882465221b5a14f2c'
    cache_path = OUT / 'canonical_reference_AVI_query_cache.sqlite'
    db = sqlite3.connect(cache_path)
    db.execute('CREATE TABLE IF NOT EXISTS scores (variant TEXT PRIMARY KEY, raw REAL, phred REAL, available INTEGER, source TEXT)')
    old_paths = [BASE / 'sleep_frailty_ml_v7/avi_exact_scores.sqlite', DATA / 'coverage_AVI_exact_scores.sqlite',
                 BASE / 'sleep_frailty_ml_v8/fresh_AVI_exact_scores.sqlite']
    old_hashes = {str(p): sha(p) for p in old_paths}
    present = {r[0] for r in db.execute('SELECT variant FROM scores WHERE available=1')}
    copied = 0
    for path in old_paths:
        old = sqlite3.connect(path.as_uri() + '?mode=ro&immutable=1', uri=True)
        for key in sorted(wanted - present):
            row = old.execute('SELECT variant,raw,phred,available,source FROM scores WHERE variant=? AND available=1', (key,)).fetchone()
            if row:
                db.execute('INSERT OR REPLACE INTO scores VALUES (?,?,?,?,?)', row)
                present.add(key)
                copied += 1
        old.close()
    db.commit()
    grouped = collections.defaultdict(set)
    for key in wanted - present:
        c, p, a, b = key.split(':')
        grouped[(c, (int(p) - 1) // 16384)].add((int(p), a, b))
    local = threading.local()
    def query(item):
        (chrom, window), keys = item
        for attempt in range(3):
            try:
                if not hasattr(local, 'tab'):
                    local.tab = pysam.TabixFile(URL, index=str(index))
                found = {}
                for line in local.tab.fetch(chrom, min(k[0] for k in keys) - 1, max(k[0] for k in keys)):
                    f = line.split('\t')
                    assert f[0] == chrom
                    key = (int(f[1]), f[2], f[3])
                    if key in keys:
                        value = float(f[4]), float(f[5])
                        assert np.isfinite(value).all() and value[1] >= 0
                        assert key not in found or found[key] == value
                        found[key] = value
                return [(f'{chrom}:{p}:{a}:{b}', *found.get((p, a, b), (None, None)),
                         int((p, a, b) in found), 'exact_public_permissive_AVI_canonical_reference_query')
                        for p, a, b in keys]
            except Exception:
                if hasattr(local, 'tab'):
                    local.tab.close()
                    del local.tab
                if attempt == 2:
                    raise
                time.sleep(1 + attempt)
    failures, completed = [], 0
    print('CANONICAL PERMITTED AVI QUERY', len(wanted), 'alleles', copied, 'reused', len(grouped), 'remote windows', flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        tasks = {pool.submit(query, item): item[0] for item in sorted(grouped.items())}
        for future in concurrent.futures.as_completed(tasks):
            completed += 1
            try:
                db.executemany('INSERT OR REPLACE INTO scores VALUES (?,?,?,?,?)', future.result())
            except Exception as error:
                failures.append({'window': tasks[future], 'error': repr(error)})
            if completed % 20 == 0 or completed == len(grouped):
                db.commit()
                progress = {'completed_windows': completed, 'total_windows': len(grouped), 'failed_windows': len(failures),
                            'elapsed_seconds': time.monotonic() - start}
                (OUT / 'canonical_AVI_query_progress.json').write_text(json.dumps(progress, indent=2) + '\n')
                print('CANONICAL PERMITTED AVI QUERY PROGRESS', json.dumps(progress), flush=True)
    db.commit()
    attempts = {'network_failures': failures, 'requested_canonical_alleles': len(wanted),
                'reused_from_older_caches': copied, 'requested_remote_windows': len(grouped)}
    (OUT / 'canonical_AVI_query_attempt_receipt.json').write_text(json.dumps(attempts, indent=2) + '\n')
    assert not failures, failures
    ready = np.full((initial['current_physical_bank_slots'], 3), np.nan, dtype='float32')
    ready[:, 2] = 0
    records, counts = [], collections.Counter()
    for row in identities:
        i = int(row['bank_id'])
        original_key = row['variant_GRCh38']
        role = roles.get(i)
        if role:
            canonical = role['canonical_genomic_variant']
            cached = db.execute('SELECT raw,phred,available,source FROM scores WHERE variant=?', (canonical,)).fetchone()
            assert cached is not None, canonical
            raw, phred, available, source = cached
            counts['canonical_reference_roles_repaired'] += 1
        else:
            canonical = original_key
            available = int(row['available'])
            raw = float(row['raw']) if available else None
            phred = float(row['phred']) if available else None
            source = 'unchanged_exact_allele_cache_values_from_initial_identity_audit'
            assert row['reference_allele_matches_official_GRCh38'] == 'True'
        # Availability decides missingness. A real score of 0.0 stays 0.0.
        if available:
            assert np.isfinite([raw, phred]).all() and phred >= 0
            ready[i] = [raw, phred, 1]
            counts['available_used_pairs'] += 1
            counts['valid_raw_score_zero'] += int(raw == 0.)
        else:
            counts['unavailable_used_pairs'] += 1
        c, p, a, b = canonical.split(':')
        assert a == row['official_GRCh38_reference_base']
        records.append({'bank_id': i, 'original_assay_variant_tag': original_key,
            'canonical_GRCh38_score_variant': canonical, 'source_assay_reference_is_genomic_reference': role is None,
            'full_construct_role_resolution_assembly': role['full_reference_match_assembly'] if role else 'unchanged_original_tag_reference_base_verified',
            'raw': raw if available else '', 'phred': phred if available else '', 'available': available,
            'score_source': source, 'source_DNA_and_signed_targets_changed': False})
    db.close()
    used = np.array(sorted(int(r['bank_id']) for r in identities), dtype='int64')
    known = ready[used, 2] == 1
    assert np.isfinite(ready[used[known]]).all()
    assert np.isnan(ready[used[~known], :2]).all()
    original_files = [BASE / 'sleep_frailty_ml_v8/avi_training_bank.npy', BASE / 'sleep_frailty_ml_v9/extra_AVI.npy', DATA / 'coverage_appended_AVI.npy']
    original = np.concatenate([np.load(p, allow_pickle=False) for p in original_files])
    for p in original_files:
        assert sha(p) == initial['source_score_array_sha256'][str(p)]
    changed = ~np.isclose(ready[used], original[used], rtol=0, atol=0, equal_nan=True).all(axis=1)
    assert np.all(np.isin(used[changed], np.r_[list(roles), 279177]))
    assert ready[279177, 0] == 0. and ready[279177, 2] == 1.
    for path, expected in old_hashes.items():
        assert sha(Path(path)) == expected
    assert sha(DATA / 'coverage_training_arrays.npz') == initial['training_data_sha256']
    array_path = OUT / 'corrected_permitted_AVI_physical_bank.npy'
    np.save(array_path, ready)
    table_path = OUT / 'corrected_permitted_AVI_exact_identities.csv.gz'
    with gzip.open(table_path, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda r: r['bank_id']))
    receipt = {'prepared': True, 'original_failure_audit_preserved': True,
        'used_physical_pairs': len(used), 'distinct_canonical_variants': len({r['canonical_GRCh38_score_variant'] for r in records}),
        'counts': dict(counts), 'changed_used_input_slots': int(changed.sum()),
        'legitimate_zero_score_restored': True, 'canonicalization_reference_assembly_counts': resolution['full_construct_match_counts'],
        'correction_scope': 'Only the auxiliary genomic score key/values change. Original paired assay DNA, all measured magnitude and signed/activity targets, split/group identities and every original model result remain unchanged.',
        'score_interpretation': 'Genomic reference-to-alternate composite AVI severity used as an auxiliary magnitude-prediction input; not an experimental effect, source-contrast direction, causal risk or cell-specific gene prediction.',
        'source': URL, 'source_index_MD5': hashlib.md5(index.read_bytes()).hexdigest(),
        'network_attempts': attempts, 'reference_role_resolution_receipt_sha256': sha(resolution_path),
        'initial_alignment_receipt_sha256': sha(initial_path),
        'original_cache_hashes': old_hashes, 'new_canonical_cache_sha256': sha(cache_path),
        'corrected_array_sha256': sha(array_path), 'corrected_identity_table_sha256': sha(table_path),
        'training_data_sha256': initial['training_data_sha256'], 'current_model_fit_changed': False,
        'new_models_trained': False, 'reserved_study_outcomes_read': False,
        'AVI_retained_for_project_predictions': False,
        'remaining': ['Independent raw-row/key/zero-value replay of the corrected bank',
                      'A prespecified, matched training comparison after a justified sequence predictor exists',
                      'Untouched-study score/no-score benefit check; keep only if it passes'],
        'code_sha256': sha(Path(__file__)), 'elapsed_seconds': time.monotonic() - start}
    (OUT / 'corrected_permitted_AVI_preparation_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['prepared', 'used_physical_pairs', 'distinct_canonical_variants', 'counts', 'changed_used_input_slots', 'legitimate_zero_score_restored', 'elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
