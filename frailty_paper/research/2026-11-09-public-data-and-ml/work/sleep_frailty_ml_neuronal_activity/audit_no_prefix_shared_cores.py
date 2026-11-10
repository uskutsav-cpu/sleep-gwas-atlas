"""Independently check every actual 145bp window without prefix acceleration."""
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'
CHANNEL = ROOT / 'outputs/sleep_frailty_ml_channels/source_candidates'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def rows(path, delimiter=','):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        yield from csv.DictReader(f, delimiter=delimiter)


def main():
    start = time.monotonic()
    proof = OUT / 'additional_activity_closed_core_receipt.json'
    receipt = json.loads(proof.read_text())
    candidates_path = OUT / 'additional_activity_closed_core_candidate_screen.csv.gz'
    components_path = OUT / 'all_neuronal_source_identity_components.csv.gz'
    assert sha(candidates_path) == receipt['candidate_screen_sha256']
    assert sha(components_path) == receipt['all_source_components_sha256']
    all_components = list(rows(components_path))
    component = {r['name']: int(r['component_after_shared_145bp_closure']) for r in all_components}
    kept = [r for r in rows(candidates_path) if r['passes_shared_core_and_duplicate_screen'] == 'True']
    wanted = {int(r['source_component_after_shared_145bp_closure']) for r in kept}
    original = json.loads((OUT / 'additional_neuronal_activity_inventory_receipt.json').read_text())
    metadata = OUT / Path(original['sources'][0]['path']).name
    assert sha(metadata) == receipt['metadata_sha256']
    identities = [(r['name'], r['sequence'].upper()) for r in rows(metadata, '\t')]
    oriented_counts = collections.Counter(s for name, s in identities)
    lookup = {}
    for name, s in identities:
        if component[name] not in wanted:
            continue
        reverse = rc(s)
        for a in range(len(s) - 144):
            # Canonical key construction; the scanner below instead visits
            # both orientations explicitly, with no prefix selection.
            word = s[a:a + 145]
            other = reverse[len(s) - a - 145:len(s) - a]
            key = word if word < other else other
            if key in lookup:
                assert lookup[key] == component[name]
            lookup[key] = component[name]
    source_windows, source_hits = 0, 0
    for name, sequence in identities:
        for s in (sequence, rc(sequence)):
            for j in range(len(s) - 144):
                source_windows += 1
                root = lookup.get(s[j:j + 145])
                if root is not None:
                    source_hits += 1
                    assert component[name] == root, ('Source component closure incomplete', name, root)
    print('NO PREFIX ENTIRE SOURCE PASS', source_windows, source_hits, flush=True)
    spatial = json.loads((OUT / 'additional_activity_spatial_and_reference_receipt.json').read_text())
    profiles = []
    for entry in spatial['historical_sources']:
        assert sha(entry['path']) == entry['sha256']
        windows, constructs = 0, 0
        for row in rows(entry['path']):
            for field in ['reference', 'alternate']:
                sequence = row[field].upper()
                constructs += 1
                # Deliberately scan repeated full constructs again. This
                # independently bypasses the first guard's SHA deduplication.
                for s in (sequence, rc(sequence)):
                    for j in range(len(s) - 144):
                        windows += 1
                        assert s[j:j + 145] not in lookup, ('Historical shared core retained', entry['path'], row.get('variant_GRCh38'), j)
        profiles.append({'path': entry['path'], 'original_constructs_without_deduplication': constructs,
                         'full_145bp_windows_both_orientations': windows, 'retained_core_matches': 0})
        print('NO PREFIX HISTORY PASS', json.dumps(profiles[-1]), flush=True)
    reserved_path = CHANNEL / 'GSE87711_design.txt.gz'
    assert sha(reserved_path) == receipt['reserved_original_design_sha256']
    reserved_constructs, reserved_windows = 0, 0
    with gzip.open(reserved_path, 'rt') as f:
        for identifier, sequence in csv.reader(f, delimiter='\t'):
            reserved_constructs += 1
            for s in (sequence, rc(sequence)):
                for j in range(len(s) - 144):
                    reserved_windows += 1
                    assert s[j:j + 145] not in lookup, ('Reserved shared core retained', identifier)
    by_name = dict(identities)
    for row in kept:
        s = by_name[row['name']]
        canonical = min(s, rc(s))
        assert hashlib.sha256(canonical.encode()).hexdigest() == row['physical_sequence_sha256_both_orientations']
        assert oriented_counts[s] == 1
        assert row['full_actual_reference_match'] == 'True'
        assert row['source_logFC'] not in ['', 'nan', 'NA']
    r_path = OUT / 'neuronal_activity_independent_R_receipt.tsv'
    with r_path.open() as f:
        r_check = next(csv.DictReader(f, delimiter='\t'))
    # This external R check covers every original finite coefficient, including
    # the retained subset; it does not constitute a raw-count BCalm refit.
    assert r_check['passed'] == 'TRUE'
    result = {'passed': True, 'scope': 'All actual source rows and all historical/reserved DNA, including repeated constructs; no prefix acceleration and no sequence-digest deduplication in scanner.',
        'candidate_rows': len(kept), 'candidate_components': len(wanted), 'canonical_145bp_core_keys': len(lookup),
        'entire_source_rows': len(identities), 'entire_source_windows_both_orientations': source_windows,
        'source_core_matches_all_within_expected_component': source_hits,
        'history_profiles': profiles, 'reserved_original_constructs': reserved_constructs,
        'reserved_windows_both_orientations': reserved_windows, 'retained_history_or_reserved_core_matches': 0,
        'retained_constructs_have_unique_original_oriented_DNA_identifiers': True,
        'candidate_hashes_and_original_reference_match_flags_verified': True,
        'original_coefficient_independent_R_check_passed': True, 'original_R_receipt_sha256': sha(r_path),
        'raw_counts_refitted': False, 'new_independent_biological_study': False, 'reserved_outcomes_read': False,
        'models_fitted_or_changed': False, 'training_experiment_still_to_specify': True,
        'source_core_receipt_sha256': sha(proof), 'candidate_file_sha256': sha(candidates_path),
        'all_source_components_sha256': sha(components_path), 'code_sha256': sha(Path(__file__)),
        'elapsed_seconds': time.monotonic() - start}
    (OUT / 'additional_activity_independent_no_prefix_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
