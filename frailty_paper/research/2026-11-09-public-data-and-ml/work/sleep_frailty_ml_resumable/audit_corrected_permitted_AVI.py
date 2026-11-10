"""Independently replay corrected score identities, assay roles and partitions."""
import ast
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import struct
import numpy as np
from audit_permitted_AVI_alignment import ROOT, BASE, DATA, OUT, sha


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def main():
    receipt_path = OUT / 'corrected_permitted_AVI_preparation_receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['prepared']
    matrix_path = OUT / 'corrected_permitted_AVI_physical_bank.npy'
    table_path = OUT / 'corrected_permitted_AVI_exact_identities.csv.gz'
    assert sha(matrix_path) == receipt['corrected_array_sha256']
    assert sha(table_path) == receipt['corrected_identity_table_sha256']
    matrix = np.load(matrix_path, allow_pickle=False)
    with gzip.open(table_path, 'rt') as f:
        rows = list(csv.DictReader(f))
    ids = np.array([int(r['bank_id']) for r in rows], dtype='int64')
    assert len(ids) == len(set(map(int, ids))) == receipt['used_physical_pairs']
    initial_path = OUT / 'permitted_AVI_alignment_receipt.json'
    assert sha(initial_path) == receipt['initial_alignment_receipt_sha256']
    initial = json.loads(initial_path.read_text())
    resolution_path = OUT / 'canonical_reference_role_resolution_receipt.json'
    assert sha(resolution_path) == receipt['reference_role_resolution_receipt_sha256']
    resolution = json.loads(resolution_path.read_text())
    role_table = OUT / 'canonical_reference_role_resolution.csv'
    assert sha(role_table) == resolution['resolution_table_sha256']
    with role_table.open() as f:
        roles = {int(r['bank_id']): r for r in csv.DictReader(f)}
    # Sources are read-only. The newly prepared query cache is fully committed.
    caches = [BASE / 'sleep_frailty_ml_v7/avi_exact_scores.sqlite', DATA / 'coverage_AVI_exact_scores.sqlite',
              OUT / 'canonical_reference_AVI_query_cache.sqlite']
    assert not Path(str(caches[-1]) + '-journal').exists()
    for path in caches[:2]:
        assert sha(path) == initial['source_cache_sha256'][str(path)]
    assert sha(caches[-1]) == receipt['new_canonical_cache_sha256']
    connections = [sqlite3.connect(p.as_uri() + '?mode=ro&immutable=1', uri=True) for p in caches]
    zero, repaired, values_by_key = 0, 0, {}
    for row in rows:
        i = int(row['bank_id'])
        canonical = row['canonical_GRCh38_score_variant']
        c, p, ref, alt = canonical.split(':')
        old = row['original_assay_variant_tag'].split(':')
        assert (c, p) == (old[0], old[1])
        reversed_role = i in roles
        assert row['source_assay_reference_is_genomic_reference'] == str(not reversed_role)
        if reversed_role:
            assert (ref, alt) == (old[3], old[2])
            assert canonical == roles[i]['canonical_genomic_variant']
            connection = connections[2]
            repaired += 1
        else:
            assert canonical == row['original_assay_variant_tag']
            connection = connections[0 if i < 277740 else 1]
        value = connection.execute('SELECT raw,phred,available FROM scores WHERE variant=?', (canonical,)).fetchone()
        assert value is not None
        available = int(row['available'])
        assert available == value[2]
        wanted = np.array([value[0], value[1], 1] if available else [np.nan, np.nan, 0], dtype='float32')
        assert np.array_equal(matrix[i], wanted, equal_nan=True)
        if available:
            assert float(row['raw']) == value[0] and float(row['phred']) == value[1]
            assert np.isfinite(matrix[i]).all()
            zero += int(value[0] == 0)
        else:
            assert row['raw'] == row['phred'] == '' and np.isnan(matrix[i, :2]).all()
        if canonical in values_by_key:
            assert np.array_equal(values_by_key[canonical], matrix[i], equal_nan=True)
        else:
            values_by_key[canonical] = matrix[i].copy()
    for connection in connections:
        connection.close()
    assert repaired == 601 and zero == 10
    assert matrix[279177, 0] == 0 and matrix[279177, 2] == 1
    # Re-read complete original source DNA for every repaired reference role.
    metadata = json.loads((ROOT / 'outputs/sleep_frailty_ml_resumable/development/frozen_reader_none_42_started.json').read_text())
    actual = {}
    for entry in metadata['physical_source_receipt']['sources']:
        assert sha(Path(entry['path'])) == entry['sha256']
        with gzip.open(entry['path'], 'rt') as f:
            for r in csv.DictReader(f):
                i = int(r['bank_id'])
                if i in roles:
                    actual[i] = r
    assert actual.keys() == roles.keys()
    reader_path = ROOT / 'work/sleep_frailty_ml_v7/reference.py'
    tree = ast.parse(reader_path.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'TwoBit')
    ns = {'struct': struct}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(reader_path), 'exec'), ns)
    genome = ns['TwoBit'](BASE / 'sleep_frailty_ml_v6/sources/official_hg38.2bit')
    payloads = {}
    for entry in resolution['public_GRCh37_reference_payloads']:
        path = ROOT / entry['path']
        assert sha(path) == entry['sha256']
        data = json.loads(path.read_text())
        assert data['genome'] == 'hg19'
        payloads[entry['sha256']] = data
    assembly_counts = collections.Counter()
    for i, role in roles.items():
        source = actual[i]
        a, b, j = source['reference'], source['alternate'], int(source['changed_position'])
        assert a[:j] == b[:j] and a[j + 1:] == b[j + 1:] and a[j] != b[j]
        if role['full_reference_match_assembly'] == 'GRCh38':
            chrom, p, ref, alt = role['canonical_genomic_variant'].split(':')
            begin = int(p) - 1 - (j if role['source_assay_orientation_to_GRCh38'] == '+' else len(b) - 1 - j)
            truth = genome.sequence(chrom, begin, begin + len(b))
            expected = b if role['source_assay_orientation_to_GRCh38'] == '+' else rc(b)
            assert truth == expected
        else:
            assert role['full_reference_match_assembly'] == 'GRCh37'
            payload = payloads[role['reference_payload_sha256']]
            assert payload['chrom'] == role['GRCh37_chromosome']
            begin = int(role['GRCh37_span_start_0_based'])
            off = begin - payload['start']
            truth = payload['dna'][off:off + len(a)].upper()
            matches = [(label, strand) for label, s in [('reference', a), ('alternate', b)]
                       for strand, candidate in [('+', s), ('-', rc(s))] if candidate == truth]
            assert matches == [(role['matched_assay_allele'], role['GRCh37_source_assay_orientation'])]
            position = begin + (j if role['GRCh37_source_assay_orientation'] == '+' else len(a) - 1 - j)
            assert position == int(role['GRCh37_variant_position_0_based'])
        assembly_counts[role['full_reference_match_assembly']] += 1
    # Canonicalizing an allele name must not reveal cross-partition copies.
    assert sha(DATA / 'coverage_training_arrays.npz') == receipt['training_data_sha256']
    with np.load(DATA / 'coverage_training_arrays.npz', allow_pickle=False) as z:
        banks, groups, splits = z['bank_index'], z['group'], z['split']
    assert np.array_equal(ids, np.unique(banks))
    gmin = np.full(len(matrix), np.iinfo('int32').max, dtype='int32')
    gmax = np.full(len(matrix), -1, dtype='int32')
    mask = np.zeros(len(matrix), dtype='uint8')
    np.minimum.at(gmin, banks, groups)
    np.maximum.at(gmax, banks, groups)
    np.bitwise_or.at(mask, banks, (1 << splits).astype('uint8'))
    assert np.all(gmin[ids] == gmax[ids])
    assert np.all(np.isin(mask[ids], [1, 2, 4]))
    by_key = collections.defaultdict(set)
    by_partition = collections.defaultdict(set)
    for row in rows:
        i = int(row['bank_id'])
        key = row['canonical_GRCh38_score_variant']
        by_key[key].add(int(gmin[i]))
        by_partition[key].add(int(mask[i]))
    assert max(map(len, by_key.values())) == 1
    assert max(map(len, by_partition.values())) == 1
    result = {'passed': True, 'used_physical_pairs': len(ids), 'exact_canonical_cache_rows_and_array_values_checked': len(ids),
        'whole_construct_reference_roles_independently_replayed': repaired, 'reference_assembly_counts': dict(assembly_counts),
        'valid_raw_zero_scores_preserved': zero, 'formerly_corrupted_zero_input_restored': True,
        'canonical_variants_have_consistent_score_values': True, 'canonical_variants_stay_in_one_existing_region_and_partition': True,
        'all_available_input_scores_are_finite': bool(np.isfinite(matrix[ids]).all()),
        'old_source_DNA_targets_arrays_and_models_modified': False, 'reserved_outcomes_read': False,
        'scope': 'Exact cached score values, corrected canonical keys, original full source construct/reference matches, zero-value handling and pre-existing region/split identity. Does not establish score predictive benefit or full remote-archive byte equivalence.',
        'preparation_receipt_sha256': sha(receipt_path), 'array_sha256': sha(matrix_path), 'identity_table_sha256': sha(table_path),
        'code_sha256': sha(Path(__file__)), 'scientific_prediction_goal_complete': False}
    (OUT / 'corrected_permitted_AVI_independent_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
