"""Prepare a separately audited exact-allele AVI input bank; no model fitting.

Only the designated downloadable composite SNV raw/PHRED scores are used.
Ordinary API/Atlas predictions and experimental outcomes are never read.
"""
import ast
import bisect
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import struct
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable/AVI_preparation'
DATA = ROOT / 'outputs/sleep_frailty_ml_next'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


class ReferenceBases:
    """Read official twoBit pages and check N blocks with binary search."""
    def __init__(self, path):
        self.f = path.open('rb')
        signature = self.f.read(4)
        self.order = '>' if signature == b'\x1aA\'C' else '<'
        self.f.seek(0)
        magic, version, n, reserved = self.ints(4)
        assert magic == 0x1A412743 and version == reserved == 0
        self.index, self.layout = {}, {}
        for i in range(n):
            name = self.f.read(self.f.read(1)[0]).decode()
            self.index[name] = self.ints(1)[0]

    def ints(self, n):
        return struct.unpack(self.order + str(n) + 'I', self.f.read(4 * n)) if n else ()

    def bases(self, chrom, positions):
        self.f.seek(self.index[chrom])
        length, count = self.ints(2)
        starts, lengths = self.ints(count), self.ints(count)
        masks = self.ints(1)[0]
        self.f.seek(8 * masks, 1)
        assert self.ints(1)[0] == 0
        offset = self.f.tell()
        assert tuple(sorted(starts)) == starts
        result, page, raw = {}, -1, b''
        for p in sorted(positions):
            assert 0 <= p < length
            wanted = (p // 4) // 4096
            if wanted != page:
                page = wanted
                self.f.seek(offset + page * 4096)
                raw = self.f.read(min(4096, (length + 3) // 4 - page * 4096))
            value = 'TCAG'[(raw[p // 4 - page * 4096] >> (6 - 2 * (p % 4))) & 3]
            i = bisect.bisect_right(starts, p) - 1
            if i >= 0 and p < starts[i] + lengths[i]:
                value = 'N'
            result[p] = value
        return result


def main():
    started = time.monotonic()
    OUT.mkdir(exist_ok=True)
    metadata_path = ROOT / 'outputs/sleep_frailty_ml_resumable/development/frozen_reader_none_42_started.json'
    metadata = json.loads(metadata_path.read_text())
    plan = json.loads((ROOT / 'outputs/sleep_frailty_ml_resumable/finite_stable_native_plan.json').read_text())
    data_path = DATA / 'coverage_training_arrays.npz'
    assert sha(data_path) == plan['data_sha256']
    with np.load(data_path, allow_pickle=False) as z:
        bank_index = z['bank_index']
    used = set(map(int, np.unique(bank_index)))
    n = int(bank_index.max()) + 1
    first_receipt = json.loads((BASE / 'sleep_frailty_ml_v8/outputs/AVI_comparison_input_receipt.json').read_text())
    second_receipt = json.loads((BASE / 'sleep_frailty_ml_v9/outputs/extra_feature_receipt.json').read_text())
    third_receipt = json.loads((DATA / 'coverage_AVI_input_receipt.json').read_text())
    files = [BASE / 'sleep_frailty_ml_v8/avi_training_bank.npy',
             BASE / 'sleep_frailty_ml_v9/extra_AVI.npy', DATA / 'coverage_appended_AVI.npy']
    expected = [first_receipt['V8_training_bank_sha256'], second_receipt['outputs']['extra_AVI.npy'], third_receipt['appended_bank_sha256']]
    for path, value in zip(files, expected):
        assert sha(path) == value
    banks = [np.load(p, allow_pickle=False) for p in files]
    assert [len(b) for b in banks] == [271904, 5836, 13817]
    values = np.concatenate(banks)
    assert values.shape == (n, 3)
    designation = BASE / 'sleep_frailty_ml_v7/outputs/AVI_permissive_artifact_receipt.json'
    assert sha(designation) == third_receipt['official_permissive_artifact_receipt_sha256']
    index = BASE / 'sleep_frailty_ml_v7/sources/AVI_SNVS.tsv.gz.tbi'
    assert hashlib.md5(index.read_bytes()).hexdigest() == third_receipt['index_MD5']
    db_paths = [BASE / 'sleep_frailty_ml_v7/avi_exact_scores.sqlite', DATA / 'coverage_AVI_exact_scores.sqlite']
    db_hashes = {str(p): sha(p) for p in db_paths}
    dbs = [sqlite3.connect(p.as_uri() + '?mode=ro&immutable=1', uri=True) for p in db_paths]
    seen = set()
    records = []
    requested = collections.defaultdict(set)
    counts, sources = collections.Counter(), collections.Counter()
    failures = []
    complement = str.maketrans('ACGT', 'TGCA')
    for source in metadata['physical_source_receipt']['sources']:
        path = Path(source['path'])
        assert sha(path) == source['sha256']
        with gzip.open(path, 'rt') as f:
            for row in csv.DictReader(f):
                i = int(row['bank_id'])
                if i not in used:
                    continue
                assert i not in seen
                seen.add(i)
                key = row['variant_GRCh38']
                c, ps, ref, alt = key.split(':')
                p = int(ps)
                assert len(ref) == len(alt) == 1 and ref != alt and set(ref + alt) <= set('ACGT')
                assert row['chromosome'] == c and int(row['position_hg38']) == p
                a, b, j = row['reference'], row['alternate'], int(row['changed_position'])
                assert len(a) == len(b) and 0 <= j < len(a)
                assert a[:j] == b[:j] and a[j + 1:] == b[j + 1:] and a[j] != b[j]
                pair = a[j] + b[j]
                if pair == ref + alt:
                    strand = '+'
                elif pair.translate(complement) == ref + alt:
                    strand = '-'
                else:
                    strand = 'unresolved'
                    failures.append({'bank_id': i, 'variant': key, 'reason': 'source_construct_alleles_disagree_with_tag'})
                counts['source_assay_orientation_' + strand] += 1
                cached = dbs[0 if i < 277740 else 1].execute('SELECT raw,phred,available,source FROM scores WHERE variant=?', (key,)).fetchone()
                assert cached is not None, key
                raw, phred, available, provenance = cached
                assert available in [0, 1]
                expected_values = np.array([raw, phred, 1] if available else [np.nan, np.nan, 0], dtype='float32')
                if not np.array_equal(values[i], expected_values, equal_nan=True):
                    failures.append({'bank_id': i, 'variant': key, 'reason': 'raw_score_bank_not_equal_to_exact_allele_cache'})
                if available:
                    assert np.isfinite(expected_values).all() and phred >= 0
                counts['available' if available else 'unavailable'] += 1
                sources[provenance] += 1
                requested[c].add(p - 1)
                records.append({'bank_id': i, 'variant_GRCh38': key, 'chromosome': c, 'position_0_based': p - 1,
                                'ref': ref, 'source_assay_strand_from_actual_changed_bases': strand,
                                'raw': raw if available else '', 'phred': phred if available else '',
                                'available': available, 'source_cache': str(db_paths[0 if i < 277740 else 1])})
    assert seen == used
    print('AVI BANK EXACT-CACHE AND PHYSICAL-ALLELE CHECKED', len(records), len(failures), flush=True)
    genome_path = BASE / 'sleep_frailty_ml_v6/sources/official_hg38.2bit'
    genome_receipt = json.loads((BASE / 'sleep_frailty_ml_v6/outputs/additional_reference_audit.json').read_text())
    assert sha(genome_path) == genome_receipt['sha256']
    reference = ReferenceBases(genome_path)
    bases = {c: reference.bases(c, points) for c, points in requested.items()}
    # Independently compare the new batched decoder with the existing checked
    # twoBit reader on a fixed position sample spanning all chromosomes.
    reader_path = ROOT / 'work/sleep_frailty_ml_v7/reference.py'
    tree = ast.parse(reader_path.read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'TwoBit')
    namespace = {'struct': struct}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(reader_path), 'exec'), namespace)
    independent = namespace['TwoBit'](genome_path)
    crosschecks = 0
    for chrom, positions in sorted(requested.items()):
        ordered = sorted(positions)
        for index in np.unique(np.linspace(0, len(ordered) - 1, min(24, len(ordered)), dtype='int64')):
            p = ordered[index]
            assert independent.sequence(chrom, p, p + 1) == bases[chrom][p]
            crosschecks += 1
    for row in records:
        base = bases[row['chromosome']][row['position_0_based']]
        row['official_GRCh38_reference_base'] = base
        row['reference_allele_matches_official_GRCh38'] = base == row['ref']
        if base != row['ref']:
            failures.append({'bank_id': row['bank_id'], 'variant': row['variant_GRCh38'], 'reason': 'reference_allele_disagrees_with_official_GRCh38', 'genome_base': base})
    for db in dbs:
        db.close()
    for path, value in db_hashes.items():
        assert sha(Path(path)) == value
    for path, value in zip(files, expected):
        assert sha(path) == value
    table = OUT / 'permitted_AVI_exact_bank_identity_audit.csv.gz'
    with gzip.open(table, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda r: r['bank_id']))
    ready_path = OUT / 'permitted_AVI_used_physical_bank.npy'
    if not failures:
        ready = np.full((n, 3), np.nan, dtype='float32')
        ready[:, 2] = 0
        selected = np.array(sorted(used), dtype='int64')
        ready[selected] = values[selected]
        np.save(ready_path, ready)
    receipt = {'passed': not failures, 'current_physical_bank_slots': n,
        'used_physical_pairs_audited': len(used), 'model_observation_rows': len(bank_index),
        'distinct_genomic_variants': len({r['variant_GRCh38'] for r in records}),
        'exact_genomic_positions_checked': sum(len(v) for v in requested.values()),
        'source_orientation_and_availability_counts': dict(counts), 'cached_score_provenance_counts': dict(sources),
        'failures': failures, 'reference_decoder_independent_position_crosschecks': crosschecks,
        'genome_sha256': genome_receipt['sha256'], 'reference_reader_sha256': sha(reader_path),
        'source_score_array_sha256': {str(p): h for p, h in zip(files, expected)},
        'source_cache_sha256': db_hashes, 'official_permissive_designation_receipt_sha256': sha(designation),
        'scope': 'Exact changed bases in each actual paired source construct, complete GRCh38 variant keys, official-genome reference bases, original bank order and raw/PHRED/availability cache values. This is not an original GWAS effect-allele convention audit or a whole-genomic-background equivalence claim.',
        'permitted_inputs_only': 'Designated downloadable composite AVI SNV raw and PHRED scores plus availability; never ordinary API tracks, embeddings, attributions or experimental labels.',
        'mirror_scope': 'Previously recorded public GRR mirror of the designated artifact; this audit does not compare every byte of the full official archive.',
        'missing_is_zero_biological_effect': False, 'ordinary_AlphaGenome_outputs_used': False,
        'model_targets_or_reserved_outcomes_read': False, 'current_running_comparison_changed': False,
        'new_models_trained': False, 'new_independent_study_success_established': False,
        'future_requirement': 'Only after a justified frozen sequence predictor: matched rows, common shared initialization, batches, epochs and update counts; zero-initialized added input columns; train-only imputation/scaling. Keep AVI only after a locked matched outside-study improvement test.',
        'training_data_sha256': sha(data_path), 'identity_audit_table_sha256': sha(table),
        'prepared_used_bank_sha256': sha(ready_path) if not failures else None,
        'code_sha256': sha(Path(__file__)), 'elapsed_seconds': time.monotonic() - started}
    (OUT / 'permitted_AVI_alignment_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['passed', 'used_physical_pairs_audited', 'distinct_genomic_variants', 'exact_genomic_positions_checked', 'source_orientation_and_availability_counts', 'reference_decoder_independent_position_crosschecks', 'failures', 'elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
