"""Compact actual historical constructs for a full native runtime replay.

No reserved new-study observations are accessed. This reduces the memory needed
to check the new score connection while the original comparison is running.
"""
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from readonly_native_inputs import ROOT, load_readonly_inputs, sha

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'
AVI = OUT / 'AVI_preparation'


def main():
    destination = AVI / 'native_runtime_diagnostic_cases.npz'
    assert not destination.exists()
    d, pairs, bank, physical, channels, sources = load_readonly_inputs()
    del pairs
    previous = json.loads((OUT / 'adapt_last_stage_none_42_epoch1_independent_audit.json').read_text())
    rows = previous['saved_validation_prediction_rows_replayed']
    extra_banks = [279177]
    with (AVI / 'canonical_reference_role_resolution.csv').open() as handle:
        extra_banks.append(int(next(csv.DictReader(handle))['bank_id']))
    rows = np.array(rows + [int(np.flatnonzero(d['bank_index'] == b)[0]) for b in extra_banks], dtype='int64')
    original_banks = d['bank_index'][rows]
    meta = json.loads((OUT / 'development/frozen_reader_none_42_started.json').read_text())
    wanted = set(map(int, original_banks))
    constructs = {}
    sequence_sources = []
    for source in meta['physical_source_receipt']['sources']:
        path = Path(source['path'])
        assert sha(path) == source['sha256']
        sequence_sources.append(source)
        with gzip.open(path, 'rt') as handle:
            for record in csv.DictReader(handle):
                i = int(record['bank_id'])
                if i in wanted:
                    assert i not in constructs
                    a, b, pos = record['reference'], record['alternate'], int(record['changed_position'])
                    assert len(a) == len(b) and [j for j, (x, y) in enumerate(zip(a, b)) if x != y] == [pos]
                    constructs[i] = (a, b, pos)
    assert set(constructs) == wanted
    lut = {base: i for i, base in enumerate('AGCT')}
    reference, length, position, alternate = [], [], [], []
    for i in original_banks:
        a, b, pos = constructs[int(i)]
        ref = np.array([lut[base] for base in a], dtype='uint8')
        alt = ref.copy(); alt[pos] = lut[b[pos]]
        assert ''.join('AGCT'[j] for j in alt) == b
        reference.append(ref); length.append(len(ref)); position.append(pos); alternate.append(lut[b[pos]])
    assert sorted(set(length)) == [150, 199, 200, 270]
    np.savez_compressed(destination,
        observation_ids=rows, original_bank_ids=original_banks,
        reference_tokens=np.concatenate(reference), length=np.array(length, dtype='int16'),
        position=np.array(position, dtype='int16'), alt=np.array(alternate, dtype='uint8'),
        bank=bank[original_banks], physical=physical[original_banks, 48:50],
        channels=channels[original_banks],
        cell=d['cell'][rows], line=d['line'][rows], effect=d['effect'][rows],
        avi=np.load(AVI / 'corrected_permitted_AVI_physical_bank.npy', mmap_mode='r')[original_banks])
    receipt = {'prepared': True, 'actual_observations': len(rows),
               'historical_observation_ids': rows.tolist(), 'physical_bank_ids': original_banks.tolist(),
               'actual_source_DNA_lengths': sorted(set(length)),
               'complete_actual_REF_and_ALT_checked': True, 'invented_or_padded_bases': False,
               'original_sequence_sources': sequence_sources,
               'original_model_metadata_sha256': sha(OUT / 'development/frozen_reader_none_42_started.json'),
               'diagnostic_data_sha256': sha(destination),
               'reserved_new_study_outcomes_read': False,
               'input_preparation_not_model_training': True,
               'loader_source_hashes': sources['source_sha256'], 'code_sha256': sha(Path(__file__))}
    (AVI / 'native_runtime_diagnostic_cases_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ['prepared', 'actual_observations', 'actual_source_DNA_lengths', 'complete_actual_REF_and_ALT_checked']}))


if __name__ == '__main__':
    main()
