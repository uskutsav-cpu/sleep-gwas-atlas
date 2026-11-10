"""Independently compare every used physical row and training normalization."""
import datetime
import json
from pathlib import Path
import resource

import numpy as np

from readonly_native_inputs import ROOT, DATA, CHANNEL, load_readonly_inputs, sha

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def main():
    destination = OUT / 'AVI_preparation/readonly_native_input_equivalence_audit.json'
    assert not destination.exists()
    d, pairs, bank, physical, channels, metadata = load_readonly_inputs()
    first = np.load(DATA / 'cell_matched_native_development_features.npy', mmap_mode='r')
    second = np.load(DATA / 'coverage_appended_reader_features.npy', mmap_mode='r')
    pfirst = np.load(DATA / 'native_development_physical_metadata.npy', mmap_mode='r')
    psecond = np.load(DATA / 'coverage_appended_physical_metadata.npy', mmap_mode='r')
    used = np.unique(d['bank_index'])
    maximum = dict(reader=0., physical=0.)
    for start in range(0, len(used), 4096):
        ids = used[start:start + 4096]
        # Independently fetch each old or appended source row; do not use the
        # new concatenation class to construct the expected result.
        expected_reader = np.stack([first[i] if i < len(first) else second[i - len(first)] for i in ids])
        pm = np.stack([pfirst[i] if i < len(pfirst) else psecond[i - len(pfirst)] for i in ids])
        expected_physical = np.empty((len(ids), 2), dtype='float32')
        expected_physical[:, 0] = pm[:, 0] / 230.
        expected_physical[:, 1] = pm[:, 1] / np.maximum(pm[:, 0] - 1, 1)
        assert np.array_equal(bank[ids], expected_reader)
        assert np.array_equal(physical[ids, 48:50], expected_physical)
        assert np.isfinite(expected_reader).all() and np.isfinite(channels[ids]).all()
        maximum['reader'] = max(maximum['reader'], float(abs(bank[ids] - expected_reader).max()))
        maximum['physical'] = max(maximum['physical'], float(abs(physical[ids, 48:50] - expected_physical).max()))
    ids = np.array([used[-1], used[0], used[-1], used[len(used) // 2]])
    assert np.array_equal(bank[ids], np.stack([bank[int(i)] for i in ids]))
    assert np.array_equal(bank[:4], np.stack([bank[i] for i in range(4)]))
    try:
        physical[ids, 49:50]
    except IndexError:
        pass
    else:
        raise AssertionError('An unsupported feature access did not fail closed.')
    train = np.flatnonzero(d['split'] == 0)
    training_banks = np.unique(d['bank_index'][train])
    meta = json.loads((OUT / 'development/frozen_reader_none_42_started.json').read_text())
    assert __import__('hashlib').sha256(train.tobytes()).hexdigest() == meta['training_index_sha256']
    total, squares, count = np.zeros(236), np.zeros(236), 0
    for start in range(0, len(training_banks), 4096):
        ix = training_banks[start:start + 4096]
        raw = bank[ix]
        scalar = np.concatenate([np.column_stack([(raw[:, r, 0] + raw[:, r, 1]) * .5, raw[:, r, 3:16]]) for r in range(3)], axis=1)
        projected = channels[ix].astype('float32')
        projected[:, :, [1, 3], :] = abs(projected[:, :, [1, 3], :])
        values = np.column_stack([scalar, physical[ix, 48:50], projected.reshape(len(ix), 192)]).astype('float32').astype('float64')
        total += values.sum(axis=0); squares += (values * values).sum(axis=0); count += len(ix)
    mean = total / count
    sd = np.sqrt(np.maximum(squares / count - mean * mean, 0))
    sd[sd < 1e-5] = 1
    mean, sd = mean.astype('float32'), sd.astype('float32')
    assert np.array_equal(mean, np.asarray(meta['normalization']['mean'], dtype='float32'))
    assert np.array_equal(sd, np.asarray(meta['normalization']['sd'], dtype='float32'))
    result = {'passed': True, 'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'all_used_physical_rows_checked': len(used), 'observations': len(d['bank_index']),
              'training_physical_rows_for_normalization': len(training_banks),
              'maximum_reader_and_physical_column_differences': maximum,
              'full_training_mean_and_SD_bitwise_equal_to_original_fit': True,
              'out_of_order_duplicate_scalar_and_slice_access_checked': True,
              'unsupported_physical_columns_raise': True,
              'avoided_whole_array_allocation_bytes': int(np.prod(bank.shape) * bank.dtype.itemsize + len(bank) * 66 * 4),
              'audit_peak_RSS_bytes_on_macOS': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'inputs': metadata,
              'current_training_code_or_model_changed': False, 'new_scientific_models_fitted': False,
              'reserved_new_study_outcomes_read': False,
              'remaining': 'Use only in a separately specified future experiment after native runtime and training checks; do not replace the running loader.',
              'code_sha256': {p.name: sha(p) for p in [Path(__file__), Path(__file__).with_name('readonly_native_inputs.py')]}}
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['passed', 'all_used_physical_rows_checked', 'full_training_mean_and_SD_bitwise_equal_to_original_fit', 'avoided_whole_array_allocation_bytes', 'audit_peak_RSS_bytes_on_macOS']}), flush=True)


if __name__ == '__main__':
    main()
