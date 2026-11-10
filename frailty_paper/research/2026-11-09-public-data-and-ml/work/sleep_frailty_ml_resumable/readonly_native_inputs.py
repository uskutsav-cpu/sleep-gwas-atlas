"""Read existing physical features in batches instead of copying whole banks.

For a future separately specified AVI comparison only. The running scientific
comparison keeps its original loader and numerical history unchanged.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'outputs/sleep_frailty_ml_next'
CHANNEL = ROOT / 'outputs/sleep_frailty_ml_channels'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


class ConcatenatedReadOnlyBank:
    """The indexing used by the pinned native model, backed by two memmaps."""
    def __init__(self, first, second):
        self.first, self.second = first, second
        if first.ndim < 1 or first.shape[1:] != second.shape[1:] or first.dtype != second.dtype:
            raise ValueError('Physical feature banks must have identical trailing shapes and dtypes.')
        self.shape = (len(first) + len(second), *first.shape[1:])
        self.dtype = first.dtype

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, indices):
        if isinstance(indices, slice):
            indices = np.arange(len(self), dtype='int64')[indices]
        ids = np.asarray(indices)
        if ids.dtype == bool:
            if ids.ndim != 1 or len(ids) != len(self):
                raise IndexError('A boolean mask must cover the physical bank.')
            ids = np.flatnonzero(ids)
        if not np.issubdtype(ids.dtype, np.integer):
            raise IndexError('Physical bank indices must be integers.')
        if np.any(ids < 0) or np.any(ids >= len(self)):
            raise IndexError('Physical bank index out of bounds.')
        flat = ids.reshape(-1)
        result = np.empty((len(flat), *self.shape[1:]), dtype=self.dtype)
        first = flat < len(self.first)
        result[first] = self.first[flat[first]]
        result[~first] = self.second[flat[~first] - len(self.first)]
        return result.reshape((*ids.shape, *self.shape[1:]))


class NativePhysicalColumns:
    """Compute exactly the two original physical columns consumed by the model."""
    def __init__(self, metadata):
        self.metadata = metadata
        self.shape = (len(metadata), 66)

    def __getitem__(self, key):
        if not isinstance(key, tuple) or len(key) != 2:
            raise IndexError('Only the audited native [physical_ids, 48:50] access is supported.')
        indices, columns = key
        if not isinstance(columns, slice) or (columns.start, columns.stop, columns.step) != (48, 50, None):
            raise IndexError('No unimplemented physical or nucleotide columns may be substituted.')
        values = self.metadata[indices]
        if np.any(values[..., 0] <= 0):
            raise ValueError('Actual construct lengths must be positive.')
        return np.stack([values[..., 0] / 230.,
                         values[..., 1] / np.maximum(values[..., 0] - 1, 1)], axis=-1).astype('float32')


def load_readonly_inputs():
    expansion = json.loads((DATA / 'coverage_expansion_receipt.json').read_text())
    appended = json.loads((DATA / 'coverage_feature_receipt.json').read_text())
    original = json.loads((DATA / 'cell_matched_feature_receipt.json').read_text())
    channel = json.loads((CHANNEL / 'channel_feature_receipt.json').read_text())
    assert expansion['passed'] and appended['passed'] and channel['passed']
    files = {
        DATA / 'coverage_training_arrays.npz': expansion['training_arrays_sha256'],
        DATA / 'coverage_measured_pairs.npz': expansion['pairs_sha256'],
        DATA / 'cell_matched_native_development_features.npy': original['features_sha256'],
        DATA / 'coverage_appended_reader_features.npy': appended['features_sha256'],
        DATA / 'native_development_physical_metadata.npy': original['physical_metadata_sha256'],
        DATA / 'coverage_appended_physical_metadata.npy': appended['physical_metadata_sha256'],
        CHANNEL / 'projected_reader_channels.npy': channel['features_sha256']}
    for path, expected in files.items():
        assert sha(path) == expected, str(path)
    with np.load(DATA / 'coverage_training_arrays.npz', allow_pickle=False) as archive:
        observations = {key: archive[key] for key in archive.files}
    with np.load(DATA / 'coverage_measured_pairs.npz', allow_pickle=False) as archive:
        pairs = {key: archive[key] for key in archive.files}
    mmap = lambda path: np.load(path, mmap_mode='r', allow_pickle=False)
    bank = ConcatenatedReadOnlyBank(mmap(DATA / 'cell_matched_native_development_features.npy'),
                                   mmap(DATA / 'coverage_appended_reader_features.npy'))
    metadata = ConcatenatedReadOnlyBank(mmap(DATA / 'native_development_physical_metadata.npy'),
                                       mmap(DATA / 'coverage_appended_physical_metadata.npy'))
    physical = NativePhysicalColumns(metadata)
    channels = mmap(CHANNEL / 'projected_reader_channels.npy')
    assert len(bank) == len(metadata) == len(channels)
    return observations, pairs, bank, physical, channels, {
        'source_sha256': {str(path): expected for path, expected in files.items()},
        'whole_reader_bank_concatenation_allocated': False,
        'whole_66_column_physical_matrix_allocated': False,
        'physical_and_reader_values_changed': False,
        'not_used_by_running_native_comparison': True}
