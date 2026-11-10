"""Train-only score scaling for a matched composite-score input experiment."""
import numpy as np


def validate(values):
    x = np.asarray(values)
    if x.ndim != 2 or x.shape[1] != 3 or not np.all(np.isin(x[:, 2], [0, 1])):
        raise ValueError('Expected raw AVI, PHRED and a 0/1 availability flag.')
    known = x[:, 2] == 1
    if not np.isfinite(x[known, :2]).all():
        raise ValueError('An available score cannot contain a missing value, including a lost numerical zero.')
    return known


def fit_scaling(bank, training_physical_ids):
    ids = np.unique(np.asarray(training_physical_ids, dtype='int64'))
    if not len(ids) or np.any(ids < 0) or np.any(ids >= len(bank)):
        raise ValueError('Nonempty, valid training-only physical IDs are required.')
    x = np.asarray(bank[ids], dtype='float64')
    known = validate(x)
    if not known.any():
        raise ValueError('No measured composite score inputs available in this training partition.')
    median = np.median(x[known, :2], axis=0)
    filled = np.where(known[:, None], x[:, :2], median)
    mean, sd = filled.mean(axis=0), filled.std(axis=0)
    sd[sd < 1e-5] = 1.
    return {'median': median.tolist(), 'mean': mean.tolist(), 'sd': sd.tolist(),
            'training_availability_mean': float(known.mean()),
            'training_unique_physical_pairs': len(ids), 'known_training_physical_pairs': int(known.sum()),
            'physical_ids_sha256': __import__('hashlib').sha256(ids.tobytes()).hexdigest(),
            'scope': 'Only training physical IDs. No phenotype, calibration or final-study labels.'}


def transform(values, scaling, use_scores):
    x = np.asarray(values, dtype='float64')
    known = validate(x)
    filled = np.where(known[:, None], x[:, :2], np.asarray(scaling['median']))
    raw = np.clip((filled - scaling['mean']) / scaling['sd'], -12., 12.)
    raw[~known] = 0.  # Neutral model feature, never a zero biological-effect label.
    if not use_scores:
        raw[:] = 0.
    # The availability information is IDENTICAL in both arms, so only score
    # values can explain an added-input difference. Constant flags center to 0.
    availability = known.astype('float64') - scaling['training_availability_mean']
    return np.column_stack([raw, availability]).astype('float32')
