#!/usr/bin/env python3
"""Exploratory group-held-out genetic-correlation predictor; not clinical software."""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np

FEATURES = ('adhd', 'mdd', 'scz', 'bipolar', 'alz', 'parkinson')
TARGETS = ('FI', 'FFS', 'HFRS')
GRID = (.1, 1., 10., 100.)


def group(name):
    if name in ('sleepdur', 'shortsleep', 'longsleep'):
        return 'dashti_2019_duration'
    if name in ('sleep_efficiency', 'accel_sleep_duration', 'sleep_timing'):
        return 'jones_accelerometer'
    return name


def load_matrix(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    names = list(dict.fromkeys(r['sleep_trait'] for r in rows))
    table = {}
    for r in rows:
        key = (r['sleep_trait'], r['outcome'])
        if key in table:
            raise ValueError('duplicate sleep trait / outcome pair')
        value = float(r['rg'])
        if not np.isfinite(value) or abs(value) >= 1:
            raise ValueError('matrix correlations must be finite and strictly inside (-1, 1)')
        table[key] = value
    if set(table) != {(n, k) for n in names for k in FEATURES + TARGETS}:
        raise ValueError('matrix must have exactly six features and three targets per sleep trait')
    return names, np.array([[table[n, k] for k in FEATURES] for n in names]), np.array([[table[n, k] for k in TARGETS] for n in names])


def ridge(train_x, train_y, test_x, lam):
    mu = train_x.mean(axis=0)
    sd = train_x.std(axis=0)
    sd = np.where(sd > 1e-12, sd, 1.)
    a = (train_x - mu) / sd
    v = (test_x - mu) / sd
    z = np.arctanh(train_y)
    coef = np.linalg.solve(a.T @ a + lam * np.eye(len(FEATURES)), a.T @ (z - z.mean()))
    return np.tanh(z.mean() + v @ coef)


def input_box(train_x, train_y, test_features, feature_ses, lam):
    """Extrema over a simultaneous one-SE feature box for a fixed model."""
    if feature_ses.shape != (len(FEATURES),) or not np.isfinite(feature_ses).all() or np.any(feature_ses < 0):
        raise ValueError('feature SEs must be six finite nonnegative numbers')
    if np.any(np.abs(test_features) + feature_ses > 1):
        raise ValueError('one-SE feature box must stay inside the correlation range')
    mu = train_x.mean(axis=0)
    sd = train_x.std(axis=0)
    sd = np.where(sd > 1e-12, sd, 1.)
    a = (train_x - mu) / sd
    z = np.arctanh(train_y)
    coef = np.linalg.solve(a.T @ a + lam * np.eye(len(FEATURES)), a.T @ (z - z.mean()))
    center = z.mean() + ((test_features - mu) / sd) @ coef
    radius = np.sum(np.abs(coef * feature_ses / sd))
    low, high = np.tanh([center - radius, center + radius])
    return {'low': float(low), 'high': float(high),
            'qualification': 'Simultaneous input stress range, not a confidence or prediction interval.'}


def feature_range_diagnostic(names, x, test_features, excluded_group):
    """Compare each input to observed ranges after removing the held-out group."""
    train = np.array([i for i, name in enumerate(names) if group(name) != excluded_group])
    if not len(train):
        raise ValueError('no training rows remain')
    low, high = x[train].min(axis=0), x[train].max(axis=0)
    return {
        'training_rows': int(len(train)),
        'features': {name: {'training_min': float(low[i]), 'training_max': float(high[i]),
                            'outside_training_range': bool(test_features[i] < low[i] or test_features[i] > high[i])}
                     for i, name in enumerate(FEATURES)},
        'qualification': 'Marginal ranges only: inside all ranges does not establish joint support or reliability.'
    }


def predict(names, x, y, test_features, excluded_group, feature_ses=None):
    train = np.array([i for i, n in enumerate(names) if group(n) != excluded_group])
    groups = list(dict.fromkeys(group(names[i]) for i in train))
    if len(groups) < 3 or len(train) < 4:
        raise ValueError('at least three training groups and four rows are required')
    mu, sd = x[train].mean(axis=0), x[train].std(axis=0)
    sd = np.where(sd > 1e-12, sd, 1.)
    nearest = train[np.argsort(np.sum(((x[train] - test_features) / sd) ** 2, axis=1), kind='stable')[:3]]
    outputs = {}
    for j, target in enumerate(TARGETS):
        scores = []
        for lam in GRID:
            errors = []
            for g in groups:
                sub = np.array([i for i in train if group(names[i]) != g])
                val = np.array([i for i in train if group(names[i]) == g])
                errors.extend((ridge(x[sub], y[sub, j], x[val], lam) - y[val, j]) ** 2)
            scores.append(float(np.mean(errors)))
        best = min(range(len(GRID)), key=lambda i: (scores[i], -GRID[i]))
        outputs[target] = {
            'bounded_ridge': float(ridge(x[train], y[train, j], test_features, GRID[best])),
            'lambda': GRID[best],
            'inner_group_mse': dict(zip(map(str, GRID), scores)),
            'training_mean': float(y[train, j].mean()),
            'knn3': float(y[nearest, j].mean()),
            'knn3_neighbors': [names[i] for i in nearest],
        }
        if feature_ses is not None:
            outputs[target]['one_se_input_box'] = input_box(x[train], y[train, j], test_features, feature_ses, GRID[best])
    return outputs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matrix', type=Path, default=Path(__file__).with_name('source_matrix.tsv'))
    ap.add_argument('--exclude-group', required=True, help='Hold out this sleep trait or related-measure group from training.')
    ap.add_argument('--features', required=True, type=Path, help='JSON object with six explicitly supplied genetic correlations.')
    ap.add_argument('--feature-ses', type=Path, help='Optional JSON object of six input standard errors for a fixed-model stress range.')
    args = ap.parse_args()
    names, x, y = load_matrix(args.matrix)
    if args.exclude_group not in {group(n) for n in names}:
        raise ValueError('excluded group not present in training matrix')
    values = json.loads(args.features.read_text())
    if set(values) != set(FEATURES):
        raise ValueError('features JSON must contain exactly: ' + ', '.join(FEATURES))
    test = np.array([float(values[k]) for k in FEATURES])
    if not np.isfinite(test).all() or np.any(np.abs(test) >= 1):
        raise ValueError('features must be finite correlations strictly inside (-1, 1)')
    ses = None
    if args.feature_ses is not None:
        raw_ses = json.loads(args.feature_ses.read_text())
        if set(raw_ses) != set(FEATURES):
            raise ValueError('feature SE JSON must contain exactly: ' + ', '.join(FEATURES))
        ses = np.array([float(raw_ses[k]) for k in FEATURES])
    result = {'excluded_group': args.exclude_group, 'features': dict(zip(FEATURES, test.tolist())),
              'feature_range_diagnostic': feature_range_diagnostic(names, x, test, args.exclude_group),
              'predictions': predict(names, x, y, test, args.exclude_group, ses),
              'warning': 'Exploratory correlations only; not individual or clinical predictions. See README.'}
    json.dump(result, sys.stdout, indent=2, allow_nan=False)
    print()


if __name__ == '__main__':
    main()
