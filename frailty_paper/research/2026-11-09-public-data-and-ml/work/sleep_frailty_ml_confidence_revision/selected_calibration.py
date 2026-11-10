"""Calibrate a fixed, outcome-independent answer-selection rule.

Empirical research intervals, not a cross-laboratory coverage guarantee.
The error critics are shared with the unchanged archived control.
"""
import numpy as np


def eligibility(package, features, baseline, supported_cell, supported_line, study_count, policy):
    critics = package['critics']
    raw = np.maximum(critics['quantile_error'].predict(features), 1e-4)
    estimated = np.maximum(critics['mean_error'].predict(features), 0)
    base = np.maximum(critics['baseline_error'].predict(features), 1e-4)
    z = np.max(abs((features - package['feature_center']) / package['feature_sd']), axis=1)
    selected = (supported_cell & supported_line
        & (study_count >= policy['minimum_training_studies_same_cell'])
        & (z <= policy['maximum_standardized_input'])
        & ((base - estimated) / base >= policy['minimum_estimated_relative_gain'])
        & (base - estimated >= policy['minimum_estimated_absolute_gain'])
        & (raw <= np.maximum(policy['maximum_halfwidth_absolute_floor'], baseline)))
    return selected, raw


def region_quantile(values, groups, level=.9):
    values = np.asarray(values, dtype='float64')
    groups = np.asarray(groups)
    assert values.size and np.isfinite(values).all() and np.all(values >= 0)
    _, inverse, counts = np.unique(groups, return_inverse=True, return_counts=True)
    weights = 1. / counts[inverse]
    order = np.argsort(values, kind='stable')
    cumulative = np.cumsum(weights[order])
    position = np.searchsorted(cumulative, level * cumulative[-1], side='left')
    return float(values[order[min(position, len(order) - 1)]])


def calibrate_selected(package, features, y, prediction, baseline, cells, sources,
                       groups, supported_cell, supported_line, study_count, policy):
    selected, raw = eligibility(package, features, baseline, supported_cell, supported_line, study_count, policy)
    ratio = abs(y - prediction) / raw
    scales = {}
    for cell in np.unique(cells):
        mask = selected & (cells == cell)
        count = int(mask.sum())
        nregions = len(np.unique(groups[mask]))
        # No unseen/sparse-cell fallback to measurements from different cells.
        if count < policy['minimum_scale_rows'] or nregions < policy['minimum_scale_regions']:
            scales[int(cell)] = {'supported': False, 'selected_rows': count, 'selected_regions': nregions}
            continue
        factor = region_quantile(ratio[mask], groups[mask])
        details = []
        for source in np.unique(sources[mask]):
            subset = mask & (sources == source)
            n, nr = int(subset.sum()), len(np.unique(groups[subset]))
            if n >= policy['minimum_scale_rows'] and nr >= policy['minimum_scale_regions']:
                value = region_quantile(ratio[subset], groups[subset])
                factor = max(factor, value)
                details.append({'source': int(source), 'selected_rows': n, 'selected_regions': nr, 'factor': value})
        scales[int(cell)] = {'supported': True, 'factor': max(factor, 1e-4), 'selected_rows': count,
                             'selected_regions': nregions, 'source_details': details}
    return {'scales': scales, 'selected_scale_rows': int(selected.sum()),
            'selected_scale_regions': len(np.unique(groups[selected])),
            'selection_used_calibration_labels': False,
            'absolute_unconditional_error_floor_used': False,
            'claim': 'Empirical region-balanced selected-residual calibration; no finite-sample or new-source coverage guarantee.'}


def apply_selected(package, calibrated, features, baseline, cells, supported_cell,
                   supported_line, study_count, policy):
    chosen, raw = eligibility(package, features, baseline, supported_cell, supported_line, study_count, policy)
    supported = np.array([calibrated['scales'].get(int(c), {}).get('supported', False) for c in cells])
    accepted = chosen & supported
    width = np.full(len(cells), np.nan, dtype='float64')
    for cell, settings in calibrated['scales'].items():
        if settings['supported']:
            mask = accepted & (cells == cell)
            width[mask] = raw[mask] * settings['factor']
    # No post-calibration interval-width selection: that would change the
    # selected population after its error distribution had been calibrated.
    assert np.isfinite(width[accepted]).all() and np.isnan(width[~accepted]).all()
    return {'precalibration_eligible': chosen, 'accepted': accepted, 'interval_halfwidth': width,
            'raw_error_scale': raw, 'unsupported_cell_calibration': ~supported}
